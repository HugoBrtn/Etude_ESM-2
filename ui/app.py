#!/usr/bin/env python3
"""Mini Flask server for UI and data access."""
from __future__ import annotations

from pathlib import Path
import csv
import gzip
import json
import sys
import numpy as np
from flask import Flask, Response, jsonify, request, send_file, send_from_directory, abort
from scipy.stats import pearsonr, spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess
import numpy.random as npr

ROOT = Path(__file__).parent
DATA_DIR = ROOT.parent / "data"
PAIRS_CSV = DATA_DIR / "GLOBAL" / "pairs_global.csv"
PROTEINS_CSV = DATA_DIR / "GLOBAL" / "proteins_global.csv"

app = Flask(__name__, static_folder=str(ROOT), static_url_path="")


@app.after_request
def disable_cache(response):
    if 'app_data.js' in request.path or request.path == '/':
        response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response.headers['Pragma'] = 'no-cache'
        response.headers['Expires'] = '0'
    return response


def _clean_xy(x: list, y: list) -> tuple[np.ndarray, np.ndarray]:
    """Return finite-only (x, y) arrays sorted by x."""
    x_arr = np.array(x, dtype=float)
    y_arr = np.array(y, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    x_c = x_arr[mask]
    y_c = y_arr[mask]
    order = np.argsort(x_c)
    return x_c[order], y_c[order]


def _lowess_on_sorted(x: np.ndarray, y: np.ndarray, frac: float, it: int) -> tuple[np.ndarray, np.ndarray]:
    """Run LOWESS on pre-sorted unique-x data and return (x_smooth, y_smooth)."""
    # Deduplicate x by averaging y values — prevents degenerate frac windows
    x_u, inv = np.unique(x, return_inverse=True)
    y_u = np.bincount(inv, weights=y) / np.maximum(np.bincount(inv), 1)

    if len(x_u) < 3:
        return x_u, y_u

    x_range = float(x_u[-1] - x_u[0]) if len(x_u) > 1 else 0.0
    delta = 0.01 * x_range if x_range > 0 else 0.0

    smoothed = lowess(
        endog=y_u,
        exog=x_u,
        frac=frac,
        it=it,
        delta=delta,
        return_sorted=True,
    )
    return smoothed[:, 0], smoothed[:, 1]


def compute_loess_and_pearson(
    x: list,
    y: list,
    frac: float = 0.3,
    iterations: int = 3,
) -> dict:
    x_c, y_c = _clean_xy(x, y)

    if len(x_c) < 3:
        return {
            "error": "Not enough valid data points",
            "loess_x": [], "loess_y": [],
            "pearson": None, "spearman": None,
            "n": int(len(x_c)),
        }

    try:
        r, _ = pearsonr(x_c, y_c)
    except Exception:
        r = np.nan
    try:
        rho, _ = spearmanr(x_c, y_c)
    except Exception:
        rho = np.nan

    try:
        lx, ly = _lowess_on_sorted(x_c, y_c, frac, iterations)
        if len(lx) < 2:
            raise ValueError("Too few unique x values")
    except Exception as exc:
        return {
            "error": f"LOESS failed: {exc}",
            "loess_x": [], "loess_y": [],
            "pearson": float(r) if np.isfinite(r) else None,
            "spearman": float(rho) if np.isfinite(rho) else None,
            "n": int(len(x_c)),
        }

    return {
        "loess_x": lx.tolist(),
        "loess_y": ly.tolist(),
        "pearson": float(r) if np.isfinite(r) else None,
        "spearman": float(rho) if np.isfinite(rho) else None,
        "n": int(len(x_c)),
        "error": None,
    }


def lowess_with_confidence_bounds(
    x: list,
    y: list,
    eval_x: list | None = None,
    N: int = 200,
    conf_interval: float = 0.95,
    frac: float = 0.3,
    iterations: int = 3,
) -> tuple[list, list, list]:
    """LOWESS + bootstrap confidence band.

    Fix vs. original
    ----------------
    The original called ``lowess(..., xvals=eval_x)`` inside the bootstrap
    loop.  When a bootstrap resample does not cover the full range of
    ``eval_x`` (very common with small N or edge points), statsmodels
    returns NaN for those positions.  ``_filter_finite_ci`` then strips those
    positions, but if NaN appears consistently at the same index across many
    bootstrap samples the reconstructed band has gaps or becomes empty.

    The fix: fit LOWESS on the bootstrap sample with ``return_sorted=True``
    to get the native (sorted x, smoothed y) pairs, then use ``np.interp``
    to evaluate at ``eval_x``.  ``np.interp`` never produces NaN — it clamps
    to boundary values for out-of-range queries, which is the correct
    behaviour for a confidence band (we keep the last known value rather
    than introducing a hole).
    """
    x_c, y_c = _clean_xy(x, y)
    n = len(x_c)

    if n < 3:
        return [], [], []

    if eval_x is None:
        eval_x = np.linspace(float(x_c[0]), float(x_c[-1]), max(50, n)).tolist()

    eval_arr = np.array(eval_x, dtype=float)

    # ---- Baseline LOWESS at eval_x -----------------------------------------
    lx, ly = _lowess_on_sorted(x_c, y_c, frac, iterations)
    # Interpolate baseline onto eval_arr (no NaN: np.interp clamps at edges)
    smoothed_eval = np.interp(eval_arr, lx, ly)

    # ---- Bootstrap -----------------------------------------------------------
    boot_matrix = np.empty((N, len(eval_arr)), dtype=float)

    for i in range(N):
        idx = npr.randint(0, n, n)
        bx, by = x_c[idx], y_c[idx]
        try:
            blx, bly = _lowess_on_sorted(bx, by, frac, iterations)
            if len(blx) < 2:
                # Degenerate sample: fall back to flat line at mean
                boot_matrix[i, :] = np.full(len(eval_arr), np.mean(by))
            else:
                # np.interp: never NaN, clamps at boundaries
                boot_matrix[i, :] = np.interp(eval_arr, blx, bly)
        except Exception:
            boot_matrix[i, :] = smoothed_eval  # neutral fallback

    # ---- Percentile bands ----------------------------------------------------
    alpha = (1.0 - conf_interval) / 2.0
    lo_pct = max(0.0, alpha) * 100.0
    hi_pct = min(1.0, 1.0 - alpha) * 100.0
    bottom = np.percentile(boot_matrix, lo_pct, axis=0)
    top    = np.percentile(boot_matrix, hi_pct, axis=0)

    return (
        smoothed_eval.tolist(),
        bottom.tolist(),
        top.tolist(),
    )


def _load_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _coerce_numeric_values(rows: list[dict]) -> list[dict]:
    numeric_keys: set[str] = set()
    for row in rows:
        for k, v in row.items():
            if v in (None, ""):
                continue
            try:
                float(v)
                numeric_keys.add(k)
            except (ValueError, TypeError):
                continue
    result = []
    for row in rows:
        coerced = {}
        for k, v in row.items():
            if k in numeric_keys and v not in (None, ""):
                try:
                    coerced[k] = float(v)
                except (ValueError, TypeError):
                    coerced[k] = v
            else:
                coerced[k] = v
        result.append(coerced)
    return result


def _safe_resolve(path_str: str) -> Path:
    path_obj = Path(path_str)
    if path_str.startswith("data/"):
        path = (ROOT.parent / path_str).resolve()
    else:
        path = path_obj.expanduser().resolve()
    data_root = DATA_DIR.resolve()
    try:
        path.relative_to(data_root)
    except ValueError:
        raise ValueError("Path outside data directory")
    return path


@app.route("/")
def index():
    return send_from_directory(str(ROOT), "index.html")


@app.route("/api/pairs_global", methods=["GET"])
def api_pairs_global():
    rows = _coerce_numeric_values(_load_csv(PAIRS_CSV))
    return jsonify({"pairs": rows})


def prediction_interval_bounds(
    x: list,
    y: list,
    eval_x: list,
    conf_interval: float = 0.95,
    frac: float = 0.3,
    iterations: int = 3,
) -> tuple[list, list, list]:
    """Prediction interval based on LOESS residuals.

    Fits LOESS, computes residuals, estimates local std via a second LOWESS
    pass on |residuals|, then constructs ±k*sigma(x) bands where k is the
    normal quantile for the desired coverage.  This reflects the actual spread
    of data around the curve, unlike a bootstrap CI which only captures
    uncertainty in the smooth mean.
    """
    from scipy.stats import norm

    x_c, y_c = _clean_xy(x, y)
    n = len(x_c)
    if n < 3:
        return [], [], []

    eval_arr = np.array(eval_x, dtype=float)

    # Baseline LOESS
    lx, ly = _lowess_on_sorted(x_c, y_c, frac, iterations)
    smoothed_eval = np.interp(eval_arr, lx, ly)

    # Residuals at data points
    fitted_at_data = np.interp(x_c, lx, ly)
    residuals = y_c - fitted_at_data

    # Estimate local spread: LOWESS on |residuals|
    abs_res = np.abs(residuals)
    try:
        rx, ry = _lowess_on_sorted(x_c, abs_res, min(frac * 1.5, 1.0), 1)
        local_sigma = np.interp(eval_arr, rx, ry)
        # Use MAD-based estimate: E[|X|] = sigma * sqrt(2/pi) for normal
        local_sigma = local_sigma / np.sqrt(2.0 / np.pi)
        # Floor at global residual std / 4 to avoid tiny bands in flat regions
        global_sigma = float(np.std(residuals)) if len(residuals) > 1 else 0.0
        local_sigma  = np.maximum(local_sigma, global_sigma / 4.0)
    except Exception:
        global_sigma = float(np.std(residuals)) if len(residuals) > 1 else 0.0
        local_sigma  = np.full(len(eval_arr), global_sigma)

    alpha = (1.0 - conf_interval) / 2.0
    k = norm.ppf(1.0 - alpha)

    bottom = smoothed_eval - k * local_sigma
    top    = smoothed_eval + k * local_sigma

    return smoothed_eval.tolist(), bottom.tolist(), top.tolist()


def combined_interval_bounds(
    x: list,
    y: list,
    eval_x: list,
    N: int = 200,
    conf_interval: float = 0.95,
    frac: float = 0.3,
    iterations: int = 3,
) -> tuple[list, list, list]:
    """Combined CI: quadratic sum of bootstrap band + prediction residuals.

    sqrt(boot_half_width^2 + pred_half_width^2)
    """
    smoothed_pred, bot_pred, top_pred = prediction_interval_bounds(
        x, y, eval_x, conf_interval, frac, iterations
    )
    smoothed_boot, bot_boot, top_boot = lowess_with_confidence_bounds(
        x, y, eval_x=eval_x, N=N, conf_interval=conf_interval, frac=frac, iterations=iterations
    )
    if not smoothed_pred or not smoothed_boot:
        return smoothed_pred, bot_pred, top_pred

    sp = np.array(smoothed_pred)
    hw_pred = sp - np.array(bot_pred)
    hw_boot = sp - np.array(bot_boot)
    hw_combined = np.sqrt(hw_pred**2 + hw_boot**2)
    return sp.tolist(), (sp - hw_combined).tolist(), (sp + hw_combined).tolist()


@app.route("/api/loess", methods=["POST"])
def api_loess():
    data = request.get_json() or {}
    x    = data.get("x", [])
    y    = data.get("y", [])
    frac = float(data.get("frac", 0.3))
    it   = int(data.get("iterations", 3))

    result = compute_loess_and_pearson(x, y, frac, it)

    do_conf = bool(data.get("confidence", False))
    if do_conf and result.get("loess_x"):
        try:
            conf_n     = int(data.get("conf_n", 200))
            raw_level  = float(data.get("conf_level", 0.95))
            conf_level = raw_level / 100.0 if raw_level > 1.0 else raw_level
            ci_method  = data.get("ci_method", "pred")  # pred | boot | combined

            eval_x = result["loess_x"]

            if ci_method == "boot":
                smoothed, bottom, top = lowess_with_confidence_bounds(
                    x, y, eval_x=eval_x, N=conf_n, conf_interval=conf_level, frac=frac, iterations=it,
                )
            elif ci_method == "combined":
                smoothed, bottom, top = combined_interval_bounds(
                    x, y, eval_x=eval_x, N=conf_n, conf_interval=conf_level, frac=frac, iterations=it,
                )
            else:  # "pred" — default, recommended
                smoothed, bottom, top = prediction_interval_bounds(
                    x, y, eval_x=eval_x, conf_interval=conf_level, frac=frac, iterations=it,
                )

            result["ci_x"]        = eval_x
            result["ci_smoothed"] = smoothed
            result["ci_bottom"]   = bottom
            result["ci_top"]      = top

        except Exception as exc:
            result["ci_error"] = str(exc)

    return jsonify(result)


@app.route("/api/proteins_global", methods=["GET"])
def api_proteins_global():
    rows = _coerce_numeric_values(_load_csv(PROTEINS_CSV))
    return jsonify({"proteins": rows})


@app.route("/api/alignment", methods=["GET"])
def api_alignment():
    path_str = request.args.get("path", "")
    if not path_str:
        abort(400, "Missing path")
    try:
        path = _safe_resolve(path_str)
    except ValueError:
        abort(403, "Invalid path")
    if not path.exists():
        abort(404, "Alignment not found")
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        return Response(text, mimetype="text/plain")
    return send_file(str(path))


@app.route("/api/pdb", methods=["GET"])
def api_pdb():
    path_str = request.args.get("path", "")
    if not path_str:
        abort(400, "Missing path")
    try:
        path = _safe_resolve(path_str)
    except ValueError:
        abort(403, "Invalid path")
    if not path.exists():
        abort(404, "PDB not found")
    return send_file(str(path))


@app.route("/api/fasta", methods=["GET"])
def api_fasta():
    path_str = request.args.get("path", "")
    if not path_str:
        abort(400, "Missing path")
    try:
        path = _safe_resolve(path_str)
    except ValueError:
        abort(403, "Invalid path")
    if not path.exists():
        abort(404, "FASTA not found")
    return send_file(str(path))


@app.errorhandler(404)
def not_found(e):
    return send_from_directory(str(ROOT), "index.html")


if __name__ == "__main__":
    print("[INFO] Starting comparaison_seq_struc_emb_new UI", file=sys.stderr, flush=True)
    print("[INFO] Launch command: conda run --no-capture-output -n etude_esm2 python -u ui/app.py", file=sys.stderr, flush=True)
    print("[INFO] Working directory: <repo>/ui", file=sys.stderr, flush=True)
    print("[INFO] Server URL: http://127.0.0.1:5000", file=sys.stderr, flush=True)
    app.run(debug=False, host="127.0.0.1", port=5000)