#!/usr/bin/env python3
"""Global aggregation helpers for comparaison_seq_struc_emb_new.

Builds two global files under data/GLOBAL:
- proteins_global.csv
- pairs_global.csv

Compact storage: no inline alignment blocks; keep file paths instead.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime
import os
from pathlib import Path
from typing import Dict, Tuple


POOLINGS = ["max", "bos", "mean", "sum", "mahalanobis"]
EMBED_CONDITIONS = [
    "raw",
    "normalized",
    "mean_outliers_filtered",
    "mean_outliers_only",
    "std_outliers_filtered",
    "std_outliers_only",
]

PAIR_PRIORITY_COLUMNS = [
    "mmseq2_pident",
    "mmseq2_qcov",
    "mmseq2_tcov",
    "mmseq2_alnlen",
    "mmseq2_evalue",
    "mmseq2_bits",
    "nw_pident",
    "nw_qcov",
    "nw_tcov",
    "nw_alnlen",
    "nw_score",
    "tm_score",
    "tm_score_reverse",
    "tm_rmsd",
    "tm_l_align",
    "tm_l_align_reverse",
    "tm_id_align_percent",
    "tm_id_align_percent_reverse",
    "foldseek_tm_score",
    "foldseek_alnlen",
    "foldseek_fident",
    "foldseek_evalue",
    "foldseek_bits",
    "foldseek_tm_score_reverse",
    "foldseek_alnlen_reverse",
    "foldseek_fident_reverse",
    "foldseek_evalue_reverse",
    "foldseek_bits_reverse",
]
PAIR_PRIORITY_COLUMNS.extend(f"emb_{pooling}_{condition}" for pooling in POOLINGS for condition in EMBED_CONDITIONS)


def _safe_read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _fasta_length(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        seq = "".join(line.strip() for line in lines if line and not line.startswith(">"))
        return len(seq) if seq else None
    except Exception:
        return None


def _to_float(value: str | float | int | None) -> float | None:
    if value is None:
        return None
    s = str(value).strip()
    if s == "" or s.lower() in {"nan", "none", "na", "n/a"}:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _to_int(value: str | float | int | None) -> int | None:
    f = _to_float(value)
    if f is None:
        return None
    try:
        return int(round(f))
    except Exception:
        return None


def _read_first_tsv_row(path: Path) -> dict | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return row
    return None


def _read_foldseek_row(path: Path, columns: list[str]) -> dict | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return None
    header = rows[0]
    if any(col in header for col in ("query", "target", "qlen", "tlen")):
        reader = csv.DictReader(["\t".join(row) for row in rows], delimiter="\t")
        for row in reader:
            return row
        return None
    row = rows[0]
    return {columns[i]: row[i] if i < len(row) else "" for i in range(len(columns))}


def _pair_key(query: str, target: str) -> Tuple[str, str]:
    return query.strip(), target.strip()


def _canonical_pair_key(query: str, target: str) -> Tuple[str, str]:
    pair = sorted((_pair_key(query, target)))
    return pair[0], pair[1]


def _is_nonempty(value: object) -> bool:
    return value not in (None, "")


def _row_priority(row: dict) -> int:
    return sum(1 for column in PAIR_PRIORITY_COLUMNS if _is_nonempty(row.get(column)))


def _merge_pair_orientations(orientation_rows: list[dict]) -> dict:
    primary = orientation_rows[0]
    row = _default_pair_row(primary.get("query", ""), primary.get("target", ""))
    sources = [row_source for row_source in orientation_rows if row_source is not None]

    def take(column: str, *sources: dict | None) -> None:
        for source in sources:
            if source is None:
                continue
            value = source.get(column, "")
            if _is_nonempty(value):
                row[column] = value
                return

    take("query", primary)
    take("target", primary)
    take("query_species", *sources)
    take("target_species", *sources)
    take("query_length", *sources)
    take("target_length", *sources)
    take("has_mmseq2", *sources)
    take("has_needleman", *sources)
    take("has_tmalign", *sources)
    take("has_foldseek", *sources)
    take("has_embedding", *sources)
    take("updated_at", *reversed(sources))

    for column in [
        "mmseq2_pident",
        "mmseq2_qcov",
        "mmseq2_tcov",
        "mmseq2_alnlen",
        "mmseq2_evalue",
        "mmseq2_bits",
        "mmseq2_alignment_file",
        "nw_pident",
        "nw_qcov",
        "nw_tcov",
        "nw_alnlen",
        "nw_score",
        "nw_alignment_file",
        # reverse Needleman fields intentionally omitted (symmetric)
        "tm_score",
        "tm_score_reverse",
        "tm_rmsd",
        "tm_l_align",
        "tm_l_align_reverse",
        "tm_id_align_percent",
        "tm_id_align_percent_reverse",
        "tm_alignment_file",
        "foldseek_tm_score",
        
        "foldseek_alnlen",
        "foldseek_fident",
        "foldseek_evalue",
        "foldseek_bits",
        "foldseek_tm_score_reverse",
        
        "foldseek_alnlen_reverse",
        "foldseek_fident_reverse",
        "foldseek_evalue_reverse",
        "foldseek_bits_reverse",
        "foldseek_alignment_file",
    ]:
        take(column, *sources)

    for pooling in POOLINGS:
        for condition in EMBED_CONDITIONS:
            take(f"emb_{pooling}_{condition}", *sources)

    return row


def _default_pair_row(query: str, target: str) -> dict:
    row = {
        "pair_id": f"{query}__{target}",
        "query": query,
        "target": target,
        "query_species": "",
        "target_species": "",
        "query_protein_name": "",
        "target_protein_name": "",
        "query_gene_names": "",
        "target_gene_names": "",
        "query_subunit_type": "",
        "target_subunit_type": "",
        "query_length": "",
        "target_length": "",
        "mmseq2_pident": "",
        "mmseq2_qcov": "",
        "mmseq2_tcov": "",
        "mmseq2_alnlen": "",
        "mmseq2_evalue": "",
        "mmseq2_bits": "",
        "mmseq2_alignment_file": "",
        "nw_pident": "",
        "nw_qcov": "",
        "nw_tcov": "",
        "nw_alnlen": "",
        "nw_score": "",
        "nw_alignment_file": "",
        "tm_score": "",
        "tm_score_reverse": "",
        "tm_rmsd": "",
        "tm_l_align": "",
        "tm_l_align_reverse": "",
        "tm_id_align_percent": "",
        "tm_id_align_percent_reverse": "",
        "tm_alignment_file": "",
        "foldseek_tm_score": "",
        "foldseek_alnlen": "",
        "foldseek_fident": "",
        "foldseek_evalue": "",
        "foldseek_bits": "",
        "foldseek_tm_score_reverse": "",
        "foldseek_alnlen_reverse": "",
        "foldseek_fident_reverse": "",
        "foldseek_evalue_reverse": "",
        "foldseek_bits_reverse": "",
        "foldseek_alignment_file": "",
        "has_mmseq2": "0",
        "has_needleman": "0",
        "has_tmalign": "0",
        "has_foldseek": "0",
        "has_embedding": "0",
        "updated_at": datetime.now().isoformat(),
    }
    for pooling in POOLINGS:
        for condition in EMBED_CONDITIONS:
            row[f"emb_{pooling}_{condition}"] = ""
    return row


def _build_proteins_global_csv(data_dir: Path, global_dir: Path) -> Path:
    inputs_dir = data_dir / "inputs"
    out_path = global_dir / "proteins_global.csv"

    rows = []
    for species_dir in sorted([p for p in inputs_dir.iterdir() if p.is_dir()]):
        species_key = species_dir.name
        for protein_dir in sorted([p for p in species_dir.iterdir() if p.is_dir()]):
            accession = protein_dir.name
            metadata = _safe_read_json(protein_dir / "metadata.json")
            fasta_path = protein_dir / "sequence.fasta"
            pdb_path = protein_dir / "structure.pdb"
            colabfold_path = protein_dir / "structure_colabfold.pdb"
            emb_path = protein_dir / "esm2_multipooling.pt"

            seq_len = _fasta_length(fasta_path)
            if seq_len is None:
                seq_len = _to_int(metadata.get("length"))

            annotations = metadata.get("annotations") or {}
            function_text = "; ".join(annotations.get("function", []) if isinstance(annotations.get("function"), list) else [])
            keywords = "; ".join(annotations.get("keywords", []) if isinstance(annotations.get("keywords"), list) else [])
            go_terms = "; ".join(annotations.get("go_terms", []) if isinstance(annotations.get("go_terms"), list) else [])

            rows.append(
                {
                    "accession": accession,
                    "species_key": species_key,
                    "species_label": metadata.get("species_label", ""),
                    "taxon_id": metadata.get("taxon_id", ""),
                    "uniprot_id": metadata.get("uniprot_id", ""),
                    "protein_name": metadata.get("protein_name", ""),
                    "gene_names": metadata.get("gene_names", ""),
                    "sequence_length": seq_len if seq_len is not None else "",
                    "mean_plddt": metadata.get("mean_plddt", ""),
                    "plddt_threshold": metadata.get("plddt_threshold", ""),
                    "subunit_type": metadata.get("subunit_type", ""),
                    "annotations_function": function_text,
                    "annotations_keywords": keywords,
                    "annotations_go_terms": go_terms,
                    "alphafold_pdb_url": metadata.get("alphafold_pdb_url", ""),
                    "has_sequence": "1" if fasta_path.exists() else "0",
                    "has_structure": "1" if pdb_path.exists() or colabfold_path.exists() else "0",
                    "has_embedding": "1" if emb_path.exists() else "0",
                    "sequence_fasta": str(fasta_path) if fasta_path.exists() else "",
                    "structure_pdb": str(pdb_path) if pdb_path.exists() else "",
                    "structure_colabfold_pdb": str(colabfold_path) if colabfold_path.exists() else "",
                    "embedding_file": str(emb_path) if emb_path.exists() else "",
                    "metadata_file": str(protein_dir / "metadata.json") if (protein_dir / "metadata.json").exists() else "",
                }
            )

    fieldnames = [
        "accession",
        "species_key",
        "species_label",
        "taxon_id",
        "uniprot_id",
        "protein_name",
        "gene_names",
        "sequence_length",
        "mean_plddt",
        "plddt_threshold",
        "subunit_type",
        "annotations_function",
        "annotations_keywords",
        "annotations_go_terms",
        "alphafold_pdb_url",
        "has_sequence",
        "has_structure",
        "has_embedding",
        "sequence_fasta",
        "structure_pdb",
        "structure_colabfold_pdb",
        "embedding_file",
        "metadata_file",
    ]

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return out_path


def _protein_lookup(proteins_csv_path: Path) -> Dict[str, dict]:
    lookup: Dict[str, dict] = {}
    if not proteins_csv_path.exists():
        return lookup
    with proteins_csv_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            accession = (row.get("accession") or "").strip()
            if accession:
                lookup[accession] = row
    return lookup

def _coverage_passes(qcov: float, tcov: float, and_threshold: float, or_threshold: float, mode: str) -> bool:
    if mode == "none":
        return True
    if mode not in {"and", "or", "both"}:
        raise ValueError(f"Unknown coverage mode: {mode}")
    ok = True
    if mode in {"and", "both"}:
        ok = ok and (qcov >= and_threshold and tcov >= and_threshold)
    if mode in {"or", "both"}:
        ok = ok and (qcov >= or_threshold or tcov >= or_threshold)
    return ok


def _load_mmseqs_pass_pairs(
    mmseqs_dir: Path,
    and_threshold: float,
    or_threshold: float,
    mode: str,
) -> set[Tuple[str, str]]:
    summary_path = mmseqs_dir / "pairwise_summary.csv"
    passed: set[Tuple[str, str]] = set()
    if not summary_path.exists():
        return passed

    with summary_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            query = (row.get("query") or "").strip()
            target = (row.get("target") or "").strip()
            if not query or not target:
                continue
            try:
                qcov = float(row.get("qcov") or 0)
                tcov = float(row.get("tcov") or 0)
            except ValueError:
                continue
            if _coverage_passes(qcov, tcov, and_threshold, or_threshold, mode):
                passed.add(_canonical_pair_key(query, target))
    return passed


def _build_pairs_global_csv(data_dir: Path, global_dir: Path, protein_by_accession: Dict[str, dict]) -> Path:
    out_path = global_dir / "pairs_global.csv"
    rows: Dict[Tuple[str, str], dict] = {}
    and_threshold = float(os.environ.get("MMSEQ2_COVERAGE_AND_THRESHOLD", "0.1"))
    or_threshold = float(os.environ.get("MMSEQ2_COVERAGE_OR_THRESHOLD", "0.1"))
    coverage_mode = os.environ.get("MMSEQ2_COVERAGE_MODE", "both")
    pass_pairs = _load_mmseqs_pass_pairs(data_dir / "alignment_mmseq2", and_threshold, or_threshold, coverage_mode)

    def get_row(query: str, target: str) -> dict:
        key = _pair_key(query, target)
        if key not in rows:
            row = _default_pair_row(query, target)
            q_meta = protein_by_accession.get(query, {})
            t_meta = protein_by_accession.get(target, {})
            row["query_species"] = q_meta.get("species_key", "")
            row["target_species"] = t_meta.get("species_key", "")
            row["query_protein_name"] = q_meta.get("protein_name", "")
            row["target_protein_name"] = t_meta.get("protein_name", "")
            row["query_gene_names"] = q_meta.get("gene_names", "")
            row["target_gene_names"] = t_meta.get("gene_names", "")
            row["query_subunit_type"] = q_meta.get("subunit_type", "")
            row["target_subunit_type"] = t_meta.get("subunit_type", "")
            row["query_length"] = q_meta.get("sequence_length", "")
            row["target_length"] = t_meta.get("sequence_length", "")
            rows[key] = row
        return rows[key]

    mmseq_dir = data_dir / "alignment_mmseq2"
    summary_path = mmseq_dir / "pairwise_summary.csv"
    if summary_path.exists():
        with summary_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row0 in reader:
                query = (row0.get("query") or "").strip()
                target = (row0.get("target") or "").strip()
                if not query or not target:
                    continue
                if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                    continue
                row = get_row(query, target)
                row["mmseq2_pident"] = row0.get("pident", "")
                row["mmseq2_qcov"] = row0.get("qcov", "")
                row["mmseq2_tcov"] = row0.get("tcov", "")
                row["mmseq2_alnlen"] = row0.get("alnlen", "")
                row["mmseq2_evalue"] = row0.get("evalue", "")
                row["mmseq2_bits"] = row0.get("bits", "")

                kept_alignment = (row0.get("kept_alignment") or "").strip().lower() == "yes"
                if kept_alignment:
                    pair_folder = (row0.get("pair_folder") or f"{query}-{target}").strip()
                    align_txt = mmseq_dir / pair_folder / "alignment.txt.gz"
                    if not align_txt.exists():
                        align_txt = mmseq_dir / pair_folder / "alignment.txt"
                    row["mmseq2_alignment_file"] = str(align_txt) if align_txt.exists() else ""
                    row["has_mmseq2"] = "1"
                else:
                    row["mmseq2_alignment_file"] = ""
                    row["has_mmseq2"] = "0"

                row["updated_at"] = datetime.now().isoformat()
    elif mmseq_dir.exists():
        for pair_dir in sorted([p for p in mmseq_dir.iterdir() if p.is_dir()]):
            row0 = _read_first_tsv_row(pair_dir / "alignment.tsv")
            if not row0:
                continue
            query = (row0.get("query") or "").strip()
            target = (row0.get("target") or "").strip()
            if not query or not target:
                continue
            if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                continue
            row = get_row(query, target)
            row["mmseq2_pident"] = row0.get("pident", "")
            row["mmseq2_qcov"] = row0.get("qcov", "")
            row["mmseq2_tcov"] = row0.get("tcov", "")
            row["mmseq2_alnlen"] = row0.get("alnlen", "")
            row["mmseq2_evalue"] = row0.get("evalue", "")
            row["mmseq2_bits"] = row0.get("bits", "")
            align_txt = pair_dir / "alignment.txt.gz"
            if not align_txt.exists():
                align_txt = pair_dir / "alignment.txt"
            row["mmseq2_alignment_file"] = str(align_txt) if align_txt.exists() else ""
            row["has_mmseq2"] = "1"
            row["updated_at"] = datetime.now().isoformat()

    nw_dir = data_dir / "alignment_needleman_wunsh"
    summary_path = nw_dir / "pairwise_summary.csv"
    if summary_path.exists():
        with summary_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row0 in reader:
                query = (row0.get("query") or "").strip()
                target = (row0.get("target") or "").strip()
                if not query or not target:
                    continue
                if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                    continue
                row = get_row(query, target)
                row["nw_pident"] = row0.get("pident", "")
                row["nw_qcov"] = row0.get("qcov", "")
                row["nw_tcov"] = row0.get("tcov", "")
                row["nw_alnlen"] = row0.get("alnlen", "")
                row["nw_score"] = row0.get("score", "")
                pair_folder = (row0.get("pair_folder") or f"{query}-{target}").strip()
                align_txt = nw_dir / pair_folder / "alignment.txt.gz"
                if not align_txt.exists():
                    align_txt = nw_dir / pair_folder / "alignment.txt"
                row["nw_alignment_file"] = str(align_txt) if align_txt.exists() else ""
                row["has_needleman"] = "1"
                row["updated_at"] = datetime.now().isoformat()
    elif nw_dir.exists():
        for pair_dir in sorted([p for p in nw_dir.iterdir() if p.is_dir()]):
            row0 = _read_first_tsv_row(pair_dir / "alignment.tsv")
            if not row0:
                continue
            query = (row0.get("query") or "").strip()
            target = (row0.get("target") or "").strip()
            if not query or not target:
                continue
            if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                continue
            row = get_row(query, target)
            row["nw_pident"] = row0.get("pident", "")
            row["nw_qcov"] = row0.get("qcov", "")
            row["nw_tcov"] = row0.get("tcov", "")
            row["nw_alnlen"] = row0.get("alnlen", "")
            row["nw_score"] = row0.get("score", "")
            align_txt = pair_dir / "alignment.txt.gz"
            if not align_txt.exists():
                align_txt = pair_dir / "alignment.txt"
            row["nw_alignment_file"] = str(align_txt) if align_txt.exists() else ""
            row["has_needleman"] = "1"
            row["updated_at"] = datetime.now().isoformat()

    tm_dir = data_dir / "alignment_structure_tmscore"
    foldseek_dir = data_dir / "alignment_structure_foldseek"
    summary_path = tm_dir / "pairwise_summary.csv"
    if summary_path.exists():
        with summary_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row0 in reader:
                query = (row0.get("query") or "").strip()
                target = (row0.get("target") or "").strip()
                if not query or not target:
                    continue
                if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                    continue
                row = get_row(query, target)
                row["tm_score"] = row0.get("tm_score", "")
                row["tm_score_reverse"] = row0.get("tm_score_reverse", "")
                row["tm_rmsd"] = row0.get("rmsd", "")
                row["tm_l_align"] = row0.get("l_align", "")
                row["tm_l_align_reverse"] = row0.get("l_align_reverse", "")
                row["tm_id_align_percent"] = row0.get("id_align_percent", "")
                row["tm_id_align_percent_reverse"] = row0.get("id_align_percent_reverse", "")
                pair_folder = (row0.get("pair_folder") or f"{query}-{target}").strip()
                align_txt = tm_dir / pair_folder / "alignment.txt.gz"
                if not align_txt.exists():
                    align_txt = tm_dir / pair_folder / "alignment.txt"
                row["tm_alignment_file"] = str(align_txt) if align_txt.exists() else ""
                row["has_tmalign"] = "1"
                row["updated_at"] = datetime.now().isoformat()
    elif tm_dir.exists():
        for pair_dir in sorted([p for p in tm_dir.iterdir() if p.is_dir()]):
            row0 = _read_first_tsv_row(pair_dir / "alignment.tsv")
            if not row0:
                continue
            query = (row0.get("query") or "").strip()
            target = (row0.get("target") or "").strip()
            if not query or not target:
                continue
            if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                continue
            row = get_row(query, target)
            row["tm_score"] = row0.get("tm_score", "")
            row["tm_score_reverse"] = row0.get("tm_score_reverse", "")
            row["tm_rmsd"] = row0.get("rmsd", "")
            row["tm_l_align"] = row0.get("l_align", "")
            row["tm_l_align_reverse"] = row0.get("l_align_reverse", "")
            row["tm_id_align_percent"] = row0.get("id_align_percent", "")
            row["tm_id_align_percent_reverse"] = row0.get("id_align_percent_reverse", "")
            align_txt = pair_dir / "alignment.txt.gz"
            if not align_txt.exists():
                align_txt = pair_dir / "alignment.txt"
            row["tm_alignment_file"] = str(align_txt) if align_txt.exists() else ""
            row["has_tmalign"] = "1"
            row["updated_at"] = datetime.now().isoformat()

    # Additionally ingest Foldseek results if present
    if foldseek_dir.exists():
        for pair_dir in sorted([p for p in foldseek_dir.iterdir() if p.is_dir()]):
            tsv_path = pair_dir / "alignment.tsv"
            tsv_path_reverse = pair_dir / "alignment_reverse.tsv"
            metadata = _safe_read_json(pair_dir / "metadata.json")
            query = (metadata.get("query_accession") or metadata.get("query") or "").strip()
            target = (metadata.get("target_accession") or metadata.get("target") or "").strip()
            if not query or not target:
                name_parts = pair_dir.name.split("-", 1)
                if len(name_parts) == 2:
                    query, target = name_parts[0].strip(), name_parts[1].strip()
            if not query or not target:
                continue
            if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                continue
            row = get_row(query, target)

            columns = metadata.get("alignment_columns") or []
            row0 = _read_foldseek_row(tsv_path, columns) if tsv_path.exists() else None
            if row0:
                row["foldseek_tm_score"] = (
                    row0.get("alntmscore") or row0.get("aln_tmscore") or row0.get("tmscore") or row0.get("tm_score", "")
                )
                row["foldseek_alnlen"] = row0.get("alnlen", "")
                row["foldseek_fident"] = row0.get("fident", "")
                row["foldseek_evalue"] = row0.get("evalue", "")
                row["foldseek_bits"] = row0.get("bits", "")

            row_rev = _read_foldseek_row(tsv_path_reverse, columns) if tsv_path_reverse.exists() else None
            if row_rev:
                row["foldseek_tm_score_reverse"] = (
                    row_rev.get("alntmscore") or row_rev.get("aln_tmscore") or row_rev.get("tmscore") or row_rev.get("tm_score", "")
                )
                row["foldseek_alnlen_reverse"] = row_rev.get("alnlen", "")
                row["foldseek_fident_reverse"] = row_rev.get("fident", "")
                row["foldseek_evalue_reverse"] = row_rev.get("evalue", "")
                row["foldseek_bits_reverse"] = row_rev.get("bits", "")

            align_txt = pair_dir / "alignment.txt"
            if align_txt.exists():
                row["foldseek_alignment_file"] = str(align_txt)
                if row0 or row_rev:
                    row["has_foldseek"] = "1"
            row["updated_at"] = datetime.now().isoformat()

    emb_dir = data_dir / "embedding_similarity"
    summary_path = emb_dir / "pairwise_summary.csv"
    if summary_path.exists():
        with summary_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row0 in reader:
                query = (row0.get("query") or "").strip()
                target = (row0.get("target") or "").strip()
                pooling = (row0.get("pooling") or "").strip().lower()
                condition = (row0.get("condition") or "").strip()
                if not query or not target:
                    continue
                if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                    continue
                row = get_row(query, target)
                if pooling in POOLINGS and condition in EMBED_CONDITIONS:
                    col = f"emb_{pooling}_{condition}"
                    row[col] = row0.get("cosine_similarity", "")
                    row["has_embedding"] = "1"
                    row["updated_at"] = datetime.now().isoformat()
    elif emb_dir.exists():
        for pair_dir in sorted([p for p in emb_dir.iterdir() if p.is_dir()]):
            tsv_path = pair_dir / "similarity.tsv"
            if not tsv_path.exists():
                continue
            with tsv_path.open("r", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row0 in reader:
                    query = (row0.get("query") or "").strip()
                    target = (row0.get("target") or "").strip()
                    pooling = (row0.get("pooling") or "").strip().lower()
                    condition = (row0.get("condition") or "").strip()
                    if not query or not target:
                        continue
                    if pass_pairs and _canonical_pair_key(query, target) not in pass_pairs:
                        continue
                    row = get_row(query, target)
                    if pooling in POOLINGS and condition in EMBED_CONDITIONS:
                        col = f"emb_{pooling}_{condition}"
                        row[col] = row0.get("cosine_similarity", "")
                        row["has_embedding"] = "1"
                        row["updated_at"] = datetime.now().isoformat()

    canonical_rows: Dict[Tuple[str, str], Dict[Tuple[str, str], dict]] = {}
    for (query, target), row in rows.items():
        canonical_key = _canonical_pair_key(query, target)
        bucket = canonical_rows.setdefault(canonical_key, {})
        bucket[(query, target)] = row

    merged_rows: list[dict] = []
    for canonical_key in sorted(canonical_rows.keys()):
        orientation_rows = canonical_rows[canonical_key]
        if not orientation_rows:
            continue
        ordered_orientations = sorted(
            orientation_rows.items(),
            key=lambda item: (-_row_priority(item[1]), item[0][0], item[0][1]),
        )
        primary_row = ordered_orientations[0][1]
        merged = _merge_pair_orientations([item[1] for item in ordered_orientations])
        
        # Only include pairs that have at least one alignment type completed
        has_any_alignment = (
            merged.get("has_mmseq2") == "1" or
            merged.get("has_needleman") == "1" or
            merged.get("has_tmalign") == "1" or
            merged.get("has_foldseek") == "1" or
            merged.get("has_embedding") == "1"
        )
        if has_any_alignment:
            merged_rows.append(merged)

    fieldnames = list(_default_pair_row("QUERY", "TARGET").keys())
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in merged_rows:
            writer.writerow(row)

    return out_path


def update_global_files(root_dir: Path) -> dict:
    root_dir = Path(root_dir)
    data_dir = root_dir / "data"
    global_dir = data_dir / "GLOBAL"
    global_dir.mkdir(parents=True, exist_ok=True)

    proteins_csv = _build_proteins_global_csv(data_dir, global_dir)
    protein_by_accession = _protein_lookup(proteins_csv)
    pairs_csv = _build_pairs_global_csv(data_dir, global_dir, protein_by_accession)

    return {
        "proteins_csv": str(proteins_csv),
        "pairs_csv": str(pairs_csv),
        "global_dir": str(global_dir),
        "updated_at": datetime.now().isoformat(),
    }


def _coerce_numeric_values(rows: list[dict]) -> list[dict]:
    numeric_columns = {
        "sequence_length",
        "mean_plddt",
        "plddt_threshold",
        "taxon_id",
        "mmseq2_pident",
        "mmseq2_qcov",
        "mmseq2_tcov",
        "mmseq2_alnlen",
        "mmseq2_evalue",
        "mmseq2_bits",
        "nw_pident",
        "nw_qcov",
        "nw_tcov",
        # Needleman reverse numeric columns omitted (symmetric)
        "nw_score",
        "tm_score",
        "tm_score_reverse",
        "tm_rmsd",
        "tm_l_align",
        "tm_l_align_reverse",
        "tm_id_align_percent",
        "tm_id_align_percent_reverse",
        "foldseek_tm_score",
        "foldseek_alnlen",
        "foldseek_fident",
        "foldseek_evalue",
        "foldseek_bits",
        "foldseek_tm_score_reverse",
        "foldseek_alnlen_reverse",
        "foldseek_fident_reverse",
        "foldseek_evalue_reverse",
        "foldseek_bits_reverse",
        "has_mmseq2",
        "has_needleman",
        "has_tmalign",
        "has_foldseek",
        "has_embedding",
    }

    for pooling in POOLINGS:
        for condition in EMBED_CONDITIONS:
            numeric_columns.add(f"emb_{pooling}_{condition}")

    result = []
    for row in rows:
        coerced_row = {}
        for key, value in row.items():
            if key in numeric_columns and value:
                try:
                    coerced_row[key] = float(value)
                except (ValueError, TypeError):
                    coerced_row[key] = value
            else:
                coerced_row[key] = value
        result.append(coerced_row)
    return result


def export_ui_data_js(root_dir: Path) -> Path:
    root_dir = Path(root_dir)
    data_dir = root_dir / "data"
    global_dir = data_dir / "GLOBAL"
    ui_dir = root_dir / "ui"
    ui_dir.mkdir(parents=True, exist_ok=True)

    proteins_csv = global_dir / "proteins_global.csv"
    pairs_csv = global_dir / "pairs_global.csv"
    out_js = ui_dir / "app_data.js"

    proteins_rows = []
    if proteins_csv.exists():
        with proteins_csv.open("r", encoding="utf-8") as handle:
            proteins_rows = list(csv.DictReader(handle))
        proteins_rows = _coerce_numeric_values(proteins_rows)

    pairs_rows = []
    if pairs_csv.exists():
        with pairs_csv.open("r", encoding="utf-8") as handle:
            pairs_rows = list(csv.DictReader(handle))
        pairs_rows = _coerce_numeric_values(pairs_rows)

    payload = {
        "metadata": {
            "proteins_count": len(proteins_rows),
            "pairs_count": len(pairs_rows),
            "source_proteins_csv": str(proteins_csv),
            "source_pairs_csv": str(pairs_csv),
            "mmseq2_coverage_mode": os.environ.get("MMSEQ2_COVERAGE_MODE", "both"),
            "mmseq2_coverage_and_threshold": os.environ.get("MMSEQ2_COVERAGE_AND_THRESHOLD", "0.1"),
            "mmseq2_coverage_or_threshold": os.environ.get("MMSEQ2_COVERAGE_OR_THRESHOLD", "0.1"),
            "updated_at": datetime.now().isoformat(),
        },
        "proteins": proteins_rows,
        "pairs": pairs_rows,
    }

    out_js.write_text("window.APP_DATA = " + json.dumps(payload, ensure_ascii=False) + ";\n", encoding="utf-8")
    return out_js
