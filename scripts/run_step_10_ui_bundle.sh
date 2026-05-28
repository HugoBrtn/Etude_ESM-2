#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="comparaison_emb"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="$ROOT_DIR/scripts"
OUTPUTS_DIR="$ROOT_DIR/data/outputs"
LOG_DIR="$ROOT_DIR/data/.logs"
source "$SCRIPTS_DIR/pipeline_config.env"
source "$SCRIPTS_DIR/progress_utils.sh"

BACKGROUND_MODE=false
DEBUG_MODE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --background)
      BACKGROUND_MODE=true
      shift
      ;;
    --debug)
      DEBUG_MODE=true
      shift
      ;;
    -h|--help)
      cat <<'EOF'
Usage:
  bash run_step_10_ui_bundle.sh [options]

Options:
  --background    Run in background with nohup
  --debug         Verbose logging
  -h, --help      Show this help
EOF
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1"
      exit 1
      ;;
  esac
done

if [[ "$BACKGROUND_MODE" == true ]]; then
  mkdir -p "$LOG_DIR"
  MAIN_LOG="$LOG_DIR/ui_bundle.nohup.log"
  PID_FILE="$LOG_DIR/ui_bundle.pid"
  
  nohup "$0" \
    $(if [[ "$DEBUG_MODE" == true ]]; then echo "--debug"; fi) \
    > "$MAIN_LOG" 2>&1 &
  CHILD_PID=$!
  echo "$CHILD_PID" > "$PID_FILE"
  
  echo "[INFO] UI bundle step started in background (PID: $CHILD_PID)"
  echo "[INFO] Main log: $MAIN_LOG"
  echo "[INFO] Monitor:"
  echo "       tail -f $MAIN_LOG"
  echo "[INFO] Stop:"
  echo "       kill $CHILD_PID"
  exit 0
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" || -z "${CONDA_PREFIX:-}" ]]; then
  echo "[ERROR] Activate conda env $ENV_NAME first."
  exit 1
fi

PYTHON_BIN="$CONDA_PREFIX/bin/python"
if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] Python not found in $PYTHON_BIN"
  exit 1
fi

unset PYTHONPATH
unset PYTHONHOME
export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/step_10_ui_bundle.log"
PROGRESS_FILE="$LOG_DIR/progress/step_10_ui_bundle.progress"

echo "[STEP 10] Bundle UI data"
echo "[INFO] Output dir: $OUTPUTS_DIR"
echo "[INFO] Logging to: $LOG_FILE"
echo "[INFO] Progress file: $PROGRESS_FILE"

if [[ "$DEBUG_MODE" == true ]]; then
  echo "[INFO] Debug mode: ON (logs in $LOG_FILE)"
fi

run_with_progress "Build UI bundle" "$PROGRESS_FILE" "$LOG_FILE" "$DEBUG_MODE" \
  "$PYTHON_BIN" "$SCRIPTS_DIR/build_ui_data.py"

echo "[OK] UI data bundling complete. Check: $LOG_FILE"
