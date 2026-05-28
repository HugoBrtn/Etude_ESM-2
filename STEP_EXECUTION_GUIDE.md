# Step-by-Step Pipeline Execution Guide

After running the initial `run_full_alphafold.sh` to collect proteins (which is the step 1-2 pipeline), you can execute the remaining analysis steps individually using the provided wrapper scripts.

## Quick Start

Once you've completed the collection phase:

```bash
cd ~/comparaison_seq_struc_emb_new

# Activate the environment
conda activate comparaison_emb

# Run steps 3-10 sequentially
bash scripts/run_step_3_embeddings.sh --debug
bash scripts/run_step_4_mmseqs2.sh --debug
bash scripts/run_step_5_needleman.sh --debug
bash scripts/run_step_6_tmalign.sh --debug
bash scripts/run_step_7_foldseek.sh --debug
bash scripts/run_step_8_similarities.sh --debug
bash scripts/run_step_9_global_csvs.sh --debug
bash scripts/run_step_10_ui_bundle.sh --debug
```

## Step-by-Step Details

### Prerequisites
- Conda environment activated: `conda activate comparaison_emb`
- Collection already completed (proteins present in `data/inputs/`)
- `pipeline_config.env` configured in `scripts/` directory

### Step 3: ESM-2 Embeddings

Compute ESM-2 (650M) protein embeddings for all collected sequences.

**Usage:**
```bash
bash scripts/run_step_3_embeddings.sh [options]

Options:
  --device cuda|cpu     Compute device (default: cuda)
  --batch-size <int>    Batch size (default: 32, reduce if CUDA OOM)
  --background          Run in background with nohup
  --debug               Verbose logging
  -h, --help            Show this help
```

**Examples:**
```bash
# Interactive with debug output
bash scripts/run_step_3_embeddings.sh --debug

# Background execution (check log with: tail -f data/.logs/embeddings.nohup.log)
bash scripts/run_step_3_embeddings.sh --background --debug

# CPU-only mode (slow, ~12-24h for 8500 proteins)
bash scripts/run_step_3_embeddings.sh --device cpu --debug

# GPU with reduced batch for limited VRAM
bash scripts/run_step_3_embeddings.sh --batch-size 8 --debug
```

**Output:**
- `data/inputs/esm_embeddings.pt` — Cached embeddings (updated incrementally)
- `data/.logs/step_3_embeddings.log` — Execution log

**Expected Duration:** 30min–2h (GPU) | 12–24h (CPU)

---

### Step 4: MMseqs2 Alignments

Fast sequence-based homology searches using the MMseqs2 database workflow.

**Usage:**
```bash
bash scripts/run_step_4_mmseqs2.sh [options]

Options:
  --threads-mmseq2 <int>         MMseqs2 threads (default: 32)
  --mmseqs-use-db                Use DB workflow (default: true, recommended)
  --mmseqs-no-db                 Use pairwise per-protein mode (slow)
  --ultra-fast                   Aggressive settings for huge datasets
  --use-gpu                      Enable GPU mode if available
  --mmseq2-coverage-and <f>      AND threshold (default: 0.1)
  --mmseq2-coverage-or <f>       OR threshold (default: 0.1)
  --mmseq2-coverage-mode <m>     one of: none, and, or, both (default: both)
  --background                   Run in background
  --debug                        Verbose logging
```

**Examples:**
```bash
# Default DB workflow (recommended)
bash scripts/run_step_4_mmseqs2.sh --debug

# Ultra-fast on huge dataset
bash scripts/run_step_4_mmseqs2.sh --ultra-fast --use-gpu --debug

# Background with custom coverage thresholds
bash scripts/run_step_4_mmseqs2.sh \
  --mmseq2-coverage-and 0.3 \
  --mmseq2-coverage-or 0.2 \
  --background --debug
```

**Output:**
- `data/inputs/mmseqs2_hits.csv` — Alignment results
- `data/.logs/step_4_mmseqs2.log` — Execution log

**Expected Duration:** 30min–2h

---

### Step 5: Needleman-Wunsch Alignments

Global sequence alignment refinement using Needleman-Wunsch.

**Usage:**
```bash
bash scripts/run_step_5_needleman.sh [options]

Options:
  --needleman-coverage-and <f>   AND threshold (default: 0.1)
  --needleman-coverage-or <f>    OR threshold (default: 0.1)
  --needleman-coverage-mode <m>  one of: none, and, or, both (default: both)
  --background                   Run in background
  --debug                        Verbose logging
```

**Examples:**
```bash
# Standard execution
bash scripts/run_step_5_needleman.sh --debug

# Background with strict coverage
bash scripts/run_step_5_needleman.sh \
  --needleman-coverage-and 0.5 \
  --needleman-coverage-or 0.3 \
  --background --debug
```

**Output:**
- `data/inputs/needleman_alignments.csv` — Alignment results
- `data/.logs/step_5_needleman.log` — Execution log

**Expected Duration:** 1–4h (depends on coverage thresholds)

---

### Step 6: TM-align Structure Alignment

Structure-based comparison using TM-align (requires AlphaFold2 structures).

**Usage:**
```bash
bash scripts/run_step_6_tmalign.sh [options]

Options:
  --tmalign-coverage-and <f>   AND threshold (default: 0.1)
  --tmalign-coverage-or <f>    OR threshold (default: 0.1)
  --tmalign-coverage-mode <m>  one of: none, and, or, both (default: both)
  --background                 Run in background
  --debug                      Verbose logging
```

**Examples:**
```bash
bash scripts/run_step_6_tmalign.sh --debug

bash scripts/run_step_6_tmalign.sh \
  --tmalign-coverage-and 0.2 \
  --background --debug
```

**Output:**
- `data/inputs/tmalign_results.csv` — Structure alignment scores
- `data/.logs/step_6_tmalign.log` — Execution log

**Expected Duration:** 2–8h

---

### Step 7: Foldseek Searches

Fast structure search using Foldseek.

**Usage:**
```bash
bash scripts/run_step_7_foldseek.sh [options]

Options:
  --threads-foldseek <int>      Foldseek threads (default: 32)
  --foldseek-coverage-and <f>   AND threshold (default: 0.1)
  --foldseek-coverage-or <f>    OR threshold (default: 0.1)
  --foldseek-coverage-mode <m>  one of: none, and, or, both (default: both)
  --background                  Run in background
  --debug                       Verbose logging
```

**Examples:**
```bash
bash scripts/run_step_7_foldseek.sh --debug

bash scripts/run_step_7_foldseek.sh --threads-foldseek 16 --debug
```

**Output:**
- `data/inputs/foldseek_hits.csv` — Foldseek search results
- `data/.logs/step_7_foldseek.log` — Execution log

**Expected Duration:** 30min–2h

---

### Step 8: Compute Similarities

Merge and compute combined similarity scores from all alignment methods.

**Usage:**
```bash
bash scripts/run_step_8_similarities.sh [options]

Options:
  --normalize-by-length   Normalize similarity scores by protein length
  --background            Run in background
  --debug                 Verbose logging
```

**Examples:**
```bash
bash scripts/run_step_8_similarities.sh --debug

bash scripts/run_step_8_similarities.sh --normalize-by-length --debug
```

**Output:**
- `data/outputs/similarities.csv` — Merged similarity matrix
- `data/.logs/step_8_similarities.log` — Execution log

**Expected Duration:** 5–15min

---

### Step 9: Build Global CSV Outputs

Aggregate similarity data into comprehensive CSV files for analysis.

**Usage:**
```bash
bash scripts/run_step_9_global_csvs.sh [options]

Options:
  --background    Run in background
  --debug         Verbose logging
```

**Examples:**
```bash
bash scripts/run_step_9_global_csvs.sh --debug

bash scripts/run_step_9_global_csvs.sh --background --debug
```

**Output:**
- `data/outputs/global_*.csv` — Aggregated results tables
- `data/.logs/step_9_global_csvs.log` — Execution log

**Expected Duration:** 1–5min

---

### Step 10: Bundle UI Data

Generate JSON data for interactive visualization dashboard.

**Usage:**
```bash
bash scripts/run_step_10_ui_bundle.sh [options]

Options:
  --background    Run in background
  --debug         Verbose logging
```

**Examples:**
```bash
bash scripts/run_step_10_ui_bundle.sh --debug

bash scripts/run_step_10_ui_bundle.sh --background --debug
```

**Output:**
- `data/outputs/ui_data.json` — Data for visualization
- `data/outputs/index.html` — Interactive dashboard
- `data/.logs/step_10_ui_bundle.log` — Execution log

**Expected Duration:** <1min

---

## Running Steps 3-10 in Sequence

Create a simple script to run all remaining steps:

```bash
cat > run_all_steps.sh << 'EOF'
#!/bin/bash
set -e

echo "[INFO] Running steps 3-10 sequentially..."

for step in 3 4 5 6 7 8 9 10; do
  echo ""
  echo "========================================"
  echo "Running step $step..."
  echo "========================================"
  bash scripts/run_step_${step}_*.sh --debug || {
    echo "[ERROR] Step $step failed!"
    exit 1
  }
  echo "[OK] Step $step complete."
done

echo ""
echo "[SUCCESS] All steps complete!"
echo "[INFO] Results in: data/outputs/"
EOF

chmod +x run_all_steps.sh
bash run_all_steps.sh
```

## Running Steps in Background

To run multiple steps simultaneously (not recommended without sufficient cores/memory):

```bash
# Start all steps in background
bash scripts/run_step_3_embeddings.sh --background
bash scripts/run_step_4_mmseqs2.sh --background
bash scripts/run_step_5_needleman.sh --background
# ... etc

# Monitor all logs
tail -f data/.logs/*.nohup.log

# Kill all background processes if needed
pkill -f "run_step.*nohup"
```

## Monitoring and Debugging

### View Log Files
```bash
# View step logs
tail -f data/.logs/step_*.log

# View background process logs
tail -f data/.logs/*.nohup.log

# View all at once (tmux or similar recommended)
for f in data/.logs/step_*.log; do tail -f "$f" & done; wait
```

### Check Progress
```bash
# View progress files (updated during execution)
ls -lh data/inputs/progress_* 2>/dev/null || echo "No progress files yet"

# Monitor process
watch ps aux | grep python
```

### Kill a Running Step
```bash
# Kill by PID (stored in logs/*.pid files)
kill $(cat data/.logs/step_N.pid)

# Or kill by script name
pkill -f run_step_N
```

## Troubleshooting

### CUDA Out of Memory
```bash
# Reduce batch size
bash scripts/run_step_3_embeddings.sh --batch-size 8
```

### Slow Performance
- For steps 3-4: Use `--use-gpu` flags
- For step 5-7: Reduce coverage thresholds with `--coverage-mode and`
- Check available CPU cores: `nproc` (adjust `--threads-*` accordingly)

### Missing Dependencies
```bash
# Validate environment
conda activate comparaison_emb
python -c "import esm, mmseqs, foldseek; print('OK')"
```

### Corrupted Intermediate Files
```bash
# Clean and restart from scratch
rm -rf data/inputs/*.csv data/inputs/*.pt
bash run_full_alphafold.sh  # Re-collect
bash scripts/run_step_3_embeddings.sh  # Resume pipeline
```

## Directory Structure

```
~/comparaison_seq_struc_emb_new/
├── data/
│   ├── inputs/           # Collected sequences + intermediate alignments
│   ├── outputs/          # Final results
│   └── .logs/            # Execution logs
├── scripts/
│   ├── run_step_*.sh     # Individual step wrappers (this guide)
│   ├── run_full_alphafold.sh  # Initial collection pipeline
│   ├── pipeline_config.env    # Configuration
│   ├── *.py              # Core pipeline scripts
│   └── ...
└── index.html            # UI dashboard (after step 10)
```

## Next Steps

After completing step 10:
1. Open `index.html` in a browser for interactive visualization
2. Explore `data/outputs/` for CSV results
3. Modify visualization parameters in UI and re-run steps as needed

---

**Questions?** Check individual script help:
```bash
bash scripts/run_step_*.sh --help
```
