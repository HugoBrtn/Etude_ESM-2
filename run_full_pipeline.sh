#!/bin/bash
# Pipeline execution with GPU support and 32 threads for MMseqs2
# Run this script or copy-paste the commands below one by one

set -e

cd ~/Bureau/Detecting-genomic-modules-across-multiple-pangenomes-using-MultiGraph-Neural-Network/Code/comparaison_seq_struc_emb_new

echo "=========================================="
echo "Pipeline Execution with GPU + 32 threads"
echo "=========================================="
echo ""

# Ensure environment is active
source /home/hugo/miniconda3/etc/profile.d/conda.sh
conda activate comparaison_emb

echo "[INFO] Environment: comparaison_emb"
echo "[INFO] GPU: enabled"
echo "[INFO] Threads: 32 (MMseqs2)"
echo ""

# Step 3: ESM-2 Embeddings (GPU)
echo "========================================"
echo "STEP 3: ESM-2 Embeddings (GPU)"
echo "========================================"
bash scripts/run_step_3_embeddings.sh --device cuda --batch-size 32 --debug
echo "[OK] Step 3 complete"
echo ""

# Step 4: MMseqs2 Alignments (GPU, 32 threads, DB mode)
echo "========================================"
echo "STEP 4: MMseqs2 Alignments (GPU, 32 threads)"
echo "========================================"
bash scripts/run_step_4_mmseqs2.sh --use-gpu --threads-mmseq2 32 --mmseqs-use-db --debug
echo "[OK] Step 4 complete"
echo ""

# Step 5: Needleman-Wunsch Alignments
echo "========================================"
echo "STEP 5: Needleman-Wunsch Alignments"
echo "========================================"
bash scripts/run_step_5_needleman.sh --debug
echo "[OK] Step 5 complete"
echo ""

# Step 6: TM-align Structure Alignment
echo "========================================"
echo "STEP 6: TM-align Structure Alignment"
echo "========================================"
bash scripts/run_step_6_tmalign.sh --debug
echo "[OK] Step 6 complete"
echo ""

# Step 7: Foldseek Structure Search
echo "========================================"
echo "STEP 7: Foldseek Structure Search (32 threads)"
echo "========================================"
bash scripts/run_step_7_foldseek.sh --threads-foldseek 32 --debug
echo "[OK] Step 7 complete"
echo ""

# Step 8: Compute Similarities
echo "========================================"
echo "STEP 8: Compute Similarities"
echo "========================================"
bash scripts/run_step_8_similarities.sh --debug
echo "[OK] Step 8 complete"
echo ""

# Step 9: Build Global CSV Outputs
echo "========================================"
echo "STEP 9: Build Global CSV Outputs"
echo "========================================"
bash scripts/run_step_9_global_csvs.sh --debug
echo "[OK] Step 9 complete"
echo ""

# Step 10: Bundle UI Data
echo "========================================"
echo "STEP 10: Bundle UI Data"
echo "========================================"
bash scripts/run_step_10_ui_bundle.sh --debug
echo "[OK] Step 10 complete"
echo ""

echo "=========================================="
echo "[SUCCESS] All steps complete!"
echo "=========================================="
echo ""
echo "Results:"
echo "  CSV outputs:    data/outputs/"
echo "  UI dashboard:   data/outputs/index.html"
echo "  Execution logs: data/.logs/"
echo ""
