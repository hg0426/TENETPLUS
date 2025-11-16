#!/bin/bash
## specify 
# $1 for input *.csv file
# $2 for number of parallel jobs
# $3 for trajectory_file
# $4 for cell_select_file
# $5 for historyLength
# $6 for species
# $7 0: TENET_TF(Only RNA), 1: TENET_Plus(RNA + ATAC) all matrix, 2: TENET_Plus only rowTF colGN, 3: TENET_Plus only rowTF colPK,  4: TENET_Plus rowTF colGN+PK, 5:TENET_Plus only rowPeak(cis-peaksource), 6: peak->peak

set -euo pipefail

LOCAL_TE_CHUNK_SIZE=${LOCAL_TE_CHUNK_SIZE:-300}
LOCAL_TE_BUFFER_EDGES=${LOCAL_TE_BUFFER_EDGES:-2000}
# Export tuning for localTE chunking
LOCAL_TE_EXPORT_WORKERS=${LOCAL_TE_EXPORT_WORKERS:-0}      # 0/1 = single process; >1 = multiprocessing
LOCAL_TE_READ_BATCH_ROWS=${LOCAL_TE_READ_BATCH_ROWS:-8192}
LOCAL_TE_USE_THREADS=${LOCAL_TE_USE_THREADS:-on}
LOCAL_TE_VALUES_DTYPE=${LOCAL_TE_VALUES_DTYPE:-float16}
LOCAL_TE_MERGE_WORKERS=${LOCAL_TE_MERGE_WORKERS:-}
# Optional: export per-selector chunk dirs (payload filtered by selector files)
LOCAL_TE_SPLIT_EXPORT=${LOCAL_TE_SPLIT_EXPORT:-on}
LOCAL_TE_SPLIT_OUTPUT_DIR=${LOCAL_TE_SPLIT_OUTPUT_DIR:-local_te_split_chunks}

TENET_GENE_FILTER=${TENET_GENE_FILTER:-none}

# Global default for pair mode (overridable via interactive prompt or TENET_PAIR_MODE env)
DEFAULT_PAIR_MODE=${TENET_PAIR_MODE:-default}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
CODE_DIR="$REPO_ROOT/code"
INPUT_DIR="$REPO_ROOT/input"
OUTPUT_DIR="$REPO_ROOT/output"
export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"
PYTHON=${PYTHON:-python3}

get_config_list() {
    local attr="$1"
    "$PYTHON" - <<PY
import importlib
cfg = importlib.import_module("code.config_options")
vals = getattr(cfg, "${attr}", [])
print(", ".join(vals))
PY
}

trim_whitespace() {
    "$PYTHON" - <<'PY' "$1"
import sys
print(sys.argv[1].strip())
PY
}

abs_path() {
    "$PYTHON" - <<'PY' "$1"
import os, sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PY
}

prompt_with_default() {
    local label="$1"
    local default_value="$2"
    local input
    if [ -n "$default_value" ]; then
        read -e -r -i "$default_value" -p "$label [$default_value]> " input
    else
        read -e -r -p "$label> " input
    fi
    input=$(trim_whitespace "${input:-$default_value}")
    echo "$input"
}

print_step() {
    local step="$1"
    local title="$2"
    shift 2
    echo ""
    echo "=== Step $step: $title ==="
    while [ $# -gt 0 ]; do
        echo "  • $1"
        shift
    done
}

prompt_input_matrix() {
    local default="${INPUT_FILE:-input/}"
    while true; do
        print_step "1" "Select Input Matrix" \
            "Expected: TE-ready expression matrix" \
            "Typical location: input/" \
            "Supports CSV/TSV/TXT/Parquet"
        local value=$(prompt_with_default "Path" "$default")
        if [ -n "$value" ]; then
            INPUT_FILE="$value"
            break
        fi
        echo "[!] Input matrix file path is required."
        default="$value"
    done
}

prompt_gene_filter() {
    local default="${TENET_GENE_FILTER:-none}"
    print_step "1b" "Housekeeping Gene Filter" \
        "none      = keep all genes" \
        "ribo_mito = drop ribosomal (RPS/RPL/MRPS/MRPL) and mitochondrial (MT-) genes" \
        "Can also be set via TENET_GENE_FILTER env."
    TENET_GENE_FILTER=$(prompt_with_default "Gene filter (none/ribo_mito)" "$default")
}

prompt_jobs() {
    local suggested=$(command -v nproc >/dev/null && nproc || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)
    local default="${NUM_JOBS:-}" 
    while true; do
        print_step "2" "Set Parallel Jobs" \
            "Defines worker processes for TE screening" \
            "Use available CPU cores (e.g. 8, 16)" \
            "Typical on this machine: $suggested"
        local value=$(prompt_with_default "Jobs" "$default")
        if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -gt 0 ]; then
            NUM_JOBS="$value"
            break
        fi
        echo "[!] Please enter a positive integer for the number of jobs."
        default="$value"
    done
}

prompt_trajectory_file() {
    local default="${TRAJECTORY_FILE:-input/}"
    while true; do
        print_step "3" "Select Trajectory File" \
            "Filtered pseudotime ordering (one value per cell)" \
            "Ensure values are aligned with the matrix rows"
        local value=$(prompt_with_default "Path" "$default")
        if [ -n "$value" ]; then
            TRAJECTORY_FILE="$value"
            break
        fi
        echo "[!] Trajectory file path is required."
        default="$value"
    done
}

prompt_cell_select_file() {
    local default="${CELL_SELECT_FILE:-input/}"
    while true; do
        print_step "4" "Select Cell-Select File" \
            "Binary mask (1/0) indicating cells to include" \
            "Should match the order of the expression matrix"
        local value=$(prompt_with_default "Path" "$default")
        if [ -n "$value" ]; then
            CELL_SELECT_FILE="$value"
            break
        fi
        echo "[!] Cell select file path is required."
        default="$value"
    done
}

prompt_history_length() {
    local default="${HISTORY_LENGTH:-1}"
    while true; do
        print_step "5" "Choose History Length k" \
            "Number of past timepoints to consider" \
            "Common choices: 1 or 2" \
            "Default: 1 (or type 'auto' to select based on series length)"
        local value
        value=$(prompt_with_default "k (or 'auto')" "$default")
        if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -gt 0 ]; then
            HISTORY_LENGTH="$value"
            break
        fi
        if [[ "${value,,}" == "auto" ]]; then
            HISTORY_LENGTH="auto"
            break
        fi
        echo "[!] Please enter a positive integer or 'auto' for history length."
        default="$value"
    done
}

prompt_pair_mode() {
    local default="${PAIR_MODE:-$DEFAULT_PAIR_MODE}"
    print_step "7-1" "Pair Set" \
        "default  = TF/peak-based pairs (original behaviour; honours mode 0-6)" \
        "all_pair = use full all-pairs TE (select variant in step 7-2)" \
        "Warning: all_pair can create O(F^2) edges and be very large."
    PAIR_MODE=$(prompt_with_default "Pair mode (default/all_pair)" "$default")
}

prompt_species() {
    local default="${SPECIES:-}" 
    while true; do
        print_step "6" "Declare Species" \
            "Used to resolve TF lists (e.g. human, mouse)"
        local value=$(prompt_with_default "Species" "$default")
        if [ -n "$value" ]; then
            SPECIES="$value"
            break
        fi
        echo "[!] Species is required."
        default="$value"
    done
}

prompt_mode() {
    while true; do
        if [ "${PAIR_MODE:-$DEFAULT_PAIR_MODE}" = "all_pair" ]; then
            print_step "7-2" "Select All-Pair TE Variant" \
                "0 = RNA-only, all gene×gene pairs" \
                "1 = TENET_Plus, all gene×gene pairs" \
                "2 = TENET_Plus, all feature×feature pairs (genes+peaks)"
            local default="${MODE_CODE:-}"
            local value
            value=$(prompt_with_default "All-pair mode (0-2)" "$default")
            if [[ "$value" =~ ^[0-2]$ ]]; then
                MODE_CODE="$value"
                break
            fi
            echo "[!] Please enter 0, 1, or 2 for all-pair mode."
        else
            print_step "7-2" "Select TENET Mode" \
                "0 = RNA-only (TENET_TF)" \
                "1-6 = TENET_Plus variations"
            "$PYTHON" - <<'PY'
import importlib
details = importlib.import_module("code.config_options").TENET_MODE_DETAILS
print("  • Mode legend:")
for key in sorted(details, key=int):
    info = details[key]
    print(f"    [{key}] {info['name']}")
    print(f"        {info['description']}")
PY
            local default="${MODE_CODE:-}" 
            local value
            value=$(prompt_with_default "Mode (0-6)" "$default")
            if [[ "$value" =~ ^[0-6]$ ]]; then
                MODE_CODE="$value"
                break
            fi
            echo "[!] Please enter an integer between 0 and 6."
        fi
    done
}

update_defaults_after_mode() {
    DEFAULT_MODALITY="none"
    DEFAULT_SCREEN_MODE="kernel"
    DEFAULT_REFINE_METHOD="none"
    DEFAULT_REFINE_TOPK="0"
    DEFAULT_REFINE_TOP_PCT="0"
    DEFAULT_PERMUTE="off"
    DEFAULT_PERM_N="100"
    DEFAULT_PERM_TOPK="0"
    DEFAULT_PERM_TOP_PCT="0"
    DEFAULT_PERM_FDR="off"
    DEFAULT_PERM_Q_ALPHA="0.05"
    DEFAULT_PERM_ALPHA="0.01"
    DEFAULT_LOCAL_TE="off"
    DEFAULT_LOCAL_EXPORT_WORKERS="${LOCAL_TE_EXPORT_WORKERS:-4}"
    DEFAULT_LOCAL_CHUNK_SIZE="$LOCAL_TE_CHUNK_SIZE"
    DEFAULT_LOCAL_BUFFER="$LOCAL_TE_BUFFER_EDGES"
    DEFAULT_LOCAL_READ_BATCH="$LOCAL_TE_READ_BATCH_ROWS"
    DEFAULT_LOCAL_THREADS="$LOCAL_TE_USE_THREADS"
    DEFAULT_LOCAL_DTYPE="$LOCAL_TE_VALUES_DTYPE"
    DEFAULT_LOCAL_MERGE_WORKERS="${LOCAL_TE_MERGE_WORKERS:-0}"
    DEFAULT_LOCAL_SPLIT_EXPORT="${LOCAL_TE_SPLIT_EXPORT:-on}"
    DEFAULT_LOCAL_SPLIT_DIR="${LOCAL_TE_SPLIT_OUTPUT_DIR:-local_te_split_chunks}"
    DEFAULT_RESULTS_BUFFER_ROWS="${TE_RESULTS_BUFFER_ROWS:-}"
    DEFAULT_INTERMEDIATE_SAVE="${TENET_INTERMEDIATE_SAVE:-on}"
    DEFAULT_BATCH_SIZE="${TE_BATCH_SIZE:-100}"
    DEFAULT_PAIR_MODE="${TENET_PAIR_MODE:-default}"
    DEFAULT_TIME_STRIDE="1"
    DEFAULT_TIME_PCT="100"
    DEFAULT_TIME_SEED="42"
}

prompt_modality() {
    local default="${MODALITY_CHOICE:-$DEFAULT_MODALITY}"
    print_step "8" "Choose Modality" \
        "Options: $MODALITY_MENU" \
        "Default: none (skip modality preprocessing)" \
        "Use rna or atac to enable modality-specific preprocessing" \
        "Current default: $default"
    MODALITY_CHOICE=$(prompt_with_default "Modality" "$default")
}

prompt_screen_mode() {
    local default="${SCREEN_MODE:-$DEFAULT_SCREEN_MODE}"
    print_step "9" "Screening Estimator" \
        "Options: $SCREEN_MENU" \
        "Default: kernel (fast, robust)"
    SCREEN_MODE=$(prompt_with_default "Screen mode" "$default")
}

prompt_refine_method() {
    local default="${REFINE_METHOD:-$DEFAULT_REFINE_METHOD}"
    print_step "10" "Refinement Strategy" \
        "Options: $REFINE_MENU" \
        "Set to none to keep screening scores"
    REFINE_METHOD=$(prompt_with_default "Refine method" "$default")
}

prompt_refine_thresholds() {
    if [ "${REFINE_METHOD,,}" = "none" ]; then
        REFINE_TOPK="$DEFAULT_REFINE_TOPK"
        REFINE_TOP_PCT="$DEFAULT_REFINE_TOP_PCT"
        return
    fi
    local default_k="${REFINE_TOPK:-$DEFAULT_REFINE_TOPK}"
    local default_pct="${REFINE_TOP_PCT:-$DEFAULT_REFINE_TOP_PCT}"
    print_step "11" "Refinement Thresholds" \
        "Top-K per target (0 = disable)" \
        "Percentile is applied after top-K filtering"
    REFINE_TOPK=$(prompt_with_default "Top-K" "$default_k")
    REFINE_TOP_PCT=$(prompt_with_default "Top percentile%" "$default_pct")
}

prompt_permutation_detail() {
    local default_n="${PERM_N:-$DEFAULT_PERM_N}"
    local default_topk="${PERM_TOPK:-$DEFAULT_PERM_TOPK}"
    local default_pct="${PERM_TOP_PCT:-$DEFAULT_PERM_TOP_PCT}"
    print_step "13" "Permutation Detail" \
        "Specify run size and target selection" \
        "Higher counts improve stability but add runtime"
    PERM_N=$(prompt_with_default "Perm count" "$default_n")
    PERM_TOPK=$(prompt_with_default "Perm top-K" "$default_topk")
    PERM_TOP_PCT=$(prompt_with_default "Perm top percentile%" "$default_pct")
}

prompt_permutation_fdr() {
    local default="${PERM_FDR:-$DEFAULT_PERM_FDR}"
    print_step "14" "Permutation FDR" \
        "Toggle: $PERM_MENU" \
        "on = control q-values (BH)" \
        "off = threshold by raw p-values"
    PERM_FDR=$(prompt_with_default "Perm FDR?" "$default")
}

prompt_permutation_thresholds() {
    local default_q="${PERM_Q_ALPHA:-$DEFAULT_PERM_Q_ALPHA}"
    local default_alpha="${PERM_ALPHA:-$DEFAULT_PERM_ALPHA}"
    print_step "15" "Permutation Thresholds" \
        "q-alpha applies when FDR is on" \
        "alpha applies to raw p-values when FDR is off"
    PERM_Q_ALPHA=$(prompt_with_default "q-alpha" "$default_q")
    PERM_ALPHA=$(prompt_with_default "p-alpha" "$default_alpha")
}

prompt_permutation_toggle() {
    local default="${PERMUTE_CHOICE:-$DEFAULT_PERMUTE}"
    print_step "12" "Permutation Testing" \
        "Toggle: $PERM_MENU" \
        "Enable for significance via permutations" \
        "Note: increases runtime considerably"
    PERMUTE_CHOICE=$(prompt_with_default "Permute?" "$default")
    if [ "${PERMUTE_CHOICE,,}" = "on" ]; then
        prompt_permutation_detail
        prompt_permutation_fdr
        prompt_permutation_thresholds
    else
        PERM_N="$DEFAULT_PERM_N"
        PERM_TOPK="$DEFAULT_PERM_TOPK"
        PERM_TOP_PCT="$DEFAULT_PERM_TOP_PCT"
        PERM_FDR="$DEFAULT_PERM_FDR"
        PERM_Q_ALPHA="$DEFAULT_PERM_Q_ALPHA"
        PERM_ALPHA="$DEFAULT_PERM_ALPHA"
    fi
}

prompt_local_te_details() {
    local default_chunk="${LOCAL_CHUNK_SIZE:-$DEFAULT_LOCAL_CHUNK_SIZE}"
    local default_buffer="${LOCAL_BUFFER_EDGES:-$DEFAULT_LOCAL_BUFFER}"
    local default_export="${LOCAL_EXPORT_WORKERS:-$DEFAULT_LOCAL_EXPORT_WORKERS}"
    local default_merge="${LOCAL_MERGE_WORKERS:-$DEFAULT_LOCAL_MERGE_WORKERS}"
    local default_read="${LOCAL_READ_BATCH:-$DEFAULT_LOCAL_READ_BATCH}"
    local default_threads="${LOCAL_THREADS_FLAG:-$DEFAULT_LOCAL_THREADS}"
    local default_dtype="${LOCAL_VALUES_DTYPE:-$DEFAULT_LOCAL_DTYPE}"
    local default_split="${LOCAL_SPLIT_EXPORT:-$DEFAULT_LOCAL_SPLIT_EXPORT}"
    local default_split_dir="${LOCAL_SPLIT_OUTPUT_DIR:-$DEFAULT_LOCAL_SPLIT_DIR}"
    print_step "17" "Local TE Parameters" \
        "Chunk size controls timepoints per parquet part" \
        "Workers tune parallel export/merge" \
        "Values dtype affects disk footprint"
    LOCAL_CHUNK_SIZE=$(prompt_with_default "Chunk size" "$default_chunk")
    LOCAL_BUFFER_EDGES=$(prompt_with_default "Buffer edges" "$default_buffer")
    LOCAL_EXPORT_WORKERS=$(prompt_with_default "Export workers" "$default_export")
    LOCAL_MERGE_WORKERS=$(prompt_with_default "Merge workers" "$default_merge")
    LOCAL_READ_BATCH=$(prompt_with_default "Read batch rows" "$default_read")
    LOCAL_THREADS_FLAG=$(prompt_with_default "Use threads (on/off)" "$default_threads")
    LOCAL_VALUES_DTYPE=$(prompt_with_default "Values dtype" "$default_dtype")
    LOCAL_SPLIT_EXPORT=$(prompt_with_default "Split export by selectors? (on/off)" "$default_split")
    LOCAL_SPLIT_OUTPUT_DIR=$(prompt_with_default "Split output directory" "$default_split_dir")
}

prompt_local_te_choice() {
    local default="${STORE_LOCAL_TE_CHOICE:-$DEFAULT_LOCAL_TE}"
    print_step "16" "Local TE Storage" \
        "Toggle: $LOCAL_MENU" \
        "on = store per-timepoint TE arrays (large output)"
    STORE_LOCAL_TE_CHOICE=$(prompt_with_default "Store local TE?" "$default")
    if [ "${STORE_LOCAL_TE_CHOICE,,}" = "on" ]; then
        prompt_local_te_details
    else
        LOCAL_CHUNK_SIZE="$DEFAULT_LOCAL_CHUNK_SIZE"
        LOCAL_BUFFER_EDGES="$DEFAULT_LOCAL_BUFFER"
        LOCAL_EXPORT_WORKERS="$DEFAULT_LOCAL_EXPORT_WORKERS"
        LOCAL_MERGE_WORKERS="$DEFAULT_LOCAL_MERGE_WORKERS"
        LOCAL_READ_BATCH="$DEFAULT_LOCAL_READ_BATCH"
        LOCAL_THREADS_FLAG="$DEFAULT_LOCAL_THREADS"
        LOCAL_VALUES_DTYPE="$DEFAULT_LOCAL_DTYPE"
    fi
}

prompt_time_options() {
    local default_stride="${TIME_STRIDE:-$DEFAULT_TIME_STRIDE}"
    local default_pct="${TIME_PCT:-$DEFAULT_TIME_PCT}"
    local default_seed="${TIME_SEED:-$DEFAULT_TIME_SEED}"
    print_step "18" "Time Subsampling" \
        "Stride > 1 = use every Nth timepoint" \
        "Otherwise use percentile sampling"
    TIME_STRIDE=$(prompt_with_default "Stride" "$default_stride")
    TIME_PCT=$(prompt_with_default "Time sampling %" "$default_pct")
    TIME_SEED=$(prompt_with_default "Time seed" "$default_seed")
}

prompt_results_buffer_rows() {
    local default="${RESULTS_BUFFER_ROWS:-$DEFAULT_RESULTS_BUFFER_ROWS}"
    local default_bs="${BATCH_SIZE:-$DEFAULT_BATCH_SIZE}"
    local default_is="${INTERMEDIATE_SAVE_CHOICE:-$DEFAULT_INTERMEDIATE_SAVE}"
    print_step "19" "TE Results Buffering" \
        "Intermediate save: on/off controls whether TE writes incremental Parquet batches for resume/restart" \
        "Rows to buffer before writing batch Parquet files (results_buffer_rows)" \
        "Leave results_buffer_rows blank for auto (200000, or 5000 when storing LocalTE)" \
        "Batch size controls how many sources per target are processed together"
    INTERMEDIATE_SAVE_CHOICE=$(prompt_with_default "Intermediate save (on/off)" "$default_is")
    RESULTS_BUFFER_ROWS=$(prompt_with_default "results_buffer_rows" "$default")
    BATCH_SIZE=$(prompt_with_default "TE batch size" "$default_bs")
}

show_summary() {
    echo ""
    echo "--- Parameter Summary ---"
    echo " 1) Input matrix         : $INPUT_FILE"
    echo " 1b) Gene filter         : ${TENET_GENE_FILTER:-none}"
    echo " 2) Parallel jobs        : $NUM_JOBS"
    echo " 3) Trajectory file      : $TRAJECTORY_FILE"
    echo " 4) Cell-select file     : $CELL_SELECT_FILE"
    echo " 5) History length (k)   : $HISTORY_LENGTH"
    echo " 6) Species              : $SPECIES"
    echo " 7-1) Pair mode          : ${PAIR_MODE:-$DEFAULT_PAIR_MODE}"
    if [ "${PAIR_MODE:-$DEFAULT_PAIR_MODE}" = "all_pair" ]; then
        case "${MODE_CODE:-0}" in
            0) echo " 7-2) All-pair mode      : 0 (RNA-only, all gene×gene)" ;;
            1) echo " 7-2) All-pair mode      : 1 (TENET_Plus, all gene×gene)" ;;
            2) echo " 7-2) All-pair mode      : 2 (TENET_Plus, all feature×feature)" ;;
            *) echo " 7-2) All-pair mode      : ${MODE_CODE:-0}" ;;
        esac
    else
        echo " 7-2) Mode               : $MODE_CODE"
    fi
    echo " 8) Modality             : $MODALITY_CHOICE"
    echo " 9) Screening estimator  : $SCREEN_MODE"
    echo "10) Refinement method    : $REFINE_METHOD"
    if [ "${REFINE_METHOD,,}" = "none" ]; then
        echo "11) Refinement thresholds: (not applicable)"
    else
        echo "11) Refinement thresholds: Top-K=$REFINE_TOPK, Percentile=$REFINE_TOP_PCT"
    fi
    echo "12) Permutation toggle   : $PERMUTE_CHOICE"
    if [ "${PERMUTE_CHOICE,,}" = "on" ]; then
        echo "13) Permutation detail   : count=$PERM_N, topK=$PERM_TOPK, pct=$PERM_TOP_PCT"
        echo "14) Permutation FDR      : $PERM_FDR"
        echo "15) Permutation thresholds: q=$PERM_Q_ALPHA, p=$PERM_ALPHA"
    else
        echo "13) Permutation detail   : (not applicable)"
        echo "14) Permutation FDR      : (not applicable)"
        echo "15) Permutation thresholds: (not applicable)"
    fi
    echo "16) Local TE storage     : $STORE_LOCAL_TE_CHOICE"
    if [ "${STORE_LOCAL_TE_CHOICE,,}" = "on" ]; then
        echo "17) Local TE parameters  : chunk=$LOCAL_CHUNK_SIZE, buffer=$LOCAL_BUFFER_EDGES, export=$LOCAL_EXPORT_WORKERS, merge=$LOCAL_MERGE_WORKERS, read=$LOCAL_READ_BATCH, threads=$LOCAL_THREADS_FLAG, dtype=$LOCAL_VALUES_DTYPE, split=$LOCAL_SPLIT_EXPORT, split_dir=$LOCAL_SPLIT_OUTPUT_DIR"
    else
        echo "17) Local TE parameters  : (not applicable)"
    fi
    if [ -n "${RESULTS_BUFFER_ROWS:-}" ]; then
        echo "19) Results buffer rows  : $RESULTS_BUFFER_ROWS"
    else
        echo "19) Results buffer rows  : <auto>"
    fi
    echo "    Intermediate save    : ${INTERMEDIATE_SAVE_CHOICE:-$DEFAULT_INTERMEDIATE_SAVE}"
    if [ -n "${BATCH_SIZE:-}" ]; then
        echo "    TE batch size        : $BATCH_SIZE"
    else
        echo "    TE batch size        : 100 (default)"
    fi
    echo "18) Time subsampling     : stride=$TIME_STRIDE, pct=$TIME_PCT, seed=$TIME_SEED"
    echo "---------------------------"
}

mkdir -p "$OUTPUT_DIR"

echo "Running memory monitor..."
chmod +x "$CODE_DIR/memory_check.sh" || true
CURRENT_STAGE_FILE="$OUTPUT_DIR/current_stage.txt"
echo "INIT" > "$CURRENT_STAGE_FILE"
# GPU monitor toggle (env override: MONITOR_GPU=on|off); default off for lower overhead
MONITOR_GPU="${MONITOR_GPU:-off}"
"$CODE_DIR/memory_check.sh" $$ "$CURRENT_STAGE_FILE" "$MONITOR_GPU" &
echo $$

# Helper to mark pipeline stage (also echoed for visibility)
set_stage() {
  local label="$1"
  echo "$label" | tee "$CURRENT_STAGE_FILE" >/dev/null
}

# Run a command while continuously stamping stage label so the monitor always sees it
run_with_stage() {
  local label="$1"; shift
  set_stage "$label"
  "$@" &
  local cmd_pid=$!
  (
    while kill -0 "$cmd_pid" >/dev/null 2>&1; do
      echo "$label" > "$CURRENT_STAGE_FILE"
      sleep 1
    done
  ) &
  local writer_pid=$!
  wait "$cmd_pid" || { kill "$writer_pid" >/dev/null 2>&1 || true; wait "$writer_pid" 2>/dev/null || true; return 1; }
  kill "$writer_pid" >/dev/null 2>&1 || true
  wait "$writer_pid" 2>/dev/null || true
}

build_replay_command() {
    local cmd=""
    local var val
    for var in PYTHON MONITOR_GPU TENET_PAIR_MODE TENET_GENE_FILTER TENET_INTERMEDIATE_SAVE TE_BATCH_SIZE LOCAL_CHUNK_SIZE LOCAL_BUFFER_EDGES LOCAL_EXPORT_WORKERS \
               LOCAL_MERGE_WORKERS LOCAL_READ_BATCH LOCAL_THREADS_FLAG LOCAL_VALUES_DTYPE \
               LOCAL_SPLIT_EXPORT LOCAL_SPLIT_OUTPUT_DIR RESULTS_BUFFER_ROWS; do
        val="${!var-}"
        if [ -n "${val:-}" ]; then
            cmd+="${var}=$(printf '%q' "$val") "
        fi
    done
    cmd+="\"$0\""
    for arg in "$@"; do
        cmd+=" $(printf '%q' "$arg")"
    done
    printf '%s\n' "$cmd"
}

write_replay_logs() {
    local cmd
    cmd=$(build_replay_command "$@")
    local txt="$OUTPUT_DIR/run_invocation.txt"
    local sh="$OUTPUT_DIR/run_invocation.sh"
    {
        echo "# Invocation (copy-paste to re-run):"
        echo "$cmd"
    } > "$txt"
    {
        echo "#!/usr/bin/env bash"
        echo "# Auto-generated by TENET_Plus_for_py.sh on $(date '+%Y-%m-%d %H:%M:%S')"
        echo ""
        echo "$cmd"
    } > "$sh"
    chmod +x "$sh" 2>/dev/null || true
    echo "Final command (replayable) logged to:"
    echo "  $txt"
    echo "  $sh"
    echo "Replay command:"
    echo "  $cmd"
}

# Function to display usage
usage() {
    cat <<USAGE
Usage: $0 <input_file> <num_jobs> <trajectory_file> <cell_select_file> <history_k> <species> <mode_code> [modality] [screen_mode] [refine_method] [refine_topk] [refine_top_pct] [permute] [perm_n] [perm_topk] [perm_top_pct] [perm_fdr] [perm_q_alpha] [perm_alpha] [store_local_te] [time_stride] [time_pct] [time_seed] [results_buffer_rows]

  <mode_code>: 0 = TENET_TF (RNA only)
               1 = TENET_Plus (RNA+ATAC all)
               2 = TENET_Plus only rowTF colGN
               3 = TENET_Plus only rowTF colPK
               4 = TENET_Plus rowTF colGN+PK
               5 = TENET_Plus only rowPeak (cis-peaksource)
               6 = TENET_Plus peak->peak (cis)

Optional:
  Pair mode (interactive step 7-1):
    default  = TF/peak-based pairs (above mode_code 0-6 semantics)
    all_pair = use full all-pairs TE, then:
                 mode_code 0 -> RNA-only, all gene×gene
                 mode_code 1 -> TENET_Plus, all gene×gene
                 mode_code 2 -> TENET_Plus, all feature×feature (genes+peaks)

  modality:    rna | atac | auto | none  (default auto for Plus modes, rna for TF mode; 'none' skips modality preprocessing)
  screen_mode: linear | poly | ksg | kernel | gcmi | disc | ordinal | kernel_grid  (default: kernel)
  refine_method: kernel | ksg | none        (default: none)
  refine_topk:  integer K per target        (default: 0)
  refine_top_pct:  percentile [0-100]       (default: 0)
  permute:     on | off                     (default: off)
  perm_n:      number of permutations       (default: 100)
  perm_topk:   K per target for perms       (default: 0 = all)
  perm_top_pct:percentile for perms         (default: 0)
  perm_fdr:    on | off                     (default: off)
  perm_q_alpha:FDR q threshold              (default: 0.05)
  perm_alpha:  p-value threshold            (default: 0.01)
  store_local_te: on | off                 (default: off)
  time_stride: sample every Nth timepoint   (default: 1 = no stride)
  time_pct:    randomly keep P% timepoints  (default: 100; used when stride==1)
  time_seed:   RNG seed for time_pct        (default: 42)
  results_buffer_rows: rows to buffer in TE before flushing batch parquet files (default auto: 200000, or 5000 when storing LocalTE)

Environment:
  TENET_GENE_FILTER   : none | ribo_mito (default: none). When 'ribo_mito',
                        drop genes whose names start with RPS, RPL, MRPS,
                        MRPL, or MT- from the input matrix before analysis.
USAGE
    exit 0
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
fi

if [ $# -eq 0 ]; then
    echo "Interactive mode: press Enter to accept defaults shown in brackets."
    echo "Path prompts start in the repository root ($REPO_ROOT); try Tab for completion."

    MODALITY_MENU=$(get_config_list MODALITY_CHOICES)
    SCREEN_MENU=$(get_config_list SCREEN_CHOICES)
    REFINE_MENU=$(get_config_list REFINE_CHOICES)
    PERM_MENU=$(get_config_list PERM_TOGGLE)
    LOCAL_MENU=$(get_config_list LOCAL_TE_TOGGLE)

    LOCAL_CHUNK_SIZE="${LOCAL_CHUNK_SIZE:-}"
    LOCAL_BUFFER_EDGES="${LOCAL_BUFFER_EDGES:-}"
    LOCAL_EXPORT_WORKERS="${LOCAL_EXPORT_WORKERS:-}"
    LOCAL_MERGE_WORKERS="${LOCAL_MERGE_WORKERS:-}"
    LOCAL_READ_BATCH="${LOCAL_READ_BATCH:-}"
    LOCAL_THREADS_FLAG="${LOCAL_THREADS_FLAG:-}"
    LOCAL_VALUES_DTYPE="${LOCAL_VALUES_DTYPE:-}"
    RESULTS_BUFFER_ROWS="${RESULTS_BUFFER_ROWS:-}"

    # Initial prompts
    prompt_input_matrix
    prompt_gene_filter
    prompt_jobs
    prompt_trajectory_file
    prompt_cell_select_file
    prompt_history_length
    prompt_species
    prompt_pair_mode
    prompt_mode
    update_defaults_after_mode
    prompt_modality
    prompt_screen_mode
    prompt_refine_method
    prompt_refine_thresholds
    prompt_permutation_toggle
    prompt_local_te_choice
    prompt_time_options
    prompt_results_buffer_rows

    while true; do
        show_summary
        read -rp "Edit step (e.g. 1, 6, 7-1, 7-2, 19) or press Enter to continue: " EDIT_CHOICE
        EDIT_CHOICE=$(trim_whitespace "$EDIT_CHOICE")
        case "$EDIT_CHOICE" in
            "")
                break
                ;;
            1)
                prompt_input_matrix
                ;;
            1b|1B)
                prompt_gene_filter
                ;;
            2)
                prompt_jobs
                ;;
            3)
                prompt_trajectory_file
                ;;
            4)
                prompt_cell_select_file
                ;;
            5)
                prompt_history_length
                ;;
            6)
                prompt_species
                ;;
            7|7-2|"7-2")
                prompt_mode
                update_defaults_after_mode
                ;;
            20|7-1|"7-1")
                prompt_pair_mode
                ;;
            8)
                prompt_modality
                ;;
            9)
                prompt_screen_mode
                ;;
            10)
                prompt_refine_method
                prompt_refine_thresholds
                ;;
            11)
                prompt_refine_thresholds
                ;;
            12)
                prompt_permutation_toggle
                ;;
            13)
                if [ "${PERMUTE_CHOICE,,}" = "on" ]; then
                    prompt_permutation_detail
                else
                    echo "[i] Permutation is off; nothing to edit for step 13."
                fi
                ;;
            14)
                if [ "${PERMUTE_CHOICE,,}" = "on" ]; then
                    prompt_permutation_fdr
                else
                    echo "[i] Permutation is off; nothing to edit for step 14."
                fi
                ;;
            15)
                if [ "${PERMUTE_CHOICE,,}" = "on" ]; then
                    prompt_permutation_thresholds
                else
                    echo "[i] Permutation is off; nothing to edit for step 15."
                fi
                ;;
            16)
                prompt_local_te_choice
                ;;
            17)
                if [ "${STORE_LOCAL_TE_CHOICE,,}" = "on" ]; then
                    prompt_local_te_details
                else
                    echo "[i] Local TE storage is off; nothing to edit for step 17."
                fi
                ;;
            18)
                prompt_time_options
                ;;
            19)
                prompt_results_buffer_rows
                ;;
            *)
                echo "[!] Unknown step selection."
                ;;
        esac
    done

    set -- "$INPUT_FILE" "$NUM_JOBS" "$TRAJECTORY_FILE" "$CELL_SELECT_FILE" \
        "$HISTORY_LENGTH" "$SPECIES" "$MODE_CODE" "$MODALITY_CHOICE" "$SCREEN_MODE" \
        "$REFINE_METHOD" "$REFINE_TOPK" "$REFINE_TOP_PCT" "$PERMUTE_CHOICE" "$PERM_N" \
        "$PERM_TOPK" "$PERM_TOP_PCT" "$PERM_FDR" "$PERM_Q_ALPHA" "$PERM_ALPHA" \
        "$STORE_LOCAL_TE_CHOICE" "$TIME_STRIDE" "$TIME_PCT" "$TIME_SEED" \
        "$RESULTS_BUFFER_ROWS"

    # After final confirmation, print summary once more and snapshot invocation
    show_summary "$@"
    cp -f "$0" "$OUTPUT_DIR/run_script_snapshot.sh" 2>/dev/null || true
    # Enable command trace to a file (keeps stdout clean)
    TRACE_LOG="$OUTPUT_DIR/command_trace.log"
    exec 3>"$TRACE_LOG"
    export BASH_XTRACEFD=3
    export PS4='+ $(date "+%F %T") [$LINENO] '
    set -x
    echo "Command trace logging to $TRACE_LOG"
fi

if [ -n "${LOCAL_CHUNK_SIZE:-}" ]; then
    LOCAL_TE_CHUNK_SIZE="$LOCAL_CHUNK_SIZE"
fi
if [ -n "${LOCAL_BUFFER_EDGES:-}" ]; then
    LOCAL_TE_BUFFER_EDGES="$LOCAL_BUFFER_EDGES"
fi
if [ -n "${LOCAL_EXPORT_WORKERS:-}" ]; then
    LOCAL_TE_EXPORT_WORKERS="$LOCAL_EXPORT_WORKERS"
fi
if [ -n "${LOCAL_MERGE_WORKERS:-}" ]; then
    LOCAL_TE_MERGE_WORKERS="$LOCAL_MERGE_WORKERS"
fi
if [ -n "${LOCAL_READ_BATCH:-}" ]; then
    LOCAL_TE_READ_BATCH_ROWS="$LOCAL_READ_BATCH"
fi
if [ -n "${LOCAL_THREADS_FLAG:-}" ]; then
    LOCAL_TE_USE_THREADS="$LOCAL_THREADS_FLAG"
fi
if [ -n "${LOCAL_VALUES_DTYPE:-}" ]; then
    LOCAL_TE_VALUES_DTYPE="$LOCAL_VALUES_DTYPE"
fi
if [ -n "${LOCAL_SPLIT_EXPORT:-}" ]; then
    LOCAL_TE_SPLIT_EXPORT="$LOCAL_SPLIT_EXPORT"
fi
if [ -n "${LOCAL_SPLIT_OUTPUT_DIR:-}" ]; then
    LOCAL_TE_SPLIT_OUTPUT_DIR="$LOCAL_SPLIT_OUTPUT_DIR"
fi

# Check if at least 7 arguments are passed
if [ $# -lt 7 ]; then
    echo "Insufficient arguments supplied."
    usage
fi

ARG1_TRIM=$(trim_whitespace "$1")
ARG2_TRIM=$(trim_whitespace "$2")
ARG3_TRIM=$(trim_whitespace "$3")
ARG4_TRIM=$(trim_whitespace "$4")
ABS_INPUT=$(abs_path "$ARG1_TRIM")
ABS_TRAJECTORY=$(abs_path "$ARG3_TRIM")
ABS_CELL_SELECT=$(abs_path "$ARG4_TRIM")
set -- "$ABS_INPUT" "$ARG2_TRIM" "$ABS_TRAJECTORY" "$ABS_CELL_SELECT" "${@:5}"
# Default HISTORY_LENGTH from 5th argument when not set via interactive mode
HISTORY_LENGTH="${HISTORY_LENGTH:-$5}"

# Normalise pair mode for downstream Python preprocess scripts
if [ -n "${PAIR_MODE:-}" ]; then
  if [ "$PAIR_MODE" = "all_pair" ]; then
    case "${MODE_CODE:-0}" in
      0) TENET_PAIR_MODE="gene_only" ;;      # RNA-only, all gene×gene
      1) TENET_PAIR_MODE="gene_only" ;;      # TENET_Plus, all gene×gene
      2) TENET_PAIR_MODE="all_feature" ;;    # TENET_Plus, all features×features
      *) TENET_PAIR_MODE="gene_only" ;;
    esac
  else
    TENET_PAIR_MODE="default"
  fi
elif [ -z "${TENET_PAIR_MODE:-}" ]; then
  TENET_PAIR_MODE="default"
fi

# Normalise intermediate save choice for logging/replay
TENET_INTERMEDIATE_SAVE="${INTERMEDIATE_SAVE_CHOICE:-${TENET_INTERMEDIATE_SAVE:-$DEFAULT_INTERMEDIATE_SAVE}}"

# Normalise TE batch size for logging/replay
if [ -n "${BATCH_SIZE:-}" ]; then
  TE_BATCH_SIZE="$BATCH_SIZE"
elif [ -n "${TE_BATCH_SIZE:-}" ]; then
  TE_BATCH_SIZE="$TE_BATCH_SIZE"
else
  TE_BATCH_SIZE="100"
fi

write_replay_logs "$@"

cd "$OUTPUT_DIR"

INPUT_FILE="$1"
ARG7="$7"
TEMP_CSV="temp_input.csv"

# Function to convert only the necessary part of Parquet to CSV (header)
convert_parquet_to_csv() {
    if ! command -v "$PYTHON" &> /dev/null; then
        echo "Python interpreter ($PYTHON) could not be found. Please install it to handle Parquet files."
        exit 1
    fi

    # Check if pandas is installed
    if ! "$PYTHON" -c "import pandas" 2>/dev/null; then
        echo "Pandas is required to handle Parquet files. Please install it in the current environment."
        exit 1
    fi

    echo "Extracting header from Parquet file to CSV..."
    "$PYTHON" - <<END
import pandas as pd
import sys

try:
    df = pd.read_parquet("$INPUT_FILE", columns=None)
    df.head(0).to_csv("$TEMP_CSV", index=False)
except Exception as e:
    print(f"Error converting Parquet to CSV: {e}")
    sys.exit(1)
END
}

# Determine file type based on extension
FILE_EXT="${INPUT_FILE##*.}"

if [[ "$FILE_EXT" == "parquet" ]]; then
    convert_parquet_to_csv
    DATA_FILE="$TEMP_CSV"
elif [[ -f "$INPUT_FILE" ]]; then
    DATA_FILE="$INPUT_FILE"
else
    echo "Input file does not exist or is not a regular file."
    exit 1
fi

# Validate ARG7
if [ -z "$ARG7" ]; then
    echo "No arguments were passed. Please check the command."
    [ "$FILE_EXT" == "parquet" ] && rm -f "$TEMP_CSV"
    exit 1
fi

if ! [[ "$ARG7" =~ ^[0-6]$ ]]; then
    echo "An invalid argument was received. Please check the command."
    [ "$FILE_EXT" == "parquet" ] && rm -f "$TEMP_CSV"
    exit 1
fi

# Process the first line of the data file
FIRST_LINE=$(head -n 1 "$DATA_FILE" | sed 's/\r//g' | sed 's/,/\n/g')

if [ "$ARG7" -eq 0 ]; then
    RESULT=$(echo "$FIRST_LINE" | grep -P "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" || true)
    if [ -n "$RESULT" ]; then
        echo "--- Chromosome information was found from the matrix in TENET_TF setting. Please check the input matrix or adjust setting($ARG7). ---"
        [ "$FILE_EXT" == "parquet" ] && rm -f "$TEMP_CSV"
        exit 1
    fi
else
    RESULT=$(echo "$FIRST_LINE" | grep -P "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" || true)
    if [ -z "$RESULT" ]; then
        echo "--- Chromosome information was not found from the matrix in TENET Plus setting. Please check the input matrix or adjust setting($ARG7). ---"
        [ "$FILE_EXT" == "parquet" ] && rm -f "$TEMP_CSV"
        exit 1
    fi
fi

if [ "$ARG7" -ne 0 ]; then
    echo 'Starting TENET_Plus'
    echo "--- Make TE_peak_list.txt from input Matrix--- "
    # If INPUT_FILE is Parquet, extract header columns; otherwise, use the original file
    if [[ "$FILE_EXT" == "parquet" ]]; then
        if ! "$PYTHON" - <<END > TE_peak_list.txt
import pandas as pd
import sys
import re
try:
    df = pd.read_parquet("$INPUT_FILE", columns=None)
    pattern = re.compile(r"^chr[0-9XY]+-([0-9]+)-([0-9]+)$")
    filtered_columns = [col for col in df.columns if pattern.match(str(col))]
    print("\n".join(filtered_columns))
except Exception as e:
    print(f"Error reading Parquet file: {e}")
    sys.exit(1)
END
        then
            echo "Error generating TE_peak_list.txt from Parquet file."
            exit 1
        fi
    else
        head -n 1 "$INPUT_FILE" | sed 's/\r//g' | sed 's/,/\n/g' | grep -P "^chr[0-9XY]+-([0-9]+)-([0-9]+)$" > TE_peak_list.txt
    fi
else
    echo "Starting TENET_TF"
fi

# Clean up temporary CSV file if it was created
if [ "$FILE_EXT" == "parquet" ]; then
    rm -f "$TEMP_CSV"
fi

# Check if cell_select_file contains only one cell or cells are sorted by pseudotime order
if grep -qv '^1$' "$4" || ! sort -gc "$3" || ! file "$1" | grep -q 'parquet'; then
    set_stage "MATRIX_PROCESSING"
    echo "--- Matrix processing, for converting matrix into parquet format, selecting cell_select_file with only 1 cell, and soring cells by pseudotime order ---"
    if ! "$PYTHON" -m code.process_matrix "$1" "$3" "$4"; then
        echo "Error in process_matrix.py. Exiting."
        exit 1
    fi
    # Move filtered artifacts into input directory and reference them there
    if [ -f filtered_matrix.parquet ]; then
      mv -f filtered_matrix.parquet "$INPUT_DIR/filtered_matrix.parquet"
    fi
    if [ -f filtered_trajectory.txt ]; then
      mv -f filtered_trajectory.txt "$INPUT_DIR/filtered_trajectory.txt"
    fi
    if [ -f filtered_cellselect.txt ]; then
      mv -f filtered_cellselect.txt "$INPUT_DIR/filtered_cellselect.txt"
    fi
    matrix="$INPUT_DIR/filtered_matrix.parquet"     # New matrix file
    trajectory="$INPUT_DIR/filtered_trajectory.txt" # New trajectory file
    cell_select="$INPUT_DIR/filtered_cellselect.txt"    # New cell select file
else
    set_stage "MATRIX_PROCESSING_SKIP"
    echo "--- Skipping Matrix processing, as cell_select_file contains only 1 cell and cells are sorted by pseudotime order."
    matrix=$1
    trajectory=$3
    cell_select=$4
fi


GENE_FILTER_MODE="${TENET_GENE_FILTER:-none}"
if [ "$GENE_FILTER_MODE" != "none" ]; then
    set_stage "GENE_FILTER"
    echo "--- Applying gene filter (${GENE_FILTER_MODE}) to matrix: ${matrix} ---"
    case "$GENE_FILTER_MODE" in
        ribo_mito|RIBO_MITO)
            GENE_FILTER_PREFIXES="RPS,RPL,MRPS,MRPL,MT-"
            ;;
        *)
            GENE_FILTER_PREFIXES=""
            ;;
    esac
    if [ -n "$GENE_FILTER_PREFIXES" ]; then
        FILTERED_MATRIX_PATH="$INPUT_DIR/filtered_matrix_gene_filtered.parquet"
        if ! "$PYTHON" -m code.filter_genes --input "${matrix}" --output "$FILTERED_MATRIX_PATH" --exclude_prefixes "$GENE_FILTER_PREFIXES"; then
            echo "Error in filter_genes.py. Exiting."
            exit 1
        fi
        matrix="$FILTERED_MATRIX_PATH"
    else
        echo "TENET_GENE_FILTER='${GENE_FILTER_MODE}' not recognised; skipping gene filtering."
    fi
fi


# Always regenerate gene_names from the current matrix
set_stage "GENE_NAMES_REGEN"
echo "--- Regenerating gene_names from matrix header ---"
if ! "$PYTHON" - <<PY
import sys
import pandas as pd
import pathlib
path = pathlib.Path("${matrix}")
ext = path.suffix.lower()
try:
    if ext == '.parquet':
        df = pd.read_parquet(path)
    elif ext == '.csv':
        df = pd.read_csv(path, index_col=0)
    elif ext in ('.tsv', '.txt'):
        df = pd.read_table(path, index_col=0)
    else:
        # fallback try csv
        df = pd.read_csv(path, index_col=0)
    cols = list(map(str, df.columns.tolist()))
    # Write gene_names into both input and output directories
    for dst in [r'${INPUT_DIR}/gene_names', r'${OUTPUT_DIR}/gene_names']:
        with open(dst, 'w') as f:
            for c in cols:
                f.write(f"{c}\n")
except Exception as e:
    print(f"Error extracting gene_names: {e}")
    sys.exit(1)
PY
then
  echo "Failed to write gene_names from matrix header"
  exit 1
fi
# Modality-specific preprocessing (produces cell_gene_trsps.parquet)
MODALITY_DEFAULT="auto"
if [ "$ARG7" -eq 0 ]; then MODALITY_DEFAULT="rna"; fi
MODALITY="${8:-$MODALITY_DEFAULT}"

# If a precomputed transposed matrix exists, remove it to regenerate fresh
if [ -f cell_gene_trsps.parquet ]; then
  echo "--- Removing existing cell_gene_trsps.parquet to regenerate ---"
  rm -f cell_gene_trsps.parquet
fi

# Run modality preprocessing unless explicitly skipped; otherwise fallback transpose runs later
if [[ "$MODALITY" == "none" || "$MODALITY" == "skip" ]]; then
  set_stage "MODALITY_SKIP"
  echo "--- Skipping modality preprocessing as requested (${MODALITY}) ---"
else
  echo "--- Modality preprocessing (${MODALITY}) -> cell_gene_trsps.parquet ---"
  if ! run_with_stage "MODALITY_PREPROCESS(${MODALITY})" "$PYTHON" -m code.modality_preprocess --input "${matrix}" --output cell_gene_trsps.parquet --modality "${MODALITY}"; then
    echo "Error in modality_preprocess.py. Exiting."
    exit 1
  fi
fi

if [ "$7" -ne 0 ]; then
  echo 'Preprocessing(TENET_Plus)'
    if [ "${TENET_PAIR_MODE:-default}" = "default" ]; then
      if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "PREPROCESS_TENET_PLUS" "$PYTHON" -m code.PreProcessScript_TE_Plus $6 $7 "${matrix}"; then
          echo "Error in PreProcessScript_TE_Plus.py. Exiting."
          exit 1
      fi
    else
      echo "--- Skipping PreProcessScript_TE_Plus (pair_mode=${TENET_PAIR_MODE}); using implicit all-pairs in TE core ---"
    fi
else
  echo 'Preprocessing(TENET_TF)'
  if [ "${TENET_PAIR_MODE:-default}" = "default" ]; then
    if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "PREPROCESS_TENET_TF" "$PYTHON" -m code.PreProcessScript_TENET_TF $6; then
        echo "Error in PreProcessScript_TENET_TF.py. Exiting."
        exit 1
    fi
  else
    echo "--- Skipping PreProcessScript_TENET_TF (pair_mode=${TENET_PAIR_MODE}); using implicit all-pairs in TE core ---"
  fi
fi

# Ensure transposed matrix exists (for TF mode or when modality preprocessing was skipped and precomputed not found)
if [ ! -f cell_gene_trsps.parquet ]; then
  echo "--- Generating cell_gene_trsps.parquet by transposing matrix ---"
  if ! "$PYTHON" - <<PY
import sys, pathlib
import pandas as pd
path = pathlib.Path("${matrix}")
ext = path.suffix.lower()
try:
    if ext == '.parquet':
        df = pd.read_parquet(path)
    elif ext == '.csv':
        df = pd.read_csv(path, index_col=0)
    elif ext in ('.tsv', '.txt'):
        df = pd.read_table(path, index_col=0)
    else:
        df = pd.read_csv(path, index_col=0)
    df.T.to_parquet('cell_gene_trsps.parquet')
except Exception as e:
    print(f"Error generating transpose: {e}")
    sys.exit(1)
PY
  then
    echo "Error creating cell_gene_trsps.parquet by transpose. Exiting."
    exit 1
  fi
fi

# If requested, automatically select history length k based on series length
if [[ "${HISTORY_LENGTH,,}" == "auto" ]]; then
  set_stage "HISTORY_AUTO_SELECT"
  echo "--- Auto-selecting history length k from cell_gene_trsps.parquet (per-series AR-BIC) ---"
  AUTO_K=$("$PYTHON" - <<'PY'
import math

import numpy as np
import pandas as pd


def per_series_best_k(series_matrix, max_k, max_series=128):
    """
    For each sampled series, pick the AR order k (1..max_k) that minimises BIC,
    then aggregate these per-series orders into a global k via the median.
    """
    n_genes, n_time = series_matrix.shape
    if n_time < 4:
        return 1

    use_genes = min(n_genes, max_series)
    idx = np.linspace(0, n_genes - 1, use_genes, dtype=int)
    sub = series_matrix[idx]

    best_orders = []

    for row in sub:
        y = np.asarray(row, dtype=float)
        if not np.all(np.isfinite(y)):
            continue
        mu = y.mean()
        sig = y.std()
        if not np.isfinite(sig) or sig <= 0:
            continue
        y = (y - mu) / sig
        N = y.size
        bics = []
        ks = []
        for k in range(1, max_k + 1):
            if N <= k + 1:
                break
            n_obs = N - k
            Y = y[k:]
            X = np.empty((n_obs, k), dtype=float)
            for i in range(k):
                X[:, i] = y[k - i - 1 : N - i - 1]
            try:
                beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
            except np.linalg.LinAlgError:
                continue
            resid = Y - X @ beta
            rss = float(resid @ resid)
            if not np.isfinite(rss) or rss <= 0:
                continue
            bic = n_obs * math.log(rss / n_obs) + k * math.log(n_obs)
            if math.isfinite(bic):
                bics.append(bic)
                ks.append(k)
        if bics:
            # series-specific best order
            k_best = ks[int(np.argmin(bics))]
            best_orders.append(int(k_best))

    if not best_orders:
        return 1

    # Use the median order across series for robustness
    orders = np.asarray(best_orders, dtype=int)
    k_med = int(np.median(orders))
    k_med = max(1, min(int(max_k), k_med))
    return k_med


try:
    df = pd.read_parquet("cell_gene_trsps.parquet")
    values = df.to_numpy()
    n_timepoints = int(values.shape[1])
except Exception:
    print(1)
else:
    # Upper bound for k: keep at least ~min_effective observations and cap at 10
    min_effective = 50
    max_k = min(10, max(1, n_timepoints - min_effective))
    if max_k <= 1:
        print(1)
    else:
        k_sel = per_series_best_k(values, max_k=max_k)
        if not isinstance(k_sel, int) or k_sel <= 0:
            k_sel = 1
        print(int(k_sel))
PY
  ) || AUTO_K=1
  if ! [[ "$AUTO_K" =~ ^[0-9]+$ ]] || [ "$AUTO_K" -le 0 ]; then
    AUTO_K=1
  fi
  HISTORY_LENGTH="$AUTO_K"
  echo "--- Auto-selected history length k=$HISTORY_LENGTH ---"
fi

# TE computation: screen and refine options
SCREEN_MODE="${9:-kernel}"
REFINE_METHOD="${10:-none}"
REFINE_TOPK="${11:-0}"
REFINE_TOP_PCT="${12:-0}"
echo "--- TE stage: screen=${SCREEN_MODE} refine=${REFINE_METHOD} topk=${REFINE_TOPK} pct=${REFINE_TOP_PCT} ---"

# Permutation options (optional)
PERMUTE="${13:-off}"
PERM_N="${14:-100}"
PERM_TOPK="${15:-0}"
PERM_TOP_PCT="${16:-0}"
PERM_FDR="${17:-off}"
PERM_Q_ALPHA="${18:-0.05}"
PERM_ALPHA="${19:-0.01}"
STORE_LOCAL_TE="${20:-off}"
LOCAL_ARGS=()
if [[ "$STORE_LOCAL_TE" == "on" || "$STORE_LOCAL_TE" == "ON" || "$STORE_LOCAL_TE" == "1" ]]; then
  LOCAL_ARGS=(--store_local_te)
  echo "--- Local TE storage enabled ---"
fi

# Time subsampling options (optional)
TIME_STRIDE="${21:-1}"
TIME_PCT="${22:-100}"
TIME_SEED="${23:-42}"
TIME_ARGS=()
if [[ "$TIME_STRIDE" != "1" ]]; then
  TIME_ARGS+=(--time_stride "$TIME_STRIDE")
else
  if [[ "$TIME_PCT" != "100" ]]; then
    TIME_ARGS+=(--time_pct "$TIME_PCT" --time_seed "$TIME_SEED")
  fi
fi
if [[ ${#TIME_ARGS[@]} -gt 0 ]]; then
  echo "--- Time subsampling: ${TIME_ARGS[@]} ---"
fi

# TE results buffer rows (optional; controls when intermediate batches flush to disk)
RESULTS_BUFFER_ROWS="${24:-}"
BUFFER_ARGS=()
if [[ -n "${RESULTS_BUFFER_ROWS:-}" ]]; then
  if [[ "$RESULTS_BUFFER_ROWS" =~ ^[0-9]+$ ]] && [ "$RESULTS_BUFFER_ROWS" -gt 0 ]; then
    BUFFER_ARGS+=(--results_buffer_rows "$RESULTS_BUFFER_ROWS")
    echo "--- TE buffering: results_buffer_rows=$RESULTS_BUFFER_ROWS ---"
  else
    echo "Warning: Ignoring invalid RESULTS_BUFFER_ROWS='$RESULTS_BUFFER_ROWS' (expected positive integer)."
  fi
fi

# Decide intermediate saving behaviour
INTERMEDIATE_SAVE_CHOICE="${INTERMEDIATE_SAVE_CHOICE:-${TENET_INTERMEDIATE_SAVE:-$DEFAULT_INTERMEDIATE_SAVE}}"
ENABLE_INTERMEDIATE_ARGS=()
case "${INTERMEDIATE_SAVE_CHOICE,,}" in
  on|1|yes)
    ENABLE_INTERMEDIATE_ARGS=(--enable_intermediate_save)
    echo "--- Intermediate save: ON ---"
    ;;
  off|0|no)
    echo "--- Intermediate save: OFF ---"
    ;;
  *)
    # Fallback to default
    if [[ "${DEFAULT_INTERMEDIATE_SAVE,,}" == "on" || "${DEFAULT_INTERMEDIATE_SAVE,,}" == "yes" ]]; then
      ENABLE_INTERMEDIATE_ARGS=(--enable_intermediate_save)
      echo "--- Intermediate save: ON (default) ---"
    else
      echo "--- Intermediate save: OFF (default) ---"
    fi
    ;;
esac

# TE batch size (optional; controls sources per target chunk)
TE_BATCH_SIZE_VALUE="${BATCH_SIZE:-${TE_BATCH_SIZE:-}}"
TE_BATCH_ARGS=()
if [[ -n "${TE_BATCH_SIZE_VALUE:-}" ]]; then
  if [[ "$TE_BATCH_SIZE_VALUE" =~ ^[0-9]+$ ]] && [ "$TE_BATCH_SIZE_VALUE" -gt 0 ]; then
    TE_BATCH_ARGS+=(--batch_size "$TE_BATCH_SIZE_VALUE")
    echo "--- TE batch size: $TE_BATCH_SIZE_VALUE ---"
  else
    echo "Warning: Ignoring invalid TE batch size '$TE_BATCH_SIZE_VALUE' (expected positive integer)."
  fi
fi

if [[ "$SCREEN_MODE" == "linear" || "$SCREEN_MODE" == "poly" || "$SCREEN_MODE" == "gcmi" || "$SCREEN_MODE" == "disc" || "$SCREEN_MODE" == "ordinal" || "$SCREEN_MODE" == "kernel_grid" ]]; then
  # Screening with fast estimator, optional refinement
  if [[ "$PERMUTE" == "on" || "$PERMUTE" == "ON" || "$PERMUTE" == "1" ]]; then
    echo "--- Permutation enabled for ${SCREEN_MODE}: n=${PERM_N} topk=${PERM_TOPK} pct=${PERM_TOP_PCT} fdr=${PERM_FDR} ---"
    PERM_ARGS=(--permute --perm_n "$PERM_N")
    if [[ "$PERM_TOPK" != "0" ]]; then PERM_ARGS+=(--permute_topk_per_target "$PERM_TOPK"); fi
    if [[ "$PERM_TOP_PCT" != "0" ]]; then PERM_ARGS+=(--permute_top_pct "$PERM_TOP_PCT"); fi
    if [[ "$PERM_FDR" == "on" || "$PERM_FDR" == "ON" || "$PERM_FDR" == "1" ]]; then
      PERM_ARGS+=(--use_fdr --perm_q_alpha "$PERM_Q_ALPHA")
    else
      PERM_ARGS+=(--perm_alpha "$PERM_ALPHA")
    fi
    if ! run_with_stage "TE_SCREEN(${SCREEN_MODE})_PERMUTE" "$PYTHON" -m code.runTE_for_py_python_batch all_pairs.csv "$2" "$HISTORY_LENGTH" \
        "${ENABLE_INTERMEDIATE_ARGS[@]}" \
        --mode "$SCREEN_MODE" \
        "${TE_BATCH_ARGS[@]}" \
        --pair_mode "$TENET_PAIR_MODE" \
        "${PERM_ARGS[@]}" \
        "${BUFFER_ARGS[@]}" \
        "${LOCAL_ARGS[@]}" \
        "${TIME_ARGS[@]}"; then
      echo "Error in runTE_for_py_python_batch.py (poly/linear permute). Exiting."
      exit 1
    fi
  else
    # Build optional refinement args only when valid
    REFINE_ARGS=()
    if [[ "$REFINE_METHOD" == "kernel" || "$REFINE_METHOD" == "ksg" ]]; then
      if [[ "$REFINE_TOPK" != "0" || "$REFINE_TOP_PCT" != "0" ]]; then
        REFINE_ARGS=(
          --hybrid_refine_method "$REFINE_METHOD" \
          --hybrid_refine_topk_per_target "$REFINE_TOPK" \
          --hybrid_refine_top_pct "$REFINE_TOP_PCT" \
          --hybrid_output TE_result_all.parquet
        )
      fi
    elif [[ "$REFINE_METHOD" == "none" ]]; then
      echo "--- Refinement disabled (refine=none); running screen-only. ---"
    else
      echo "Warning: Unknown refine method '$REFINE_METHOD'; skipping refinement."
    fi

    if ! run_with_stage "TE_SCREEN(${SCREEN_MODE})_REFINE_${REFINE_METHOD}" "$PYTHON" -m code.runTE_for_py_python_batch all_pairs.csv "$2" "$HISTORY_LENGTH" \
        "${ENABLE_INTERMEDIATE_ARGS[@]}" \
        --mode "$SCREEN_MODE" \
        "${TE_BATCH_ARGS[@]}" \
        --pair_mode "$TENET_PAIR_MODE" \
        "${REFINE_ARGS[@]}" \
        "${BUFFER_ARGS[@]}" \
        "${LOCAL_ARGS[@]}" \
        "${TIME_ARGS[@]}"; then
        echo "Error in runTE_for_py_python_batch.py (screen/refine). Exiting."
        exit 1
    fi
  fi
else
  # Direct method (ksg/kernel)
  if [[ "$PERMUTE" == "on" || "$PERMUTE" == "ON" || "$PERMUTE" == "1" ]]; then
    echo "--- Permutation enabled for ${SCREEN_MODE}: n=${PERM_N} topk=${PERM_TOPK} pct=${PERM_TOP_PCT} fdr=${PERM_FDR} ---"
    PERM_ARGS=(--permute --perm_n "$PERM_N")
    if [[ "$PERM_TOPK" != "0" ]]; then PERM_ARGS+=(--permute_topk_per_target "$PERM_TOPK"); fi
    if [[ "$PERM_TOP_PCT" != "0" ]]; then PERM_ARGS+=(--permute_top_pct "$PERM_TOP_PCT"); fi
    if [[ "$PERM_FDR" == "on" || "$PERM_FDR" == "ON" || "$PERM_FDR" == "1" ]]; then
      PERM_ARGS+=(--use_fdr --perm_q_alpha "$PERM_Q_ALPHA")
    else
      PERM_ARGS+=(--perm_alpha "$PERM_ALPHA")
    fi
    if ! run_with_stage "TE_DIRECT(${SCREEN_MODE})_PERMUTE" "$PYTHON" -m code.runTE_for_py_python_batch all_pairs.csv "$2" "$HISTORY_LENGTH" \
        "${ENABLE_INTERMEDIATE_ARGS[@]}" \
        --mode "$SCREEN_MODE" \
        "${TE_BATCH_ARGS[@]}" \
        --pair_mode "$TENET_PAIR_MODE" \
        "${PERM_ARGS[@]}" \
        "${BUFFER_ARGS[@]}" \
        "${LOCAL_ARGS[@]}" \
        "${TIME_ARGS[@]}"; then
      echo "Error in runTE_for_py_python_batch.py (direct permute). Exiting."
      exit 1
    fi
  else
    if ! run_with_stage "TE_DIRECT(${SCREEN_MODE})" "$PYTHON" -m code.runTE_for_py_python_batch all_pairs.csv "$2" "$HISTORY_LENGTH" \
        "${ENABLE_INTERMEDIATE_ARGS[@]}" \
        --mode "$SCREEN_MODE" \
        "${TE_BATCH_ARGS[@]}" \
        --pair_mode "$TENET_PAIR_MODE" \
        "${BUFFER_ARGS[@]}" \
        "${LOCAL_ARGS[@]}" \
        "${TIME_ARGS[@]}"; then
      echo "Error in runTE_for_py_python_batch.py (direct). Exiting."
      exit 1
    fi
  fi
fi

if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "MATRIX_GENERATE" "$PYTHON" -m code.Matrix_generate TE_result_all.parquet $6 $7; then
    echo "Error in Matrix_generate.py. Exiting."
    exit 1
fi

if [[ "$STORE_LOCAL_TE" == "on" || "$STORE_LOCAL_TE" == "ON" || "$STORE_LOCAL_TE" == "1" ]]; then
    echo "--- Chunking LocalTE payloads (${LOCAL_TE_CHUNK_SIZE} timepoints per chunk) ---"
    set_stage "LOCALTE_EXPORT"
    if [[ "${LOCAL_TE_SPLIT_EXPORT,,}" == "on" ]]; then
        # Build selector-based inputs if selector files exist
        SEL_INPUTS=()
        if [[ -f TE_TF_GN.parquet ]]; then
            SEL_INPUTS+=("TE_result_all.parquet=TE_TF_GN.parquet")
        fi
        if [[ -f TE_TF_PK.parquet ]]; then
            SEL_INPUTS+=("TE_result_all.parquet=TE_TF_PK.parquet")
        fi
        if [[ -f TE_PK_GN.parquet ]]; then
            SEL_INPUTS+=("TE_result_all.parquet=TE_PK_GN.parquet")
        fi
        # If no selectors found, fall back to single export
        if [[ ${#SEL_INPUTS[@]} -gt 0 ]]; then
            echo "[localte_chunk_export] Split export using selectors -> $LOCAL_TE_SPLIT_OUTPUT_DIR"
            set_stage "LOCALTE_EXPORT_SPLIT"
            if ! run_with_stage "LOCALTE_EXPORT_SPLIT" "$PYTHON" -m code.localte_chunk_export \
                --input "${SEL_INPUTS[@]}" \
                --output_dir "$LOCAL_TE_SPLIT_OUTPUT_DIR" \
                --chunk_size "$LOCAL_TE_CHUNK_SIZE" \
                --buffer_edges "$LOCAL_TE_BUFFER_EDGES" \
                --read_batch_rows "$LOCAL_TE_READ_BATCH_ROWS" \
                --use_threads "$LOCAL_TE_USE_THREADS" \
                --values_dtype "$LOCAL_TE_VALUES_DTYPE" \
                --workers "$LOCAL_TE_EXPORT_WORKERS" \
                ${LOCAL_TE_MERGE_WORKERS:+--merge_workers "$LOCAL_TE_MERGE_WORKERS"}; then
                echo "Error in localte_chunk_export.py (split mode). Exiting."
                exit 1
            fi
        else
            echo "[localte_chunk_export] No selector files found; exporting all edges to local_te_chunks"
            if ! run_with_stage "LOCALTE_EXPORT_ALL" "$PYTHON" -m code.localte_chunk_export \
                --input TE_result_all.parquet \
                --output_dir local_te_chunks \
                --chunk_size "$LOCAL_TE_CHUNK_SIZE" \
                --buffer_edges "$LOCAL_TE_BUFFER_EDGES" \
                --read_batch_rows "$LOCAL_TE_READ_BATCH_ROWS" \
                --use_threads "$LOCAL_TE_USE_THREADS" \
                --values_dtype "$LOCAL_TE_VALUES_DTYPE" \
                --workers "$LOCAL_TE_EXPORT_WORKERS" \
                ${LOCAL_TE_MERGE_WORKERS:+--merge_workers "$LOCAL_TE_MERGE_WORKERS"} \
                --scores_output TE_result_scores.parquet; then
                echo "Error in localte_chunk_export.py. Exiting."
                exit 1
            fi
        fi
    else
        if ! run_with_stage "LOCALTE_EXPORT_ALL" "$PYTHON" -m code.localte_chunk_export \
            --input TE_result_all.parquet \
            --output_dir local_te_chunks \
            --chunk_size "$LOCAL_TE_CHUNK_SIZE" \
            --buffer_edges "$LOCAL_TE_BUFFER_EDGES" \
            --read_batch_rows "$LOCAL_TE_READ_BATCH_ROWS" \
            --use_threads "$LOCAL_TE_USE_THREADS" \
            --values_dtype "$LOCAL_TE_VALUES_DTYPE" \
            --workers "$LOCAL_TE_EXPORT_WORKERS" \
            ${LOCAL_TE_MERGE_WORKERS:+--merge_workers "$LOCAL_TE_MERGE_WORKERS"} \
            --scores_output TE_result_scores.parquet; then
            echo "Error in localte_chunk_export.py. Exiting."
            exit 1
        fi
    fi
fi
