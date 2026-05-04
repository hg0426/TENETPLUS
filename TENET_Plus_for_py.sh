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
LOCAL_TE_BUFFER_EDGES=${LOCAL_TE_BUFFER_EDGES:-10000}
# Export tuning for localTE chunking
LOCAL_TE_EXPORT_WORKERS=${LOCAL_TE_EXPORT_WORKERS:-auto}      # auto = min(jobs,16); 0/1 = single process; >1 = multiprocessing
LOCAL_TE_READ_BATCH_ROWS=${LOCAL_TE_READ_BATCH_ROWS:-8192}
LOCAL_TE_USE_THREADS=${LOCAL_TE_USE_THREADS:-on}
LOCAL_TE_VALUES_DTYPE=${LOCAL_TE_VALUES_DTYPE:-float16}
LOCAL_TE_MERGE_WORKERS=${LOCAL_TE_MERGE_WORKERS:-}
LOCAL_TE_MERGE_PARTS=${LOCAL_TE_MERGE_PARTS:-off}
LOCAL_TE_CODEC=${LOCAL_TE_CODEC:-${TE_LOCALTE_CODEC:-zlib}}
LOCAL_TE_ADVANCED=${LOCAL_TE_ADVANCED:-${TENET_LOCALTE_ADVANCED:-off}}
# Optional: export per-selector chunk dirs (payload filtered by selector files)
LOCAL_TE_SPLIT_EXPORT=${LOCAL_TE_SPLIT_EXPORT:-on}
LOCAL_TE_SPLIT_OUTPUT_DIR=${LOCAL_TE_SPLIT_OUTPUT_DIR:-local_te_split_chunks}

TENET_GENE_FILTER=${TENET_GENE_FILTER:-none}

# Global default for pair mode (overridable via interactive prompt or TENET_PAIR_MODE env)
DEFAULT_PAIR_MODE=${TENET_PAIR_MODE:-default}
DEFAULT_INTERMEDIATE_SAVE=${TENET_INTERMEDIATE_SAVE:-off}
DEFAULT_MAKE_GRN=${TENET_MAKE_GRN:-on}
DEFAULT_GRN_FDR=${TENET_GRN_FDR:-0.01}
DEFAULT_GRN_TE_CUTOFF=${TENET_GRN_TE_CUTOFF:-0.0}
DEFAULT_TRIM_INDIRECT=${TENET_TRIM_INDIRECT:-on}
DEFAULT_TRIM_THRESHOLD=${TENET_TRIM_THRESHOLD:--0.01}
DEFAULT_GRN_THREADS=${TENET_GRN_THREADS:-}

# Kernel mode defaults. These make the exact grouped-state implementation the
# package default while still allowing advanced users to opt out via env vars.
export TE_KERNEL_GROUPED_STATE=${TE_KERNEL_GROUPED_STATE:-on}
export TE_KERNEL_GROUPED_STATE_CODED=${TE_KERNEL_GROUPED_STATE_CODED:-on}
export TE_KERNEL_GROUPED_STATE_FACTORIZED=${TE_KERNEL_GROUPED_STATE_FACTORIZED:-auto}
export TE_KERNEL_GROUPED_STATE_SOURCE_PREFIX=${TE_KERNEL_GROUPED_STATE_SOURCE_PREFIX:-auto}
export TE_KERNEL_SOURCE_CODE_INTEGER=${TE_KERNEL_SOURCE_CODE_INTEGER:-auto}
export TE_KERNEL_GROUPED_STATE_SPARSE=${TE_KERNEL_GROUPED_STATE_SPARSE:-off}
export TE_KERNEL_GROUPED_STATE_2D_RANGE=${TE_KERNEL_GROUPED_STATE_2D_RANGE:-on}
export TE_KERNEL_PERM_TABLE_SAMPLER=${TE_KERNEL_PERM_TABLE_SAMPLER:-off}
export TE_PERM_KERNEL_GROUPED_SOURCES_PER_CHUNK=${TE_PERM_KERNEL_GROUPED_SOURCES_PER_CHUNK:-2}
export TE_PERM_SORT_WORK=${TE_PERM_SORT_WORK:-on}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
CODE_DIR="$REPO_ROOT/code"
INPUT_DIR="$REPO_ROOT/input"
OUTPUT_DIR="${TENET_OUTPUT_DIR:-$REPO_ROOT/output}"
if [[ "$OUTPUT_DIR" != /* ]]; then
    OUTPUT_DIR="$REPO_ROOT/$OUTPUT_DIR"
fi
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


set_analysis_mode_from_choice() {
    local choice="${1,,}"
    case "$choice" in
        rna_tf|tf_rna|tf|rna|tenet_tf|0)
            PAIR_MODE="default"; MODE_CODE="0"; ANALYSIS_MODE_CHOICE="rna_tf" ;;
        plus_full|full|tenet_plus_full|1)
            PAIR_MODE="default"; MODE_CODE="1"; ANALYSIS_MODE_CHOICE="plus_full" ;;
        tf_gene|tf_gn|tf-gene|gene_grn|grn|2)
            PAIR_MODE="default"; MODE_CODE="2"; ANALYSIS_MODE_CHOICE="tf_gene" ;;
        tf_peak|tf_pk|tf-peak|3)
            PAIR_MODE="default"; MODE_CODE="3"; ANALYSIS_MODE_CHOICE="tf_peak" ;;
        tf_targets|tf_gene_peak|tf_gn_pk|tf-gene-peak|tf_gene+peak|4)
            PAIR_MODE="default"; MODE_CODE="4"; ANALYSIS_MODE_CHOICE="tf_targets" ;;
        peak_gene|peak_gn|peak-gene|cis_peak_gene|5)
            PAIR_MODE="default"; MODE_CODE="5"; ANALYSIS_MODE_CHOICE="peak_gene" ;;
        peak_peak|peak_pk|peak-peak|cis_peak_peak|6)
            PAIR_MODE="default"; MODE_CODE="6"; ANALYSIS_MODE_CHOICE="peak_peak" ;;
        all_gene_rna|rna_all_gene|a0|7)
            PAIR_MODE="all_pair"; MODE_CODE="0"; ANALYSIS_MODE_CHOICE="7" ;;
        all_gene|all_gene_plus|gene_all_gene|a1|8)
            PAIR_MODE="all_pair"; MODE_CODE="1"; ANALYSIS_MODE_CHOICE="8" ;;
        all_feature|all_features|feature_all_feature|a2|9)
            PAIR_MODE="all_pair"; MODE_CODE="2"; ANALYSIS_MODE_CHOICE="9" ;;
        *)
            return 1 ;;
    esac
}

analysis_mode_default_choice() {
    if [ -n "${ANALYSIS_MODE_CHOICE:-}" ]; then
        echo "$ANALYSIS_MODE_CHOICE"
        return
    fi
    local pair="${PAIR_MODE:-$DEFAULT_PAIR_MODE}"
    local mode="${MODE_CODE:-1}"
    if [ "$pair" = "all_pair" ]; then
        case "$mode" in
            0) echo "7" ;;
            1) echo "8" ;;
            2) echo "9" ;;
            *) echo "8" ;;
        esac
    else
        case "$mode" in
            0) echo "0" ;;
            1) echo "1" ;;
            2) echo "2" ;;
            3) echo "3" ;;
            4) echo "4" ;;
            5) echo "5" ;;
            6) echo "6" ;;
            *) echo "1" ;;
        esac
    fi
}

analysis_mode_label() {
    local pair="${PAIR_MODE:-$DEFAULT_PAIR_MODE}"
    local mode="${MODE_CODE:-1}"
    if [ "$pair" = "all_pair" ]; then
        case "$mode" in
            0) echo "7: RNA-only every gene -> every gene" ;;
            1) echo "8: every gene -> every gene" ;;
            2) echo "9: every gene/peak -> every gene/peak" ;;
            *) echo "all_pair mode=$mode" ;;
        esac
    else
        case "$mode" in
            0) echo "0: RNA-only TF->gene" ;;
            1) echo "1: TENET_Plus full (TF->gene, TF->peak, peak->gene)" ;;
            2) echo "2: TF->gene" ;;
            3) echo "3: TF->peak" ;;
            4) echo "4: TF->gene + TF->peak" ;;
            5) echo "5: peak->gene (cis)" ;;
            6) echo "6: peak->peak (cis)" ;;
            *) echo "default mode=$mode" ;;
        esac
    fi
}

prompt_analysis_mode() {
    local default
    default="$(analysis_mode_default_choice)"
    while true; do
        print_step "7" "Select Output Type" \
            "Enter a number. Press Enter for [1]." \
            "Default: [1] TENET_Plus full" \
            "All modes:" \
            "  [0] RNA-only TENET_TF" \
            "  [1] TENET_Plus full: TF->gene, TF->peak, peak->gene" \
            "  [2] TF->gene" \
            "  [3] TF->peak" \
            "  [4] TF->gene + TF->peak" \
            "  [5] peak->gene (cis)" \
            "  [6] peak->peak (cis)" \
            "  [7] RNA-only every gene -> every gene" \
            "  [8] every gene -> every gene" \
            "  [9] every gene/peak -> every gene/peak"
        local value
        value=$(prompt_with_default "Network mode number" "$default")
        if set_analysis_mode_from_choice "$value"; then
            TENET_PAIR_MODE="$PAIR_MODE"
            break
        fi
        echo "[!] Please enter a number: 0,1,2,3,4,5,6,7,8,9."
        echo "    Backward-compatible aliases also work: a0/a1/a2 or names such as tf_gene/all_feature."
        default="$value"
    done
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
    DEFAULT_PERM_SCOPE="grn_fdr"
    DEFAULT_PERM_GRN_FDR="${TENET_PERM_CANDIDATE_GRN_FDR:-0.01}"
    DEFAULT_PERM_GRN_TE_CUTOFF="${TENET_PERM_CANDIDATE_TE_CUTOFF:-0.0}"
    DEFAULT_LOCAL_TE="off"
    DEFAULT_LOCAL_EXPORT_WORKERS="${LOCAL_TE_EXPORT_WORKERS:-auto}"
    DEFAULT_LOCAL_CHUNK_SIZE="$LOCAL_TE_CHUNK_SIZE"
    DEFAULT_LOCAL_BUFFER="$LOCAL_TE_BUFFER_EDGES"
    DEFAULT_LOCAL_READ_BATCH="$LOCAL_TE_READ_BATCH_ROWS"
    DEFAULT_LOCAL_THREADS="$LOCAL_TE_USE_THREADS"
    DEFAULT_LOCAL_DTYPE="$LOCAL_TE_VALUES_DTYPE"
    DEFAULT_LOCAL_MERGE_WORKERS="${LOCAL_TE_MERGE_WORKERS:-auto}"
    DEFAULT_LOCAL_SPLIT_EXPORT="${LOCAL_TE_SPLIT_EXPORT:-on}"
    DEFAULT_LOCAL_SPLIT_DIR="${LOCAL_TE_SPLIT_OUTPUT_DIR:-local_te_split_chunks}"
    DEFAULT_RESULTS_BUFFER_ROWS="${TE_RESULTS_BUFFER_ROWS:-}"
    DEFAULT_INTERMEDIATE_SAVE="${TENET_INTERMEDIATE_SAVE:-off}"
    DEFAULT_BATCH_SIZE="${TE_BATCH_SIZE:-100}"
    DEFAULT_PAIR_MODE="${TENET_PAIR_MODE:-default}"
    DEFAULT_TIME_STRIDE="1"
    DEFAULT_TIME_PCT="100"
    DEFAULT_TIME_SEED="42"
    DEFAULT_MAKE_GRN="${TENET_MAKE_GRN:-on}"
    DEFAULT_GRN_FDR="${TENET_GRN_FDR:-0.01}"
    DEFAULT_GRN_TE_CUTOFF="${TENET_GRN_TE_CUTOFF:-0.0}"
    DEFAULT_TRIM_INDIRECT="${TENET_TRIM_INDIRECT:-on}"
    DEFAULT_TRIM_THRESHOLD="${TENET_TRIM_THRESHOLD:--0.01}"
    DEFAULT_GRN_THREADS="${TENET_GRN_THREADS:-$NUM_JOBS}"
}

set_permutation_off_defaults() {
    PERMUTE_CHOICE="off"
    PERM_N="$DEFAULT_PERM_N"
    PERM_TOPK="$DEFAULT_PERM_TOPK"
    PERM_TOP_PCT="$DEFAULT_PERM_TOP_PCT"
    PERM_FDR="$DEFAULT_PERM_FDR"
    PERM_Q_ALPHA="$DEFAULT_PERM_Q_ALPHA"
    PERM_ALPHA="$DEFAULT_PERM_ALPHA"
    PERM_SCOPE="$DEFAULT_PERM_SCOPE"
    PERM_GRN_FDR="$DEFAULT_PERM_GRN_FDR"
    PERM_GRN_TE_CUTOFF="$DEFAULT_PERM_GRN_TE_CUTOFF"
    unset TENET_PERM_CANDIDATE_GRN_FDR
    unset TENET_PERM_CANDIDATE_TE_CUTOFF
}

prompt_permutation_detail() {
    local default_n="${PERM_N:-$DEFAULT_PERM_N}"
    local default_grn_fdr="${PERM_GRN_FDR:-$DEFAULT_PERM_GRN_FDR}"
    if [ "${PERMUTE_CHOICE,,}" != "on" ]; then
        echo "[i] Permutation is off; no detail to edit."
        return
    fi
    print_step "8" "Permutation Settings" \
        "Only the necessary settings are shown." \
        "Advanced p/q thresholds keep defaults: raw p <= ${DEFAULT_PERM_ALPHA}, q <= ${DEFAULT_PERM_Q_ALPHA} if FDR is enabled via CLI/env."
    PERM_N=$(prompt_with_default "Permutation count" "$default_n")
    PERM_TOPK="0"
    PERM_TOP_PCT="0"
    PERM_FDR="$DEFAULT_PERM_FDR"
    PERM_Q_ALPHA="$DEFAULT_PERM_Q_ALPHA"
    PERM_ALPHA="$DEFAULT_PERM_ALPHA"
    case "${PERM_SCOPE:-$DEFAULT_PERM_SCOPE}" in
        grn_fdr)
            PERM_GRN_FDR=$(prompt_with_default "Candidate GRN FDR alpha" "$default_grn_fdr")
            PERM_GRN_TE_CUTOFF="$DEFAULT_PERM_GRN_TE_CUTOFF"
            export TENET_PERM_CANDIDATE_GRN_FDR="$PERM_GRN_FDR"
            export TENET_PERM_CANDIDATE_TE_CUTOFF="$PERM_GRN_TE_CUTOFF"
            ;;
        all)
            unset TENET_PERM_CANDIDATE_GRN_FDR
            unset TENET_PERM_CANDIDATE_TE_CUTOFF
            ;;
    esac
}

permutation_mode_label() {
    if [ "${PERMUTE_CHOICE:-off}" != "on" ]; then
        echo "off"
        return
    fi
    case "${PERM_SCOPE:-$DEFAULT_PERM_SCOPE}" in
        grn_fdr) echo "grn_edges" ;;
        all) echo "all_pairs" ;;
        *) echo "${PERM_SCOPE:-on}" ;;
    esac
}

prompt_permutation_toggle() {
    local default_mode
    if [ "${PERMUTE_CHOICE:-$DEFAULT_PERMUTE}" = "on" ]; then
        case "${PERM_SCOPE:-$DEFAULT_PERM_SCOPE}" in
            all) default_mode="all_pairs" ;;
            *) default_mode="grn_edges" ;;
        esac
    else
        default_mode="off"
    fi
    print_step "8" "Permutation Mode" \
        "off       = skip permutation testing" \
        "grn_edges = after kernel TE, first make FDR-filtered GRN candidates, then permute those edges" \
        "all_pairs = permute every selected TE pair (much slower)"
    while true; do
        local mode
        mode=$(prompt_with_default "Permutation mode (off/grn_edges/all_pairs)" "$default_mode")
        case "${mode,,}" in
            off|no|none|0)
                set_permutation_off_defaults
                break
                ;;
            grn_edges|grn|grn_fdr|fdr|1|on|yes)
                PERMUTE_CHOICE="on"
                PERM_SCOPE="grn_fdr"
                prompt_permutation_detail
                break
                ;;
            all_pairs|all|full|2)
                PERMUTE_CHOICE="on"
                PERM_SCOPE="all"
                prompt_permutation_detail
                break
                ;;
            *)
                echo "[!] Please enter off, grn_edges, or all_pairs."
                default_mode="$mode"
                ;;
        esac
    done
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
    print_step "9" "Local TE Parameters" \
        "Uses Step 2 jobs automatically for TE calculation and LocalTE export" \
        "Only the output-shape settings are shown here" \
        "Advanced worker knobs remain available through environment variables"
    LOCAL_CHUNK_SIZE=$(prompt_with_default "Chunk size" "$default_chunk")
    LOCAL_VALUES_DTYPE=$(prompt_with_default "Values dtype" "$default_dtype")
    LOCAL_SPLIT_EXPORT=$(prompt_with_default "Split export by selectors? (on/off)" "$default_split")
    LOCAL_SPLIT_OUTPUT_DIR=$(prompt_with_default "Split output directory" "$default_split_dir")
    LOCAL_BUFFER_EDGES="$default_buffer"
    LOCAL_EXPORT_WORKERS="$default_export"
    LOCAL_MERGE_WORKERS="$default_merge"
    LOCAL_READ_BATCH="$default_read"
    LOCAL_THREADS_FLAG="$default_threads"
    if is_on "$LOCAL_TE_ADVANCED"; then
        echo ""
        echo "Advanced LocalTE tuning (LOCAL_TE_ADVANCED=on)"
        LOCAL_BUFFER_EDGES=$(prompt_with_default "Buffer edges" "$default_buffer")
        LOCAL_EXPORT_WORKERS=$(prompt_with_default "Export workers (auto uses Step 2 jobs)" "$default_export")
        LOCAL_MERGE_WORKERS=$(prompt_with_default "Merge workers (auto follows export workers)" "$default_merge")
        LOCAL_READ_BATCH=$(prompt_with_default "Read batch rows" "$default_read")
        LOCAL_THREADS_FLAG=$(prompt_with_default "Use threads (on/off)" "$default_threads")
    fi
}

prompt_local_te_choice() {
    local default="${STORE_LOCAL_TE_CHOICE:-$DEFAULT_LOCAL_TE}"
    print_step "9" "Local TE Storage" \
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

prompt_grn_choice() {
    local default="${MAKE_GRN_CHOICE:-$DEFAULT_MAKE_GRN}"
    local default_fdr="${GRN_FDR:-$DEFAULT_GRN_FDR}"
    local default_trim="${TRIM_INDIRECT_CHOICE:-$DEFAULT_TRIM_INDIRECT}"
    local default_trim_threshold="${TRIM_THRESHOLD:-$DEFAULT_TRIM_THRESHOLD}"
    case "${default_trim,,}" in
        tf_gene|tf_gn|tf-gn|all|default|auto) default_trim="on" ;;
    esac
    print_step "10" "GRN / FDR Network" \
        "on = after TE, create FDR-filtered .sif networks for the selected Step 7 outputs" \
        "Indirect trimming is applied only to TF->gene; untrimmed FDR networks are always kept."
    MAKE_GRN_CHOICE=$(prompt_with_default "Create FDR-filtered network? (on/off)" "$default")
    if is_on "$MAKE_GRN_CHOICE"; then
        GRN_FDR=$(prompt_with_default "Network FDR alpha" "$default_fdr")
        TRIM_INDIRECT_CHOICE=$(prompt_with_default "Indirect trimming? (on/off)" "$default_trim")
        case "${TRIM_INDIRECT_CHOICE,,}" in
            tf_gene|tf_gn|tf-gn|all|default|auto) TRIM_INDIRECT_CHOICE="on" ;;
        esac
        if is_on "$TRIM_INDIRECT_CHOICE"; then
            TRIM_THRESHOLD=$(prompt_with_default "Indirect trim cutoff" "$default_trim_threshold")
        else
            TRIM_THRESHOLD="$default_trim_threshold"
        fi
    else
        GRN_FDR="$DEFAULT_GRN_FDR"
        TRIM_INDIRECT_CHOICE="off"
        TRIM_THRESHOLD="$DEFAULT_TRIM_THRESHOLD"
    fi
}

prompt_results_buffer_rows() {
    local default="${RESULTS_BUFFER_ROWS:-$DEFAULT_RESULTS_BUFFER_ROWS}"
    local default_bs="${BATCH_SIZE:-$DEFAULT_BATCH_SIZE}"
    local default_is="${INTERMEDIATE_SAVE_CHOICE:-$DEFAULT_INTERMEDIATE_SAVE}"
    print_step "11" "Runtime / Output" \
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
    echo " 2) Parallel jobs        : $NUM_JOBS"
    echo " 3) Trajectory file      : $TRAJECTORY_FILE"
    echo " 4) Cell-select file     : $CELL_SELECT_FILE"
    echo " 5) History length (k)   : $HISTORY_LENGTH"
    echo " 6) Species              : $SPECIES"
    echo " 7) Network mode         : $(analysis_mode_label)"
    echo "    Kernel estimator     : on (fixed default line)"
    if [ "${PERMUTE_CHOICE,,}" = "on" ]; then
        if [[ "${PERM_SCOPE:-$DEFAULT_PERM_SCOPE}" = "grn_fdr" ]]; then
            echo " 8) Permutation          : grn_edges, count=$PERM_N, candidate_GRN_FDR=${PERM_GRN_FDR:-$DEFAULT_PERM_GRN_FDR}"
        else
            echo " 8) Permutation          : all_pairs, count=$PERM_N"
        fi
    else
        echo " 8) Permutation          : off"
    fi
    echo " 9) Local TE storage     : $STORE_LOCAL_TE_CHOICE"
    if [ "${STORE_LOCAL_TE_CHOICE,,}" = "on" ]; then
        echo "    Local TE export      : chunk=$LOCAL_CHUNK_SIZE, dtype=$LOCAL_VALUES_DTYPE, split=$LOCAL_SPLIT_EXPORT, split_dir=$LOCAL_SPLIT_OUTPUT_DIR"
        echo "                           workers=auto from Step 2 jobs unless overridden"
        if is_on "$LOCAL_TE_ADVANCED" || [ "${LOCAL_BUFFER_EDGES:-$DEFAULT_LOCAL_BUFFER}" != "$DEFAULT_LOCAL_BUFFER" ] || [ "${LOCAL_EXPORT_WORKERS:-$DEFAULT_LOCAL_EXPORT_WORKERS}" != "$DEFAULT_LOCAL_EXPORT_WORKERS" ] || [ "${LOCAL_MERGE_WORKERS:-$DEFAULT_LOCAL_MERGE_WORKERS}" != "$DEFAULT_LOCAL_MERGE_WORKERS" ] || [ "${LOCAL_READ_BATCH:-$DEFAULT_LOCAL_READ_BATCH}" != "$DEFAULT_LOCAL_READ_BATCH" ] || [ "${LOCAL_THREADS_FLAG:-$DEFAULT_LOCAL_THREADS}" != "$DEFAULT_LOCAL_THREADS" ]; then
            echo "                           advanced: buffer=$LOCAL_BUFFER_EDGES, export=$LOCAL_EXPORT_WORKERS, merge=$LOCAL_MERGE_WORKERS, read=$LOCAL_READ_BATCH, threads=$LOCAL_THREADS_FLAG"
        fi
    fi
    if is_on "${MAKE_GRN_CHOICE:-$DEFAULT_MAKE_GRN}"; then
        echo "10) GRN / FDR network    : on, FDR=${GRN_FDR:-$DEFAULT_GRN_FDR}, trim=${TRIM_INDIRECT_CHOICE:-$DEFAULT_TRIM_INDIRECT}, trim_cutoff=${TRIM_THRESHOLD:-$DEFAULT_TRIM_THRESHOLD}"
    else
        echo "10) GRN / FDR network    : off"
    fi
    if [ -n "${RESULTS_BUFFER_ROWS:-}" ] || [ "${INTERMEDIATE_SAVE_CHOICE:-$DEFAULT_INTERMEDIATE_SAVE}" != "$DEFAULT_INTERMEDIATE_SAVE" ] || [ "${BATCH_SIZE:-$DEFAULT_BATCH_SIZE}" != "$DEFAULT_BATCH_SIZE" ]; then
        echo "11) Runtime/output       : buffer=${RESULTS_BUFFER_ROWS:-auto}, intermediate=${INTERMEDIATE_SAVE_CHOICE:-$DEFAULT_INTERMEDIATE_SAVE}, batch=${BATCH_SIZE:-$DEFAULT_BATCH_SIZE}"
    else
        echo "11) Runtime/output       : defaults"
    fi
    echo "---------------------------"
}

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

is_on() {
    case "${1,,}" in
        on|1|yes|true) return 0 ;;
        *) return 1 ;;
    esac
}

should_trim_te_file() {
    local trim_mode="${1,,}"
    local te_file="$2"
    case "$trim_mode" in
        on|1|yes|true|tf_gene|tf_gn|tf-gn|default|auto)
            [[ "$te_file" == "TE_TF_GN.parquet" ]]
            return
            ;;
        all)
            echo "Warning: TENET_TRIM_INDIRECT=all is no longer supported; trimming only TF->gene." >&2
            [[ "$te_file" == "TE_TF_GN.parquet" ]]
            return
            ;;
        off|0|no|false|none|"")
            return 1
            ;;
        *)
            echo "Warning: unknown TENET_TRIM_INDIRECT=${trim_mode}; skipping indirect trimming." >&2
            return 1
            ;;
    esac
}

build_grn_networks() {
    if ! is_on "${TENET_MAKE_GRN:-$DEFAULT_MAKE_GRN}"; then
        echo "--- GRN/FDR network generation: OFF ---"
        return 0
    fi

    local grn_fdr="${TENET_GRN_FDR:-$DEFAULT_GRN_FDR}"
    local grn_te_cutoff="${TENET_GRN_TE_CUTOFF:-$DEFAULT_GRN_TE_CUTOFF}"
    local trim_mode="${TENET_TRIM_INDIRECT:-$DEFAULT_TRIM_INDIRECT}"
    local trim_threshold="${TENET_TRIM_THRESHOLD:-$DEFAULT_TRIM_THRESHOLD}"
    local grn_threads="${TENET_GRN_THREADS:-${2:-1}}"
    local target_files=()
    local te_file base outdir sif_file trimmed_file

    for te_file in \
        TE_TF_GN.parquet \
        TE_TF_PK.parquet \
        TE_PK_GN.parquet \
        TE_PK_PK.parquet \
        TE_GN_GN.parquet \
        TE_all_features.parquet; do
        if [[ -f "$te_file" ]]; then
            target_files+=("$te_file")
        fi
    done

    if [[ ${#target_files[@]} -eq 0 ]]; then
        echo "[GRN] Warning: network generation requested, but no split TE parquet files were found."
        return 0
    fi

    echo "[GRN] FDR network generation: ON (alpha=${grn_fdr}, indirect_trim=${trim_mode})"
    printf "te_parquet\tnetwork_dir\tsif_before_trim\ttrimmed_sif\n" > grn_network_outputs.tsv
    for te_file in "${target_files[@]}"; do
        base="${te_file%.parquet}"
        outdir="network_${base}_fdr${grn_fdr}"
        echo "[GRN] Building FDR network for ${te_file}; untrimmed .sif will be kept."
        if ! run_with_stage "MAKE_GRN_${base}" "$PYTHON" -m code.network_from_te_parquet "$te_file" \
            --outdir "$outdir" \
            --output-prefix "$base" \
            --fdr "$grn_fdr" \
            --te-cutoff "$grn_te_cutoff" \
            --threads "$grn_threads"; then
            echo "Error in network_from_te_parquet.py for ${te_file}. Exiting."
            exit 1
        fi

        sif_file=$(find "$outdir" -maxdepth 1 -type f -name "${base}_fdr*.sif" ! -name '*.trimIndirect*.sif' | sort | head -n 1 || true)
        trimmed_file=""
        if [[ -n "$sif_file" && -f "$sif_file" ]] && should_trim_te_file "$trim_mode" "$te_file"; then
            echo "[GRN] Running indirect trimming for ${te_file}; source .sif remains: ${sif_file}"
            if ! run_with_stage "TRIM_INDIRECT_${base}" "$PYTHON" "$CODE_DIR/trim_indirect.py" "$sif_file" "$trim_threshold"; then
                echo "Error in trim_indirect.py for ${sif_file}. Exiting."
                exit 1
            fi
            trimmed_file="${sif_file/.sif/.trimIndirect${trim_threshold}.sif}"
        elif [[ -n "$sif_file" && -f "$sif_file" ]]; then
            echo "[GRN] Skipping indirect trimming for ${te_file} (TENET_TRIM_INDIRECT=${trim_mode}); kept untrimmed .sif: ${sif_file}"
        fi
        printf "%s\t%s\t%s\t%s\n" "$te_file" "$outdir" "$sif_file" "$trimmed_file" >> grn_network_outputs.tsv
    done
    echo "[GRN] Output manifest: grn_network_outputs.tsv"
}

build_replay_command() {
    local cmd=""
    local var val default

    append_env_if_nondefault() {
        var="$1"
        default="$2"
        val="${!var-}"
        if [[ -n "${val:-}" && "$val" != "$default" ]]; then
            cmd+="${var}=$(printf '%q' "$val") "
        fi
    }

    append_env_if_set() {
        var="$1"
        val="${!var-}"
        if [[ -n "${val:-}" ]]; then
            cmd+="${var}=$(printf '%q' "$val") "
        fi
    }

    local mode="${MODE_CODE:-${7:-1}}"
    local replay_mode="$mode"
    case "${PAIR_MODE:-}" in
        all_pair)
            case "$mode" in
                0) replay_mode="7" ;;
                1) replay_mode="8" ;;
                2) replay_mode="9" ;;
            esac
            ;;
        *)
            case "${TENET_PAIR_MODE:-default}:$mode" in
                gene_only:0) replay_mode="7" ;;
                gene_only:1) replay_mode="8" ;;
                all_feature:2) replay_mode="9" ;;
            esac
            ;;
    esac

    append_env_if_nondefault PYTHON "python3"
    append_env_if_nondefault MONITOR_GPU "off"
    append_env_if_nondefault TENET_OUTPUT_DIR "$REPO_ROOT/output"
    append_env_if_nondefault TENET_INTERMEDIATE_SAVE "off"
    append_env_if_nondefault TE_BATCH_SIZE "100"
    append_env_if_nondefault TENET_MAKE_GRN "on"
    append_env_if_nondefault TENET_GRN_FDR "0.01"
    append_env_if_nondefault TENET_GRN_TE_CUTOFF "0.0"
    append_env_if_nondefault TENET_TRIM_INDIRECT "on"
    append_env_if_nondefault TENET_TRIM_THRESHOLD "-0.01"
    if [[ -n "${TENET_GRN_THREADS:-}" && "${TENET_GRN_THREADS}" != "${2:-}" ]]; then
        cmd+="TENET_GRN_THREADS=$(printf '%q' "$TENET_GRN_THREADS") "
    fi
    if [[ -n "${TENET_PERM_CANDIDATE_GRN_FDR:-}" ]]; then
        append_env_if_set TENET_PERM_CANDIDATE_GRN_FDR
    fi
    append_env_if_set TENET_PERM_CANDIDATE_TE_CUTOFF

    append_env_if_nondefault TE_KERNEL_GROUPED_STATE "on"
    append_env_if_nondefault TE_KERNEL_GROUPED_STATE_CODED "on"
    append_env_if_nondefault TE_KERNEL_GROUPED_STATE_FACTORIZED "auto"
    append_env_if_nondefault TE_KERNEL_GROUPED_STATE_SOURCE_PREFIX "auto"
    append_env_if_nondefault TE_KERNEL_SOURCE_CODE_INTEGER "auto"
    append_env_if_nondefault TE_KERNEL_GROUPED_STATE_SPARSE "off"
    append_env_if_nondefault TE_KERNEL_GROUPED_STATE_2D_RANGE "on"
    append_env_if_nondefault TE_KERNEL_PERM_TABLE_SAMPLER "off"
    append_env_if_nondefault TE_KERNEL_PERM_TABLE_MAX_CELLS ""
    append_env_if_nondefault TE_PERM_KERNEL_GROUPED_SOURCES_PER_CHUNK "2"
    append_env_if_nondefault TE_PERM_SORT_WORK "on"

    local store_local="${20:-off}"
    if [[ "${store_local,,}" == "on" || "$store_local" == "1" ]]; then
        append_env_if_nondefault LOCAL_CHUNK_SIZE "300"
        append_env_if_nondefault LOCAL_BUFFER_EDGES "10000"
        append_env_if_nondefault LOCAL_EXPORT_WORKERS "auto"
        append_env_if_nondefault LOCAL_MERGE_WORKERS "auto"
        append_env_if_nondefault LOCAL_READ_BATCH "8192"
        append_env_if_nondefault LOCAL_THREADS_FLAG "on"
        append_env_if_nondefault LOCAL_VALUES_DTYPE "float16"
        append_env_if_nondefault LOCAL_SPLIT_EXPORT "on"
        append_env_if_nondefault LOCAL_SPLIT_OUTPUT_DIR "local_te_split_chunks"
        append_env_if_nondefault LOCAL_TE_ADVANCED "off"
    fi

    local opt_args=(
        "${8:-none}"
        "${9:-kernel}"
        "${10:-none}"
        "${11:-0}"
        "${12:-0}"
        "${13:-off}"
        "${14:-100}"
        "${15:-0}"
        "${16:-0}"
        "${17:-off}"
        "${18:-0.05}"
        "${19:-0.01}"
        "${20:-off}"
        "${21:-1}"
        "${22:-100}"
        "${23:-42}"
        "${24:-}"
    )
    local last_opt=0
    local permute="${13:-off}"
    if [[ "${permute,,}" == "on" || "$permute" == "1" ]]; then
        last_opt=7
        [[ "${15:-0}" != "0" ]] && last_opt=8
        [[ "${16:-0}" != "0" ]] && last_opt=9
        if [[ "${17:-off}" != "off" ]]; then
            last_opt=11
        elif [[ "${19:-0.01}" != "0.01" ]]; then
            last_opt=12
        fi
    fi
    if [[ "${store_local,,}" == "on" || "$store_local" == "1" ]]; then
        (( last_opt < 13 )) && last_opt=13
    fi
    if [[ -n "${24:-}" ]]; then
        last_opt=17
    fi

    cmd+="$(printf '%q' "$0")"
    local replay_args=("${1:-}" "${2:-}" "${3:-}" "${4:-}" "${5:-}" "${6:-}" "$replay_mode")
    local arg
    for arg in "${replay_args[@]}"; do
        cmd+=" $(printf '%q' "$arg")"
    done
    local i
    for ((i = 0; i < last_opt; i++)); do
        cmd+=" $(printf '%q' "${opt_args[$i]}")"
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
    echo "[Run] Replay files saved:"
    echo "  $txt"
    echo "  $sh"
    echo "[Run] Replay command:"
    echo "  $cmd"
}

# Function to display usage
usage() {
    cat <<USAGE
Usage: $0 <input_file> <num_jobs> <trajectory_file> <cell_select_file> <history_k> <species> <network_mode> [none] [kernel] [none] [0] [0] [permute] [perm_n] [0] [0] [perm_fdr] [perm_q_alpha] [perm_alpha] [store_local_te] [1] [100] [42] [results_buffer_rows]

  <network_mode>:
    Default: 1 = TENET_Plus full
    All modes:
      0 = RNA-only TENET_TF
      1 = TENET_Plus full: TF->gene, TF->peak, peak->gene
      2 = TF->gene
      3 = TF->peak
      4 = TF->gene + TF->peak
      5 = peak->gene (cis)
      6 = peak->peak (cis)
      7 = RNA-only every gene -> every gene
      8 = every gene -> every gene
      9 = every gene/peak -> every gene/peak
    Name aliases and a0/a1/a2 remain backward-compatible.

Optional:
  screen_mode: kernel only in this default line.
  modality/refine/gene-filter/time-sampling:
               not part of this kernel-only package. Prepare inputs before running.
  permute:     on | off                     (default: off)
  perm_n:      number of permutations       (default: 100)
  perm_fdr:    on | off                     (default: off)
  perm_q_alpha:FDR q threshold              (default: 0.05)
  perm_alpha:  p-value threshold            (default: 0.01)
  store_local_te: on | off                 (default: off)
  results_buffer_rows: rows to buffer in TE before flushing batch parquet files (default auto: 200000, or 5000 when storing LocalTE)

Environment:
  TENET_PAIR_MODE     : default | gene_only | all_feature.
                        Kept for backward-compatible non-interactive all-pair runs.
  TENET_INTERMEDIATE_SAVE: on | off (default: off).
  TENET_MAKE_GRN    : on | off (default: on). Create FDR-filtered network files
                      from the Step 7 TE outputs after Matrix_generate.
  TENET_GRN_FDR     : FDR alpha for generated networks (default: 0.01).
  TENET_TRIM_INDIRECT: on | off (default: on). When on, run trim_indirect only
                      on TF->gene. Untrimmed FDR networks are always kept.
  TENET_TRIM_THRESHOLD: indirect-trimming cutoff (default: -0.01).
  TENET_PERM_CANDIDATE_GRN_FDR:
                        when permutation is on, restrict candidates to GRN FDR
                        edges at this alpha.
  TENET_PERM_CANDIDATE_TE_CUTOFF:
                        optional TE cutoff for GRN-FDR permutation candidates.
USAGE
    exit 0
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
fi

mkdir -p "$OUTPUT_DIR"

if [ $# -eq 0 ]; then
    echo "Interactive mode: press Enter to accept defaults shown in brackets."
    echo "Path prompts start in the repository root ($REPO_ROOT); try Tab for completion."

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
    prompt_jobs
    prompt_trajectory_file
    prompt_cell_select_file
    prompt_history_length
    prompt_species
    prompt_analysis_mode
    update_defaults_after_mode
    TENET_GENE_FILTER="none"
    MODALITY_CHOICE="none"
    SCREEN_MODE="kernel"
    REFINE_METHOD="none"
    REFINE_TOPK="0"
    REFINE_TOP_PCT="0"
    TIME_STRIDE="1"
    TIME_PCT="100"
    TIME_SEED="42"
    prompt_permutation_toggle
    prompt_local_te_choice
    prompt_grn_choice
    INTERMEDIATE_SAVE_CHOICE="${INTERMEDIATE_SAVE_CHOICE:-$DEFAULT_INTERMEDIATE_SAVE}"
    RESULTS_BUFFER_ROWS="${RESULTS_BUFFER_ROWS:-$DEFAULT_RESULTS_BUFFER_ROWS}"
    BATCH_SIZE="${BATCH_SIZE:-$DEFAULT_BATCH_SIZE}"

    while true; do
        show_summary
        read -rp "Edit step (1-11) or press Enter to continue: " EDIT_CHOICE
        EDIT_CHOICE=$(trim_whitespace "$EDIT_CHOICE")
        case "$EDIT_CHOICE" in
            "")
                break
                ;;
            1)
                prompt_input_matrix
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
            7)
                prompt_analysis_mode
                update_defaults_after_mode
                ;;
            8)
                prompt_permutation_toggle
                ;;
            9)
                prompt_local_te_choice
                ;;
            10)
                prompt_grn_choice
                ;;
            11)
                prompt_results_buffer_rows
                ;;
            *)
                echo "[!] Unknown step selection."
                ;;
        esac
    done

    TENET_MAKE_GRN="${MAKE_GRN_CHOICE:-$DEFAULT_MAKE_GRN}"
    TENET_GRN_FDR="${GRN_FDR:-$DEFAULT_GRN_FDR}"
    TENET_GRN_TE_CUTOFF="${DEFAULT_GRN_TE_CUTOFF}"
    TENET_TRIM_INDIRECT="${TRIM_INDIRECT_CHOICE:-$DEFAULT_TRIM_INDIRECT}"
    TENET_TRIM_THRESHOLD="${TRIM_THRESHOLD:-$DEFAULT_TRIM_THRESHOLD}"
    TENET_GRN_THREADS="${DEFAULT_GRN_THREADS:-$NUM_JOBS}"

    set -- "$INPUT_FILE" "$NUM_JOBS" "$TRAJECTORY_FILE" "$CELL_SELECT_FILE" \
        "$HISTORY_LENGTH" "$SPECIES" "$MODE_CODE" "none" "kernel" \
        "none" "0" "0" "$PERMUTE_CHOICE" "$PERM_N" \
        "$PERM_TOPK" "$PERM_TOP_PCT" "$PERM_FDR" "$PERM_Q_ALPHA" "$PERM_ALPHA" \
        "$STORE_LOCAL_TE_CHOICE" "1" "100" "42" \
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

# Accept the compact analysis-mode aliases used by the interactive prompt while
# preserving the older numeric $7 + TENET_PAIR_MODE environment contract.
ANALYSIS_ARG="${7,,}"
case "$ANALYSIS_ARG" in
  rna_tf|tf_rna|tf|rna|tenet_tf)
    PAIR_MODE="default"
    MODE_CODE="0"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  plus_full|full|tenet_plus_full)
    PAIR_MODE="default"
    MODE_CODE="1"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  tf_gene|tf_gn|tf-gene|gene_grn|grn)
    PAIR_MODE="default"
    MODE_CODE="2"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  tf_peak|tf_pk|tf-peak)
    PAIR_MODE="default"
    MODE_CODE="3"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  tf_targets|tf_gene_peak|tf_gn_pk|tf-gene-peak|tf_gene+peak)
    PAIR_MODE="default"
    MODE_CODE="4"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  peak_gene|peak_gn|peak-gene|cis_peak_gene)
    PAIR_MODE="default"
    MODE_CODE="5"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  peak_peak|peak_pk|peak-peak|cis_peak_peak)
    PAIR_MODE="default"
    MODE_CODE="6"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  7|a0|all_gene_rna|rna_all_gene)
    PAIR_MODE="all_pair"
    MODE_CODE="0"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  8|a1|all_gene|all_gene_plus|gene_all_gene)
    PAIR_MODE="all_pair"
    MODE_CODE="1"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  9|a2|all_feature|all_features|feature_all_feature)
    PAIR_MODE="all_pair"
    MODE_CODE="2"
    set -- "${@:1:6}" "$MODE_CODE" "${@:8}"
    ;;
  *)
    MODE_CODE="${MODE_CODE:-$7}"
    ;;
esac

# Normalise pair mode for downstream Python preprocess scripts
REQUESTED_PAIR_MODE="${PAIR_MODE:-${TENET_PAIR_MODE:-}}"
if [ -n "$REQUESTED_PAIR_MODE" ]; then
  if [ "$REQUESTED_PAIR_MODE" = "all_pair" ]; then
    case "${MODE_CODE:-${7:-0}}" in
      0) TENET_PAIR_MODE="gene_only" ;;      # RNA-only, all gene x gene
      1) TENET_PAIR_MODE="gene_only" ;;      # TENET_Plus, all gene x gene
      2) TENET_PAIR_MODE="all_feature" ;;    # TENET_Plus, all features x features
      *) TENET_PAIR_MODE="gene_only" ;;
    esac
  else
    case "$REQUESTED_PAIR_MODE" in
      default|gene_only|all_feature) TENET_PAIR_MODE="$REQUESTED_PAIR_MODE" ;;
      *) TENET_PAIR_MODE="default" ;;
    esac
  fi
elif [ -z "${TENET_PAIR_MODE:-}" ]; then
  TENET_PAIR_MODE="default"
fi

if [[ "${TENET_GENE_FILTER:-none}" != "none" ]]; then
  echo "Error: gene filtering is separated from the kernel-only default line (requested TENET_GENE_FILTER=${TENET_GENE_FILTER})."
  echo "Pre-filter the matrix before running this kernel package."
  exit 1
fi
REQUESTED_MODALITY="${8:-none}"
if [[ "$REQUESTED_MODALITY" != "none" && "$REQUESTED_MODALITY" != "skip" ]]; then
  echo "Error: modality preprocessing is separated from the kernel-only default line (requested modality=${REQUESTED_MODALITY})."
  echo "Provide an already prepared matrix before running this kernel package."
  exit 1
fi
REQUESTED_SCREEN_MODE="${9:-kernel}"
if [[ "${REQUESTED_SCREEN_MODE,,}" != "kernel" ]]; then
  echo "Error: this default TENETPLUS line is kernel-only (requested screen_mode=${REQUESTED_SCREEN_MODE})."
  echo "This standalone package only supports kernel TE."
  exit 1
fi
REQUESTED_REFINE_METHOD="${10:-none}"
if [[ "${REQUESTED_REFINE_METHOD,,}" != "none" ]]; then
  echo "Error: refinement is not part of the kernel-only default line (requested refine_method=${REQUESTED_REFINE_METHOD})."
  echo "Run refinement outside this standalone kernel package if needed."
  exit 1
fi
REQUESTED_TIME_STRIDE="${21:-1}"
REQUESTED_TIME_PCT="${22:-100}"
REQUESTED_TIME_SEED="${23:-42}"
if [[ "$REQUESTED_TIME_STRIDE" != "1" || "$REQUESTED_TIME_PCT" != "100" ]]; then
  echo "Error: time sampling is separated from the kernel-only default line (requested stride=${REQUESTED_TIME_STRIDE}, pct=${REQUESTED_TIME_PCT}, seed=${REQUESTED_TIME_SEED})."
  echo "Subsample timepoints before running this standalone kernel package."
  exit 1
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

TENET_MAKE_GRN="${MAKE_GRN_CHOICE:-${TENET_MAKE_GRN:-$DEFAULT_MAKE_GRN}}"
TENET_GRN_FDR="${GRN_FDR:-${TENET_GRN_FDR:-$DEFAULT_GRN_FDR}}"
TENET_GRN_TE_CUTOFF="${TENET_GRN_TE_CUTOFF:-$DEFAULT_GRN_TE_CUTOFF}"
TENET_TRIM_INDIRECT="${TRIM_INDIRECT_CHOICE:-${TENET_TRIM_INDIRECT:-$DEFAULT_TRIM_INDIRECT}}"
TENET_TRIM_THRESHOLD="${TRIM_THRESHOLD:-${TENET_TRIM_THRESHOLD:-$DEFAULT_TRIM_THRESHOLD}}"
TENET_GRN_THREADS="${TENET_GRN_THREADS:-$2}"

write_replay_logs "$@"

chmod +x "$CODE_DIR/memory_check.sh" || true
CURRENT_STAGE_FILE="$OUTPUT_DIR/current_stage.txt"
echo "INIT" > "$CURRENT_STAGE_FILE"
export TENET_STAGE_FILE="$CURRENT_STAGE_FILE"
# GPU monitor toggle (env override: MONITOR_GPU=on|off); default off for lower overhead
MONITOR_GPU="${MONITOR_GPU:-off}"
"$CODE_DIR/memory_check.sh" $$ "$CURRENT_STAGE_FILE" "$MONITOR_GPU" &
MONITOR_PID=$!
echo "[Run] Memory monitor started: pid=${MONITOR_PID}, gpu=${MONITOR_GPU}"

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
    echo "=== TENET+ Kernel Run ==="
    echo "[Setup] Extracting peak names from matrix header."
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
    echo "=== TENET Kernel Run: RNA-only TF mode ==="
fi

# Clean up temporary CSV file if it was created
if [ "$FILE_EXT" == "parquet" ]; then
    rm -f "$TEMP_CSV"
fi

# Check if cell_select_file contains only one cell or cells are sorted by pseudotime order
if grep -qv '^1$' "$4" || ! sort -gc "$3" || ! file "$1" | grep -q 'parquet'; then
    set_stage "MATRIX_PROCESSING"
    echo "[Matrix] Filtering selected cells, sorting by pseudotime, and writing Parquet input."
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
    echo "[Matrix] Input already appears selected/sorted; using it as-is."
    matrix=$1
    trajectory=$3
    cell_select=$4
fi


GENE_FILTER_MODE="${TENET_GENE_FILTER:-none}"
if [ "$GENE_FILTER_MODE" != "none" ]; then
    echo "Error: gene filtering is separated from the kernel-only default line (requested TENET_GENE_FILTER=${GENE_FILTER_MODE})."
    echo "Pre-filter the matrix before running this kernel package."
    exit 1
fi


# Always regenerate gene_names from the current matrix
set_stage "GENE_NAMES_REGEN"
echo "[Setup] Writing gene_names from matrix header."
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
# Modality-specific preprocessing has been separated from the kernel-only line.
MODALITY="${8:-none}"
if [[ "$MODALITY" != "none" && "$MODALITY" != "skip" ]]; then
  echo "Error: modality preprocessing is separated from the kernel-only default line (requested modality=${MODALITY})."
  echo "Provide an already prepared matrix before running this kernel package."
  exit 1
fi

# If a precomputed transposed matrix exists, remove it to regenerate fresh
if [ -f cell_gene_trsps.parquet ]; then
  echo "[Matrix] Removing stale cell_gene_trsps.parquet before regeneration."
  rm -f cell_gene_trsps.parquet
fi

set_stage "MODALITY_SKIP"
echo "[Setup] Modality preprocessing skipped; kernel package expects a prepared matrix."

if [ "$7" -ne 0 ]; then
  echo "[Pairs] Building selected source-target pairs for mode ${7}."
    if [ "${TENET_PAIR_MODE:-default}" = "default" ]; then
      if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "PREPROCESS_TENET_PLUS" "$PYTHON" -m code.PreProcessScript_TE_Plus $6 $7 "${matrix}"; then
          echo "Error in PreProcessScript_TE_Plus.py. Exiting."
          exit 1
      fi
    else
      echo "[Pairs] Using implicit all-pair mode in TE core (${TENET_PAIR_MODE}); no all_pairs.csv materialization."
    fi
else
  echo "[Pairs] Building RNA-only TF source-target pairs."
  if [ "${TENET_PAIR_MODE:-default}" = "default" ]; then
    if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "PREPROCESS_TENET_TF" "$PYTHON" -m code.PreProcessScript_TENET_TF $6; then
        echo "Error in PreProcessScript_TENET_TF.py. Exiting."
        exit 1
    fi
  else
    echo "[Pairs] Using implicit all-pair mode in TE core (${TENET_PAIR_MODE}); no all_pairs.csv materialization."
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
if [[ "${SCREEN_MODE,,}" != "kernel" ]]; then
  echo "Error: this default TENETPLUS line is kernel-only (requested screen_mode=${SCREEN_MODE})."
  echo "This standalone package only supports kernel TE."
  exit 1
fi
if [[ "${REFINE_METHOD,,}" != "none" ]]; then
  echo "Error: refinement is not part of the kernel-only default line (requested refine_method=${REFINE_METHOD})."
  echo "Run refinement outside this standalone kernel package if needed."
  exit 1
fi
echo "[TE] Kernel estimator: exact grouped-state fast path enabled."

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
  LOCAL_ARGS=(--store_local_te --localte_codec "$LOCAL_TE_CODEC")
  echo "[TE] LocalTE storage: ON (codec=${LOCAL_TE_CODEC}). This adds per-pair timepoint arrays and can be I/O/compression bound."
else
  echo "[TE] LocalTE storage: OFF."
fi

# Reject time-sampling placeholders in the kernel-only line.
TIME_STRIDE="${21:-1}"
TIME_PCT="${22:-100}"
TIME_SEED="${23:-42}"
TIME_ARGS=()
if [[ "$TIME_STRIDE" != "1" || "$TIME_PCT" != "100" ]]; then
  echo "Error: time sampling is separated from the kernel-only default line (requested stride=${TIME_STRIDE}, pct=${TIME_PCT}, seed=${TIME_SEED})."
  echo "Subsample timepoints before running this standalone kernel package."
  exit 1
fi

# TE results buffer rows (optional; controls when intermediate batches flush to disk)
RESULTS_BUFFER_ROWS="${24:-}"
BUFFER_ARGS=()
if [[ -n "${RESULTS_BUFFER_ROWS:-}" ]]; then
  if [[ "$RESULTS_BUFFER_ROWS" =~ ^[0-9]+$ ]] && [ "$RESULTS_BUFFER_ROWS" -gt 0 ]; then
	    BUFFER_ARGS+=(--results_buffer_rows "$RESULTS_BUFFER_ROWS")
	    echo "[TE] Output buffer rows: $RESULTS_BUFFER_ROWS"
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
	    echo "[TE] Intermediate save: ON"
	    ;;
	  off|0|no)
	    echo "[TE] Intermediate save: OFF"
	    ;;
  *)
    # Fallback to default
    if [[ "${DEFAULT_INTERMEDIATE_SAVE,,}" == "on" || "${DEFAULT_INTERMEDIATE_SAVE,,}" == "yes" ]]; then
	      ENABLE_INTERMEDIATE_ARGS=(--enable_intermediate_save)
	      echo "[TE] Intermediate save: ON (default)"
	    else
	      echo "[TE] Intermediate save: OFF (default)"
	    fi
    ;;
esac

# TE batch size (optional; controls sources per target chunk)
TE_BATCH_SIZE_VALUE="${BATCH_SIZE:-${TE_BATCH_SIZE:-}}"
TE_BATCH_ARGS=()
if [[ -n "${TE_BATCH_SIZE_VALUE:-}" ]]; then
  if [[ "$TE_BATCH_SIZE_VALUE" =~ ^[0-9]+$ ]] && [ "$TE_BATCH_SIZE_VALUE" -gt 0 ]; then
	    TE_BATCH_ARGS+=(--batch_size "$TE_BATCH_SIZE_VALUE")
	    echo "[TE] Fallback batch size: $TE_BATCH_SIZE_VALUE"
  else
    echo "Warning: Ignoring invalid TE batch size '$TE_BATCH_SIZE_VALUE' (expected positive integer)."
  fi
fi

# Kernel-only TE execution path.
if [[ "$PERMUTE" == "on" || "$PERMUTE" == "ON" || "$PERMUTE" == "1" ]]; then
  echo "[Permutation] ON: n=${PERM_N}, topk=${PERM_TOPK}, pct=${PERM_TOP_PCT}, fdr=${PERM_FDR}"
  PERM_ARGS=(--permute --perm_n "$PERM_N")
  PERM_CANDIDATE_GRN_FDR="${TENET_PERM_CANDIDATE_GRN_FDR:-}"
  if [[ -n "${PERM_CANDIDATE_GRN_FDR:-}" && "${PERM_CANDIDATE_GRN_FDR:-0}" != "0" ]]; then
    PERM_ARGS+=(--perm_candidate_grn_fdr "$PERM_CANDIDATE_GRN_FDR")
    PERM_ARGS+=(--perm_candidate_te_cutoff "${TENET_PERM_CANDIDATE_TE_CUTOFF:-0.0}")
    echo "--- Permutation candidate filter: grn_fdr=${PERM_CANDIDATE_GRN_FDR}, te_cutoff=${TENET_PERM_CANDIDATE_TE_CUTOFF:-0.0} ---"
  fi
  if [[ "$PERM_TOPK" != "0" ]]; then PERM_ARGS+=(--permute_topk_per_target "$PERM_TOPK"); fi
  if [[ "$PERM_TOP_PCT" != "0" ]]; then PERM_ARGS+=(--permute_top_pct "$PERM_TOP_PCT"); fi
  if [[ "$PERM_FDR" == "on" || "$PERM_FDR" == "ON" || "$PERM_FDR" == "1" ]]; then
    PERM_ARGS+=(--use_fdr --perm_q_alpha "$PERM_Q_ALPHA")
  else
    PERM_ARGS+=(--perm_alpha "$PERM_ALPHA")
  fi
  if ! run_with_stage "TE_KERNEL_PERMUTE" "$PYTHON" -m code.runTE_for_py_python_batch all_pairs.csv "$2" "$HISTORY_LENGTH" \
      "${ENABLE_INTERMEDIATE_ARGS[@]}" \
      --mode kernel \
      "${TE_BATCH_ARGS[@]}" \
      --pair_mode "$TENET_PAIR_MODE" \
      "${PERM_ARGS[@]}" \
      "${BUFFER_ARGS[@]}" \
      "${LOCAL_ARGS[@]}" \
      "${TIME_ARGS[@]}"; then
    echo "Error in kernel TE permutation stage. Exiting."
    exit 1
  fi
else
  if ! run_with_stage "TE_KERNEL" "$PYTHON" -m code.runTE_for_py_python_batch all_pairs.csv "$2" "$HISTORY_LENGTH" \
      "${ENABLE_INTERMEDIATE_ARGS[@]}" \
      --mode kernel \
      "${TE_BATCH_ARGS[@]}" \
      --pair_mode "$TENET_PAIR_MODE" \
      "${BUFFER_ARGS[@]}" \
      "${LOCAL_ARGS[@]}" \
      "${TIME_ARGS[@]}"; then
    echo "Error in kernel TE stage. Exiting."
    exit 1
  fi
fi

rm -f TE_TF_GN.parquet TE_TF_PK.parquet TE_PK_GN.parquet TE_PK_PK.parquet TE_GN_GN.parquet TE_all_features.parquet
if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "MATRIX_GENERATE" "$PYTHON" -m code.Matrix_generate TE_result_all.parquet $6 $7; then
    echo "Error in Matrix_generate.py. Exiting."
    exit 1
fi

build_grn_networks

if [[ "$STORE_LOCAL_TE" == "on" || "$STORE_LOCAL_TE" == "ON" || "$STORE_LOCAL_TE" == "1" ]]; then
    echo "[LocalTE] Exporting LocalTE chunks (${LOCAL_TE_CHUNK_SIZE} timepoints per chunk)."
    set_stage "LOCALTE_EXPORT"
    LOCAL_TE_EFFECTIVE_WORKERS="$LOCAL_TE_EXPORT_WORKERS"
    if [[ "${LOCAL_TE_EFFECTIVE_WORKERS,,}" == "auto" ]]; then
        LOCAL_TE_EFFECTIVE_WORKERS=$("$PYTHON" - <<PY
jobs = int("${2:-1}")
print(max(1, min(jobs, 16)))
PY
)
    fi
    if ! [[ "$LOCAL_TE_EFFECTIVE_WORKERS" =~ ^[0-9]+$ ]]; then
        LOCAL_TE_EFFECTIVE_WORKERS=1
    fi
    if [[ -z "${LOCAL_TE_MERGE_WORKERS:-}" || "${LOCAL_TE_MERGE_WORKERS,,}" == "auto" ]]; then
        LOCAL_TE_EFFECTIVE_MERGE_WORKERS=$("$PYTHON" - <<PY
jobs = int("${2:-1}")
print(max(1, min(jobs, 32)))
PY
)
    else
        LOCAL_TE_EFFECTIVE_MERGE_WORKERS="$LOCAL_TE_MERGE_WORKERS"
    fi
    echo "[LocalTE] Export workers=${LOCAL_TE_EFFECTIVE_WORKERS}, merge_workers=${LOCAL_TE_EFFECTIVE_MERGE_WORKERS}, merge_parts=${LOCAL_TE_MERGE_PARTS}"
    LOCAL_TE_INPUT_PARQUET="TE_result_all.parquet"
    if [[ ! -f "$LOCAL_TE_INPUT_PARQUET" && -f "TE_fast.parquet" ]]; then
        LOCAL_TE_INPUT_PARQUET="TE_fast.parquet"
    fi
    if [[ "${LOCAL_TE_SPLIT_EXPORT,,}" == "on" ]]; then
        # Build selector-based inputs if selector files exist
        SEL_INPUTS=()
        if [[ -f TE_TF_GN.parquet ]]; then
            SEL_INPUTS+=("${LOCAL_TE_INPUT_PARQUET}=TE_TF_GN.parquet")
        fi
        if [[ -f TE_TF_PK.parquet ]]; then
            SEL_INPUTS+=("${LOCAL_TE_INPUT_PARQUET}=TE_TF_PK.parquet")
        fi
        if [[ -f TE_PK_GN.parquet ]]; then
            SEL_INPUTS+=("${LOCAL_TE_INPUT_PARQUET}=TE_PK_GN.parquet")
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
                --workers "$LOCAL_TE_EFFECTIVE_WORKERS" \
                --merge_workers "$LOCAL_TE_EFFECTIVE_MERGE_WORKERS" \
                --merge_parts "$LOCAL_TE_MERGE_PARTS"; then
                echo "Error in localte_chunk_export.py (split mode). Exiting."
                exit 1
            fi
        else
            echo "[localte_chunk_export] No selector files found; exporting all edges to local_te_chunks"
            if ! run_with_stage "LOCALTE_EXPORT_ALL" "$PYTHON" -m code.localte_chunk_export \
                --input "$LOCAL_TE_INPUT_PARQUET" \
                --output_dir local_te_chunks \
                --chunk_size "$LOCAL_TE_CHUNK_SIZE" \
                --buffer_edges "$LOCAL_TE_BUFFER_EDGES" \
                --read_batch_rows "$LOCAL_TE_READ_BATCH_ROWS" \
                --use_threads "$LOCAL_TE_USE_THREADS" \
                --values_dtype "$LOCAL_TE_VALUES_DTYPE" \
                --workers "$LOCAL_TE_EFFECTIVE_WORKERS" \
                --merge_workers "$LOCAL_TE_EFFECTIVE_MERGE_WORKERS" \
                --merge_parts "$LOCAL_TE_MERGE_PARTS" \
                --scores_output TE_result_scores.parquet; then
                echo "Error in localte_chunk_export.py. Exiting."
                exit 1
            fi
        fi
    else
        if ! run_with_stage "LOCALTE_EXPORT_ALL" "$PYTHON" -m code.localte_chunk_export \
            --input "$LOCAL_TE_INPUT_PARQUET" \
            --output_dir local_te_chunks \
            --chunk_size "$LOCAL_TE_CHUNK_SIZE" \
            --buffer_edges "$LOCAL_TE_BUFFER_EDGES" \
            --read_batch_rows "$LOCAL_TE_READ_BATCH_ROWS" \
            --use_threads "$LOCAL_TE_USE_THREADS" \
            --values_dtype "$LOCAL_TE_VALUES_DTYPE" \
            --workers "$LOCAL_TE_EFFECTIVE_WORKERS" \
            --merge_workers "$LOCAL_TE_EFFECTIVE_MERGE_WORKERS" \
            --merge_parts "$LOCAL_TE_MERGE_PARTS" \
            --scores_output TE_result_scores.parquet; then
            echo "Error in localte_chunk_export.py. Exiting."
            exit 1
        fi
    fi
fi
