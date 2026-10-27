# Interactive input, mode labels and parameter summary.

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

set_permutation_off_defaults() {
    PERMUTE_CHOICE="off"
    PERM_N="$DEFAULT_PERM_N"
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
    echo "    Kernel estimator     : exact grouped-state kernel"
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
        if [[ "${LOCAL_TE_VIRTUAL,,}" == "on" || "${LOCAL_TE_VIRTUAL,,}" == "1" || "${LOCAL_TE_VIRTUAL,,}" == "yes" ]]; then
            echo "    Local TE export      : virtual exact array, out=$LOCAL_TE_VIRTUAL_OUTPUT_DIR"
            echo "                           dtype=$LOCAL_TE_VIRTUAL_VALUE_DTYPE, state_code_bits=$LOCAL_TE_VIRTUAL_STATE_CODE_BITS, compression=$LOCAL_TE_VIRTUAL_COMPRESSION, buffer=$LOCAL_TE_VIRTUAL_BUFFER_EDGES"
            echo "                           code_index=$LOCAL_TE_VIRTUAL_CODE_INDEX, index_dir=$LOCAL_TE_CODE_INDEX_OUTPUT_DIR"
            echo "                           export materialized chunks later with code.localte_virtual_array export"
        elif [[ "${LOCAL_TE_LAZY,,}" == "on" || "${LOCAL_TE_LAZY,,}" == "1" || "${LOCAL_TE_LAZY,,}" == "yes" ]]; then
            echo "    Local TE export      : lazy exact recipe, out=$LOCAL_TE_RECIPE_DIR"
            echo "                           export selected edges later with code.localte_lazy_export"
        elif [[ "${LOCAL_TE_VIDEO_TILES,,}" == "on" ]]; then
            echo "    Local TE export      : video tiles, chunk=$LOCAL_CHUNK_SIZE, dtype=$LOCAL_VALUES_DTYPE, out=$LOCAL_TE_VIDEO_OUTPUT_DIR"
            echo "                           codec=${LOCAL_TE_VIDEO_CODEC}-${LOCAL_TE_VIDEO_LEVEL}, edges=$LOCAL_TE_VIDEO_EDGE_FORMAT, factorize=$LOCAL_TE_VIDEO_FACTORIZE"
            echo "                           direct_buffer=$LOCAL_TE_DIRECT_BUFFER_EDGES"
        else
            echo "    Local TE export      : parquet chunks, chunk=$LOCAL_CHUNK_SIZE, dtype=$LOCAL_VALUES_DTYPE, split=$LOCAL_SPLIT_EXPORT, split_dir=$LOCAL_SPLIT_OUTPUT_DIR"
            echo "                           direct_stream=$LOCAL_TE_DIRECT_CHUNKS, worker_parts=yes, fallback_export_workers=auto"
        fi
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
    if [[ "${TENET_REUSE_PREPARED,,}" != "off" && "${TENET_REUSE_PREPARED}" != "0" && "${TENET_REUSE_PREPARED,,}" != "no" ]]; then
        echo "    Prepared reuse       : ${TENET_REUSE_PREPARED}, dir=${TENET_PREPARED_DIR:-<unset>}, link=${TENET_PREPARED_LINK_MODE}"
    fi
    echo "12) Output directory     : $OUTPUT_DIR"
    echo "---------------------------"
}


prompt_output_directory() {
    print_step "12" "Output Directory" \
        "Use a separate directory for each run to keep previous results."
    local value
    value=$(prompt_with_default "Output directory" "$OUTPUT_DIR")
    if [[ "$value" != /* ]]; then value="$REPO_ROOT/$value"; fi
    OUTPUT_DIR=$(abs_path "$value")
    TENET_OUTPUT_DIR="$OUTPUT_DIR"
}
