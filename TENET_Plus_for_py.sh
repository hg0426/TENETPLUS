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
case "${1:-}" in
    --config|--check-config|--write-config)
        exec "$PYTHON" "$CODE_DIR/run_config.py" "$@" ;;
esac


source "$CODE_DIR/runner/defaults.sh"
source "$CODE_DIR/runner/interactive.sh"
source "$CODE_DIR/runner/replay.sh"

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

link_or_copy_prepared_file() {
  local src="$1"
  local dst="$2"
  if [[ ! -f "$src" ]]; then
    echo "[Prepared] Missing required prepared file: $src" >&2
    return 1
  fi
  if [[ "$src" -ef "$dst" ]]; then
    return 0
  fi
  rm -f "$dst"
  if [[ "${TENET_PREPARED_LINK_MODE,,}" == "copy" ]]; then
    cp -f "$src" "$dst"
  else
    ln -s "$(abs_path "$src")" "$dst"
  fi
}

validate_prepared_manifest() {
  local prepared_dir="$1"
  local manifest="$prepared_dir/prepared_input_manifest.json"
  local reuse_mode="${TENET_REUSE_PREPARED,,}"
  if [[ "$reuse_mode" == "force" ]]; then
    return 0
  fi
  if [[ ! -f "$manifest" ]]; then
    echo "[Prepared] No prepared_input_manifest.json found. Use TENET_REUSE_PREPARED=force to trust this directory." >&2
    return 1
  fi
  "$PYTHON" - "$manifest" "$RUN_INPUT_FILE" "$RUN_TRAJECTORY_FILE" "$RUN_CELL_SELECT_FILE" "$HISTORY_LENGTH" "$RUN_SPECIES" "$RUN_NETWORK_MODE" "${TENET_PAIR_MODE:-default}" <<'PY'
import json, os, sys
manifest_path, input_matrix, trajectory_file, cell_select_file, history, species, mode, pair_mode = sys.argv[1:]
manifest = json.load(open(manifest_path, encoding="utf-8"))
expected = {
    "input_matrix": os.path.abspath(os.path.expanduser(input_matrix)),
    "trajectory_file": os.path.abspath(os.path.expanduser(trajectory_file)),
    "cell_select_file": os.path.abspath(os.path.expanduser(cell_select_file)),
    "history_length": history,
    "species": species,
    "network_mode": mode,
    "pair_mode": pair_mode,
}
errors = []
for key, expected_value in expected.items():
    observed = str(manifest.get(key, ""))
    if observed != str(expected_value):
        errors.append(f"{key}: expected {expected_value!r}, observed {observed!r}")
if errors:
    print("[Prepared] Manifest mismatch:")
    for item in errors:
        print("  - " + item)
    sys.exit(1)
PY
}

try_reuse_prepared_inputs() {
  PREPARED_REUSED=off
  local reuse_mode="${TENET_REUSE_PREPARED,,}"
  if [[ "$reuse_mode" == "off" || "$reuse_mode" == "0" || "$reuse_mode" == "no" || -z "$reuse_mode" ]]; then
    return 1
  fi
  if [[ -z "${TENET_PREPARED_DIR:-}" ]]; then
    echo "[Prepared] TENET_REUSE_PREPARED=${TENET_REUSE_PREPARED}, but TENET_PREPARED_DIR is empty." >&2
    return 1
  fi
  local prepared_dir
  if [[ "$TENET_PREPARED_DIR" == /* ]]; then
    prepared_dir="$(abs_path "$TENET_PREPARED_DIR")"
  else
    prepared_dir="$(abs_path "$REPO_ROOT/$TENET_PREPARED_DIR")"
  fi
  if [[ ! -d "$prepared_dir" ]]; then
    echo "[Prepared] Prepared directory does not exist: $prepared_dir" >&2
    return 1
  fi
  validate_prepared_manifest "$prepared_dir" || return 1
  if [[ ! -f "$prepared_dir/feature_manifest.json" ]]; then
    echo "[Prepared] Missing reference-aware feature manifest; regenerate prepared inputs." >&2
    return 1
  fi
  link_or_copy_prepared_file "$prepared_dir/feature_manifest.json" "feature_manifest.json" || return 1

  set_stage "PREPARED_REUSE"
  echo "[Prepared] Reusing prepared inputs from: $prepared_dir"
  echo "[Prepared] Link mode: ${TENET_PREPARED_LINK_MODE}"

  link_or_copy_prepared_file "$prepared_dir/cell_gene_trsps.parquet" "cell_gene_trsps.parquet" || return 1
  link_or_copy_prepared_file "$prepared_dir/gene_names" "gene_names" || return 1
  link_or_copy_prepared_file "$prepared_dir/TE_peak_list.txt" "TE_peak_list.txt" || return 1
  if [[ "${TENET_PAIR_MODE:-default}" == "default" ]]; then
    link_or_copy_prepared_file "$prepared_dir/all_pairs.csv" "all_pairs.csv" || return 1
  fi
  if [[ -f "$prepared_dir/time_index_map.parquet" ]]; then
    link_or_copy_prepared_file "$prepared_dir/time_index_map.parquet" "time_index_map.parquet" || true
  fi
  matrix="cell_gene_trsps.parquet"
  trajectory="$RUN_TRAJECTORY_FILE"
  cell_select="$RUN_CELL_SELECT_FILE"
  PREPARED_REUSED=on
  echo "[Prepared] Reuse ready: matrix processing, transpose, and pair generation will be skipped."
  return 0
}

write_prepared_manifest() {
  "$PYTHON" - "$RUN_INPUT_FILE" "$RUN_TRAJECTORY_FILE" "$RUN_CELL_SELECT_FILE" "$HISTORY_LENGTH" "$RUN_SPECIES" "$RUN_NETWORK_MODE" "${TENET_PAIR_MODE:-default}" <<'PY'
import json, os, time, sys
input_matrix, trajectory_file, cell_select_file, history, species, mode, pair_mode = sys.argv[1:]
payload = {
    "format": "tenet_prepared_input_manifest_v1",
    "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    "input_matrix": os.path.abspath(os.path.expanduser(input_matrix)),
    "trajectory_file": os.path.abspath(os.path.expanduser(trajectory_file)),
    "cell_select_file": os.path.abspath(os.path.expanduser(cell_select_file)),
    "history_length": history,
    "species": species,
    "network_mode": mode,
    "pair_mode": pair_mode,
    "prepared_files": {
        "cell_gene_trsps": os.path.abspath("cell_gene_trsps.parquet") if os.path.exists("cell_gene_trsps.parquet") else None,
        "gene_names": os.path.abspath("gene_names") if os.path.exists("gene_names") else None,
        "te_peak_list": os.path.abspath("TE_peak_list.txt") if os.path.exists("TE_peak_list.txt") else None,
        "all_pairs": os.path.abspath("all_pairs.csv") if os.path.exists("all_pairs.csv") else None,
    },
}
with open("prepared_input_manifest.json", "w", encoding="utf-8") as fh:
    json.dump(payload, fh, indent=2)
PY
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

# Function to display usage
usage() {
    cat <<USAGE
Usage:
  $0
  $0 --write-config analysis.ini
  $0 --check-config analysis.ini
  $0 --config analysis.ini [--output NEW_DIRECTORY]
  $0 <input_file> <num_jobs> <trajectory_file> <cell_select_file> <history_k> <species> <network_mode>

Use --config for a non-interactive run with named settings in an INI file.
Every run saves run_config.ini and run_invocation.sh in its output directory.
Config paths are relative to the config file; --output is relative to the current
directory. The seven-argument interface and environment settings remain supported.

Species resources:
  Presets: human, mouse, zebrafish, fly, worm, rat.
  TENET_REFERENCE=/path/to/bundle overrides the species reference.
  The verified bundled human reference defaults to human_legacy cis mapping.
  Other references default to legacy (skip unmapped RNA; last row for multi-contig RNA).
  TENET_CIS_MAPPING_POLICY=strict requires unique RNA contigs for cis modes.
  TENET_TF_LIST=/path/to/TFs.txt supplies matrix-matching TF IDs (modes 0-4).
  Only human/mouse TF lists are bundled. Modes 5/6 and all-pair modes need no TF list.
  Run python -m code.species_support --species SPECIES --mode MODE --matrix MATRIX
  before a large run. See README.md and docs/OTHER_SPECIES.md.

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
    Name aliases and a0/a1/a2 are accepted.

Optional:
  This package runs kernel TE on prepared inputs. Optional behavior is
  controlled with the environment variables below.

Environment:
  TENET_PAIR_MODE     : default | gene_only | all_feature.
                        Optional pair-scope override for non-interactive runs.
  TENET_OUTPUT_DIR    : output directory. Relative paths are resolved from
                        the package root (default: output).
  TENET_INTERMEDIATE_SAVE: on | off (default: off).
  TENET_PERMUTE       : on | off (default: off).
  TENET_PERM_N        : number of permutations (default: 100).
  TENET_PERM_FDR      : on | off. Use permutation q-value threshold
                        instead of p-value threshold (default: off).
  TENET_PERM_Q_ALPHA  : permutation q-value threshold (default: 0.05).
  TENET_PERM_ALPHA    : permutation p-value threshold (default: 0.01).
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
  TENET_STORE_LOCAL_TE: on | off (default: off).
  LOCAL_TE_VIRTUAL   : on | off (default: off). When LocalTE is on, store
                       exact compact state tables for all edge x timepoint
                       values instead of materializing dense LocalTE vectors.
  TENET_RESULTS_BUFFER_ROWS:
                        rows to buffer in TE before flushing batch parquet
                        files (default auto: 200000, or 5000 with LocalTE).
  TENET_REUSE_PREPARED:
                        off | on | force (default: off). Reuse a prepared
                        output directory and skip CSV filtering/transpose/pair
                        generation. 'on' requires prepared_input_manifest.json
                        to match current inputs; 'force' still requires valid feature/reference metadata.
  TENET_PREPARED_DIR  : output directory containing cell_gene_trsps.parquet,
                        gene_names, TE_peak_list.txt, and all_pairs.csv.
  TENET_PREPARED_LINK_MODE:
                        symlink | copy (default: symlink).
  LOCAL_TE_*          : LocalTE export tuning; most users can leave these unset.
USAGE
    return 0
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
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
    prompt_permutation_toggle
    prompt_local_te_choice
    prompt_grn_choice
    prompt_output_directory
    INTERMEDIATE_SAVE_CHOICE="${INTERMEDIATE_SAVE_CHOICE:-$DEFAULT_INTERMEDIATE_SAVE}"
    RESULTS_BUFFER_ROWS="${RESULTS_BUFFER_ROWS:-$DEFAULT_RESULTS_BUFFER_ROWS}"
    BATCH_SIZE="${BATCH_SIZE:-$DEFAULT_BATCH_SIZE}"

    while true; do
        show_summary
        read -rp "Edit step (1-12) or press Enter to continue: " EDIT_CHOICE
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
            12)
                prompt_output_directory
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
        "$HISTORY_LENGTH" "$SPECIES" "$MODE_CODE"
    export TENET_PERMUTE="$PERMUTE_CHOICE"
    export TENET_PERM_N="$PERM_N"
    export TENET_PERM_FDR="$PERM_FDR"
    export TENET_PERM_Q_ALPHA="$PERM_Q_ALPHA"
    export TENET_PERM_ALPHA="$PERM_ALPHA"
    export TENET_STORE_LOCAL_TE="$STORE_LOCAL_TE_CHOICE"
    if [ -n "${RESULTS_BUFFER_ROWS:-}" ]; then
        export TENET_RESULTS_BUFFER_ROWS="$RESULTS_BUFFER_ROWS"
    else
        unset TENET_RESULTS_BUFFER_ROWS
    fi

    mkdir -p "$OUTPUT_DIR"

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

# Check the compact seven-argument non-interactive interface.
if [ $# -lt 7 ]; then
    echo "Insufficient arguments supplied."
    usage
    exit 2
fi
if [ $# -gt 7 ]; then
    echo "Too many positional arguments supplied."
    echo "Use the seven core arguments plus TENET_* / LOCAL_TE_* environment variables for optional settings."
    usage
    exit 2
fi

ARG1_TRIM=$(trim_whitespace "$1")
ARG2_TRIM=$(trim_whitespace "$2")
ARG3_TRIM=$(trim_whitespace "$3")
ARG4_TRIM=$(trim_whitespace "$4")
ABS_INPUT=$(abs_path "$ARG1_TRIM")
ABS_TRAJECTORY=$(abs_path "$ARG3_TRIM")
ABS_CELL_SELECT=$(abs_path "$ARG4_TRIM")
set -- "$ABS_INPUT" "$ARG2_TRIM" "$ABS_TRAJECTORY" "$ABS_CELL_SELECT" "${@:5}"
RUN_INPUT_FILE="$1"
RUN_TRAJECTORY_FILE="$3"
RUN_CELL_SELECT_FILE="$4"
RUN_SPECIES="$6"
RUN_NETWORK_MODE="$7"
# Default HISTORY_LENGTH from 5th argument when not set via interactive mode
HISTORY_LENGTH="${HISTORY_LENGTH:-$5}"

# Accept the compact analysis-mode aliases used by the interactive prompt.
ANALYSIS_ARG="${7,,}"
case "$ANALYSIS_ARG" in
  rna_tf|tf_rna|tf|rna|tenet_tf)
    PAIR_MODE="default"
    MODE_CODE="0"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  plus_full|full|tenet_plus_full)
    PAIR_MODE="default"
    MODE_CODE="1"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  tf_gene|tf_gn|tf-gene|gene_grn|grn)
    PAIR_MODE="default"
    MODE_CODE="2"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  tf_peak|tf_pk|tf-peak)
    PAIR_MODE="default"
    MODE_CODE="3"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  tf_targets|tf_gene_peak|tf_gn_pk|tf-gene-peak|tf_gene+peak)
    PAIR_MODE="default"
    MODE_CODE="4"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  peak_gene|peak_gn|peak-gene|cis_peak_gene)
    PAIR_MODE="default"
    MODE_CODE="5"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  peak_peak|peak_pk|peak-peak|cis_peak_peak)
    PAIR_MODE="default"
    MODE_CODE="6"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  7|a0|all_gene_rna|rna_all_gene)
    PAIR_MODE="all_pair"
    MODE_CODE="0"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  8|a1|all_gene|all_gene_plus|gene_all_gene)
    PAIR_MODE="all_pair"
    MODE_CODE="1"
    set -- "${@:1:6}" "$MODE_CODE"
    ;;
  9|a2|all_feature|all_features|feature_all_feature)
    PAIR_MODE="all_pair"
    MODE_CODE="2"
    set -- "${@:1:6}" "$MODE_CODE"
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
export TENET_PAIR_MODE

INPUT_FILE="$1"
NUM_JOBS="$2"
TRAJECTORY_FILE="$3"
CELL_SELECT_FILE="$4"
HISTORY_LENGTH="$5"
SPECIES="$6"
MODE_CODE="${MODE_CODE:-$7}"
update_defaults_after_mode

PERMUTE_CHOICE="${TENET_PERMUTE:-$DEFAULT_PERMUTE}"
PERM_N="${TENET_PERM_N:-$DEFAULT_PERM_N}"
PERM_FDR="${TENET_PERM_FDR:-$DEFAULT_PERM_FDR}"
PERM_Q_ALPHA="${TENET_PERM_Q_ALPHA:-$DEFAULT_PERM_Q_ALPHA}"
PERM_ALPHA="${TENET_PERM_ALPHA:-$DEFAULT_PERM_ALPHA}"
PERM_SCOPE="$DEFAULT_PERM_SCOPE"
if [[ -n "${TENET_PERM_CANDIDATE_GRN_FDR:-}" && "${TENET_PERM_CANDIDATE_GRN_FDR:-0}" != "0" ]]; then
    PERM_SCOPE="grn_fdr"
    PERM_GRN_FDR="${TENET_PERM_CANDIDATE_GRN_FDR}"
else
    PERM_SCOPE="all_pairs"
    PERM_GRN_FDR="$DEFAULT_PERM_GRN_FDR"
fi
PERM_GRN_TE_CUTOFF="${TENET_PERM_CANDIDATE_TE_CUTOFF:-$DEFAULT_PERM_GRN_TE_CUTOFF}"

STORE_LOCAL_TE_CHOICE="${TENET_STORE_LOCAL_TE:-$DEFAULT_LOCAL_TE}"
LOCAL_CHUNK_SIZE="${LOCAL_TE_CHUNK_SIZE}"
LOCAL_BUFFER_EDGES="${LOCAL_TE_BUFFER_EDGES}"
LOCAL_EXPORT_WORKERS="${LOCAL_TE_EXPORT_WORKERS}"
LOCAL_MERGE_WORKERS="${LOCAL_TE_MERGE_WORKERS:-$DEFAULT_LOCAL_MERGE_WORKERS}"
LOCAL_READ_BATCH="${LOCAL_TE_READ_BATCH_ROWS}"
LOCAL_THREADS_FLAG="${LOCAL_TE_USE_THREADS}"
LOCAL_VALUES_DTYPE="${LOCAL_TE_VALUES_DTYPE}"
LOCAL_SPLIT_EXPORT="${LOCAL_TE_SPLIT_EXPORT}"
LOCAL_SPLIT_OUTPUT_DIR="${LOCAL_TE_SPLIT_OUTPUT_DIR}"

MAKE_GRN_CHOICE="${TENET_MAKE_GRN:-$DEFAULT_MAKE_GRN}"
GRN_FDR="${TENET_GRN_FDR:-$DEFAULT_GRN_FDR}"
TRIM_INDIRECT_CHOICE="${TENET_TRIM_INDIRECT:-$DEFAULT_TRIM_INDIRECT}"
TRIM_THRESHOLD="${TENET_TRIM_THRESHOLD:-$DEFAULT_TRIM_THRESHOLD}"
INTERMEDIATE_SAVE_CHOICE="${TENET_INTERMEDIATE_SAVE:-$DEFAULT_INTERMEDIATE_SAVE}"
RESULTS_BUFFER_ROWS="${TENET_RESULTS_BUFFER_ROWS:-}"
BATCH_SIZE="${TE_BATCH_SIZE:-$DEFAULT_BATCH_SIZE}"

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

export TENET_OUTPUT_DIR="$OUTPUT_DIR"
# The Python payload writer must receive the dtype selected in the prompt too.
export LOCAL_TE_VALUES_DTYPE
show_summary "$@"
write_replay_logs "$@"

# Fail before matrix conversion, monitoring or TE if species resources are absent.
if [[ -n "${TENET_TF_LIST:-}" ]]; then
  export TENET_TF_LIST="$(abs_path "$TENET_TF_LIST")"
fi
"$PYTHON" -m code.species_support --species "$6" --mode "$7" \
  --matrix "$1" --output "$OUTPUT_DIR/species_preflight.json"

chmod +x "$CODE_DIR/memory_check.sh" || true
CURRENT_STAGE_FILE="$OUTPUT_DIR/current_stage.txt"
echo "INIT" > "$CURRENT_STAGE_FILE"
export TENET_STAGE_FILE="$CURRENT_STAGE_FILE"
"$CODE_DIR/memory_check.sh" $$ "$CURRENT_STAGE_FILE" &
MONITOR_PID=$!
cleanup_memory_monitor() {
    local exit_status=$?
    trap - EXIT
    kill "$MONITOR_PID" 2>/dev/null || true
    wait "$MONITOR_PID" 2>/dev/null || true
    exit "$exit_status"
}
trap cleanup_memory_monitor EXIT
echo "[Run] RAM monitor started: pid=${MONITOR_PID}"

if [[ -n "${TENET_REFERENCE:-}" && -d "$TENET_REFERENCE" ]]; then
  export TENET_REFERENCE="$(abs_path "$TENET_REFERENCE")"
fi
if [[ -n "${TENET_FEATURE_TYPES:-}" ]]; then
  export TENET_FEATURE_TYPES="$(abs_path "$TENET_FEATURE_TYPES")"
fi
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
    "$PYTHON" - "$INPUT_FILE" "$TEMP_CSV" <<'END'
import pandas as pd
import sys

try:
    df = pd.read_parquet(sys.argv[1], columns=None)
    df.head(0).to_csv(sys.argv[2], index=False)
except Exception as e:
    print(f"Error converting Parquet to CSV: {e}")
    sys.exit(1)
END
}

PREPARED_REUSED=off
if ! try_reuse_prepared_inputs; then
  case "${TENET_REUSE_PREPARED,,}" in
    on|1|yes|true|force)
      echo "Error: prepared input reuse was requested but could not be activated. Exiting instead of falling back to slow preprocessing." >&2
      exit 1
      ;;
  esac

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

# Classify features using the shared candidate parser.
"$PYTHON" -m code.feature_utils "$INPUT_FILE" "$ARG7" "$6"

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
    # Keep generated inputs isolated within this run's output directory.
    matrix="$OUTPUT_DIR/filtered_matrix.parquet"
    trajectory="$OUTPUT_DIR/filtered_trajectory.txt"
    cell_select="$OUTPUT_DIR/filtered_cellselect.txt"
else
    set_stage "MATRIX_PROCESSING_SKIP"
    echo "[Matrix] Input already appears selected/sorted; using it as-is."
    matrix=$1
    trajectory=$3
    cell_select=$4
fi

# Always regenerate gene_names from the current matrix
set_stage "GENE_NAMES_REGEN"
echo "[Setup] Writing gene_names from matrix header."
if ! "$PYTHON" - "$matrix" "$OUTPUT_DIR/gene_names" <<'PY'
import sys
import pandas as pd
import pathlib
path = pathlib.Path(sys.argv[1])
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
    for dst in [sys.argv[2]]:
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

# If a precomputed transposed matrix exists, remove it to regenerate fresh
if [ -f cell_gene_trsps.parquet ]; then
  echo "[Matrix] Removing stale cell_gene_trsps.parquet before regeneration."
  rm -f cell_gene_trsps.parquet
fi

"$PYTHON" -m code.feature_utils "${matrix}" "$ARG7" "$6"

set_stage "PREPARED_MATRIX"
echo "[Setup] Prepared matrix accepted."

if [ "$7" -ne 0 ]; then
  echo "[Pairs] Building selected source-target pairs for mode ${7}."
    if [ "${TENET_PAIR_MODE:-default}" = "default" ]; then
      if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "PREPROCESS_TENET_PLUS" "$PYTHON" -m code.PreProcessScript_TE_Plus "$6" "$7" "${matrix}"; then
          echo "Error in PreProcessScript_TE_Plus.py. Exiting."
          exit 1
      fi
    else
      echo "[Pairs] Using implicit all-pair mode in TE core (${TENET_PAIR_MODE}); no all_pairs.csv materialization."
    fi
else
  echo "[Pairs] Building RNA-only TF source-target pairs."
  if [ "${TENET_PAIR_MODE:-default}" = "default" ]; then
    if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "PREPROCESS_TENET_TF" "$PYTHON" -m code.PreProcessScript_TENET_TF "$6"; then
        echo "Error in PreProcessScript_TENET_TF.py. Exiting."
        exit 1
    fi
  else
    echo "[Pairs] Using implicit all-pair mode in TE core (${TENET_PAIR_MODE}); no all_pairs.csv materialization."
  fi
fi

# Ensure transposed matrix exists for TE calculation.
if [ ! -f cell_gene_trsps.parquet ]; then
  echo "--- Generating cell_gene_trsps.parquet by transposing matrix ---"
  if ! "$PYTHON" - "$matrix" <<'PY'
import sys, pathlib
import pandas as pd
path = pathlib.Path(sys.argv[1])
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

fi

# Validate copied manifests before any TE calculation (including reuse).
"$PYTHON" -m code.feature_utils --validate-prepared "$ARG7" "$6"

# A valid feature manifest may still yield no edges in a selected network mode.
if [[ "${TENET_PAIR_MODE:-default}" == "default" && ! -s all_pairs.csv ]]; then
  echo "Error: No candidate pairs for network mode ${ARG7}; TE calculation was not started." >&2
  echo "Check the selected features, TF identifiers, and cis mappings. Mode 6 needs at least two peaks on the same contig." >&2
  exit 2
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

write_prepared_manifest

# TE computation
echo "[TE] Kernel estimator: exact grouped-state fast path enabled."

# Permutation options (optional)
PERMUTE="${TENET_PERMUTE:-off}"
PERM_N="${TENET_PERM_N:-100}"
PERM_FDR="${TENET_PERM_FDR:-off}"
PERM_Q_ALPHA="${TENET_PERM_Q_ALPHA:-0.05}"
PERM_ALPHA="${TENET_PERM_ALPHA:-0.01}"
STORE_LOCAL_TE="${TENET_STORE_LOCAL_TE:-off}"
LOCAL_TE_LAZY_ON=off
if [[ "${LOCAL_TE_LAZY,,}" == "on" || "${LOCAL_TE_LAZY,,}" == "1" || "${LOCAL_TE_LAZY,,}" == "yes" ]]; then
  LOCAL_TE_LAZY_ON=on
fi
LOCAL_TE_VIRTUAL_ON=off
if [[ "${LOCAL_TE_VIRTUAL,,}" == "on" || "${LOCAL_TE_VIRTUAL,,}" == "1" || "${LOCAL_TE_VIRTUAL,,}" == "yes" ]]; then
  LOCAL_TE_VIRTUAL_ON=on
fi
LOCAL_ARGS=()
if [[ "$STORE_LOCAL_TE" == "on" || "$STORE_LOCAL_TE" == "ON" || "$STORE_LOCAL_TE" == "1" ]]; then
  if [[ "$LOCAL_TE_VIRTUAL_ON" == "on" ]]; then
    rm -f "${LOCAL_TE_VIRTUAL_OUTPUT_DIR}/direct_virtual_localte_done.json" 2>/dev/null || true
    echo "[TE] LocalTE storage: VIRTUAL exact array. TE scores are kept separately; timepoint values are reconstructed from compact state tables."
  elif [[ "$LOCAL_TE_LAZY_ON" == "on" ]]; then
    echo "[TE] LocalTE storage: LAZY exact recipe. Full LocalTE arrays are not materialized during TE."
  else
    LOCAL_ARGS=(--store_local_te --localte_codec "$LOCAL_TE_CODEC")
    rm -f "${LOCAL_TE_SPLIT_OUTPUT_DIR}/direct_localte_done.json" 2>/dev/null || true
    rm -f "${LOCAL_TE_VIDEO_OUTPUT_DIR}/direct_localte_done.json" 2>/dev/null || true
    echo "[TE] LocalTE storage: ON (codec=${LOCAL_TE_CODEC}). This adds per-pair timepoint arrays and can be I/O/compression bound."
  fi
else
  echo "[TE] LocalTE storage: OFF."
fi

TIME_ARGS=()

# TE results buffer rows (optional; controls when intermediate batches flush to disk)
RESULTS_BUFFER_ROWS="${TENET_RESULTS_BUFFER_ROWS:-}"
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
  echo "[Permutation] ON: n=${PERM_N}, fdr=${PERM_FDR}"
  PERM_ARGS=(--permute --perm_n "$PERM_N")
  PERM_CANDIDATE_GRN_FDR="${TENET_PERM_CANDIDATE_GRN_FDR:-}"
  if [[ -n "${PERM_CANDIDATE_GRN_FDR:-}" && "${PERM_CANDIDATE_GRN_FDR:-0}" != "0" ]]; then
    PERM_ARGS+=(--perm_candidate_grn_fdr "$PERM_CANDIDATE_GRN_FDR")
    PERM_ARGS+=(--perm_candidate_te_cutoff "${TENET_PERM_CANDIDATE_TE_CUTOFF:-0.0}")
    echo "--- Permutation candidate filter: grn_fdr=${PERM_CANDIDATE_GRN_FDR}, te_cutoff=${TENET_PERM_CANDIDATE_TE_CUTOFF:-0.0} ---"
  fi
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
      --species "$6" \
      --split_type "$7" \
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
      --species "$6" \
      --split_type "$7" \
      "${BUFFER_ARGS[@]}" \
      "${LOCAL_ARGS[@]}" \
      "${TIME_ARGS[@]}"; then
    echo "Error in kernel TE stage. Exiting."
    exit 1
  fi
fi

rm -f TE_TF_GN.parquet TE_TF_PK.parquet TE_PK_GN.parquet TE_PK_PK.parquet TE_GN_GN.parquet TE_all_features.parquet
if ! TENET_PAIR_MODE="$TENET_PAIR_MODE" run_with_stage "MATRIX_GENERATE" "$PYTHON" -m code.Matrix_generate TE_result_all.parquet "$6" "$7"; then
    echo "Error in Matrix_generate.py. Exiting."
    exit 1
fi

build_grn_networks

if [[ "$STORE_LOCAL_TE" == "on" || "$STORE_LOCAL_TE" == "ON" || "$STORE_LOCAL_TE" == "1" ]]; then
    if [[ "$LOCAL_TE_VIRTUAL_ON" == "on" ]]; then
        mkdir -p "$LOCAL_TE_RECIPE_DIR"
        "$PYTHON" - "$LOCAL_TE_RECIPE_DIR" "$LOCAL_TE_VIRTUAL_OUTPUT_DIR" "$LOCAL_TE_CODE_INDEX_OUTPUT_DIR" "$HISTORY_LENGTH" "$SPECIES" "$7" "$TENET_PAIR_MODE" <<'PY'
import json, os, time, sys
recipe_path, virtual_path, code_index_path, history, species, mode, pair_mode = sys.argv[1:]
recipe_dir = os.path.abspath(recipe_path)
recipe = {
    "format": "tenet_localte_lazy_recipe_v1",
    "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    "work_dir": os.getcwd(),
    "matrix_parquet": os.path.abspath("cell_gene_trsps.parquet"),
    "gene_names": os.path.abspath("gene_names"),
    "time_index_map": os.path.abspath("time_index_map.parquet") if os.path.exists("time_index_map.parquet") else None,
    "all_pairs": os.path.abspath("all_pairs.csv") if os.path.exists("all_pairs.csv") else None,
    "te_result_all": os.path.abspath("TE_result_all.parquet"),
    "virtual_array": os.path.abspath(virtual_path),
    "code_index": os.path.abspath(code_index_path) if os.path.exists(os.path.join(code_index_path, "metadata.json")) else None,
    "history_length": int(history),
    "kernel_width": 0.5,
    "species": species,
    "network_mode": mode,
    "pair_mode": pair_mode,
    "note": "Exact LocalTE timepoint values are stored as compact virtual state tables. Use python -m code.localte_virtual_array export to materialize chunks.",
}
with open(os.path.join(recipe_dir, "recipe.json"), "w", encoding="utf-8") as fh:
    json.dump(recipe, fh, indent=2)
PY
        echo "[LocalTE] Virtual recipe written: ${LOCAL_TE_RECIPE_DIR}/recipe.json"
        if [[ -f "${LOCAL_TE_VIRTUAL_OUTPUT_DIR}/direct_virtual_localte_done.json" ]]; then
            echo "[LocalTE] Virtual array ready: ${LOCAL_TE_VIRTUAL_OUTPUT_DIR}"
        else
            echo "Error: virtual LocalTE marker not found at ${LOCAL_TE_VIRTUAL_OUTPUT_DIR}/direct_virtual_localte_done.json"
            exit 1
        fi
        set_stage "LOCALTE_VIRTUAL_READY"
    elif [[ "$LOCAL_TE_LAZY_ON" == "on" ]]; then
        mkdir -p "$LOCAL_TE_RECIPE_DIR"
        "$PYTHON" - "$LOCAL_TE_RECIPE_DIR" "$LOCAL_TE_VIRTUAL_OUTPUT_DIR" "$LOCAL_TE_CODE_INDEX_OUTPUT_DIR" "$HISTORY_LENGTH" "$SPECIES" "$7" "$TENET_PAIR_MODE" <<'PY'
import json, os, time, sys
recipe_path, virtual_path, code_index_path, history, species, mode, pair_mode = sys.argv[1:]
recipe_dir = os.path.abspath(recipe_path)
recipe = {
    "format": "tenet_localte_lazy_recipe_v1",
    "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    "work_dir": os.getcwd(),
    "matrix_parquet": os.path.abspath("cell_gene_trsps.parquet"),
    "gene_names": os.path.abspath("gene_names"),
    "time_index_map": os.path.abspath("time_index_map.parquet") if os.path.exists("time_index_map.parquet") else None,
    "all_pairs": os.path.abspath("all_pairs.csv") if os.path.exists("all_pairs.csv") else None,
    "te_result_all": os.path.abspath("TE_result_all.parquet"),
    "history_length": int(history),
    "kernel_width": 0.5,
    "species": species,
    "network_mode": mode,
    "pair_mode": pair_mode,
    "note": "Exact LocalTE arrays were not materialized during TE. Use python -m code.localte_lazy_export to export selected edges on demand.",
}
with open(os.path.join(recipe_dir, "recipe.json"), "w", encoding="utf-8") as fh:
    json.dump(recipe, fh, indent=2)
PY
        echo "[LocalTE] Lazy recipe written: ${LOCAL_TE_RECIPE_DIR}/recipe.json"
        echo "[LocalTE] Full LocalTE export skipped. Export selected edges later with code.localte_lazy_export."
        set_stage "LOCALTE_LAZY_READY"
    else
    DIRECT_LOCALTE_MARKER="${LOCAL_TE_SPLIT_OUTPUT_DIR}/direct_localte_done.json"
    if [[ "${LOCAL_TE_VIDEO_TILES,,}" == "on" ]]; then
        DIRECT_LOCALTE_MARKER="${LOCAL_TE_VIDEO_OUTPUT_DIR}/direct_localte_done.json"
    fi
    if [[ "${LOCAL_TE_DIRECT_CHUNKS,,}" == "on" && "${LOCAL_TE_SPLIT_EXPORT,,}" == "on" && -f "$DIRECT_LOCALTE_MARKER" ]]; then
        echo "[LocalTE] Direct LocalTE streaming already completed; skipping secondary LocalTE export."
        echo "[LocalTE] Marker: $DIRECT_LOCALTE_MARKER"
        set_stage "LOCALTE_DIRECT_DONE"
    else
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
    if [[ ! -f "$LOCAL_TE_INPUT_PARQUET" ]]; then
        echo "Error: ${LOCAL_TE_INPUT_PARQUET} not found for LocalTE export."
        exit 1
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
    fi
fi
