#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"
PYTHON=${PYTHON:-python3}

prompt_with_default() {
    local label="$1"
    local default_value="$2"
    local input=""
    if [ -n "$default_value" ]; then
        read -e -r -i "$default_value" -p "$label [$default_value]> " input
    else
        read -e -r -p "$label> " input
    fi
    input="${input:-$default_value}"
    echo "${input#"${input%%[![:space:]]*}"}" | sed 's/[[:space:]]*$//'
}

show_summary() {
    echo ""
    echo "--- make_GRN_localTE configuration ---"
    echo " 1) Input file          : $INPUT_FILE"
    echo " 2) Output file         : ${OUTPUT_FILE:-<not set>}"
    echo " 3) History length (k)  : $HISTORY_LENGTH"
    echo " 4) FDR threshold       : $FDR_THRESHOLD"
    echo " 5) TE cutoff           : $TE_CUTOFF"
    echo " 6) Chunk mode          : $USE_CHUNK_MODE"
    if [ "${USE_CHUNK_MODE,,}" = "yes" ]; then
        echo "    • chunk_dir         : $CHUNK_DIR"
        echo "    • chunk_size        : ${CHUNK_SIZE:-<metadata>}"
        echo "    • chunk_workers     : $CHUNK_WORKERS"
    else
        echo "    • workers           : $WORKERS"
    fi
    echo " 7) time_labels         : ${TIME_LABELS:-<not set>}"
    echo " 8) time_map            : ${TIME_MAP:-<not set>}"
    echo " 9) max_timepoints      : ${MAX_TIMEPOINTS:-<not set>}"
    echo "--------------------------------------"
}

usage() {
    cat <<USAGE
Usage (recommended interactive): $0

Interactive mode (no arguments):
  - Prompts for:
      1) Input TE parquet path
      2) Output file path (optional; blank = skip write)
      3) History length k
      4) FDR threshold
      5) TE cutoff
      6) Chunk mode (yes/no) and parameters:
           - chunk_dir (default: output/local_te_split_chunks)
           - chunk_size (optional; auto from metadata when blank)
           - chunk_workers
      7) time_labels file (optional)
      8) time_map parquet (optional)
      9) max_timepoints (optional)

Non-interactive shorthand:
  $0 <input> <output> <history_k> <fdr> <te_cutoff>
    - Uses chunk_dir=output/local_te_split_chunks, chunk_size=300,
      chunk_workers = nproc (or 4 if unavailable).
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [ $# -eq 0 ]; then
    # Auto-detect a reasonable default input from recent Matrix_generate outputs.
    if [ -z "${INPUT_FILE:-}" ]; then
        if   [ -f "output/TE_TF_GN.parquet" ]; then
            INPUT_FILE="output/TE_TF_GN.parquet"
        elif [ -f "output/TE_GN_GN.parquet" ]; then
            INPUT_FILE="output/TE_GN_GN.parquet"
        elif [ -f "output/TE_all_features.parquet" ]; then
            INPUT_FILE="output/TE_all_features.parquet"
        else
            INPUT_FILE="output/TE_TF_GN.parquet"
        fi
    fi

    echo "Interactive mode for make_GRN_localTE"
    echo "Press Enter to accept defaults shown in brackets."

    OUTPUT_FILE="${OUTPUT_FILE:-output/Local_TE_TF_GN_fdr_0_01.parquet}"
    HISTORY_LENGTH="${HISTORY_LENGTH:-1}"
    FDR_THRESHOLD="${FDR_THRESHOLD:-0.01}"
    TE_CUTOFF="${TE_CUTOFF:-0.0}"
    USE_CHUNK_MODE="${USE_CHUNK_MODE:-yes}"
    CHUNK_DIR="${CHUNK_DIR:-output/local_te_split_chunks}"
    CHUNK_SIZE="${CHUNK_SIZE:-}"
    CHUNK_WORKERS="${CHUNK_WORKERS:-1}"
    WORKERS="${WORKERS:-0}"
    TIME_LABELS="${TIME_LABELS:-}"
    TIME_MAP="${TIME_MAP:-}"
    MAX_TIMEPOINTS="${MAX_TIMEPOINTS:-}"

    prompt_all() {
        INPUT_FILE="$(prompt_with_default "Input TE parquet path" "$INPUT_FILE")"
        OUTPUT_FILE="$(prompt_with_default "Output file path (leave blank to skip)" "$OUTPUT_FILE")"
        HISTORY_LENGTH="$(prompt_with_default "History length k" "$HISTORY_LENGTH")"
        FDR_THRESHOLD="$(prompt_with_default "FDR threshold" "$FDR_THRESHOLD")"
        TE_CUTOFF="$(prompt_with_default "TE cutoff" "$TE_CUTOFF")"
        USE_CHUNK_MODE="$(prompt_with_default "Use chunk_dir processing? (yes/no)" "$USE_CHUNK_MODE")"
        if [ "${USE_CHUNK_MODE,,}" = "yes" ]; then
            CHUNK_DIR="$(prompt_with_default "Chunk directory" "$CHUNK_DIR")"
            CHUNK_SIZE="$(prompt_with_default "Chunk size (timesteps per chunk, blank=auto)" "$CHUNK_SIZE")"
            CHUNK_WORKERS="$(prompt_with_default "Chunk workers (0=single)" "$CHUNK_WORKERS")"
        else
            WORKERS="$(prompt_with_default "Workers (0/1=single)" "$WORKERS")"
        fi
        TIME_LABELS="$(prompt_with_default "time_labels file (optional)" "$TIME_LABELS")"
        TIME_MAP="$(prompt_with_default "time_map parquet (optional)" "$TIME_MAP")"
        MAX_TIMEPOINTS="$(prompt_with_default "max_timepoints (optional)" "$MAX_TIMEPOINTS")"
    }

    prompt_all

    while true; do
        show_summary
        read -rp "Edit step number (1-9) or press Enter to run: " EDIT_CHOICE
        EDIT_CHOICE="$(echo "$EDIT_CHOICE" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        case "$EDIT_CHOICE" in
            "")
                break
                ;;
            1) INPUT_FILE="$(prompt_with_default "Input TE parquet path" "$INPUT_FILE")" ;;
            2) OUTPUT_FILE="$(prompt_with_default "Output file path (leave blank to skip)" "$OUTPUT_FILE")" ;;
            3) HISTORY_LENGTH="$(prompt_with_default "History length k" "$HISTORY_LENGTH")" ;;
            4) FDR_THRESHOLD="$(prompt_with_default "FDR threshold" "$FDR_THRESHOLD")" ;;
            5) TE_CUTOFF="$(prompt_with_default "TE cutoff" "$TE_CUTOFF")" ;;
            6)
                USE_CHUNK_MODE="$(prompt_with_default "Use chunk_dir processing? (yes/no)" "$USE_CHUNK_MODE")"
                if [ "${USE_CHUNK_MODE,,}" = "yes" ]; then
                    CHUNK_DIR="$(prompt_with_default "Chunk directory" "$CHUNK_DIR")"
                    CHUNK_SIZE="$(prompt_with_default "Chunk size (timesteps per chunk, blank=auto)" "$CHUNK_SIZE")"
                    CHUNK_WORKERS="$(prompt_with_default "Chunk workers (0=single)" "$CHUNK_WORKERS")"
                else
                    WORKERS="$(prompt_with_default "Workers (0/1=single)" "$WORKERS")"
                fi
                ;;
            7) TIME_LABELS="$(prompt_with_default "time_labels file (optional)" "$TIME_LABELS")" ;;
            8) TIME_MAP="$(prompt_with_default "time_map parquet (optional)" "$TIME_MAP")" ;;
            9) MAX_TIMEPOINTS="$(prompt_with_default "max_timepoints (optional)" "$MAX_TIMEPOINTS")" ;;
            *) echo "[!] Unknown selection." ;;
        esac
    done
else
    echo "Usage (interactive preferred): $0"
    echo "Or provide explicit arguments: $0 <input> <output> <history_k> <fdr> <te_cutoff>"
    if [ $# -lt 5 ]; then
        exit 1
    fi
    INPUT_FILE="$1"
    OUTPUT_FILE="$2"
    HISTORY_LENGTH="$3"
    FDR_THRESHOLD="$4"
    TE_CUTOFF="$5"
    USE_CHUNK_MODE="yes"
    CHUNK_DIR="${CHUNK_DIR:-output/local_te_split_chunks}"
    CHUNK_SIZE="${CHUNK_SIZE:-300}"
    CHUNK_WORKERS="${CHUNK_WORKERS:-$(command -v nproc >/dev/null && nproc || echo 4)}"
    TIME_LABELS="${TIME_LABELS:-}"
    TIME_MAP="${TIME_MAP:-}"
    MAX_TIMEPOINTS="${MAX_TIMEPOINTS:-}"
fi

CMD_ARGS=(
    --input "$INPUT_FILE"
    --history_length "$HISTORY_LENGTH"
    --fdr "$FDR_THRESHOLD"
    --te_cutoff "$TE_CUTOFF"
)

if [ -n "${OUTPUT_FILE:-}" ]; then
    CMD_ARGS+=(--output_file "$OUTPUT_FILE")
fi

if [ "${USE_CHUNK_MODE,,}" = "yes" ]; then
    CMD_ARGS+=(--chunk_dir "$CHUNK_DIR")
    if [ -n "${CHUNK_SIZE:-}" ]; then
        CMD_ARGS+=(--chunk_size "$CHUNK_SIZE")
    fi
    if [ -n "${CHUNK_WORKERS:-}" ]; then
        CMD_ARGS+=(--chunk_workers "$CHUNK_WORKERS")
    fi
else
    if [ -n "${WORKERS:-}" ]; then
        CMD_ARGS+=(--workers "$WORKERS")
    fi
fi

if [ -n "${TIME_LABELS:-}" ]; then
    CMD_ARGS+=(--time_labels "$TIME_LABELS")
fi
if [ -n "${TIME_MAP:-}" ]; then
    CMD_ARGS+=(--time_map "$TIME_MAP")
fi
if [ -n "${MAX_TIMEPOINTS:-}" ]; then
    CMD_ARGS+=(--max_timepoints "$MAX_TIMEPOINTS")
fi

# Apply chunk parameter overrides for the exporter to pick up from environment
echo ""
echo "Running: $PYTHON -m code.make_GRN_localTE ${CMD_ARGS[*]}"
"$PYTHON" -m code.make_GRN_localTE "${CMD_ARGS[@]}"
