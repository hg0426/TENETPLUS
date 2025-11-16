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
    echo "--- make_GRN_averageTE configuration ---"
    echo " 1) Input file     : $INPUT_FILE"
    echo " 2) Output file    : $OUTPUT_FILE"
    echo " 3) FDR threshold  : $FDR_THRESHOLD"
    echo " 4) TE cutoff      : $TE_CUTOFF"
    echo "----------------------------------------"
}

usage() {
    cat <<USAGE
Usage: $0 <input_file> <output_file> <fdr_threshold> <te_cutoff>

Interactive mode (no arguments):
  - Prompts for:
      1) Input TE table path
      2) Output file (default: output/average_GRN.fdr0.01.parquet)
      3) FDR threshold (default: 0.01)
      4) TE cutoff (default: 0.0)

Non-interactive mode:
  - Provide all four arguments explicitly as shown above.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [ $# -eq 0 ]; then
    # Prefer mode-aware defaults from the newest Matrix_generate outputs.
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

    echo "Interactive mode for make_GRN_averageTE"
    echo "Press Enter to accept defaults shown in brackets."

    OUTPUT_FILE="${OUTPUT_FILE:-output/average_GRN.fdr0.01.parquet}"
    FDR_THRESHOLD="${FDR_THRESHOLD:-0.01}"
    TE_CUTOFF="${TE_CUTOFF:-0.0}"

    prompt_sequence() {
        INPUT_FILE="$(prompt_with_default "Input TE table path" "$INPUT_FILE")"
        OUTPUT_FILE="$(prompt_with_default "Output file path" "$OUTPUT_FILE")"
        FDR_THRESHOLD="$(prompt_with_default "FDR threshold" "$FDR_THRESHOLD")"
        TE_CUTOFF="$(prompt_with_default "TE cutoff" "$TE_CUTOFF")"
    }

    prompt_sequence

    while true; do
        show_summary
        read -rp "Edit step number (1-4) or press Enter to run: " EDIT_CHOICE
        EDIT_CHOICE="$(echo "$EDIT_CHOICE" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        case "$EDIT_CHOICE" in
            "")
                break
                ;;
            1|2|3|4)
                case "$EDIT_CHOICE" in
                    1) INPUT_FILE="$(prompt_with_default "Input TE table path" "$INPUT_FILE")" ;;
                    2) OUTPUT_FILE="$(prompt_with_default "Output file path" "$OUTPUT_FILE")" ;;
                    3) FDR_THRESHOLD="$(prompt_with_default "FDR threshold" "$FDR_THRESHOLD")" ;;
                    4) TE_CUTOFF="$(prompt_with_default "TE cutoff" "$TE_CUTOFF")" ;;
                esac
                ;;
            *)
                echo "[!] Unknown selection."
                ;;
        esac
    done
else
    if [ $# -lt 4 ]; then
        usage
        exit 1
    fi
    INPUT_FILE="$1"
    OUTPUT_FILE="$2"
    FDR_THRESHOLD="$3"
    TE_CUTOFF="$4"
fi

CMD_ARGS=(
    --input "$INPUT_FILE"
    --output_file "$OUTPUT_FILE"
    --fdr "$FDR_THRESHOLD"
    --te_cutoff "$TE_CUTOFF"
)

echo ""
echo "Running: $PYTHON -m code.make_GRN_averageTE ${CMD_ARGS[*]}"
"$PYTHON" -m code.make_GRN_averageTE "${CMD_ARGS[@]}"
