# Save resolved settings for interactive and batch runs through the same path.
write_replay_logs() {
    local replay_mode="$7"
    case "$TENET_PAIR_MODE:$7" in
        gene_only:0) replay_mode=7 ;;
        gene_only:1) replay_mode=8 ;;
        all_feature:2) replay_mode=9 ;;
    esac
    local name
    local settings=()
    while IFS= read -r name; do
        [[ -v "$name" ]] || continue
        case "$name" in
            TENET_STAGE_FILE|MEMORY_MONITOR_GPU) continue ;;
            TENET_*|TE_*|LOCAL_TE_*|MEMORY_MONITOR_*|PYTHON|PYTHONNOUSERSITE|OMP_NUM_THREADS|MKL_NUM_THREADS|OPENBLAS_NUM_THREADS|NUMEXPR_NUM_THREADS)
                settings+=("$name=${!name}") ;;
        esac
    # Inherited/exported advanced settings already reach Python through env.
    # Add the named user settings, but do not turn shell-only internal defaults
    # into exported overrides on replay (Python may have different defaults).
    done < <("$PYTHON" "$CODE_DIR/run_config.py" setting-names)
    env "${settings[@]}" TENET_OUTPUT_DIR="$OUTPUT_DIR" \
        TENET_PERMUTE="$PERMUTE_CHOICE" TENET_PERM_N="$PERM_N" \
        TENET_PERM_FDR="$PERM_FDR" TENET_PERM_ALPHA="$PERM_ALPHA" \
        TENET_PERM_Q_ALPHA="$PERM_Q_ALPHA" TENET_STORE_LOCAL_TE="$STORE_LOCAL_TE_CHOICE" \
        "$PYTHON" "$CODE_DIR/run_config.py" snapshot "$OUTPUT_DIR/run_config.ini" \
        "$1" "$2" "$3" "$4" "$5" "$6" "$replay_mode"
    echo "[Run] Settings saved: $OUTPUT_DIR/run_config.ini"
    echo "[Run] Replay script: $OUTPUT_DIR/run_invocation.sh"
    echo "[Run] Rerun into a new directory with:"
    printf '  bash %q --config %q --output NEW_DIRECTORY\n' "$REPO_ROOT/TENET_Plus_for_py.sh" "$OUTPUT_DIR/run_config.ini"
}
