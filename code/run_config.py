"""Read, validate and save runner settings without evaluating shell code."""

from __future__ import annotations

import argparse
import configparser
import math
import os
from pathlib import Path
import shlex
import sys

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "TENET_Plus_for_py.sh"
CORE = ("matrix", "jobs", "trajectory", "cells", "history", "species", "mode")
# Defaults belong to runner/defaults.sh and the runner's normalisation step.
# Omitted config entries use those defaults, never an ambient TENET setting.
FIELDS = {
    "run": {"output": "TENET_OUTPUT_DIR"},
    "reference": {"bundle": "TENET_REFERENCE", "tf_list": "TENET_TF_LIST",
                  "feature_types": "TENET_FEATURE_TYPES", "cis_policy": "TENET_CIS_MAPPING_POLICY"},
    "permutation": {"count": "TENET_PERM_N", "fdr": "TENET_PERM_FDR",
                    "p_alpha": "TENET_PERM_ALPHA", "q_alpha": "TENET_PERM_Q_ALPHA",
                    "candidate_fdr": "TENET_PERM_CANDIDATE_GRN_FDR",
                    "candidate_te_cutoff": "TENET_PERM_CANDIDATE_TE_CUTOFF"},
    "localte": {"enabled": "TENET_STORE_LOCAL_TE", "chunk_size": "LOCAL_TE_CHUNK_SIZE",
                "values_dtype": "LOCAL_TE_VALUES_DTYPE", "split_export": "LOCAL_TE_SPLIT_EXPORT",
                "split_output": "LOCAL_TE_SPLIT_OUTPUT_DIR", "virtual": "LOCAL_TE_VIRTUAL",
                "lazy": "LOCAL_TE_LAZY", "buffer_edges": "LOCAL_TE_BUFFER_EDGES",
                "export_workers": "LOCAL_TE_EXPORT_WORKERS", "merge_workers": "LOCAL_TE_MERGE_WORKERS",
                "read_batch_rows": "LOCAL_TE_READ_BATCH_ROWS", "use_threads": "LOCAL_TE_USE_THREADS"},
    "network": {"enabled": "TENET_MAKE_GRN", "fdr": "TENET_GRN_FDR",
                "te_cutoff": "TENET_GRN_TE_CUTOFF", "trim_indirect": "TENET_TRIM_INDIRECT",
                "trim_threshold": "TENET_TRIM_THRESHOLD", "threads": "TENET_GRN_THREADS"},
    "reuse": {"mode": "TENET_REUSE_PREPARED", "directory": "TENET_PREPARED_DIR",
              "link_mode": "TENET_PREPARED_LINK_MODE"},
    "runtime": {"python": "PYTHON", "intermediate_save": "TENET_INTERMEDIATE_SAVE",
                "batch_size": "TE_BATCH_SIZE", "buffer_rows": "TENET_RESULTS_BUFFER_ROWS"},
}
PATH_VARS = {"TENET_OUTPUT_DIR", "TENET_REFERENCE", "TENET_TF_LIST",
             "TENET_FEATURE_TYPES", "TENET_PREPARED_DIR"}
PREFIXES = ("TENET_", "LOCAL_TE_", "TE_", "MEMORY_MONITOR_")
EXTRA_ENV = {"PYTHON", "PYTHONNOUSERSITE", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"}
EXCLUDE = {"TENET_STAGE_FILE", "MEMORY_MONITOR_GPU"}


def is_setting(key):
    return (key.startswith(PREFIXES) or key in EXTRA_ENV) and key not in EXCLUDE


def parser():
    cfg = configparser.ConfigParser(interpolation=None)
    cfg.optionxform = str
    return cfg


def absolute(value, base):
    path = Path(value).expanduser()
    return str((base / path).resolve()) if not path.is_absolute() else str(path.resolve())


def validate_env(env):
    for key in ("TENET_PERM_N", "LOCAL_TE_CHUNK_SIZE", "TENET_GRN_THREADS",
                "TE_BATCH_SIZE", "TENET_RESULTS_BUFFER_ROWS"):
        if env.get(key) and (not env[key].isdigit() or int(env[key]) < 1):
            raise ValueError(f"{key} must be a positive integer")
    for key in ("TENET_PERM_ALPHA", "TENET_PERM_Q_ALPHA", "TENET_GRN_FDR",
                "TENET_PERM_CANDIDATE_GRN_FDR"):
        if env.get(key):
            value = float(env[key])
            if not 0 < value <= 1:
                raise ValueError(f"{key} must be in (0, 1]")
    for key in ("TENET_GRN_TE_CUTOFF", "TENET_TRIM_THRESHOLD", "TENET_PERM_CANDIDATE_TE_CUTOFF"):
        if env.get(key) and not math.isfinite(float(env[key])):
            raise ValueError(f"{key} must be finite")
    for key in ("TENET_PERM_FDR", "TENET_STORE_LOCAL_TE", "TENET_MAKE_GRN",
                "LOCAL_TE_SPLIT_EXPORT", "LOCAL_TE_VIRTUAL", "LOCAL_TE_LAZY",
                "TENET_INTERMEDIATE_SAVE"):
        if key in env and env[key].lower() not in {"on", "off", "1", "0", "yes", "no", "true", "false"}:
            raise ValueError(f"{key} must be on or off")
    choices = {"TENET_REUSE_PREPARED": {"off", "on", "force"},
               "TENET_PREPARED_LINK_MODE": {"symlink", "copy"},
               "TENET_CIS_MAPPING_POLICY": {"legacy", "human_legacy", "strict"},
               "LOCAL_TE_VALUES_DTYPE": {"float16", "float32"},
               "TENET_TRIM_INDIRECT": {"on", "off", "all", "tf_gene", "tf_gn", "tf-gn", "default", "auto", "1", "0", "yes", "no", "true", "false", "none"}}
    for key, allowed in choices.items():
        if env.get(key) and env[key].lower() not in allowed:
            raise ValueError(f"{key} must be one of {', '.join(sorted(allowed))}")


def load_config(path, output=None):
    path = Path(path).resolve()
    cfg = parser()
    with path.open() as handle:
        cfg.read_file(handle)
    if cfg.defaults():
        raise ValueError("[DEFAULT] is not supported; put settings in their named sections")
    if "run" not in cfg:
        raise ValueError("Missing [run] section")
    env = {k: v for k, v in os.environ.items() if not is_setting(k)}
    # Use the interpreter that launched this reader unless the config selects one.
    env["PYTHON"] = sys.executable
    explicit = {}
    for section in cfg.sections():
        if section not in FIELDS and section != "advanced":
            raise ValueError(f"Unknown section [{section}]")
        for key, value in cfg[section].items():
            if section == "run" and key in CORE:
                continue
            if section == "permutation" and key == "mode":
                continue
            if section == "advanced":
                if not is_setting(key) or key in {"TENET_PERMUTE", "TENET_PAIR_MODE"}:
                    raise ValueError(f"Unsupported advanced setting: {key}")
                variable = key
            else:
                if key not in FIELDS[section]:
                    raise ValueError(f"Unknown setting [{section}] {key}")
                variable = FIELDS[section][key]
            if variable in explicit:
                raise ValueError(f"Setting specified twice: {variable}")
            explicit[variable] = value
    explicit = {k: v for k, v in explicit.items() if v}
    for key in CORE:
        if not cfg["run"].get(key, "").strip():
            raise ValueError(f"Missing [run] {key}")
    args = [cfg["run"][key] for key in CORE]
    for index in (0, 2, 3):
        args[index] = absolute(args[index], path.parent)
        if not Path(args[index]).is_file():
            raise ValueError(f"Input file does not exist: {args[index]}")
    if not args[1].isdigit() or int(args[1]) < 1:
        raise ValueError("[run] jobs must be a positive integer")
    if args[4] != "auto" and (not args[4].isdigit() or int(args[4]) < 1):
        raise ValueError("[run] history must be a positive integer or auto")
    if args[6] not in {str(i) for i in range(10)}:
        raise ValueError("[run] mode must be a number from 0 to 9 (see README)")
    if any("\n" in v or "\r" in v for v in args + list(explicit.values())):
        raise ValueError("Settings must fit on a single line")
    for key in PATH_VARS & explicit.keys():
        if explicit[key]:
            explicit[key] = absolute(explicit[key], path.parent)
    if explicit.get("PYTHON") and "/" in explicit["PYTHON"]:
        explicit["PYTHON"] = absolute(explicit["PYTHON"], path.parent)
    mode = cfg.get("permutation", "mode", fallback="off")
    if mode not in {"off", "grn_edges", "all_pairs"}:
        raise ValueError("[permutation] mode must be off, grn_edges or all_pairs")
    explicit["TENET_PERMUTE"] = "off" if mode == "off" else "on"
    if mode == "grn_edges":
        if not explicit.get("TENET_PERM_CANDIDATE_GRN_FDR"):
            raise ValueError("grn_edges requires [permutation] candidate_fdr")
    else:
        explicit.pop("TENET_PERM_CANDIDATE_GRN_FDR", None)
        explicit.pop("TENET_PERM_CANDIDATE_TE_CUTOFF", None)
    if output is not None:
        explicit["TENET_OUTPUT_DIR"] = absolute(output, Path.cwd())
    elif "TENET_OUTPUT_DIR" not in explicit:
        explicit["TENET_OUTPUT_DIR"] = str(path.parent / "output")
    validate_env(explicit)
    env.update(explicit)
    return args, env


def save_config(path, args):
    cfg = parser()
    cfg["run"] = dict(zip(CORE, args))
    settings = {k: v for k, v in os.environ.items() if is_setting(k) and v}
    # The canonical network mode already encodes the pair scope.
    settings.pop("TENET_PAIR_MODE", None)
    settings["PYTHON"] = sys.executable
    for key in PATH_VARS & settings.keys():
        settings[key] = absolute(settings[key], Path.cwd())
    for section, fields in FIELDS.items():
        if section not in cfg:
            cfg[section] = {}
        for name, variable in fields.items():
            if variable in settings:
                cfg[section][name] = settings.pop(variable)
    enabled = settings.pop("TENET_PERMUTE", "off") in {"on", "ON", "1"}
    candidate = cfg["permutation"].get("candidate_fdr", "0")
    cfg["permutation"]["mode"] = ("grn_edges" if float(candidate) > 0 else "all_pairs") if enabled else "off"
    cfg["advanced"] = dict(sorted(settings.items()))
    path = Path(path)
    with path.open("w") as handle:
        handle.write("# Resolved TENET+ settings. Paths are absolute for reruns on this machine.\n")
        handle.write("# Use --output NEW_DIRECTORY to keep the original results.\n")
        cfg.write(handle)
    # Check that the saved form can actually be read before advertising replay.
    load_config(path)
    command = ["bash", str(RUNNER), "--config", str(path.resolve())]
    (path.parent / "run_invocation.txt").write_text(shlex.join(command) + "\n")
    replay = path.parent / "run_invocation.sh"
    replay.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' + shlex.join(command) + ' "$@"\n')
    replay.chmod(0o755)


def write_template(path):
    path = Path(path).resolve()
    cfg = parser()
    cfg["run"] = {
        "matrix": os.path.relpath(ROOT / "input/merged_expression_data.csv", path.parent),
        "jobs": "16", "trajectory": os.path.relpath(ROOT / "input/trajectory.txt", path.parent),
        "cells": os.path.relpath(ROOT / "input/cell_select.txt", path.parent),
        "history": "1", "species": "human", "mode": "1", "output": "output/my_run"}
    cfg["permutation"] = {"mode": "off", "count": "100"}
    cfg["localte"] = {"enabled": "off"}
    cfg["network"] = {"enabled": "on", "fdr": "0.01", "trim_indirect": "on"}
    with path.open("x") as handle:
        handle.write("# Edit this file, then run TENET_Plus_for_py.sh --config FILE.\n")
        handle.write("# Relative input/reference/output paths are relative to this file.\n")
        handle.write("# Paths are literal text: do not quote them or use shell variables.\n")
        cfg.write(handle)
    print(f"Example configuration written: {path}")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "setting-names":
        print("\n".join(sorted({variable for fields in FIELDS.values() for variable in fields.values()})))
        return 0
    if len(sys.argv) > 1 and sys.argv[1] == "snapshot":
        save_config(sys.argv[2], sys.argv[3:])
        return 0
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--config", type=Path, help="Run without prompts using an INI file")
    group.add_argument("--check-config", type=Path, help="Validate settings and input paths without running TE")
    group.add_argument("--write-config", type=Path, help="Write an editable example; does not overwrite")
    ap.add_argument("--output", help="Override output directory (relative to the current directory)")
    ns = ap.parse_args()
    if ns.write_config:
        if ns.output:
            ap.error("--output requires --config or --check-config")
        write_template(ns.write_config)
        return 0
    args, env = load_config(ns.config or ns.check_config, ns.output)
    if ns.check_config:
        print("Settings and input paths are valid. Matrix contents and species resources are checked at run time.")
        print("Output:", env["TENET_OUTPUT_DIR"])
        print("Command:", shlex.join(["bash", str(RUNNER), *args]))
        return 0
    # Replace this process so signals and exit status belong to the actual runner.
    os.execvpe("bash", ["bash", str(RUNNER), *args], env)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, configparser.Error) as error:
        print(f"Configuration error: {error}", file=sys.stderr)
        sys.exit(2)
