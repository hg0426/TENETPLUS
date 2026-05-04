# TENETPLUS

TENETPLUS reconstructs Transfer Entropy (TE)-based causal gene networks from pseudotime-ordered single-cell transcriptomic and epigenetic data, with an emphasis on RNA + ATAC multiome assays. This package provides the current TENETPLUS workflow for kernel-based TE analysis.

## Method

![tenetplus_workflow](docs/tenetplus_workflow.jpg)

## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014  
Original TENET repository: https://github.com/neocaleb/TENET

## Highlights

- **Multiome-aware outputs**: supports TF->gene, TF->peak, peak->gene, peak->peak, RNA-only, and all-pair modes.
- **Fast post-processing**: TE parquet outputs are converted directly into FDR-filtered `.sif` and parquet edge tables.
- **Clear network files**: FDR network filenames include the split name, e.g. `TE_TF_GN_fdr0.01.sif`.
- **Focused indirect trimming**: trimming is applied only to TF->gene when enabled; untrimmed FDR networks are always kept.

## Requirements

- `python3`
- Python packages listed in `requirements.txt`

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Directory Layout

- `TENET_Plus_for_py.sh`: main runner.
- `code/`: Python utilities, kernel TE engine, matrix splitting, LocalTE export, FDR network export, and indirect trimming.
- `input/`: optional place for input matrices, pseudotime files, and cell-select masks.
- `output/`: default output directory created and used by the runner.
- `docs/`: workflow figure and documentation assets.
- `requirements.txt`: Python runtime dependencies.

## Step 1 - Run `TENET_Plus_for_py.sh`

The shell wrapper is the entry point for both multiome and RNA-only kernel TE modes. Launch it without arguments for the guided prompt sequence, or supply the arguments below to run non-interactively.

### Interactive Mode

```bash
./TENET_Plus_for_py.sh
```

The prompt asks for the input matrix, CPU jobs, trajectory file, cell-select file, history length (`k`), species, output mode, permutation setting, LocalTE storage, and GRN/FDR network export.

### Batch Usage

```bash
./TENET_Plus_for_py.sh \
  <input_matrix> \
  <num_jobs> \
  <trajectory_file> \
  <cell_select_file> \
  <history_k> \
  <species> \
  <network_mode>
```

Example:

```bash
./TENET_Plus_for_py.sh \
  input/merged_expression_data.csv \
  32 \
  input/trajectory.txt \
  input/cell_select.txt \
  1 human 1
```

The compact command above runs the included example dataset with the current defaults: kernel TE, no permutation, no LocalTE storage, FDR network export on, and TF->gene indirect trimming on.

### Optional Toggles for Scripts

The kernel package keeps the positional command short. Use environment
variables for optional behavior.

Common examples:

```bash
# Write outputs to a separate run directory.
TENET_OUTPUT_DIR=output/my_run \
./TENET_Plus_for_py.sh input/merged_expression_data.csv 32 input/trajectory.txt input/cell_select.txt 1 human 1

# Skip automatic FDR network generation.
TENET_MAKE_GRN=off \
./TENET_Plus_for_py.sh input/merged_expression_data.csv 32 input/trajectory.txt input/cell_select.txt 1 human 1

# Enable LocalTE storage.
TENET_STORE_LOCAL_TE=on \
./TENET_Plus_for_py.sh input/merged_expression_data.csv 32 input/trajectory.txt input/cell_select.txt 1 human 1

# Enable permutation on GRN-FDR candidates.
TENET_PERMUTE=on \
TENET_PERM_N=100 \
TENET_PERM_CANDIDATE_GRN_FDR=0.01 \
./TENET_Plus_for_py.sh input/merged_expression_data.csv 32 input/trajectory.txt input/cell_select.txt 1 human 1
```

Useful environment variables:

- `TENET_OUTPUT_DIR`: output directory; relative paths are resolved from the package root.
- `TENET_MAKE_GRN`: `on` or `off` (default `on`).
- `TENET_GRN_FDR`: FDR alpha for generated networks (default `0.01`).
- `TENET_TRIM_INDIRECT`: `on` or `off`; trims TF->gene only when enabled.
- `TENET_STORE_LOCAL_TE`: `on` or `off` (default `off`).
- `TENET_PERMUTE`: `on` or `off` (default `off`).
- `TENET_PERM_N`: number of permutations (default `100`).
- `TENET_PERM_FDR`: `on` or `off`; use permutation q-value threshold when enabled.
- `TENET_PERM_Q_ALPHA`: permutation q-value threshold (default `0.05`).
- `TENET_PERM_ALPHA`: permutation p-value threshold (default `0.01`).
- `TENET_RESULTS_BUFFER_ROWS`: TE result write buffer size; leave unset for auto.

## Network Modes

```text
0 - RNA-only TENET_TF.
1 - TENET_Plus full: TF->gene, TF->peak, peak->gene.
2 - TF->gene only.
3 - TF->peak only.
4 - TF->gene + TF->peak.
5 - peak->gene (cis).
6 - peak->peak (cis).
7 - RNA-only every gene -> every gene.
8 - every gene -> every gene.
9 - every gene/peak -> every gene/peak.
```

Modes `7`, `8`, and `9` are exhaustive all-pair modes and can be very large.

## Step 2 - TE Outputs

By default, TE outputs are written under `output/`. To isolate large runs or keep multiple runs side by side, set `TENET_OUTPUT_DIR`:

```bash
TENET_OUTPUT_DIR=/path/to/run_output ./TENET_Plus_for_py.sh ...
```

Relative `TENET_OUTPUT_DIR` values are resolved from the package root. Depending on the selected mode, the runner writes parquet tables such as:

- `TE_result_all.parquet`: full TE result table before split export.
- `TE_TF_GN.parquet`: TF->gene TE table.
- `TE_TF_PK.parquet`: TF->peak TE table.
- `TE_PK_GN.parquet`: peak->gene TE table.
- `TE_PK_PK.parquet`: peak->peak TE table.
- `TE_GN_GN.parquet`: gene->gene TE table.
- `TE_all_features.parquet`: all-feature TE table.

Each TE parquet table uses `Source`, `Target`, and `TE` columns. Permutation and LocalTE runs add their own output files when enabled.

## Step 3 - Reconstruct FDR Networks

By default, `TENET_Plus_for_py.sh` converts the selected TE outputs into FDR-filtered network files immediately after TE calculation.

For mode `1`, typical outputs are:

```text
output/network_TE_TF_GN_fdr0.01/TE_TF_GN_fdr0.01.parquet
output/network_TE_TF_GN_fdr0.01/TE_TF_GN_fdr0.01.sif
output/network_TE_TF_GN_fdr0.01/TE_TF_GN_fdr0.01.trimIndirect-0.01.sif
output/network_TE_TF_PK_fdr0.01/TE_TF_PK_fdr0.01.parquet
output/network_TE_TF_PK_fdr0.01/TE_TF_PK_fdr0.01.sif
output/network_TE_PK_GN_fdr0.01/TE_PK_GN_fdr0.01.parquet
output/network_TE_PK_GN_fdr0.01/TE_PK_GN_fdr0.01.sif
output/grn_network_outputs.tsv
```

The FDR edge parquet files contain only:

```text
Source, Target, TE
```

To skip automatic network generation:

```bash
TENET_MAKE_GRN=off ./TENET_Plus_for_py.sh ...
```

To change the FDR threshold:

```bash
TENET_GRN_FDR=0.05 ./TENET_Plus_for_py.sh ...
```

## Step 4 - Manual FDR Network Export

If you already have a TE parquet table, you can create an FDR network manually:

```bash
python -m code.network_from_te_parquet \
  output/TE_TF_GN.parquet \
  --outdir output/network_TE_TF_GN_fdr0.01 \
  --output-prefix TE_TF_GN \
  --fdr 0.01 \
  --te-cutoff 0 \
  --threads 32
```

This writes:

```text
TE_TF_GN_fdr0.01.parquet
TE_TF_GN_fdr0.01.sif
network_summary.txt
```

## Step 5 - Indirect Trimming

Indirect trimming is optional and controlled by `TENET_TRIM_INDIRECT`.

```bash
TENET_TRIM_INDIRECT=on  ./TENET_Plus_for_py.sh ...  # default; trim TF->gene only
TENET_TRIM_INDIRECT=off ./TENET_Plus_for_py.sh ...  # keep only untrimmed FDR networks
```

Manual trimming:

```bash
python code/trim_indirect.py \
  output/network_TE_TF_GN_fdr0.01/TE_TF_GN_fdr0.01.sif \
  -0.01
```

This writes:

```text
TE_TF_GN_fdr0.01.trimIndirect-0.01.sif
```

The untrimmed `TE_TF_GN_fdr0.01.sif` is always kept.

## Step 6 - Summarize Outdegree

```bash
python code/countOutdegree.py output/network_TE_TF_GN_fdr0.01/TE_TF_GN_fdr0.01.sif
python code/countOutdegree.py output/network_TE_TF_GN_fdr0.01/TE_TF_GN_fdr0.01.trimIndirect-0.01.sif
```

Each command writes `<network>.outdegree.txt`.

## Optional - Permutation

Permutation can be enabled through environment variables or the interactive prompt. For faster candidate-restricted permutation, use GRN-FDR candidates:

```bash
TENET_PERMUTE=on \
TENET_PERM_N=100 \
TENET_PERM_CANDIDATE_GRN_FDR=0.01 \
./TENET_Plus_for_py.sh \
  <input_matrix> <num_jobs> <trajectory_file> <cell_select_file> <history_k> <species> <network_mode>
```

## Optional - LocalTE Storage

LocalTE storage can be enabled interactively or with `TENET_STORE_LOCAL_TE=on`. It creates larger outputs and can be tuned with:

- `LOCAL_TE_CHUNK_SIZE` (default `300`)
- `LOCAL_TE_VALUES_DTYPE` (`float16` by default)
- `LOCAL_TE_SPLIT_EXPORT` (`on` by default)
- `LOCAL_TE_SPLIT_OUTPUT_DIR` (default `local_te_split_chunks`)
- `LOCAL_TE_MERGE_PARTS` (`off` by default; keeps equivalent worker-part datasets instead of merging them into one chunk file)

Most users only need to set `<num_jobs>` in the main command. LocalTE export automatically uses that value:

- `LOCAL_TE_EXPORT_WORKERS=auto` means `min(<num_jobs>, 16)`
- `LOCAL_TE_MERGE_WORKERS=auto` means `min(<num_jobs>, 32)`

Set `LOCAL_TE_MERGE_PARTS=on` only when you need the older single-file-per-chunk layout (`chunk_0000.parquet`, etc.). The default `off` avoids the expensive merge step and is faster/lower-memory for large LocalTE outputs.

Advanced tuning is still available through environment variables when needed:

- `LOCAL_TE_ADVANCED=on` exposes buffer/worker/read settings in interactive mode
- `LOCAL_TE_BUFFER_EDGES` (default `10000`)
- `LOCAL_TE_READ_BATCH_ROWS` (default `8192`)
- `LOCAL_TE_USE_THREADS` (`on` by default)
