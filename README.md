# TENETPLUS

TENETPLUS reconstructs Transfer Entropy (TE)-based causal gene networks from pseudotime-ordered single-cell transcriptomic and epigenetic data, with an emphasis on RNA + ATAC multiome assays. This package provides the current TENETPLUS workflow for kernel-based TE analysis.

## Method

![tenetplus_workflow](docs/tenetplus_workflow.jpg)

The figure summarizes the study workflow. The runner implements the TE and
network steps described below; it does not automatically perform every
downstream analysis in the figure. In particular, cis candidate generation
uses the same contig, not a 1 Mb peak-to-gene distance cutoff. Distance-based
filtering and triplet reconstruction use the optional R post-processing described
in [advanced usage](docs/ADVANCED_USAGE.md).

## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014  
Original TENET repository: https://github.com/neocaleb/TENET

## Highlights

- **Multiome-aware outputs**: supports TF->gene, TF->peak, peak->gene, peak->peak, RNA-only, and all-pair modes.
- **Species-aware inputs**: bundled references for six species, with annotation-based RNA/peak classification and checks before TE calculation.
- **Consistent LocalTE scores**: storing LocalTE uses the same scalar TE and permutation calculation as LocalTE OFF in the validated paths.
- **Fast post-processing**: TE parquet outputs are converted directly into FDR-filtered `.sif` and parquet edge tables.
- **Clear network files**: FDR network filenames include the split name, e.g. `TE_TF_GN_fdr0.01.sif`.
- **Focused indirect trimming**: trimming is applied only to TF->gene when enabled; untrimmed FDR networks are always kept.

## Requirements

- Python 3.10 recommended (the release checks use Python 3.10.15)
- Standard Linux shell tools used by `TENET_Plus_for_py.sh`
- Python packages listed in `requirements.txt`
- Optional R post-processing: R and the packages listed in [advanced usage](docs/ADVANCED_USAGE.md#r-post-processing)

Install dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

`requirements.txt` pins the direct dependency versions used in the release
checks. It is not a complete lockfile for transitive dependencies. Bitwise
reproducibility across different Python/NumPy/library versions has not been
established.

## Directory Layout

- `TENET_Plus_for_py.sh`: main runner.
- `code/runner/`: shared defaults, interactive prompts, and replay recording.
- `code/run_config.py`: configuration-file reading and validation.
- `code/`: Python utilities, kernel TE engine, matrix splitting, LocalTE export, FDR network export, and indirect trimming.
- `input/`: example matrix, trajectory, selection mask, human/mouse TF lists, and human annotation files.
- `output/`: default output directory created and used by the runner.
- `references/`: species-specific annotation bundles and the reference catalog.
- `scripts/rebuild_references.py`: rebuilds the bundled references from their recorded sources.
- `docs/`: workflow figure, other-species setup, and optional analysis commands.
- `requirements.txt`: Python runtime dependencies.

## Input Files and Species References

Run commands from the package root. Provide a cells-by-features matrix, a
trajectory file, and a cell-selection file. CSV/TSV matrices use the first
column for cell IDs; Parquet matrices can use pandas index metadata. Feature
IDs must be unique strings. The included example has 1,448 cells and 73 features;
its matrix, trajectory, and selection mask are under `input/`.

Trajectory and cell-selection files are headerless, one value per line,
in the original matrix row order. Trajectory values determine ascending
pseudotime order; selection value 1 retains a cell and 0 excludes it. Both
files must have one row per input cell.

Mixed RNA/ATAC modes use a reference bundle matching the input species,
genome assembly, and RNA identifier format. The following bundles are included:

| Species argument | Bundled reference | TF list |
|---|---|---|
| human | GENCODE v44 / GRCh38.p14; preserves the original human mapping | Included |
| mouse | GENCODE M39 / GRCm39 primary assembly | Included |
| zebrafish | NCBI RefSeq GRCz13ab | User-provided |
| fly | FlyBase Release 6.54 / NCBI Release 6 plus ISO1 MT | User-provided |
| worm | WormBase WS298 / WBcel235 | User-provided |
| rat | NCBI RefSeq GRCr8, RS_2026_08 | User-provided |

The species argument selects the corresponding bundled reference. These are
fixed releases, not downloads of the latest annotation at run time. Use a
custom bundle if your input was built with another assembly or annotation.
Human/mouse TF lists are preserved from the existing package; replace them
when they do not match the identifiers in your matrix.

The bundled human reference uses `human_legacy`; other references use the same
`legacy` cis mapping rules. Unmapped RNA is omitted from peak->gene candidates;
multiple contig mappings use the last row. Both cases produce warnings and are
recorded in `feature_manifest.json`. TF->gene candidates are not removed for
these reasons. Set `TENET_CIS_MAPPING_POLICY=strict` to stop on either case in
cis modes 1/5.

RNA/peak classification uses annotation membership and contig coordinates.
Peak IDs may be `contig-start-end` or `contig:start-end`, with 0-based, half-open
coordinates. Ambiguous names require `TENET_FEATURE_TYPES` with explicit types.

The runner checks required references and TF lists before TE calculation.
For another species, assembly, or identifier scheme, follow
[other-species and reference setup](docs/OTHER_SPECIES.md).

## Step 1 - Run `TENET_Plus_for_py.sh`

The same runner handles both multiome and RNA-only kernel TE modes. Launch it
without arguments for guided prompts, or use a configuration file for a batch
run. Both paths use the same calculation workflow.

### Interactive Mode

```bash
./TENET_Plus_for_py.sh
```

The prompt asks for the input matrix, CPU jobs, trajectory file, cell-select
file, history length (`k`), species, network mode, permutation setting, LocalTE
storage, GRN/FDR export, and output directory. Review the summary and enter a
step number to edit it before starting. Choose a separate output directory for
each run.

### Run from a Configuration File

Generate an editable example, change its settings, and run:

```bash
./TENET_Plus_for_py.sh --write-config analysis.ini
./TENET_Plus_for_py.sh --check-config analysis.ini
./TENET_Plus_for_py.sh --config analysis.ini
```

The generated INI file points to the included example inputs. Its `[run]`
section contains the matrix, jobs, trajectory, cells, history, species, mode,
and output directory. Other sections control permutation, LocalTE, and network
export. No environment variables are needed for these settings.

`--check-config` checks settings and input-file existence without computing TE.
Matrix contents and species resources are checked when the run starts. Relative
input, reference, and output paths in an INI file are relative to that file.
Paths are literal text: do not add shell quotes or `$VARIABLE` expressions.

### Repeat or Adjust a Run

Every run saves its resolved settings to `run_config.ini` in its output directory.
Use that file to repeat the analysis with a new output directory:

```bash
./TENET_Plus_for_py.sh --config output/my_run/run_config.ini \
  --output output/repeated_run
```

To change the analysis, copy `run_config.ini`, edit it, and pass the copy to
`--config`. The saved `run_invocation.sh` is a shortcut that also accepts
`--output`. This repeats the workflow from the beginning; it is not a resume
command. Omitting `--output` reuses the saved output directory and may replace
previous files. Use a new directory for changed settings and concurrent runs.

For file-based runs, the configuration supplies TENET settings; exported TENET
variables from the current terminal are ignored. `--output` takes precedence
and resolves relative to the current working directory. Saved configurations
contain absolute input paths and the Python interpreter used by the run; update
these when moving to another machine. Reproducing results also requires the
same input contents, references, package, and software environment.

The [configuration guide](docs/ADVANCED_USAGE.md#configuration-files-and-reruns)
explains all sections, permutation modes, and prepared-input reuse.

### Existing Batch Commands

The seven-argument interface remains supported:

```bash
./TENET_Plus_for_py.sh input/merged_expression_data.csv \
  16 input/trajectory.txt input/cell_select.txt 1 human 1
```

The arguments are matrix, jobs, trajectory, cell selection, history length,
species, and network mode. Existing `TENET_*` / `LOCAL_TE_*` environment settings
remain supported for this interface and interactive runs. Defaults are kernel
TE, no permutation, no LocalTE storage, FDR network export on, and TF->gene
indirect trimming on.

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

By default, TE outputs are written under `output/`. Set the output directory
in the interactive prompt or the configuration file. For a file-based rerun,
use `--output` as shown above. In the existing environment-variable interface,
relative `TENET_OUTPUT_DIR` values are resolved from the package root. Generated
preprocessing files stay in that run's output directory. Use a different output
directory for each concurrent run. Depending on the selected mode, the runner
writes parquet tables such as:

- `TE_result_all.parquet`: full TE result table before split export.
- `TE_TF_GN.parquet`: TF->gene TE table.
- `TE_TF_PK.parquet`: TF->peak TE table.
- `TE_PK_GN.parquet`: peak->gene TE table.
- `TE_PK_PK.parquet`: peak->peak TE table.
- `TE_GN_GN.parquet`: gene->gene TE table.
- `TE_all_features.parquet`: all-feature TE table.

Each TE parquet table uses `Source`, `Target`, and `TE` columns. Permutation and LocalTE runs add their own output files when enabled.

`TE_result_all.parquet` uses integer 1-based feature IDs; the split TE tables
use feature names. Preserve `gene_names` and the feature manifest alongside
these outputs when exporting LocalTE or running exact permutation.

## Step 3 - Reconstruct FDR Networks

By default, `TENET_Plus_for_py.sh` converts the selected TE outputs into FDR-filtered network files immediately after TE calculation.

These networks use one-sided p-values from global TE z-scores within each
split table, followed by BH correction. They are distinct from per-edge
shuffle-based permutation p/q values in the optional permutation workflow.

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

Configure automatic network export in the INI file:

```ini
[network]
enabled = on
fdr = 0.01
```

Set `enabled = off` to skip this step, or change `fdr` to the desired threshold.

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

Indirect trimming is optional. In the INI file, use:

```ini
[network]
trim_indirect = on
trim_threshold = -0.01
```

Set `trim_indirect = off` to keep only untrimmed FDR networks. When enabled,
automatic trimming applies only to TF->gene.

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

Choose permutation in the interactive prompt or edit the INI file:

```ini
[permutation]
mode = grn_edges
count = 100
candidate_fdr = 0.01
```

This tests GRN-FDR candidates. Use `mode = all_pairs` to test all selected
pairs, or `off` to skip permutation. See the configuration guide for p/q
thresholds and the distinction between these testing scopes.

### Optimized Exact Permutation

For source-wise shuffle reuse, run the separate optimized runner on prepared
inputs and observed TE from a matching run:

```bash
TE_KERNEL_PERM_TABLE_SAMPLER=off \
python -m code.run_exact_permutation \
  --matrix output/my_run/cell_gene_trsps.parquet \
  --observed output/my_run/TE_result_all.parquet \
  --gene-names output/my_run/gene_names \
  --output output/my_exact_permutation \
  --mode optimized --perms 1000 --jobs 16 --seed 42
```

Permutation enabled through the main runner (interactive, INI, or environment)
uses the shell permutation path. The optimized runner above must be invoked
separately. Both paths include the LocalTE ON/OFF consistency fix.

The optimized runner requires finite prepared feature-by-cell data, 3–8193
cells, history length 1, kernel width 0.5, integer 1-based endpoints, and
matching feature manifest/peak list. Keep any explicitly selected reference
and mapping-policy settings from preprocessing. It checks observed TE for
bitwise agreement before calculating permutation counts.

The runner preserves the source-specific shuffle sequence and writes
`TE_kernel_perm.parquet` with exceedance counts, p-values, and BH q-values,
plus q-filtered edge tables. BH correction applies to the provided pair set.
Candidate-restricted and full-pair runs therefore use different testing sets.
Permutation p/q values also differ in definition from the score-based network
filtering in Steps 3–4.

Use a new output directory for a new analysis. Reusing the same directory
resumes only when code, input hashes, and settings match. Up to 32 workers
are supported; `--mode baseline` provides a comparison path with the same
inputs and shuffle settings.

## Optional - LocalTE Storage

Enable LocalTE storage interactively or edit the INI file:

```ini
[localte]
enabled = on
```

LocalTE storage uses the same scalar TE/permutation path as storage OFF.
Stored arrays default to float16; their average must not be substituted for
scalar TE. Outputs are written under `local_te_split_chunks/`, split by network
type. Worker parts are kept by default; set `LOCAL_TE_MERGE_PARTS=on` to merge
them into single chunk files. Additional formats and TF activity analysis are
covered in [advanced usage](docs/ADVANCED_USAGE.md).

### Export LocalTE for Selected Edges

Lazy mode first calculates scalar TE/permutation results and saves a recipe.
Export LocalTE later for a selected edge table:

In the existing INI file, set:

```ini
[localte]
enabled = on
lazy = on
```

Run it with `--config` and a separate output directory. Then export selected
edges (adjust the paths to that run):

```bash
python -m code.localte_lazy_export \
  --recipe output/my_lazy_run/local_te_recipe/recipe.json \
  --edges output/my_lazy_run/network_edges.parquet \
  --output_dir output/my_selected_localte
```

Lazy export requires history length 1 and validates the feature manifest and
canonical TE scores. CSV/TSV selectors use feature names; use integer
1-based Source/Target columns in Parquet for numeric endpoints. Exported
scalar scores retain float64 precision; LocalTE arrays use float16 by default.
Use a new or empty export directory to avoid mixing outputs from different selections.
Do not substitute the average of a quantized LocalTE array for the scalar TE.

Lazy export saves only the selected edges, with an additional recomputation
step. Use it when full edge-by-timepoint output is unnecessary.

## Optional - Peak–Gene Distance and Triplets

`input/gene_anno.bed` contains the original human gene intervals, with columns
`gene_chr`, `gene_start`, `gene_end`, and `gene`. It is used for R distance
filtering and TF–peak–gene triplets, not for TE calculation. `gene_chr.txt`
contains only contig/name mappings and cannot replace this coordinate table.
Use matching gene coordinates for other species or genome builds; see
[the R post-processing instructions](docs/ADVANCED_USAGE.md#r-post-processing).
