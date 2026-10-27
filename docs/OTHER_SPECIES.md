# Other species and reference setup

Use the species, genome assembly, annotation release, and RNA identifiers that
match your input. A newer annotation does not make coordinates from an older
assembly compatible. The bundled releases are listed in the main
[README](../README.md#input-files-and-species-references); source URLs and
checksums are recorded in each bundle's `reference.json`.

## Required files

| Network mode | Reference bundle | TF list |
|---|---|---|
| 0: RNA-only TF->RNA | Optional | Required |
| 1/2/3/4: TF modes with mixed inputs | Required | Required |
| 5/6: cis peak modes | Required | Not required |
| 7: RNA-only all pairs | Optional | Not required |
| 8/9: mixed-input all pairs | Required | Not required |

Each reference bundle contains:

- `gene_chr.txt`: headerless `contig<TAB>RNA identifier` mappings.
- `contigs.json`: contig names, explicit aliases, and lengths where available.
- `reference.json`: species, assembly, annotation sources, and checksums.

`gene_chr.txt` is generated when building the bundle, not during each TE run.
The copy in `input/` preserves the original human mapping; the runner uses the
selected bundle. Human preserves the original 62,700 gene-name rows. Other
bundles include gene IDs, names, and supplied aliases, so their row counts are
not gene counts. Human/mouse TF lists are preserved from the original package.
Zebrafish, fly, worm, and rat need a user-provided TF list for TF-based modes.

## Build from GTF or GFF3

```bash
python -m code.reference_bundle \
  --annotation /data/reference/genes.gtf.gz \
  --format gtf \
  --output references/my_species_my_build \
  --species 'My species' \
  --assembly my_build \
  --release my_annotation_release \
  --source-url 'https://provider.example/exact_release/genes.gtf.gz'

export TENET_REFERENCE=/absolute/path/to/references/my_species_my_build
export TENET_TF_LIST=/absolute/path/to/my_species_TFs.txt
```

Use `--format gff3` for GFF3. GTF gene IDs/names and exon-only GTF are supported;
GFF3 uses gene/pseudogene records and their declared IDs/names/aliases. Supply
`--fai genome.fa.fai` or `--assembly-report assembly_report.txt` to add contig
lengths and explicit aliases. Without lengths, the peak-end upper bound cannot
be checked. Contig aliases are taken from the source, not inferred by adding
or removing `chr`. `--assembly-unit` selects exact assembly-report unit names.

The TF list has one RNA ID per line, without a header or duplicates. IDs must
match the matrix: use symbols for a symbol matrix, or the same stable-ID scheme
for an ID matrix. Keep the TF source/version with your analysis.

```bash
python -m code.species_support \
  --species 'My species' --mode 1 --matrix /data/cells_by_features.csv

TENET_OUTPUT_DIR=/data/my_run \
bash TENET_Plus_for_py.sh \
  /data/cells_by_features.csv 8 /data/trajectory.txt /data/cell_select.txt \
  1 'My species' 1
```

The species must match the bundle; the six bundled species also accept the
short names in the README. Missing required files, species mismatches, and
TF lists without matrix overlap stop execution before TE calculation.
Do not edit bundle files by hand: checksums are verified during loading.

To rebuild a bundled reference from its recorded sources:

```bash
python scripts/rebuild_references.py \
  --species mouse --cache /data/reference_sources --output /data/rebuilt_references
# --species all rebuilds all six bundles.
```

## Feature types and cis mapping

Peak IDs use `contig-start-end` or `contig:start-end`, a known contig, and
0 <= start < end in 0-based, half-open coordinates. Labels such as `peak_123`
do not contain enough information. If an ID could be RNA or a peak, provide
`TENET_FEATURE_TYPES=/data/types.tsv` with every matrix feature:

```text
feature_id	feature_type
chr1-10-20	RNA
chr1-30-40	peak
```

Types are `RNA` or `peak`. Duplicate feature IDs must first be distinguished
using their original annotation; a type table cannot distinguish two columns
with the same ID.

Default `legacy` mapping excludes unmapped RNA from cis peak->gene candidates
and selects the last row for multi-contig RNA. Warnings and
`feature_manifest.json` record the choices. This does not remove these RNA
features from TF->gene candidates. The last row is a compatibility rule, not
a preferred biological locus. `human_legacy` names the same policy for the
verified original human bundle.

Set `TENET_CIS_MAPPING_POLICY=strict` to reject missing/multiple cis mappings
in modes 1/5. Missing references, invalid peaks, and unresolved feature types
are errors under either policy. After changing reference or mapping settings,
start a new output directory and rerun preprocessing.

## Coordinates for downstream distance filtering

The bundles above supply contig mappings, not gene intervals. The optional R
post-processing needs a separate coordinate table matching your species,
assembly, and matrix IDs. Only the original human table is bundled as
`input/gene_anno.bed`; do not use it for mouse or another species. See
[R post-processing](ADVANCED_USAGE.md#r-post-processing) for its format.
