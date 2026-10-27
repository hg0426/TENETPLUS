"""Resolve feature types once and share a checked manifest across all consumers."""
import csv
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from code.reference_bundle import selected_reference, sha256, cis_mapping_policy
from code.path_utils import resolve_output

PEAK_PATTERN = re.compile(r'(?P<contig>[^\s:]+)[:-](?P<start>[0-9]+)-(?P<end>[0-9]+)')


def parse_peak(name):
    match = PEAK_PATTERN.fullmatch(str(name))
    if match is None or max(len(match['start']), len(match['end'])) > 19:
        return None
    start, end = int(match['start']), int(match['end'])
    if not 0 <= start < end <= 2**63-1:
        return None
    return match['contig'], start, end


def read_names(path):
    with open(path, encoding='utf-8') as handle:
        return [line.rstrip('\r\n') for line in handle]


def validate_names(names):
    if not names or any(not isinstance(n, str) or not n or any(c in n for c in '\t\r\n') for n in names):
        raise ValueError('Feature IDs must be nonempty strings without tab/newline characters')
    if len(names) != len(set(names)):
        raise ValueError('Duplicate feature IDs: RNA and peaks must have distinct IDs')


def explicit_types(path, names):
    if not path:
        return {}
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or not {'feature_id', 'feature_type'} <= set(reader.fieldnames):
            raise ValueError('Feature metadata requires feature_id and feature_type TSV columns')
        result = {}
        for row in reader:
            name, kind = row['feature_id'], row['feature_type']
            if name in result or kind not in {'RNA', 'peak'}:
                raise ValueError(f'Duplicate ID or invalid type in feature metadata: {name!r}')
            result[name] = kind
    if set(result) != set(names):
        raise ValueError('Feature metadata must cover exactly the current matrix feature IDs')
    return result


def classify(names, reference=None, types=None, mode='1'):
    validate_names(names)
    types = types or {}
    if types and set(types) != set(names):
        raise ValueError('Incomplete feature types')
    rows, errors = [], []
    for name in names:
        parsed = parse_peak(name)
        chromosomes = sorted(reference.genes.get(name, set())) if reference else []
        canonical = reference.aliases.get(parsed[0]) if reference and parsed else None
        kind = types.get(name)
        if kind is None:
            if mode == '0':
                kind = 'RNA'  # RNA-only explicitly declares the modality.
            elif chromosomes and canonical:
                errors.append({'feature_id': name, 'reason': 'Both annotated RNA and reference-contig coordinate; supply explicit feature types'})
                continue
            elif chromosomes:
                kind = 'RNA'
            elif parsed and canonical:
                kind = 'peak'
            elif parsed:
                errors.append({'feature_id': name, 'reason': 'Coordinate-shaped name on unknown contig; supply annotation or explicit RNA type'})
                continue
            elif re.search(r'[:-][+-]?[0-9]+-[+-]?[0-9]+$', name):
                errors.append({'feature_id': name, 'reason': 'Invalid or unsupported coordinate-shaped name; supply explicit RNA type if this is RNA'})
                continue
            else:
                kind = 'RNA'
        row = {'feature_id': name, 'feature_type': kind, 'gene_contigs': chromosomes}
        if kind == 'peak':
            if not parsed or not canonical:
                errors.append({'feature_id': name, 'reason': 'Peak needs valid coordinates on a reference contig'})
                continue
            if parsed[2] > reference.lengths.get(canonical, 2**63-1):
                errors.append({'feature_id': name, 'reason': 'Peak end exceeds reference contig length'})
                continue
            row.update(contig=canonical, start=parsed[1], end=parsed[2])
        rows.append(row)
    if mode == '0' and any(r['feature_type'] == 'peak' for r in rows):
        errors.append({'reason': 'Explicit peaks are incompatible with RNA-only mode'})
    if mode != '0' and not any(r['feature_type'] == 'peak' for r in rows):
        errors.append({'reason': 'No validated peaks in mixed input'})
    return rows, errors


def matrix_header(path):
    path = Path(path)
    if path.suffix.lower() == '.parquet':
        import pyarrow.parquet as pq
        schema = pq.read_schema(path)
        meta = json.loads((schema.metadata or {}).get(b'pandas', b'{}'))
        excluded = {x for x in meta.get('index_columns', []) if isinstance(x, str)}
        names = [n for n in schema.names if n not in excluded]
    else:
        with path.open(newline='', encoding='utf-8-sig') as handle:
            row = next(csv.reader(handle, delimiter='\t' if path.suffix.lower() in {'.tsv', '.txt'} else ','))
        names = row[1:]
    validate_names(names)
    return names


def prepare(names, mode, species=None):
    reference = selected_reference(required=mode != '0', species=species)
    policy=cis_mapping_policy(reference)
    from code.species_support import needs_tf,tf_list_path,load_tf_names
    tf_hash=None
    if needs_tf(mode):
        load_tf_names(species)
        tf_hash=sha256(tf_list_path(species))
    type_path = os.getenv('TENET_FEATURE_TYPES')
    types = explicit_types(type_path, names)
    rows, errors = classify(names, reference, types, mode)
    for row in rows:
        if row['feature_type']=='RNA':
            row['cis_contig']=(reference.last_contig.get(row['feature_id']) if policy in {'human_legacy', 'legacy'}
                               else row['gene_contigs'][0] if len(row['gene_contigs'])==1 else None)
    if mode in {'1', '5'} and os.getenv('TENET_PAIR_MODE', 'default') == 'default':
        for row in rows:
            if policy=='strict' and row['feature_type'] == 'RNA' and len(row['gene_contigs']) != 1:
                errors.append({'feature_id': row['feature_id'], 'reason': 'cis peak→RNA needs exactly one annotated contig', 'contigs': row['gene_contigs']})
    out = resolve_output('feature_manifest.json')
    out.parent.mkdir(parents=True, exist_ok=True)
    report = {'schema_version': 1, 'mode': mode, 'species': species,
              'reference_sha256': reference.fingerprint if reference else None,
              'reference_metadata': reference.metadata if reference else None,
              'feature_types_sha256': sha256(type_path) if type_path else None,
              'tf_list_sha256':tf_hash,'cis_mapping_policy':policy,
              'features': rows, 'errors': errors,
              'unmapped_RNA': [r['feature_id'] for r in rows if r['feature_type'] == 'RNA' and not r['gene_contigs']]}
    report['multi_contig_RNA']=[{'feature_id':r['feature_id'],'contigs':r['gene_contigs'],'selected_cis_contig':r['cis_contig']}
                               for r in rows if r['feature_type']=='RNA' and len(r['gene_contigs'])>1]
    if policy in {'human_legacy', 'legacy'} and mode in {'1','5'} and os.getenv('TENET_PAIR_MODE','default')=='default':
        print(f"[Legacy cis mapping WARNING] {len(report['unmapped_RNA'])} unmapped RNA receive no cis peak→RNA candidates; {len(report['multi_contig_RNA'])} multi-contig RNA use the last mapping row. Names and choices are recorded in feature_manifest.json.",file=sys.stderr)
    if errors:
        resolve_output('feature_validation_errors.json').write_text(json.dumps(report, indent=2, ensure_ascii=False))
        raise ValueError(f'{len(errors)} feature/reference errors; see feature_validation_errors.json. First: {errors[0]}')
    peaks = [r['feature_id'] for r in rows if r['feature_type'] == 'peak']
    peak_text = ''.join(n+'\n' for n in peaks)
    report['peak_list_sha256'] = hashlib.sha256(peak_text.encode()).hexdigest()
    peak_path = resolve_output('TE_peak_list.txt')
    tmp = peak_path.with_suffix('.tmp'); tmp.write_text(peak_text); tmp.replace(peak_path)
    tmp = out.with_suffix('.tmp'); tmp.write_text(json.dumps(report, indent=2, ensure_ascii=False)); tmp.replace(out)
    return report


def load_manifest(names):
    validate_names(names)
    path = resolve_output('feature_manifest.json')
    if not path.exists():
        raise ValueError('Missing feature_manifest.json; prepare features before TE or matrix export')
    data = json.loads(path.read_text())
    if data.get('schema_version') != 1 or data.get('errors'):
        raise ValueError('Invalid feature manifest')
    if [r['feature_id'] for r in data['features']] != list(names):
        raise ValueError('Feature manifest does not match current feature IDs/order')
    reference = selected_reference(required=data.get('reference_sha256') is not None, species=data.get('species'))
    if cis_mapping_policy(reference)!=data.get('cis_mapping_policy','strict'):
        raise ValueError('Prepared cis mapping policy differs from current selection')
    if reference and reference.fingerprint != data['reference_sha256']:
        raise ValueError('Prepared reference differs from selected reference')
    peaks = resolve_output('TE_peak_list.txt')
    if not peaks.exists() or sha256(peaks) != data['peak_list_sha256']:
        raise ValueError('Peak list checksum differs from feature manifest')
    expected = [r['feature_id'] for r in data['features'] if r['feature_type'] == 'peak']
    if read_names(peaks) != expected:
        raise ValueError('Peak list differs from feature manifest')
    return data


def load_peak_names(names, path=None):
    data = load_manifest(names)
    return {r['feature_id'] for r in data['features'] if r['feature_type'] == 'peak'}


def main():
    if sys.argv[1] == '--validate-prepared':
        data = load_manifest(read_names(resolve_output('gene_names')))
        mode, species = sys.argv[2:4]
        selected_reference(required=mode != '0', species=species)
        if data['mode'] != mode or data['species'] != species:
            raise ValueError('Prepared feature mode/species mismatch')
        path = os.getenv('TENET_FEATURE_TYPES')
        if (sha256(path) if path else None) != data['feature_types_sha256']:
            raise ValueError('Prepared feature metadata differs from current selection')
        from code.species_support import needs_tf,tf_list_path
        tf_hash=sha256(tf_list_path(species)) if needs_tf(mode) else None
        if tf_hash!=data.get('tf_list_sha256'):
            raise ValueError('Prepared TF list differs from current selection; regenerate candidate pairs')
        return
    path, mode, *species = sys.argv[1:]
    prepare(matrix_header(path), mode, species[0] if species else None)


if __name__ == '__main__':
    main()
