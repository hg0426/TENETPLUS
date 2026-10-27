"""Build and load explicit, versioned annotation bundles (stdlib only)."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from pathlib import Path
from urllib.parse import unquote


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def text_open(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def attributes(text, fmt):
    if fmt == 'gtf':
        return {k: [v] for k, v in re.findall(r'(\S+)\s+"([^"]*)"\s*;?', text)}
    result = {}
    for part in text.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            # Split before percent-decoding: %2C is a literal comma in a name.
            result[key] = [unquote(v) for v in value.split(',')]
    return result


def build(annotation, output, species, assembly, release, source_url,
          fmt=None, assembly_report=None, fai=None, assembly_unit=None,
          mapping_style='identifiers'):
    annotation, output = Path(annotation), Path(output)
    if output.exists():
        raise FileExistsError(f'Reference destination already exists: {output}')
    fmt = fmt or ('gtf' if '.gtf' in annotation.name.lower() else 'gff3')
    if mapping_style not in {'identifiers', 'gene-symbol-rows'}:
        raise ValueError('Unknown mapping style')
    if mapping_style == 'gene-symbol-rows' and fmt != 'gtf':
        raise ValueError('gene-symbol-rows requires GTF gene records')
    symbol_rows = []
    allowed = None
    if assembly_unit:
        if not assembly_report:
            raise ValueError('Assembly-unit selection requires an NCBI assembly report')
        allowed = set()
        units = {unit.strip() for unit in assembly_unit.split(',')}
        found_units = set()
        with text_open(assembly_report) as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():
                    continue
                row = line.rstrip().split('\t')
                if len(row) >= 10 and row[7] in units:
                    found_units.add(row[7])
                    allowed.update(row[i] for i in (0, 4, 6, 9) if row[i] != 'na')
        if units != found_units:
            raise ValueError(f'Assembly units not found: {units - found_units}')
    genes, contigs, lengths, aliases = {}, set(), {}, {}
    def alias(name, canonical):
        if name in aliases and aliases[name] != canonical:
            raise ValueError(f'Ambiguous contig alias: {name}')
        aliases[name] = canonical
    headers = []
    records = 0
    with text_open(annotation) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                if len(headers) < 100:
                    headers.append(line.rstrip())
                if line.startswith('##sequence-region '):
                    _, chrom, start, end = line.split()
                    if allowed is not None and unquote(chrom) not in allowed:
                        continue
                    contigs.add(unquote(chrom))
                    lengths[unquote(chrom)] = int(end)
                continue
            if not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'Annotation line {number}: expected 9 columns')
            chrom = unquote(fields[0]) if fmt == 'gff3' else fields[0]
            if allowed is not None and chrom not in allowed:
                continue
            contigs.add(chrom)
            attrs = attributes(fields[8], fmt)
            if mapping_style == 'gene-symbol-rows':
                if fields[2] != 'gene':
                    continue
                if not attrs.get('gene_name'):
                    raise ValueError(f'GTF gene record lacks gene_name on line {number}')
                symbol_rows.append((chrom, attrs['gene_name'][0]))
            if fmt == 'gff3' and fields[2] not in {'gene', 'pseudogene'}:
                continue
            keys = ['gene_id', 'gene_name'] if fmt == 'gtf' else ['ID', 'Name', 'gene', 'gene_id', 'locus_tag', 'gene_synonym', 'Alias']
            if mapping_style == 'gene-symbol-rows':
                keys = ['gene_name']
            names = {n for key in keys for n in attrs.get(key, []) if n}
            if not names:
                if fmt == 'gtf':
                    continue
                raise ValueError(f'Gene record without an identifier on line {number}')
            records += 1
            for name in names:
                if any(c in name for c in '\t\r\n'):
                    raise ValueError(f'Unsupported control character in gene name: {name!r}')
                genes.setdefault(name, set()).add(chrom)
    if not genes:
        raise ValueError('No gene identifiers found; GFF3 needs gene/pseudogene records')
    for chrom in contigs:
        alias(chrom, chrom)
    sources = {'annotation': {'sha256': sha256(annotation), 'url': source_url, 'filename': annotation.name}}
    if assembly_report:
        sources['assembly_report'] = {'sha256': sha256(assembly_report), 'filename': Path(assembly_report).name}
        with text_open(assembly_report) as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():
                    continue
                row = line.rstrip().split('\t')
                if len(row) < 10:
                    raise ValueError('Malformed NCBI assembly report')
                if assembly_unit and row[7] not in units:
                    continue
                names = [row[i] for i in (0, 4, 6, 9) if row[i] != 'na']
                hits = set(names) & contigs
                if len(hits) > 1:
                    raise ValueError(f'Annotation contains multiple aliases of one contig: {hits}')
                canonical = next(iter(hits)) if hits else (row[6] if row[6] != 'na' else row[4])
                contigs.add(canonical)
                length = int(row[8])
                if canonical in lengths and lengths[canonical] != length:
                    raise ValueError(f'Contig length mismatch: {canonical}')
                lengths[canonical] = length
                for name in names:
                    alias(name, canonical)
    if fai:
        sources['fai'] = {'sha256': sha256(fai), 'filename': Path(fai).name}
        with open(fai, encoding='utf-8') as handle:
            for line in handle:
                name, length, *_ = line.split('\t')
                canonical = aliases.get(name, name)
                length = int(length)
                if canonical in lengths and lengths[canonical] != length:
                    raise ValueError(f'FAI length mismatch: {name}')
                contigs.add(canonical)
                lengths[canonical] = length
                alias(name, canonical)
    output.mkdir(parents=True)
    with (output/'gene_chr.txt').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        if mapping_style == 'gene-symbol-rows':
            writer.writerows(symbol_rows)
        else:
            for name in sorted(genes):
                for chrom in sorted(genes[name]):
                    writer.writerow([chrom, name])
    (output/'contigs.json').write_text(json.dumps({'aliases': aliases, 'lengths': lengths}, sort_keys=True), encoding='utf-8')
    metadata = {'schema_version': 1, 'species': species, 'assembly': assembly,
                'assembly_unit': assembly_unit or 'all',
                'annotation_release': release, 'format': fmt, 'sources': sources,
                'mapping_style': mapping_style,
                'annotation_headers': headers, 'gene_records_read': records,
                'gene_identifiers': len(genes), 'contigs': len(contigs),
                'multi_contig_identifiers': sum(len(v) > 1 for v in genes.values()),
                'files': {f: sha256(output/f) for f in ['gene_chr.txt', 'contigs.json']}}
    (output/'reference.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    return metadata


class Reference:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.metadata = json.loads((self.path/'reference.json').read_text())
        if self.metadata.get('schema_version') != 1:
            raise ValueError('Unsupported reference schema')
        if not {'gene_chr.txt', 'contigs.json'} <= self.metadata['files'].keys():
            raise ValueError('Reference is missing required files')
        for name in self.metadata['files']:
            if Path(name).name != name:
                raise ValueError('Reference file names must not contain directory components')
            if sha256(self.path/name) != self.metadata['files'][name]:
                raise ValueError(f'Reference checksum mismatch: {name}')
        self.fingerprint = sha256(self.path/'reference.json')
        self.genes = {}
        self.last_contig = {}
        with (self.path/'gene_chr.txt').open(newline='', encoding='utf-8') as handle:
            for chrom, name in csv.reader(handle, delimiter='\t'):
                self.genes.setdefault(name, set()).add(chrom)
                self.last_contig[name] = chrom
        data = json.loads((self.path/'contigs.json').read_text())
        self.aliases, self.lengths = data['aliases'], data['lengths']


def cis_mapping_policy(reference):
    expected='23050441eba9dd4a9e5632e4062926dd77fabe321a9fcdbb1b6f1bcade22bab3'
    verified_legacy = bool(
        reference and reference.metadata.get('species') == 'Homo sapiens'
        and reference.metadata.get('legacy_mapping_policy_supported')
        and reference.metadata.get('files', {}).get('gene_chr.txt') == expected
    )
    policy=os.getenv('TENET_CIS_MAPPING_POLICY', 'human_legacy' if verified_legacy else 'legacy' if reference else 'strict')
    if policy not in {'strict','human_legacy','legacy'}:
        raise ValueError('TENET_CIS_MAPPING_POLICY must be strict, legacy or human_legacy')
    if policy in {'legacy', 'human_legacy'} and not reference:
        raise ValueError('Legacy cis mapping requires a reference bundle; leave the policy unset for RNA-only modes')
    if policy=='human_legacy':
        if not reference or reference.metadata.get('species')!='Homo sapiens' or not reference.metadata.get('legacy_mapping_policy_supported'):
            raise ValueError('human_legacy requires the byte-identical human GENCODE v44 legacy bundle')
        if reference.metadata['files']['gene_chr.txt']!=expected:
            raise ValueError('human_legacy requires the verified original mapping checksum')
    return policy


def selected_reference(required=False, species=None):
    from code.species_support import species_record, canonical_species
    value = os.getenv('TENET_REFERENCE', '')
    if not value and required and species:
        value = species_record(species).get('reference', '')
    if not value:
        if required:
            raise ValueError(f'No reference for species {species!r}. Set TENET_REFERENCE to a bundle containing gene_chr.txt, contigs.json and reference.json. Build it from your GTF/GFF3 with python -m code.reference_bundle; see docs/OTHER_SPECIES.md.')
        return None
    path = Path(value).expanduser()
    if not path.is_dir():
        path = Path(__file__).resolve().parents[1]/'references'/value
    missing = [name for name in ['gene_chr.txt', 'contigs.json', 'reference.json'] if not (path/name).is_file()]
    if missing:
        raise ValueError(f'Reference files missing at {path}: {", ".join(missing)}. Supply a matching bundle via TENET_REFERENCE; see docs/OTHER_SPECIES.md.')
    ref = Reference(path)
    if species and canonical_species(species).casefold() != ref.metadata['species'].casefold():
        raise ValueError(f'Species/reference mismatch: {species} vs {ref.metadata["species"]}')
    return ref


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotation', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--species', required=True)
    parser.add_argument('--assembly', required=True)
    parser.add_argument('--release', required=True)
    parser.add_argument('--source-url', required=True)
    parser.add_argument('--format', choices=['gtf', 'gff3'])
    parser.add_argument('--assembly-report')
    parser.add_argument('--assembly-unit', help='Comma-separated exact NCBI assembly units; default includes all units')
    parser.add_argument('--fai')
    parser.add_argument('--mapping-style', choices=['identifiers', 'gene-symbol-rows'], default='identifiers',
                        help='gene-symbol-rows preserves GTF gene_name row order and duplicates for legacy file reproduction')
    args = parser.parse_args()
    print(json.dumps(build(args.annotation, args.output, args.species, args.assembly,
                           args.release, args.source_url, args.format,
                           args.assembly_report, args.fai, args.assembly_unit,
                           args.mapping_style), indent=2))


if __name__ == '__main__':
    main()
