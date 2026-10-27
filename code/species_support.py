"""Resolve species presets and check required resources before computation."""
import argparse
import json
import os
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def species_record(species):
    catalog=json.loads((ROOT/'references/catalog.json').read_text())
    value=str(species).casefold().replace('_',' ').strip()
    for key,record in catalog['species'].items():
        names=[key,record['scientific_name']]+record.get('aliases',[])
        if value in [n.casefold().replace('_',' ') for n in names]:
            return dict(record,key=key)
    return {}

def canonical_species(species):
    return species_record(species).get('scientific_name',species)

def needs_tf(mode,pair_mode=None):
    pair_mode=(pair_mode or os.getenv('TENET_PAIR_MODE','default')).strip().lower()
    return pair_mode=='default' and str(mode) in {'0','1','2','3','4'}

def tf_list_path(species):
    explicit=os.getenv('TENET_TF_LIST')
    if explicit:
        path=Path(explicit).expanduser().resolve()
    else:
        record=species_record(species)
        path=ROOT/'input'/record.get('tf_file',f'GO_symbol_{record.get("key",species)}_regulation_of_transcription+sequence-specific_DNA_binding_list.txt')
    if not path.is_file():
        raise ValueError(f'TF list missing for {species!r}: {path}. TF modes 0-4 require a species-specific list, one matrix-matching RNA ID per line. Set TENET_TF_LIST=/absolute/path/to/TFs.txt. Modes 5/6 and all-pair modes do not require it. See docs/OTHER_SPECIES.md.')
    return path

def load_tf_names(species):
    path=tf_list_path(species)
    names=path.read_text(encoding='utf-8-sig').splitlines()
    if not names or any(not n or '\t' in n or not n.strip() for n in names):
        raise ValueError(f'TF list must contain one nonempty ID per line, without a header or tabs: {path}')
    if len(names)!=len(set(names)):
        raise ValueError(f'TF list contains duplicate IDs: {path}')
    return names

def preflight(species,mode,matrix=None,pair_mode=None):
    from code.reference_bundle import selected_reference,sha256,cis_mapping_policy
    reference=selected_reference(required=str(mode)!='0',species=species)
    policy=cis_mapping_policy(reference)
    result={'species':canonical_species(species),'mode':str(mode),
            'reference':str(reference.path) if reference else None,
            'reference_metadata':reference.metadata if reference else None,
            'tf_required':needs_tf(mode,pair_mode),'tf_list':None,'cis_mapping_policy':policy}
    if result['tf_required']:
        names=load_tf_names(species);path=tf_list_path(species)
        result.update(tf_list=str(path),tf_list_sha256=sha256(path),TF_identifiers=len(names))
        if matrix:
            from code.feature_utils import matrix_header
            present=set(matrix_header(matrix)) & set(names)
            if not present:
                raise ValueError('No TF list identifiers match the input matrix. Use the same species and naming convention (symbols vs stable IDs). See docs/OTHER_SPECIES.md.')
            result['TF_identifiers_in_matrix']=len(present)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species',required=True)
    parser.add_argument('--mode',required=True,choices=[str(i) for i in range(10)])
    parser.add_argument('--matrix')
    parser.add_argument('--output')
    args=parser.parse_args()
    mode=args.mode;pair_mode=None
    if mode in {'7','8','9'}:
        pair_mode='all_feature' if mode=='9' else 'gene_only'
        mode='0' if mode=='7' else '1'
    try:
        result=preflight(args.species,mode,args.matrix,pair_mode)
        if args.output:Path(args.output).write_text(json.dumps(result,indent=2))
        display={k:v for k,v in result.items() if k!='reference_metadata'}
        if result['reference_metadata']:
            display.update(assembly=result['reference_metadata']['assembly'],annotation=result['reference_metadata']['annotation_release'])
        print('[Species preflight] '+json.dumps(display,ensure_ascii=False))
    except (ValueError,FileNotFoundError,KeyError) as error:
        print(f'[Species preflight ERROR] {error}',file=sys.stderr)
        raise SystemExit(2)

if __name__=='__main__':main()
