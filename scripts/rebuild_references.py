"""Rebuild pinned bundled references from checksum-verified official sources.

Example: python scripts/rebuild_references.py --species human --cache /data/sources
         --output /data/rebuilt_references
"""
import argparse
import json
import sys
import urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from code.reference_bundle import build,sha256

def fetch(source,cache):
    path=cache/source['filename']
    if not path.exists():
        req=urllib.request.Request(source['url'],headers={'User-Agent':'TENETPLUS-reference-rebuild/1.0'})
        with urllib.request.urlopen(req,timeout=180) as response,path.with_suffix(path.suffix+'.part').open('wb') as out:
            while chunk:=response.read(1024*1024):out.write(chunk)
        path.with_suffix(path.suffix+'.part').replace(path)
    if sha256(path)!=source['sha256']:
        raise ValueError(f'Pinned source checksum differs: {path}; do not silently substitute another release')
    return path

def main():
    specs=json.loads((ROOT/'references/build_sources.json').read_text())['species']
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species',choices=list(specs)+['all'],required=True)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.cache.mkdir(parents=True,exist_ok=True)
    selected=list(specs) if args.species=='all' else [args.species]
    for key in selected:
        item=specs[key];destination=args.output/item['bundle']
        if destination.exists():raise FileExistsError(f'Output already exists: {destination}')
        files={role:fetch(source,args.cache) for role,source in item['sources'].items()}
        metadata=build(files['annotation'],destination,item['species'],item['assembly'],item['release'],
                       item['sources']['annotation']['url'],fmt=item['format'],
                       assembly_report=files.get('assembly_report'),mapping_style=item['mapping_style'])
        if key=='human':
            expected='23050441eba9dd4a9e5632e4062926dd77fabe321a9fcdbb1b6f1bcade22bab3'
            if sha256(destination/'gene_chr.txt')!=expected:raise ValueError('Human mapping does not reproduce the original file')
            metadata['legacy_mapping_policy_supported']=True
            metadata['legacy_runtime_policy']='human_legacy is the default for this verified bundle; TENET_CIS_MAPPING_POLICY=strict opts into strict validation'
            (destination/'reference.json').write_text(json.dumps(metadata,indent=2))
        print(f'{key}: {destination}',flush=True)

if __name__=='__main__':main()
