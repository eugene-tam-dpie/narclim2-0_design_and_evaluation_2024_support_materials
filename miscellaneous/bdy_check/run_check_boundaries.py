#!/usr/bin/env python3
from pathlib import Path
import sys
from string import Template
import argparse
import json
import os

import check_boundaries as cb

runpath=Path('.').resolve()
script_path=Path(__file__).resolve().parent
file_lists_path=runpath.joinpath('file_lists')
split_paths=file_lists_path.joinpath('chunked')
templates_path=script_path.joinpath(
    'monitor',
    'file_checks',
    'templates'
    )
pbs_path=runpath.joinpath('pbs')

num_chunks=64

def replace_control_strings(job_dict,tup):
    key,val=tup
    val0=val
    if type(val)==str:
        if '$' in val:
            template=Template(val)
            val0=template.substitute(os.environ[val])
        if '&' in val:
            val0=val0.replace('&','$')
            val0=Template(val0)
            val0=val0.substitute(**job_dict)
    return key,val0

def parse_json(json_file):
    with open(json_file,'rt') as f:
        job_dict=json.load(f)
    jd_copy=job_dict.copy()
    job_dict=map(
        lambda x: replace_control_strings(jd_copy,x),
        job_dict.items()
    )
    return dict(job_dict)

def parse_args():
    parser = argparse.ArgumentParser(
        description = (
            'Reads inputs and extracts grid information'
            )
    )
    parser.add_argument(
        '--definition_json',
        dest = 'definition_json',
        default = runpath.joinpath('inputs.json').as_posix(),
        type=Path,
        help = (
            'Job definition file.'
            )
        )
    parser.add_argument(
        '--run_mode',
        dest = 'run_mode',
        default = 0,
        help = (
            'Set for dry run. Effects depends on script'
            )
        )
    args = parser.parse_args()
    job_json=args.definition_json
    job_dict=parse_json(job_json)
    job_dict['run_mode']=args.run_mode
    return job_dict

def check_missing(file_df):
    missing=file_df[file_df['exists']==False]
    pbs=missing['pbs_script']
    pbs=set(pbs)
    pbs=list(pbs)
    pbs.sort()
    print('\n'.join(pbs))

def pre_pbs():
    split_paths.mkdir(parents=True,exist_ok=True)
    pbs_path.mkdir(parents=True,exist_ok=True)
    outnames=cb.find_shuffle_chunk(num_chunks,split_paths)
    with templates_path.joinpath('pbs_template.pbs').open('rt') as f:
        pbs_template=f.read().strip()
    pbs_template=Template(pbs_template)
    for zfill, split_names in outnames:
        job_string=f'{zfill} {split_names}'
        outstring=pbs_template.substitute(job_num=job_string).strip()
        with pbs_path.joinpath(f'run_{zfill}.pbs').open('wt') as f:
            f.write(outstring)

def pbs_tasks(flist):
    with flist.open('rt') as f:
        all_files=f.read().strip()
    all_files=all_files.split('\n')
    var_vals=cb.check_vars_in_files(all_files)
    return var_vals

def main():
    zfill=sys.argv[1]
    job_id=int(zfill)
    if job_id<0:
        pre_pbs()
    else:
        parquet_path=runpath.parent.joinpath('parquet')
        parquet_path.mkdir(parents=True,exist_ok=True)
        flist=Path(sys.argv[2])
        var_vals=pbs_tasks(flist)
        fout=parquet_path.joinpath(f'variable_checks_{zfill}.parquet')
        var_vals.to_parquet(fout,engine='pyarrow')

if __name__=='__main__':
    main()
