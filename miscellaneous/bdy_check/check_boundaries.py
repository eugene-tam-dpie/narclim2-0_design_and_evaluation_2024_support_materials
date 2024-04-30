from pathlib import Path
import datetime as dt
import pandas as pd
from collections import namedtuple
import os
import itertools
import xarray as xr

import common_functions as cf

project_base=Path("base-path-of-your-project")

bdy_prefixes=(
    'bdy',
    'input',
    'lowinp'
)
single_domain='bdy'

MissingFile=namedtuple(
    'MissingFile',
    (
        'file',
        'pbs_script',
        'namelist',
        'exists',
        'size'
    )
)

ignore_zeros=(
    'LAI',
    'SST',
    'VAR',
    'VAR_SSO',
    'VEGFRAC',
)

try:
    num_cpus=int(os.environ['PBS_NCPUS'])
except:
    num_cpus=1

def name_to_date(f):
    stem=f.stem
    datestring=stem.split('_')
    datestring=datestring[-2:]
    datestring='_'.join(datestring)
    return dt.datetime.strptime(datestring,'%Y-%m-%d_%H:%M:%S')

def check_pbs(leaf_path):
    pbs_files=set(leaf_path.glob('*.pbs'))
    output_files=map(
        lambda x: leaf_path.joinpath(x.stem),
        leaf_path.glob('*.o*')
    )
    output_files=set(output_files)
    missing_files=pbs_files.difference(output_files)
    missing_files=list(missing_files)
    missing_files.sort()
    return missing_files

def check_existance(pobj,pattern):
    files=pobj.glob(pattern)
    files=tuple(files)
    lf=len(files)
    assert lf<=1, f'Pattern too general - {pobj.as_posix()}/{pattern} - {files}'
    if lf==1:
        exists=True
    else:
        exists=False
    return exists

def make_missing_tup(test_file,file_part,num):
    exists=test_file.is_file()
    if exists:
        size=test_file.stat().st_size
    else:
        size=0
    leaf_path=test_file.parent
    year=leaf_path.name
    numz=str(num).zfill(3)
    wps_pbs_equiv=leaf_path.as_posix().replace('bdy','pbs')
    wps_pbs_equiv=f'{wps_pbs_equiv}/{file_part}_Real_{year}_{numz}.pbs'
    nml_equiv=leaf_path.as_posix().replace('bdy','namelists_wps')
    nml_equiv=f'{nml_equiv}/namelist_{year}_{numz}.wrf'
    out_tup=MissingFile(
        test_file.as_posix(),
        wps_pbs_equiv,
        nml_equiv,
        exists,
        size
    )
    return out_tup

def check_bdy(lp):
    leaf_path=Path(lp)
    year=leaf_path.name
    dates=cf.split_dates(60,int(year))
    parts=leaf_path.parts
    # gcm,ens_mem,rcm,exp=parts[0:4]
    file_part='_'.join(parts[-6:-2])
    out_tup=tuple()
    for num, date in enumerate(dates):
        string_date=date.strftime('%Y-%m-%d_%H:%M:%S')
        for p in bdy_prefixes:
            prefix_base=f'wrf{p}'
            if p==single_domain:
                glob_pattern=f'{prefix_base}_d01_{string_date}'
                test_file=leaf_path.joinpath(glob_pattern)
                out_tup+=make_missing_tup(test_file,file_part,num),
            else:
                for dom in range(2):
                    dom_string=str(dom+1).zfill(2)
                    glob_pattern=f'{prefix_base}_d{dom_string}_{string_date}'
                    test_file=leaf_path.joinpath(glob_pattern)
                    out_tup+=make_missing_tup(test_file,file_part,num),
    return out_tup

def calc_diffs(da):
    dims=da.dims
    diff_dict={}
    for dim in dims:
        diff=da.diff(dim)
        try:
            diff_dict[f'{dim}_min']=float(diff.min())
            diff_dict[f'{dim}_max']=float(diff.max())
        except:
            pass
    return diff_dict

def check_ds(ds,fname):
    out_tup=tuple()
    for var in ds.data_vars:
        if (
            str(ds[var].dtype).startswith('float') or
            str(ds[var].dtype).startswith('int')
        ):
            da=ds[var]
            if var in ignore_zeros:
                if float(da.min())!=float(da.max()):
                    da=da.where(da!=0)
            dom=os.path.basename(fname)
            dom=dom.split('_')
            dom=dom[1]
            key=f'{var}_{dom}'
            main_dict={
                'file':fname,
                'variable':var,
                'var_key':key,
                'min':float(da.min()),
                'mean':float(da.mean()),
                'max':float(da.max()),
                'nan_count':float((~da.notnull()).astype(float).mean())
            }
            diff_dict=calc_diffs(da)
            out_tup+={**main_dict,**diff_dict},
    return out_tup

def check_vars(fname):
    ds=xr.open_dataset(fname,engine='netcdf4')
    tups=check_ds(ds,fname)
    return tups

def find_all_leaf_bdy():
    all_bdy_base=[
        'find',
        project_base.as_posix(),
        '-maxdepth',
        '5',
        '-name',
        'bdy'
    ]
    all_bdy_base=cf.sub_no_dir(all_bdy_base)
    leaf_dirs=[
        'xargs',
        '-I',
        '%',
        'find',
        '%',
        '-type',
        'd',
        '-links',
        '2'
    ]
    leaf_dirs=cf.sub_pipe_stdin(leaf_dirs,all_bdy_base)
    leaf_dirs=cf.stdout_to_list(leaf_dirs)
    return leaf_dirs

def leaf_dir_to_all_files(ld):
    cmd=[
        'find',
        ld,
        '-type',
        'f'
    ]
    all_f=cf.sub_no_dir(cmd)
    all_f=cf.stdout_to_list(all_f)
    return all_f

def find_missing_files(leaf_dirs):
    if num_cpus>1:
        output_tups=cf.pool_function(check_bdy,leaf_dirs,num_cpus)
        mising_files=tuple(itertools.chain.from_iterable(output_tups))
    else:
        mising_files=tuple()
        for ld in leaf_dirs:
            mising_files+=check_bdy(ld)
    mising_files=pd.DataFrame(mising_files)
    return mising_files

def check_vars_in_files(all_files):
    if num_cpus>1:
        output_tups=cf.pool_function(check_vars,all_files,num_cpus)
        var_vals=tuple(itertools.chain.from_iterable(output_tups))
        print(len(var_vals))
    else:
        for f in all_files:
            var_vals+=check_vars(f)
    var_vals=pd.DataFrame(var_vals)
    return var_vals

def check_all_bdy():
    leaf_dirs=find_all_leaf_bdy()
    mising_files=find_missing_files(leaf_dirs)
    if num_cpus>1:
        output_tups=cf.pool_function(leaf_dir_to_all_files,leaf_dirs,num_cpus)
        all_files=tuple(itertools.chain.from_iterable(output_tups))
        print(len(all_files))
    else:
        for ld in leaf_dirs:
            all_files+=tuple(leaf_dir_to_all_files(ld))
    var_vals=check_vars_in_files(all_files)
    mising_files=pd.DataFrame(mising_files)
    var_vals=pd.DataFrame(var_vals)
    return mising_files,var_vals

def find_shuffle_chunk(chunks,chunk_output_dir):
    all_bdy_base=[
        'find',
        project_base.as_posix(),
        '-maxdepth',
        '5',
        '-name',
        'bdy'
    ]
    all_bdy_base=cf.sub_no_dir(all_bdy_base)
    all_files=[
        'xargs',
        '-I',
        '%',
        'find',
        '%',
        '-type',
        'f',
        '-name',
        'wrf*'
    ]
    all_files=cf.sub_pipe_stdin(all_files,all_bdy_base)
    all_files=cf.stdout_to_list(all_files)
    fnames=tuple()
    for i in range(chunks):
        zfill=str(i).zfill(2)
        tmp=all_files[i::chunks]
        outstring='\n'.join(tmp)
        outfile=chunk_output_dir.joinpath(f'file_split.{zfill}')
        fnames+=(zfill,outfile),
        with outfile.open('wt') as f:
            f.write(outstring)
    return fnames