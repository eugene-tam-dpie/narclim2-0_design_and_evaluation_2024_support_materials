import logging
from logging.handlers import RotatingFileHandler
import functools
import json
import subprocess as sp
import multiprocessing as mp
import cftime as cft
from string import Template
from pathlib import Path
import time

#COMMON FUNCTIONS
#Setup logger (functionally the same for every script, differing only in log file name)
def logger_setup(log_name,log_id):
    Path(log_name).resolve().parent.mkdir(exist_ok=True,parents=True)
    log_file = RotatingFileHandler(log_name)
    log_file.setLevel(logging.DEBUG)
    logging.captureWarnings(True)
    logger = logging.getLogger(log_id)
    warnings_logger = logging.getLogger("py.warnings")
    logger.setLevel(logging.DEBUG)
    # create file handler
    log_cons = logging.StreamHandler()
    log_cons.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    log_form = logging.Formatter('%(asctime)s | %(filename)s:%(lineno)s | %(levelname)s | %(message)s')
    log_cons.setFormatter(log_form)
    log_file.setFormatter(log_form)
    # add the handlers to the logger
    logger.addHandler(log_file)
    logger.addHandler(log_cons)
    warnings_logger.addHandler(log_file)
    warnings_logger.addHandler(log_cons)
    return logger

# Run function with potential errors
def try_except_wrapper(log_id):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger=logging.getLogger(log_id)
            try:
                return func(*args, **kwargs)
            except Exception as e:
                logger.exception(f'Error raise in {func.__name__}')
            # re-raise the exception
                raise SystemExit(e)
        return wrapper
    return decorator

def type_to_namelist_type(string):
    if type(string)==str:
        return f"'{string}'"
    elif type(string)==bool:
        return f'.{string}.'
    else:
        return str(string)

def duplicate(string,num):
    out=map(
        lambda x: type_to_namelist_type(string),
        range(num)
    )
    out=tuple(out)
    out=', '.join(out)
    return out

def tuple_to_namelist_string(tup,num):
    assert len(tup)==num, f'Input tuple does not match number of domains {tup} : {num}'
    string_tup=map(
        lambda x: type_to_namelist_type(x),
        tup
    )
    string=', '.join(string_tup)
    return string

def base_ensemble_output_path(base,ensemble_key):
    ensemble=ensemble_key.split('/')
    return base.joinpath(*ensemble)

def dict_pretty_print(in_dict):
    string=json.dumps(in_dict,indent=4)
    print(string)

def parts_split_year(days_in_year,parts):
    days_per_part=days_in_year/parts
    split_days=map(
        lambda x: int(days_per_part*x),
        range(parts)
    )
    return split_days

def cftime_relativetime(
    reference,
    delta,
    units,
    calendar
):
    iso=reference.isoformat()
    timeunits=f'{units} since {iso}'
    num_ref=cft.date2num(reference,timeunits,calendar)
    out_num=num_ref+delta
    return cft.num2date(out_num,timeunits,calendar)

def split_dates(parts,year,calendar):
    year_start=cft.datetime(year,1,1,0,0,calendar=calendar)
    year_end=cft.datetime(year+1,1,1,0,0,calendar=calendar)
    days_in_year=year_end-year_start
    days_in_year=days_in_year.days
    split_days=tuple(parts_split_year(days_in_year,parts))
    split_days=map(
        lambda x: cftime_relativetime(
            year_start,
            x,
            'days',
            calendar
        ),
        split_days
    )
    return tuple(split_days)

def parts_split_calendar(
    ensemble_start,
    ensemble_end,
    parts,
    calendar
):
    splits={}
    for year in range(ensemble_start, ensemble_end+1):
        split_days=split_dates(parts,year,calendar)
        splits[year]=split_days
    return splits

#########################################################
# File and bash routines
def sub_with_dir(cmd, wdir):
    shell = sp.Popen(cmd, cwd = wdir, stdout = sp.PIPE)
    return shell.stdout
    
def sub_no_dir(cmd):
    shell = sp.Popen(cmd, stdout = sp.PIPE)
    return shell.stdout

def sub_pipe_stdin(cmd,stdin):
    shell = sp.Popen(cmd, stdout = sp.PIPE, stdin=stdin)
    return shell.stdout

def stdout_to_utf8(stdout):
    try:
        raw = stdout.read()
    except AttributeError:
        raw = stdout
    utf = raw.decode('utf-8')
    return utf.strip()

def string_to_list(string):
    string = string.strip()
    if string == "":
        return []
    return string.split('\n')

def stdout_to_list(stdout):
    foo = stdout_to_utf8(stdout)
    return string_to_list(foo)

def sub_with_retry(cmd,timeout,wdir=False):
    logger=logging.getLogger('__main__')
    not_finished=True
    retries=0
    max_retries=10
    max_sleep_retries=60
    seconds_per_sleep=1
    while not_finished and retries<max_retries:
        if wdir:
            shell = sp.Popen(cmd, cwd = wdir, stdout = sp.PIPE, stderr = sp.PIPE)
        else:
            shell = sp.Popen(cmd, stdout = sp.PIPE, stderr = sp.PIPE)
        try:
            stdout,stderr=shell.communicate(timeout=timeout)
            not_finished=False
        except:
            shell.kill()
            logger.warning(f'timeout error for {cmd}, retry {retries}')
            sleep_retry=0
            while shell.returncode is None:
                time.sleep(seconds_per_sleep)
                sleep_retry+=1
                if sleep_retry>=(max_sleep_retries-2):
                    shell.terminate()
                if sleep_retry>=max_sleep_retries:
                    logger.warning(f'unable to kill stuck process {cmd}')
                    break
            retries+=1
    assert retries<max_retries, f'Too many retries for subprocess {cmd}'
    stdout=stdout_to_utf8(stdout)
    stderr=stdout_to_utf8(stderr)
    return stdout,stderr,shell.returncode

#########################################################
# parallelization
def pool_function(func, in_list, num_procs):
    tlen=len(in_list)
    chunksize=int(tlen/(num_procs*256))
    chunksize=max(chunksize, 1)
    p = mp.Pool(processes = num_procs)
    out = p.imap_unordered(func, in_list, chunksize)
    return out

#########################################################
# String template
def create_template_from_file(pobj):
    with pobj.open('rt') as f:
        string=f.read()
    return Template(string)

def string_to_file(pobj,string):
    pobj.parent.mkdir(parents=True,exist_ok=True)
    with pobj.open('wt') as f:
        f.write(string)

