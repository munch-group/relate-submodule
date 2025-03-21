

import os, re
from collections import defaultdict
from pathlib import Path
import pandas as pd

from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect

def modpath(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


def groupby_chrom(files, pattern='chr_?([XYxy\d+]+)'):
    groups = defaultdict(list)
    for f in files:
        chrom = re.search(pattern, os.path.splitext(os.path.basename(f))[0]).group(1)
        groups[chrom].append(f)
    return groups


def state_segments(posterior_file):

    stepsdir = 'steps/state_segments'
    # if not os.path.exists(stepsdir):
    #     os.makedirs(stepsdir)

    segment_file = modpath(posterior_file, parent=stepsdir, suffix='.h5')

    inputs = {'posterior_file': posterior_file}
    outputs = {'segment_file': segment_file}

    options = {'memory': '40g', 'walltime': '01:00:00'} 
    spec = f"""
    mkdir -p {stepsdir}    
    python scripts/state_segments.py {posterior_file} {segment_file}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def trio_segments(state_segment_files, trio_segment_file):

    # python scripts/trio_segments.py gene_trees.h5 steps/state_segments/*_chr_22.h5

    stepsdir = 'steps/trio_segments'
    # if not os.path.exists(stepsdir):
    #     os.makedirs(stepsdir)

    inputs = {'state_segment_files': state_segment_files}
    outputs = {'trio_segment_file': trio_segment_file}

    options = {'memory': '40g', 'walltime': '01:00:00'} 
    spec = f"""
    mkdir -p {stepsdir}    
    python scripts/trio_segments.py {trio_segment_file} {state_segment_files}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def ils_in_windows(segment_file):

    stepsdir = 'steps/ils_in_windows'
    # if not os.path.exists(stepsdir):
    #     os.makedirs(stepsdir)

    window_file = modpath(segment_file, parent=stepsdir, suffix='.h5')

    inputs = {'segment_file': segment_file}
    outputs = {'window_file': window_file}

    options = {'memory': '8g', 'walltime': '01:00:00'} 
    spec = f"""
    mkdir -p {stepsdir}
    python scripts/ils_in_windows.py {segment_file} {window_file}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def low_ils_regions(window_file):

    stepsdir = 'steps/low_ils_regions'
    # if not os.path.exists(stepsdir):
    #     os.makedirs(stepsdir)

    low_ils_file = modpath(window_file, parent=stepsdir, suffix='.csv')

    inputs = {'window_file': window_file}
    outputs = {'low_ils_file': low_ils_file}

    options = {'memory': '8g', 'walltime': '01:00:00'} 
    spec = f"""
    mkdir -p {stepsdir}
    python scripts/low_ils_regions.py {window_file} {low_ils_file}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def workflow(working_dir=os.getcwd(), defaults={}, input_files=[]):

    gwf = Workflow(working_dir=working_dir, defaults=defaults)

    # dict of targets as info for other workflows
    targets = defaultdict(list)

    # compute state segments
    targets_state_segments = gwf.map(state_segments, input_files)
    targets['segments'] = targets_state_segments

    # compute ils in windows
    targets_ils_in_windows = gwf.map(ils_in_windows, targets_state_segments.outputs)
    targets['windows'] = targets_ils_in_windows

    stepsdir = 'steps/merge_ils_data'
    # if not os.path.exists(stepsdir):
    #     os.makedirs(stepsdir)

    ils_window_files = collect(targets_ils_in_windows.outputs, ['window_file'])['window_files']
    merged_ils_file = os.path.join(stepsdir, 'merged_ils_data.h5')
    input_args = ' '.join(ils_window_files)
    target = gwf.target('merge_ils_data', memory='36g', walltime='01:00:00', inputs=ils_window_files, outputs=[merged_ils_file]) << f"""
    mkdir -p {stepsdir}
    python scripts/merge_hdf_files.py {input_args} {merged_ils_file}
    """
    targets['merge_ils_data'] = [merged_ils_file]

    # compute low ils regions
    targets_low_ils_regions = gwf.map(low_ils_regions, targets_ils_in_windows.outputs)
    targets['regions'] = targets_low_ils_regions

    stepsdir = 'steps/merge_low_data'
    # if not os.path.exists(stepsdir):
    #     os.makedirs(stepsdir)

    low_ils_region_files = collect(targets_low_ils_regions.outputs, ['low_ils_file'])['low_ils_files']

    merged_low_region_file = os.path.join(stepsdir, 'merged_low_ils_regions.csv')

    input_args = ' '.join(low_ils_region_files)
    target = gwf.target('merge_low_ils_regions', memory='36g', walltime='01:00:00', inputs=low_ils_region_files, outputs=[merged_low_region_file]) << f"""
    mkdir -p {stepsdir}
    python scripts/merge_csv_files.py {input_args} {merged_low_region_file}
    """
    targets['merge_regions'] = targets_low_ils_regions

    state_segment_files = collect(targets_state_segments.outputs, ['segment_file'])['segment_files']



    # for chrom, input_files in groupby_chrom(state_segment_files).items():
    #     output_file_name = f'XXXXXXXXXX'
    #     gwf.target_from_template(
    #     name=f'tree_segments_{chrom}',
    #     trio_segments(
    #         state_segment_files=input_files,
    #         trio_segment_file=output_file_name
    #     )
    # )

    return gwf, targets


####################################################################
# Use code like this to run this as standalone workflow: 
####################################################################

# data_dir = '/home/kmt/Primategenomes/data/final_tables'
# state_posterior_files = sorted(Path(data_dir).glob('**/*.HDF'))
# gwf, codeml_targets  = workflow(working_dir=working_dir, 
#                                     defaults={'account': 'xy-drive'},
#                                     input_files=state_posterior_files)

####################################################################
# Use code like this to run this as a submodule workflow: 
####################################################################

# data_dir = '/home/kmt/Primategenomes/data/final_tables'
# state_posterior_files = sorted(Path(data_dir).glob('**/*.HDF'))
# ils = importlib.import_module('primate-ils.workflow')
# gwf, codeml_targets  = ils.workflow(working_dir=working_dir, 
#                                     defaults={'account': 'xy-drive'},
#                                     input_files=state_posterior_files)
# globals()['ils'] = gwf
