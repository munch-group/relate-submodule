

# %% [markdown]
# ---
# title: GWF workflow
# execute:
#   eval: false
# ---

# %%

import os
os.environ['NUMEXPR_MAX_THREADS'] = '16'

from gwf import Workflow
import re
from collections import defaultdict
from pathlib import Path
import pandas as pd
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect
import os, re
from collections import defaultdict



# function that modifies file path
def modify_path(p, parent=None, base=None, suffix=None):
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


# function to combine 2 target outputs as an input
def combine(*args, only=None):
    assert all(len(args[0]) == len(args[i]) for i in range(len(args)))
    combined = []
    for j in range(len(args[0])):
        output_group = {}
        for i in range(len(args)):
            if only:
                output_group.update({k: v for k, v in args[j].items() if k in only})
            else:
                output_group.update(args[i][j])
        combined.append(output_group)
    return combined


def download_data(config):
    inputs = []
    outputs = {'ancestral_vcf': f"steps/data/{os.path.basename(config['ancestral_vcf'])}",
               'sample_vcf': f"steps/data/{os.path.basename(config['sample_vcf'])}",
               'sample_vcf_index': f"steps/data/{os.path.basename(config['sample_vcf_index'])}",
               '1000G_2504_seq_index': f"steps/data/{os.path.basename(config['1000G_2504_seq_index'])}",
               '1000G_698_seq_index': f"steps/data/{os.path.basename(config['1000G_698_seq_index'])}",
               'mask': f"steps/data/{os.path.basename(config['mask']).replace('.gz', '')}", # remove gz suffix because it is unpacked
               'ancestral_fa': 'steps/data/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_X.fa',
    }
    options = {'memory': '8g', 'walltime': '10:00:00'}
    spec = f'''
    mkdir -p steps/data
    wget --directory-prefix steps/data {config['sample_vcf']}
    wget --directory-prefix steps/data {config['sample_vcf_index']}
    wget --directory-prefix steps/data {config['1000G_2504_seq_index']}
    wget --directory-prefix steps/data {config['1000G_698_seq_index']}

    wget --directory-prefix steps/data {config['mask']}
    wget --directory-prefix steps/data {config['ancestral_vcf']}
    cd steps/data/
    gzip -d {os.path.basename(config['mask'])}
    tar xfvz {os.path.basename(config['ancestral_vcf'])}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# map of recombination rate across the X chromosome made by DECODE genetics
def decode_genetic_maps(decode_hg38_sexavg_per_gen, genetic_map_chrX):
    inputs = [decode_hg38_sexavg_per_gen]
    outputs = [genetic_map_chrX]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    mkdir -p {os.path.dirname(genetic_map_chrX)}
    cat {decode_hg38_sexavg_per_gen} | tail -n +2 | grep chrX | cut -f 2,4,5 | (echo pos COMBINED_rate Genetic_Map ; cat - ; ) > {genetic_map_chrX}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# turn diploid females (XX) into two individual haplotypes (haploid individuals) like males
def female_haploid(chrX_filtered_eagle2_phased):
    phased_haplotypes = 'steps/1000g_phased_haplotypes.vcf.gz'
    inputs = [chrX_filtered_eagle2_phased]
    outputs = {'haplotypes': phased_haplotypes}
    options = {'memory': '10g', 'walltime': '01:20:00'}
    spec = f'''
    python scripts/haploid_vcf.py {chrX_filtered_eagle2_phased} | gzip > {phased_haplotypes}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct files with haplotype IDs
def haplotype_id(phased_haplotypes):
    phased_haplotypes_id = 'steps/1000g_phased_haplotypes_ids.txt'
    inputs = [phased_haplotypes]
    outputs = {'ids': phased_haplotypes_id}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    # conda install -c bioconda bcftools
    # conda install openssl   ## to install libcrypto.so.1.0.0 library
    bcftools query -l {phased_haplotypes} > {phased_haplotypes_id}
    sleep 5
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct populations labels mapping each haplotype to a population
# (group haplotypes according to the population to which the individuals carrying those haplotypes belong)
def all_pop_labels(phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index):
    phased_haplotypes_poplabels = 'steps/1000g_phased_haplotypes_poplabels.txt'
    inputs = [phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index]
    outputs = {'poplabels': phased_haplotypes_poplabels}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    python scripts/make_poplabels.py {phased_haplotypes_id} {high_coverage_seq_index} {related_high_coverage_seq_index} > {phased_haplotypes_poplabels} 
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Define the function to convert VCF to haps/sample format
def convert_vcf(phased_haplotypes, phased_haplotypes_poplabels):
    phased_haplotypes_haps = modify_path(phased_haplotypes, suffix='.haps')
    phased_haplotypes_sample = modify_path(phased_haplotypes, suffix='.sample')
    inputs = [phased_haplotypes_poplabels, phased_haplotypes]
    outputs = {'haps': phased_haplotypes_haps, 'sample': phased_haplotypes_sample}
    options = {'memory': '10g', 'walltime': '01:00:00'}
    spec = f'''
    {config['relate_dist_dir']}/bin/RelateFileFormats --mode ConvertFromVcf --haps {phased_haplotypes_haps} --sample {phased_haplotypes_sample} -i {phased_haplotypes.replace('.vcf.gz', '')} --poplabels {phased_haplotypes_poplabels}
    sleep 20
    touch {phased_haplotypes_haps}
    touch {phased_haplotypes_sample}
    '''
    # Returning outputs as well
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


## start with specific population ##

# exclude related individuals to avoid biases arising from shared genetic material
def exclude_related(path, population):
    output_dir = f'steps/{population}/excluded'
    output_path = modify_path(path, parent=output_dir, suffix='_related.txt')
    inputs = {'path' : path}
    outputs = {'path' : output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v '#' {path} | cut -f 10 > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# find IDs of haplotypes from all other populations so we can exclude them
def ids_other_ppl(path, population):
    output_dir = f'steps/{population}/excluded'
    output_path = modify_path(path, parent=output_dir, suffix='_non_ppl.txt')
    inputs = {'path' : path}
    outputs = {'path' : output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v {population} {path} | cut -f 1 -d ' ' > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# combine excluded files: both related and non population individuals
def combine_files(path, population=None, related=None):
    output_dir = f'steps/{population}/combined'
    output_path = modify_path(path, parent=output_dir, base='', suffix='excluded_combined.txt')
    # output_dir = modify_path(output_path, base='', suffix='')
    inputs = {'path': path, 'related': related}
    outputs = {'path': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    cat {path} {related} | sort | uniq > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct a list of excluded individuals
def excluded_list(path, population=None, haplotype_id=None):
    output_dir = f'steps/{population}/combined'
    output_path = modify_path(path, parent=output_dir, base='', suffix='excluded_list.txt')
    # output_dir = modify_path(output_path, base='', suffix='')
    inputs = {'path': path, 'haplotype_id': haplotype_id}
    outputs = {'exclude_list': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -f {path} {haplotype_id} > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct a list of only individuals from the population of interest
def pop_labels(exclude_list, population=None, poplabels=None):
    output_dir = f'steps/{population}/included'
    output_path = os.path.join(output_dir, 'included_pop_labels.txt')
    inputs = {'exclude_list': exclude_list, 'poplabels': poplabels}
    outputs = {'pop_label_list': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {output_dir}
    grep -v -f {exclude_list} {poplabels} > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# prepare input files for RELATE
def prepare_files(exclude_list, population=None, haps=None, sample=None, ancestor=None, mask=None, poplabels=None):
    output_dir = f'steps/{population}/relate'
    inputs = {'haps': haps, 'sample': sample, 'ancestor': ancestor, 'mask':mask, 'poplabels':poplabels, 'exclude_list':exclude_list}
    output_path = os.path.join(output_dir, '1000g_ppl_phased_haplotypes')
    outputs = {'haps': output_path + '.haps', 'sample': output_path + '.sample', 'dist': output_path + '.dist', 'poplabels': output_path + '.poplabels', 'annot': output_path + '.annot'} 
    options = {'memory': '20g', 'walltime': '10:00:00'}
    spec = f'''
    mkdir -p {output_dir}
    {config['relate_dist_dir']}/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps {haps} --sample {sample} --ancestor {ancestor} --mask {mask} --remove_ids {exclude_list} --poplabels {poplabels} -o {output_path}
    sleep 20
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# compute sfs to make sure singletons are not missing (sanity check)
# zcat 1000g_LWK_phased_haplotypes.haps.gz | cut -d ' ' -f 4- | tr -d -c '1\n' | awk '{ print length; }' | sort -n | uniq -c

# run the inference of tree sequences using RELATE
def relate(genetic_map, population=None, sample_relate=None, haps_relate=None, annot_relate=None, dist_relate=None):
    output_dir = f'steps/{population}/relate'
    file_base_name = '1000g_ppl_phased_haplotypes'
    output_path = os.path.join(output_dir, file_base_name)
    inputs = {'sample_relate': sample_relate, 'haps_relate': haps_relate, 'annot_relate': annot_relate, 'dist_relate': dist_relate}
    outputs = {'anc': output_path + '.anc', 'mut': output_path + '.mut'}
    options = {'memory': '24g', 'walltime': '10:00:00'}
    # program creates a temporary folder for temporary files and if it already exists relate won't run
    spec= f'''
    mkdir -p {output_dir}
    cd {output_dir}
    rm -rf {file_base_name}
    /bin/Relate --mode All -m 1.25e-8 -N 20000 --sample {sample_relate} --haps {haps_relate} --map {genetic_map} --annot {annot_relate} --dist {dist_relate} --memory 20 -o {file_base_name}
    sleep 90
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# estimate historical population size trajectory from initially inferred tree sequences
# setting --threshold 0. This is so that the branch lengths in all trees are updated for the estimated population size history. 
def estimate_ppl_size(anc_size, population=None, mut_size=None, poplabels_size=None):
    output_dir = f'steps/{population}/relate'
    file_name_input = '1000g_ppl_phased_haplotypes'
    file_name_output = '1000g_ppl_phased_haplotypes_demog'
    output_path = os.path.join(output_dir, file_name_output)
    # inputs: inferred .anc/.mut files and a .poplabels file
    inputs = {'anc_size': anc_size, 'mut_size': mut_size, 'poplabels_size': poplabels_size}
    # outputs: two versions of coalescence rates/population sizes are outputted
    ## .coal --> contains coalescence rates and cross-coalescence rates, treating all samples as one population
    ## *.pairwise.coal/.bin --> coalescence rate file and corresponding binary file containing coalescence rates between pairs of samples
    outputs = {'anc': output_path + '.anc', 'mut': output_path + '.mut', 'coal': output_path + '.coal', 'pairwise_coal': output_path + '.pairwise.coal', 'pairwise_bin': output_path + '.pairwise.bin'}
    options = {'memory': '8g', 'walltime': '08:00:00'}
    spec = f'''
    mkdir -p {output_dir}
    cd {output_dir}
    rm -rf {file_name_output}
    {config['relate_dist_dir']}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i {file_name_input} --poplabels {poplabels_size} -o {file_name_output} --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
    sleep 20
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# detect selection using RELATEs builtin statistic
def detect_selection(anc_selection, population=None, mut_selection=None, poplabels_selection=None):
    output_dir = f'steps/{population}/relate'
    file_name_input = '1000g_ppl_phased_haplotypes_demog'
    file_name_output = '1000g_ppl_phased_haplotypes_selection'
    output_path = os.path.join(output_dir, file_name_output)
    inputs = {'anc_selection': anc_selection, 'mut_selection': mut_selection, 'poplabels_selection': poplabels_selection}
    # .freq --> Records the frequency of the derived allele at generations genN .. gen1
    # .lin --> Records the number of lineages in the tree at generations genN .. gen1 as well as the number of lineages when the mutation had frequency 2
    # .sele --> Records the log10 p-value for selection evidence at generations genN .. gen1 as well as the log10 p-value when the
    # mutation had frequency 2. Log10 p-value is set to 1 if mutation had frequency <= 1 at a generation. 
    outputs = {'freq_selection': output_path + '.freq', 'lin_selection': output_path + '.lin', 'sele_selection': output_path + '.sele'}
    options = {'memory': '20g', 'walltime': '10:00:00'}
    spec = f'''
    mkdir -p {output_dir}
    cd {output_dir}
    rm -rf {file_name_output}
    {config['relate_dist_dir']}/scripts/DetectSelection/DetectSelection.sh -i {file_name_input} -m 1.25e-8 --poplabels steps/{population}/relate -o {file_name_output}
    sleep 80
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# convert to tree sequence file format (tskit)
# this function converts anc/mut files inferred by Relate into the tree sequence file format used by tskit. In the current
# implementation, each tree is stored with new nodes in the tree sequence file format, leading to no compression. In addition,
# information about how long branches persist, and how many mutations map to a branch are lost by this conversion.
def tree_seq(anc_convert=None, population=None, mut_convert=None):
    output_dir = f'steps/{population}/relate'
    file_name_input = '1000g_ppl_phased_haplotypes'
    file_name_output = '1000g_ppl_phased_haplotypes'
    output_path = os.path.join(output_dir, file_name_output)
    # inputs: .anc (genealogical relationships info) and .mut (mutations info)
    inputs = {'anc_convert': anc_convert, 'mut_convert': mut_convert}
    # outputs: .trees (combination of both inputs)
    outputs = {'trees_convert': output_path + '.trees'}
    options = {'memory': '8g', 'walltime': '04:00:00'}
    spec = f'''
    mkdir -p {output_dir}
    cd {output_dir}
    rm -rf {file_name_output}
    {config['relate_dist_dir']}/bin/RelateFileFormats --mode ConvertToTreeSequence -i {file_name_input} -o {file_name_output}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def workflow(working_dir=os.getcwd(), defaults={}, config=None):

    gwf = Workflow(working_dir=working_dir, defaults=defaults)

    targets = defaultdict(list)

    targets['download'] = gwf.target_from_template(f'download',
        download_data(config))

    targets['decode_genetic_maps'] = gwf.target_from_template(
        f'decode_genetic_maps',
        decode_genetic_maps(
            config['decode_hg38_sexavg_per_gen'], 
            'steps/X/genetic_map_chrX.tsv'
            )
        )
    targets['female_haploid'] = gwf.target_from_template(
        f'female_haploid',
        female_haploid(
            targets['download'].outputs['sample_vcf_index']
            )
        )
    targets['haplotype_id'] = gwf.target_from_template(
        f'haplotype_id', 
        haplotype_id(
            targets['female_haploid'].outputs['haplotypes']
            )
        )
    targets['all_pop_labels'] = gwf.target_from_template(
        f'all_pop_labels',
        all_pop_labels(
            targets['haplotype_id'].outputs['ids'], 
            targets['download'].outputs['1000G_2504_seq_index'], 
            targets['download'].outputs['1000G_698_seq_index']
            )
        )
    targets['convert_vcf'] = gwf.target_from_template(
        'convert_vcf', 
        convert_vcf(
            targets['female_haploid'].outputs['haplotypes'], 
            targets['all_pop_labels'].outputs['poplabels']
            )
        )

    populations = ['KHV']
    # populations = ['LWK', 'GWD', 'ESN', 'MSL', 'YRI', 'GBR', 'FIN', 'IBS', 'TSI', 'CDX', 'CHB', 'CHS', 'JPT', 'KHV']

    for population in populations:
        
        # exlcude related
        related_target = gwf.map(exclude_related,
                                 [(targets['download'].outputs['1000G_698_seq_index'], population)],
                                 name=f"exclude_related_{population}")
        targets[f'exclude_related_{population}'] = related_target
        related = related_target.outputs[0]  # list

        # get ids for other populations
        input_other_ppl = [(targets['all_pop_labels'].outputs['poplabels'], population)]
        targets['ids_other_ppl_{population}'] = gwf.map(ids_other_ppl, 
                                   input_other_ppl, 
                                   name=f"ids_other_ppl_{population}")

        # combine related and other populations
        targets[f'combine_files_{population}'] = gwf.map(combine_files, 
                                 targets['ids_other_ppl_{population}'].outputs, 
                                 extra = {'population': population, 
                                          'related':related
                                          }, 
                                 name=f"combine_files_{population}")

        # list of excluded
        targets[f'excluded_list_{population}'] = gwf.map(excluded_list, 
                                      targets[f'combine_files_{population}'].outputs, 
                                      extra = {'population': population, 
                                               'haplotype_id':targets['haplotype_id'].outputs['ids']
                                               }, 
                                      name=f"excluded_list_{population}")

        # list of included
        pop_labels_target = gwf.map(pop_labels, 
                                    targets[f'excluded_list_{population}'].outputs, 
                                    extra = {'population': population, 
                                             'poplabels':targets['all_pop_labels'].outputs['poplabels']
                                             }, 
                                    name=f"pop_labels_{population}")
        targets[f"pop_labels_{population}"] = pop_labels_target

        # prepare input for relate
        targets[f"prepare_files_{population}"] = gwf.map(prepare_files, 
                                 targets[f'excluded_list_{population}'].outputs, 
                                 extra = {'population': population,                                           
                                          'haps': targets['convert_vcf'].outputs['haps'],
                                          'sample': targets['convert_vcf'].outputs['sample'],
                                          'ancestor': targets['download'].outputs['ancestral_fa'], 
                                          'mask':targets['download'].outputs['mask'],  
                                          'poplabels': f'steps/1000g_phased_haplotypes_poplabels.txt'
                                          },
                                name=f"prepare_files_{population}")

        # run relate
        genetic_map = 'steps/genetic_map_chrX.tsv'
        targets[f"relate_{population}/relate"] = gwf.map(relate, 
                                    [genetic_map], 
                                    extra = {'population': population, 
                                             'haps_relate': targets[f"prepare_files_{population}"].outputs[0]['haps'],
                                             'sample_relate': targets[f"prepare_files_{population}"].outputs[0]['sample'], 
                                             'annot_relate': targets[f"prepare_files_{population}"].outputs[0]['annot'], 
                                             'dist_relate': targets[f"prepare_files_{population}"].outputs[0]['dist'],
                                             }, 
                                    name=f"relate_{population}")

        # estimate population sizes
        ppl_size_target = gwf.map(estimate_ppl_size, 
                                  [f'steps/{population}/relate/1000g_ppl_phased_haplotypes.anc'], 
                                  extra = {'population': population, 
                                           'mut_size': targets[f"relate_{population}/relate"].outputs[0]['mut'],
                                           'poplabels_size':targets[f"relate_{population}/relate"].outputs[0]['anc'],
                                           }, 
                                  name=f"estimate_ppl_size_{population}")
        targets[f"estimate_ppl_size_{population}"] = ppl_size_target

        # detect selection
        detect_selection_target = gwf.map(detect_selection, 
                                          [f'steps/{population}/relate/1000g_ppl_phased_haplotypes_demog.anc'], 
                                          extra = {'population': population, 
                                                   'mut_selection': f'steps/{population}/relate/1000g_ppl_phased_haplotypes_demog.mut', 
                                             'poplabels_selection': f'steps/{population}/relate/1000g_ppl_phased_haplotypes_demog.mut'
                                                   }, 
                                          name=f"detect_selection_{population}")
        targets[f"detect_selection_{population}"] = detect_selection_target

        # convert to tree sequence file format (tskit)
        tree_seq_target = gwf.map(tree_seq, 
                                  [f'steps/{population}/relate/1000g_ppl_phased_haplotypes.anc'], 
                                  extra = {'population': population, 
                                           'mut_convert': f'steps/{population}/relate/1000g_ppl_phased_haplotypes.mut'
                                           }, 
                                  name=f"tree_convert_{population}")
        targets[f"tree_convert_{population}"] = tree_seq_target

    return gwf, targets



####################################################################
# Use code like this to run this as standalone workflow: 
####################################################################

import yaml
with open('config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

gwf, targets  = workflow(working_dir=os.getcwd(), 
                        defaults={'account': 'xy-drive'},
                        config=config)

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
