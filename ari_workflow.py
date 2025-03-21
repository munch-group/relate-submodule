
# %% [markdown]
# ---
# title: GWF workflow
# execute:
#   eval: false
# ---

# %%
 
from gwf import Workflow
import sys, os, re
from collections import defaultdict
from pathlib import Path
import pandas as pd
from gwf import Workflow, AnonymousTarget
from gwf.workflow import collect

gwf = Workflow(defaults={'account': 'ari-intern'})


# directories
out_dir = '/home/ari/ari-intern/people/ari/ariadna-intern/steps'
data_big = '/faststorage/project/ari-intern/data'
script_dir = '/faststorage/project/ari-intern/people/ari/ariadna-intern/scripts'
data_dir = '/faststorage/project/ari-intern/people/ari/ariadna-intern/steps/1000genome'

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


# map of recombination rate across the X chromosome made by DECODE genetics
def decode_genetic_maps(decode_hg38_sexavg_per_gen, genetic_map_chrX):
    inputs = [decode_hg38_sexavg_per_gen]
    outputs = [genetic_map_chrX]
    options = {'memory': '1g', 'walltime': '00:10:00'}
    spec = f'''
    cat {decode_hg38_sexavg_per_gen} | tail -n +2 | grep chrX | cut -f 2,4,5 | (echo pos COMBINED_rate Genetic_Map ; cat - ; ) > {genetic_map_chrX}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


decode_hg38_sexavg_per_gen=f'{data_big}/decode_hg38_sexavg_per_gen.tsv'
genetic_map_chrX=f'{out_dir}/genetic_map_chrX.tsv'

gwf.target_from_template(f'decode_genetic_maps',
    decode_genetic_maps(decode_hg38_sexavg_per_gen, genetic_map_chrX))


# turn diploid females (XX) into two individual haplotypes (haploid individuals) like males
def female_haploid(haploid_vcf, chrX_filtered_eagle2_phased, phased_haplotypes):
    inputs = [haploid_vcf, chrX_filtered_eagle2_phased]
    outputs = [phased_haplotypes]
    options = {'memory': '10g', 'walltime': '01:20:00'}
    spec = f'''
    python {haploid_vcf} {chrX_filtered_eagle2_phased} | gzip > {phased_haplotypes}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

haploid_vcf=f'{script_dir}/haploid_vcf.py'
chrX_filtered_eagle2_phased='/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000genome/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz'
phased_haplotypes=f'{out_dir}/1000g_phased_haplotypes.vcf.gz'

gwf.target_from_template(f'female_haploid',
    female_haploid(haploid_vcf, chrX_filtered_eagle2_phased, phased_haplotypes))


# construct files with haplotype IDs
def haplotype_id(phased_haplotypes, phased_haplotypes_id):
    inputs = [phased_haplotypes]
    outputs = [phased_haplotypes_id]
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    # conda install -c bioconda bcftools
    # conda install openssl   ## to install libcrypto.so.1.0.0 library
    bcftools query -l {phased_haplotypes} > {phased_haplotypes_id}
    sleep 5
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

phased_haplotypes=f'{out_dir}/1000g_phased_haplotypes.vcf.gz'
phased_haplotypes_id=f'{out_dir}/1000g_phased_haplotypes_ids.txt'

gwf.target_from_template(f'haplotype_id', haplotype_id(phased_haplotypes, phased_haplotypes_id))


# construct populations labels mapping each haplotype to a population
# (group haplotypes according to the population to which the individuals carrying those haplotypes belong)
def pop_labels(make_poplabels, phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index, phased_haplotypes_poplabels):
    inputs = [make_poplabels, phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index]
    outputs = [phased_haplotypes_poplabels]
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    python {make_poplabels} {phased_haplotypes_id} {high_coverage_seq_index} {related_high_coverage_seq_index} > {phased_haplotypes_poplabels} 
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

make_poplabels=f'{script_dir}/make_poplabels.py'
phased_haplotypes_id=f'{out_dir}/1000g_phased_haplotypes_ids.txt'
high_coverage_seq_index=f'{data_dir}/seq_index/1000G_2504_high_coverage.sequence.index'
related_high_coverage_seq_index=f'{data_dir}/seq_index/1000G_698_related_high_coverage.sequence.index'
phased_haplotypes_poplabels=f'{out_dir}/1000g_phased_haplotypes_poplabels.txt'

gwf.target_from_template(f'pop_labels',
    pop_labels(make_poplabels, phased_haplotypes_id, high_coverage_seq_index, related_high_coverage_seq_index, phased_haplotypes_poplabels))


# Define the function to convert VCF to haps/sample format
def convert_vcf(RelateFileFormats, phased_haplotypes_haps, phased_haplotypes_sample, phased_haplotypes, phased_haplotypes_poplabels):
    inputs = [RelateFileFormats, phased_haplotypes_poplabels, phased_haplotypes]
    outputs = [phased_haplotypes_haps, phased_haplotypes_sample]
    options = {'memory': '10g', 'walltime': '01:00:00'}
    spec = f'''
    {RelateFileFormats} --mode ConvertFromVcf --haps {phased_haplotypes_haps} --sample {phased_haplotypes_sample} -i {phased_haplotypes.replace('vcf.gz', '')} --poplabels {phased_haplotypes_poplabels}
    sleep 20
    touch {phased_haplotypes_haps}
    touch {phased_haplotypes_sample}
    '''
    # Returning outputs as well
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Define the file paths for input and output
RelateFileFormats = '/faststorage/project/ari-intern/people/ari/relate/bin/RelateFileFormats'
phased_haplotypes_poplabels = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000g_phased_haplotypes_poplabels.txt'
phased_haplotypes = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000g_phased_haplotypes.vcf.gz'
phased_haplotypes_haps = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000g_phased_haplotypes.haps'
phased_haplotypes_sample = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000g_phased_haplotypes.sample'

# Creating the target using the function
convert_vcf_target = convert_vcf(RelateFileFormats, phased_haplotypes_haps, phased_haplotypes_sample, phased_haplotypes, phased_haplotypes_poplabels)

# Adding the target to the workflow
gwf.target_from_template(f'convert_vcf', convert_vcf_target)



## start with specific population ##

# exclude related individuals to avoid biases arising from shared genetic material
def exclude_related(path, population):
    output_dir = f'{out_dir}/{population}/excluded'
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
    output_dir = f'{out_dir}/{population}/excluded'
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
def combine_files(path, related=None):
    output_path = modify_path(path, base='', suffix='excluded_combined.txt')
    out_dir = modify_path(output_path, base='', suffix='')
    inputs = {'path': path, 'related': related}
    outputs = {'path': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {out_dir}
    cat {path} {related} | sort | uniq > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct a list of excluded individuals
def excluded_list(path, haplotype_id=None):
    output_path = modify_path(path, base='', suffix='excluded_list.txt')
    out_dir = modify_path(output_path, base='', suffix='')
    inputs = {'path': path, 'haplotype_id': haplotype_id}
    outputs = {'exclude_list': output_path}
    options = {'memory': '10g', 'walltime': '00:60:00'}
    spec = f'''
    mkdir -p {out_dir}
    grep -f {path} {haplotype_id} > {output_path}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# construct a list of only individuals from the population of interest
def pop_labels(exclude_list, poplabels=None):
    output_dir = f'{out_dir}/{population}/included'
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
def prepare_files(exclude_list, haps=None, sample=None, ancestor=None, mask=None, poplabels=None):
    directory = '/home/ari/ari-intern/people/ari/ariadna-intern/steps'
    output_dir = f'{directory}/{population}/relate'
    inputs = {'haps': haps, 'sample': sample, 'ancestor': ancestor, 'mask':mask, 'poplabels':poplabels, 'exclude_list':exclude_list}
    output_path = os.path.join(output_dir, '1000g_ppl_phased_haplotypes')
    # outputs: .haps, .sample, .dist (if --mask specified), .poplabels (if remove_ids & poplabels specified), .annot (if poplabels specified)
    outputs = {'haps': output_path + '.haps', 'sample': output_path + '.sample', 'dist': output_path + '.dist', 'poplabels': output_path + '.poplabels', 'annot': output_path + '.annot'} 
    options = {'memory': '20g', 'walltime': '10:00:00'}
    spec = f'''
    mkdir -p {output_dir}
    /home/ari/ari-intern/people/ari/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps {haps} --sample {sample} --ancestor {ancestor} --mask {mask} --remove_ids {exclude_list} --poplabels {poplabels} -o {output_path}
    sleep 20
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# compute sfs to make sure singletons are not missing (sanity check)
# zcat 1000g_LWK_phased_haplotypes.haps.gz | cut -d ' ' -f 4- | tr -d -c '1\n' | awk '{ print length; }' | sort -n | uniq -c

# run the inference of tree sequences using RELATE
def relate(genetic_map, sample_relate=None, haps_relate=None, annot_relate=None, dist_relate=None):
    output_dir = f'/home/ari/ari-intern/people/ari/ariadna-intern/steps/{population}/relate/run_relate'
    file_name = '1000g_ppl_phased_haplotypes'
    output_path = os.path.join(output_dir, file_name)
    inputs = {'sample_relate': sample_relate, 'haps_relate': haps_relate, 'annot_relate': annot_relate, 'dist_relate': dist_relate}
    outputs = {'anc': output_path + '.anc', 'mut': output_path + '.mut'}
    options = {'memory': '24g', 'walltime': '10:00:00'}
    # program creates a temporary folder for temporary files and if it already exists relate won't run
    spec= f'''
    mkdir -p {output_dir}
    cd {output_dir}
    rm -rf {file_name}
    /home/ari/ari-intern/people/ari/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample {sample_relate} --haps {haps_relate} --map {genetic_map} --annot {annot_relate} --dist {dist_relate} --memory 20 -o {file_name}
    sleep 90
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# estimate historical population size trajectory from initially inferred tree sequences
# setting --threshold 0. This is so that the branch lengths in all trees are updated for the estimated population size history. 
def estimate_ppl_size(anc_size=None, mut_size=None, poplabels_size=None):
    output_dir = f'/home/ari/ari-intern/people/ari/ariadna-intern/steps/{population}/relate/run_relate'
    file_name_input = '1000g_ppl_phased_haplotypes'
    file_name_output = '1000g_ppl_phased_haplotypes_demog'
    output_path = os.path.join(output_dir, file_name_output)
    # inputs: inferred .anc/.mut files and a .poplabels file
    inputs = {'anc_size': anc_size, 'mut_size': mut_size, 'poplabels_size': poplabels_size}
    # outputs: two versions of coalescence rates/population sizes are outputted
    ## .coal --> contains coalescence rates and cross-coalescence rates, treating all samples as one population
    ## *.pairwise.coal/.bin --> coalescence rate file and corresponding binary file containing coalescence rates between pairs of samples
    outputs = {'coal': output_path + '.coal', 'pairwise_coal': output_path + '.pairwise.coal', 'pairwise_bin': output_path + '.pairwise.bin'}
    options = {'memory': '8g', 'walltime': '08:00:00'}
    spec = f'''
    mkdir -p {output_dir}
    cd {output_dir}
    rm -rf {file_name_output}
    /home/ari/ari-intern/people/ari/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i {file_name_input} --poplabels {poplabels_size} -o {file_name_output} --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0
    sleep 20
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



# detect selection using RELATEs builtin statistic
def detect_selection(anc_selection=None, mut_selection=None, poplabels_selection=None):
    output_dir = f'/home/ari/ari-intern/people/ari/ariadna-intern/steps/{population}/relate/run_relate'
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
    /home/ari/ari-intern/people/ari/relate/scripts/DetectSelection/DetectSelection.sh -i {file_name_input} -m 1.25e-8 --poplabels {poplabels_size} -o {file_name_output}
    sleep 80
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



# convert to tree sequence file format (tskit)
# this function converts anc/mut files inferred by Relate into the tree sequence file format used by tskit. In the current
# implementation, each tree is stored with new nodes in the tree sequence file format, leading to no compression. In addition,
# information about how long branches persist, and how many mutations map to a branch are lost by this conversion.
def tree_seq(anc_convert=None, mut_convert=None):
    output_dir = f'/home/ari/ari-intern/people/ari/ariadna-intern/steps/{population}/relate/run_relate'
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
    /home/ari/ari-intern/people/ari/relate/bin/RelateFileFormats --mode ConvertToTreeSequence -i {file_name_input} -o {file_name_output}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
                    



#populations = ['GWD', 'ESN', 'MSL', 'YRI', 'LWK'] # african ancestry --> DONE
populations = ['KHV']
#populations = ['GBR', 'FIN', 'IBS',# 'TSI'] + PUR (puerto ricans)# european ancestry --> DONE
#populations = ['CDX', 'CHB', 'CHS', #'JPT', 'KHV'] # east asian ancestry
#populations = ['LWK', 'GWD', 'ESN', 'MSL', 'YRI', 'GBR', 'FIN', 'IBS', 'TSI', 'CDX', 'CHB', 'CHS', 'JPT', 'KHV']


# append a unique identifier to each target name to ensure they are distinct (name = ...{population})
for population in populations:
    # exlcude related
    input_related = [(f'{data_dir}/seq_index/1000G_698_related_high_coverage.sequence.index', population)]
    related_target = gwf.map(exclude_related, input_related, name=f"exclude_related_{population}")
    related = related_target.outputs[0]  # list


    # get ids for other populations
    input_other_ppl = [(f'{out_dir}/1000g_phased_haplotypes_poplabels.txt', population)]
    other_ppl_target = gwf.map(ids_other_ppl, input_other_ppl, name=f"ids_other_ppl_{population}")

    # combine related and other populations
    combine_target = gwf.map(combine_files, other_ppl_target.outputs, extra = {'related':related}, name=f"combine_files_{population}")

    # list of excluded
    haplotype_ids = f'{out_dir}/1000g_phased_haplotypes_ids.txt'
    exclude_list_target = gwf.map(excluded_list, combine_target.outputs, extra = {'haplotype_id':haplotype_ids}, name=f"excluded_list_{population}")

    # list of included
    poplabels = f'{out_dir}/1000g_phased_haplotypes_poplabels.txt'
    pop_labels_target = gwf.map(pop_labels, exclude_list_target.outputs, extra = {'poplabels':poplabels}, name=f"pop_labels_{population}")


    # RELATE DIRECTORY !
    relate_dir = f'/home/ari/ari-intern/people/ari/ariadna-intern/steps/{population}/relate' # relate directory


    # PREPARE INPUT
    haps = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000g_phased_haplotypes.haps'
    sample = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/1000g_phased_haplotypes.sample'
    ancestor = f'{data_dir}/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_X.fa'
    mask = f'{data_dir}/20160622.chrX.mask.fasta'
    poplabels = f'{out_dir}/1000g_phased_haplotypes_poplabels.txt'

    prepare_target = gwf.map(prepare_files, exclude_list_target.outputs, 
                            extra = {'haps': haps, 'sample': sample, 'ancestor': ancestor, 'mask':mask, 'poplabels':poplabels}, name=f"prepare_files_{population}")


    # RUN RELATE
    sample_relate = f'{relate_dir}/1000g_ppl_phased_haplotypes.sample.gz'
    haps_relate = f'{relate_dir}/1000g_ppl_phased_haplotypes.haps.gz'
    genetic_map = '/home/ari/ari-intern/people/ari/ariadna-intern/steps/genetic_map_chrX.tsv'
    annot_relate = f'{relate_dir}/1000g_ppl_phased_haplotypes.annot'
    dist_relate = f'{relate_dir}/1000g_ppl_phased_haplotypes.dist.gz'

    run_relate_target = gwf.map(relate, [genetic_map], extra = {'haps_relate': haps_relate, 'sample_relate': sample_relate, 'annot_relate': annot_relate, 'dist_relate': dist_relate}, name=f"relate_{population}")


    # ESTIMATE POPULATION SIZES
    anc_size = f'{relate_dir}/run_relate/1000g_ppl_phased_haplotypes.anc'
    mut_size  = f'{relate_dir}/run_relate/1000g_ppl_phased_haplotypes.mut'
    poplabels_size = f'{relate_dir}/1000g_ppl_phased_haplotypes.poplabels'
    ppl_size_target = gwf.map(estimate_ppl_size, [anc_size], extra = {'mut_size': mut_size, 'poplabels_size': poplabels_size}, name=f"estimate_ppl_size_{population}")


    # DETECT SELECTION
    anc_selection = f'{relate_dir}/run_relate/1000g_ppl_phased_haplotypes_demog.anc.gz'
    mut_selection  = f'{relate_dir}/run_relate/1000g_ppl_phased_haplotypes_demog.mut.gz'
    poplabels_selection = f'{relate_dir}/1000g_ppl_phased_haplotypes.poplabels'
    detect_selection_target = gwf.map(detect_selection, [anc_selection], extra = {'mut_selection': mut_selection, 'poplabels_selection': poplabels_selection}, name=f"detect_selection_{population}")


   # CONVERT TO TREE SEQUENCE FILE FORMAT (tskit)
    anc_convert = f'{relate_dir}/run_relate/1000g_ppl_phased_haplotypes.anc'
    mut_convert  = f'{relate_dir}/run_relate/1000g_ppl_phased_haplotypes.mut'
    tree_seq_target = gwf.map(tree_seq, [anc_convert], extra = {'mut_convert': mut_convert}, name=f"tree_convert_{population}")