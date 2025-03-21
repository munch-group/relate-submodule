# Get raw input data

## get ancestral sequence

Download the human genome as it were at the common ancestor of all humans:

    wget http://ftp.ensembl.org/pub/release-109/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz

Unpack it:

    tar xfvz homo_sapiens_ancestor_GRCh38.tar.gz

## get decode genetic maps 

Get map of recombination rate accross the X chromosome made by DECODE genetics:

    cat ~/ari-intern/data/decode_hg38_sexavg_per_gen.tsv | tail -n +2 | grep chrX | cut -f 2,4,5 | (echo pos COMBINED_rate Genetic_Map ; cat - ; ) > genetic_map_chrX.tsv

## Download VCF files for phased vcf 1000g

Get file with SNPs for individuals from the 1000 genomes project:

    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz

Get index for VCF file:

    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz.tbi

Get documentation files for the phases version of the 1000 genomes data set:

    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/README_SNV_INDEL_phasing_111822.pdf
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/1000G_NYGC_phasing_CHANGE_LOG.pdf

Get sequence index files:

    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index 
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index

Get strict mask of unreliably called bases:

    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.chrX.mask.fasta.gz
gunzip 20160622.chrX.mask.fasta.gz

# Construct input files for RELATE:

## Make females haploid

Turn diploid females into two individual haplotypes (haploid individuals) like males:

    python haploid_vcf.py CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz | gzip > 1000g_phased_haplotypes.vcf.gz

## ID file

Construct files with haplotype IDs:

    bcftools query -l 1000g_phased_haplotypes.vcf.gz > 1000g_phased_haplotypes_ids.txt

## Poplabels file

Construct populations labels mapping each haplotype to a population:

    python make_poplabels.py 1000g_phased_haplotypes_ids.txt 1000G_2504_high_coverage.sequence.index 1000G_698_related_high_coverage.sequence.index > 1000g_phased_haplotypes_poplabels.txt

## Convert X chromosome VCF for all samples to haps/sample format

Convert VCF format to haps/sample format required by RELATE:

    ~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps 1000g_phased_haplotypes.haps --sample 1000g_phased_haplotypes.sample -i 1000g_phased_haplotypes --poplabels 1000g_phased_haplotypes_poplabels.txt

## Excluded haplotypes (related and non-LWK)

Always exclude related individuals:

    grep -v '#' 1000G_698_related_high_coverage.sequence.index | cut -f 10 > 1000g_related_ids.txt

Here I analyze only individuals from the African LWK population. So we find IDs fo haplotypes from all other populations so we can exclude them:

    grep -v 'LWK' 1000g_phased_haplotypes_poplabels.txt | cut -f 1 -d ' ' > 1000g_excluded_pop_ids.txt

Combine the files:

    cat 1000g_related_ids.txt 1000g_excluded_pop_ids.txt | sort | uniq > all_excluded.txt

Construct a list of excluded individuals:

    grep -f all_excluded.txt 1000g_phased_haplotypes_ids.txt > 1000g_excluded_non_LWK_haplotype_ids.txt

Construct a list of excluded individuals In this case leaves only indivuals from the LWK population:

    grep -v -f 1000g_excluded_haplotype_ids.txt 1000g_phased_haplotypes_poplabels.txt > 1000g_phased_haplotypes_LWK_poplabels.txt

## Prepare input files

Prepare input files for RELATE:

    srun --mem-per-cpu=8g --time=04:00:00 --account=xy-drive ~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_phased_haplotypes.haps --sample 1000g_phased_haplotypes.sample --ancestor homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_X.fa --mask 20160622.chrX.mask.fasta --poplabels 1000g_phased_haplotypes_poplabels.txt --remove_ids 1000g_excluded_non_LWK_haplotype_ids.txt -o 1000g_LWK_phased_haplotypes

## Compute sfs to make sure singletons are not missing:

Sanity check (Leo said singletons were missing in one of the data set versions):

    zcat 1000g_LWK_phased_haplotypes.haps.gz | cut -d ' ' -f 4- | tr -d -c '1\n' | awk '{ print length; }' | sort -n | uniq -c

# Run RELATE

Run the inference of tree sequences using RELATE:

    srun --mem-per-cpu=24g --time=08:00:00 --account=xy-drive ~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 --sample 1000g_LWK_phased_haplotypes.sample.gz --haps 1000g_LWK_phased_haplotypes.haps.gz --map genetic_map_chrX.tsv --annot 1000g_LWK_phased_haplotypes.annot --dist 1000g_LWK_phased_haplotypes.dist.gz --memory 20 -o 1000g_LWK_phased_haplotypes

# Estimate population sizes

Estimate historical population size trajectory from initially inferred tree sequences

    srun --mem-per-cpu=8g --cpus-per-task=15 --time=08:00:00 --account=xy-drive ~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_LWK_phased_haplotypes --poplabels 1000g_LWK_phased_haplotypes.poplabels -o 1000g_LWK_phased_haplotypes_demog --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14 --threshhold 0

# Detect selection

Detect selection using RELATEs builtin statistic

    srun --mem-per-cpu=8g --time=04:00:00 --account=xy-drive ~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i 1000g_LWK_phased_haplotypes_demog -m 1.25e-8 --poplabels 1000g_LWK_phased_haplotypes.poplabels -o 1000g_LWK_phased_haplotypes_selection




################





# standard 1000g data set on hg38
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz












# poplabels file
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index 
cat 1000G_2504_high_coverage.sequence.index | grep -v '#' | awk -F'\t' '{ print $10 " " $11 " "$11 " NA"  }' | (echo sample population group sex ; cat - ; ) > 1000g_unrelated_poplabels.txt



# id file
tail -n +2 1000g_unrelated_poplabels.txt | cut -d ' ' -f 1 > 1000g_unrelated_ids.txt

# subset VCF to get only unrelated indivs
bcftools view --samples-file 1000g_unrelated_ids.txt -Oz CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz > 1000g_chrX_phased_unrelated.vcf.gz

<!-- 
# lists of ids
bcftools query -l CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz > 100g_all_ids.txt
grep -v '#' 1000G_698_related_high_coverage.sequence.index | cut -f 10 > 1000g_related_ids.txt
grep -v '#' 1000G_2504_high_coverage.sequence.index | cut -f 10 > 1000g_unrelated_ids.txt 
-->



# convert chromosome vcf for all samples to haps/sample format
~/populationgenomics/software/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps 1000g_all.haps --sample 1000g_all.sample -i 1000g_chrX_phased_unrelated --poplabels 1000g_unrelated_poplabels.txt

# make a list of the all the samples
tail -n +3 1000g_all.sample | cut -f 1 > 1000g_all_ids.txt

# make poplabels file for all the samples extracted from the meta info that goes with the vcf
cat ~/simons/faststorage/data/1000Genomes/metainfo/sample_info.csv | tail -n +2 | awk -F';' '{ print $1 " " $3 " "$3 " NA"  }' | sort | (echo sample population group sex ; cat - ; ) > 1000g_all_poplabels.txt

# make a poplabels file with only the samples in the samples file (in case the vcf meta info file has samples that are not in the vcf)
grep -f 1000g_all_ids.txt 1000g_all_poplabels.txt | (echo sample population group sex ; cat - ; ) > 1000g_poplabels.txt

# extract the subset of samples to run relate on and make a poplabels file
cat ~/simons/faststorage/data/1000Genomes/metainfo/sample_info.csv | tail -n +2 | grep female | awk -F';' '{ print $1 " " $3 " "$3 " NA"  }' | grep -E 'YRI' | sort | grep -f 1000g_poplabels.txt | (echo sample population group sex ; cat - ; ) > 1000g_eur_females_poplabels.txt

# make a list of ids of those samples
tail -n +2 1000g_eur_females_poplabels.txt | cut -f 1 -d ' ' > 1000g_eur_females_ids.txt

# make a list of the samples you want to exclude (all but the above ones)
cat 1000g_all_ids.txt | grep -v -w -f 1000g_eur_females_ids.txt | awk -F';' '{ print $1 }' > 1000g_excluded_ids.txt





# prepare input files
 ~/populationgenomics/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps 1000g_all.haps --sample 1000g_all.sample --ancestor homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_X.fa --mask 20160622.chrX.mask.fasta --poplabels 1000g_poplabels.txt --remove_ids 1000g_excluded_ids.txt -o 1000g_eur_females

# run relate
~/populationgenomics/software/relate/bin/Relate --mode All -m 1.25e-8 -N 20000 -haps 1000g_eur_females.sample.gz --sample 1000g_eur_females.haps.gz --map genetic_map_chrX.tsv --annot 1000g_eur_females.annot --dist 1000g_eur_females.dist.gz --memory 20 -o 1000g_eur_females

# estimate population sizes
#~/populationgenomics/software/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -m 1.25e-8 -N 20000 -i 1000g_eur_females --poplabels 1000g_poplabels.txt -o 1000g_eur_females_demog --threshold 0 --num_iter 5 --years_per_gen 29 --threads 14
sbatch sbatch_demog.sh

# detect selection
~/populationgenomics/software/relate/scripts/DetectSelection/DetectSelection.sh -i eur_phased_demog -m 1.25e-8 --poplabels eur_poplabels.txt -o eur_selection

