import pandas as pd
import sys
import gzip

_, vcf_file_name = sys.argv

males = []
females = []

print('sexing samples', file=sys.stderr)
with gzip.open(vcf_file_name, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
        else:
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            if any('|' not in call for call in calls):
                first_non_par_pos = int(POS)
                print('first non-par pos', first_non_par_pos, file=sys.stderr)
                with open('sexes.txt', 'w') as sex_file:
                    for i, call in enumerate(calls):
                        if '|' in call:
                            females.append(all_samples[i])
                            print(all_samples[i], 'F', file=sex_file)
                        else:
                            males.append(all_samples[i])
                            print(all_samples[i], 'M', file=sex_file)
                assert len(all_samples) == len(males) + len(females)
                break

print('writing haploids', file=sys.stderr)
with gzip.open(vcf_file_name, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *all_samples = line.split()
                ids = []
                for sample in all_samples:
                    if sample in males:
                        ids.append(sample)
                    else:
                        ids.append(sample + "_1")
                        ids.append(sample + "_2")

                start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
                start_line.extend(ids)
                print('\t'.join(start_line))
            else:
                print(line, end='')
        else:
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *calls = line.split()
            if int(POS) < first_non_par_pos:
                continue
            haploid_calls = []
            for call in calls:
                for x in call.split('|'):
                    haploid_calls.append(x)
            if len(haploid_calls) > 2* len(females) + len(males):
                print('aborting at par2:', POS,file=sys.stderr)
                break
            start_line = [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT]
            start_line.extend(haploid_calls)
            print('\t'.join(start_line))
