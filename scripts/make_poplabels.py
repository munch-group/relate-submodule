
import sys

_, id_file, unrelated_meta_info_file, related_meta_info_file = sys.argv

id_pop_pairs = []
with open(unrelated_meta_info_file) as f:
    for line in f:      
        if line.startswith("#"):
            continue
        sample_id, pop_id = line.split('\t')[9:11]
        id_pop_pairs.append((sample_id.strip(), pop_id.strip()))
with open(related_meta_info_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        sample_id, pop_id = line.split('\t')[9:11]
        id_pop_pairs.append((sample_id.strip(), pop_id.strip()))      

pops = dict(id_pop_pairs)

print('sample population group sex')
with open(id_file) as f:
    for line in f:
        hap_id = line.strip()
        if hap_id.endswith('_1') or hap_id.endswith('_2'):
            pop = pops[hap_id[:-2]]
        else:
            pop = pops[hap_id]
        print(hap_id, pop, pop, 1)