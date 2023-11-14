'''
unlike get_prevalences, which retrieves the prevalences of specific mutations,
this program retrieves all mutations (either within a specific gene or all
genes) that meet a threshold coverage and a threshold alternate count. This
version imposes an additional filter by requiring at least x samples to contain
the mutation of interest in order to report it.
'''
import pickle
filtered_mutations=snakemake.output.filtered_mutations
search_term, cov_threshold, alt_threshold, count_threshold=filtered_mutations.split('/')[-1].split('_')[:4]
cov_threshold, alt_threshold, count_threshold=map(int, [cov_threshold, alt_threshold, count_threshold])
filtered_mutations=open(filtered_mutations, 'w')
metadata=snakemake.input.metadata
valid_samples=set([line.strip().split(',')[11] for line in open(metadata)][1:])
#print('valid samples are', valid_samples)

UMI_counts=snakemake.input.UMI_counts
UMI_counts=pickle.load(open(UMI_counts, 'rb'))

def check_thresholds(valid_samples, search_term, cov_threshold, alt_threshold, count_threshold):
	cov_dict={}
	if search_term=='all':
		search_term=''
	for sample in UMI_counts:
		if sample in valid_samples:
			for mutation in UMI_counts[sample]:
#				print(search_term, mutation)
				if search_term in mutation:
					cov_dict.setdefault(mutation, [set([]), set([])])
					cov, ref, alt=UMI_counts[sample][mutation]
					if cov>=cov_threshold:
						cov_dict[mutation][0].add(sample)
					if alt>=alt_threshold and cov>=cov_threshold:
						cov_dict[mutation][1].add(sample)
	for mutation in list(cov_dict.keys()):
		if len(cov_dict[mutation][1])<count_threshold:
			cov_dict.pop(mutation)
	return cov_dict



cov_dict=check_thresholds(valid_samples, search_term, cov_threshold, alt_threshold, count_threshold)
filtered_mutations.write('\n'.join(list(cov_dict.keys())))
cov_threshold, alt_threshold, count_threshold=map(str, [cov_threshold, alt_threshold, count_threshold])
for mutation in cov_dict:
	cov_file=open('prevalences_by_threshold/'+'_'.join([mutation, cov_threshold, alt_threshold, count_threshold, 'cov.txt']), 'w')
	alt_file=open('prevalences_by_threshold/'+'_'.join([mutation, cov_threshold, alt_threshold, count_threshold, 'alt.txt']), 'w')
	cov_file.write('\n'.join(list(cov_dict[mutation][0])))
	alt_file.write('\n'.join(list(cov_dict[mutation][1])))
