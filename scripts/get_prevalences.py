import pickle
UMI_counts=snakemake.input.UMI_counts
UMI_counts=pickle.load(open(UMI_counts, 'rb'))
metadata=snakemake.input.metadata
mutation=snakemake.wildcards.mutation
cov_threshold=int(snakemake.wildcards.cov)
alt_threshold=int(snakemake.wildcards.alt)
region=snakemake.wildcards.region.split(':')
region_type, region_value=region
coverage_samples=open(snakemake.output.coverage_samples, 'w')
alternate_samples=open(snakemake.output.alternate_samples, 'w')

def check_thresholds(sample):
	num, denom=False, False
	if sample in UMI_counts:
		for UMI_mutation in UMI_counts[sample]:
			if mutation in UMI_mutation:
				cov, ref, alt=UMI_counts[sample][UMI_mutation]
				if cov>=cov_threshold:
					denom=True
				if alt>=alt_threshold and cov>=cov_threshold:
					num=True
	return num, denom

h={}
for line_number, line in enumerate(open(metadata)):
	line=line.strip().split(',')
	if line_number==0:
		for column_number, column in enumerate(line):
			h[column]=column_number
	else:
		sample=line[h['Sample_ID']]
		if region_value=='all' or line[h[region_type]]==region_value:
			num, denom=check_thresholds(sample)
			if denom:
				coverage_samples.write(sample+'\n')
			if num:
				alternate_samples.write(sample+'\n')
