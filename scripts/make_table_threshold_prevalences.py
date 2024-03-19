output_folder=snakemake.params.output_folder
threshold_mutations=snakemake.input.threshold_mutations
prevalence_files=[]
for mutation_file in threshold_mutations:
	print(mutation_file)
	gene, cov, alt, count=mutation_file.split('/')[-1].split('_')[:4]
	mutations=[line.strip() for line in open(mutation_file)]
	mutations=sorted([[int(mutation.split('-')[-1][3:-3]), mutation] for mutation in mutations])
	for mutation in mutations:
		mutation=mutation[1]
		mutation=mutation.replace(' ', '_').strip()
		prevalence_files.append(output_folder+'/prevalences_by_threshold/'+'_'.join([mutation, cov, alt, count, 'cov.txt']))

metadata_sheet=snakemake.input.metadata_sheet
summary=snakemake.output.summary
heirarchy=snakemake.params.heirarchy
heirarchy=heirarchy.split(':')

def make_metadata_dict(metadata_sheet):
	metadata_dict={}
	h={}
	for line_number, line in enumerate(open(metadata_sheet)):
		line=line.strip().split(',')
		if line_number==0:
			for column_number, column in enumerate(line):
				h[column]=column_number
		else:
			sample=line[h['Sample_ID']]
			metadata_dict[sample]=line
	return metadata_dict, h

metadata_dict, header_dict=make_metadata_dict(metadata_sheet)

def get_counts(sample_file, category, filter_type):
	parsed_counts={}
	samples=[line.strip() for line in open(sample_file)]
	category_index=heirarchy.index(category)
	subsidiary=heirarchy[category_index+1]
	for sample in samples:
		metadata=metadata_dict[sample]
		category_value=metadata[header_dict[category]]
		subsidiary_value=metadata[header_dict[subsidiary]]
		if filter_type=='all':
			parsed_counts[category_value]=parsed_counts.setdefault(category_value, 0)+1
		else:
			parsed_counts[subsidiary_value]=parsed_counts.setdefault(subsidiary_value, 0)+1
		parsed_counts['overall']=parsed_counts.setdefault('overall', 0)+1
	return parsed_counts

def format_line(mutation, prevalence_file, category_values, filtered_files, output_table):
		from statsmodels.stats.proportion import proportion_confint
		output_lines=[[mutation], [''], ['']]
		for category_value in category_values:
			cov_dict, alt_dict=prevalence_dict[prevalence_file]
			cov_count=cov_dict.setdefault(category_value, 0)
			alt_count=alt_dict.setdefault(category_value, 0)
			if cov_count>0:
				lower_bound, upper_bound=proportion_confint(count=alt_count, nobs=cov_count, alpha=0.1)
				lower_bound, upper_bound=round(lower_bound*100, 1), round(upper_bound*100, 1)
				output_lines[0].append(f'{alt_count}/{cov_count}')
				output_lines[1].append(f'{round(alt_count/cov_count*100, 1)}%')
				output_lines[2].append(f'{lower_bound}%-{upper_bound}%')
			else:
				output_lines[0].append(f'no coverage')
				output_lines[1].append(f'no coverage')
				output_lines[2].append(f'no coverage')
		for output_line in output_lines:
			output_table.write('\t'.join(output_line)+'\n')

def format_table(prevalence_dict, category_values, output_header, filtered_files):
	output_table=open(summary, 'w')
	output_header=output_header+category_values
	output_table.write('\t'.join(output_header)+'\n')
	for mutation_name, prevalence_file in sorted(filtered_files):
		format_line(mutation_name, prevalence_file, category_values, filtered_files, output_table)

prevalence_dict={}
category_type, filter_type='Region', 'all'
categories=set([])
output_header=['mutation']
filtered_files=[]
sorting_table=[]
for prevalence_file in prevalence_files:
	prefix='_'.join(prevalence_file.split('_')[:-1])
	mutation, cov, alt, count=prevalence_file.split('/')[-1].split('_')[1:5]
	obs_category, obs_filter='Region', 'all'
	if obs_filter==filter_type and obs_category==category_type:
		sorting_table.append([mutation, prevalence_file])
		cov_counts=get_counts(prevalence_file, category_type, filter_type)
		alt_counts=get_counts(prefix+'_alt.txt', category_type, filter_type)
		prevalence_dict[prevalence_file]=[cov_counts, alt_counts]
		categories=categories | set(cov_counts.keys()) | set(alt_counts.keys())
		filtered_files=sorted(sorting_table) #sorting filtered_files by mutation name

categories=sorted(list(categories-set(['overall'])))+['overall']
format_table(prevalence_dict, categories, output_header, filtered_files)
