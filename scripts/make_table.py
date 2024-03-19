prevalence_files=snakemake.input.desired_files
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

#prevalences/Region:all_3_1_summary.tsv

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

def format_line(prevalence_dict, category_value, filtered_files, output_table):
		from statsmodels.stats.proportion import proportion_confint
		output_lines=[[category_value], [category_value], [category_value]]
		for prevalence_file in filtered_files:
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
	output_table.write('\t'.join(output_header)+'\n')
	for category_value in category_values:
		format_line(prevalence_dict, category_value, filtered_files, output_table)
	format_line(prevalence_dict, 'overall', filtered_files, output_table)
#k13-Arg561His_Region:all_3_1_cov_samples.txt

prevalence_dict={}
category_type, filter_type=summary.split('/')[-1].split('_')[0].split(':')
categories=set([])
output_header=[category_type]
filtered_files=[]
for prevalence_file in prevalence_files:
	print('prevalence file is', prevalence_file)
	prefix='_'.join(prevalence_file.split('_')[:-2])
	file_name=prevalence_file.split('/')[-1]
	mutation, region, cov, alt=file_name.split('_')[:4]
	print('region is', region)
	obs_category, obs_filter=region.split(':')
	if obs_filter==filter_type and obs_category==category_type:
		filtered_files.append(prevalence_file)
		output_header.append(mutation.split('/')[-1])
		cov_counts=get_counts(prevalence_file, category_type, filter_type)
		alt_counts=get_counts(prefix+'_alt_samples.txt', category_type, filter_type)
		prevalence_dict[prevalence_file]=[cov_counts, alt_counts]
		categories=categories | set(cov_counts.keys()) | set(alt_counts.keys())
categories=categories-set(['overall'])
format_table(prevalence_dict, categories, output_header, filtered_files)
