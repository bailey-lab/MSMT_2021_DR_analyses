'''
merges my metadata (metadata originally sent by Dativa plus two additional
metadata files she emailed me) with metadata that Filbert fixed up. I later used
email correspondence with Filbert to manually further clean this output file.
'''

filbert_metadata='/nfs/jbailey5/baileyweb/asimkin/github_pipelines/MSMT_2021_DR_analyses/scripts/NEW_MIP_2021_Revised_Metadata_Alfred_combined_preprint_metada.csv'
my_metadata='/nfs/jbailey5/baileyweb/asimkin/github_pipelines/MSMT_2021_DR_analyses/scripts/preprint_metadata-all_Dativa_metadata.csv'
def parse_metadata(input_file):
	metadata_dict={}
	for line_number, line in enumerate(open(input_file)):
		if line_number>0:
			split_line=line.split(',')
			sample=split_line[11]
			#print('sample is', sample)
			metadata_dict.setdefault(sample, []).append(line)	
		else:
			header=line
	return metadata_dict, header

def merge_dicts(my_dict, filbert_dict):
	merged_dict={}
	my_samples, filbert_samples=list(my_dict.keys()), list(filbert_dict.keys())
	extra=set(my_dict.keys())-set(filbert_dict.keys())
	print(extra)
	for key in filbert_dict:
		merged_dict[key]=filbert_dict[key]
	for key in extra:
		merged_dict[key]=my_dict[key]
	return merged_dict

def print_merged_metadata(merged_dict, header, output_path):
	output_file=open(output_path, 'w')
	output_file.write(header)
	for key in sorted(list(merged_dict.keys())):
		for line in merged_dict[key]:
			output_file.write(line)


my_dict, header=parse_metadata(my_metadata)
print('header is', header)
filbert_dict, header=parse_metadata(filbert_metadata)
print('header is', header)
merged_dict=merge_dicts(my_dict, filbert_dict)
print_merged_metadata(merged_dict, header, 'manuscript_metadata.csv')
