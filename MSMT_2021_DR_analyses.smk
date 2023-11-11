configfile: 'MSMT_2021_DR_analyses.yaml'
variant_folder=config['variant_folder']
#rule all:
#	input:
#		IRNL='IRNL_prevalences_by_region.tsv',
#		IRNGEG='IRNGEG_prevalences_by_region.tsv',
#		k13='k13_prevalences_by_region.tsv',
#		dhfr='dhfr_prevalences_by_region.tsv',
#		dhps='dhps_prevalences_by_region.tsv',
#		k13_kagera='k13_kagera_prevalences_by_district.tsv'

rule all:
	input:
		summary='prevalences/Region:all_3_1_summary.tsv'

rule get_UMIs:
	input:
		cov_counts=variant_folder+'/coverage_AA_table.csv',
		alt_counts=variant_folder+'/alternate_AA_table.csv',
		ref_counts=variant_folder+'/reference_AA_table.csv'
	output:
		UMI_counts='counts/all_AA_counts.pkl'
	script:
		'scripts/get_UMIs.py'


rule get_prevalences:
	input:
		UMI_counts='counts/all_AA_counts.pkl',
		metadata=config['metadata_sheet'],
	output:
		coverage_samples='prevalences/{mutation}_{region}_{cov}_{alt}_cov_samples.txt',
		alternate_samples='prevalences/{mutation}_{region}_{cov}_{alt}_alt_samples.txt'
	script:
		'scripts/get_prevalences.py'


rule make_table:
	input:
		desired_files=expand('prevalences/{file}', file=config['prevalence_files']),
		metadata_sheet=config['metadata_sheet']
	output:
		summary='prevalences/{region}_{cov}_{alt}_summary.tsv'
	script:
		'scripts/make_table.py'

'''
rule intersect_samples:
	input:
		file_list=expand(config['intersections']['wildcard'])
	params:
		file_list=config['intersections']
		
	output:
'''
