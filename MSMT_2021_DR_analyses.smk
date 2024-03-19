configfile: 'MSMT_2021_DR_analyses.yaml'
output_folder=config['output_folder']
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
		all_summary=output_folder+'/summaries/Region:all_3_1_summary.tsv',
		Kagera_summary=output_folder+'/summaries/Region:Kagera_3_1_summary.tsv',
		threshold_summary=output_folder+'/summaries/10_3_3_threshold-summary.tsv',
		background_mutations=output_folder+'/background_mutations/561_on_Asian_backgrounds.tsv',
		D_samples=output_folder+'/prevalences/mdr1-Asp1246_3_1_cov.txt'


rule get_UMIs:
	'''
	generates a pickled table of UMI counts
	'''
	input:
		cov_counts=variant_folder+'/coverage_AA_table.csv',
		alt_counts=variant_folder+'/alternate_AA_table.csv',
		ref_counts=variant_folder+'/reference_AA_table.csv'
	output:
		UMI_counts=output_folder+'/counts/all_AA_counts.pkl'
	script:
		'scripts/get_UMIs.py'


rule prevalences_by_mutation:
	'''
	generates a list of which samples have alternate and coverage levels greater
	than some threshold for a given mutation (thresholds are baked into
	mutation names). Thresholds are coverage UMI count and alternate UMI count.
	Format of mutation names is mutation_region_cov_alt
	'''
	input:
		UMI_counts=output_folder+'/counts/all_AA_counts.pkl',
		metadata=config['metadata_sheet'],
	output:
		coverage_files=expand(output_folder+'/prevalences/{mutation}_cov_samples.txt', mutation=config['mutations']),
		alternate_files=expand(output_folder+'/prevalences/{mutation}_alt_samples.txt', mutation=config['mutations'])
	script:
		'scripts/get_prevalences.py'

rule get_threshold_prevalences:
	'''
	gets all mutations that pass a given set of thresholds, and the samples that
	pass these thresholds for each mutation. Adds an additional threshold - in
	addition to a minimum coverage and alternate UMI count, this code introduces
	a minimum count of number of samples that the mutation is observed in.
	Format is searchterm_cov_alt_count (where searchterm is e.g. a gene to
	search within.)
	'''
	input:
		UMI_counts=output_folder+'/counts/all_AA_counts.pkl',
		metadata=config['metadata_sheet']
	params:
		output_folder=output_folder
	output:
		filtered_mutations=output_folder+'/prevalences_by_threshold/{threshold}_mutation_list.txt'
	script:
		'scripts/get_threshold_prevalences.py'

rule check_background_mutations:
	'''
	checks to see if any Asian background K13 mutations exist
	'''
	input:
		background_mutations=expand(output_folder+'/prevalences/{mutation}_alt_samples.txt', mutation=config['background_mutations']),
		foreground_mutation=output_folder+'/prevalences/k13-Arg561His_Region:all_3_1_alt_samples.txt'
	output:
		background_mutations=output_folder+'/background_mutations/561_on_Asian_backgrounds.tsv'
	script:
		'scripts/check_background_mutations.py'

rule make_table_named_mutations:
	input:
		desired_files=expand(output_folder+'/prevalences/{file}_cov_samples.txt', file=config['mutations']),
		metadata_sheet=config['metadata_sheet']
	params:
		heirarchy=config['heirarchy']
	output:
		summary=output_folder+'/summaries/{region}:{region_type}_{cov}_{alt}_summary.tsv'
	script:
		'scripts/make_table.py'

rule parse_NFD:
	'''
	NFD haplotype is a mixture of reference alleles (N and D) and mutant (F).
	Looking for samples that have non-zero reference counts for N and D, and
	non-zero alternate counts for F
	'''
	input:
		UMI_counts=output_folder+'/counts/all_AA_counts.pkl',
		metadata=config['metadata_sheet'],
		F_alt=output_folder+'/prevalences/mdr1-Tyr184Phe_Region:all_3_1_alt_samples.txt',
		F_cov=output_folder+'/prevalences/mdr1-Tyr184Phe_Region:all_3_1_cov_samples.txt'
	params:
		N_mutation='mdr1-Asn86_3_1',
		D_mutation='mdr1-Asp1246_3_1'
	output:
		N_ref=output_folder+'/prevalences/mdr1-Asn86_3_1_ref.txt',
		D_ref=output_folder+'/prevalences/mdr1-Asp1246_3_1_ref.txt',
		N_cov=output_folder+'/prevalences/mdr1-Asn86_3_1_cov.txt',
		D_cov=output_folder+'/prevalences/mdr1-Asp1246_3_1_cov.txt',
		NFD_cov=output_folder+'/prevalences/mdr1-NFD_3_1_cov.txt',
		NFD_alt=output_folder+'/prevalences/mdr1-NFD_3_1_alt.txt'
	script:
		'scripts/parse_NFD.py'


rule make_table_threshold_prevalences:
	input:
		threshold_mutations=expand(output_folder+'/prevalences_by_threshold/{threshold}_mutation_list.txt', threshold=config['prevalence_thresholds']),
		metadata_sheet=config['metadata_sheet']
	params:
		heirarchy=config['heirarchy'],
		output_folder=config['output_folder']
	output:
		summary=output_folder+'/summaries/{cov}_{alt}_{sample_count}_threshold-summary.tsv'
	script:
		'scripts/make_table_threshold_prevalences.py'
'''

rule intersect_samples:
	input:
		file_list=expand('prevalences/{file}_cov_samples.txt', file=config['mutations'])
	params:
		file_list=config['intersections']
		
	output:
'''
