rule all:
	input:
		IRNL='IRNL_prevalences_by_region.tsv',
		IRNGEG='IRNGEG_prevalences_by_region.tsv',
		k13='k13_prevalences_by_region.tsv',
		dhfr='dhfr_prevalences_by_region.tsv',
		dhps='dhps_prevalences_by_region.tsv',
		k13_kagera='k13_kagera_prevalences_by_district.tsv'

rule get_UMIs:
	input:
		cov_counts=
		alt_counts=
		ref_counts=
	output:
		
	script:
		get_UMIs.py

rule get_metadata:
	input:
		metadata_sheet=
	output:
		metadata_file=
	script:
		'scripts/get_metadata'

rule get_prevalences:
	input:
	params:
		metadata=
		thresholds=config['thresholds'] #any way to make this a wildcard across multiple named searches?
		input_samples=
	output:
		coverage_samples=
		alternate_samples=
		formatted_prevalences=
	script:
		'scripts/get_prevalences.py'

rule intersect_samples:
	input:
		file_list=expand(config['intersections']['wildcard'])
	params:
		file_list=config['intersections']
		
	output:
