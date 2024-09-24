'''
special purpose script to pull metadata for all non-control non-Geita samples
'''

metadata='/nfs/jbailey5/baileyweb/asimkin/MSMT_561H_paper/preprint_metadata-all_Dativa_metadata.tsv'
IRNL_samples=[line.strip()+'-2021-MSMT-1' for line in open('IRNL_samples.txt')]
#R561H_samples=set([line.strip().split('\t')[0] for line in open('/nfs/jbailey5/baileyweb/asimkin/MSMT_561H_paper/R561_cov_3.tsv')][6:])
#attempted_samples=[line.rstrip().split('\t') for line in open('/nfs/jbailey5/baileyweb/asimkin/MSMT_561H_paper/DR2_wrangler_output/mip_ids/allMipsSamplesNames.tab.txt')][1:]
#attempted_samples=set([line[1] for line in attempted_samples if len(line)>1])
#successful_samples=set([line.strip().split(',')[0] for line in open('/nfs/jbailey5/baileyweb/asimkin/MSMT_561H_paper/variant/k13_results/alternate_AA_table.csv')][6:])

regions={'DS': 'Dar es salaam', 'DO': 'Dodoma', 'KG': 'Kagera', 'TA': 'Tabora', 'KL': 'Kilimanjaro', 'MY': 'Manyara', 'MA': 'Mara', 'KI': 'Kigoma', 'RU': 'Ruvuma', 'MT': 'Mtwara', 'NJ': 'Njombe', 'SO': 'Songwe', 'TB': 'Tabora'}
community_regions={'LUN-': 'Ruvuma', 'MAB-': 'Tanga', 'MGD-': 'Tanga', 'MPY-': 'Tanga', 'NYA-': 'Kigoma'}

#email from Filbert, 11-4-23
#KMW should be KLMW, Kilimanjaro, Mwanga district
#DMP should be DOMP, Dodoma Mpwapwa district
#DOM should be DOMP, Dodoma Mpwapwa district
#MYK should be MYKI, Manyara Kiteto Dc district
typo_regions={'KMW-': 'Kilimanjaro', 'DMP-':'Dodoma', 'DOM-':'Dodoma', 'MYK-':'Manyara'}

def clean_samples(samples):
	community_samples=set([])
	facility_samples=set([])
	for sample in samples:
		prefix=sample.split('-')[0]
		if not prefix.isdigit() and (len(prefix)==3 or len(prefix)==4):
			if sample[:4] in community_regions:
				community_samples.add(sample)
			elif sample[:4] in typo_regions:
				facility_samples.add(sample)
			elif sample[:2] in regions:
				facility_samples.add(sample)
	return facility_samples, community_samples

def get_info(sample_set, output_path):
	output_file=open(output_path, 'w')
	found_samples=set([])
	for line_number, line in enumerate(open(metadata)):
		split_line=line.strip().split('\t')
		if line_number==0:
			output_file.write(line)
		if split_line[8] in sample_set or split_line[11] in sample_set:
			found_samples.add(split_line[11])
			output_file.write(line)
	print('original was', len(sample_set), len(found_samples))
	missing_data=sample_set-found_samples
	print('missing is', len(missing_data))


#attempted_facility, attempted_community=clean_samples(attempted_samples)
#successful_facility, successful_community=clean_samples(successful_samples)
IRNL_facility, IRNL_community=clean_samples(IRNL_samples)
print(len(IRNL_facility), len(IRNL_community))
get_info(IRNL_facility|IRNL_community, 'IRNL_coverage_3_alt_1_metadata.tsv')

#R561_facility, R561_community=clean_samples(R561H_samples)
#print(len(R561_facility), len(R561_community))
#get_info(R561_facility, 'R561_coverage_3_facility_metadata')
#get_info(R561_community, 'R561_coverage_3_community_metadata')
#get_info(R561_community|R561_facility, 'R561_coverage_3_metadata')


#print(len(attempted_facility), len(attempted_community))
#print(len(successful_facility), len(successful_community))
#get_info(attempted_facility, 'attempted_facility_metadata.tsv')
#get_info(attempted_community, 'attempted_community_metadata.tsv')
#get_info(successful_facility, 'successful_facility_metadata.tsv')
#get_info(successful_community, 'successful_community_metadata.tsv')
