'''
filters out metadata of samples that were not part of this study (Geita samples
and samples that we didn't attempt to wrangle)
'''
valid_regions=set(['Njombe', 'Dar Es Salaam', 'Ruvuma', 'Kagera', 'Mara', 'Dodoma', 'Kilimanjaro', 'Manyara', 'Tanga', 'Mtwara', 'Kigoma', 'Songwe', 'Tabora'])
valid_metadata=open('../manuscript_metadata_attempted_samples.csv', 'w')
attempted_samples=[line.rstrip().split('\t') for line in open('/nfs/jbailey5/baileyweb/asimkin/MSMT_561H_paper/DR2_wrangler_output/mip_ids/allMipsSamplesNames.tab.txt')][1:]
attempted_samples=set([line[1] for line in attempted_samples if len(line)>1])
attempted_valid_samples=open('attempted_samples.txt', 'w')

meta_samples=set([])
for line_number, line in enumerate(open('../../../manuscript_metadata.csv')):
	split_line=line.strip().split(',')
	region=split_line[0]
	sample=split_line[11]
	if line_number==0:
		valid_metadata.write(line)
	if region in valid_regions:
		if sample in attempted_samples:
			attempted_valid_samples.write(sample+'\n')
			valid_metadata.write(line)
			
