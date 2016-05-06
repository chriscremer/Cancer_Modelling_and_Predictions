

import csv

clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt' 




with open(clinical_file, 'rU') as f:
	reader = csv.reader(f, delimiter='\t')
	count = 0
	header = []
	for row in reader:

		#print count
		if count == 0:
			header = row
		elif count > 2:
			for i in range(len(row)):
				print str(i) + '  ' + str(header[i]) + '  ' + str(row[i])

		print

		count += 1