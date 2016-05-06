

import csv

import sys


if len(sys.argv) > 1:

	if sys.argv[1] == 'cqcf':
		#quality form
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_cqcf_brca.txt'
		print 'File is ', clinical_file

	elif sys.argv[1] == 'nte':
		#new tumor event 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_nte_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'drug':
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_drug_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'patient': 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'radiation': 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_radiation_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'omf': 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_omf_v4.0_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'nte_follow': 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_follow_up_v4.0_nte_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'follow_v4': 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'follow_v2': 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt'
		print 'File is ', clinical_file 

	elif sys.argv[1] == 'follow_v1': 
		clinical_file = '/data1/morrislab/ccremer/TCGA_data/BRCA_batch61/Clinical/Biotab/nationwidechildrens.org_clinical_follow_up_v1.5_brca.txt'
		print 'File is ', clinical_file 

	else:
		print 'Not sure what file you want to read. Options are: cqcf, nte, patient, drug, radiation, omf, nte_follow, follow_v4, follow_v2, follow_v1'


print 'Options are: cqcf, nte, patient, drug, radiation, omf, nte_follow, follow_v4, follow_v2, follow_v1'

with open(clinical_file, 'rU') as f:
	reader = csv.reader(f, delimiter='\t')
	count = 0
	header = []
	for row in reader:

		print count
		if count == 0:
			header = row
		elif count > 2:
			for i in range(len(row)):
				print str(i) + '  ' + str(header[i]) + '  ' + str(row[i])

		print

		count += 1