

#This will go through batch folders and record some clinical variables and make matrices based on the variables.


from os.path import expanduser
home = expanduser("~")

# import os,sys,inspect
# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#sys.path.insert(0,parentdir) 
# sys.path.insert(0,currentdir+'/simulated_data')

import numpy as np
import csv
import random



def assemble():

	tcga_batch_directories = ['BRCA_batch47', 'BRCA_batch56', 'BRCA_batch61', 'BRCA_batch72', 'BRCA_batch74', 'BRCA_batch93', 'BRCA_batch96']
	# tcga_batch_directories = ['BLCA_batch86', 'BLCA_batch113', 'BLCA_batch128', 'BLCA_batch150_170_175', 'BLCA_batch192_199_207']

	clinical_data = {}


	for directory in tcga_batch_directories:

		print 'Going through ' + str(directory) + ' ...'

		with open(home + '/TCGA_data/' + directory + '/file_manifest.txt') as f:
			reader = csv.reader(f, delimiter='\t')
			for row in reader:

				if 'genes.normalized_results' in row[6]:
					sample = row[4]
					clinical_data[sample] = {}

					clinical_data[sample]['batch'] = directory
					clinical_data[sample]['gene_exps_file'] = row[6]
					clinical_data[sample]['random_value'] = random.randint(0, 1)

					#Get clinical data for this sample
					#HOW DO I MATCH PATIENT TO SAMPLE?
					#manifest has sample name and clinical has patient name
					#patient barcode is beginning of sample name
					cancer_type = directory[:4].lower()
					with open(home + '/TCGA_data/' + directory + '/Clinical/Biotab/nationwidechildrens.org_clinical_patient_' + cancer_type + '.txt') as f2:
							reader2 = csv.reader(f2, delimiter='\t')
							for row2 in reader2:
								#if patient barcode in sample name
								if row2[1] in sample:
									clinical_data[sample]['vital_status'] = row2[13]
									clinical_data[sample]['gender'] = row2[6]
									clinical_data[sample]['birth_days_to'] = row2[5]
									clinical_data[sample]['tumor_status'] = row2[12]
									clinical_data[sample]['race'] = row2[8]
									clinical_data[sample]['last_contact_days_to'] = row2[14]
									break

	for samp in clinical_data:
		print 'Example. Sample ' + samp
		for key in clinical_data[samp]:
			print key, clinical_data[samp][key]
		print
		break

	print 'Number of total samples is ' + str(len(clinical_data))

	return clinical_data



def select_samples(clinical_data):
	'''
	Return a list of sample names that fit some criteria.
	'''

	samples = []
	for samp in clinical_data:
		# if clinical_data[samp]['batch'] == 'BRCA_batch47':
		# 	samples.append(samp)

		#if its a solid tumour
		if samp.endswith('01'):
			if clinical_data[samp]['vital_status'] == 'Dead':
				samples.append(samp)

	print 'Number of samples after selection criteria is ' + str(len(samples))

	return samples



def get_expressions(clinical_data, samples):
	'''
	Given the batch and names of samples that we want for this study, this function
	will get their gene expressions and return a matrix.
	'''

	print 'Assembling X ... '

	X = []
	for samp in samples:

		batch = clinical_data[samp]['batch']
		exp_file = clinical_data[samp]['gene_exps_file']

		gene_expressions = []
		with open(home + '/TCGA_data/' + batch + '/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/' + exp_file) as f2:
			b = 0
			for line in f2:
				if b == 0:
					b+=1
					continue
				firstSplit = line.split()
				secondSplit = firstSplit[0].split('|')

				if secondSplit[0] == '?':
					continue

				gene_expressions.append(float(firstSplit[1]))
		X.append(gene_expressions)

	X = np.array(X)
	print 'X shape ' + str(X.shape)

	return X



def make_target(clinical_data, sample_order):

	y = []
	for samp in sample_order:


		#RACE
		# race = clinical_data[samp]['race']
		# if race == 'WHITE':
		# 	y.append(1)
		# else:
		# 	y.append(0)


		#GENDER
		# gend = clinical_data[samp]['gender']
		# if gend == 'MALE':
		# 	y.append(1)
		# else:
		# 	y.append(0)

		#RANDOM VALUE
		y.append(clinical_data[samp]['random_value'])


	y = np.array(y)

	print 'y shape ' + str(y.shape)

	return y






def main():

	clinical_data = assemble()
	samples = select_samples(clinical_data)
	X = get_expressions(clinical_data, samples)
	y = make_target(clinical_data, samples)
	

	return X, y


if __name__ == "__main__":

	main()
