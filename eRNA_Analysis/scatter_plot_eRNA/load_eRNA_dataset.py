
import csv
import numpy as np
# this file has one method, it returns the dataset


def return_dataset_matrix(file_name):

	dataset = []
	sample_names = []
	gene_names = []
	row_numb = 0 
	with open(file_name, 'rb') as file1:
			reader = csv.reader(file1, delimiter='\t')
			header = True

			for row in reader:

				if header == True:
					sample_names = row
					header = False
					continue

				gene_names.append(row[0])

				this_gene = []
				for i in range(1,len(row)):
					this_gene.append(float(row[i]))
				dataset.append(this_gene)



				#if row_numb > 5:
				#	break
				row_numb += 1

	dataset = np.array(dataset)
	dataset = dataset.T

	#print dataset.shape
	#print len(sample_names)
	#print len(gene_names)

	#rows are samples, columns are expressions
	return dataset, sample_names, gene_names