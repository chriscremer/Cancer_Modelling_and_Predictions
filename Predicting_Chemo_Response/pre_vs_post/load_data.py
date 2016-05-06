

import xlrd
import csv
import pprint
import numpy as np
import random


def get_all_expressions(expr_file, first_row, samples):

	columns = []
	name_order = []
	cur_column = 0
	for samp in first_row:
		if samp in samples:
			columns.append(cur_column)
			name_order.append(samp)
		cur_column += 1

	#gene are rows, samp are columns
	data = []
	with open(expr_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		row_count = 0
		for row in reader:

			if row_count == 0:
				row_count += 1
				continue

			#if row[len(row)-1] == 'BCL2':
			this_gene = []
			for i in columns:
				this_gene.append(float(row[i]))
			data.append(this_gene)

			#for j in row:
			#	print row


			row_count += 1
			#if row_count > 40000:
			#	break

	#pprint.pprint(training_data)

	print 'Numb of genes ' + str(len(data))
	print 'Numb of samples ' + str(len(data[0]))


	data = np.array(data)
	#samples will be rows and genes are columns
	data = data.T

	return data, name_order



def assemble_data(expr_file, clinical_file, translation_file):

	#################################################
	#extracting data from xlsx files into lists
	################################################
	sample_translation = []
	workbook = xlrd.open_workbook(translation_file)
	worksheet = workbook.sheet_by_name('Sheet1')
	#numb of rows in the file
	num_rows = worksheet.nrows - 1
	curr_row = -1
	num_cells = worksheet.ncols - 1
	#go through all rows
	while curr_row < num_rows:
		curr_row += 1
		row = worksheet.row(curr_row)
		#print 'Row:', curr_row
		curr_cell = -1
		this_row = []
		#go through all cells of the row
		while curr_cell < num_cells:
			curr_cell += 1
			# Cell Types: 0=Empty, 1=Text, 2=Number, 3=Date, 4=Boolean, 5=Error, 6=Blank
			cell_type = worksheet.cell_type(curr_row, curr_cell)
			cell_value = worksheet.cell_value(curr_row, curr_cell)
			#print '	', cell_type, ':', cell_value
			this_row.append(str(cell_value))
			#print cell_value
		#this has all the data in the file
		sample_translation.append(this_row)

	#now do the same as translation file
	#save all data into this list
	clinpathdata = []
	workbook = xlrd.open_workbook(clinical_file)
	worksheet = workbook.sheet_by_name('cohort.data')
	num_rows = worksheet.nrows - 1
	curr_row = 4
	num_cells = worksheet.ncols - 1
	while curr_row < num_rows:
		curr_row += 1
		row = worksheet.row(curr_row)
		#print 'Row:', curr_row
		curr_cell = -1
		this_row = []
		while curr_cell < num_cells:
			curr_cell += 1
			# Cell Types: 0=Empty, 1=Text, 2=Number, 3=Date, 4=Boolean, 5=Error, 6=Blank
			cell_type = worksheet.cell_type(curr_row, curr_cell)
			cell_value = worksheet.cell_value(curr_row, curr_cell)
			#print '	', cell_type, ':', cell_value
			this_row.append(str(cell_value))
			#print cell_value

		clinpathdata.append(this_row)

	#print clinpathdata


	#################################################
	#list of the order of the samps in the raw data file
	################################################
	first_row = []
	with open(expr_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in reader:
			first_row = row
			break

	#################################################
	#selecting samples that have TURBT and that arent red and have RC <-resection I think
	#################################################
	red_samp_ids = ['BER104', 'BER096', 'BER062', 'BER027', 'BER020', 'BER011', 'BER010','BER026','BER028','BER077','BER092','BER105','BER115',]

	TURBT_sample_ids = []
	RC_sample_ids = []

	header = True
	for row in sample_translation:
		#skip header
		if header:
			header = False
			continue
		#if red, skip it
		if row[0] in red_samp_ids:
			continue

		patient_numb = int(row[1].split('-')[1])

		#see if the patient has a TURBT and RC
		TURBT = False
		RC = False
		header2 = True
		for row2 in sample_translation:
			if header2:
				header2 = False
				continue
			if int(row2[1].split('-')[1]) == patient_numb and row2[0] not in red_samp_ids and row2[0] not in TURBT_sample_ids and row2[0] not in RC_sample_ids:
				if row2[2] == 'TURBT':
					TURBT = True
					TURBT_samp_id = row2[0]
				if row2[2] == 'RC':
					RC = True
					RC_samp_id = row2[0]
		if TURBT and RC:
			TURBT_sample_ids.append(TURBT_samp_id)
			RC_sample_ids.append(RC_samp_id)



	#print TURBT_sample_ids
	#print RC_sample_ids

	TURBT_data, TURBT_sample_order = get_all_expressions(expr_file, first_row, TURBT_sample_ids)
	RC_data, RC_sample_order = get_all_expressions(expr_file, first_row, RC_sample_ids)


	return TURBT_data, TURBT_sample_order, RC_data, RC_sample_order



def get_gene_names(expr_file):

	gene_names = []

	with open(expr_file, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',', quotechar='|')
		row_count = 0
		for row in reader:

			if row_count == 0:
				row_count += 1
				continue

			#gene_names.append((row[0], row[len(row)-1]))
			gene_names.append(row[len(row)-1])

			row_count += 1

	#print len(gene_names)
	#print gene_names[:1000]

	return gene_names


if __name__ == "__main__":

	TURBT_data, TURBT_sample_order, RC_data, RC_sample_order = assemble_data('/data1/morrislab/ccremer/chemo_response/rawdata_main_genes_only.csv', '/data1/morrislab/ccremer/chemo_response/ClinpathDataBernCohort.xlsx', '/data1/morrislab/ccremer/chemo_response/SampleTranslation.xlsx')

	print TURBT_data.shape
	print RC_data.shape