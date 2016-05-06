

import numpy as np

#add parent directory to path to get my packages there
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

#code for loading dataset
from load_eRNA_dataset import return_dataset_matrix

from calc_p_value import p_value

def get_labels(sample_names):


	set3_LG = ['BC060','BC005','BC030','BC026','BC050',
			'BC006','BC019','BC052','BC040','BC037',
			'BC032','BC039','BC041','BC029','BC024',
			'BC020','BC038','BC054','BC053','BC014',
			'BC059','BC051', 'BC063']

	set3_HG = ['BC048','BC046','BC012','BC061','BC010',
				'BC013','BC035','BC027','BC036','BC008',
				'BC031','BC064','BC034','BC042','BC044',
				'BC002','BC017','BC022']


	labels = []
	for samp_name in sample_names:
		if samp_name in set3_LG:
			labels.append(0.0)
		else:
			labels.append(1.0)

	labels = np.array(labels)

	return labels


if __name__ == "__main__":

	dataset, sample_names, gene_names = return_dataset_matrix('/data1/morrislab/ccremer/eRNA/bc_eRNA_1k_noOverlap_41Samples_09.txt')
	
	labels = get_labels(sample_names)

	p_value(dataset, labels, gene_names)

