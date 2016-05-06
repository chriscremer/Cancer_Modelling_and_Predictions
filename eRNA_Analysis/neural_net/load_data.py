# We'll need the csv module to read the file
import csv
# We'll need numpy to manage arrays of data
import numpy as np

# We'll need the DenseDesignMatrix class to return the data
from pylearn2.datasets.dense_design_matrix import DenseDesignMatrix

def load_data(start, stop):
    """
    Loads the eRNA dataset

    The dataset contains 41 samples (18 LG, 23 HG).  4320 eRNAs per sample.

    Parameters
    ----------
    start: int
    stop: int

    Returns
    -------

    dataset : DenseDesignMatrix
        A dataset include examples start (inclusive) through stop (exclusive).
        The start and stop parameters are useful for splitting the data into
        train, validation, and test data.
    """
    X = []
    y = []
    sample_names = []
    gene_names = []
    with open('/data1/morrislab/ccremer/eRNA/bc_eRNA_1k_noOverlap_41Samples_09.txt', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = True
        for row in reader:
            # Skip the first row containing the string names of each attribute
            if header:
                sample_names = row
                header = False
                continue
            gene_names.append(row[0])
            # Convert the row into numbers
            row_values = [float(elem) for elem in row[1:]]
            X.append(row_values)

    X = np.asarray(X)
    X = X.T


    set3_LG = ['BC060','BC005','BC030','BC026','BC050',
            'BC006','BC019','BC052','BC040','BC037',
            'BC032','BC039','BC041','BC029','BC024',
            'BC020','BC038','BC054','BC053','BC014',
            'BC059','BC051', 'BC063']

    set3_HG = ['BC048','BC046','BC012','BC061','BC010',
                'BC013','BC035','BC027','BC036','BC008',
                'BC031','BC064','BC034','BC042','BC044',
                'BC002','BC017','BC022']

    y = []
    for samp_name in sample_names:
        if samp_name in set3_LG:
            y.append(0.0)
        else:
            y.append(1.0)

    #need to order the samples in alternating order
    #so that the training is balanced
    y_ordered = []
    index_order = []
    used_index = []
    current_value = 0.0
    for i in range(len(y)):
        while len(y) != len(y_ordered):
            for j in range(len(y)):
                if y[j] == current_value and j not in used_index:
                    #add it to new list and record the index in the original list
                    y_ordered.append(current_value)
                    index_order.append(j)
                    #switch current value to the other 
                    if current_value == 0.0:
                        current_value = 1.0
                    else:
                        current_value = 0.0
                    #take it out of the original list
                    #y.pop(j)
                    used_index.append(j)
                    break
                #if last element, just take rest
                if j == (len(y)-1):
                    for k in range(len(y)):
                        if k not in used_index:
                            #add it to new list and record the index in the original list
                            y_ordered.append(y[k])
                            index_order.append(k)
                            #take it out of the original list
                            #y.pop(j)
                            used_index.append(k)
                    break


    X = X[index_order]
    y = np.asarray(y_ordered)
    y = y.reshape(y.shape[0], 1)

    # print X[0][0]
    # print y[0]

    # print X.shape
    # print y.shape


    X = X[start:stop, :]
    y = y[start:stop, :]


    return DenseDesignMatrix(X=X, y=y)


if __name__ == "__main__":

    load_data(0, 10)