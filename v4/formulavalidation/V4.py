# This script contains the validation of the formula given
import numpy as np
from scipy.io import mmread
# Importing matrix
coo_matrix =  mmread('s12.mtx')
# Conver coo to dense np array
matrix = coo_matrix.todense()
# c3 formula
nodes = np.matmul(np.multiply(matrix,np.matmul(matrix,matrix)),np.ones(np.shape(matrix)[1])*0.5)
# Print result
print (nodes)
