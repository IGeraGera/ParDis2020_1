import numpy as np
from scipy.io import mmwrite
from scipy.sparse import coo_matrix


# Make a symmetric upper triangle matrix
a = [[ 0,1,0,0,1,0,1],
    [ 0,0,1,0,1,1,0 ],
    [ 0,0,0,1,0,1,1],
    [ 0,0,0,0,1,0,1],
    [ 0,0,0,0,0,1,1],
    [ 0,0,0,0,0,0,1],
    [ 0,0,0,0,0,0,0]]
a=np.array(a)
a=a.T+a
# Convert it to COO
acoo=coo_matrix(a)
# Write it to a file
mmwrite('test.mtx',acoo,field='pattern')

# TEST --- TEST --- TEST
# Create the v vector

#for i in range(a.shape[0]):
#    v[i] = i+1

v = np.ones(a.shape[0])*2
ans = np.matmul(a,v)

print(ans)
