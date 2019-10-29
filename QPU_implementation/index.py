import numpy as np

# To keep track of Lie algebra. Let P represent some Pauli operator. We want to know PiPj = cijPij.

idx = np.zeros([4,4],dtype=int)
idx[0,0] = 0
idx[0,1] = 1
idx[0,2] = 2
idx[0,3] = 3
idx[1,0] = 1
idx[1,1] = 0
idx[1,2] = 3
idx[1,3] = 2
idx[2,0] = 2
idx[2,1] = 3
idx[2,2] = 0
idx[2,3] = 1
idx[3,0] = 3
idx[3,1] = 2
idx[3,2] = 1
idx[3,3] = 0

coeff = np.zeros([4,4],dtype=complex)
coeff[0,0] = 1
coeff[0,1] = 1
coeff[0,2] = 1
coeff[0,3] = 1
coeff[1,0] = 1
coeff[1,1] = 1
coeff[1,2] = 1j
coeff[1,3] = -1j
coeff[2,0] = 1
coeff[2,1] = -1j
coeff[2,2] = 1
coeff[2,3] = 1j
coeff[3,0] = 1
coeff[3,1] = 1j
coeff[3,2] = -1j
coeff[3,3] = 1