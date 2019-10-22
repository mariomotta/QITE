import numpy as np
import scipy
from   binary_functions  import Int2Bas,Bas2Int,Opp2Str,Psi2Str,Lst2Str

# ---------------------------------------------------------- #

pauli_product=[np.zeros((4,4),dtype=int),np.zeros((4,4),dtype=complex)]

pauli_product[0][0,:]=[0,1,2,3]
pauli_product[0][1,:]=[1,0,3,2]
pauli_product[0][2,:]=[2,3,0,1]
pauli_product[0][3,:]=[3,2,1,0]
pauli_product[1][0,:]=[1,  1,  1,  1]
pauli_product[1][1,:]=[1,  1, 1j,-1j]
pauli_product[1][2,:]=[1,-1j,  1, 1j]
pauli_product[1][3,:]=[1, 1j,-1j,  1]

sigma_matrices=np.zeros((2,2,4),dtype=complex)
for i in range(2): 
 j = (i+1)%2
 sigma_matrices[i,i,0] = 1.0
 sigma_matrices[i,j,1] = 1.0
 sigma_matrices[i,j,2] = 1.0j*(-1.0)**(i+1.0)
 sigma_matrices[i,i,3] = (-1.0)**i

d12  = lambda t: 1 if t%3>0 else 0
d12f = np.vectorize(d12)
d2   = lambda t: 1 if t==2 else 0
d2f  = np.vectorize(d2)
d23  = lambda t: 1 if t>1 else 0
d23f = np.vectorize(d23)

# ---------------------------------------------------------- #

def computational_basis(nbit_):
 N=2**nbit_
 for i in range(N):
  print i,Psi2Str(Int2Bas(i,2,nbit_))

def pauli_basis(nbit_):
 M=4**nbit_
 for i in range(M):
  print i,Opp2Str(Int2Bas(i,4,nbit_))

# ---------------------------------------------------------- #

def pauli_action(active_,nbit_,verbose=False):
 nact = len(active_)
 N    = 2**nbit_
 M    = 4**nact

 dot    = [2**(nbit_-1-i) for i in range(nbit_)]
 ind_sx = np.zeros((M,N),dtype=int)
 gmm_sx = np.zeros((M,N),dtype=complex)+1

 svec = np.zeros((M,nbit_),dtype=int) 
 for mu in range(M):
  svec[mu,active_]=Int2Bas(mu,4,nact)
 sxyvec = d12f(svec)
 nyvec  = d2f(svec)
 syzvec = d23f(svec)
 nyvec  = np.einsum('ab->a',nyvec)

 xvec = np.zeros((N,nbit_),dtype=int)
 for xi in range(N):
  xvec[xi,:] = np.asarray(Int2Bas(xi,2,nbit_))

 gmm_sx=np.einsum('am,bm->ba',xvec,syzvec)+0j
 gmm_sx[:,:]=(-1)**gmm_sx[:,:]
 for mu in range(M):
  gmm_sx[mu,:] *= 1j**nyvec[mu]
  yvec          = (xvec[:,:]+sxyvec[mu,:])%2
  ind_sx[mu,:]  = np.einsum('a,ba->b',dot,yvec)

 return ind_sx,gmm_sx 

