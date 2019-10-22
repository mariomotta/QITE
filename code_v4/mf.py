import numpy as np
import itertools
from   scipy import optimize as opt
from   pauli            import pauli_action
from   binary_functions import Bas2Int,Int2Bas
from   numpy            import linalg as LA
from   scipy            import linalg as SciLA
from   tools            import print_state,fidelity,dgr
from   pauli            import sigma_matrices

def mf_state(theta_):
 nbit = theta_.shape[0]
 chi  = np.zeros((nbit,2))
 for i in range(nbit):
  chi[i,0] = np.cos(theta_[i])
  chi[i,1] = np.sin(theta_[i])
 # -----
 N   = 2**nbit
 phi = np.zeros(N,dtype=complex)
 for i in range(N):
  x      = Int2Bas(i,2,nbit)
  phi[i] = np.prod([chi[k,x[k]] for k in range(nbit)])
 # -----
 return phi

def mf_energy(theta_,H_):
 nbit = theta_.shape[0]
 chi  = np.zeros((nbit,2))
 for i in range(nbit):
  chi[i,0] = np.cos(theta_[i])
  chi[i,1] = np.sin(theta_[i])
 # -----
 xjm = np.einsum('ja,abc,jb->jc',chi,sigma_matrices,chi)
 xjm = np.real(xjm)
 ea = 0.0
 for (A,h,imp,gmp) in H_:
  nact = len(A)
  for m in np.where(np.abs(h)>1e-8)[0]:
   xm  = Int2Bas(m,4,nact)
   ea += h[m]*np.prod([xjm[A[k],xm[k]] for k in range(nact)])  
 return ea

def mf_solution(theta0_,H_):
 res = opt.minimize(mf_energy,theta0_,args=H_,method='SLSQP')
 return res.x,mf_energy(res.x,H_)

# ------------------------------------------------- #

def hom_mf_state(theta_,nbit_):
 chi  = np.zeros(2)
 chi[0] = np.cos(theta_)
 chi[1] = np.sin(theta_)
 N   = 2**nbit_
 phi = np.zeros(N,dtype=complex)
 for i in range(N):
  x      = Int2Bas(i,2,nbit_)
  phi[i] = np.prod([chi[x[k]] for k in range(nbit_)])
 return phi

def hom_mf_energy(theta,nbit,H_):
 chi  = np.zeros(2)
 chi[0] = np.cos(theta)
 chi[1] = np.sin(theta)
 xjm = np.einsum('a,abc,b->c',chi,sigma_matrices,chi)
 xjm = np.real(xjm)
 ea = 0.0
 for (A,h,imp,gmp) in H_:
  nact = len(A)
  for m in np.where(np.abs(h)>1e-8)[0]:
   xm  = Int2Bas(m,4,nact)
   ea += h[m]*np.prod([xjm[xm[k]] for k in range(nact)])
 return ea

def hom_mf_solution(theta0_,nbit_,H_):
 res = opt.minimize(hom_mf_energy,theta0_,args=(nbit_,H_),method='SLSQP')
 return res.x,hom_mf_energy(res.x,nbit_,H_)
