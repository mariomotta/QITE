import numpy as np
import itertools
from   scipy import optimize as opt
from   numpy import linalg   as LA
from   scipy import linalg   as SciLA

from   pauli            import pauli_action,sigma_matrices
from   binary_functions import Bas2Int,Int2Bas
from   tools            import print_state,fidelity,dump_state,read_state,dump_lanz_rte
from   hamiltonian      import Hmat,Hmoms,Hpsi
from   pyscf.lib        import davidson

# ------------------------------------------------- #

def ExpitH(H_,psi_,dT):
 phi = psi_.copy()
 chi = psi_.copy()
 i   = 0
 while(LA.norm(chi)>1e-6):
  chi  = (1j*dT/float(i+1))*Hpsi(H_,chi)
  phi += chi.copy()
  i   += 1
 nu = LA.norm(phi)
 return phi.copy()/nu,nu

# ------------------------------------------------- #

def RTE(H_,dT,Tmax,lanczos=False,psi0=None):
 # --- diagonalization ---
 N     = H_[0][2].shape[1]
 nbit  = int(np.log2(N))
 hdiag = np.zeros(N,dtype=complex)
 ei    = np.zeros(N,dtype=complex)
 for i in range(N):
  xi       = Int2Bas(i,2,nbit)
  for (A,h,imp,gmp) in H_:
   nact=len(A)
   for m in np.where(np.abs(h)>1e-8)[0]:
    sm  = Int2Bas(m,4,nact)
    #print A,sm,[xi[A[i]] for i in range(nact)]
    smx = [ sigma_matrices[xi[A[w]],xi[A[w]],sm[w]] for w in range(nact)]
    hdiag[i] += h[m]*np.prod(smx)
  if(i%1000==0): print i,N,hdiag[i]
 
 precond = lambda x,e, *args: x/(hdiag-e+1e-4)
 def hop(c_):
  return Hpsi(H_,c_)

 epsm0,Um0 = davidson(hop,psi0,precond)

 fout = open('RTE_davidson.out','w')
 fout.write("gs energy %.6f \n" % epsm0)
 
 # --- initial state ---
 if(psi0 is None):
  i0       = np.argmin(hdiag)
  psi0     = np.zeros(N,dtype=complex)
  psi0[i0] = 1.0

 # --- real-time evolution ---
 bra_RTE  = psi0[:]
 braH_RTE = Hpsi(H_,psi0[:])[:]
 ket_RTE  = psi0[:]

 nbeta = int(Tmax/dT)+1
 hvect_LANZ = np.zeros(nbeta+1,dtype=complex)
 svect_LANZ = np.zeros(nbeta+1,dtype=complex)

 fout.write("ITE\n")
 for ib in range(nbeta):
  hvect_LANZ[ib] = np.einsum('a,a',np.conj(braH_RTE),ket_RTE)
  svect_LANZ[ib] = np.einsum('a,a',np.conj(bra_RTE),ket_RTE)
  ket_RTE        = ExpitH(H_,ket_RTE,dT)[0]
  print ib,hvect_LANZ[ib]

 dump_lanz_rte(hvect_LANZ[:nbeta],svect_LANZ[:nbeta],'qlanz.vecs')

 fout.close()

# ------------------------------------------------- #
