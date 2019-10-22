import numpy as np
import itertools
from   scipy import optimize as opt
from   numpy import linalg   as LA
from   scipy import linalg   as SciLA

from   pauli            import pauli_action
from   binary_functions import Bas2Int,Int2Bas
from   tools            import print_state,fidelity,dump_state,read_state,dump_lanz_vecs
from   hamiltonian      import Hmat,Hmoms,Hpsi,Hii

def ITE_FCI(H_,db,bmax,psi0=None,omega=None):

 Hm    = Hmat(H_)
 N     = Hm.shape[0]
 nbit  = int(np.log2(N))
 eps,U = SciLA.eigh(Hm)
 m0    = np.argmin(eps)
 zeta  = np.exp(-db*(eps-eps[m0]))

 fout = open('ITE_FCI.out','w')
 fout.write("FCI gs energy %.6f \n" % eps[m0])
 fout.write("FCI gs wfn \n")
 print_state(U[:,m0],nbit,fout)

 if(psi0 is None):
  i0          = np.argmin(np.diag(Hm))
  psi_FCI     = np.zeros(N,dtype=complex)
  psi_FCI[i0] = 1.0
 else:
  psi_FCI     = psi0.copy()

 nbeta    = int(bmax/db)+1
 fout.write("FCI ITE\n")
 for ib in range(nbeta):
  ea,ev    = Hmoms(H_,psi_FCI)
  print ib*db,ea,ev
  psi_FCI  = np.dot(np.conj(U.T),psi_FCI)
  psi_FCI  = zeta*psi_FCI
  psi_FCI  = np.dot(U,psi_FCI)
  psi_FCI /= LA.norm(psi_FCI)
  if(omega is None): fide = fidelity(psi_FCI,U[:,m0])
  else:              fide = LA.norm(psi_FCI[omega])**2
  fout.write("%.6f %.6f %.6f %.6f \n" % (ib*db,ea,ev,fide))

 fout.write("FCI ITE gs wfn \n")
 print_state(psi_FCI,nbit,fout)

 fout.close()

# ------------------------------------------------- #

def ExpmbH(H_,psi_,db):
 phi = psi_.copy()
 chi = psi_.copy()
 i   = 0
 while(LA.norm(chi)>1e-6):
  chi  = (-db/float(i+1))*Hpsi(H_,chi)
  phi += chi.copy()
  i   += 1
 nu = LA.norm(phi)
 return phi.copy()/nu,nu

def Lanczos_kernel(Hmat_,chi_):
 # to match the QITE implementation
 chi_    = chi_[:,::2]
 Sact    = np.einsum('ic,id->cd'   ,np.conj(chi_),chi_)
 # regularization
 Sact   += 1e-8*np.eye(chi_.shape[1])
 Hact    = np.einsum('ic,ij,jd->cd',np.conj(chi_),Hmat_,chi_)
 eps_,c_ = SciLA.eig(Hact,Sact)
 m0      = np.argmin(eps_)
 eta_    = np.einsum('c,ic->i',c_[:,m0],chi_)
 eta_   /= LA.norm(eta_)
 return eta_

# ------------------------------------------------- #

def ITE(H_,db,bmax,lanczos=False,psi0=None):

 N     = H_[0][2].shape[1]
 nbit  = int(np.log2(N))
 hdiag = np.zeros(N,dtype=complex)
 for i in range(N):
  hdiag[i] = Hii(H_,i)

 precond = lambda x,e, *args: x/(hdiag-e+1e-4)

 def hop(c_):
  return Hpsi(H_,c_)

 if(psi0 is None):
  i0       = np.argmin(hdiag)
  psi0     = np.zeros(N,dtype=complex)
  psi0[i0] = 1.0
 
 from pyscf.lib import davidson
 epsm0,Um0 = davidson(hop,psi0,precond)
 
 fout = open('ITE.out','w')
 fout.write("FCI gs energy %.6f \n" % epsm0)
 fout.write("FCI gs wfn \n")
 print_state(Um0,nbit,fout)

 if(psi0 is None):
  i0          = np.argmin(hdiag)
  psi_ITE     = np.zeros(N,dtype=complex)
  psi_ITE[i0] = 1.0
 else:
  psi_ITE     = psi0.copy()

 nbeta = int(bmax/db)+1
 hvect_LANZ = np.zeros(nbeta+1)
 svect_LANZ = np.zeros(nbeta+1)

 if(lanczos):
  space_LANZ = np.zeros((N,nbeta),dtype=complex)

 fout.write("ITE\n")
 for ib in range(nbeta):
  ea,ev          = Hmoms(H_,psi_ITE)
  hvect_LANZ[ib] = ea
  fide           = fidelity(psi_ITE,Um0)

  if(lanczos):
   space_LANZ[:,ib] = psi_ITE.copy()
   psi_LANZ         = Lanczos_kernel(Hm,space_LANZ[:,:ib+1])
   ea_,ev_          = Hmoms(H_,psi_LANZ)
   fide_            = fidelity(psi_LANZ,Um0)
   fout.write("%.6f %.6f %.6f %.6f %.6f %.6f %.6f \n" % (ib*db,ea,ev,fide,ea_,ev_,fide_))
  else:
   fout.write("%.6f %.6f %.6f %.6f \n" % (ib*db,ea,ev,fide)) 

  psi_ITE,dn = ExpmbH(H_,psi_ITE,db)
  svect_LANZ[ib+1]  = svect_LANZ[ib]+np.log(dn)

 fout.write("ITE gs wfn \n")
 print_state(psi_ITE,nbit,fout)
 dump_state(psi_ITE,nbit,'ite.psi')
 dump_lanz_vecs(hvect_LANZ[:nbeta],svect_LANZ[:nbeta],'qlanz.vecs')

 fout.close()

# ------------------------------------------------- #
