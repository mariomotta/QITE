import numpy as np
import itertools
from   pauli            import pauli_action
from   binary_functions import Bas2Int,Int2Bas
from   numpy            import linalg as LA
from   scipy            import linalg as SciLA
from   tools            import print_state,fidelity,dump_state,read_state,dump_lanz_vecs
from   hamiltonian      import Hmat,Hmoms

# ----- term-by-term ITE

def H_alpha_psi(H_,psi_,alpha_):
 (A,h,imp,gmp) = H_[alpha_] 
 phi = np.zeros(psi_.shape,dtype=complex)
 for m in np.where(np.abs(h)>1e-8)[0]:
  phi += h[m]*gmp[m,imp[m,:]]*psi_[imp[m,:]]
 return phi.copy()

def ExpmbH_alpha(H_,psi_,alpha_,db):
 phi = psi_.copy()
 chi = psi_.copy()
 i   = 0
 while(LA.norm(chi)>1e-8):
  chi  = (-db/float(i+1))*H_alpha_psi(H_,chi,alpha_)
  phi += chi.copy()
  i   += 1
 nu = LA.norm(phi)
 return phi.copy()/nu,nu

# ----- unitary evolution

def xP_psi(x_,psi_,imp_,gmp_):
 phi = np.zeros(psi_.shape,dtype=complex)
 for m in np.where(np.abs(x_)>1e-8)[0]:
  phi += x_[m]*gmp_[m,imp_[m,:]]*psi_[imp_[m,:]]
 return phi.copy()

def Exp_ixP(x_,psi_,imp_,gmp_):
 phi = psi_.copy()
 chi = psi_.copy()
 i   = 0
 while(LA.norm(chi)>1e-8):
  chi  = (1j/float(i+1))*xP_psi(x_,chi,imp_,gmp_)
  phi += chi.copy()
  i   += 1
 nu = LA.norm(phi)
 return phi.copy()/nu

def QITE_step(H_,psi_,db,xv,check):

 import time

 nalpha = len(H_)
 dn_    = 1.0

 if(xv is None): xv = [0]*nalpha

 for alpha in range(nalpha):
  # ----- target state

  #t0  = time.time()
  delta_alpha,dnalpha_ = ExpmbH_alpha(H_,psi_,alpha,db)
  delta_alpha         -= psi_.copy()
  dn_                 *= dnalpha_

  # ----- pauli action
  (A,h,imp,gmp) = H_[alpha]
  nact = imp.shape[0]
  Pmu_psi = np.zeros(imp.shape,dtype=complex)
  for m in range(nact):
   Pmu_psi[m,:] = gmp[m,imp[m,:]]*psi_[imp[m,:]]

  #t1  = time.time()
 
  # ----- set linear system 
  Amat = [[],[],[]]
  for mu in range(nact):
   Acol    = 2*np.real(np.einsum('a,na->n',np.conj(Pmu_psi[mu,:]),Pmu_psi))
   idx     = np.where(np.abs(Acol)>1e-8)[0]
   Amat[0] = Amat[0]+list(Acol[idx])
   Amat[1] = Amat[1]+[mu]*len(idx)
   Amat[2] = Amat[2]+list(idx)

  import scipy
  from scipy.sparse import csr_matrix
  from scipy.sparse import linalg as SPLA
  Amat_sparse = scipy.sparse.coo_matrix((Amat[0],(Amat[1],Amat[2])))
  Amat_T      = scipy.sparse.coo_matrix.transpose(Amat_sparse)

  #t2  = time.time()
  bvec = np.dot(Pmu_psi,np.conj(delta_alpha))
  bvec = -2.0*np.imag(bvec)
  #t3  = time.time()

  if(check):
   x = SPLA.lsqr(Amat_sparse,bvec)[0]
  else:
   zct = Amat_T.dot(bvec)
   def cost_fun(vct):
    kvc = Amat_sparse.dot(vct)
    return LA.norm(kvc-bvec)**2
   def J_cost_fun(vct):
    wct = Amat_sparse.dot(vct)
    wct = Amat_T.dot(wct)
    return 2.0*(wct-zct)
   import scipy
   x    = scipy.optimize.minimize(cost_fun,x0=xv[alpha],method='Newton-CG',jac=J_cost_fun).x
  xv[alpha] = x.copy()

  #t4  = time.time()
  psi_ = Exp_ixP(x,psi_,imp,gmp)
  #t5  = time.time()

 return psi_,dn_,xv

def Lanczos_QITE(hv,sv,db):
 nv = len(range(0,len(hv),2))
 hm = np.zeros((nv,nv),dtype=complex)
 sm = np.zeros((nv,nv),dtype=complex)
 for jr in range(0,len(hv),2):
  for js in range(0,len(hv),2):
   jk = (jr+js)//2
   sm[jr//2,js//2] = np.exp(2*sv[jk]-sv[jr]-sv[js])
   hm[jr//2,js//2] = hv[jk]*sm[jr//2,js//2]

 # rarefied sampling
 idx = []
 for l in range(nv):
  if(int(np.sqrt(2.0)**l)<nv+1): idx.append(l)
 if(nv-1 not in idx): idx.append(nv-1)
 sm = sm[idx,:]
 sm = sm[:,idx]
 hm = hm[idx,:]
 hm = hm[:,idx]
 nv = sm.shape[0]

 # regularization
 for jk in range(nv):
  sm[jk,jk] *= 1.0+2*db
  hm[jk,jk] *= 1.0+2*db

 eps,U = SciLA.eigh(hm,sm)
 eps   = np.real(eps)
 return np.min(eps)

def QITE(H_,db,bmax,lanczos=False,psi0=None,omega=None,ncheck=1):
 Hm    = Hmat(H_)
 N     = Hm.shape[0]
 nbit  = int(np.log2(N))
 eps,U = SciLA.eigh(Hm)
 m0    = np.argmin(eps)
 zeta  = np.exp(-db*(eps-eps[m0]))
 fide  = 1.0

 fout = open('QITE.out','w')
 fout.write("FCI gs energy %.6f \n" % eps[m0])
 fout.write("FCI gs wfn \n")
 print_state(U[:,m0],nbit,fout)

 psi_QITE = psi0[:]

 nbeta = int(bmax/db)+1
 hvect_LANZ = np.zeros(nbeta+1)
 svect_LANZ = np.zeros(nbeta+1)

 xv = None
 fout.write("QITE\n")
 for ib in range(nbeta):
  ea,ev           = Hmoms(H_,psi_QITE)
  hvect_LANZ[ib]  = ea

  if(omega is None): fide = fidelity(psi_QITE,U[:,m0])
  else:              fide = LA.norm(psi_QITE[omega])**2

  if(lanczos):
   ea_            = Lanczos_QITE(hvect_LANZ[:ib+1],svect_LANZ[:ib+1],db) 
   fout.write("%.6f %.6f %.6f %.6f %.6f \n" % (ib*db,ea,ev,fide,ea_))
  else:
   fout.write("%.6f %.6f %.6f %.6f \n" % (ib*db,ea,ev,fide))
  fout.flush()

  psi_QITE,dnorm,xv = QITE_step(H_,psi_QITE,db,xv,ib%ncheck==0)
  svect_LANZ[ib+1]  = svect_LANZ[ib]+np.log(dnorm)

 fout.write("QITE gs wfn \n")
 print_state(psi_QITE,nbit,fout)
 dump_state(psi_QITE,nbit,'qite.psi')
 dump_lanz_vecs(hvect_LANZ[:nbeta],svect_LANZ[:nbeta],'qlanz.vecs')

 fout.close()
