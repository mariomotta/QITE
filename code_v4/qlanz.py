import numpy as np
import sys
from   numpy import linalg as LA
from   scipy import linalg as SciLA

# ---------

def read_lanz_vecs(fname):
 V     = np.loadtxt(fname)
 nbeta = V.shape[0]
 hv    = V[:,1]
 sv    = V[:,2]
 return nbeta,hv,sv

def rarefied_sampling(hm_,sm_,idx):
 sm_ = sm_[idx,:]
 sm_ = sm_[:,idx]
 hm_ = hm_[idx,:]
 hm_ = hm_[:,idx]
 nv_ = sm_.shape[0]
 return nv_,hm_,sm_

def diagonal_regularization(hm_,sm_,f_=1):
 nv = hm_.shape[0]
 for jk in range(nv):
  hm_[jk,jk] *= f_
  sm_[jk,jk] *= f_
 return hm_,sm_

# ---------

def Lanczos_kernel(hm,sm,idx,Eref):
 nv,hm_,sm_ = rarefied_sampling(hm[:,:],sm[:,:],idx) 
 sigma,V    = SciLA.eigh(sm_)
 jdx        = np.where(sigma>1e-14)[0]

 sm_ = np.dot(V.T,np.dot(sm_,V))
 hm_ = np.dot(V.T,np.dot(hm_,V))
 sm_ = sm_[jdx,:]
 sm_ = sm_[:,jdx]
 hm_ = hm_[jdx,:]
 hm_ = hm_[:,jdx]

 eps,U = SciLA.eigh(hm_,sm_)
 return np.min(np.real(eps))

# ---------

def Lanczos_QITE(fname_,nmax_,db_,ds_,nnn=None):

 nbeta,hv,sv = read_lanz_vecs(fname_)
 Et = (hv[0]+hv[nbeta-1])/2
 for ib in range(nbeta):
  sv[ib] += db_*ib*Et
  hv[ib] -= Et

 vlst = range(0,nbeta,2)
 nv   = len(vlst)
 hm   = np.zeros((nv,nv))
 sm   = np.zeros((nv,nv))
 for jr in vlst:
  for js in vlst:
   jk = (jr+js)//2
   sm[jr//2,js//2] = np.exp(2*sv[jk]-sv[jr]-sv[js])
   hm[jr//2,js//2] = hv[jk]*sm[jr//2,js//2]

 nv = nmax_
 sm = sm[:nmax_,:nmax_]
 hm = hm[:nmax_,:nmax_]

 idx = [0]
 ii  = 0
 jj  = 0
 while(ii<nv and jj<nv-1):
  for jj in range(ii+1,nv):
   if(sm[ii,jj]<ds_):
    idx.append(jj)
    break
  ii=idx[len(idx)-1]

 idx += [nv-1]
 ea_ = 0
 ea_ = Lanczos_kernel(hm,sm,idx,ea_)
 print nmax_,ds_,ea_+Et,hm[nmax_-1,nmax_-1]+Et
 sys.stdout.flush()
