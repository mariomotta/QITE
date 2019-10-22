import numpy as np
import pickle
from binary_functions import Int2Bas

def fidelity(psi_,phi_):
 return np.abs(np.vdot(psi_,phi_))

def print_state(psi_,nbit,outf):
 for i in range(psi_.shape[0]):
  if(np.abs(psi_[i])>1e-4):
   for x in Int2Bas(i,2,nbit): outf.write(str(x))
   outf.write(" %.12f %.12f I \n" % (np.real(psi_[i]),np.imag(psi_[i])))

def dump_state(psi_,nbit,fname):
 outf = open(fname,'w')
 for i in np.where(np.abs(psi_)>1e-4)[0]:
  outf.write("%d %.12f %.12f \n" % (i,np.real(psi_[i]),np.imag(psi_[i])))
 outf.close()

def read_state(nbit,fname):
 psi = np.zeros(2**nbit,dtype=complex)
 V   = np.loadtxt(fname)
 idx = V[:,0].astype(int)
 psi[idx] = V[:,1]+1j*V[:,2]
 return psi

def dump_lanz_vecs(hv,sv,fname):
 nv = len(hv)
 outf = open(fname,'w')
 for i in range(nv):
  outf.write("%d %.12f %.12f \n" % (i,hv[i],sv[i]))
 outf.close() 

def dump_lanz_rte(hv,sv,fname):
 nv = len(hv)
 outf = open(fname,'w')
 for i in range(nv):
  outf.write("%d %.8f %.8f %.8f %.8f\n" % (i,np.real(hv[i]),np.imag(hv[i]),np.real(sv[i]),np.imag(sv[i])))
 outf.close()

# ------------------------------------------------- #

def dpbc(a,b,n):
 return np.min([(a-b)%n,(b-a)%n])

def dobc(a,b,n):
 return np.abs(a-b)

def dgr(graph,i,j):
 VV,EE = graph
 nbit  = len(VV)
 paths = []
 if(i==j): return 0

 for (a,b) in EE:
  if(i==a and j==b or i==b and j==a): return 1

 for (a,b) in EE:
  if(i==a):   paths.append([a,b])
  if(i==b):   paths.append([b,a])

 while(True):
  new_paths = []
  for p in paths:
   end = p[len(p)-1]
   for (a,b) in EE:
    if(end==a):
     if(b==j): return len(p)
     else:     new_paths.append(p+[b])
    if(end==b):
     if(a==j): return len(p)
     else:     new_paths.append(p+[a])
  paths = [x for x in new_paths]
