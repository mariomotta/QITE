import numpy as np
import itertools
from   scipy import optimize as opt
from   pauli            import pauli_action
from   binary_functions import Bas2Int,Int2Bas,Opp2Str
from   numpy            import linalg as LA
from   scipy            import linalg as SciLA
from   tools            import print_state,fidelity,dgr,dpbc,dobc
from   pauli            import sigma_matrices

# ------------------------------------------------- #

def Hpsi(H_,psi_):
 phi = np.zeros(psi_.shape,dtype=complex)
 for (A,h,imp,gmp) in H_:
  for m in np.where(np.abs(h)>1e-8)[0]:
   phi += h[m]*gmp[m,imp[m,:]]*psi_[imp[m,:]]
 return phi.copy()

def Hmat(H_):
 N  = H_[0][2].shape[1]
 Hm = np.zeros((N,N),dtype=complex)
 for i in range(N):
  ei      = np.zeros(N,dtype=complex)
  ei[i]   = 1.0
  Hm[:,i] = Hpsi(H_,ei).copy()
 return Hm

def Hmoms(H_,psi_):
 phi_ = Hpsi(H_,psi_)
 ea   = np.vdot(psi_,phi_)
 ev   = np.vdot(phi_,phi_)
 return np.real(ea),np.real(ev-ea**2)

def print_Hamiltonian(H_):
 mu = 0
 for (A,h,imp,gmp) in H_:
  nact = len(A)
  print "term ",mu
  print "active qubits ",A
  print "operators: "
  for m in np.where(np.abs(h)>1e-8)[0]:
   print Opp2Str(Int2Bas(m,4,nact)),h[m]
  mu += 1

def Hii(H_,i):
 N     = H_[0][2].shape[1]
 nbit  = int(np.log2(N))
 hii = 0.0
 xi  = Int2Bas(i,2,nbit)
 for (A,h,imp,gmp) in H_:
  nact = len(A)
  for m in np.where(np.abs(h)>1e-8)[0]:
   sm   = Int2Bas(m,4,nact)
   smx  = [ sigma_matrices[xi[A[w]],xi[A[w]],sm[w]] for w in range(nact)]
   hii += np.real(h[m]*np.prod(smx))
 return hii

# ------------------------------------------------- #

def Heisenberg_SR(nspin,R):
 H = []
 imax = nspin
 if(nspin==2): imax = nspin-1
 for i in range(imax):
  j=(i+1)%nspin
  active = [k for k in range(nspin) if dpbc(i,k,nspin)<R or dpbc(j,k,nspin)<R]
  active = np.asarray(active)
  nact   = len(active)
  # -----
  h_alpha = np.zeros(4**nact)
  ii = np.where(active==i)[0][0]
  jj = np.where(active==j)[0][0]
  for alpha in range(1,4):
   idx     = [0]*nact
   idx[ii] = alpha
   idx[jj] = alpha
   h_alpha[Bas2Int(idx,4)] = 1.0
  # -----
  imap,gmap = pauli_action(active,nspin)
  H.append((active,h_alpha,imap,gmap))
 return H

def Heisenberg_LR(nspin,R):
 H = []
 for i in range(nspin):
  for j in range(i+1,nspin):
   print i,j
   # -----
   active = [k for k in range(nspin) if dobc(i,k,nspin)<R or dobc(j,k,nspin)<R]
   active = np.asarray(active)
   nact   = len(active)      
   print active
   # -----
   h_alpha = np.zeros(4**nact)
   ii = np.where(active==i)[0][0]
   jj = np.where(active==j)[0][0]
   for alpha in range(1,4):
    idx     = [0]*nact
    idx[ii] = alpha
    idx[jj] = alpha
    h_alpha[Bas2Int(idx,4)] = 1.0/(dobc(i,j,nspin)+1.0)
   # -----
   #print h_alpha
   imap,gmap = pauli_action(active,nspin)
   H.append((active,h_alpha,imap,gmap))
 return H

def Ising(nspin,R,psi):
 H = []
 for i in range(nspin):
  j = (i+1)%nspin
  # -----
  active = [k for k in range(nspin) if dpbc(i,k,nspin)<R or dpbc(j,k,nspin)<R]
  active = np.asarray(active)
  nact   = len(active)
  # -----
  h_alpha = np.zeros(4**nact)
  ii = np.where(active==i)[0][0]
  jj = np.where(active==j)[0][0]

  idx     = [0]*nact
  idx[ii] = 3
  h_alpha[Bas2Int(idx,4)] = np.sin(psi)/2.0
  idx     = [0]*nact
  idx[jj] = 3
  h_alpha[Bas2Int(idx,4)] = np.sin(psi)/2.0
  idx     = [0]*nact
  idx[ii] = 1
  idx[jj] = 1
  h_alpha[Bas2Int(idx,4)] = np.cos(psi)

  # -----
  imap,gmap = pauli_action(active,nspin)
  H.append((active,h_alpha,imap,gmap))
 return H

def MaxCut(graph,R):
 VV,EE = graph
 nbit  = len(VV)
 H     = []
 for (i,j) in EE:
  # -----
  active = [k for k in range(nbit) if dgr(graph,i,k)<R or dgr(graph,j,k)<R]
  active = np.asarray(active)
  nact   = len(active)
  # -----
  h_alpha = np.zeros(4**nact)
  ii = np.where(active==i)[0][0]
  jj = np.where(active==j)[0][0]
  
  idx = [0]*nact
  h_alpha[Bas2Int(idx,4)] = -0.5
  idx[ii] = 3
  idx[jj] = 3
  h_alpha[Bas2Int(idx,4)] =  0.5

  # -----
  imap,gmap = pauli_action(active,nbit)
  H.append((active,h_alpha,imap,gmap))
 return H

def Hubbard(norb,R,U):
 H     = []
 nspin = 2*norb

 # ----- encoding:
 # 0u 0d 1u 1d 2u 2d 3u 3d ... (n-1)u (n-1)d
 #  0  1  2  3  4  5  6  7 ...   2n-2   2n-1

 # ----- formulas:
 # n0 => (1-Zp)/2
 # (a*_p a_q + a*_q a_p) =  (1/2) X_p X_q (prod_{k=q+1}^{p-1} Z_k ) (1- Z_p Z_q)   with p>q

 # ----- potential energy:
 # n_iu n_id = (1-Z_{2i}) (1-Z_{2i+1})/4      for i = 0 ... n-1

 # ----- kinetic energy (open boundary conditions):
 # a*_{iu} a_{i+1u} + hc => (a*_p a_q + a*_q a_p)                                  with p=2i   q=2i+2, i=0 ... n-2
 # a*_{id} a_{i+1d} + hc => (a*_p a_q + a*_q a_p)                                  with p=2i+1 q=2i+3, i=1 ... n-2

 # ----- potential energy
 for i in range(norb-1):
  print ">>>>>>> sites ",i,i+1
  print "neighborhood"
  dij   = np.asarray([min(np.abs(i-j),np.abs(i+1-j)) for j in range(norb)])
  idx   = np.where(dij<R)[0]
  print idx
  pmin  = 2*min(idx)
  pmax  = 2*max(idx)+1
  act   = range(pmin,pmax+1)
  nact  = len(act)
  print act

  print "-----"
  h_alpha = np.zeros(4**nact)
  for k in (i,i+1):  
   pk = 2*k-pmin
   print "interaction on site ",k
   wk = 0.5

   #idx     = [0]*nact
   #h_alpha[Bas2Int(idx,4)] += U*wk/4.0+(U/2)*wk/2.0

   if(k==0 or k==norb-1): wk = 1.0
 
   print "indices",pk,pk+1,nact,wk
   # -----
   idx     = [0]*nact
   idx[pk] = 3
   print idx
   h_alpha[Bas2Int(idx,4)] = -U*wk/4.0+(U/2)*wk/2.0 # half-filling Hubbard
   # -----
   idx       = [0]*nact
   idx[pk+1] = 3
   print idx
   h_alpha[Bas2Int(idx,4)] = -U*wk/4.0+(U/2)*wk/2.0 # half-filling Hubbard
   # -----
   idx       = [0]*nact
   idx[pk]   = 3
   idx[pk+1] = 3
   print idx
   h_alpha[Bas2Int(idx,4)] =  U*wk/4.0
   # -----

  print "kinetic"
  for sigma in (0,1):
   print "spin ",sigma
   p = (2*i+sigma)-pmin
   q = p+2
   # -----
   idx      = [0]*nact
   idx[p]   = 1
   idx[p+1] = 3
   idx[q]   = 1
   print idx
   h_alpha[Bas2Int(idx,4)] = -0.5
   # -----
   idx      = [0]*nact
   idx[p]   = 2
   idx[p+1] = 3
   idx[q]   = 2
   print idx
   h_alpha[Bas2Int(idx,4)] = -0.5

  imap,gmap = pauli_action(act,nspin)
  H.append((act,h_alpha,imap,gmap))

  print "remember to add the constant correction ",U*norb/4.0-(U/2)*2*norb/2.0

 return H

def H_molecule(line):
 # from O'Malley et al, https://arxiv.org/pdf/1512.06860.pdf
 nspin = 2
 V = np.loadtxt('../../code_v4/h2.dat').T
 V = V[:,line]
 H = []
 active = [0,1]
 active = np.asarray(active)
 nact   = len(active)
 # -----
 h_alpha = np.zeros(4**nact)
 h_alpha[Bas2Int([0,0],4)] = V[1]
 h_alpha[Bas2Int([3,0],4)] = V[2]
 h_alpha[Bas2Int([0,3],4)] = V[3]
 h_alpha[Bas2Int([3,3],4)] = V[4]
 h_alpha[Bas2Int([1,1],4)] = V[5]
 h_alpha[Bas2Int([2,2],4)] = V[6]
 # -----
 imap,gmap = pauli_action(active,nspin)
 H.append((active,h_alpha,imap,gmap))
 return H

