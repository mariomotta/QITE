import numpy as np
import itertools
from   pyscf import lib
import sys
sys.path.append("./../code_v3/")
import lattice as lat
import pauli

def main(Nx, Ny, nmetts,beta,db):

 # set up system 
 lattice = [x for x in itertools.product(range(Nx),range(Ny))]
 nbit    = len(lattice)
 Nbasis  = 2**nbit
 nn      = [] # nearest neighbors
 dlst    = [] # list of distances between points
 d       = np.zeros((nbit,nbit))
 for i in range(nbit):
  for j in range(nbit):
   d[i,j]=np.sqrt(lat.pbc_dist(lattice[i][0],lattice[j][0],Nx)**2+lat.pbc_dist(lattice[i][1],lattice[j][1],Ny)**2)
   if(d[i,j] not in dlst): dlst.append(d[i,j])
  nn.append([j for j in range(nbit) if j!=i and d[i,j]<1+1e-6])
 dlst = sorted(dlst)
 
 # parameters for time evolution
 nt    = beta/(2*db)
 nexp  = (5,2) 
 Rstar = 1.2
 if (abs(nt-int(nt))>1e-9):
  raise ValueError("beta/2 cannot be divided by the time step!")
 else:
  nt   = int(nt)
 
 # construct Hamiltonian
 # 3 denotes the three spin components: x,y,z
 v = np.zeros(3*nbit) #one body
 w = np.zeros((3*nbit,3*nbit)) #two body
 for i in range(nbit):
  #v[3*i+2]=-1.0 # z direction
  v[3*i]=-1.0 # z direction transformed by Hadamard
 for i in range(nbit):
  for j in range(nbit):
   if(j in nn[i]):
    #w[3*i,3*j]=0.5 # x direction
    w[3*i+2,3*j+2]=0.5 # x direction transformed by Hadamard
 
 sxind,sxcof = pauli.pauli_action(nbit)
 ci0,hdiag,e,c = sub.kernel(nbit,Nbasis,v,w,sxind,sxcof)
 
 # get initial guess
 ci0   = np.zeros(Nbasis,dtype=np.complex128)
 ci0[np.random.randint(0,Nbasis)] = 1
 
 # loop for 1) imaginary time evolution 2) measurement 3) collapse
 Elst = []
 for k in range(nmetts):
  print "####Constructing METTS # %d ..."%(k+1)
  print "The chosen CPS is number: ",np.where(ci0!=0)[0][0]
  #print "\o/"*10
  for i in range(nt):
   ci0 = sub.expmbH_approximate(nbit,Nbasis,v,w,sxind,sxcof,db,nexp,d,Rstar,ci0,'it')
   Hci0 = sub.Hpsi(nbit,Nbasis,v,w,sxind,sxcof,ci0)
   E0   = np.real(np.dot(ci0.conj(),Hci0))
   #print "Time step: %d     E = %0.6f"%(i+1, E0)
  #print "o_O "*10
  Emetts = E0 # measure
  print "METTS # %d      E = %0.6f"%(k+1,E0/(Nx*Ny))
  Elst.append(Emetts)
  # collapse 
  ci0 = collapse_metts(ci0)

 Elst = np.asarray(Elst)/nbit
 err  = np.std(Elst)
 E    = np.average(Elst)
 return E, err
 

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def collapse_metts(ci0):
 # randomly pick a product state (from the basis set) based on the 
 # probabilitys
 ci0 = ci0/np.linalg.norm(ci0)
 pval = (ci0*ci0.conj()).real
 #print "Probability distribution for CPS basis:"
 #print pval
 #pval = pval/np.sum(pval)
 cps = np.random.multinomial(1,pval).astype(np.complex128)
 return cps  

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#def exact_ftsol(Nx, Ny, beta):
 '''
 Explicitly calculate the partition function and observable
 averages by constructing the Hamiltonian matrix and diagonalize it.
 '''
# # set up system 
# lattice = [x for x in itertools.product(range(Nx),range(Ny))]
# nbit    = len(lattice)
# Nbasis  = 2**nbit
# nn      = [] # nearest neighbors
# dlst    = [] # list of distances between points
# d       = np.zeros((nbit,nbit))
# for i in range(nbit):
#  for j in range(nbit):
#   d[i,j]=np.sqrt(lat.pbc_dist(lattice[i][0],lattice[j][0],Nx)**2+lat.pbc_dist(lattice[i][1],lattice[j][1],Ny)**2)
#   if(d[i,j] not in dlst): dlst.append(d[i,j])
#  nn.append([j for j in range(nbit) if j!=i and d[i,j]<1+1e-6])
# dlst = sorted(dlst)
#
# # construct Hamiltonian
# # 3 denotes the three spin components: x,y,z
# v = np.zeros(3*nbit) #one body
# w = np.zeros((3*nbit,3*nbit)) #two body
# for i in range(nbit):
#  v[3*i+2]=-1.0 # z direction
# for i in range(nbit):
#  for j in range(nbit):
#   if(j in nn[i]):
#    w[3*i,3*j]=0.5 # x direction
# 
# sxind,sxcof = sub.pauli_action(nbit)
#
# # Construct Hamiltonian matrix explicitly
# H = np.zeros((Nbasis,Nbasis),np.complex128)
# for i in range(Nbasis):
#  ci0 = np.zeros(Nbasis, dtype=np.complex128)
#  ci0[i] = 1.
#  H[i,:] = sub.Hpsi(nbit,Nbasis,v,w,sxind,sxcof,ci0)
# #H = (H + H.conj().T)/2. #symmetrize
# # solve Hamiltonian
# ew,ev = np.linalg.eigh(H)
# Z = np.sum(np.exp(-beta*ew))
# E = np.sum(ew*np.exp(-beta*ew))/Z
# return E/(Nx*Ny)
#
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#def test(Nx, Ny, beta,db):
 '''
 Test if the time evolution part is correct by doing
 the trace explicitly. (using the natural basis)
 '''
# # set up system 
# lattice = [x for x in itertools.product(range(Nx),range(Ny))]
# nbit    = len(lattice)
# Nbasis  = 2**nbit
# nn      = [] # nearest neighbors
# dlst    = [] # list of distances between points
# d       = np.zeros((nbit,nbit))
# for i in range(nbit):
#  for j in range(nbit):
#   d[i,j]=np.sqrt(lat.pbc_dist(lattice[i][0],lattice[j][0],Nx)**2+lat.pbc_dist(lattice[i][1],lattice[j][1],Ny)**2)
#   if(d[i,j] not in dlst): dlst.append(d[i,j])
#  nn.append([j for j in range(nbit) if j!=i and d[i,j]<1+1e-6])
# dlst = sorted(dlst)
# 
# # parameters for time evolution
# nt    = beta/(2*db)
# nexp  = (5,2) 
# Rstar = 1.2
# if (abs(nt-int(nt))>1e-9):
#  raise ValueError("beta/2 cannot be divided by the time step!")
# else:
#  nt   = int(nt)
# 
# # construct Hamiltonian
# # 3 denotes the three spin components: x,y,z
# v = np.zeros(3*nbit) #one body
# w = np.zeros((3*nbit,3*nbit)) #two body
# for i in range(nbit):
#  v[3*i+2]=-1.0 # z direction
# for i in range(nbit):
#  for j in range(nbit):
#   if(j in nn[i]):
#    w[3*i,3*j]=0.5 # x direction
# 
# sxind,sxcof = sub.pauli_action(nbit)
# ci0,hdiag,e,c = sub.kernel(nbit,Nbasis,v,w,sxind,sxcof)
# 
# E = 0.
# Z = 0.
# for i in range(Nbasis):
#  ci0 = np.zeros(Nbasis, dtype=ci0.dtype)
#  ci0[i] = 1.0
#  for t in range(2*nt):
#   ci0  = sub.expmbH_approximate(nbit,Nbasis,v,w,sxind,sxcof,db,nexp,d,Rstar,ci0,'it')
#   Hci0 = sub.Hpsi(nbit,Nbasis,v,w,sxind,sxcof,ci0)
#  E   += np.real(np.dot(ci0.conj(),Hci0))
#  Z   += np.real(np.dot(ci0,ci0.conj()))
# E = E/(Z*Nx*Ny)
# return E
   
  
 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

if __name__ == "__main__":

 Nx      = 3
 Ny      = 2
 nmetts  = 200
 beta    = 2.0 
 db      = 0.01 # time step 
 
 Efci = exact_ftsol(Nx, Ny, beta)
 E,err = main(Nx, Ny, nmetts,beta,db)
 #Etest = test(Nx, Ny, beta,db)
 print "fci result: ", Efci
 print "METTS result: %0.6f   standard deviation: %0.6f"%(E,err)
 #print "test time evolve: ", Etest
