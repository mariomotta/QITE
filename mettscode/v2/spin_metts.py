import numpy as np
import itertools
from   pyscf import lib
from   scipy import linalg as LA

from   mooncake.code_v3.lattice          import torus
from   mooncake.code_v3.binary_functions import Int2Bas,Bas2Int,Opp2Str,Psi2Str
from   mooncake.code_v3.pauli            import sigma_matrices,pauli_action
from   mooncake.code_v3.ising            import Exp_mbH,H_average,H_ising_transverse,define_neighbors,MF_ave,MF_to_CI,HPsi
from   mooncake.code_v3.it2rt            import rt2it


def metts(Nx,Ny,nmetts,beta,db,realtime=False,eta=0.25*np.pi,outdir='./',tofile=True):

    # set up system 
    nbit,lattice,d = torus(Nx, Ny)
    Nbasis = 2**nbit
    # Sx->Sz and Sz->Sx so that |0>,|1> -> |+>,|->
    h_terms = H_ising_transverse(lattice,d,nbit,eta,4,hadamard=False)
    h_terms = define_neighbors(h_terms,lattice,d,2.2)
   
    # parameters for time evolution
    nt    = beta/(2*db)
    if (abs(nt-int(nt))>1e-9):
        raise ValueError("beta/2 cannot be divided by the time step!")
    else:
        nt   = int(nt)
   
    # get initial guess
    ci0   = np.zeros(Nbasis,dtype=np.complex128)
    ci0[np.random.randint(0,Nbasis)] = 1
    ind_sx,gmm_sx = pauli_action(nbit)
   
    # loop for 1) imaginary time evolution 2) measurement 3) collapse
    if tofile:
        fin  = open(outdir+"X%dY%d_b%0.1f.txt"%(Nx,Ny,beta), "w")
        fin.write("# metts#            E\n")
    Elst = []
    for k in range(nmetts):
        print "####Constructing METTS # %d ..."%(k+1)
        #print "The chosen CPS is number: ",np.where(ci0!=0)[0][0]
        # time evolution
        # collapse onto X 
        for i in range(nt):
            for (omega,act,hlist) in h_terms:
                if realtime:
                    ci0_it, nc = Exp_mbH(db, [(omega,act,hlist)], ind_sx, gmm_sx,ci0)
                    ci0        = rt2it(ci0,ci0_it,act,ind_sx,gmm_sx)
                else:
                    ci0, nc    = Exp_mbH(db, [(omega,act,hlist)], ind_sx, gmm_sx,ci0)
        # done time evolution
        ea, ev = H_average(h_terms,ci0,ind_sx,gmm_sx) 
        Emetts = ea # measure
        if tofile:
         fin.write("%d          %0.6f\n"%(k+1, Emetts))
        print "METTS # %d      E = %0.6f"%(k+1,Emetts/(Nx*Ny))
        Elst.append(Emetts)
        # collapse 
        if k%2 == 0:
            ci0 = collapse_metts(ci0, nbit, basis='X')
        else:
            ci0 = collapse_metts(ci0, nbit, basis='Z')
        # done collapsing onto X
           
            
    if tofile:
        fin.close()
   
    Elst = np.asarray(Elst)/nbit #energy per site
    err  = np.std(Elst)/np.sqrt(nmetts)
    E    = np.average(Elst)
    return E, err
    

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def collapse_metts(ci0, nbit, basis):
 # randomly pick a product state (from the basis set) based on the 
 # probabilitys
    if basis=='X':
        Hadamard = np.array([[1.,1.],[1.,-1.]])/np.sqrt(2)
        U = Hadamard.copy()
        for i in range(nbit-1):
            U = np.kron(U, Hadamard)
        ci0 = np.dot(U, ci0)

    ci0 = ci0/np.linalg.norm(ci0)
    pval = (ci0*ci0.conj()).real
        
 #print "Probability distribution for CPS basis:"
 #print pval
 #pval = pval/np.sum(pval)
    cps = np.random.multinomial(1,pval).astype(np.complex128)
    print "The chosen CPS is number: ",np.where(cps!=0)[0][0]
    if basis=='X':
        cps = np.dot(U,cps)
    return cps  

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def finiteT_ED(Nx, Ny, beta, eta=0.25*np.pi):
 '''
 Explicitly calculate the partition function and observable
 averages by constructing the Hamiltonian matrix and diagonalize it.
 '''
 # set up system 
 nbit,lattice,d = torus(Nx, Ny)
 Nbasis = 2**nbit
 # Sx->Sz and Sz->Sx so that |0>,|1> -> |+>,|->
 h_terms = H_ising_transverse(lattice,d,nbit,eta,4,hadamard=False)
 h_terms = define_neighbors(h_terms,lattice,d,2.2)
 ind_sx,gmm_sx = pauli_action(nbit)

 # Construct Hamiltonian matrix explicitly
 H = np.zeros((Nbasis,Nbasis),np.complex128)
 for i in range(Nbasis):
  ci0 = np.zeros(Nbasis, dtype=np.complex128)
  ci0[i] = 1.
  H[i,:] = HPsi(h_terms,ci0,ind_sx,gmm_sx)

 H = (H + H.conj().T)/2. #symmetrize
 # solve Hamiltonian
 ew,ev = np.linalg.eigh(H)
 Z = np.sum(np.exp(-beta*ew))
 E = np.sum(ew*np.exp(-beta*ew))/Z
 return E/(Nx*Ny)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
def full_trace(Nx, Ny, beta, db, realtime=False, eta=0.25*np.pi):
 '''
 Test if the time evolution part is correct by doing
 the trace explicitly. (using the natural basis)
 '''
 # set up system 
 nbit,lattice,d = torus(Nx, Ny)
 Nbasis = 2**nbit
 # Sx->Sz and Sz->Sx so that |0>,|1> -> |+>,|->
 h_terms = H_ising_transverse(lattice,d,nbit,eta,4,hadamard=False)
 h_terms = define_neighbors(h_terms,lattice,d,2.2)
 ind_sx,gmm_sx = pauli_action(nbit)

 # parameters for time evolution
 nt    = beta/(2*db)
 if (abs(nt-int(nt))>1e-9):
  raise ValueError("beta/2 cannot be divided by the time step!")
 else:
  nt   = int(nt)

 E = 0.
 Z = 0.
 for i in range(Nbasis):
  ci0 = np.zeros(Nbasis, dtype=np.complex128)
  ci0[i] = 1.0
  prefact = 1.0
  for t in range(nt):
   for (omega,act,hlist) in h_terms:
    if realtime:
     ci0_it, nc = Exp_mbH(db, [(omega,act,hlist)], ind_sx, gmm_sx,ci0)
     ci0    = rt2it(ci0,ci0_it,act,ind_sx,gmm_sx)
    else:
     ci0, nc  = Exp_mbH(db, [(omega,act,hlist)], ind_sx, gmm_sx,ci0)
    prefact *= nc
  #print prefact

  ci0  /= np.linalg.norm(ci0) 
  ci0 *= prefact
  ea,ev = H_average(h_terms,ci0,ind_sx,gmm_sx)
  E    += ea
  Z    += np.real(np.dot(ci0,ci0.conj()))

 E = E/(Z*Nx*Ny)
 return E
   
  
 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

if __name__ == "__main__":

 Nx      = 2
 Ny      = 1
 nmetts  = 50
 beta    = 2.0
 db      = 0.1 # time step 
 eta     = 0.25*np.pi
 # FCI solution 
 Efci = finiteT_ED(Nx, Ny, beta)
 # Test real time evolution
 #Etest = full_trace(Nx, Ny, beta,db)
 # Metts
 E,err = metts(Nx, Ny, nmetts,beta,db,realtime=False,eta=eta)
 print "fci result: ", Efci
 #print "test time evolve: ", Etest
 print "METTS result: %0.6f   standard deviation: %0.6f"%(E,err)
