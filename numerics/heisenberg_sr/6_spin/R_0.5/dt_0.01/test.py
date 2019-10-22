import numpy as np
import sys
sys.path.append('../../../../../code_v4/')

from hamiltonian      import Heisenberg_SR,print_Hamiltonian
from mf               import hom_mf_solution,hom_mf_state,hom_mf_energy,mf_solution,mf_state,mf_energy
from ite              import ITE_FCI
from qite             import QITE
from binary_functions import Bas2Int
nspin =  6
R     =  0.50
db    =  0.01
bmax  =  4.00

H = Heisenberg_SR(nspin,R)
print_Hamiltonian(H)

# AFM initial guess

psi_0       = np.zeros(2**nspin,dtype=complex)
xvec        = [0,1]*(nspin/2)
xind        = Bas2Int(xvec,2)
psi_0[xind] = 1.0

ITE_FCI(H,db,bmax,psi0=psi_0)
QITE(H,db,bmax,lanczos=False,psi0=psi_0,ncheck=10)
