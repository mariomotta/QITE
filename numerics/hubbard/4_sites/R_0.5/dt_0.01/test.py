import numpy as np
import sys
sys.path.append('../../../../../code_v4/')

from hamiltonian      import Hubbard,print_Hamiltonian
from ite              import ITE
from qite             import QITE
from binary_functions import Bas2Int

norb =  4
R    =  0.50
db   =  0.05
bmax =   8.0
U    =   1.0

H = Hubbard(norb,R,U)
print_Hamiltonian(H)

# AFM initial guess

nspin = 2*norb
psi_0 = np.zeros(2**nspin,dtype=complex)
xvec  = [0]*nspin
for ii in range(0,nspin,4):
 xvec[ii]   = 1
 xvec[ii+3] = 1
xind        = Bas2Int(xvec,2)
psi_0[xind] = 1.0

ITE(H,db,bmax,psi0=psi_0)
QITE(H,db,bmax,lanczos=False,psi0=psi_0,ncheck=10)
