import numpy as np
import sys,shutil
sys.path.append('../../code_v4/')

from hamiltonian      import H_molecule,print_Hamiltonian
from ite              import ITE
from qite             import QITE
from binary_functions import Bas2Int

db   = 0.01
bmax = 8.00

# number of lines in file code_v4/h2.dat, from https://arxiv.org/pdf/1512.06860.pdf
for line in range(54):
 H   = H_molecule(line)
 print_Hamiltonian(H)
 # RHF initial guess
 nspin    = 2
 psi_0    = np.zeros(2**nspin,dtype=complex)
 psi_0[2] = 1.0
 #ITE(H,db,bmax,psi0=psi_0)
 QITE(H,db,bmax,lanczos=False,psi0=psi_0,ncheck=10)
 shutil.move('QITE.out','QITE.'+str(line))
