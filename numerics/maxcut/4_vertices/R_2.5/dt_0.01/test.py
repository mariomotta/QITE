import numpy as np
import sys

sys.path.append('../../../../../code_v4/')
from hamiltonian      import MaxCut,print_Hamiltonian,Hmat
from mf               import hom_mf_solution,hom_mf_state,hom_mf_energy,mf_solution,mf_state,mf_energy
from ite              import ITE_FCI
from qite             import QITE
from binary_functions import Bas2Int

nspin = 4
R     = 2.50
db    = 0.01
bmax  =  8.00

vrtex = [x for x in range(nspin)]
if(nspin==4): np.random.seed(1953)
else:         np.random.seed(14657)
links = []
for i in range(nspin):
 for j in range(i+1,nspin):
  x = np.random.randint(2,size=1)[0]
  if(x==1): links.append((i,j))

#if(nspin==4): links=[(0,1),(0,2),(0,3),(1,3),(2,3)] # IBM graph

gamma = (vrtex,links)

H  = MaxCut(gamma,R)
print_Hamiltonian(H)
Hm = np.diag(Hmat(H))
E0 = np.min(Hm)
omega = np.where(np.abs(Hm-E0)<1e-6)[0]

theta0    = np.random.random()*2.0*np.pi
theta1,e1 = hom_mf_solution(theta0,nspin,H)
psi_1     = hom_mf_state(theta1,nspin)
print "HOM mf-energy ",e1

ITE_FCI(H,db,bmax,psi0=psi_1,omega=omega)
QITE(H,db,bmax,lanczos=False,psi0=psi_1,omega=omega,ncheck=10)
