import numpy
from pyscf import gto, scf, ao2mo, mcscf

'''
User-defined Hamiltonian for CASSCF module.

Defining Hamiltonian once for SCF object, the derivate post-HF method get the
Hamiltonian automatically.
'''

n   = 4
mol = gto.M()
mol.nelectron = n
mol.incore_anyway = True
U  = 1.0
h1 = numpy.zeros((n,n))
for i in range(n-1):
 h1[i,i+1] = h1[i+1,i] = -1.0

for i in range(n):
 h1[i,i]   = -U/2

for i in range(n):
 print h1[i,:]

eri = numpy.zeros((n,n,n,n))
for i in range(n):
 eri[i,i,i,i] = U

mf            = scf.RHF(mol)
mf.get_hcore  = lambda *args: h1
mf.get_ovlp   = lambda *args: numpy.eye(n)
mf._eri       = ao2mo.restore(8, eri, n)
mf.init_guess = '1e'
mf.kernel()

mycas = mcscf.CASSCF(mf,ncas=n,nelecas=(n//2,n//2))
mycas.kernel()
