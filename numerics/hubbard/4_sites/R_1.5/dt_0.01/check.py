import numpy
from pyscf import gto, scf, ao2mo, mcscf

'''
User-defined Hamiltonian for CASSCF module.

Defining Hamiltonian once for SCF object, the derivate post-HF method get the
Hamiltonian automatically.
'''

mol = gto.M()
mol.nelectron = 6
mol.incore_anyway = True

n  = 6
h1 = numpy.zeros((n,n))
for i in range(n-1):
 h1[i,i+1] = h1[i+1,i] = -1.0

for i in range(n):
 print h1[i,:]

eri = numpy.zeros((n,n,n,n))
for i in range(n):
 eri[i,i,i,i] = 10.0

mf            = scf.RHF(mol)
mf.get_hcore  = lambda *args: h1
mf.get_ovlp   = lambda *args: numpy.eye(n)
mf._eri       = ao2mo.restore(8, eri, n)
mf.init_guess = '1e'
mf.kernel()

mycas = mcscf.CASSCF(mf,ncas=n,nelecas=(3,3))
mycas.kernel()
