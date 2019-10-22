The source (code_v4) is organized in several python files

binary_functions.py functions to convert a string of n integers modulo 2 into an integer between 0 and 2**n-1 and viceversa, eg 101 = 5 
pauli.py            functions to deal with the algebra of Pauli operators, written as integer strings modulo 4, eg 3102 = ZXIY, and to apply a Pauli operator to an element of the computational basis
hamiltonian.py      defines the Hamiltonians we studied in the QITE paper, the Heisenberg short-range, Maxcut etc, as linear combinations of Pauli operators; 
                    contains functions to apply H to a wavefunction, construct the matrix associated to the Hamiltonian in the computational basis, 
                    and compute the moments <x|H**k|x> of H on a wavefunction x
mf.py               produces the lowest-energy "symmetry-adapted" ie |psi,psi,psi,...> or "broken-symmetry", i.e. |psi1,psi2,psi3...> factorized state, 
                    to serve as starting point for the imaginary-time evolution
ite.py	            performs the imaginary-time evolution (not the QITE) in various ways:
                    it either (1) constructs the matrix associate to H in the computational basis, diagonalizes it, and applies its exponential to a wavefunction x
                    or (2) applies exp(-dt*H) to x by Taylor series
                    it also prepares the matrix elements of the Hamiltonian and overlap matrices needed by the QLanczos algorithm
qite.py             performs the quantum imaginary time evolution -- constructs the matrix A and the vector b (Amat,bvec) of the linear system, 
                    solves it by least-squares methods, and applies the corresponding unitary transformation to the current wavefunction
                    it also prepares the approximate matrix elements of the Hamiltonian and overlap matrices actually used by the QLanczos algorithm
qlanz.py            performs the quantum Lanczos calculations, and implements their numeric stabilization
tools.py            readout and computation of distances in OBC, PBC
style.py            conventions for the figures

