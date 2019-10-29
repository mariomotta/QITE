import numpy as np
import matplotlib.pyplot as plt
from pyquil.gates import *
from pyquil.noise import estimate_bitstring_probs
from pyquil import Program,get_qc
from helper import measure,propagate,estimate_bitstring_probs,estimate_assignment_probs,update_alist

def ansatz(p,qbits):
	None

def measure_energy(alist,shots,qc,qbits,hm_list,correction_matrix):
	# Measure the energy at the end of each time step
	Energy = 0
	Nterms = len(hm_list)
	for i in range(len(hm_list)):
		for hm in hm_list[i]:
			for j in range(len(hm[0])):
				p = Program()
				p.wrap_in_numshots_loop(shots)
				ro = p.declare('ro','BIT',1)
				ansatz(p,qbits)
				propagate(p,alist,qbits)
				Energy += hm[1][j]*measure(p,ro,hm[0][j],qc,qbits,correction_matrix)
	return Energy

def get_expectation(alist,shots,qc,qbits,correction_matrix):
	# Obtain the expectation values of the Pauli string at each time step

	sigma_expectation = np.zeros([4],dtype=complex)
	for j in range(4):
		p = Program()
		p.wrap_in_numshots_loop(shots)
		ro = p.declare('ro','BIT',1)
		ansatz(p,qbits)
		propagate(p,alist,qbits)
		sigma_expectation[j] = measure(p,ro,j,qc,qbits,correction_matrix)
	return sigma_expectation

def qite_step(alist,shots,qc,qbits,correction_matrix,db,delta,hm_list):
	for j in range(len(hm_list)):
		sigma_expectation = get_expectation(alist,shots,qc,qbits,correction_matrix)
		norm = update_alist(sigma_expectation,alist,db,delta,hm_list[j])
	return alist

def qite(qc,qbits,shots,db,delta,N,hm_list):
	E = np.zeros([N+1],dtype=complex)
	alist = []

	# Readout error mitigation, estimate p(0|0),p(0|1),p(1|0),p(1|1)
	correction_matrix = estimate_assignment_probs(qbits[0],20000,qc)
	E[0] = measure_energy(alist,shots,qc,qbits,hm_list,correction_matrix)

	# Qite main loop
	for i in range(1,N+1):
		correction_matrix = estimate_assignment_probs(qbits[0],20000,qc)
		alist = qite_step(alist,shots,qc,qbits,correction_matrix,db,delta,hm_list)
		E[i] = measure_energy(alist,shots,qc,qbits,hm_list,correction_matrix)
	return E

if __name__ == '__main__':

	# ---- input parameters for qite
	# Produces Figure 2(e) of https://arxiv.org/pdf/1901.07653.pdf
	N = 25
	shots = 1000
	db = 0.1
	qc = '1q-qvm'
	qbits = [0]
	hm_list = []
	hm_list.append([])
	hm_list[0].append([[1],[1/np.sqrt(2)]])
	hm_list.append([])
	hm_list[1].append([[3],[1/np.sqrt(2)]])
	delta = 0.1
	
	E = qite(qc,qbits,shots,db,delta,N,hm_list)
	plt.plot(np.arange(0,N+1)*db,E,'ro',label='QITE')
	plt.grid()
	plt.show()

