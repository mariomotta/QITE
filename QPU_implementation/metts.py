import numpy as np
import matplotlib.pyplot as plt
from pyquil.gates import *
from pyquil.noise import estimate_bitstring_probs
from pyquil import Program,get_qc
from helper import measure,propagate,estimate_bitstring_probs,estimate_assignment_probs,update_alist
from scipy.linalg import expm

def ansatz(p,qbits,ci0):
	# Ansatz
	if ci0 == 0:
		None
	elif ci0 == 1:
		p += X(qbits[0])

def measure_energy(alist,ci0,shots,qc,qbits,hm_list,correction_matrix):
	# Measure the energy at the end of each time step
	Energy = 0
	Nterms = len(hm_list)
	for i in range(len(hm_list)):
		for hm in hm_list[i]:
			for j in range(len(hm[0])):
				p = Program()
				p.wrap_in_numshots_loop(shots)
				ro = p.declare('ro','BIT',1)
				ansatz(p,qbits,ci0)
				propagate(p,alist,qbits)
				Energy += hm[1][j]*measure(p,ro,hm[0][j],qc,qbits,correction_matrix)
	return Energy

def get_expectation(alist,shots,ci0,qc,qbits,correction_matrix):
	# Obtain the expectation values of the Pauli string at each time step

	sigma_expectation = np.zeros([4],dtype=complex)
	for j in range(4):
		p = Program()
		p.wrap_in_numshots_loop(shots)
		ro = p.declare('ro','BIT',1)
		ansatz(p,qbits,ci0)
		propagate(p,alist,qbits)
		sigma_expectation[j] = measure(p,ro,j,qc,qbits,correction_matrix)
	return sigma_expectation

def qite_step(alist,shots,qc,qbits,correction_matrix,db,delta,ci0,hm_list):
	for j in range(len(hm_list)):
		sigma_expectation = get_expectation(alist,shots,ci0,qc,qbits,correction_matrix)
		norm = update_alist(sigma_expectation,alist,db,delta,hm_list[j])
	return alist

def metts(qc,qbits,shots,beta,nt,nmetts,hm_list,lblock=10):
	#nt = int(beta/(2*db))
	db = (beta/2)/nt
	finqc = open('test.dat',"w")
	ci0 = np.random.randint(0,2)
	Elst_qc = []
	for n in range(0,nmetts):
		alist = []
		for _ in range(nt):
			correction_matrix = estimate_assignment_probs(qbits[0],20000,qc)
			alist = qite_step(alist,shots,qc,qbits,correction_matrix,db,delta,ci0,hm_list)
		correction_matrix = estimate_assignment_probs(qbits[0],20000,qc)
		Elst_qc.append(measure_energy(alist,ci0,shots,qc,qbits,hm_list,correction_matrix))
		ci0 = collapse(ci0,alist,qc,qbits)
	Elst_qc = np.asarray(Elst_qc)
	err_qc = block(lblock,Elst_qc)
	E_qc = np.average(Elst_qc)
	print('QMETTs value ',E_qc,'+/- ',err_qc)
	return(Elst_qc)


def collapse(ci0,alist,qc,qbits):
	p = Program()
	ro = p.declare('ro','BIT',1)
	ansatz(p,qbits,ci0)
	propagate(p,alist,qbits)
	p += MEASURE(qbits[0],ro[0])
	qc = get_qc(qc)
	exe = qc.compile(p)
	res = qc.run(exe)
	if res[0][0] == 0:
		return 0
	elif res[0][0] == 1:
		return 1
	else:
		raise ValueError('Oppps')

def block(lblock,arr):
	l = len(arr)
	ndata = int(l/lblock)
	arr_ = (arr.reshape(ndata,lblock)).copy()
	eblock = np.average(arr_,axis=1)
	err = np.std(eblock)/float(np.sqrt(ndata))
	return err




if __name__ == '__main__':

	# ---- input parameters for METTs.
	# Simulation will loop through different beta and perform METTs simulation to estimate thermal averages.
	# Produces Figure 4(b) of https://arxiv.org/pdf/1901.07653.pdf
	shots = 1500
	beta_list = [1.0,2.0,3.0,4.0,5.0]
	nt = 10
	nmetts = 200
	qc = '1q-qvm'
	qbits = [0]
	hm_list = []
	hm_list.append([])
	hm_list[0].append([[1,3],[1/np.sqrt(2),1/np.sqrt(2)]])
	delta = 0.1
	
	for beta in beta_list:
		print('beta: ', beta)
		E = metts(qc,qbits,shots,beta,nt,nmetts,hm_list,lblock=10)

		Xgate = np.array([[0,1],[1,0]],dtype=complex)
		Ygate = np.array([[0,-1j],[1j,0]],dtype=complex)
		Zgate= np.array([[1,0],[0,-1]],dtype=complex)
		Hamiltonian = (Xgate+Zgate)/np.sqrt(2)
		ew,ev = np.linalg.eigh(Hamiltonian)
		Partition_fnc = np.sum(np.exp(-beta*ew))
		E = np.sum(ew*np.exp(-beta*ew))/Partition_fnc
		print('Exact ' ,E)

