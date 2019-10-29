import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as SciLA
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
		c = update_alist(sigma_expectation,alist,db,delta,hm_list[j])
	return alist,c

def qlanz(i,E,norm,ds=0.7,dlanz=1e-2):
	vect = list(np.arange(1,i+2,2))
	n = len(vect)
	H = np.zeros([n,n],dtype=complex)
	S = np.zeros([n,n],dtype=complex)
	j = 0
	k = 0
	for l in vect:
		k = 0
		for lp in vect:
			S[j,k] = norm[(l+lp)//2]**2/(norm[l]*norm[lp])
			H[j,k] = E[(l+lp)//2]*S[j,k]
			k += 1
		j += 1
	# Determine which vectors to keep to stabilize lanczos scheme
	idx = [0]
	ii = 0
	jj = 0
	while(ii<n and jj<n-1):
		for jj in range(ii+1,n):
			if(S[ii,jj]<ds):
				idx.append(jj)
				break
		ii = idx[len(idx)-1]
	if S.shape[0]-1 not in idx:
		idx += [S.shape[0]-1]
	H_ = H[idx,:]
	H_ = H_[:,idx]
	S_ = S[idx,:]
	S_ = S_[:,idx]
	n = S_.shape[0]

	# Regularize the overlap matrix
	sigma,V = np.linalg.eigh(S_)
	jdx = np.where(sigma>dlanz)[0]
	S_ = np.dot(V.T,np.dot(S_,V))
	H_ = np.dot(V.T,np.dot(H_,V))
	S_ = S_[jdx,:]
	S_ = S_[:,jdx]
	H_ = H_[jdx,:]
	H_ = H_[:,jdx]

	eps,c = SciLA.eig(H_,S_)
	return np.min(eps)


def qite(qc,qbits,shots,db,delta,N,hm_list):
	E = np.zeros([N+1],dtype=complex)
	norm = np.zeros([N+1],dtype=complex)
	EQlanz  = np.zeros([(N+1)//2],dtype=complex)
	alist = []

	# Readout error mitigation, estimate p(0|0),p(0|1),p(1|0),p(1|1)
	correction_matrix = estimate_assignment_probs(qbits[0],20000,qc)
	E[0] = measure_energy(alist,shots,qc,qbits,hm_list,correction_matrix)
	norm[0] = 1.
	# Qite main loop
	flag = 0
	for i in range(1,N+1):
		correction_matrix = estimate_assignment_probs(qbits[0],20000,qc)
		alist,norm_c = qite_step(alist,shots,qc,qbits,correction_matrix,db,delta,hm_list)
		norm[i] = norm[i-1]*(norm_c)
		E[i] = measure_energy(alist,shots,qc,qbits,hm_list,correction_matrix)
		if np.mod(i,2):
			EQlanz[flag] = qlanz(i,E,norm)
			flag += 1

	return E,EQlanz

if __name__ == '__main__':

	# ---- input parameters for qite
	# Produces Figure 2(e) of https://arxiv.org/pdf/1901.07653.pdf
	
	N = 50
	shots = 20000
	db = 0.1
	qc = '1q-qvm'
	qbits = [0]
	hm_list = []
	hm_list.append([])
	hm_list[0].append([[1,3],[1/np.sqrt(2),1/np.sqrt(2)]])
	delta = 0.1

	E,EQlanz = qite(qc,qbits,shots,db,delta,N,hm_list)
	plt.plot(np.arange(0,N+1)*db,E,'ro',label='QITE')
	plt.plot(np.arange(1,N+1,2)*db,EQlanz,'-ko',label='QLanz')
	plt.legend()
	plt.grid()
	plt.show()

