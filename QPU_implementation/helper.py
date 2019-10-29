import numpy as np
from pyquil.gates import *
from pyquil.noise import estimate_bitstring_probs,correct_bitstring_probs
from pyquil import Program,get_qc
from index import idx,coeff


def measure(p,ro,idx,qc,qbits,correction_matrix):
	# Circuit to measure the expectation value of any Pauli string

	if idx == 0:
		return 1
	elif idx == 1:
		qc = get_qc(qc)
		p += H(qbits[0])
		p += MEASURE(qbits[0],ro[0])
		exe = qc.compile(p)
		res = qc.run(exe)
		probs = estimate_bitstring_probs(res)
		probs = correct_bitstring_probs(probs,[correction_matrix])
		return probs[0] - probs[1]
	elif idx == 2:
		qc = get_qc(qc)
		p += RX(np.pi/2,qbits[0])
		p += MEASURE(qbits[0],ro[0])
		exe = qc.compile(p)
		res = qc.run(exe)
		probs = estimate_bitstring_probs(res)
		probs = correct_bitstring_probs(probs,[correction_matrix])
		return probs[0] - probs[1]
	elif idx == 3:
		qc = get_qc(qc)
		p += MEASURE(qbits[0],ro[0])
		exe = qc.compile(p)
		res = qc.run(exe)
		probs = estimate_bitstring_probs(res)
		probs = correct_bitstring_probs(probs,[correction_matrix])
		return probs[0] - probs[1]


def propagate(p,alist,qbits):
	# Circuit to propagate the state
	if len(alist) == 0:
		None
	else:
		for t in range(len(alist)):
			for gate in range(1,4):
				angle = np.real(alist[t][gate])
				if gate == 1:
					p += RX(angle,qbits[0])
				elif gate == 2:
					p += RY(angle,qbits[0])
				elif gate == 3:
					p += RZ(angle,qbits[0])
				else:
					raise ValueError

def update_alist(sigma_expectation,alist,db,delta,hm):
	# Obtain A[m]

	# Step 1: Obtain S matrix
	S = np.zeros([4,4],dtype=complex)
	for i in range(4):
		for j in range(4):
			S[i,j] = sigma_expectation[idx[i,j]]*coeff[i,j]

	# Step 2: Obtain b vector
	b = np.zeros([4],dtype=complex)
	c = 1
	for i in range(len(hm[0][0])):
		c -= 2*db*hm[0][1][i]*sigma_expectation[hm[0][0][i]]
	c = np.sqrt(c)
	for i in range(4):
		b[i] += (sigma_expectation[i]/c-sigma_expectation[i])/(db)
		for j in range(len(hm[0][0])):
			b[i] -= hm[0][1][j]*coeff[i,hm[0][0][j]]*sigma_expectation[idx[i,hm[0][0][j]]]/c 
		b[i] = 1j*b[i] - 1j*np.conj(b[i])

	# Step 3: Add regularizer
	dalpha = np.eye(4)*delta

	# Step 4: Solve for linear equation, the solution is multiplied by -2 because of the definition of unitary rotation gates is exp(-i theta/2)
	x = np.linalg.lstsq(S+np.transpose(S)+dalpha,-b,rcond=-1)[0]
	alist.append([])
	for i in range(len(x)):
		alist[-1].append(-x[i]*2*db)
	return c

def estimate_assignment_probs(q, shots, qc_name,p00=None,p11=None):
    """
    Estimate the readout assignment probabilities for a given qubit ``q``.
    The returned matrix is of the form::

            [[p00 p01]
             [p10 p11]]

    :param int q: The index of the qubit.
    :param int trials: The number of samples for each state preparation.
    :param qc: The quantum computer we wish to sample.
    :return: The assignment probability matrix
    :rtype: np.array
    """
    qc = get_qc(qc_name)
    if p00 is None:
        p = Program()
        p.wrap_in_numshots_loop(shots)
        ro = p.declare('ro','BIT',1)
        p += I(q)
        p += MEASURE(q,ro[0])
        exe = qc.compile(p)
        results_i = np.sum(qc.run(exe))

        p = Program()
        p.wrap_in_numshots_loop(shots)
        ro = p.declare('ro','BIT',1)
        p += X(q)
        p += MEASURE(q,ro[0])
        exe = qc.compile(p)
        results_x = np.sum(qc.run(exe))
    else:
        p = Program()
        p.define_noisy_readout(q,p00=p00,p11=p11)
        ro = p.declare('ro','BIT',1)
        p += I(q)
        p += MEASURE(q,ro[0])
        cxn = QVMConnection()
        results_i = np.sum(cxn.run(p,[0],trials=shots))

        p = Program()
        p.define_noisy_readout(q,p00=p00,p11=p11)
        ro = p.declare('ro','BIT',1)
        p += X(q)
        p += MEASURE(q,ro[0])
        cxn = QVMConnection()
        results_x = np.sum(cxn.run(p,[0],trials=shots))


    p00 = 1. - results_i / float(shots)
    p11 = results_x / float(shots)
    return np.array([[p00, 1 - p11],
                     [1 - p00, p11]])