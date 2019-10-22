import numpy as np
import os,sys,subprocess

for line in range(54):
 V = np.loadtxt('/home/mmotta/QITE/imag_time_for_qcomp/code_v4/h2.dat').T
 R = V[0,line]
 subprocess.call(["sed -n 6,806p QITE."+str(line)+" > eofbeta.out"],shell=True)
 W = np.loadtxt('./eofbeta.out')
 for i in range(W.shape[0]//2):
  print R,W[i,0],W[i,1],0
 print " " 
