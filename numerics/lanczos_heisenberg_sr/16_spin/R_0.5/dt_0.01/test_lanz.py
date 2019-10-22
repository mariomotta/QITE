import numpy as np
import sys
sys.path.append('../../../../../code_v4/')
from   qlanz import Lanczos_QITE

db   = 0.01
k    = 1
klst = list(np.arange(k,200,k))
if(k>1): klst = [1]+klst
nnn  = 9

for ds in [0.95,1.00]:
 for n in klst:
  Lanczos_QITE('qlanz.vecs',n,db,ds,nnn)
 print " "
 print " "
