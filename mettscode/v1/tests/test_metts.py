import numpy as np
import sys
sys.path.append("../")
sys.path.append("/Users/sunchong/WorkStation/Mario/mooncake/code_v1")
from spin_metts import *

Nx      = 3
Ny      = 2
nmetts  = 200
db      = 0.01 # time step 
blist   = [0.5, 1.0, 2.0, 3.0, 4.0]
result = open("energy_%dx%d.txt"%(Nx, Ny),"w")
result.write("# Nx = %d ; Ny = %d \n"%(Nx, Ny))
result.write("# Nmetts = %d\n"%nmetts)
result.write("# time step = %f\n"%db)
result.write("#==========================================\n")
result.write("#beta      Efci         Emetts       err\n")
for beta in blist:
 Efci = exact_ftsol(Nx, Ny, beta)
 E,err = main(Nx, Ny, nmetts,beta,db)
 result.write(" %.2f    %.6f    %.6f    %.6f\n"%(beta,Efci,E,err))
result.close()
