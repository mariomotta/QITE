import numpy as np
import os,subprocess
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/mmotta/QITE/imag_time_for_qcomp/code_v4/')
import style as st

def qite_plot(RR,TT,EE,FF,FCI,xrng,yrng,fname):
 x0        = 7.00
 xfig      = 1.00*1*x0
 yfig      = 0.75*1*x0

 fig,pan = plt.subplots(1,1,figsize=(xfig,yfig))

 pan.get_xaxis().set_tick_params(which='both',direction='in',labelsize=20)
 pan.get_yaxis().set_tick_params(which='both',direction='in',labelsize=20) 

 pan.set_xlabel(r'$\displaystyle{\beta}$',fontsize=20)
 pan.set_ylabel(r'$\displaystyle{E(\beta)}$',fontsize=20)
 pan.set_xlim(xrng[0])
 pan.set_ylim(yrng[0])

 ip=0
 u = np.arange(round(np.min(xrng[ip])),round(np.max(xrng[ip])),0.05)
 v = [FCI[ip]]*len(u)
 st.Palette[st.fci].dataplot(pan,u,v,spline=True,connected=False)

 for jR in range(len(Rlst)):
  Lx = np.max(xrng[0])-np.min(xrng[0])
  Ly = np.max(yrng[0])-np.min(yrng[0])
  mask = st.equidistant_points(TT[:,jR],Lx,EE[:,jR],Ly*3.0,0.015)
  st.Palette[jR].dataplot(pan,TT[mask,jR],EE[mask,jR]    ,spline=True ,connected=False)
  st.Palette[jR].dataplot(pan,TT[mask,jR],EE[mask,jR]    ,spline=False,connected=False)
  st.Palette[jR].dataplot(pan,TT[mask,jR],EE[mask,jR]-5e3,spline=False,connected=True)

 plt.text(3.55,-2.05,r'$(a)$',fontsize=21)

 x0 = -0.04
 dx =  1.08
 y0 =  1.00
 dy =  0.10
 order = [0]+[3*(x+1) for x in range(len(Rlst))]
 lgd   = st.make_legend(pan,order,ncol=range(len(Rlst)),bbox_to_anchor=(x0,y0,dx,dy),\
         loc=3,mode="expand",borderaxespad=1,fontsize=14,handlelength=3,numpoints=2)

 fig.savefig(fname,format='eps',bbox_inches='tight')
 os.system('ps2epsi '+fname)
 os.system('mv '+fname+'i '+fname)

# ------------------------------------------------- #

Nlst  = [2,4,6] #,8]
strt  = [6,10,24,74,256]
dtlst = [0.01]
fnsh  = [400,800]
Efci  = [-1.500000,-3.333333,-4.664214,-6.266279]

xrng  = [[[ -0.2, 4.2],[-0.20,4.20]],\
         [[ -0.2, 4.2],[-0.20,4.20]],\
         [[ -0.2, 4.2],[-0.20,4.20]],\
         [[ -0.2, 4.2],[-0.20,4.20]],\
         [[ -0.2, 4.2],[-0.20,4.20]]]
yrng  = [[[ -2.0, 0.0],[ 0.65,1.05]],\
         [[ -8.4,-1.0],[ 0.45,1.05]],\
         [[ -5.0,-1.5],[ 0.35,1.05]],\
         [[ -7.0,-2.0],[ 0.35,1.05]],\
         [[-18.4,-9.6],[ 0.25,1.05]]]

owd = os.getcwd()
jn  = 0
for N in Nlst:
 jt   = 0
 Rlst = [str(k)+'.5' for k in range(N//2-1+1)]
 if(N==8): Rlst=['0.5','1.5','3.5']
 for dt in dtlst:                               
  Times      = np.zeros((fnsh[jt]+1,len(Rlst)))
  Energies   = np.zeros((fnsh[jt]+1,len(Rlst)))
  Fidelities = np.zeros((fnsh[jt]+1,len(Rlst)))
  jR = 0
  for R in Rlst:
   os.chdir('./'+str(N)+'_spin/R_'+R+'/dt_'+str(dt))
   subprocess.call(["sed -n "+str(strt[jn])+","+str(strt[jn]+fnsh[jt])+"p QITE.out > eofbeta.out"],shell=True)
   V = np.loadtxt('./eofbeta.out')
   print './'+str(N)+'_spin/R_'+R+'/dt_'+str(dt)
   print "sed -n "+str(strt[jn])+","+str(strt[jn]+fnsh[jt])+"p QITE.out > eofbeta.out"
   print V.shape
   Times[:,jR]      = V[:,0]
   Energies[:,jR]   = V[:,1]
   Fidelities[:,jR] = V[:,3]
   jR += 1
   os.chdir(owd)
  qite_plot(Rlst,Times,Energies,Fidelities,[Efci[jn],1],xrng[jn],yrng[jn],'qite_N_'+str(N)+'_dt_'+str(dt)+'.eps')
  jt += 1
 jn += 1
