import numpy as np
import os,subprocess
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/mmotta/QITE/imag_time_for_qcomp/code_v4/')
import style as st

def qite_plot(RR,TT,EE,FF,FCI,xrng,yrng,fname):
 x0        = 7.00
 xfig      = 1.00*2*x0
 yfig      = 0.75*1*x0

 fig,pan = plt.subplots(1,2,figsize=(xfig,yfig))

 for js in range(2):
  pan[js].get_xaxis().set_tick_params(which='both',direction='in',labelsize=18)
  pan[js].get_yaxis().set_tick_params(which='both',direction='in',labelsize=18)

 pan[0].set_xlabel(r'$\displaystyle{\beta}$',fontsize=18)
 pan[0].set_ylabel(r'$\displaystyle{E_{QITE}(\beta)}, \, E_{Q-Lanz}(\beta)$',fontsize=18)
 pan[0].set_xlim(xrng[0])
 pan[0].set_ylim(yrng[0])

 pan[1].set_xlabel(r'$\displaystyle{\beta}$',fontsize=18)
 pan[1].set_ylabel(r'$\displaystyle{E_{QITE}(\beta)}, \, E_{Q-Lanz}(\beta)$',fontsize=18)
 pan[1].set_xlim(xrng[1])
 pan[1].set_ylim(yrng[1])

 for ip in range(2):
  u = np.arange(round(np.min(xrng[ip])),round(np.max(xrng[ip])),0.05)
  v = [FCI[ip]]*len(u)
  st.Palette[st.fci].dataplot(pan[ip],u,v,spline=True,connected=False)

 ipan = 0
 for jR in range(len(Rlst)):
  kR = jR + 10
  Lx = np.max(xrng[0])-np.min(xrng[0])
  Ly = np.max(yrng[0])-np.min(yrng[0])
  mask = st.equidistant_points(TT[:,jR],Lx,EE[:,jR],Ly*3.0,0.015)
  st.Palette[jR].dataplot(pan[ipan],TT[mask,jR],EE[mask,jR]    ,spline=True ,connected=False)
  st.Palette[jR].dataplot(pan[ipan],TT[mask,jR],EE[mask,jR]    ,spline=False,connected=False)
  st.Palette[jR].dataplot(pan[ipan],TT[mask,jR],EE[mask,jR]-5e3,spline=False,connected=True)
  # for the legend
  st.Palette[jR].dataplot(pan[0],TT[mask,jR],EE[mask,jR]-5e3,spline=True ,connected=False)
  st.Palette[jR].dataplot(pan[0],TT[mask,jR],EE[mask,jR]-5e3,spline=False,connected=False)
  st.Palette[jR].dataplot(pan[0],TT[mask,jR],EE[mask,jR]-5e3,spline=False,connected=True)

  Lx = np.max(xrng[1])-np.min(xrng[1])
  Ly = np.max(yrng[1])-np.min(yrng[1])
  SS   = [TT[2*i,jR] for i in range(FF.shape[0])]
  SS   = np.asarray(SS)
  mask = st.equidistant_points(SS,Lx,FF[:,jR],Ly*3.0,0.015)
  pt0,ps0                             = st.Palette[kR].pt,st.Palette[kR].ps
  st.Palette[kR].pt,st.Palette[kR].ps = '.',0
  st.Palette[kR].dataplot(pan[ipan],SS,FF[:,jR]             ,spline=False,connected=True)
  st.Palette[kR].pt,st.Palette[kR].ps = pt0,ps0
  st.Palette[kR].dataplot(pan[ipan],SS[mask],FF[mask,jR]    ,spline=False,connected=False)
  st.Palette[kR].dataplot(pan[ipan],SS[mask],FF[mask,jR]-5e3,spline=False,connected=True)

  # for the legend
  st.Palette[kR].dataplot(pan[0],SS,FF[:,jR]-5e3         ,spline=False,connected=True)
  st.Palette[kR].dataplot(pan[0],SS[mask],FF[mask,jR]-5e3,spline=False,connected=False)
  st.Palette[kR].dataplot(pan[0],SS[mask],FF[mask,jR]-5e3,spline=False,connected=True)

  ipan += 1

 x0 = -0.04
 dx =  2.28
 y0 =  1.002
 dy =  0.10
 order = [0,3,9,15,16]
 lgd   = st.make_legend(pan[0],order,ncol=5,bbox_to_anchor=(x0,y0,dx,dy),
         loc=3,mode="expand",borderaxespad=1,fontsize=17,handlelength=3,numpoints=2)

 fig.savefig(fname,format='eps',bbox_inches='tight')
 os.system('ps2epsi '+fname)
 os.system('mv '+fname+'i '+fname)

# ------------------------------------------------- #

Nlst  = [16]
strt  = [10370]
dtlst = [0.01]
fnsh  = [400]
lanz  = [199]
Efci  = [-28.569185]

xrng  = [[[ -0.2,  4.2],[-0.20, 4.20]]]
yrng  = [[[-30.0,-14.0],[-30.0,-14.0]]]

owd = os.getcwd()
jn  = 0
for N in Nlst:
 jt   = 0
 Rlst = [0.5,1.5]
 for dt in dtlst:                               
  Times      = np.zeros((fnsh[jt]+1,len(Rlst)))
  Energies   = np.zeros((fnsh[jt]+1,len(Rlst)))
  Fidelities = np.zeros((lanz[jt],len(Rlst)))
  jR = 0
  for R in Rlst:
   os.chdir('./R_'+str(R)+'/dt_'+str(dt))
   subprocess.call(["sed -n "+str(strt[jn])+","+str(strt[jn]+fnsh[jt])+"p QITE.out > eofbeta.out"],shell=True)
   V = np.loadtxt('./eofbeta.out')
   print './'+str(N)+'_spin/R_'+str(R)+'/dt_'+str(dt)
   print "sed -n "+str(strt[jn])+","+str(strt[jn]+fnsh[jt])+"p QITE.out > eofbeta.out"
   print V.shape
   Times[:,jR]      = V[:,0]
   Energies[:,jR]   = V[:,1]
   V = np.loadtxt('./test_lanz.out')
   print V.shape
   Fidelities[:,jR] = V[:,2]
   jR += 1
   os.chdir(owd)
  qite_plot(Rlst,Times,Energies,Fidelities,[Efci[jn],Efci[jn]],xrng[jn],yrng[jn],'qite_N_'+str(N)+'_dt_'+str(dt)+'.eps')
  #exit()
  jt += 1
 jn += 1
