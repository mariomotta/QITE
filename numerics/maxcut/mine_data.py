import numpy as np
import os,subprocess
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/mmotta/QITE/imag_time_for_qcomp/code_v4/')
import style as st

def qite_plot(RR,TT,EE,FF,FCI,xrng,yrng,gamma,fname):
 x0        = 7.00
 xfig      = 1.00*1*x0
 yfig      = 0.75*1*x0

 fig,pan = plt.subplots(1,1,figsize=(xfig,yfig))

 pan.get_xaxis().set_tick_params(which='both',direction='in',labelsize=18)
 pan.get_yaxis().set_tick_params(which='both',direction='in',labelsize=18)

 '''
 pan.set_xlabel(r'$\displaystyle{\beta}$',fontsize=14)
 pan.set_ylabel(r'$\displaystyle{C(\beta)}$',fontsize=14)
 pan.set_xlim(xrng[0])
 pan.set_ylim(yrng[0])
 '''

 pan.set_xlabel(r'$\displaystyle{\beta}$',fontsize=18)
 pan.set_ylabel(r'$P(C=C_{max})$',fontsize=18) 
 pan.set_xlim(xrng[1])
 pan.set_ylim(yrng[1])
 pan.set_yticks([0,0.2,0.4,0.6,0.8,1.0,1.2000000001],[r'$\displaystyle{\phantom{-}0.0}$',
                                             r'$\displaystyle{\phantom{-}0.2}$',
                                             r'$\displaystyle{\phantom{-}0.4}$',
                                             r'$\displaystyle{\phantom{-}0.6}$',
                                             r'$\displaystyle{\phantom{-}0.8}$',
                                             r'$\displaystyle{\phantom{-}1.0}$',
                                             r'$\displaystyle{\phantom{-}1.2]$'])

 #for ip in range(2):
 ip=1
 u = np.arange(0,16,0.1) #round(np.min(xrng[ip])),round(np.max(xrng[ip])),0.05)
 v = [FCI[ip]]*len(u)
 st.Palette[st.fci].dataplot(pan,u,v,spline=False,connected=True)

 for jR in range(len(Rlst)):
  '''
  Lx = np.max(xrng[0])-np.min(xrng[0])
  Ly = np.max(yrng[0])-np.min(yrng[0])
  mask = st.equidistant_points(TT[:,jR],Lx,EE[:,jR],Ly*3.0,0.01)
  st.Palette[jR].dataplot(pan[0],TT[   :,jR],EE[   :,jR]    ,spline=True ,connected=False)
  st.Palette[jR].dataplot(pan[0],TT[mask,jR],EE[mask,jR]    ,spline=False,connected=False)
  st.Palette[jR].dataplot(pan[0],TT[mask,jR],EE[mask,jR]-5e3,spline=False,connected=True)
  '''
  Lx = np.max(xrng[1])-np.min(xrng[1])
  Ly = np.max(yrng[1])-np.min(yrng[1])
  mask = st.equidistant_points(TT[:,jR],Lx,FF[:,jR],Ly*3.0,0.006)
  st.Palette[jR].dataplot(pan,TT[   :,jR],FF[   :,jR]    ,spline=False,connected=True)
  st.Palette[jR].dataplot(pan,TT[mask,jR],FF[mask,jR]    ,spline=False,connected=False)
  st.Palette[jR].dataplot(pan,TT[mask,jR],FF[mask,jR]-5e3,spline=False,connected=True)

 # graph inset

 from matplotlib.patches import Polygon
 ins      = plt.axes([0.48,0.10,0.40,0.40],aspect='equal')
 nv,E,cut = gamma

 ins.axis("off")
 ins.set_xlim([-1.2,1.2])
 ins.set_ylim([-1.2,1.2])
 u = np.linspace(-np.pi,np.pi,100)
 v = np.sin(u)
 u = np.cos(u)
 st.Palette[st.link].dataplot(ins,u, v,spline=False,connected=True)
 st.Palette[st.link].dataplot(ins,u,-v,spline=False,connected=True)

 for (iv,jv) in E:
  phii  = 2*np.pi*float(iv)/float(nv)
  phij  = 2*np.pi*float(jv)/float(nv)
  xi,xj = np.cos(phii),np.cos(phij)
  yi,yj = np.sin(phii),np.sin(phij)
  u     = np.linspace(xi,xj,10)
  v     = np.linspace(yi,yj,10)
  idx   = st.linkin
  print iv,jv,iv in cut,jv in cut
  if((iv in cut) and (not jv in cut)): idx = st.linkbtw
  if((jv in cut) and (not iv in cut)): idx = st.linkbtw
  st.Palette[idx].dataplot(ins,u,v,spline=False,connected=True)

 for iv in range(nv):
  if(iv in cut): idx = st.vertexin
  else:          idx = st.vertexout
  phi = 2*np.pi*float(iv)/float(nv)
  st.Palette[idx].dataplot(ins,[np.cos(phi)],[np.sin(phi)],spline=False,connected=False)

 # figure label
 plt.text(1.05,3.15,r'$(b)$',fontsize=21)

 x0 = -0.04
 dx =  1.08
 y0 =  1.002
 dy =  0.10
 order = [0]+[3*(x+1) for x in range(len(Rlst))]
 lgd   = st.make_legend(pan,order,ncol=range(len(Rlst)),bbox_to_anchor=(x0,y0,dx,dy),\
         loc=3,mode="expand",borderaxespad=1,fontsize=14,handlelength=3,numpoints=2)

 fig.savefig(fname,format='eps',bbox_inches='tight')
 os.system('ps2epsi '+fname)
 os.system('mv '+fname+'i '+fname)

# ------------------------------------------------- #

Nlst  = [4,6]
strt  = [5,5]
dtlst = [0.01]
fnsh  = [800,1600]
Efci  = [3,5]
Graph = [(4,[(0,2),(1,2),(1,3)],[2,3]),\
         (6,[(0,3),(1,4),(2,3),(2,4),(2,5),(4,5)],[1,3,5])]

xrng  = [[[-0.40, 8.40],[-0.40, 8.40]],\
         [[-0.50,16.50],[-0.50,16.50]]]
yrng  = [[[ 1.40, 3.10],[-0.05, 1.05]],\
         [[ 2.85, 5.15],[-0.00, 1.00]]]

owd = os.getcwd()
jn  = 0
for N in Nlst:
 jt   = 0
 Rlst = [str(k)+'.5' for k in range(3)]
 for dt in dtlst:
  Times    = np.zeros((fnsh[jn]+1,len(Rlst)))
  Energies = np.zeros((fnsh[jn]+1,len(Rlst)))
  Probabs  = np.zeros((fnsh[jn]+1,len(Rlst)))
  jR = 0
  for R in Rlst:
   os.chdir('./'+str(N)+'_vertices/R_'+R+'/dt_'+str(dt))
   subprocess.call(["sed -n "+str(strt[jn])+","+str(strt[jn]+fnsh[jn])+"p QITE.out > eofbeta.out"],shell=True)
   V = np.loadtxt('./eofbeta.out')
   print './'+str(N)+'_vertices/R_'+R+'/dt_'+str(dt)
   print "sed -n "+str(strt[jn])+","+str(strt[jn]+fnsh[jn])+"p QITE.out > eofbeta.out"
   print V.shape
   Times[:,jR]    =  V[:,0]
   Energies[:,jR] = -V[:,1]
   Probabs[:,jR]  =  V[:,3]
   jR += 1
   os.chdir(owd)
  qite_plot(Rlst,Times,Energies,Probabs,[Efci[jn],1],xrng[jn],yrng[jn],Graph[jn],'qite_N_'+str(N)+'_dt_'+str(dt)+'.eps')
  jt += 1
 jn += 1
