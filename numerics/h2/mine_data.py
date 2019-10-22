import numpy as np
import os,subprocess
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/mmotta/QITE/imag_time_for_qcomp/code_v4/')
import style as st

fname='h2.eps'

nb=400
nr=54
my_data = np.genfromtxt('h2_surf.dat',skip_header=0)
R       = my_data[:,0]
b       = my_data[:,1]
E       = my_data[:,2]

print R[:10]
print b[:10]
print E[:10]

x0        = 7.00
xfig      = 1.00*1*x0
yfig      = 0.75*1*x0

fig,pan = plt.subplots(1,1,figsize=(xfig,yfig))

pan.get_xaxis().set_tick_params(which='both',direction='in',labelsize=18)
pan.get_yaxis().set_tick_params(which='both',direction='in',labelsize=18)

pan.set_xlabel(r'$\displaystyle{R \, [\AA]}$',fontsize=18)
pan.set_ylabel(r'$\displaystyle{E(\beta) \, \big[\mathrm{E_{Ha}}\big]}$',fontsize=18)
pan.set_xlim([-0.2,3.2])
pan.set_ylim([-1.2001,0.2])

u = [R[nb*x]      for x in range(nr)]
v = [E[nb-1+nb*x] for x in range(nr)]
st.Palette[st.fci].dataplot(pan,u,v,spline=False,connected=True)

mu = 16
for ib in [0,100,200,300]:
 u = [R[nb*x]    for x in range(nr)]
 v = [E[ib+nb*x] for x in range(nr)]
 st.Palette[mu].dataplot(pan,u,v,spline=False,connected=True)
 mu += 1

# figure label
plt.text(2.4,-0.03,r'$(d)$',fontsize=21)

x0 = -0.035
dx =  1.07
y0 =  1.00
dy =  0.10
order = [0,1,2,3,4]
lgd   = st.make_legend(pan,order,ncol=5,bbox_to_anchor=(x0,y0,dx,dy),\
        loc=3,mode="expand",borderaxespad=1,fontsize=13,handlelength=2.5,numpoints=2)

fig.savefig(fname,format='eps',bbox_inches='tight')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)

