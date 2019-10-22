import os,pylab,sys,itertools,subprocess
from   scipy.optimize    import curve_fit
from   matplotlib.ticker import FormatStrFormatter
from   scipy.linalg      import inv
sys.path.append('/home/mmotta/QITE/imag_time_for_qcomp/code_v4')
import style             as st
import numpy             as np
import matplotlib.pyplot as plt

# --------------------------------------- setting the figure

outfig='h2_surf.eps'

xfig = 7
yfig = 0.75*xfig
fig,pan = plt.subplots(1,1,figsize=(xfig,yfig))

dy   =  1.0
yoff =  0.0
ymin =  0.0
ymax =  5.0

st.set_canvas(pan,r'$\displaystyle R \left[\AA \right]$',
              [0.3,2.8],[0.3,0.8,1.3,1.8,2.3,2.8],[0.3,0.8,1.3,1.8,2.3,2.8], 
              r'$\displaystyle \beta \Big[\mathrm{E^{-1}_{Ha}} \Big]$',
              [ymin-yoff,ymax+yoff],np.arange(ymin,ymax+0.0001,dy),
              ['0.0','1.0','2.0','3.0','4.0','5.0'],20,20)

# --------------------------------------- plotting the data

st.colorplot(pan,'h2_surf.dat',[0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27],
             r'$\displaystyle E(R,\beta)-E_{exact}(R) [\mathrm{E_{Ha}}]$',fig,pan,shr=1.0,pad=0.025,smooth=3)

'''
# --- inset with cartoon

left, bottom, width, height = [0.2,0.8,0.32,0.32/2]
ax = fig.add_axes([left, bottom, width, height])

x  =  np.arange(0,1,0.01)
y1 = [ 0.008]*len(x)+0.07*(x-0.5)**2
y2 = [-0.008]*len(x)-0.07*(x-0.5)**2
ax.fill_between(x,y2,y1,facecolor='#808080')

x  =  np.arange(4,5,0.01)
y1 = [ 0.008]*len(x)+0.07*(x-4.5)**2
y2 = [-0.008]*len(x)-0.07*(x-4.5)**2
ax.fill_between(x,y2,y1,facecolor='#808080')

x  =  np.arange(1,4,0.01)
y1 = [ 0.005]*len(x)+0.003*(x-2.5)**2
y2 = [-0.005]*len(x)-0.003*(x-2.5)**2
ax.fill_between(x,y2,y1,facecolor='#DCDCDC')

x  =  np.arange(5,6,0.01)
y1 = [ 0.005]*len(x)+0.003*(x-6.5)**2
y2 = [-0.005]*len(x)-0.003*(x-6.5)**2
ax.fill_between(x,y2,y1,facecolor='#DCDCDC')

x  =  np.arange(-1,0,0.01)
y1 = [ 0.005]*len(x)+0.003*(x+1.5)**2
y2 = [-0.005]*len(x)-0.003*(x+1.5)**2
ax.fill_between(x,y2,y1,facecolor='#DCDCDC')

x = [0,1,4,5]
y = [0]*4
x = np.asarray(x)
y = np.asarray(y)

nx=50
for ix in range(0,nx,2):
 dx=0.025*float(ix)/float(nx)
 dy=0.015*float(ix)/float(nx)
 c=str(0.2+0.8*(float(ix)/float(nx))**2)
 s=250*float(nx-ix)/float(nx)
 ax.scatter(x+dx,y+dy,s=s,facecolors=c,edgecolors=None)

from matplotlib.pyplot import text

text(0.4, 0.06,r'$\displaystyle R$')
ax.arrow(0, 0.05, 1,0,head_width=0.007,head_length=0.1,fc='k',ec='k',width=0.0004)
ax.arrow(1, 0.05,-1,0,head_width=0.007,head_length=0.1,fc='k',ec='k',width=0.0004)
text(4.4, 0.06,r'$\displaystyle R$')
ax.arrow(4, 0.05, 1,0,head_width=0.007,head_length=0.1,fc='k',ec='k',width=0.0004)
ax.arrow(5, 0.05,-1,0,head_width=0.007,head_length=0.1,fc='k',ec='k',width=0.0004)
text(2.3,-0.08,r'$\displaystyle R^\prime$')
ax.arrow(1,-0.05, 3,0,head_width=0.01,head_length=0.1,fc='k',ec='k')
ax.arrow(4,-0.05,-3,0,head_width=0.01,head_length=0.1,fc='k',ec='k')

st.set_canvas(ax,'',[-1,6],[],[],
                 '',[-0.1,0.1],[],[],5,5)
'''

fig.savefig(outfig,format='eps',bbox_inches='tight')
os.system('ps2epsi '+outfig)
os.system('mv '+outfig+'i '+outfig)

exit()

