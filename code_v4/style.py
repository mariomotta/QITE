import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import pylab             as pl
import matplotlib.colors as col

import os,math

from   matplotlib           import rc
from   scipy                import interpolate
from   matplotlib.lines     import Line2D
from   mpl_toolkits.mplot3d import Axes3D
from   scipy.interpolate    import griddata,interp1d

class Style:

 def __init__(self):
  self.name  = 'empty'
  self.title = 'empty'
  self.color = 'empty'
  self.pt    = 'empty'
  self.ps    = 'empty'
  self.fs    = 'empty'
  self.lt    = 'empty'
  self.lw    = 'empty'
  self.lw2   = 'empty'

 def fill(self,name,title,color,pt,ps,fs,lt,lw,lw2):
  self.name  = name
  self.title = title
  self.color = color
  self.pt    = pt
  self.ps    = ps
  self.fs    = fs
  self.lt    = lt
  self.lw    = lw
  self.lw2   = lw2

 def dataplot(self,target,x,y,dx=None,dy=None,spline=None,connected=None,capsize=None,logscale=None):

  if(logscale=='xy'): target.set_xscale("log"); target.set_yscale("log")

  if(capsize==None): capsize=5
  mec='black'
  if(self.fs=='none'): mec=self.color

  if(dy is None):

   if(spline==False):
    if(connected==False):
     target.plot(x,y,label=self.title,marker=self.pt,markersize=self.ps,fillstyle=self.fs,linestyle=' ', \
                 color=self.color,markeredgecolor=mec,mew=self.lw,linewidth=self.lw2)
    else:
     target.plot(x,y,label=self.title,marker=self.pt,markersize=self.ps,fillstyle=self.fs,linestyle=self.lt, \
                 color=self.color,markeredgecolor=mec,mew=self.lw,linewidth=self.lw2)
   else:

    f2 = interp1d(x,y,kind='cubic')
    u  = np.linspace(min(x),max(x),2000)
    v  = f2(u)

    target.plot(u,v,label=self.title,linestyle=self.lt, \
                color=self.color,markeredgecolor=mec,mew=self.lw,linewidth=self.lw2)
    #target.plot(x,y,label=self.title,marker=self.pt,markersize=self.ps,fillstyle=self.fs,linestyle=' ', \
    #            color=self.color,markeredgecolor=mec,mew=self.lw,linewidth=self.lw2)

  else:

   if(spline==False):
    if(connected==False):
     (_, caps, _) = target.errorbar(x,y,label=self.title,marker=self.pt,markersize=self.ps,fillstyle=self.fs,linestyle=' ', \
                                    color=self.color,yerr=dy,mew=self.lw,markeredgecolor=mec,linewidth=0,capsize=capsize,elinewidth=self.lw2)
     for cap in caps:
          cap.set_markeredgewidth(self.lw2)
    else:
     (_, caps, _) = target.errorbar(x,y,label=self.title,marker=self.pt,markersize=self.ps,fillstyle=self.fs,linestyle=self.lt, \
                    color=self.color,yerr=dy,mew=self.lw,markeredgecolor=mec,linewidth=self.lw2,capsize=capsize,elinewidth=self.lw2)
     for cap in caps:
          cap.set_markeredgewidth(self.lw2)
   else:

    f2 = interp1d(x,y,kind='cubic')
    u  = np.linspace(min(x),max(x),2000)
    v  = f2(u)

    target.plot(u,v,label=self.title, \
                linestyle=self.lt,color=self.color,mew=self.lw,linewidth=self.lw2)
    (_, caps, _) = target.errorbar(x,y,label=self.title,marker=self.pt,markersize=self.ps,fillstyle=self.fs,linestyle=' ', \
                   color=self.color,yerr=dy,markeredgecolor=mec,mew=self.lw,linewidth=self.lw2,capsize=capsize,elinewidth=self.lw2)
    for cap in caps:
         cap.set_markeredgewidth(self.lw2)

 def funplot(self,target,u,v):
     target.plot(u,v,label=self.title,linestyle=self.lt,color=self.color,mew=self.lw,linewidth=self.lw2)

#---------------------------------------------------

def make_legend(pan,order,ncol,bbox_to_anchor,loc,mode,borderaxespad,fontsize,handlelength,numpoints=None):

 handles, labels = pan.get_legend_handles_labels()

 handles2 = []
 labels2  = []

 for i in range(len(order)):
  handles2.append(handles[order[i]])
  labels2.append(labels[order[i]])

 if(numpoints==None):
  return pan.legend(handles2,labels2,ncol=ncol,fancybox=True,shadow=True,bbox_to_anchor=bbox_to_anchor,\
                    loc=loc,mode=mode,borderaxespad=borderaxespad,fontsize=fontsize,handlelength=handlelength)
 else:
  return pan.legend(handles2,labels2,ncol=ncol,fancybox=True,shadow=True,bbox_to_anchor=bbox_to_anchor,\
                    loc=loc,mode=mode,borderaxespad=borderaxespad,fontsize=fontsize,handlelength=handlelength,numpoints=numpoints)

def set_canvas(pan,xlab,xlim,xtic,xticlab,ylab,ylim,ytic,yticlab,ticfontsize,labfontsize):
 pan.set_xlabel(xlab,fontsize=ticfontsize)
 pan.set_xlim(xlim)
 pan.set_xticks(xtic)
 pan.set_xticklabels(xticlab)
 pan.set_ylabel(ylab,fontsize=ticfontsize)
 pan.set_ylim(ylim)
 pan.set_yticks(ytic)
 pan.set_yticklabels(yticlab)
 pan.tick_params(axis='both',labelsize=labfontsize,direction='in')

def mask(ymin,ymax,y):
 for iy in range(len(y)):
  if(y[iy]<ymin)or(y[iy]>ymax):
     y[iy]='nan'
 return y

def colorplot(ax,fname,levels,clabel,fig,pan=None,shr=None,pad=None,smooth=None):

 import scipy.ndimage
 from mpl_toolkits.axes_grid1 import make_axes_locatable

 my_data = np.genfromtxt(fname,skip_header=0)

 X       = my_data[:,0]
 Y       = my_data[:,1]
 Z       = my_data[:,2]   

 print "extrema of Z",min(Z),max(Z)

 xi      = np.linspace(min(X),max(X),1000)
 yi      = np.linspace(min(Y),max(Y),1000)
 zi      = griddata((X,Y),Z+1e-6,(xi[None,:],yi[:,None]),method='cubic')

 print "min of Z on grid",np.min(np.min(zi, axis=1), axis=0)

 xig,yig = np.meshgrid(xi,yi)
 surf    = ax.contour(xig, yig, zi, levels, linewidths=0.5,linestyle=':', colors='k')
 cax     = ax.contourf(xig, yig, zi, np.linspace(min(levels),max(levels),100), cmap=plt.cm.Spectral)
 
 if(pan==None):
    cbar = plt.colorbar(cax,ticks=levels)
 else:
    if(pad==None):
       cbar = plt.colorbar(cax,ax=pan,ticks=levels,shrink=shr)
    else:
       cbar = plt.colorbar(cax,ax=pan,ticks=levels,shrink=shr,pad=pad)
 cbar.ax.set_ylabel(clabel,size=18,labelpad=15)
 cbar.ax.tick_params(axis='y',direction='in',labelsize=18)
# text = cbar.ax.yaxis.label
 

def equidistant_points(x,Lx,y,Ly,d):
 mask   = [0]
 lngth2 = 0.0
 for i in range(1,x.shape[0]):
  lngth2 += (x[i]-x[i-1])**2/Lx**2+(y[i]-y[i-1])**2/Ly**2
  if(lngth2>d**2):
   mask.append(i)
   lngth2 = 0.0
 mask.append(x.shape[0]-1)
 return mask

#---------------------------------------------------

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

#---------------------------------------------------   

Nstyles = 200
d       =  11
lw      = 0.2
lwb     = 0.2
lw2     = 2.0
lw2b    = 2.0
ps      =  12

red    = '#E42F0C'
orange = '#F28500'
yellow = '#fdb827'
green  = '#78bb17'
blue   = '#146eff'
purple = '#A26BA2'
black  = '#000000'
white  = '#ffffff'
gray   = '#999999'

Palette = [ Style() for i in range(Nstyles) ]

#Palette[0].fill('c0',r'$\displaystyle D=2$',  red   ,'o',8,'full','-.',lw,lw2b)
#Palette[1].fill('c1',r'$\displaystyle D=4$',  green ,'D',7,'full','--',lw,lw2b)
#Palette[2].fill('c2',r'$\displaystyle D=6$',  blue  ,'s',7,'full',':', lw,lw2b)
#Palette[3].fill('c3',r'$\displaystyle D=8$',  orange,'^',7,'full','-', lw,lw2b)
#Palette[4].fill('c4',r'$\displaystyle exact$',black ,'.',0,'full','-.',lw,lw2b)
Palette[0].fill('c0',r'$\displaystyle D=2$',  red   ,'.',0,'full','--', lw,lw2b)
Palette[1].fill('c1',r'$\displaystyle D=4$',  green ,'.',0,'full','-.',lw,lw2b+1)
Palette[2].fill('c2',r'$\displaystyle D=6$',  blue  ,'.',0,'full',':', lw,lw2b+1)
Palette[3].fill('c3',r'$\displaystyle D=8$',  yellow,'.',0,'full','--',lw,lw2b)
Palette[4].fill('c3',r'$\displaystyle D=10$', yellow,'.',0,'full','-.',lw,lw2b)
Palette[5].fill('c4',r'$\displaystyle \mathrm{exact}$',black ,'.',0,'full','-', lw,lw2b)
fci = 5

Palette[ 6].fill('c5','',black,'o',9,'full','-',  lw,lw2b)
Palette[ 7].fill('c5','',white,'o',9,'full','-',  lw,lw2b)
Palette[ 8].fill('c4','',gray ,'.',0,'full',':',  lw,  lw)
Palette[ 9].fill('c4','',black,'.',0,'full','-',  lw,  lw)
Palette[10].fill('c4','',black,'.',0,'full','-',2*lw,2*lw)
vertexin  = 6
vertexout = 7
link      = 8
linkin    = 9
linkbtw   = 10

Palette[11].fill('c0',r'$\displaystyle D=2, \mathrm{QLanczos}$',  red   ,'.',0,'none','-',lw2b,lw2b)
Palette[12].fill('c1',r'$\displaystyle D=4, \mathrm{QLanczos}$',  green ,'.',0,'none','-',lw2b,lw2b)
Palette[13].fill('c2',r'$\displaystyle D=6, \mathrm{QLanczos}$',  blue  ,'.',0,'none','-',lw2b,lw2b)
Palette[14].fill('c3',r'$\displaystyle D=8, \mathrm{QLanczos}$',  orange,'.',0,'none','-',lw2b,lw2b)
Palette[15].fill('c4',r'$\displaystyle \mathrm{exact}$',black ,'.',0,'none','-',lw2b,lw2b)

Palette[16].fill('c0',r'$\displaystyle \beta=0$',red   ,'.',0,'full','--', lw,lw2b)
Palette[17].fill('c1',r'$\displaystyle \beta=1$',green ,'.',0,'full','-.',lw,lw2b+1)
Palette[18].fill('c2',r'$\displaystyle \beta=2$',blue  ,'.',0,'full',':', lw,lw2b+1)
Palette[19].fill('c3',r'$\displaystyle \beta=3$',orange,'.',0,'full','--',lw,lw2b)
