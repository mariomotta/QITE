import numpy as np
import time

def Int2Bas(n,b,nbit):
 if(n==0): return [0]*nbit
 x=[]
 while(n):
  x.append(int(n%b))
  n//=b
 return [0]*(nbit-len(x))+x[::-1]

def Bas2Int(x,b):
 nbit=len(x)
 z = [b**(nbit-1-i) for i in range(nbit)]
 return np.dot(z,x) 

def Psi2Str(x):
 s='|'
 for i in x: s+=str(i)
 return s+'>'

def str_op(i):
 if(i==0): return 'I'
 if(i==1): return 'X'
 if(i==2): return 'Y'
 if(i==3): return 'Z'

def Opp2Str(x):
 s=''
 for i in x: s+=str_op(i)
 return s

def Lst2Str(x):
 s=''
 for i in x: s+=str(i)
 return s

if __name__ == "__main__":
 nbmax=10
 b=4
 v=[]
 for nbit in range(1,nbmax):
  t0 = time.time()
  for n0 in range(b**nbit):
   assert(n0==Bas2Int(Int2Bas(n0,b,nbit),b))
   if(nbit==4):
    print n0,Int2Bas(n0,b,nbit)
  t1 = time.time()
  dt = t1-t0
  v.append(dt)

 import matplotlib.pyplot as plt
 plt.plot(range(1,nbmax),v,'bo-')
 plt.xlabel("bits")
 plt.ylabel("t [s]")
 plt.show()
