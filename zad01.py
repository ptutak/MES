# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 11:01:53 2017

@author: PiotrTutak
"""
import numpy as np
import scipy.linalg as lg
import matplotlib.pyplot as plt
print("Podaj L1 L2 L3 L4")

L=[float(x) for x in input().strip().split()]

print('Podaj k S q alfa tInf')

k,S,q,alfa,tInf=(float(x) for x in input().strip().split())

C=[k*S/l for l in L]
L=[sum(L[:i] for i in range(1,len(L)+1)]
       
A=np.array([
   [C[0],-C[0],0,0,0],
   [-C[0],C[0]+C[1],-C[1],0,0],
   [0,-C[1],C[1]+C[2],-C[2],0],
   [0,0,-C[2],C[2]+C[3],-C[3]],
   [0,0,0,-C[3],C[3]+alfa*S]
   ])

P=np.array([
        q*S,
        0,
        0,
        0,
        -alfa*S*tInf
        ])

P=-P

#t=np.linalg.solve(A,P)
t=lg.solve(A,P)
print(t)
plt.plot(L,t)
plt.show()
