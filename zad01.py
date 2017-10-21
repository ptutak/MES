# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 11:01:53 2017

@author: PiotrTutak
"""
import numpy
import scipy

print("Podaj L1,L2,L3,L4")

L=[float(x) for x in input().strip().split()]
print("k=",end='')
k=float(input().strip())
print("S=",end='')
S=float(input().strip())
print("q=",end='')
q=float(input().strip())
print("alfa=",end='')
alfa=float(input().strip())
print("tInf=",end='')
tInf=float(input().strip())

C=[k*S/l for l in L]

A=numpy.array([
   [C[0],-C[0],0,0,0],
   [-C[0],C[0]+C[1],-C[1],0,0],
   [0,-C[1],C[1]+C[2],-C[2],0],
   [0,0,-C[2],C[2]+C[3],-C[3]],
   [0,0,0,-C[3],C[3]+alfa*S]
   ])

P=numpy.array([
        q*S,
        0,
        0,
        0,
        -alfa*S*tInf
        ])

P=-P

t=scipy.linalg.solve(A,P)
print(t)