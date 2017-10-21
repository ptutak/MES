# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 11:01:53 2017

@author: PiotrTutak
"""
import numpy

print("Podaj L1,L2,L3,L4")

L=[float(x) for x in input().strip().split()]
print("k=")
k=float(input().strip())
print("S=")
S=float(input().strip())
print("q=")
q=float(input().strip())
print("alfa=")
alfa=float(input().strip())
print("tInf=")
tInf=float(input().strip())

C=[k*S/l for l in L]

A=numpy.array([
   [C[1],-C[1],0,0,0],
   [-C[1],C[1]+C[2],-C[2],0,0],
   [0,-C[2],C[2]+C[3],-C[3],0],
   [0,0,-C[3],C[3]+C[4],-C[4]],
   [0,0,0,-C[4],C[4]+alfa*S]
   ])

B=numpy.array([
        -q*S,
        0,
        0,
        0,
        alfa*S*tInf
        ])

