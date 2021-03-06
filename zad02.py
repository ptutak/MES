# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 16:58:46 2017

@author: PiotrTutak
"""

import numpy as np

def f(x,y):
    return 2*(x**2)*(y**2)+6*x+5

def calkaPodwojna(f,ax,bx,ay,by,n,m):
    dx=(bx-ax)/n
    dy=(by-ay)/m
    x=ax
    calka=0.0
    for i in range(n):
        y=ay
        for j in range(m):
            f1=f(x,y)
            f2=f(x,y+dy)
            f3=f(x+dx,y)
            f4=f(x+dx,y+dy)
            calka+=(f1+f2+f3+f4)
            y+=dy
        x+=dx
    calka=calka*dx*dy*0.25
    return calka


def calkaGauss(f, kwadratury, wagi):
    calka=0.0
    for i in range(len(kwadratury)):
        for j in range(len(kwadratury)):
            calka+=f(kwadratury[i],kwadratury[j])*wagi[i]*wagi[j]   
    return calka


print(calkaPodwojna(f,-1.0,1.0,-1.0,1.0,100,100))
print(calkaGauss(f,[-np.sqrt(3)/3,np.sqrt(3)/3],[1.0,1.0]))
