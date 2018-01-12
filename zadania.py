#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 07:18:21 2018

@author: ptutak
"""

from meshmes import *
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib import cm as cm


n1=ShapeFunc(lambda xsi,eta:0.25*(1-xsi)*(1-eta),
             {'Xsi':lambda xsi,eta:-0.25*(1-eta),
              'Eta':lambda xsi,eta:-0.25*(1-xsi)})
n2=ShapeFunc(lambda xsi,eta:0.25*(1+xsi)*(1-eta),
             {'Xsi':lambda xsi,eta:0.25*(1-eta),
              'Eta':lambda xsi,eta:-0.25*(1+xsi)})
n3=ShapeFunc(lambda xsi,eta:0.25*(1+xsi)*(1+eta),
             {'Xsi':lambda xsi,eta:0.25*(1+eta),
              'Eta':lambda xsi,eta:0.25*(1+xsi)})
n4=ShapeFunc(lambda xsi,eta:0.25*(1-xsi)*(1+eta),
             {'Xsi':lambda xsi,eta:-0.25*(1+eta),
              'Eta':lambda xsi,eta:0.25*(1-xsi)})
#    x=Compute([n1,n2,n3,n4],dict([('X','Xsi'),('Y','Eta')]),[-0.7745966692414834,0.0,0.7745966692414834],[0.5555555555555556,0.8888888888888888,0.5555555555555556])
x=Compute([n1,n2,n3,n4],dict([('X','Xsi'),('Y','Eta')]),[-1.0/np.sqrt(3),1.0/np.sqrt(3)],[1.0,1.0])


fileName='mroz'

globalData=loadData('data-'+fileName+'.yml')
print(globalData)
g=Grid(globalData['B'],globalData['H'],globalData['nB'],globalData['nH'],globalData['t0'],globalData['edges'])
g.setElemPhysParams(k=globalData['k'],ro=globalData['ro'],c=globalData['c'])
x.compGridElemPoints(g)
tau=0.0
results=[]
start=time.clock()
while tau<globalData['tau']:
    for k,v in sorted(globalData['tInf'].items()):
        if tau<k:
            temp=v
            break;
    t1=x.compGridTempPoints(g,alfa=globalData['alfa'],tInf=temp,dTau=globalData['dTau'])
    g.updateNodes(t1,'t')
    tau+=globalData['dTau']
    if (tau%7200.0==0):
        t1=np.transpose(np.reshape(t1,(globalData['nB'],globalData['nH'])))
        results.append((t1,tau,temp[-1],time.clock()-start))
        print(tau,temp[-1],time.clock()-start)
print("time: ",time.clock()-start)
f,axarr=plt.subplots(len(results))
f.set_size_inches((12,6*len(results)))
j=0
for r,tau,temp,tauR in results:
    axarr[j].set_xlim(0,53)
    axarr[j].set_xlabel('grubość [cm]')
    axarr[j].set_ylim(0,1)
    axarr[j].axvline(1.0)
    axarr[j].axvline(26.0)
    axarr[j].axvline(41.0)
    axarr[j].axvline(53.0)
    axarr[j].set_title('tau:  '+str(int(tau/3600.0//24.0))+'d '+str(int(tau/3600.0%24.0))+'h, temp: '+str(temp))
    im=axarr[j].imshow(r,cmap=cm.jet,interpolation='bicubic',aspect='auto',origin='lower',vmax=25,vmin=-10)
    f.colorbar(im,ax=axarr[j])
    j+=1
plt.tight_layout()
plt.savefig('results-'+fileName+'.png')
plt.show()
with open('results-'+fileName+'.txt','w') as fileOut:
    for r,tau,temp,tauR in results:
        print(tau,temp,tauR,file=fileOut)
        print(printSeq2(r,5),file=fileOut)