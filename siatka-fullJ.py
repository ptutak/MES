# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:26:28 2017

@author: Piotrek
"""
import yaml
import numpy as np
import scipy.linalg as lg
from collections import abc
import time

def loadData(fileName):
    data = {}
    with open(fileName) as file:
        data = yaml.load(file)
    for x in data:
        if isinstance(data[x],str):
            data[x]=eval(data[x])
    return data

def printSeq(seq,prec):
    formatstr="{0: >"+str(prec+6)+"."+str(prec)+"f}"
    res=""
    for x in seq:
        res+=formatstr.format(float(x))
    return res

def printSeq2(seq,prec):
    formatstr="{0: >"+str(prec+6)+"."+str(prec)+"f}"
    res=""
    for x in seq:
        for y in x:
            res+=formatstr.format(float(y))
        res+="\n"
    return res

class ShapeFunc:
    def __init__(self, f, dFd):
        self.f = f
        self.dFd = dFd
        self.varNumber = len(dFd)
        self.vars = dFd.keys()

    def __call__(self, *args):
        if len(args) != self.varNumber:
            raise TypeError('Wrong number of variables')
        return self.f(*args)

    def __getitem__(self, index):
        return self.dFd[index]


class Node:
    def __init__(self,X,Y,t=0,edge=False):
        self.X=X
        self.Y=Y
        self.edge=edge
        self.t=t
    def __repr__(self):
        return "({0},{1})".format(self.X,self.Y)
    def __str__(self):
        return "({0},{1},t={2},edge={3})".format(self.X,self.Y,self.t,self.edge)
    def __setitem__(self,index,value):
        self.__dict__[index]=value
    def __getitem__(self,index):
        return self.__dict__[index]


class Element:
    def __init__(self,nodes,ids):
        self.nodes=nodes
        self.ids=ids
        self.edge=False
        self.surface=[dict() for i in range(4)]
        for i in range(4):
            self.surface[i]['edge']=all([nodes[(i+3)%4].edge,nodes[(i+4)%4].edge])
            if self.surface[i]['edge']:
                self.edge=True
            self.surface[i]['len']=((nodes[(i+3)%4].X - nodes[(i+4)%4].X)**2 + (nodes[(i+3)%4].Y - nodes[(i+4)%4].Y)**2)**0.5
        self.points=None
        self.H=None
        self.P=None
        self.C=None
#        self.Ct0=None
    def __repr__(self):
        rep="{3!r} {2!r}\n{0!r} {1!r}\n{4!r}\n".format(*self.ids,self.surface)
        return rep
    def __str__(self):
        strg="{3!s} {2!s}\n{0!s} {1!s}\n{4!r}\n".format(*self.nodes,self.surface)
        return strg
    def __getitem__(self,index):
        if type(index)==str:
            tmp=[]
            for n in self.nodes:
                tmp.append(n.__dict__[index])
            return np.array(tmp)
        return self.nodes[index]
    def __setitem__(self,index,value):
        self.nodes[index]=value
    def __len__(self):
        return len(self.nodes)


class Grid:
    def __init__(self,B,H,nB,nH,t0=0,edges=False,x=0,y=0):
        self.nB=nB
        self.nH=nH
        self.B=B
        self.H=H
        self.db=self.B/(self.nB-1)
        self.dh=self.H/(self.nH-1)
        self.nn=nB*nH
        self.ne=(nH-1)*(nB-1)
        nodes=[]
        if isinstance(t0,abc.Sequence):
            t0=iter(t0)
        else:
            t0=iter([t0 for x in range(self.nn)])
        if isinstance(edges,abc.Sequence):
            edges=iter(edges)
        for i in range(nB):
            for j in range(nH):
                if edges:
                    nextEdge=False
                    try:
                        nextEdge=next(edges)
                    except (StopIteration):
                        edges=False
                    nodes.append(Node(x+self.db*i,y+self.dh*j,next(t0),nextEdge))
                else:
                    if (i==0 or i==nB-1 or j==0 or j==nH-1):
                        nodes.append(Node(x+self.db*i,y+self.dh*j,next(t0),True))
                    else:
                        nodes.append(Node(x+self.db*i,y+self.dh*j,next(t0)))
        self.nodes=nodes
        elements=[]
        for i in range(nB-1):
            for j in range(nH-1):
                elements.append(Element([nodes[i*nH+j],nodes[(i+1)*nH+j],nodes[(i+1)*nH+j+1],nodes[i*nH+j+1]],[i*nH+j,(i+1)*nH+j,(i+1)*nH+j+1,i*nH+j+1]))
        self.elements=elements
    def __getitem__(self,index):
        if type(index)==str:
            tmp=[]
            for n in self.nodes:
                tmp.append(n.__dict__[index])
            return np.array(tmp)
        return self.elements[index]
    def __len__(self):
        return len(self.elements)
    def __call__(self,index):
        return self.nodes[index]
    def printNodeAttrs(self,attr,prec=None):
        if prec:
            formatStr="{0: >"+str(prec+6)+"."+str(prec)+"f}"
        else:
            formatStr="{0!r:^5}\t"
        res=""
        for j in range(self.nH-1,-1,-1):
            for i in range(self.nB):
                res+=formatStr.format(self.nodes[i*self.nH+j][attr])
            res+="\n\n"
        return res
    def updateNodes(self,values,attr):
        values=iter(values)
        for n in self.nodes:
            n[attr]=next(values)

class Compute:
    def __init__(self,N,variables,gaussQ,gaussW):
        self.N=N
        self.lenN=len(N)
        self.variables=variables
        self.globVar=set(variables.keys())
        self.locVar=set(variables.values())
        self.q=gaussQ
        self.w=gaussW
        self.lenGauss=len(self.q)
        self.points=[dict([('q',(gaussQ[i],gaussQ[j])),('w',gaussW[i]*gaussW[j]),('N',np.array([n(gaussQ[i],gaussQ[j]) for n in self.N]))]) for i in range(len(gaussQ)) for j in range(len(gaussQ))]
        
        for point in self.points:
            point['N^2']=np.matmul(np.transpose(np.array([point['N']])),np.array([point['N']]))
        for point in self.points:
            for v in self.locVar:
                listdNd=[]
                for n in self.N:
                    listdNd.append(n[v](*point['q']))
                point['dNd'+v]=np.array(listdNd)
        self.lenPoints=len(self.points)
        surfaces=[[dict() for y in range(self.lenGauss)] for x in range(4)]
        for j in range(self.lenGauss):
            nSurf0=[]
            nSurf1=[]
            nSurf2=[]
            nSurf3=[]
            for n in self.N:
                nSurf0.append(n(-1.0,self.q[j]))
                nSurf1.append(n(self.q[j],-1.0))
                nSurf2.append(n(1.0,self.q[j]))
                nSurf3.append(n(self.q[j],1.0))
            surfaces[0][j]['N']=np.array(nSurf0)
            surfaces[1][j]['N']=np.array(nSurf1)
            surfaces[2][j]['N']=np.array(nSurf2)
            surfaces[3][j]['N']=np.array(nSurf3)
            for i in range(4):
                surfaces[i][j]['w']=self.w[j]
                surfaces[i][j]['N^2']=np.matmul(np.transpose(np.array([surfaces[i][j]['N']])),np.array([surfaces[i][j]['N']]))
        self.surface=surfaces
    def compElementPoints(self,element):
        element.points=[]
        for i in range(self.lenPoints):
            element.points.append(dict())
            for vG in self.globVar:
                for vL in self.locVar:
                    element.points[i]['d'+vG+'d'+vL]=np.dot(self.points[i]['dNd'+vL],element[vG])
            J=np.array([
                    [element.points[i]['dXdXsi'],element.points[i]['dYdXsi']],
                    [element.points[i]['dXdEta'],element.points[i]['dYdEta']]
                    ])
            element.points[i]['J^-1']=np.array([
                    [J[1,1],-J[0,1]],
                    [-J[1,0],J[0,0]]
                    ])
            detJ=lg.det(J)
            element.points[i]['detJ']=detJ
            element.points[i]['1/detJ']=1.0/detJ
            dNdX=[]
            dNdY=[]
            for j in range(self.lenN):
                res=element.points[i]['1/detJ']*np.dot(element.points[i]['J^-1'],np.array([self.points[i]['dNdXsi'][j],self.points[i]['dNdEta'][j]]))
                dNdX.append(res[0])
                dNdY.append(res[1])
            element.points[i]['dNdX']=np.array(dNdX)
            element.points[i]['dNdY']=np.array(dNdY)
    def compFunctionalMatrices(self,element,k,alfa,c,ro,tInf):
        elemH=0
        elemC=0
#        elemCt0=0
        elemP=np.array([0.0 for x in range(self.lenN)])
        for i in range(self.lenPoints):
            H=0
            for v in self.globVar:
                H=np.add(np.matmul(np.transpose(np.array([element.points[i]['dNd'+v]])),np.array([element.points[i]['dNd'+v]])),H)
            elemH+=H*element.points[i]['detJ']*self.points[i]['w']
            C=self.points[i]['N^2']*element.points[i]['detJ']*self.points[i]['w']
#            Ct0=C*np.dot(element['t'],self.points[i]['N'])
            elemC+=C
#            elemCt0+=Ct0
        elemH*=k
        elemC*=c*ro
#        elemCt0*=c*ro
        surfaceH=0
        surfaceP=0
        if element.edge:
            for j in range(len(element.surface)):
                if element.surface[j]['edge']:
                    H=0
                    P=0
                    for point in self.surface[j]:
                        H+= point['N^2']*point['w']
                        P+= point['N']*point['w']
                    surfaceH+=H*element.surface[j]['len']
                    surfaceP+=P*element.surface[j]['len']
        surfaceH*=alfa*0.5
        surfaceP*=-alfa*tInf*0.5
        element.H=np.array(elemH+surfaceH)
        element.P=np.array(elemP+surfaceP)
        element.C=np.array(elemC)
#        element.Ct0=np.array(elemCt0)
    def compTempPoints(self,grid,dTau):
        HH=[[0 for y in range(grid.nn)] for x in range(grid.nn)]
        CC=[[0 for y in range(grid.nn)] for x in range(grid.nn)]
        PP=[0 for x in range(grid.nn)]
        for element in grid:
            C=element.C
            H=element.H
            P=element.P
            for i in range(self.lenN):
                PP[element.ids[i]]+=P[i]
                for j in range(self.lenN):
                    HH[element.ids[i]][element.ids[j]]+=H[i][j]
                    CC[element.ids[i]][element.ids[j]]+=C[i][j]
        A=np.array(HH)+np.array(CC)/dTau
        B=np.matmul(np.array(CC),np.transpose(np.array([grid['t']])))/dTau-np.transpose(np.array([PP]))
        return lg.solve(A,B)
    """
    def compTempPointsIntra(self,grid,dTau):
        HH=[[0 for y in range(grid.nn)] for x in range(grid.nn)]
        CC=[[0 for y in range(grid.nn)] for x in range(grid.nn)]
        CCt0=[[0 for y in range(grid.nn)] for x in range(grid.nn)]
        PP=[0 for x in range(grid.nn)]
        for element in grid:
            C=element.C
            Ct0=element.Ct0
            H=element.H
            P=element.P
            for i in range(self.lenN):
                PP[element.ids[i]]+=P[i]
                for j in range(self.lenN):
                    HH[element.ids[i]][element.ids[j]]+=H[i][j]
                    CC[element.ids[i]][element.ids[j]]+=C[i][j]
                    CCt0[element.ids[i]][element.ids[j]]+=Ct0[i][j]
        CCt0s=np.sum(CCt0,axis=0)
        A=np.array(HH)+np.array(CC)/dTau
        B=CCt0s/dTau-np.array(PP)
        return lg.solve(A,B)
    """
    def compGridTemp(self, grid, k, alfa, c, ro, tInf,dTau):
        gridLen=len(grid)
        if isinstance(k,abc.Sequence):
            k=iter(k)
        else:
            k=iter([k for i in range(gridLen)])
        if isinstance(alfa,abc.Sequence):
            alfa=iter(alfa)
        else:
            alfa=iter([alfa for i in range(gridLen)])
        if isinstance(c,abc.Sequence):
            c=iter(c)
        else:
            c=iter([c for i in range(gridLen)])
        if isinstance(ro,abc.Sequence):
            ro=iter(ro)
        else:
            ro=iter([ro for i in range(gridLen)])
        if isinstance(tInf,abc.Sequence):
            tInf=iter(tInf)
        else:
            tInf=iter([tInf for i in range(gridLen)])
        
        for element in grid:
            self.compElementPoints(element)
            self.compFunctionalMatrices(element,next(k),next(alfa),next(c),next(ro),next(tInf))
        t1=self.compTempPoints(grid,dTau)
        return (np.reshape(t1,len(t1)))


if __name__=='__main__':
    globalData=loadData('data.yml')
    print(globalData)
    g=Grid(globalData['B'],globalData['H'],globalData['nB'],globalData['nH'],globalData['t0'],globalData['edges'])
    print(g.ne)
    print(g.printNodeAttrs('edge'))
    input()
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

    x=Compute([n1,n2,n3,n4],dict([('X','Xsi'),('Y','Eta')]),[-0.7745966692414834,0.0,0.7745966692414834],[0.5555555555555556,0.8888888888888888,0.5555555555555556])
#    x=Compute([n1,n2,n3,n4],dict([('X','Xsi'),('Y','Eta')]),[-1.0/np.sqrt(3),1.0/np.sqrt(3)],[1.0,1.0])
    
    tau=0.0
    start=time.clock()
    while tau<globalData['tau']:
        t1=x.compGridTemp(g,k=globalData['k'],alfa=globalData['alfa'],c=globalData['c'],ro=globalData['ro'],tInf=globalData['tInf'],dTau=globalData['dTau'])
#        print(printSeq(t1,16),printSeq(t2,16),sep='\n')
        g.updateNodes(t1,'t')
        tau+=globalData['dTau']
        print(g.printNodeAttrs('t'))
    print("time: ",time.clock()-start)
    