# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:26:28 2017

@author: Piotrek
"""
import yaml
import numpy as np
import numpy.linalg as lg

def loadData(fileName):
    data={}
    with open(fileName) as file:
        data=yaml.load(file)
    return data

class ShapeFunc:
    def __init__(self,f,dFd):
        self.f=f
        self.dFd=dFd
        self.varNumber=len(dFd)
        self.var=dFd.keys()
    def __call__(self,*args):
        if len(args)!=self.varNumber:
            raise TypeError('Wrong number of variables')
        return self.f(*args)
    def __getitem__(self,index):
        return self.dFd[index]


class Node:
    def __init__(self,x,y,t=0,edge=False):
        self.x=x
        self.y=y
        self.t=t
        self.edge=edge
    def __repr__(self):
        return "({0},{1})".format(self.x,self.y)
    def __str__(self):
        return "({0},{1},t={2},edge={3})".format(self.x,self.y,self.t,self.edge)

class Element:
    def __init__(self,nodes,ids):
        self.nodes=nodes
        self.ids=ids
        self.surface=[all([nodes[(i+3)%4].edge,nodes[(i+4)%4].edge]) for i in range(4)]
        self.dXd=dict()
        self.dYd=dict()
        self.recipDetJ=[]
        self.recipJ=[]
        self.dNd={'x':[],'y':[]}
    def __repr__(self):
        rep="{3} {2}\n{0} {1}\n{4!r}\n".format(*self.ids,self.surface)
        return rep
    def __str__(self):
        strg="{3} {2}\n{0} {1}\n{4!r}\n".format(*self.nodes,self.surface)
        return strg
    def __getitem__(self,index):
        return self.nodes[index]
    def __len__(self):
        return len(self.nodes)
    def getXs(self):
        xs=[]
        for n in self.nodes:
            xs.append(n.x)
        return np.array(xs)
    def getYs(self):
        ys=[]
        for n in self.nodes:
            ys.append(n.y)
        return np.array(ys)


class Grid:
    def __init__(self,nB,nH,db,dh,x=0,y=0,t=0):
        nodes=[]
        for i in range(nB):
            for j in range(nH):
                if (i==0 or i==nB-1 or j==0 or j==nH-1):
                    nodes.append(Node(x+db*i,y+dh*j,t,True))
                else:
                    nodes.append(Node(x+db*i,y+dh*j,t))
        self.nodes=nodes
        elements=[]
        for i in range(nB-1):
            for j in range(nH-1):
                elements.append(Element([nodes[i*nH+j],nodes[(i+1)*nH+j],nodes[(i+1)*nH+j+1],nodes[i*nH+j+1]],[i*nH+j,(i+1)*nH+j,(i+1)*nH+j+1,i*nH+j+1]))
        self.elements=elements
    def __getitem__(self,index):
        return self.elements[index]
    def __len__(self):
        return len(self.elements)
    def __call__(self,index):
        return self.nodes[index]

class MesObject:
    def __init__(self,B,H,nB,nH):
        self.B=B
        self.H=H
        self.nH=nH
        self.nB=nB
        self.ne=(nH-1)*(nB-1)
        self.nn=nH*nB
    def generateGrid(self):
        self.grid=Grid(self.nB,self.nH,self.B/(self.nB-1),self.H/(self.nH-1))
    def printGrid(self):
        print(self.grid)


class Compute:
    def __init__(self,N,var,gaussQ,gaussW):
        self.N=N
        self.var=var
        self.points=[((gaussQ[i],gaussQ[j]),gaussW[i]*gaussW[j]) for i in range(len(gaussQ)) for j in range(len(gaussQ))]
        listdNd=dict()
        dNd=dict()
        for name in self.var:
            listdNd[name]=[list() for _ in range(len(N))]
            for i in range(len(N)):
                for point in self.points:
                    listdNd[name][i].append(self.N[i][name](*point[0])*point[1])
            dNd[name]=list(zip(*listdNd[name]))
            for i in range(len(dNd[name])):
                dNd[name][i]=np.array(dNd[name][i])
        self.dNd=dNd
    def compElement(self,element):
        for name in self.var:
            element.dXd[name]=[]
            element.dYd[name]=[]
            for d in self.dNd[name]:
                element.dXd[name].append(np.dot(d,element.getXs()))
                element.dYd[name].append(np.dot(d,element.getYs()))
        for name in self.var:
            for i in range(len(self.points)):
                jacobi=np.array([[element.dXd['xsi'][i],element.dYd['xsi'][i]],[element.dXd['eta'][i], element.dYd['eta'][i]]])
                recipJacobi=np.array([[jacobi[1,1],-jacobi[0,1]],[-jacobi[1,0],jacobi[0,0]]])
                element.recipDetJ.append(1.0/lg.det(jacobi))
                element.recipJ.append(recipJacobi)

        for i in range(len(self.points)):
            element.dNd['x'].append([])
            element.dNd['y'].append([])
            for j in range(len(self.N)):
                res=element.recipDetJ[i]*np.dot(element.recipJ[i],np.array([self.dNd['xsi'][i][j],self.dNd['eta'][i][j]]))
                element.dNd['x'][i].append(res[0])
                element.dNd['y'][i].append(res[1])


if __name__=='__main__':
    globalData=loadData('data.txt')
    print(globalData)
    mO=MesObject(globalData['B'],globalData['H'],globalData['nB'],globalData['nH'])
    mO.generateGrid()
    print(mO.grid(0),mO.grid(4),mO.grid(20),mO.grid(24),'\n')
    print(mO.grid[0],mO.grid[3],mO.grid[12],mO.grid[15],sep='')

    n1=ShapeFunc(lambda xsi,eta:0.25*(1-xsi)*(1-eta),
                 {'xsi':lambda xsi,eta:-0.25*(1-eta),'eta':lambda xsi,eta:-0.25*(1-xsi)})
    n2=ShapeFunc(lambda xsi,eta:0.25*(1+xsi)*(1-eta),
                 {'xsi':lambda xsi,eta:0.25*(1-eta),'eta':lambda xsi,eta:-0.25*(1+xsi)})
    n3=ShapeFunc(lambda xsi,eta:0.25*(1+xsi)*(1+eta),
                 {'xsi':lambda xsi,eta:0.25*(1+eta),'eta':lambda xsi,eta:0.25*(1+xsi)})
    n4=ShapeFunc(lambda xsi,eta:0.25*(1-xsi)*(1+eta),
                 {'xsi':lambda xsi,eta:-0.25*(1+eta),'eta':lambda xsi,eta:0.25*(1-xsi)})

    print(n1['xsi'](1.0,0.0))

    x=Compute([n1,n2,n3,n4],{'xsi','eta'},[-np.sqrt(3)/3,np.sqrt(3)/3],[1.0,1.0])
    x.compElement(mO.grid[0])
    print(mO.grid[0].dNd)
