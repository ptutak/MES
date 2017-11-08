# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:26:28 2017

@author: Piotrek
"""
import yaml


def loadData(fileName):
    data={}
    with open(fileName) as file:
        data=yaml.load(file)
    return data

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
    def __repr__(self):
        rep="{3} {2}\n{0} {1}\n{4!r}\n".format(*self.ids,self.surface)
        return rep
    def __str__(self):
        strg="{3} {2}\n{0} {1}\n{4!r}\n".format(*self.nodes,self.surface)
        return strg

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
        return self.nodes[index]
    def __call__(self,index):
        return self.elements[index]
        
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


if __name__=='__main__':
    globalData=loadData('data.txt')
    print(globalData)
    mO=MesObject(globalData['B'],globalData['H'],globalData['nB'],globalData['nH'])
    mO.generateGrid()
    print(mO.grid[0],mO.grid[4],mO.grid[20],mO.grid[24],'\n')
    print(mO.grid(0),mO.grid(3),mO.grid(12),mO.grid(15),sep='')