# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:26:28 2017

@author: Piotrek
"""



class Node:
    def __init__(self,x,y,t=0,edge=0):
        self.x=x
        self.y=y
        self.t=t
        self.edge=edge
    def __repr__(self):
        return "({0},{1})".format(self.x,self.y)

class Element:
    def __init__(self,nodes,ids):
        self.nodes=nodes
        self.ids=ids
        self.surface=[all([nodes[(i+3)%4].edge,nodes[(i+4)%4].edge] for i in range(4))]
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
        
class GlobalData:
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
    with open('data.txt') as f:
        for line in f:
            if f[0]=='#':
                continue
            else:
                inp=line.strip().split(',')
                B=float(inp[0].strip())
                H=float(inp[1].strip())
                nB=int(inp[2].strip())
                nH=int(inp[3].strip())
                break
    gD=GlobalData(B,H,nB,nH)
    gD.generateGrid()
    print(gD.grid[0],gD.grid[4],gD.grid[20],gD.grid[24])
    print(gD.grid(0),gD.grid(3),gD.grid(12),gD.grid(15),sep='')