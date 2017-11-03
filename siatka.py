# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 13:26:28 2017

@author: Piotrek
"""

H,B,nH,nB=[float(x) for x in input().strip().split()]


class Node:
    def __init__(self,x,y,t,edge=0):
        self.x=x
        self.y=y
        self.t=t
        self.edge=edge
    def __repr__(self):
        return "({0},{1})".format(self.x,self.y)

class Element:
    def __init__(self,nodes):
        id=[]
        x=nodes[0].x
        y=nodes[0].y
        id[0],id[1],id[2],id[3]=x*5+y+1,(x+1)*5+y+1,(x+1)*5+y+2,x*5+y+2
        self.id=id
        self.nodes=nodes
        self.x=x
        self.y=y
        self.surface=[all([nodes[i%4].edge,nodes[i%4+1].edge] for i in range(3,7))]
    def __repr__(self):
        repr="{3} {2}\n{0} {1}\n{4}".format(*self.id,self.surface)
        return repr
    
class Grid:
    def __init__(self,nH,nB,x=0,y=0):
        
class GlobalData:
    def __init__(self,B,H,nB,nH):
        self.B=B
        self.H=H
        self.nH=nH
        self.nB=nB
        self.ne=(nH-1)*(nB-1)
        self.nn=nH*nB
        




[Node(x,y,t,any(surface[0:2])),Node(x+1,y,t,any(surface[1:3])),Node(x+1,y+1,t,any(surface[2:4])),Node(x,y+1,t,any(surface[3]))]