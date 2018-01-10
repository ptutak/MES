#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 20:42:47 2018

@author: ptutak
"""

for x in range(20):
    if x<5:
        alfa=3.49+0.93*x
    else:
        alfa=2.32*x**0.25
    print(x,alfa)
    