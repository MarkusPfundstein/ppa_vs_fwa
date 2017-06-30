# -*- coding: utf-8 -*-
"""
Created on Thu May  4 17:16:05 2017

@author: mpfun
"""

#%%
import numpy as np

def R():
    return 1.0

def dx(ni):
    return 2.0 * (1.0 - ni) * (R() - 0.5);

def update(xj, dx, LB, UB):
    return xj + (UB - LB) * dx;

xs = np.linspace(0, 1, 1000)

UB=100
LB=-100

xmapped = []
for x in xs:
    xmapped.append(update(1.0001, dx(x), UB, LB))

xmapped.reverse()
