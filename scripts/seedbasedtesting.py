# -*- coding: utf-8 -*-
"""
Created on Tue May  2 23:45:16 2017

@author: mpfun
"""

#%%
import numpy as np

def eq4(A, l, u):
    return A - (l / u)

def g(l, u):
    print("{}, {}".format(l, u))
    return l < (u)

ls = np.linspace(1.0, 1.5, 20)
us = np.linspace(0.01, 1, 20)
As = np.linspace(0, 10)

xs = []
for l in ls:
    for u in us:
        if g(l, u):
            for A in As:
                xs.append({"v" : eq4(A, l, u), "w" : [l,u,A]})
            


print(eq4(10, 1.1, 0.1))
print(max(xs, key=lambda x : x["v"]))