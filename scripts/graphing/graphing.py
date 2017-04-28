# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 23:12:10 2017

@author: mpfun
"""

#%%
import json
import numpy as np
import matplotlib.pyplot as plt

#%%
def read_json_file(j): 
    with open(j) as json_data:
        return json.load(json_data)
        
def plot_objective_values(data, run, alg_name, function_name, to=-1):
    
    objective_values = data[run]['bestMembersInGeneration']['values']
    if to != -1:
        objective_values = objective_values[:to]
    
    plt.plot(
            objective_values, 
            "g")
    
    plt.legend(['alg_name'])
    plt.title(function_name)
    plt.ylabel('f')
    plt.xlabel('generation')
    plt.show()
    
    
def plot_ppa_vs_fwa(ppa_data, fwa_data, run, function_name):
    
    ppa_objective_values = ppa_data[run]['bestMembersInGeneration']['values']
    fwa_objective_values = fwa_data[run]['bestMembersInGeneration']['values']
    
    plt.plot(
            ppa_objective_values, 
            "g")

    plt.plot(    
            fwa_objective_values, 
            "b")
    
    plt.legend(['ppa', 'fwa'])
    plt.title(function_name)
    plt.ylabel('f')
    plt.xlabel('generation')
    plt.show()
    