# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 23:12:10 2017

@author: mpfun
"""

#%%
import json
import numpy as np
import matplotlib.pyplot as plt
import glob   
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


#%%
def read_json_file(j): 
    with open(j) as json_data:
        return json.load(json_data)
        
SHIFT = 25

#path = 'C:\\Users\\mpfun\\Documents\\Visual Studio 2017\\Projects\\gen\\Release\\experiment8\\';   

def read_experiment_dir(path, alg, fun, dim):
    path = path + '\\' + alg + '_' + fun + '_' + str(dim) + '_*.json'
    
    files=glob.glob(path)   
    
    data = {}
    for file in files:
        print(os.path.basename(file))
        with open(file) as json_data:
            data[os.path.basename(file).split('_')[0]] = json.load(json_data)
    
    return data

# plot_avg_errs([d2,d3,d5,d8,d12,d20], ['ppa', 'fwa', 'ppalevy'])
def plot_avg_errs(data, tx = [], emod=1.0):
    errs = {}
    means = {}
    ubs = {}
    lbs = {}
    for d in data:
        for k in d.keys():
            if k not in errs:
                errs[k] = []
                means[k] = []
                ubs[k] = []
                lbs[k] = []
                
            s = d[k]['stats']
            errs[k].append({'mean' : s['mean'], 'sd' : s['sd']})
            means[k].append(s['mean'] / emod)
            ubs[k].append((s['mean'] + s['sd']) / emod)
            lbs[k].append((s['mean'] + s['sd']) / emod)
    
    lgnd = []
    
    for k in means.keys():
        if len(tx) > 0:
            plt.plot(tx, means[k])
        else:
            plt.plot(means[k])
        #plt.plot(ubs[k], '.')
        #plt.plot(lbs[k], '.')
        lgnd.append(k)
        #lgnd.append(k + '-UB')
        #lgnd.append(k + '-LB')
    
    akey = [k for k in means.keys()]
    params = data[0][akey[0]]['params']
    #print(params)
    of = '{} - bounds: ({},{}) - fevals: {}'.format(
            params['objectiveFunction'], 
            params['initBounds[0]'][0], 
            params['initBounds[0]'][1],
            params['maxFevals'])
    plt.legend(lgnd)#[k for k in errs.keys()])
    plt.title('f(x) - ' + of)
    plt.xlabel('dimensions')
    plt.ylabel('mean f(x)')
    #if len(tx) > 0:
     #   plt.xticks(tx)
    plt.show()
    
    #return means

# plot_vs([d2['fwa'], d2['ppa'], d2['ppalevy']], 0, fromGen=2000)
def plot_vs(data, run, fromGen = -1, toGen=-1):
    
    for (idx, d) in enumerate(data):
        x = d['runs'][run]
        vals = x['bestMembersInGeneration']['values']
        if toGen != -1:
            vals = vals[:toGen]
        if fromGen != -1:
            vals = vals[fromGen:]
        plt.plot(vals)
   
    algs = [d['params']['algorithm'] for (idx, d) in enumerate(data)]
    of = data[0]['params']['objectiveFunction'] + ' d:' + str(data[0]['params']['dimensions'])
    of = of + ' run: ' + str(run)
    of = of + ' from: {} to: {}'.format(fromGen if fromGen != -1 else 0, toGen if toGen != -1 else 'inf')
    
    plt.legend(algs)
    plt.title(of)
    plt.ylabel('cost')
    plt.xlabel('gen')
    plt.show()
    
#%%

#d2['fwa']['runs'][0]['bestMembersInGeneration']['coordinates']
def draw_function_and_members(f, data, fname, aname, LB, UB, run=0, sh=1, azim=None, opt=None):
    
    members = data[aname]['runs'][run]['bestMembersInGeneration']['coordinates']
    
    xs = np.arange(LB, UB, 0.1)
    ys = np.arange(LB, UB, 0.1)
    
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    
    # Compute z to make the pringle surface.
    z = np.zeros([len(xs), len(ys)])
    for (ix, x) in enumerate(xs):
        for (iy, y) in enumerate(ys):
            z[ix][iy] = (f([x, y], sh))
            
    x, y = np.meshgrid(xs, ys)

    #ax.plot_surface(X, Y, zs, linewidth=0.3, antialiased=True, cmap=cm.jet)
    #ax.view_init(0, 0)
    #print(ax.azim)
    if azim is not None:
        ax.view_init(azim=azim)
    #ax.plot_wireframe(x, y, z, cmap=cm.jet)
    ax.plot_surface(x, y, z, linewidth=0.3, antialiased=True, cmap=cm.jet, alpha=0.6)
    cset = ax.contour(x, y, z, zdir='z', offset=725, cmap=cm.coolwarm)
    #cset = ax.contour(x, y, z, zdir='x', offset=40.2, cmap=cm.coolwarm)
    #cset = ax.contour(x, y, z, zdir='y', offset=40.2, cmap=cm.coolwarm)

    if opt is not None:
        ax.scatter([opt[0]], [opt[1]], [f(opt, sh)], s=50, c='#FF0000', marker='*')
    
    mxs1 = []
    mys1 = []
    mzs1 = []
    
    mxs2 = []
    mys2 = []
    mzs2 = []
    
    ll = 0
    for m in members:
        if m[0] > LB and m[0] < UB and m[1] > LB and m[1] < UB:
            ll = ll + 1
    print(ll)
    c = 0
    for m in members:
        if m[0] > LB and m[0] < UB and m[1] > LB and m[1] < UB:
            if False and c > (ll/100.0 * 10.0):
                mxs2.append(m[0])
                mys2.append(m[1])
                mzs2.append(f(m, sh))          
            else:
                mxs1.append(m[0])
                mys1.append(m[1])
                mzs1.append(f(m, sh))
 
            c = c + 1
    #print(sh)
    #print(mzs)
    
    print(len(mzs1))
    print(len(mzs2))
    print("{},{} -> {}", mxs1[-1], mys1[-1], mzs1[-1])
    ax.scatter(mxs1, mys1, mzs1, s=90, c = '#000000', marker='x')
    ax.scatter(mxs2, mys2, mzs2, s=60, c = '#000000', marker='x')
    
    ax.scatter([mxs1[-1]], [mys1[-1]], [[mzs1[-1]]], c='#000000', marker='^', s=90)
    


    plt.title(fname + ' - ' + aname + ' - ({}, {}) - run: {}'.format(LB, UB, run))
    plt.show()
    
    #return mzs1
#%%

def rosenbrock(member, shift = 0):
    sum = 0.0;

    o = shift;

    n = len(member)
    for i in range(0, n-1):
		
        x_1 = member[i + 1] - o + 1;
        x = member[i] - o + 1;

        inner1 = x_1 - x * x
        part1 = 100.0 * inner1 * inner1
        part2 = (x - 1.0) * (x - 1.0)

        sum += (part1 + part2)
    return sum

def schwefel(x, shift = 0):
    f = 0;
    for i in range(0, len(x)):
        a = 0.0;
        for j in range(0, i):
            a += x[j] + shift;
        f += a * a;
	
    return f

def schwefel7(member, shift=0):
    c = 418.98288727243369;

    s = sum([(x + shift) * np.sin(np.sqrt(abs(x + shift))) for x in member])

    return c * len(member) - s;
# quite nice rosenbrock 


#exp 8

#path='C:\\Users\\mpfun\\Documents\\Visual Studio 2017\\Projects\\gen\\Release\\experiment8\\'
#2 = read_experiment_dir(path, '*', 'rosenbrock', '2')
#draw_function_and_members(rosenbrock, d2, 'rosenbrock', 'fwa', 0, 70, sh=25, azim=-20, opt=[25, 25])



#%%

import copy

def exp_mod_ppa_1_read():
    path="C:\\Users\\mpfun\\Documents\\Visual Studio 2017\\Projects\\gen\\Release\\experiments_mod_ppa\\experiment1\\"
    dims=[2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 50, 100]
    algs=['ppa', 'ppalevy', 'fwa']
    funs=['rosenbrock', 'schwefel7', 'griewank', 'schwefel', 'ackley', 'michalewicz12']
    data = {}
    for f in funs:
        subdata = []
        for d in dims:
            subsubdata = {}
            for a in algs:
                print('read {} - {} - {}'.format(f, d, a))

                ssd = read_experiment_dir(path, a, f, d)[a]
                subsubdata[a] = ssd
                print(ssd)
            subdata.append(copy.deepcopy(subsubdata))
        data[f] = copy.deepcopy(subdata)
    return data

def exp_mod_ppa_2_read():
    path="C:\\Users\\mpfun\\Documents\\Visual Studio 2017\\Projects\\gen\\Release\\experiments_mod_ppa\\experiment2\\"
    dims=[2, 5, 10, 15, 30]# 50, 100]
    algs=['ppalevy', 'ppalevy-birds', 'ppa']
    funs=['rosenbrock', 'schwefel7', 'griewank', 'schwefel', 'ackley', 'michalewicz12']
    data = {}
    for f in funs:
        subdata = []
        for d in dims:
            subsubdata = {}
            for a in algs:
                print('read {} - {} - {}'.format(f, d, a))

                ssd = read_experiment_dir(path, a, f, d)[a]
                subsubdata[a] = ssd
                print(ssd)
            subdata.append(copy.deepcopy(subsubdata))
        data[f] = copy.deepcopy(subdata)
    return data

