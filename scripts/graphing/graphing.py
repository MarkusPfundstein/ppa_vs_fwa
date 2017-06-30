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

SHIFT = 0

#path = 'C:\\Users\\mpfun\\Documents\\Visual Studio 2017\\Projects\\gen\\Release\\experiment8\\';   

#%%
def draw_schwefel():
    f = schwefel7
    fname = 'Schwefel'
    UB=100
    LB=-100
    
    xs = np.arange(LB, UB, 0.1)    
    
    # Compute z to make the pringle surface.
    ys = np.zeros([len(xs)])
    for (ix, x) in enumerate(xs):
        ys[ix]= f([x])
        
    fig, ax = plt.subplots()
    ax.plot(xs, ys)


    bbox_props = dict(boxstyle="circle,pad=0.3", fc="white", ec="black", lw=2)
    ax.text(-100, 357, "1", ha="center", va="center", rotation=0,
            size=10,
            bbox=bbox_props)
    

    bbox_props = dict(boxstyle="circle,pad=0.3", fc="white", ec="black", lw=2)
    ax.text(-26, 388, "2", ha="center", va="center", rotation=0,
            size=10,
            bbox=bbox_props)
    
    bbox_props = dict(boxstyle="circle,pad=0.3", fc="white", ec="black", lw=2)
    ax.text(6, 408, "3", ha="center", va="center", rotation=0,
            size=10,
            bbox=bbox_props)
    
    bbox_props = dict(boxstyle="circle,pad=0.3", fc="white", ec="black", lw=2)
    ax.text(65, 365, "4", ha="center", va="center", rotation=0,
            size=10,
            bbox=bbox_props)

    plt.title('{} [{}, {}]'.format(fname, LB, UB))
    plt.show()
    


#%%

def draw_function(f, fname, LB, UB, sh=1, azim=None, opt=None):
    
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
    #ax.contour(x, y, z, zdir='z', offset=725, cmap=cm.coolwarm)
    if opt is not None:
        ax.scatter([opt[0]], [opt[1]], [f(opt, sh)], s=100, c='#FF0000', marker='*')

    plt.title('{} [{}, {}]'.format(fname, LB, UB))
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

    n = len(member)
    for i in range(0, n-1):
		
        x_1 = member[i + 1] + 1;
        x = member[i] + 1;

        inner1 = x_1 - x * x
        part1 = 100.0 * inner1 * inner1
        part2 = (x - 1.0) * (x - 1.0)

        sum += (part1 + part2)
    return sum

def ackleys_path(member, shift=0):
    a = 20
    b = 0.2
    c = 2 * np.pi
    n = len(member)
    o = shift
    s1 = 0.0
    s2 = 0.0
    for i in range(0, n):
        x = member[i]
        s1 += pow(x, 2)
        s2 += np.cos(c * x)

    e1 = np.exp(-b * np.sqrt(s1 / n))
    e2 = np.exp(s2 / n)
    
    return -a * e1 - e2 + a + np.exp(1);

def schwefel(x, shift = 0):
    f = 0;
    for i in range(0, len(x)):
        a = 0.0;
        for j in range(0, i):
            a += x[j];
        f += a * a;
	
    return f

def schwefel7(member, shift=0):
    c = 418.98289;

    s = sum([(x) * np.sin(np.sqrt(abs(x))) for x in member])

    return c * len(member) - s;
# quite nice rosenbrock 


#exp 8

#path='C:\\Users\\mpfun\\Documents\\Visual Studio 2017\\Projects\\gen\\Release\\experiment8\\'
#2 = read_experiment_dir(path, '*', 'rosenbrock', '2')
#draw_function_and_members(rosenbrock, d2, 'rosenbrock', 'fwa', 0, 70, sh=25, azim=-20, opt=[25, 25])



#%%

EXP_DIR = "/Users/markus/Documents/THESIS/experiment_paper"

import copy

def read_json_file(j): 
    with open(j) as json_data:
        return json.load(json_data)
    
def read_experiment_dir(path, alg, fun, dim, shift):
    path = path + '/' + alg + '_' + fun + '_' + str(dim) + '_' + str(shift) + '_*.json'
    
    files=glob.glob(path)   
    
    data = {}
    for file in files:
        #print(os.path.basename(file))
        with open(file) as json_data:
            data[os.path.basename(file).split('_')[0]] = json.load(json_data)
    
    return data

def exp_thesis_read(algs, funs, dims, shifts):
    path=EXP_DIR
    data = {}
    for f in funs:
        subdata = {}
        for d in dims:
            subsubdata = {}
            for a in algs:
                subsubsubdata = {}
                for s in shifts: 
                    #rint('read {} - {} - {}'.format(f, d, a))
    
                    ssd = read_experiment_dir(path, a, f, d, s)[a]
                    subsubsubdata[s] = ssd
                    #print(ssd)
                subsubdata[a] = copy.deepcopy(subsubsubdata) 
            subdata[d] = copy.deepcopy(subsubdata)
        data[f] = copy.deepcopy(subdata)
    return data


#%%


# plot_avg_error_for_function(x)
def plot_avg_error_for_function(dataset, function_name, algorithms, shifts, file_title, log=False):
    #print dataset[function_name]
    
    function_data = dataset[function_name]
    
    errs = {}
    means = {}
    ubs = {}
    lbs = {}
    
  #  algorithms = []

    for shift in shifts:  
        errs[shift] = {}
        means[shift] = {}
        ubs[shift] = {}
        lbs[shift] = {}
        
        for dimension in function_data.keys():
           # print('dimension {}'.format(dimension))
            for algorithm in algorithms:
#                if algorithm not in algorithms:
 #                   algorithms.append(algorithm)
                #print('algorithm {}'.format(algorithm))
                exp_data = function_data[dimension][algorithm][shift]
                #print(exp_data['stats'])
                
                if algorithm not in errs[shift]:
                    errs[shift][algorithm] = []
                    means[shift][algorithm] = []                    
                    ubs[shift][algorithm] = []
                    lbs[shift][algorithm] = []
                
                s = exp_data['stats']
                errs[shift][algorithm].append({'mean' : s['mean'], 'sd' : s['sd']})
                if log is True:
                    means[shift][algorithm].append(np.log(s['mean']))
                else:
                    means[shift][algorithm].append(s['mean'])
                    
                ubs[shift][algorithm].append((s['mean'] + s['sd']))
                lbs[shift][algorithm].append((s['mean'] + s['sd']))
    
    f = plt.figure()
    axes = f.add_subplot(111)
    
#    for algorithm in algorithms:
    lgn = []
    for algorithm in algorithms:
        for shift in shifts:
            axes.plot([i for i in range(0, len(function_data.keys()))], means[shift][algorithm])
#            plt.plot([d for d in function_data.keys()], means[shift][algorithm])
            lgn.append('{}+{}'.format(algorithm, shift))
    
    params = function_data[2][algorithms[0]][0]['params']
    title = '{} - mean error'.format(
            params['objectiveFunction'])
    
    a=axes.get_xticks().tolist()
    for i in range(1, len(function_data.keys()) + 1):
        a[i] = pow(2, i)
    axes.set_xticklabels(a)
    
    plt.title(title)
    plt.xlabel('dimensions')
    if log is True:
        plt.ylabel('log f(x)')
    else:
        plt.ylabel('f(x)')
    plt.legend(lgn)
    add = ''
    if log is True:
        add = '_LOG_'

    plt.savefig(EXP_DIR + '/FIGURE_' + file_title + add + function_name + '.png')
    
#%%
def latex_tables(dataset, function_name, algorithms, shifts, file_title, log=False):
    #print dataset[function_name]
    
    function_data = dataset[function_name]
    
    

    for shift in shifts:  
        for dimension in function_data.keys():
            for algorithm in algorithms:
#                if algorithm not in algorithms:
 #                   algorithms.append(algorithm)
                #print('algorithm {}'.format(algorithm))
                exp_data = function_data[dimension][algorithm][shift]
                #print(exp_data['stats'])
                
                if algorithm not in errs[shift]:
                    errs[shift][algorithm] = []
                    means[shift][algorithm] = []                    
                    ubs[shift][algorithm] = []
                    lbs[shift][algorithm] = []
                
                s = exp_data['stats']
                errs[shift][algorithm].append({'mean' : s['mean'], 'sd' : s['sd']})
                if log is True:
                    means[shift][algorithm].append(np.log(s['mean']))
                else:
                    means[shift][algorithm].append(s['mean'])
                    
                ubs[shift][algorithm].append((s['mean'] + s['sd']))
                lbs[shift][algorithm].append((s['mean'] + s['sd']))
#%%
def print_table(X, dim, shifts):
#    \emph{Ackley}$_{0}$ & $\mathbf{4.44\cdot 10^{-16} \pm 0}$
#                  & $20.99 \pm 0.06$
#                  & $20.02 \pm 0.01$\\
#    \emph{Ackley}$_{12}$ & $\mathbf{4.37 \pm 0.97}$
#                  & $20.97 \pm 0.07$
#                  & $20.02 \pm 0.00$\\
#    \emph{Ackley}$_{25}$ & $\mathbf{5.76 \pm 2.26}$
#                  & $20.98 \pm 0.07$
#                  & $20.02 \pm 0.01$\\\hline
    
    for function_name in X.keys():
        print('\emph\{' + function_name + '\} & ', end='')
        function_data = X[function_name][dim]
        for shift in shifts: 
            for algorithm in function_data.keys():
                run_data = function_data[algorithm][shift]
                
                mean = run_data['stats']['mean']
                sd = run_data['stats']['sd']
                print('${} \pm {}$ &'.format(mean, sd))

#%%
def print_all():
    DIMS=[2, 4, 8, 16, 32, 64, 128]# 50, 100]
    ALGS=['ppalevy', 'ppa', 'fwa']
    FUNS=['ackley', 
          'cigar', 
          'ellipse', 
          'griewank', 
          'michalewicz12', 
          'rastrigrin', 
          'rosenbrock', 
          'schwefel7', 
          'schwefelFWA', 
          'schwefel', 
          'sphere', 
          'tablet']
    SHIFTS=[0, 12, 25]
    ex
    X = exp_thesis_read(ALGS, FUNS, DIMS, SHIFTS)
    for function_name in FUNS:
        plot_avg_error_for_function(X, function_name, ALGS, SHIFTS, 'all')
        #plot_avg_error_for_function(X, function_name, ALGS, SHIFTS, 'all', log=True)
        
def print_shifts():
    DIMS=[2, 4, 8, 16, 32, 64, 128]# 50, 100]
    ALGS=['ppalevy', 'ppa', 'fwa']
    FUNS=['ackley', 
          'cigar', 
          'ellipse', 
          'griewank', 
          'michalewicz12', 
          'rastrigrin', 
          'rosenbrock', 
          'schwefel7', 
          'schwefelFWA', 
          'schwefel', 
          'sphere', 
          'tablet']
    SHIFTS=[0, 12, 25]
    
    X = exp_thesis_read(ALGS, FUNS, DIMS, SHIFTS)
    for function_name in FUNS:
        for shift in SHIFTS:
            plot_avg_error_for_function(X, function_name, ALGS, [shift], '_shift_' + str(shift))
        #plot_avg_error_for_function(X, function_name, ALGS, SHIFTS, 'all', log=True)

def print_tables():
    DIMS=[128]# 50, 100]
    ALGS=['ppalevy', 'ppa', 'fwa']
    FUNS=['ackley', 
          'cigar', 
          'ellipse', 
          'griewank', 
          'michalewicz12', 
          'rastrigrin', 
          'rosenbrock', 
          'schwefel7', 
          'schwefelFWA', 
          'schwefel', 
          'sphere', 
          'tablet']
    SHIFTS=[0, 12, 25]
    
    X = exp_thesis_read(ALGS, FUNS, DIMS, SHIFTS)
    
    return X