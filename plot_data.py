#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 15:50:36 2019

@author: rplazevedo
"""

#import numpy as np
import matplotlib.pyplot as plt
import os
import configparser
import numpy as np

def run():
    config = configparser.ConfigParser()
    config.read('input.ini')
    name_def = config['Parameters']['name']
    Nx = int(config['PRS']['Lattice_size'])
    t0 = float(config['PRS']['initial_t'])
    w0 = float(config['PRS']['initial_width'])
    phi_0 = float(config['PRS']['initial_phi'])
    alpha = float(config['PRS']['alpha'])
    dt = float(config['PRS']['dt'])
    delta_x = float(config['PRS']['dx'])
    delta_y = float(config['PRS']['dy'])
    alpha_1 = float(config['Parameters']['alpha_1'])
    alpha_2 = float(config['Parameters']['alpha_2'])
    alpha_e = float(config['Parameters']['alpha_e'])
    Nw_alpha = float(config['Parameters']['Nw_alpha'])
    daeta = int(config['Parameters']['dlna/dlneta'])
    last_part = int(config['PRS']['sPRS_parts'])
    Nt = int(config['PRS']['points_per_part'])
    nparts = int(config['PRS']['mPRS_parts'])
    name = name_def       
    
    # name = str(input('Name of the data to plot?\n'))
    # if name == '':
    #     name = name_def
    
    
    
    lbl = ['sPRS', r'no vac repl.', r'vac repl.']
    
    xlim_val =(dt*(Nt-1)*(last_part),dt*(Nt-1)*(nparts+last_part))
    ylim_val =  0
    line_width_val = 1.0
    xscale_val = 'log'
    yscale_val = 'linear'
    title_val = fr"$w_0={w0:G}$, $\phi_0={phi_0:G}$, $N_x={Nx:G}$, $dt={dt:G}$, $\alpha_w={Nw_alpha:G}$, $\alpha_\rho={alpha_e:G}$"
    
    
    
    # V plots 
    plt.figure()
    number = 0
    while os.path.isfile(str(name) +'_'+str(number)+ '_v.data'):
        T, V = [], []
        for line in open(str(name) +'_'+str(number)+ '_v.data',"r"):
            x, y = line.split()
            T.append(float(x))
            V.append(float(y))
        plt.plot(T, V, label=lbl[number], lw=line_width_val)         
        if number == 0:
            V0 = V
        number += 1
    plt.title(title_val)
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'$\gamma^2 v^2$')
    plt.xscale(xscale_val)
    plt.yscale(yscale_val)
    plt.xlim(xlim_val)
    plt.ylim(bottom=ylim_val)
    plt.grid()
    plt.legend()
    plt.show()
    # Comparison
    plt.figure()
    number = 0
    V0a = np.array(V0)
    while os.path.isfile(str(name) +'_'+str(number)+ '_v.data'):
        T, V = [], []
        for line in open(str(name) +'_'+str(number)+ '_v.data',"r"):
            x, y = line.split()
            T.append(float(x))
            V.append(float(y))
        Va = np.array(V)
        plt.plot(T, np.abs(Va-V0a)/V0a, label=lbl[number], lw=line_width_val)
        number += 1
    plt.title(title_val)
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'$\gamma^2 v^2$ (frac. diff.)')
    plt.xscale(xscale_val)
    plt.yscale(yscale_val)
    plt.xlim(xlim_val)
    plt.ylim(bottom=ylim_val)
    plt.grid()
    plt.legend()
    plt.show()
    
    # Number of walls
    plt.figure()
    number = 0
    while os.path.isfile(str(name) +'_'+str(number)+ '_nwalls.data'):
        T, Nw = [], []
        for line in open(str(name) +'_'+str(number)+ '_nwalls.data',"r"):
            x, y = line.split()
            T.append(float(x))
            Nw.append(float(y))
        plt.plot(T, Nw, label=lbl[number], lw=line_width_val)         
        if number == 0:
            Nw0 = Nw
        number += 1
    plt.title(title_val)
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'$N_{walls}$')
    plt.xscale(xscale_val)
    plt.yscale(yscale_val)
    plt.xlim(xlim_val)
    plt.ylim(bottom=ylim_val)
    plt.grid()
    plt.legend()
    plt.show()
    # Comparison
    plt.figure()
    number = 0
    Nw0a = np.array(Nw0)
    while os.path.isfile(str(name) +'_'+str(number)+ '_nwalls.data'):
        T, Nw = [], []
        for line in open(str(name) +'_'+str(number)+ '_nwalls.data',"r"):
            x, y = line.split()
            T.append(float(x))
            Nw.append(float(y))
        Nwa = np.array(Nw)
        plt.plot(T, np.abs(Nwa-Nw0a)/Nw0a, label=lbl[number], lw=line_width_val)
        number += 1
    plt.title(title_val)
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'$N_{walls}$ (frac. diff.)')
    plt.xscale(xscale_val)
    plt.yscale(yscale_val)
    plt.xlim(xlim_val)
    plt.ylim(bottom=ylim_val)
    plt.grid()
    plt.legend()
    plt.show()
    
    # Vacuum plots
    plt.figure()
    number = 0
    while os.path.isfile(str(name) +'_'+str(number)+ '_exc_pts.data'):
        T, P = [], []
        for line in open(str(name) +'_'+str(number)+ '_exc_pts.data',"r"):
            x, y = line.split()
            T.append(float(x))
            P.append(float(y))
            
        plt.plot(T, P, label=lbl[number], lw=line_width_val)
        number += 1
    plt.title(title_val)
    plt.xlabel(r'$\eta$')
    plt.ylabel('Frac. of vacua')
    plt.xscale(xscale_val)
    plt.yscale(yscale_val)
    plt.xlim(xlim_val)
    plt.ylim(bottom=ylim_val)
    plt.grid()
    plt.legend()
    plt.show()

# Run the function
# run()