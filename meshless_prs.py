#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:55:16 2020

@author: rplazevedo, fferreira
"""

import time
import numpy as np
import sys
import configparser

# Read paratmeters from the config file:
config = configparser.ConfigParser()
config.read('input.ini')

Nx = int(config['PRS']['Lattice_size'])
t0 = float(config['PRS']['initial_t'])
w0 = float(config['PRS']['initial_width'])
phi_0 = float(config['PRS']['initial_phi'])
alpha = float(config['PRS']['alpha'])
dt = float(config['PRS']['dt'])
delta_x = float(config['PRS']['dx'])
delta_y = float(config['PRS']['dy'])
Nw_alpha = float(config['Parameters']['Nw_alpha'])
alpha_e = float(config['Parameters']['alpha_e'])
daeta = int(config['Parameters']['dlna/dlneta'])
name = config['Parameters']['name']
all_phi = config['Parameters']['all_phi']
last_part = int(config['PRS']['sPRS_parts'])
Nt = int(config['PRS']['points_per_part'])
nparts = int(config['PRS']['mPRS_parts'])


# Set initial counters and parameters:
V0 = (np.pi**2)*(phi_0)**2/(2*w0**2)
count_op = 0
Nw = 0

start_time = time.time()    # initial time for measuring speed of calculation

orig_stdout = sys.stdout
f = open('out.txt', 'a')
sys.stdout = f

# ----- FUNCTIONS -----

def V(phi):
    """
    Returns the potential of phi.
    """
    return V0*((phi/phi_0)**2-1)**2 

def array_to_dic(phi, d_phi, V):
    """
    Reads the initial phi and d_phi arrays and applies the energy condition.
    Returns a diccionary for storing the values: phi_dic[(i,j)].
    (i,j,s) is a tuple, where i and j are the coordinates of the points.
    The value stored is a list where the indices are as following:
    0:φ_{i,j}, 1:φ_dot_{i,j}, 2:φ_{i+1,j}, 3:φ_{i−1,j}, 4:φ_{i,j+1}, 5:φ_{i,j−1}
    """
    x, y = (np.abs(d_phi**2 
                    +(0.5*(np.roll(phi,1,axis=0)-np.roll(phi,-1,axis=0))/delta_x)**2
                    +(0.5*(np.roll(phi,1,axis=1)-np.roll(phi,-1,axis=1))/delta_y)**2
                    +V) >= V0*alpha_e).nonzero()
    dic = {}
    for p in range(len(x)):
        dic[(x[p], y[p])] = [phi[x[p]], d_phi[x[p]],
                             np.roll(phi,1,axis=0)[x[p],y[p]], 
                             np.roll(phi,-1,axis=0)[x[p],y[p]],
                             np.roll(phi,1,axis=1)[x[p],y[p]], 
                             np.roll(phi,-1,axis=1)[x[p],y[p]]]
    return dic

def energy_cond(phi_list, V):
    """
    Parameters
    ----------
    phi_list : list
        List with the values of phi, d_phi, etc.
    V : float
        Value of the potential at phi.

    Returns
    -------
    True if the point is not a vacuum, False if it is a vacuum, according to
    the energy condition.
    """
    return np.abs(phi_list[1]**2 
                    +(0.5*(phi_list[2]-phi_list[3])/delta_x)**2
                    +(0.5*(phi_list[4]-phi_list[5])/delta_y)**2
                    +V) >= V0*alpha_e

# Load final data from initial PRS calculation:
phi_init = np.load(str(name)+'_phi_data'+str(last_part)+'.npy')[-1]
d_phi_init = np.load(str(name)+'_d_phi_data'+str(last_part)+'.npy')[-1]
tspan_init = np.load(str(name)+'_tspan_data'+str(last_part)+'.npy')[-1]
v_init = np.load(str(name)+'_vdata'+str(last_part)+'.npy')[-1]
N_walls_init = np.load(str(name)+'_nwalls_data'+str(last_part)+'.npy')[-1]
n_op_init = np.load(str(name)+'_nop_data'+str(last_part)+'.npy')[-1]

tspan = np.arange(tspan_init, Nt*dt+tspan_init, dt)

# Steps 1, 2 & 3: Initialize the diccionary for storing the values phi_dic[(i,j)].
# (i,j) is a tuple, where i and j are the coordinates of the points.
# The value stored is a list where the indices are as following:
# 0:φ_{i,j}, 1:φ_dot_{i,j}, 2:φ_{i+1,j}, 3:φ_{i−1,j}, 4:φ_{i,j+1}, 5:φ_{i,j−1}

phi_dic = array_to_dic(phi_init, d_phi_init, V(phi_init))


for t in tspan:
    # Step 4: Calculate the next time step for the values in the diccionary
    start_prs = time.time()   #timer for counting the time it takes to complete a prs step      
    non_vac = int(len(phi_dic))
    print("Vector size is:", non_vac)   #prints the number of non_vacuum elements 
    print("Starting PRS algorithm")     
    Nw = 0
    for point in phi_dic.values():
        nabla_phi = ( (point[2]+point[3]-2*point[0])/(delta_x**2)+
                              (point[4]+point[5]-2*point[0])/(delta_y**2))                          
        delta = 0.5*alpha*dt*(daeta)/t
        d_V = V0*4*(point[0]/phi_0)*( (point[0]**2)/(phi_0**2)-1 )
        point[1] = ((1-delta)*point[1]+dt*(nabla_phi-d_V))/(1+delta)
        point[0] += dt * point[1]
        if point[0] <= Nw_alpha:
            v = (1.0/(2.0))*((point[1])**2/(1.0))/(V0*( (point[0]**2)/(phi_0**2)-1)**2)       
        count_op += 1
    print("PRS algorithm is finished %s seconds ---" % (time.time() - start_prs)) 
    
    # Step 5: Now we must update the neighbours
    start_neig = time.time()   #timer for counting the time it takes to complete a neighbours update step  
    print("Updating neighbors")    
    for coord, point in phi_dic.items():
        # Get coordinates of the neighbours
        x,y = coord
        if x != 0:
            x_l = x-1
        else:
            x_l = Nx
        if x != Nx:
            x_r = x+1
        else:
            x_r = 0            
        if y != 0:
            y_d = y-1
        else:
            y_d = Nx
        if y != Nx:
            y_u = y+1
        else:
            y_u = 0
            
        try:
            point[2] = phi_dic[(x_r, y)][0]
        except KeyError:
            point[2] = np.sign[phi_dic[x,y]]
        try:
            point[3] = phi_dic[(x_l, y)][0]
        except KeyError:
            point[3] = np.sign[phi_dic[x,y]]
        try:
            point[4] = phi_dic[(x, y_u)][0]
        except KeyError:
            point[4] = np.sign[phi_dic[x,y]]
        try:
            point[5] = phi_dic[(x, y_d)][0]
        except KeyError:
            point[5] = np.sign[phi_dic[x,y]]
            
    print("Neighbors update is finished %s seconds ---" % (time.time() - start_neig))
    
    # Step 6: update the list of non-vacuums
    
    for coord, point in phi_dic.items():
        
        x,y = coord
        
        if x != 0:
            x_l = x-1
        else:
            x_l = Nx
        if x != Nx:
            x_r = x+1
        else:
            x_r = 0
            
        if y != 0:
            y_d = y-1
        else:
            y_d = Nx
        if y != Nx:
            y_u = y+1
        else:
            y_u = 0

        # First check if the point still belongs in the dictionary and remove
        # it if it is a vacuum
            
        if not energy_cond(point, V(point[0])):
            del phi_dic[coord]
        
        # Now check the 4 points beside it (skipping over them if already in
        # the dictionary) and add them to the list if they are not a vacuum    
        
        if (x_r, y) not in phi_dic:
            # Get the neighbouring points, apply condition at the end
            if x_r != Nx:
                x_rr = x_r+1
            else:
                x_rr = 0
            phi = point[2]
            phi_l = point[0]
            try:
                phi_r = phi_dic[(x_rr, y)][0]
            except KeyError:
                phi_r = np.sign(phi)            
            try:
                phi_u = phi_dic[(x_r, y_u)][0]
            except KeyError:
                phi_u = np.sign(phi)
            try:
                phi_d = phi_dic[(x_r, y_d)][0]
            except KeyError:
                phi_d = np.sign(phi)
            if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_l], V(phi)):
                phi_dic[x_r, y] = [phi, 0, phi_r, phi_l, phi_u, phi_l]
        
        if (x_l, y) not in phi_dic:
            # Get the neighbouring points
            if x_l != 0:
                x_ll = x_r-1
            else:
                x_ll = Nx
            phi = point[3]
            phi_r = point[0]
            try:
                phi_l = phi_dic[(x_ll, y)][0]
            except KeyError:
                phi_l = np.sign(phi)            
            try:
                phi_u = phi_dic[(x_l, y_u)][0]
            except KeyError:
                phi_u = np.sign(phi)
            try:
                phi_d = phi_dic[(x_l, y_d)][0]
            except KeyError:
                phi_d = np.sign(phi)
                
            if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_l], V(phi)):
                phi_dic[x_l, y] = [phi, 0, phi_r, phi_l, phi_u, phi_l]

        if (x, y_u) not in phi_dic:
            # Get the neighbouring points
            if y_u != Nx:
                y_uu = y_u+1
            else:
                y_uu = 0
            phi = point[4]
            phi_d = point[0]
            try:
                phi_u = phi_dic[(x, y_uu)][0]
            except KeyError:
                phi_u = np.sign(phi)            
            try:
                phi_r = phi_dic[(x_r, y_u)][0]
            except KeyError:
                phi_r = np.sign(phi)
            try:
                phi_l = phi_dic[(x_l, y_u)][0]
            except KeyError:
                phi_l = np.sign(phi)
                
            if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_l], V(phi)):
                phi_dic[x, y_u] = [phi, 0, phi_r, phi_l, phi_u, phi_l]            
            
        if (x, y_d) not in phi_dic:
            # Get the neighbouring points
            if y_d != 0:
                y_dd = y_d-1
            else:
                y_dd = Nx
            phi = point[5]
            phi_u = point[0]
            try:
                phi_d = phi_dic[(x, y_dd)][0]
            except KeyError:
                phi_d = np.sign(phi)            
            try:
                phi_r = phi_dic[(x_r, y_d)][0]
            except KeyError:
                phi_r = np.sign(phi)
            try:
                phi_l = phi_dic[(x_l, y_d)][0]
            except KeyError:
                phi_l = np.sign(phi)
                
            if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_l], V(phi)):
                phi_dic[x, y_d] = [phi, 0, phi_r, phi_l, phi_u, phi_l]
                

                        
                        
                    














