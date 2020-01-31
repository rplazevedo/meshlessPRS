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

DTYPE = np.float32

def run():
    # Read paratmeters from the config file:
    config = configparser.ConfigParser()
    config.read('input.ini')
    
    Nx = int(config['PRS']['Lattice_size'])
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
    # all_phi = config['Parameters']['all_phi']
    last_part = int(config['PRS']['sPRS_parts'])
    Nt = int(config['PRS']['points_per_part'])
    nparts = int(config['PRS']['mPRS_parts'])
    
    # Set initial counters and parameters:
    V0 = (np.pi**2)*(phi_0)**2/(2*w0**2)
    count_op = 0
    
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
    
    def get_neighbours(coord, N):
        """
        Parameters
        ----------
        coord : tuple
            Indices (i,j) of a given point in the mesh.
        N : size of the mesh
    
        Returns
        -------
        Coordinates of the neighbouring points, given periodic boundary conditions.
        Output is a tuple: (i+1, i-1, j+1, j-1)
        """
        if x != 0:
            x_l = x-1
        else:
            x_l = N
        if x != N:
            x_r = x+1
        else:
            x_r = 0            
        if y != 0:
            y_d = y-1
        else:
            y_d = N
        if y != N:
            y_u = y+1
        else:
            y_u = 0
        return (x_r, x_l, y_u, y_d)
    
    def neighbours_array(cond, N):
        """
        Parameters
        ----------
        cond : tuple of arrays
            Tuple with the coordinates of the points that satisfy the condition.
        N : integer
            Size of the grid.
        Returns
        -------
        cond_l : array
            X coordinates of the points to the left.
        cond_r : array
            X coordinates of the points to the right.
        cond_d : array
            Y coordinates of the points down.
        cond_u : array
            Y coordinates of the points up.
        """
        cond_l = cond[0]-1
        cond_r = cond[0]+1
        cond_d = cond[1]-1
        cond_u = cond[1]+1
        
        cond_l[cond_l == -1] = N-1
        cond_r[cond_r == N] = 0
        cond_d[cond_d == -1] = N-1
        cond_u[cond_u == N] = 0
        
        return cond_l, cond_r, cond_d, cond_u   
    
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
        x_l, x_r, y_d, y_u = neighbours_array((x,y),Nx)
        for p in range(len(x)):
            # dic[(x[p], y[p])] = [phi[x[p],y[p]], d_phi[x[p],y[p]],
            #                      np.roll(phi,1,axis=0)[x[p],y[p]], 
            #                      np.roll(phi,-1,axis=0)[x[p],y[p]],
            #                      np.roll(phi,1,axis=1)[x[p],y[p]], 
            #                      np.roll(phi,-1,axis=1)[x[p],y[p]]]            
            dic[(x[p], y[p])] = [phi[x[p],y[p]], d_phi[x[p],y[p]],
                                 phi[x_r[p],y[p]], phi[x_l[p],y[p]], 
                                 phi[x[p],y_u[p]], phi[x[p],y_d[p]]]
            # print(len(dic))
        return dic    
    
    v = np.zeros((Nt),dtype=DTYPE)
    N_walls = np.zeros((Nt),dtype=DTYPE)
    n_op = np.zeros((Nt),dtype=DTYPE)
    tspan = np.zeros((Nt),dtype=DTYPE)

    # Load final data from initial PRS calculation:
    phi_init = np.load(str(name)+'_phi_data'+str(last_part)+'.npy')[-1]
    d_phi_init = np.load(str(name)+'_d_phi_data'+str(last_part)+'.npy')[-1]
    tspan_init = np.load(str(name)+'_tspan_data'+str(last_part)+'.npy')[-1]
    v[0] = np.load(str(name)+'_vdata'+str(last_part)+'.npy')[-1]
    N_walls[0] = np.load(str(name)+'_nwalls_data'+str(last_part)+'.npy')[-1]
    n_op[0] = np.load(str(name)+'_nop_data'+str(last_part)+'.npy')[-1]   
    count_op = n_op[0]
    cond_pts = np.zeros((Nt),dtype=DTYPE)
    tspan[0] = tspan_init
    for i in range(1,Nt):
        tspan[i] = tspan[0]+dt*i

    
    # Steps 1, 2 & 3: Initialize the diccionary for storing the values phi_dic[(i,j)].
    # (i,j) is a tuple, where i and j are the coordinates of the points.
    # The value stored is a list where the indices are as following:
    # 0:φ_{i,j}, 1:φ_dot_{i,j}, 2:φ_{i+1,j}, 3:φ_{i−1,j}, 4:φ_{i,j+1}, 5:φ_{i,j−1}
    phi_dic = array_to_dic(phi_init, d_phi_init, V(phi_init))
    
    for part in range(1,nparts+1):
        for n in range(0,Nt-1):
            # Step 4: Calculate the next time step for the values in the diccionary
            start_prs = time.time()   #timer for counting the time it takes to complete a prs step      
            non_vac = int(len(phi_dic))
            print("Vector size is:", non_vac)   #prints the number of non_vacuum elements 
            print("Starting PRS algorithm")     
            Nwn1 = 0
            vn1 = 0
            for point in phi_dic.values():
                # print((point[4]+point[5]-2*point[0])/(delta_y**2))
                nabla_phi = ( (point[2]+point[3]-2*point[0])/(delta_x**2)+
                                      (point[4]+point[5]-2*point[0])/(delta_y**2))  
                delta = 0.5*alpha*dt*(daeta)/tspan[n]
                d_V = V0*4*(point[0]/phi_0)*( (point[0]**2)/(phi_0**2)-1 )
                point[1] = ((1-delta)*point[1]+dt*(nabla_phi-d_V))/(1+delta)
                point[0] += dt * point[1]
                if abs(point[0]) <= Nw_alpha:
                    vn1 += (point[1]**2/(1.0))/(V0*( (point[0]**2)/(phi_0**2)-1)**2)       
                    Nwn1 += 1 
                count_op += 1
            cond_pts[n+1] = 1-len(phi_dic)/Nx**2
            N_walls[n+1] = Nwn1
            v[n+1] = vn1/(2*Nwn1)    
            n_op[n+1] = count_op
            print("PRS algorithm is finished %s seconds ---" % (time.time() - start_prs)) 
            
            # Step 5: Now we must update the neighbours
            start_neig = time.time()   #timer for counting the time it takes to complete a neighbours update step  
            print("Updating neighbors")
            for coord, point in phi_dic.items():
                # Get coordinates of the neighbours
                x,y = coord
                x_r, x_l, y_u, y_d = get_neighbours(coord, Nx)
                
                phi_sign = np.sign(point[0])
                try:
                    point[2] = phi_dic[(x_r, y)][0]
                except KeyError:
                    point[2] = phi_sign
                try:
                    point[3] = phi_dic[(x_l, y)][0]
                except KeyError:
                    point[3] = phi_sign
                try:
                    point[4] = phi_dic[(x, y_u)][0]
                except KeyError:
                    point[4] = phi_sign
                try:
                    point[5] = phi_dic[(x, y_d)][0]
                except KeyError:
                    point[5] = phi_sign
                    
            print("Neighbors update is finished %s seconds ---" % (time.time() - start_neig))
            
            # Step 6: update the list of non-vacuums
            
            # We can't change the size of a dictionary during an iteration, so we
            # store the coordinates of the values to delete in a list and the points
            # to add in an auxiliary dictionary
            
            to_remove = []
            new_points = {}        
            
            for coord, point in phi_dic.items():
                # Get coordinates of the neighbours
                x,y = coord
                x_r, x_l, y_u, y_d = get_neighbours(coord, Nx)
        
                # First check if the point still belongs in the dictionary and remove
                # it if it is a vacuum
                    
                if not energy_cond(point, V(point[0])):
                    to_remove.append(coord)
                
                # Now check the 4 points beside it (skipping over them if already in
                # the dictionary) and add them to the list if they are not a vacuum    
                
                if (x_r, y) not in phi_dic:
                    # Get (i+1, j), and check condition at the end
                    if x_r != Nx:
                        x_rr = x_r+1
                    else:
                        x_rr = 0
                    phi = point[2]
                    phi_l = point[0]
                    sign_phi = np.sign(phi)
                    try:
                        phi_r = phi_dic[(x_rr, y)][0]
                    except KeyError:
                        phi_r = sign_phi            
                    try:
                        phi_u = phi_dic[(x_r, y_u)][0]
                    except KeyError:
                        phi_u = sign_phi
                    try:
                        phi_d = phi_dic[(x_r, y_d)][0]
                    except KeyError:
                        phi_d = sign_phi
                        
                    if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_d], V(phi)):
                        new_points[x_r, y] = [phi, 0, phi_r, phi_l, phi_u, phi_d]
                
                if (x_l, y) not in phi_dic:
                    # Get (i-1, j), and check condition at the end
                    if x_l != 0:
                        x_ll = x_r-1
                    else:
                        x_ll = Nx
                    phi = point[3]
                    phi_r = point[0]
                    sign_phi = np.sign(phi)
                    try:
                        phi_l = phi_dic[(x_ll, y)][0]
                    except KeyError:
                        phi_l = sign_phi           
                    try:
                        phi_u = phi_dic[(x_l, y_u)][0]
                    except KeyError:
                        phi_u = sign_phi
                    try:
                        phi_d = phi_dic[(x_l, y_d)][0]
                    except KeyError:
                        phi_d = sign_phi
                        
                    if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_d], V(phi)):
                        new_points[x_l, y] = [phi, 0, phi_r, phi_l, phi_u, phi_d]
        
                if (x, y_u) not in phi_dic:
                    # Get (i, j+1), and check condition at the end
                    if y_u != Nx:
                        y_uu = y_u+1
                    else:
                        y_uu = 0
                    phi = point[4]
                    phi_d = point[0]
                    sign_phi = np.sign(phi)
                    try:
                        phi_u = phi_dic[(x, y_uu)][0]
                    except KeyError:
                        phi_u = sign_phi            
                    try:
                        phi_r = phi_dic[(x_r, y_u)][0]
                    except KeyError:
                        phi_r = sign_phi
                    try:
                        phi_l = phi_dic[(x_l, y_u)][0]
                    except KeyError:
                        phi_l = sign_phi
                        
                    if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_d], V(phi)):
                        new_points[x, y_u] = [phi, 0, phi_r, phi_l, phi_u, phi_d]            
                    
                if (x, y_d) not in phi_dic:
                    # Get (i, j-1), and check condition at the end
                    if y_d != 0:
                        y_dd = y_d-1
                    else:
                        y_dd = Nx
                    phi = point[5]
                    phi_u = point[0]
                    sign_phi = np.sign(phi)
                    try:
                        phi_d = phi_dic[(x, y_dd)][0]
                    except KeyError:
                        phi_d = sign_phi            
                    try:
                        phi_r = phi_dic[(x_r, y_d)][0]
                    except KeyError:
                        phi_r = sign_phi
                    try:
                        phi_l = phi_dic[(x_l, y_d)][0]
                    except KeyError:
                        phi_l = sign_phi
                        
                    if energy_cond([phi, 0, phi_r, phi_l, phi_u, phi_d], V(phi)):
                        new_points[x, y_d] = [phi, 0, phi_r, phi_l, phi_u, phi_d]
            # Now we delete vacuums        
            for pt in to_remove:
                del phi_dic[pt]
            # And update the dictionary with the new points
            phi_dic.update(new_points)
            
        # save part to disk
        np.save(str(name) + '_tspan_data' + str(part+int(last_part)) + '.npy', tspan)
        np.save(str(name) + '_vdata' + str(part+int(last_part)) + '.npy', v)
        np.save(str(name) + '_nwalls_data' + str(part+int(last_part)) + '.npy', N_walls)
        np.save(str(name) + '_nop_data' + str(part+int(last_part)) + '.npy', n_op)
        np.save(str(name) + '_exc_pts_data' + str(part+int(last_part)) + '.npy', cond_pts)       
        if (part < nparts):
            #Lets initialize the t again from the last point
            tspan_init = tspan[-1] 
            tspan[0] = tspan_init
            for i in range(1,Nt):
                tspan[i] = tspan_init+dt*i
            v[0] = v[-1]
            v[1:] = 0
        
            N_walls[0] = N_walls[-1]
            N_walls[1:] = 0
        
            n_op[0] = n_op[-1]
            n_op[1:] = 0
            
            cond_pts[0] = cond_pts[-1]

    print("Seconds:" +  str(time.time() - start_time))
    sys.stdout = orig_stdout
    f.close()