#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  10 14:10:52 2019

@authors: rplazevedo, fferreira
"""

import numpy as np
import sys
import time
import configparser
    
DTYPE = np.float32

def run():
  
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
    daeta = int(config['Parameters']['dlna/dlneta'])
    name = config['Parameters']['name']
    all_phi = config['Parameters']['all_phi']
    
    
    V0 = (np.pi**2)*(phi_0)**2/(2*w0**2)
    count_op = 0
    Nw = 0
    
    Nt = int(config['PRS']['points_per_part'])
    nparts = int(config['PRS']['sPRS_parts'])
    
    tspan= np.zeros((Nt),dtype=DTYPE)
    phi = np.zeros((Nt,Nx,Nx),dtype=DTYPE)
    nabla_phi = np.zeros((Nt,Nx,Nx),dtype=DTYPE)
    d_phi = np.zeros((Nt,Nx,Nx),dtype=DTYPE)
    d_V = np.zeros((Nt,Nx,Nx),dtype=DTYPE)
    #phi_w = np.zeros((Nt,Nx,Nx),dtype=DTYPE)
    v = np.zeros((Nt),dtype=DTYPE)
    N_walls = np.zeros((Nt),dtype=DTYPE)
    n_op = np.zeros((Nt),dtype=DTYPE)

    start_time = time.time()
    
    orig_stdout = sys.stdout
    f = open('out.txt', 'a')
    sys.stdout = f
    
    print ("Running with initial sPRS")
    print("GRID POINTS on space interval (Nx=Ny): " + str(Nx))
    print("GRID POINTS on time interval (Nt): " + str(Nt))
    print("w0 = ", w0, ", V0 = ", V0, ", phi_0 = ", phi_0, "t_0 =", t0, "dt = ", dt)
    
    
    tspan[0] = t0
    
    for i in range(1,Nt):
        tspan[i] = t0+dt*i
    
    
    #phi = np.memmap('phi.mymemmap', dtype='float32', mode='w+', shape=(Nt,Nx,Nx)) 
    #nabla_phi = np.memmap('nabla_phi.mymemmap', dtype='float32', mode='w+', shape=(Nt,Nx,Nx))
    #d_phi = np.memmap('d_phi.mymemmap', dtype='float32', mode='w+', shape=(Nt,Nx,Nx))
    #d_V =  np.memmap('dV.mymemmap', dtype='float32', mode='w+', shape=(Nt,Nx,Nx))
    #Ek = np.memmap('Ek.mymemmap', dtype='float32', mode='w+', shape=(Nt))
    
    
    
    #Condition: initial value of phi_ij lies between -1 and 1

    phi[0] = np.random.uniform(-1,1,(Nx,Nx))
    wall_loc = abs(phi[0]) <= Nw_alpha
    Nw = np.sum(wall_loc)
    v[0] = (1.0/(2.0*Nw))*np.sum(np.where(wall_loc,
                         ((d_phi[0])**2/(1.0))/(V0*( (phi[0]**2)/(phi_0**2)-1)**2),
                         0))
    
    #print(((d_phi[0])**2/(8*pi)))
    #print((V0*(phi[0]**2-1))**2)
        
    #phi[0,i,j] = t0
    #print("Phi", 0, i, j, phi[0,i,j], sep="  ")
    #print(phi)
    #phi[0,-1,:]=phi[0,0,:]
    #phi[0,:,-1]=phi[0,:,0]
    
    N_walls[0] = Nw
    n_op[0] = 1
    
    #Lets calculate the initial values for nabla(phi_ij)
    nabla_phi[0] = ( (np.roll(phi[0],1,axis=0)+np.roll(phi[0],-1,axis=0)-2*phi[0])/(delta_x**2)+
                         (np.roll(phi[0],1,axis=1)+np.roll(phi[0],-1,axis=1)-2*phi[0])/(delta_y**2))      
    #print("Nabla_Phi", 0, i, j, nabla_phi[0,i,j], sep="  ")   
    
    #phi[0,-1,:]=0
    #phi[0,:,-1]=0
    
    #Initial value of d_phi is zero
       
    #Computing the next d_phi
    delta = 0.5*alpha*dt*(daeta)/(tspan[0]) 
    d_V[0]=V0*4*(phi[0]/phi_0)*( (phi[0]**2)/(phi_0**2)-1 )                                 
    d_phi[1] = ((1-delta)*d_phi[0]+dt*(nabla_phi[0]-d_V[0]))/(1+delta)
         
    #Ek[1] = Ek[1] +  1/(8*pi)*(d_phi[1])**2
    
    #Computing the next phi_ij
    phi[1] = phi[0] + dt*d_phi[1]

    # abs(phi) < Nw_alpha -> Wall, else it's radiation
    # v is the average velocity of the walls!
    wall_loc = abs(phi[1]) <= Nw_alpha
    Nw = np.sum(wall_loc)
    v[1] = (1.0/(2.0*Nw))*np.sum(np.where(wall_loc,
                         ((d_phi[1])**2/(1.0))/(V0*( (phi[1]**2)/(phi_0**2)-1)**2),
                         0))
    N_walls[1] = Nw
    n_op[1] = n_op[0]+1
        
    print("Initialization is done--- %s seconds ---" % (time.time() - start_time))
    
    
    #Loop
    for part in range(1,nparts+1):
        
        if (part ==1):
            start=1
        else:
            start=0
    
        for n in range(start,Nt-1):
            #            step_start_time = time.time()                        
            # Calculates the gradient of phi   
            nabla_phi[n] = ( (np.roll(phi[n],1,axis=0)+np.roll(phi[n],-1,axis=0)-2*phi[n])/(delta_x**2)+
                             (np.roll(phi[n],1,axis=1)+np.roll(phi[n],-1,axis=1)-2*phi[n])/(delta_y**2))                          
 
            delta = 0.5*alpha*dt*(daeta)/(tspan[n])                                  
            d_V[n]=V0*4*(phi[n]/phi_0)*( (phi[n]**2)/(phi_0**2)-1 )
            d_phi[n+1] = ((1-delta)*d_phi[n]+dt*(nabla_phi[n]-d_V[n]))/(1+delta)
                  
            #Ek[n+1] = Ek[n+1] +  1/(8*pi)*(d_phi[n+1])**2
            phi[n+1] = phi[n] + dt*d_phi[n+1]

            # abs(phi) < Nw_alpha -> Wall, else it's radiation
            # v is the average velocity of the walls!
            wall_loc = abs(phi[n+1]) <= Nw_alpha
            Nw = np.sum(wall_loc)
            v[n+1] = (1.0/(2.0*Nw))*np.sum(np.where(wall_loc,
                                 ((d_phi[n+1])**2/(1.0))/(V0*( (phi[n+1]**2)/(phi_0**2)-1)**2),
                                 0))
            

#            else:
#                phi_w[n+1,i,j]=0


            count_op = count_op+1    
            N_walls[n+1] = Nw
        
#            print("n: " , n, "count: ", count_op, "part = ", part)
#            print("--- %s seconds ---" % (time.time() - step_start_time))
            n_op[n+1] = 2+count_op
    
            #print("n_op[n+1]", n_op[n+1])
    
    #    print("Part " + str(part+int(last_part))+" is done--- %s seconds ---" % (time.time() - start_time))
    
        #numpart = part-1
    
    #    print("Part " + str(part)+" is done--- %s seconds ---" % (time.time() - start_time))
    
    
        if (all_phi in ['yes', 'YES', 'init', 'INIT']):
    
            np.save(str(name) + '_phi_data' + str (part) +'.npy', phi)
            np.save(str(name) + '_d_phi_data' + str(part) + '.npy', d_phi)

        elif part == nparts:
            np.save(str(name) + '_phi_data' + str (part) +'.npy', phi)
            np.save(str(name) + '_d_phi_data' + str(part)+ '.npy', d_phi)
            
        
    
        #np.save('d_V_data' + str(part) + '.npy', d_V)
        np.save(str(name) + '_tspan_data' + str(part) + '.npy', tspan)
        #np.save('phiw_data' + str (part) +'.npy', phi_w)
        np.save(str(name) + '_vdata' + str(part) +'.npy', v)
        np.save(str(name) + '_nwalls_data' + str(part) + '.npy', N_walls)
        np.save(str(name) + '_nop_data' + str(part) + '.npy', n_op)
    
    #    print("Part " + str(part)+" is saved--- %s seconds ---" % (time.time() - start_time))
    
        #Lets initialize the t again from the last point
        tspan[0] = tspan[-1]
        tspan[1:] = 0
    
        for i in range(1,Nt):
            tspan[i] = tspan[0]+dt*i  
    
    
        d_V[:,:] = 0
    
        d_phi[0,:,:] =  d_phi[-1,:,:]
        d_phi[1:,:,:] =  0
    
        phi[0,:,:] =  phi[-1,:,:]
        phi[1:,:,:] =  0
    
        nabla_phi[:,:] = 0
        v[0] = v[-1]
        v[1:] = 0
    
        N_walls[0] = N_walls[-1]
        N_walls[1:] = 0
    
        n_op[0] = n_op[-1]
        n_op[1:] = 0
    
    #    print("Going to Part " + str(part+1)+"--- %s seconds ---" % (time.time() - start_time))
    
    #np.save(str(name) + '_Ek_data.npy', Ek)
    
    print('Count= ', count_op)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    
    sys.stdout = orig_stdout
    f.close()
