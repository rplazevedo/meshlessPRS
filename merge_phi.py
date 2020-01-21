import numpy as np
import random
import sys
import time
import os
import configparser


DTYPE = np.float32

config = configparser.ConfigParser()
config.read('input.ini')

def run():    
    Nx = int(config['PRS']['Lattice_size'])
    Nt_int = int(config['PRS']['points_per_part'])
    parts = int(config['PRS']['sPRS_parts'])+int(config['PRS']['mPRS_parts'])
    name = config['Parameters']['name']
    
    #parts = 150
    #Nt_int = 10
    
    Nt=parts*Nt_int
    #Nx = 2**6
    
    c=0
    
    # Deletes old memory-map and creates new one
    try:
        os.remove(str(name)+'_phi.mymemmap')
    except FileNotFoundError:
        pass
    phi = np.memmap(str(name)+'_phi.mymemmap', dtype='float32', mode='w+', shape=(Nt,Nx,Nx))
    
    for p in range(1, parts+1):
        # Load the data files
        phi_fname = str(name) + '_phi_data'+str(p)+'.npy'
        phi_r = np.load(phi_fname)
    
        #print (phi.shape[0])
        try:
            if (p==1):
                for i in range (0, phi_r.shape[0]):
    #               print (c,i)
                    phi[c,:,:]=phi_r[i,:,:]
                    c = c+1
                        
            else:
                for i in range (0, phi_r.shape[0]-1):
    #               print(c,i+1)
                    phi[c,:,:]=phi_r[i+1,:,:]
                    c=c+1
        except:
            print("Error in the loop")

            

