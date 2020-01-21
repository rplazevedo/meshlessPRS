import numpy as np
import os
import configparser

def run():
    config = configparser.ConfigParser()
    config.read('input.ini')
    name = config['Parameters']['name']
    
    dt = float(config['PRS']['dt'])
    
    parts = int(config['PRS']['sPRS_parts'])+int(config['PRS']['mPRS_parts'])
    Np = int(config['PRS']['points_per_part'])
    
    N = parts*Np-parts+1
    
    y = np.empty(N)
    
    reduce_data = False
    
    reduction_number = 2
    
    for i in range(0, y.shape[0]):
        y[i] = dt*i
    
    v_data = np.zeros(y.shape[0])
    nop_data = np.zeros(y.shape[0])
    nwalls_data = np.zeros(y.shape[0])
    exc_pts_data = np.zeros(y.shape[0])
    
    number = 0
    while os.path.isfile(str(name) +'_'+str(number)+ '_v.data'):
        number += 1
    
    v_predata = np.loadtxt(str(name) +'_'+str(number)+ '_v.txt')
    nwalls_predata = np.loadtxt(str(name) +'_'+str(number)+ '_nwalls.txt')
    exc_pts_predata = np.loadtxt(str(name) +'_'+str(number)+ '_exc_pts.txt')
    try:
        nop_predata = np.loadtxt(str(name) +'_'+str(number)+ '_nop.txt')
    except FileNotFoundError:
        pass
    
    
    
    for i in range(0, y.shape[0]):
        v_data[i] = v_predata[i]
        try:
            nop_data[i] = nop_predata[i]
        except:
            pass
        nwalls_data[i] = nwalls_predata[i]
        exc_pts_data[i] = exc_pts_predata[i]
        
    if not reduce_data:
    
        np.savetxt(str(name) +'_'+str(number)+ '_v.data', np.c_[y,v_data])
        np.savetxt(str(name) +'_'+str(number)+ '_nwalls.data', np.c_[y,nwalls_data])
        np.savetxt(str(name) +'_'+str(number)+ '_exc_pts.data', np.c_[y,exc_pts_data])
        try:
            np.savetxt(str(name) +'_'+str(number)+ '_nop.data', np.c_[y,nop_data])
        except:
            pass
    
    else:
    
        old_v_data = np.c_[y,v_data]
        try:
            old_nop_data = np.c_[y,nop_data]
        except:
            pass
        old_nwalls_data = np.c_[y,nwalls_data]
        old_exc_pts_data = np.c_[y,exc_pts_data]
    
        v_new_data = np.zeros(( int(v_data.shape[0]/reduction_number),2) )
        try:
            nop_new_data = np.zeros(( int(v_data.shape[0]/reduction_number),2) )
        except:
            pass
    
        nwalls_new_data = np.zeros(( int(nwalls_data.shape[0]/reduction_number),2) )
        exc_pts_new_data = np.zeros(( int(exc_pts_data.shape[0]/reduction_number),2) )
    
        for i in range(0, v_new_data.shape[0]):
            
            try:
                v_new_data[i] = old_v_data[i*reduction_number]
                try:
                    nop_new_data[i] = old_nop_data[i*reduction_number]
                except:
                    pass
                
                nwalls_new_data[i] = old_nwalls_data[i*reduction_number]
                exc_pts_new_data[i] = old_exc_pts_data[i*reduction_number]
                
            
            except:
                print("Error at: " + str(i))
                
        np.savetxt(str(name) +'_'+str(number)+ '_v.data', v_new_data)
        np.savetxt(str(name) +'_'+str(number)+ '_nwalls.data', nwalls_new_data)
        np.savetxt(str(name) +'_'+str(number)+ '_exc_pts.data', exc_pts_new_data)
        try:
            np.savetxt(str(name) +'_'+str(number)+ '_nop.data', nop_new_data)
        except:
            pass


