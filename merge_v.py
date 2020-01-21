import numpy as np
import os
import configparser


DTYPE = np.float32

config = configparser.ConfigParser()
config.read('input.ini')


def run():
    Nt_int = int(config['PRS']['points_per_part'])
    parts = int(config['PRS']['sPRS_parts'])+int(config['PRS']['mPRS_parts'])
    name = str(config['Parameters']['name'])
    
    Nt=parts*Nt_int
    
    c=0
    
    # Creates memory-maps
    v = np.memmap(str(name) + '_v.mymemmap', dtype='float32', mode='w+', shape=(Nt))
    nwalls = np.memmap(str(name) + '_nwalls.mymemmap', dtype='float32', mode='w+', shape=(Nt))
    nop = np.memmap(str(name) + '_nop.mymemmap', dtype='float32', mode='w+', shape=(Nt))
    exc_pts = np.memmap(str(name) + '_exc_pts.mymemmap', dtype='float32', mode='w+', shape=(Nt))

    
    for p in range(1, parts+1):
        # Load the data files
        v_fname = str(name) + '_vdata' +str(p)+'.npy'
        v_r = np.load(v_fname)
        nwalls_fname = str(name) + '_nwalls_data' + str(p) + '.npy'
        nwalls_r = np.load(nwalls_fname)
        exc_pts_fname = str(name) + '_exc_pts_data' +str(p)+'.npy'
        try:
            exc_pts_r = np.load(exc_pts_fname)
        except:
            exc_pts_r = np.zeros(v_r.shape)
        try:
            nop_fname = str(name) + '_nop_data' + str(p) + '.npy'
            # print(str(name) + '_nop_data' + str(p) + '.npy')
            nop_r = np.load(nop_fname)
        except:
            pass
    
        if (p==1):
            for i in range (0, v_r.shape[0]):
                # print(c,i)
                v[c]=v_r[i]
                exc_pts[c] = exc_pts_r[i]
                nwalls[c]=nwalls_r[i]
                try:
                    nop[c] = nop_r[i]
                except: 
                    pass
                c = c+1
                        
        else:
            for i in range (0, v_r.shape[0]-1):
                # print(c,i+1)
    
                v[c]=v_r[i+1]
                nwalls[c]=nwalls_r[i]
                exc_pts[c] = exc_pts_r[i]
                try:
                    nop[c]=nop_r[i]
                except:
                    pass
    
                c=c+1
        
    # Save new .txt files
    number = 0
    while os.path.isfile(str(name) +'_'+str(number)+ '_v.txt'):
        number+=1
    np.savetxt(str(name) +'_'+str(number)+ '_v.txt', v)
    np.savetxt(str(name) +'_'+str(number)+ '_nwalls.txt', nwalls)
    np.savetxt(str(name) +'_'+str(number)+ '_exc_pts.txt', exc_pts)    
    try:
        np.savetxt(str(name) +'_'+str(number)+ '_nop.txt', nop)
    except:
        pass


