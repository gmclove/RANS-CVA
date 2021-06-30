import numpy as np


class PtsList(object):
    pass

def cost(gamma, sigma):
    # load DNS data from file
    dx, dy, ijlpts, lpts, ptset, ptsetval = np.load('dnsData.npy', allow_pickle=True)
    # load RANS data from file 
    tracers = np.load(f'../base_case/tracers_cva-{gamma}-{sigma}.npy')
    eddies = np.load(f'eddies_cva-{gamma}-{sigma}.npy')

    ijlpts.    

    x = ptset[::2,0]
    y = ptset[::2,1]
    c = ptsetval[::2]
    
    x = tracers[:,i_x]
    y = tracers[:,i_y]
    c = tracers[:,i_inj]
    
    x = eddies[:,i_x]
    y = eddies[:,i_y]
    c = eddies[:,i_G]

    # compare DNS tracers to RANS-CVA tracers 
    
    # compare location of eddies in DNS to that of RANS-CVA
    # fields = h5py.File('2D_cyl.sol000400_1.sol.h5', 'r')
    # data = fields['Data']
    
    
cost(1, 0.05)
    