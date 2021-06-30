import numpy as np

def IDW(ptsint,pts,phi):
    n = 2
    # R = 0.5
    if np.isscalar(phi):
        phiint = phi
    else:
        d = np.sqrt((ptsint[0] - pts[:,0])**2 + (ptsint[1] - pts[:,1])**2)    
        mask = np.where(d < 1e-6)
        if len(mask[0]) > 0:
            phiint = phi[mask]
        else:
            den = np.sum(1./d**n)
            num = np.sum(phi/d**2)
            phiint = num/den
#             if np.amin(d) < R :
#                 mask = np.where( d < R)
                
#                 print(mask,len(mask))
#                 den = np.sum(1./d[mask]**n)
#                 num = np.sum(phi[mask]/d[mask]**n)
#                 phi = num/den
#             else:
#                 il = np.argmin(d)
#                 phiint = phi[il]
    return phiint

def NN(ptsint,pts,phi):
    # n = 2
    # R = 0.2
    if np.isscalar(phi):
        phiint = phi
    else:
        d = np.sqrt((ptsint[0] - pts[:,0])**2 + (ptsint[1] - pts[:,1])**2)   
        il = np.argmin(d)
        phiint = phi[il]
    return phiint

