import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.tri as tri

# def main():
#     xcfd, ycfd = loadMesh()
#     u_mean, v_mean, u_rms, v_rms = loadStats()


def loadMesh():
    unstrmesh = h5py.File('dump/2D_cyl.sol000000_1.mesh.h5','r')
    coords = unstrmesh['Coordinates']
    xcfd = np.array(coords['XYZ'][:,0])
    ycfd = np.array(coords['XYZ'][:,1])
    return xcfd, ycfd
    # Nr = 21
    # r = 0.5
    # xr = np.linspace(-r,r,Nr)
    # xadded = np.array([])
    # yadded = np.array([])
    # for i in range(Nr):
    #     for j in range(Nr):
    #         if (np.sqrt(xr[i]**2 + xr[j]**2) < r):
    #             xadded = np.append(xadded,xr[i])
    #             yadded = np.append(yadded,xr[j])
    # phiadded = np.zeros_like(xadded)
    # xcfd = np.append(xcfd,xadded)
    # ycfd = np.append(ycfd,yadded)

def loadStats():
    Ncfdpts = len(xcfd)
    cfdpts = np.zeros((Ncfdpts,2),dtype='float64')
    cfdpts[:,0] = xcfd
    cfdpts[:,1] = ycfd
    fields = h5py.File('dump/2D_cyl.sol000400_1.sol.h5', 'r')
    data = fields['Data']
    u_mean = np.array(data['U_MEAN'][:,0])
    # u_mean = np.append(u_mean,phiadded)
    v_mean = np.array(data['U_MEAN'][:,1])
    # v_mean = np.append(v_mean,phiadded)
    u_rms = np.array(data['U_RMS'][:,0])
    # u_rms = np.append(u_rms,phiadded)
    v_rms = np.array(data['U_RMS'][:,1])
    # v_rms = np.append(v_rms,phiadded)
    return u_mean, v_mean, u_rms, v_rms


def visualize():
    xreg = np.linspace(-10,50,600)
    yreg = np.linspace(-15,15,300)
    Xreg,Yreg = np.meshgrid(xreg,yreg)
    Um = griddata((x,y),u_mean,(Xreg,Yreg),method='cubic')
    plt.figure()
    plt.scatter(x,y,s = 0.1,c=u_mean)
    plt.show()
    plt.figure()
    plt.contourf(Xreg,Yreg,Um)
    plt.show()
    triang = tri.Triangulation(xcfd, ycfd)
    uminterpolator = tri.CubicTriInterpolator(triang, u_mean)
    Um = uminterpolator(Xreg,Yreg)
    vminterpolator = tri.CubicTriInterpolator(triang, v_mean)
    Vm = vminterpolator(Xreg,Yreg)
    urinterpolator = tri.CubicTriInterpolator(triang, u_rms)
    Ur = urinterpolator(Xreg,Yreg)
    vrinterpolator = tri.CubicTriInterpolator(triang, v_rms)
    Vr = vrinterpolator(Xreg,Yreg)
    fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
    axs[0,0].contourf(Xreg,Yreg,Um)
    axs[1,0].contourf(Xreg,Yreg,Vm)
    axs[0,1].contourf(Xreg,Yreg,Ur)
    axs[1,1].contourf(Xreg,Yreg,Vr)
    plt.show()

def loadParticles():
    particles = h5py.File("dump/2D_cyl.sol000400.ptset1_1.sol.h5",'r')
    ptset = particles['Coordinates']['XYZ']
    # ptset[0,:]
    ptsetval = particles['Data']['INJECTOR']
    return ptset[:], ptsetval[:]

# def loadVorticity():


xcfd, ycfd = loadMesh()
u_mean, v_mean, u_rms, v_rms = loadStats()
ptset, ptsetval = loadParticles()
