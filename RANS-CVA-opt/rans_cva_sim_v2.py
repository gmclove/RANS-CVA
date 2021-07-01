import numpy as np
import os
import matplotlib.pyplot as plt
from pymooCFD.setupOpt import dataDir
# from IPython.display import display, clear_output

from interpRANS import IDW

class PtsList(object):
    pass

# load DNS data from file
dx, dy, ijlpts, lpts, ptset, ptsetval = np.load('./flowdata/dnsData.npy', allow_pickle=True)

i_x = 0; i_y = 1; i_u = 2; i_v = 3; i_G = 4; i_s = 5; i_t = 6; i_count = 7; i_D = 8; i_inj=4

xmin = -10.
xmax = 50.
ymin = -15.
ymax = 15.

def run(var):
    global caseDir, indDir
    indDir = f'ind{var[-1]}'
    caseDir = f'{dataDir}/{indDir}'
    print(f'### running {caseDir}')
    try:
        os.mkdir(caseDir)
    except OSError as err:
        print(err)
    # global caseDir
    # caseDir = path
    # if os.path.exists(f'{caseDir}/tracers_all-{var[0]}-{var[1]}'):
    try:
        tracers = np.load(f'{caseDir}/tracers_all-{var[0]}-{var[1]}.npy', allow_pickle=True)
        print(f'{indDir} simulation already complete')
    except FileNotFoundError as err:
        print(err)
        tracers = runSim(var)
        try:
            tracers = np.load(f'{caseDir}/tracers_all-{var[0]}-{var[1]}.npy', allow_pickle=True)
        except FileNotFoundError as err:
            raise err

    obj = postProc(tracers)
    saveObj(obj)

def saveObj(obj):
    try:
        prev_obj = np.load(f'{caseDir}/obj.npy')
        new_obj = prev_obj + obj
        np.save(f'{caseDir}/obj.npy', new_obj)
    except FileNotFoundError as err:
        print(err)
        print(f'### creating {caseDir}/obj.npy file')
        np.save(f'{caseDir}/obj.npy', obj)



def runSim(var):
    gamma = var[0]
    sigma = var[1]

    dims = 2
    trdims = 5
    eddims = 9

    dt = 0.1
    ymin_inj1 = -1.0
    ymax_inj1 = -0.1
    ymin_inj2 = 0.1
    ymax_inj2 = 1.0
    Npt_inj = 5
    y_inj1 = np.linspace(ymin_inj1,ymax_inj1,Npt_inj)
    y_inj2 = np.linspace(ymin_inj2,ymax_inj2,Npt_inj)
    # x_inj = xmin*np.ones(Npt_inj)
    tracers_inj1 = np.zeros((Npt_inj,trdims),dtype='float64')
    tracers_inj1[:,i_y] = y_inj1
    tracers_inj1[:,i_x] = xmin
    tracers_inj1[:,i_u] = 1.0
    tracers_inj1[:,i_inj] = 1.0
    tracers_inj2 = np.zeros((Npt_inj,trdims),dtype='float64')
    tracers_inj2[:,i_y] = y_inj2
    tracers_inj2[:,i_x] = xmin
    tracers_inj2[:,i_u] = 1.0
    tracers_inj2[:,i_inj] = 0.0
    T_tracer_inj = 0.2
    it_inj = int(T_tracer_inj / dt)
    T_eddies = 1/0.18527791687531295  #3.
    it_eddy_inj = int(T_eddies / dt)
    x_eddy_inj = 1.0
    y_eddy_inj = 0.25
    sigma_inj = sigma #0.05
    # Drag = 1.0

    # Gamma_max = 1.0
    # # c_growth = 2.0
    # t_G_max = 4.0
    # c_decay = 0.001
    # def Gfunc(t):
    #     #global c_growth,c_decay #,Gamma_max
    #     if t <= t_G_max:
    #         G = Gamma_max*t/t_G_max
    #     else:
    #         G = np.max([Gamma_max*(1-c_decay*(t-t_G_max)),0.])
    #     return G

    eddy_inj = np.zeros((1,eddims),dtype='float64')
    eddy_inj[0,i_x] = x_eddy_inj
    eddy_inj[0,i_y] = y_eddy_inj
    eddy_inj[0,i_u] = 0.0
    eddy_inj[0,i_v] = 0.0
    eddy_inj[0,i_G] = 0.0
    eddy_inj[0,i_s] = sigma_inj
    eddy_inj[0,i_t] = 0.0
    Tsimu = 100.
    Nt = int(Tsimu / dt)

    i = 0
    # iphase = 0
    ie_count = 0
    tracers_all = [] #np.dtype(object)#np.array(0, dtype=object) #[] # tracers at every iteration
    eddies_all = []
    while i < Nt:
        if (i == 0):
            tracers = np.copy(tracers_inj1)
            tracers = np.insert(tracers,np.shape(tracers)[0],tracers_inj2,axis=0)
        elif (i % it_inj) == 0:
            tracers = np.insert(tracers,np.shape(tracers)[0],tracers_inj1,axis=0)
            tracers = np.insert(tracers,np.shape(tracers)[0],tracers_inj2,axis=0)
        Ntracers = np.shape(tracers)[0]
        if (i == 0):
            ie = 0
            eddies = np.copy(eddy_inj)
            eddies[ie,i_y] *= (-1.0)**ie_count
            eddies[ie,i_count] = ie_count
            Neddies = 1
            ie_count += 1
        elif (i % it_eddy_inj) == 0:
            eddies = np.insert(eddies,np.shape(eddies)[0],eddy_inj,axis=0)
            Neddies = np.shape(eddies)[0]
            ie = Neddies - 1
            eddies[ie,i_y] *= (-1.0)**ie_count
            eddies[ie,i_count] = ie_count
            ie_count += 1
        Ntracers = np.shape(tracers)[0]
    #     print(Ntracers)
        # RANS velocity
    #     for j in range(Ntracers):
        pts = np.copy(tracers[:,:i_y+1])
        itracers = (tracers[:,i_x] - xmin)/dx
        itracers = itracers.astype(int)
    #     print(itracers)
        jtracers = (tracers[:,i_y] - ymin)/dy
        jtracers = jtracers.astype(int)
        for it in range(len(itracers)):
            il = ijlpts[itracers[it],jtracers[it]]
            Nlistpts = lpts[il].Npts
            itpts = np.zeros((Nlistpts,2))
            itpts[:,0] = lpts[il].xlist
            itpts[:,1] = lpts[il].ylist
            itu_mean = lpts[il].ulist
            itv_mean = lpts[il].vlist
    #         print(itu_mean,lpts[il].ulist,il)
            tracers[it,i_u] = IDW(tracers[it,:i_y+1],itpts,itu_mean)
            tracers[it,i_v] = IDW(tracers[it,:i_y+1],itpts,itv_mean)
    #     tracers[:,i_u] = np.ones_like(tracers[:,i_u]) #np.diag(f_Um(tracers[:,i_x],tracers[:,i_y]))
    #     tracers[:,i_v] = np.zeros_like(tracers[:,i_u]) #np.diag(f_Vm(tracers[:,i_x],tracers[:,i_y]))
        # Eddy-Induced velocity
        for ie in range(Neddies):
            rvec = np.zeros((Ntracers,dims),dtype='float64')
            for j in range(dims):
                rvec[:,j] = tracers[:,j] - eddies[ie,j]
            r2 = rvec[:,i_x]**2+rvec[:,i_y]**2
            mask = np.where(r2 != 0.0)
            r2inv = np.zeros_like(r2)
            r2inv[mask] = 1.0/r2[mask]
            GKs = eddies[ie,i_G]/(2*np.pi)*r2inv*(1.0 - np.exp(-r2/(2*eddies[ie,i_s]**2)))
            tracers[:,i_u] -= GKs*rvec[:,i_y]
            tracers[:,i_v] += GKs*rvec[:,i_x]
        tracers[:,i_x:i_y+1] += dt*tracers[:,i_u:i_v+1]
        # Boundary conditions
        maskx = np.where(tracers[:,i_x] > xmax)
        if (len(maskx[0]) > 0):
            tracers = np.delete(tracers,maskx[0],axis=0)
        masky = np.where(np.abs(tracers[:,i_y]) > ymax)
        if (len(masky[0]) > 0):
            tracers = np.delete(tracers,masky[0],axis=0)
        pts = np.copy(tracers[:,:i_y+1])
        ieddies = (eddies[:,i_x] - xmin)/dx
        ieddies = ieddies.astype(int)
    #     print(itracers)
        jeddies = (eddies[:,i_y] - ymin)/dy
        jeddies = jeddies.astype(int)
        for ie in range(Neddies):
            il = ijlpts[ieddies[ie],jeddies[ie]]
            Nlistpts = lpts[il].Npts
            itpts = np.zeros((Nlistpts,2))
            itpts[:,0] = lpts[il].xlist
            itpts[:,1] = lpts[il].ylist
            itu_mean = lpts[il].ulist
            itv_mean = lpts[il].vlist
    #         print(itu_mean,lpts[il].ulist,il)
            eddies[ie,i_u] = IDW(eddies[ie,:i_y+1],itpts,itu_mean)
            eddies[ie,i_v] = IDW(eddies[ie,:i_y+1],itpts,itv_mean)
    #     eddies[:,i_u] = 0.#bruteIDW(eddies[:,:i_y+1],u_mean)
    #     eddies[:,i_v] = 0.#bruteIDW(eddies[:,:i_y+1],v_mean)
        for ie in range(Neddies):
            for ie1 in range(Neddies):
                if (ie != ie1):
                    rvec = np.zeros(2,dtype='float64')
                    rvec = eddies[ie,i_x:i_y+1] - eddies[ie1,i_x+i_y+1]
                    r2 = np.sum(rvec**2)
                    r2inv = 1.0/r2
                    GKs = eddies[ie1,i_G]/(2*np.pi)*r2inv*(1.0 - np.exp(-r2/(2*eddies[ie1,i_s]**2)))
                    eddies[ie,i_u] -= GKs*rvec[i_y]
                    eddies[ie,i_v] += GKs*rvec[i_x]
        eddies[:,i_x:i_y+1] += dt*eddies[:,i_u:i_v+1]
        maskx = np.where(eddies[:,i_x] > xmax)
        if len(maskx[0]) > 0 :
            eddies = np.delete(eddies,maskx[0],axis=0)
        masky = np.where(np.abs(eddies[:,i_y]) > ymax)
        if len(masky[0]) > 0 :
            eddies = np.delete(eddies,masky[0],axis=0)
        Neddies = np.shape(eddies)[0]
        eddies[:,i_t] += dt
        for ie in range(Neddies):
            eddies[ie,i_G] = (-1.0)**(eddies[ie,i_count]+1)*gamma #*Gfunc(eddies[ie,i_t])
        i += 1

        # store tracer info
        # tracers_all = np.append(tracers_all, tracers)
        tracers_all.append(tracers)
        # print(tracers.shape)
        # print(tracers)
        # print(tracers_all.shape)
        # print(len(tracers_all), len(tracers_all[0]), len(tracers_all[0][0]))
        # print(tracers_all)

        eddies_all.append(eddies)

    #     print(i)
        # if (i % it_inj) == 0:
        #     print(Nt,i,Ntracers,Neddies)
            # clear_output(wait=True)

    # SIMULATION COMPLETE
    np.save(f'{caseDir}/tracers_all-{var[0]}-{var[1]}', tracers_all)
    # np.save(f'{caseDir}/tracers_cva-{var[0]}-{var[1]}', tracers)
    # np.save(f'{caseDir}/eddies_cva-{var[0]}-{var[1]}', eddies)

    plot(gamma, sigma, tracers, eddies)
    return tracers

def postProc(tracers_all):
    ############################################################################
    ###### POST-PROCESS ######
    # create bins to measure distribution of tracer particles across the domain
    binsize = 0.5
    Lx = xmax - xmin
    Ly = ymax - ymin
    nxbins = int(Lx/binsize)
    nybins = int(Ly/binsize)
    xsbins = np.linspace(xmin,xmax,nxbins+1)
    ysbins = np.linspace(ymin,ymax,nybins+1)
    xbins = (xsbins[1:] + xsbins[:-1])/2.
    ybins = (ysbins[1:] + ysbins[:-1])/2.
    Xbins,Ybins = np.meshgrid(ybins,xbins)
    # print(Xbins.shape)
    i_inj1 = 0; i_inj2 = 1; ninj = 2
    meanC = np.zeros((nxbins,nybins,2),dtype='float64')

    isamplestart = 1 - 1
    Nsamples = len(tracers_all)

    for it in range(isamplestart,Nsamples):
        # filename = "dump/2D_cyl.sol"+str(ifile).zfill(6)+".ptset1_1.sol.h5"
        # # print(filename)
        # particles = h5py.File(filename,'r')
        # ptset = np.array(particles['Coordinates']['XYZ'])
        # ptsetval = np.array(particles['Data']['INJECTOR'])
        # particles.close()
        # for n in range(len(ptsetval)):
            # i = int((ptset[n,0] - xmin)/binsize)
            # j = int((ptset[n,1] - ymin)/binsize)
            # if ptsetval[n] == 1:
            #     iinj = i_inj1
            # elif ptsetval[n] == 2:
            #     iinj = i_inj2
            # else:
            #     print("problem")
            # meanC[i,j,iinj] += 1
        # plot(gamma, sigma, tracers_all[it], eddies_all[it])
        tracers = tracers_all[it]
        # print(tracers_all)
        # tracers = np.array(tracers_all[it])
        # print(tracers.type)
        # print(tracers.shape)
        # print(tracers[:,i_x].shape)
        # print(tracers[:][i_x].shape)
        # print(tracers.shape[0])

        # I = ((tracers[:][i_x] - xmin)/binsize).astype(int)
        # J = ((tracers[:][i_y] - ymin)/binsize).astype(int)
        I = ((tracers[:,i_x] - xmin)/binsize).astype(int)
        J = ((tracers[:,i_y] - ymin)/binsize).astype(int)
        mask = np.where(I > nxbins-1)
        I[mask] = nxbins-1
        # Iinj = (tracers[:][i_inj] - 1.).astype(int)
        Iinj = (tracers[:, i_inj] - 1.).astype(int)
        meanC[I,J,Iinj] += 1./tracers.shape[0]
    meanC /= Nsamples
    meanC_cva = meanC
    np.save(f'{caseDir}/meanC-cva', meanC)
    # print(meanC.shape)
    # return meanC
    # np.save(f'{caseDir}/meanC-cva', meanC)

    # print(meanC[20,:,i_inj1].shape)
    # meanC[20,:,i_inj2]
    # meanC[30,:,i_inj1]
    # meanC[30,:,i_inj2]
    meanC_dns = np.load('flowdata/meanC-dns.npy')
    obj = []
    for i, x in enumerate([10, 30]):
        meanDif1 = np.mean(np.abs(meanC_cva[x,:,i_inj1]-meanC_dns[x,:,i_inj1]))
        meanDif2 = np.mean(np.abs(meanC_cva[x,:,i_inj2]-meanC_dns[x,:,i_inj2]))
        # print(meanDif1, meanDif2)
        meanDif = np.mean([meanDif1, meanDif2])
        obj.append(meanDif)
    print(f'{caseDir} objectives: {obj}')
    return obj

def plotBins():
    fig, ax = plt.subplots(nrows=2, ncols=1, constrained_layout=True,dpi=150)
    for iinj in range(ninj):
        ax[iinj].contourf(Ybins,Xbins,meanC[:,:,iinj])
        ax[iinj].set_aspect('equal')
    plt.savefig(f'{caseDir}/meanC-cva_contour.png')
    plt.show()

    fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True,dpi=150)
    ax.plot(ybins,meanC[20,:,i_inj1],'k-', label='x=20')
    ax.plot(ybins,meanC[20,:,i_inj2],'k--', label='x=20')
    ax.plot(ybins,meanC[30,:,i_inj1],'C0-', label='x=30')
    ax.plot(ybins,meanC[30,:,i_inj2],'C0--', label='x=30')
    ax.plot(ybins,meanC[40,:,i_inj1],'C1-', label='x=40')
    ax.plot(ybins,meanC[40,:,i_inj2],'C1--', label='x=40')
    ax.legend()
    ax.set_xlim(-6,6)
    plt.savefig(f'{caseDir}/meanC-cva.png')
    plt.show()

def plot(gamma_max, sigma, tracers, eddies):
    # tracers = np.load(f'tracers_cva-{self.gamma_max}-{self.sigma}.npy')
    # eddies = np.load(f'eddies_cva-{self.gamma_max}-{self.sigma}.npy')

    Nr = 2
    fig, ax = plt.subplots(nrows=Nr, ncols=1, constrained_layout=True,dpi=150)

    x = ptset[::2,0]
    y = ptset[::2,1]
    c = ptsetval[::2]
    ax[0].scatter(x, y, c=c, s=0.05,cmap='bwr')
    ax[0].set_title('DNS')


    x = tracers[:,i_x]
    y = tracers[:,i_y]
    c = tracers[:,i_inj]
    ax[1].scatter(x, y, c=c, s=0.05,cmap='bwr',vmin=0,vmax=1)
    ax[1].set_title(f'RANS-CVA: {gamma_max}, {sigma}')
    # fig.colorbar(scat,orientation='horizontal')
    x = eddies[:,i_x]
    y = eddies[:,i_y]
    c = eddies[:,i_G]
    # ax[1].scatter(x,y, c=c, s=20, marker = "o",cmap='RdBu_r',vmin=-0.9*Gamma_max,vmax=0.9*Gamma_max)


    # x = tracers_rans[:,i_x]
    # y = tracers_rans[:,i_y]
    # c = tracers_rans[:,i_inj]
    # ax[2].scatter(x, y, c=c, s=0.05,cmap='bwr',vmin=0,vmax=1)
    # ax[2].set_title('RANS')


    for ir in range(Nr):
        ax[ir].set_aspect('equal')
        ax[ir].set_xlim(xmin,xmax)
        ax[ir].set_ylim(-5,5)

    # plt.show()
    plt.savefig(f'{caseDir}/plt-{gamma_max}-{sigma}.png')

# run([1, 0.05])
