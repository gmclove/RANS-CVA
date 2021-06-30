import numpy as np
import matplotlib.pyplot as plt
# from IPython.display import display, clear_output

from interpRANS import IDW



class PtsList(object):
    pass

# load DNS data from file
dx, dy, ijlpts, lpts, ptset, ptsetval = np.load('../flowdata/dnsData.npy', allow_pickle=True)

i_x = 0; i_y = 1; i_u = 2; i_v = 3; i_G = 4; i_s = 5; i_t = 6; i_count = 7; i_D = 8; i_inj=4

xmin = -10
xmax = 50

class RANS_CVA:
    def __init__(self, gamma_max, sigma):
        self.gamma_max = gamma_max
        self.sigma = sigma
        # load DNS data from file
        # self.dx, self.dy, self.ijlpts, self.lpts, self.ptset, self.ptsetval = np.load('../flowdata/dnsData.npy', allow_pickle=True)
        # self.i_x = 0; self.i_y = 1; self.i_u = 2; self.i_v = 3; self.i_G = 4; self.i_s = 5; self.i_t = 6; self.i_count = 7; self.i_D = 8; self.i_inj=4
        # self.xmin = -10.
        # self.xmax = 50.
        

    def runSim(self):

        dims = 2
        trdims = 5
        eddims = 9
        
        dt = 0.1
        xmin = -10.
        xmax = 50.
        ymin = -15.
        ymax = 15.
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
        T_eddies = 3.
        it_eddy_inj = int(T_eddies / dt)
        x_eddy_inj = 1.0
        y_eddy_inj = 0.25
        # Drag = 1.0
        Gamma_max = self.gamma_max #1.0
        # c_growth = 2.0
        t_G_max = 4.0
        c_decay = 0.001
        sigma_inj = self.sigma #0.05
        def Gfunc(t):
            #global c_growth,c_decay #,Gamma_max
            if t <= t_G_max:
                G = Gamma_max*t/t_G_max
            else:
                G = np.max([Gamma_max*(1-c_decay*(t-t_G_max)),0.])
            return G
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
                eddies[ie,i_G] = (-1.0)**(eddies[ie,i_count]+1)*Gfunc(eddies[ie,i_t])
            i += 1\
            
        #     print(i)
            # if (i % it_inj) == 0:
            #     print(Nt,i,Ntracers,Neddies)
                # clear_output(wait=True)
                
        tracers_cva = np.copy(tracers)
        np.save(f'tracers_cva-{self.gamma_max}-{self.sigma}', tracers)
        np.save(f'eddies_cva-{self.gamma_max}-{self.sigma}', tracers)
        
        self.tracers = tracers_cva
        self.eddies = eddies
        self.plot()
        
        
    def plot(self):
        # tracers = np.load(f'tracers_cva-{self.gamma_max}-{self.sigma}.npy')
        # eddies = np.load(f'eddies_cva-{self.gamma_max}-{self.sigma}.npy')
        
        Nr = 2
        fig, ax = plt.subplots(nrows=Nr, ncols=1, constrained_layout=True,dpi=150)
        
        x = ptset[::2,0]
        y = ptset[::2,1]
        c = ptsetval[::2]
        ax[0].scatter(x, y, c=c, s=0.05,cmap='bwr')
        ax[0].set_title('DNS')
    
        
        x = self.tracers[:,i_x]
        y = self.tracers[:,i_y]
        c = self.tracers[:,i_inj]
        ax[1].scatter(x, y, c=c, s=0.05,cmap='bwr',vmin=0,vmax=1)
        ax[1].set_title(f'RANS-CVA: {self.gamma_max}, {self.sigma}')
        # fig.colorbar(scat,orientation='horizontal')
        x = self.eddies[:,i_x]
        y = self.eddies[:,i_y]
        c = self.eddies[:,i_G]
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
        plt.savefig(f'plt-{self.gamma_max}-{self.sigma}.png')
                    

# if __name__=='__main__':
#     runSim(1, 0.05)
#     # plot(1, 0.05)
