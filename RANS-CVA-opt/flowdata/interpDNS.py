from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline
import numpy as np 

# from extractData import 
from extractData import xcfd, ycfd, u_mean, u_rms, v_mean, v_rms, ptset, ptsetval, vort

if __name__ == "__main__":
    xmin = -10.; xmax = 50.
    ymin = -15.; ymax = 15.
    xbnds = np.linspace(xmin,xmax,61)
    ybnds = np.linspace(ymin,ymax,31)
    xcn = (xbnds[1:] + xbnds[:-1])/2.
    ycn = (ybnds[1:] + ybnds[:-1])/2.
    dx = xbnds[1] - xbnds[0]
    dy = ybnds[1] - ybnds[0]
    # print(dx,dy)
    class PtsList(object):
        def __init__(self,i,j):
            global u_mean,v_mean,u_rms,v_rms,xbnds,ybnds,xcn,ycn,xcfd,ycfd
            self.i = i
            self.j = j
            self.xmin = xbnds[i]
            self.xmax = xbnds[i+1]
            self.ymin = ybnds[j]
            self.ymax = ybnds[j+1]
            self.xcenter = xcn[i]
            self.ycenter = ycn[j]
            listpts = []
            xlist = []
            ylist = []
            ulist = []
            vlist = []
            klist = []
            Npts = 0
            for k in range(len(xcfd)):
                if (xcfd[k] >= self.xmin) and (xcfd[k] < self.xmax):
                    if (ycfd[k] >= self.ymin) and (ycfd[k] < self.ymax):
                        listpts.append(k)
                        xlist.append(xcfd[k])
                        ylist.append(xcfd[k])
                        ulist.append(u_mean[k])
                        vlist.append(v_mean[k])
                        klist.append((u_rms[k]**2 + v_rms[k]**2)/2.)
                        Npts += 1
            if Npts == 0:
                self.xlist = [xcn[j]]
                self.ylist = [ycn[j]]
                self.ulist = [1.0]
                self.vlist = [1.0]
                self.klist = [0.0]
                self.listpts = [-1]
                Npts = 1
            else:
                self.xlist = xlist
                self.ylist = ylist
                self.ulist = ulist
                self.vlist = vlist
                self.klist = klist
                self.listpts = listpts
            
            self.Npts = Npts
            
            
    Nxgrid = len(xcn)
    Nygrid = len(ycn)
    lpts = []
    il = 0
    ijlpts = np.zeros((Nxgrid,Nygrid),dtype='int')
    for i in range(Nxgrid):
        for j in range(Nygrid):
            lpts.append(PtsList(i,j))
    #         print(lpts[il].Npts)
            ijlpts[i,j] = il
            il += 1
            
    for il in range(len(lpts)):
        if lpts[il].Npts == 0:
            print(il,lpts[il].xmin,lpts[il].ymin)
    
    # Save to file
    dat = np.array([dx, dy, ijlpts, lpts, ptset, ptsetval], dtype=object)
    np.save('dnsData.npy', dat) #, ptset, ptsetval])
    



        