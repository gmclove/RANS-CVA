#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 18:51:52 2021

@author: dubief
"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
xmin = -10.
xmax = 50.
ymin = -15.
ymax = 15.
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
print(Xbins.shape)
i_inj1 = 0; i_inj2 = 1; ninj = 2
meanC = np.zeros((nxbins,nybins,2),dtype='float64')

isamplestart = 1
Nsamples = 400

for ifile in range(isamplestart,Nsamples+1):
    filename = "dump/2D_cyl.sol"+str(ifile).zfill(6)+".ptset1_1.sol.h5"
    # print(filename)
    particles = h5py.File(filename,'r')
    ptset = np.array(particles['Coordinates']['XYZ'])
    ptsetval = np.array(particles['Data']['INJECTOR'])
    particles.close()
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
    print(ptsetval.shape)
    print(ptsetval)
    I = ((ptset[:,0] - xmin)/binsize).astype(int)
    J = ((ptset[:,1] - ymin)/binsize).astype(int)
    mask = np.where(I > nxbins-1)
    I[mask] = nxbins-1
    Iinj = (ptsetval[:] - 1.).astype(int)
    meanC[I,J,Iinj] += 1./ptsetval.shape[0]
meanC /= Nsamples
print(meanC.shape)
np.save('meanC-dns', meanC)

fig, ax = plt.subplots(nrows=2, ncols=1, constrained_layout=True,dpi=150)
for iinj in range(ninj):
    ax[iinj].contourf(Ybins,Xbins,meanC[:,:,iinj])
    ax[iinj].set_aspect('equal')
plt.savefig('meanC-dns_contour.png')
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
plt.savefig('meanC-dns.png')
plt.show()
