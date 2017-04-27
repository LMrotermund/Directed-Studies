#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 16:01:57 2017

@author: Lina
"""

from numpy import *
from scipy import *
from pylab import *
import numpy.matlib as matlib
from shutil import copy
from os import mkdir
import os
import math


# GRID DIMENSIONS
    
def lininc(n,Dx,dx0):
    a=(Dx-n*dx0)*2./n/(n+1)
    DX = dx0+arange(1.,n+1.,1.)*a
    return DX
    
ny = 1          # These must match ../code/SIZE.h
nx = 80
nz = 100

g=9.81
g_reduced=2E-3*g
D = 200          #total depth     
    
nmid = 50.
dx0=300.
nleft = (nx-nmid)/2
nright = (nx-nmid)/2
dx = zeros(nx)
dxleft = flipud(lininc(nleft,100.e3,dx0))
dxright = lininc(nright,100.e3,dx0)
dx[0:nleft]=dxleft
dx[(nleft):(nleft+nmid)]=dx0
dx[(nleft+nmid):]=dxright
x=cumsum(dx)
x = x-x[nx/2]
    
dz=zeros(nz)+D/nz
z=cumsum(dz)
        
X,Z=np.meshgrid(x,z,indexing='ij')

run = dict(
    # r, Hm, Fr
    # Unidirectional
    Fr_0p25_r_0p25_Hm_1=[0.25, 1, 0.25],
    Fr_0p75_r_0p25_Hm_1=[0.25, 1, 0.75],
    Fr_1_r_0p25_Hm_0p25=[0.25, 0.25, 1],
    Fr_0p15_r_0p35_Hm_1=[0.35, 1, 0.15],
    Fr_0p55_r_0p35_Hm_1=[0.35, 1, 0.55],
    Fr_1_r_0p35_Hm_1=[0.35, 1, 1],
    Fr_0p15_r_0p55_Hm_1=[0.55, 1, 0.15],
    Fr_0p45_r_0p55_Hm_0p75=[0.55, 0.75, 0.45],
    Fr_0p75_r_0p55_Hm_1=[0.55, 1, 0.75],
    Fr_0p25_r_0p75_Hm_0p25=[0.75, 0.25, 0.25],
    Fr_0p75_r_0p75_Hm_0p25=[0.75, 0.25, 0.75])


for name, (r, Hm, Fr) in run.items():
    
    outdir='../runs/'+name
    try:
      mkdir(outdir)
    except:
      print( outdir+' Exists')
    try:
      mkdir(outdir+'/figs')
    except:
      print( outdir+'/figs Exists')
    copy('./gendata.py',outdir)

    # Model parameters
     
    print(r)                 # VARIABLE
    d_1o=r*D         #thickness of lower layer
    d_2o=(1-r)*D     #thickness of upper layer
    print(Hm)                   # VARIABLE     
    hm=Hm*d_1o  
    print(Fr)                   # VARIABLE
    u_1o=Fr*np.sqrt(r*(1-r)*g_reduced*D)                # CHANGES b/c Fr is a Variable
    u_2o=u_1o
    
    # TOPOGRAPHY
    
    sigma = 4000.
    topo = hm*exp(-x*x/(sigma**2))-D # correct sign?
    topo[topo>0.]=0.
    
    with open(outdir+"/topo.bin", "wb") as f:
    		topo.tofile(f)
    f.close()

    # bottom layer
    
    boundary=np.ones(nx)*(d_2o)
    
    Dbottom= z< d_2o
    Dbottom2= Z < boundary [:,None]
    
    # Temperature profile
    alpha = 2e-4
    dT=g_reduced/(g*alpha)
    Tref=np.ones(nz)*dT
    Tref[Dbottom]=15
    T0 = dT * np.ones((nx, nz)) #T0 = Tref[:,np.newaxis]*np.ones((nz,nx))
    T0[Dbottom2] = 15
      
    # Initial Vel
    U_initial = u_2o * np.ones(nz)
    U_initial[Dbottom] = u_1o
    U0 = u_2o * np.ones((nx, nz)) #U0=U_initial[:,np.newaxis]*np.ones((nz,nx))
    U0[Dbottom2] = u_1o
    T0=T0.transpose()
    U0=U0.transpose()
    
    with open(outdir+"/delZvar.bin", "wb") as f:
        	dz.tofile(f)
    f.close()
    
    with open(outdir+"/delXvar.bin", "wb") as f:
        	dx.tofile(f)
    f.close()
    
    with open(outdir+"/TRef.bin", "wb") as f:
        Tref.tofile(f)
    f.close()
    
    with open(outdir+"/T0.bin", "wb") as f:
        T0.tofile(f)
    f.close()
    
    with open(outdir+"/Uin.bin","wb") as f:
        U_initial.tofile(f)
    f.close()
    
    with open(outdir+"/U0.bin", "wb") as f:
        U0.tofile(f)
    f.close()
    
    with open(outdir+"/Ue.bin","wb") as f:
        U_initial.tofile(f)
    f.close()
    
    with open(outdir+"/Uw.bin", "wb") as f:
        U_initial.tofile(f)
    f.close()
    
    with open(outdir+"/Te.bin", "wb") as f:
        Tref.tofile(f)
    f.close()
    
    with open(outdir+"/Tw.bin", "wb") as f:
        Tref.tofile(f)
    f.close()  

    with open(outdir+ "/summary.txt", "wt") as f:
        summary_str = 'r: {0}, Hm: {1:2.2f}, Fr: {2}'.format(r, Hm, Fr)
        #print(summary_str, file=f) 
    
    
    ## Copy some other files
    import shutil
    shutil.copy('data', outdir+'/data')
    shutil.copy('eedata', outdir)
    shutil.copy('data.kl10', outdir)
    shutil.copy('data.mnc', outdir)
    shutil.copy('data.obcs', outdir)
    shutil.copy('data.diagnostics', outdir)
    shutil.copy('data.pkg', outdir+'/data.pkg')
    # also store these.  They are small and helpful to document what we did
    for nm in {'input','code','build_options','analysis'}:
        to_path = outdir+'/'+nm
        if os.path.exists(to_path):
            shutil.rmtree(to_path)
        shutil.copytree('../'+nm, outdir+'/'+nm)