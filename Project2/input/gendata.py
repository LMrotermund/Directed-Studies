# 1d nonlinear hydrostatic flow over a gaussian bump - simulated by reduced gravity boundary as the free surface
# this reproduces the 5 regimes shown in Figure 2.11 of Baines and two other simulations with Fr= 0.8 and 1.2
# the 5 regimes are Supercritical flow. Subcritical flow. Partially blocked no lee jump. Partially blocked with lee jump. Complete blocking.
# need to vary Fr and H_m=h_m/d_o

from numpy import *
from scipy import *
from pylab import *
import numpy.matlib as matlib
from shutil import copy
from os import mkdir
import os


Fr=2                                                               #Variable
H = 1000
d0=100
H_m=0.25                                                           #Variable
h_m=H_m*d0
g=9.81

outdir='../runs/Run_trial_dT_2'

# 5 different regimes + Fr=0.8, Fr=1.2
#Supercritical, Fr=2 H_m=0.25
#Subcritical,   Fr=0.25 H_m=0.25
#Complete_blocking,   Fr=0.25 H_m=2
#Partially_blocked_with_lee_jump, Fr=0.25 H_m=1
#Partially_blocked_no_lee_jump,   Fr=1.2 H_m=1
#Fr0p8H0p25,  Fr=0.8 H_m=0.25
#Fr1p2H0p25   Fr=1.2 H_m=0.25

try:
  mkdir(outdir)
except:
  print outdir+' Exists'
try:
  mkdir(outdir+'/figs')
except:
  print outdir+'/figs Exists'
copy('./gendata.py',outdir)

# These must match ../code/SIZE.h
ny = 1
nx = 4*20
nz = 50

# y direction:
dy = 1000
# x direction
xt = 410e3

# dx
dx=zeros(nx)+1000
p1=np.arange(nx/2)
p2=np.flipud(p1)
p3=np.concatenate((p2,p1), axis=0)
dx=dx*1.02**p3
x=cumsum(dx)
x = x-x[nx/2]

with open(outdir+"/delXvar.bin", "wb") as f:
	dx.tofile(f)
f.close()

# plot
if 1:
    plot(x/1000.,dx)
    #xlim([-10,10])
    savefig(outdir+'/figs/dx.pdf')

# topo
sigma = 4000. # m
topo = h_m*exp(-x*x/(sigma**2))-H #correct?
print shape(topo)
topo[topo>0.]=0.

with open(outdir+"/topo.bin", "wb") as f:
	topo.tofile(f)
f.close()

# plot
if 1:
    clf()
    plot(x/1.e3,topo)
    # xlim([-20.,20.])
    savefig(outdir+'/figs/topo.pdf')

# dz
dz=zeros(nz)+4.78
p4=np.arange(nz)
p5=np.flipud(p4)
dz=dz*1.05**p5
z=cumsum(dz)
z = [(z[:-1]+z[1:])/2]

with open(outdir+"/delZvar.bin", "wb") as f:
	dz.tofile(f)
f.close()

# Temperature profile ... two layer
gravity=9.81
alpha = 2e-4
g_reduced=0.4E-3*gravity
dT=g_reduced/(gravity*alpha)
u0=Fr*sqrt(g_reduced*d0)

Tref=np.ones(nz)*dT
dz_flip=np.flipud(dz)
CumlativeSum=cumsum(dz_flip)
CumlativeSum_flip=np.flipud(CumlativeSum)

for j in range(0, nz):
    if CumlativeSum_flip[j] < d0:
        Tref[j]=0

#T0 = Tref[np.newaxis, :]*np.ones((nx,nz))
T0 = Tref[:,np.newaxis]*np.ones((nz,nx))

with open(outdir+"/TRef.bin", "wb") as f:
	Tref.tofile(f)
f.close()

# save T0 over whole domain

with open(outdir+"/T0.bin", "wb") as f:
	T0.tofile(f)

#initial velocity
U_initial=u0*np.ones(nz)
for j in range(0, nz):
    if CumlativeSum_flip[j] > d0:
        U_initial[j]=0

#U0=U_initial[np.newaxis, :]*np.ones((nx, nz))
U0=U_initial[:,np.newaxis]*np.ones((nz,nx))

with open(outdir+"/Uin.bin","wb") as f:
  U_initial.tofile(f)

with open(outdir+"/U0.bin", "wb") as f:
  U0.tofile(f)

# plot:
if 1:
    clf()
    plot(Tref,CumlativeSum)
    savefig(outdir+'/figs/Tref.pdf')


with open(outdir+"/Ue.bin","wb") as f:
  U_initial.tofile(f)

with open(outdir+"/Uw.bin", "wb") as f:
  U_initial.tofile(f)

with open(outdir+"/Te.bin", "wb") as f:
	Tref.tofile(f)
f.close()
with open(outdir+"/Tw.bin", "wb") as f:
	Tref.tofile(f)
f.close()

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
