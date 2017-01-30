from numpy import *
from scipy import *
from pylab import *
import numpy.matlib as matlib
from shutil import copy
from os import mkdir
import os

def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  print a
  dx = dx0+arange(1.,n+1.,1.)*a
  return dx

Fr=0.13
H = 2000.
h0=500.
om = 2.*pi/12.42/3600.
N0=5.2e-3
u0=Fr*N0*h0


outdir='../runs/RunFr%03d' % (10000*Fr)
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
nz = 25

# y direction:
dy = 1000
# x direction
xt = 410e3

nmid = 50.
dx0=300.
nleft = (nx-nmid)/2
nright = (nx-nmid)/2
dx = zeros(nx)
dxleft = flipud(lininc(nleft,200.e3,dx0))
dxright = lininc(nright,200.e3,dx0)
dx[0:nleft]=dxleft
dx[(nleft):(nleft+nmid)]=dx0
dx[(nleft+nmid):]=dxright
x=cumsum(dx)
x = x-x[nx/2]

with open(outdir+"/delXvar.bin", "wb") as f:
	dx.tofile(f)
f.close()
# plot
if 1:
    plot(x/1000.,dx)
    xlim([-10,10])
    savefig(outdir+'/figs/dx.pdf')

# topo
sigma = 4000. # m

topo = 1500*exp(-x*x/(sigma**2))-1500+h0
#topo = h0*exp(-x*x/(3000**2))
print shape(topo)
topo[topo<0.]=0.
topo=-H+topo
topo[topo<-H]=-H

# plot
if 1:
    clf()
    plot(x/1.e3,topo)
    # xlim([-20.,20.])
    savefig(outdir+'/figs/topo.pdf')


with open(outdir+"/topo.bin", "wb") as f:
	topo.tofile(f)
f.close()
# dz:
# dz is from the surface down (right?).  Its saved as positive.

dz=zeros(nz)+H/nz

with open(outdir+"/delZvar.bin", "wb") as f:
	dz.tofile(f)
f.close()

# temperature profile...
g=9.8
alpha = 2e-4
T0 = 28+cumsum(N0**2/g/alpha*(-dz))

with open(outdir+"/TRef.bin", "wb") as f:
	T0.tofile(f)
f.close()

# save T0 over whole domain
TT0 = matlib.repmat(T0,nx,1).T
with open(outdir+"/T0.bin", "wb") as f:
	TT0.tofile(f)


z=cumsum(dz)
# plot:
if 1:
    clf()
    plot(T0,z)
    savefig(outdir+'/figs/TO.pdf')

# Forcing for boundaries
dt=3720.
time = arange(0,12.*3720.,dt)
print time/3600./12.4
om = 2*pi/12.40/3600;
uw = u0+0.*time
ue = u0+0.*time
# plot:
if 1:
    clf()
    plot(time/3600./12.4,ue,label='Ue')
    plot(time/3600/12.4,uw,label='Uw')
    legend()
    xlabel('t/T')
    ylabel('Vel')
    title('%d' % time[-1])
    savefig(outdir+'/figs/Vels.pdf')

# try time,nz,ny...

uen=zeros((shape(time)[0],nz,ny))
for j in range(0,ny):
  for i in range(0,nz):
    uen[:,i,j]=ue
#print(uen)

uwn=zeros((shape(time)[0],nz,ny))
print shape(uwn)
for j in range(0,ny):
  for i in range(0,nz):
    uwn[:,i,j]=uw
#print(uwn)

with open(outdir+"/Ue.bin","wb") as f:
  uen.tofile(f)

with open(outdir+"/Uw.bin", "wb") as f:
  uwn.tofile(f)

t=zeros((shape(time)[0],nz,ny))
for j in range(0,ny):
	for i in range(0,nz):
		for k in range(0,shape(time)[0]):
			t[k,i,j]=T0[i]
print shape(t)
with open(outdir+"/Te.bin", "wb") as f:
	t.tofile(f)
f.close()
with open(outdir+"/Tw.bin", "wb") as f:
	t.tofile(f)
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
