#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 18:43:29 2020

@author: marc

This routine will set up the size of the domain, the size of the vortex, and
other parameters needed to start a run
"""
## Initial Setup Parameters for model run

n = 301                 # number of grid points in the x-direction, must be odd
m = 301                 # number of grid points in the y-direction, must be odd
dx = 2.0e3              # x-direction grid spacing in meters
dy = 2.0e3              # y-direction grid spacing in meters
dt = 7.5                # timestep in seconds
g = 9.81                # acceleration due to gravity in m/s^2
f0 = 5.0e-5             # Coriolis parameter in s^-1
beta = 0.0              # beta parmeter in 1/(m*s)
omega = 7.292E-5        # angular rotaion rate of the earth 1/s
nu = 100.0              # viscosity
maxTimeCount = 11520    # max time count iterations
outputInc = 960         # number of increments between model output
order = 0               # 2nd order (0) or 4th order (1)
state = 0               # ring (0) or monopole (1) or cimcumpolar (2)
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
from IPython.display import display,clear_output
import time as Time
import math, os
import numpy as np
#import scipy.fftpack
import scipy
import imageio
import matplotlib
matplotlib.rcParams.update({'font.size': 22})
from IPython.core.display import HTML
import urllib.request

## Now that we have some key parameters defined

# compute angle and radius values at each array element
rows,cols = m,n
xmdpt = int((n-1)/2)
ymdpt = int((m-1)/2)
xdist = np.arange(cols)
xdist = xdist-xmdpt
x = np.zeros((m,n))
for j in range(rows):
    x[j] = xdist*dx
ydist = np.arange(rows)
ydist = ymdpt-ydist
y = np.zeros((m,n))
for i in range(cols):
    y[i] = ydist*dy
y = np.transpose(y)

radToDeg = 180.0/np.pi
radius = np.zeros((rows,cols))
radius = np.sqrt(x*x + y*y)

import pandas as pd 
pd.DataFrame(radius).to_csv("radius.csv", header=None, index=None)

theta = np.zeros((rows,cols))
for i in range(cols):
    for j in range(rows):
        if (xdist[i] != 0.0) :
            theta[j,i] = np.arctan(ydist[j]/xdist[i])
            if (ydist[j] >= 0.0) and (xdist[i] < 0.0):
                theta[j,i] = np.pi - abs(theta[j,i])
            elif (ydist[j] < 0.0):
                if (xdist[i] < 0.0):
                    theta[j,i] = np.pi + abs(theta[j,i])
                else:
                    theta[j,i] = 2.0*np.pi - abs(theta[j,i])
        else :
            if (ydist[j] >= 0.0):
                theta[j,i] = np.pi/2.0
            else :
                theta[j,i] = 3.0*np.pi/2.0
theta[ymdpt,xmdpt] = -2.0*np.pi
pd.DataFrame(theta).to_csv("theta.csv", header=None, index=None)

# build the ring vortex profile from Schubert et. al
r1 = 15.0e3
r2 = 22.5e3
r3 = 25.0e3
r4 = 32.5e3
zeta1 = 4.1825e-4
zeta2 = 7.000808e-3
zeta3 = 1.0e-5
zetaZero = 0.0
zeta = np.zeros((m,n))

print("compute basicsate for ring vortex")
# assign the basicstate vorticity values
for j in range(n):
    for i in range(m):
        sumCos = 0.0
        for k in range(1,8):
            sumCos = sumCos + np.cos(k*theta[i,j])

        r = radius[i,j]
        if (r <= r1):
            zeta[i,j] = zeta1
        elif ((r > r1) and (r <= r2)):
            fn1 = (r-r1)/(r2-r1)
            fn2 = (r2-r)/(r2-r1)
            fn3 = (r2-r+1)/(r2-r1)
            s1 = 1.0 - 3.0*(fn1*fn1) + 2.0*(fn1*fn1*fn1)
            s2 = 1.0 - 3.0*(fn2*fn2) + 2.0*(fn2*fn2*fn2)
            s3 = 1.0 - 3.0*(fn3*fn3) + 2.0*(fn3*fn3*fn3)
            zeta[i,j] = zeta1*s1+zeta2*s2+zeta3*sumCos*s3
        elif ((r > r2) and (r <= r3)):
            zeta[i,j] = zeta2 + zeta3*sumCos
        elif ((r > r3) and (r <= r4)):
            fn2 = (r-r3)/(r4-r3)
            fn3 = (r-r3)/(r4-r3)
            s2 = 1.0 - 3.0*(fn2*fn2) + 2.0*(fn2*fn2*fn2)
            s3 = 1.0 - 3.0*(fn3*fn3) + 2.0*(fn3*fn3*fn3)
            zeta[i,j] = zetaZero + zeta2*s2 + zeta3*sumCos*s3
    
# equate the boundaries so it's doubly periodic
for j in range(n):
    zeta[m-1,j] = zeta[0,j]
for i in range(m):
    zeta[i,n-1] = zeta[i,0]

# ensure net zero circulation. first compute the area
area = ((m-1)*dx)*((n-1)*dy)

# compute the area integrated vorticity divided by the total area
# to get the circulation correction factor. The integral is via 
# trapezoidal quadrature. Correct the vorticity so that there's
# zero net circulation

print("correct vorticity for zero net circulation")
integral = 0.0
for j in range (n - 1):
    for i in range (m - 1):
        integral = integral + 0.25*dx*dy*\
        (zeta[i,j]+zeta[i,j+1]+zeta[i+1,j]+zeta[i+1,j+1])
correction = integral/area
      
# Then use the correction factor to adjust the vorticity to have
# zero net circulation

zeta = zeta - correction

# Expands the margins of a matplotlib axis, 
# and so prevents arrows on boundaries from being clipped. 
def stop_clipping(ax,marg=.02): # default is 2% increase
    l,r,b,t = ax.axis()
    dx,dy = r-l, t-b
    ax.axis([l-marg*dx, r+marg*dx, b-marg*dy, t+marg*dy])

# dqdt requires a list of the time derivatives for q, stored 
# in order from present to the past
def ab_blend(dqdt,order):
    if order==1:
        return dqdt[0]
    elif order==2:
        return 1.5*dqdt[0]-.5*dqdt[1]
    elif order==3:
        return (23*dqdt[0]-16*dqdt[1]+5*dqdt[2])/12.
    else:
        print("order", order ," not supported ")

def advect(q,u,v,dx,dy): 
# third-order upwind advection
# q,u,v are co-located    
    dqdt = np.zeros(q.shape)
    
    dqmx = np.zeros(q.shape)
    dqpx = np.zeros(q.shape)
    dqmy = np.zeros(q.shape)
    dqpy = np.zeros(q.shape)
    
    dqmx[:,1]  = -q[:,0] + q[:,1] # 1st order, plus side at left wall
    dqmx[:,2:-1] = (q[:,:-3] - 6*q[:,1:-2] + 3*q[:,2:-1] + 2*q[:,3:])/6. # 3rd order, minus side
    dqpx[:,-2] = -q[:,-2] + q[:,-1] # 1st order, plus side at right wall
    dqpx[:,1:-2] = (-2*q[:,0:-3] - 3*q[:,1:-2] + 6*q[:,2:-1] -1*q[:,3:])/6. #3rd order, plus side

    dqmy[1,:]  = -q[0,:] + q[1,:] # 1st order, minus side at bottom wall
    dqmy[2:-1,:] =  (q[:-3,:] - 6*q[1:-2,:] + 3*q[2:-1,:] + 2*q[3:,:])/6. # 3rd order, minus side
    dqpy[-2,:] = -q[-2,:] + q[-1,:] # 1st order, plus side at top wall
    dqpy[1:-2,:] = ( - 2*q[0:-3,:]  - 3*q[1:-2,:] + 6*q[2:-1,:] - q[3:,:] )/6. # 3rd order, plus side


    dqdx = np.where(u>0.,dqmx,dqpx)/dx # upwind, emphasize side from where fluid is coming from
    dqdy = np.where(v>0.,dqmy,dqpy)/dy # ditto
    
    dqdt += -u*dqdx
    dqdt += -v*dqdy
    
    return dqdt

#############################################################
def divergence(u,v,dx,dy):
    # du/dx + dv/dy at p-grid
    div = .5*( u[:-1,1:] + u[1:,1:] - u[:-1,:-1] - u[1:,:-1])/dx + \
          .5*( v[1:,:-1] + v[1:,1:] - v[:-1,:-1] - v[:-1,1:])/dy
    return div
#############################################################
def vortp(u,v,dx,dy):
    # dv/dx - du/dy at p-grid
    vort = .5*( v[:-1,1:] + v[1:,1:] - v[:-1,:-1] - v[1:,:-1])/dx - \
           .5*( u[1:,:-1] + u[1:,1:] - u[:-1,:-1] - u[:-1,1:])/dy
    return vort 
#############################################################
def vortU(u,v,dx,dy):
    # dv/dx - du/dy at U-grid interior points
    vort = np.zeros(u.shape)
    vort[1:-1,1:-1] =  (v[1:-1,2:] - v[1:-1,:-2])/(2*dx) - \
                       (u[2:,1:-1] - u[:-2,1:-1])/(2*dy)
    return vort

def psi_to_uv(q,dx,dy):
# q is streamfunction (psi) on u-grid, assumed to be 0 on boundaries
# returns v = dq/dx and u= -dq/dy, on U-grid
    u = 0.*q
    v = 0.*q
    
    u[1:-1,1:-1] = -( q[2:,1:-1] - q[:-2,1:-1] )/(2*dy)
    u[0,1:-1] = -q[1,1:-1]/dy
    u[-1,1:-1] = q[-2,1:-1]/dy
  
    v[1:-1,1:-1] = +( q[1:-1,2:] - q[1:-1,:-2])/(2*dx)
    v[1:-1,0] =  q[1:-1,1]/dx
    v[1:-1,-1] = -q[1:-1,-2]/dx
    return u,v

def laplacian(p,dx,dy,il=None, ir=None, jb=None, jt=None):
# Returns Laplacian of p, d^2p/dx^2 + d^2/dy^2.
# If needed, specify how to grab the image of a point outside
# the domain.  Otherwise, the d^2p/dx^2 or d^2/dy^2 term is not included
# on the boundary.  
    rdx2 = 1./(dx*dx)
    rdy2 = 1./(dy*dy)
    lapl = np.zeros(p.shape)
    lapl[:,1:-1]  =  rdx2*( p[:,:-2] -2*p[:,1:-1] + p[:,2:] )
    lapl[1:-1,:] +=  rdy2*( p[:-2,:] -2*p[1:-1,:] + p[2:,:] ) 
    if il in [-2,-1,0,1]:    
        lapl[:,0]  +=  rdx2*( p[:,il] -2*p[:,0] + p[:,1] ) 
    if ir in [-2,-1,0,1]:    
        lapl[:,-1] +=  rdx2*( p[:,-2] -2*p[:,-1] + p[:,ir] )
    if jb in [-2,-1,0,1]:
        lapl[0,:]  +=  rdy2*( p[jb,: ] -2*p[0,:] + p[1,:] ) 
    if jt in [-2,-1,0,1]:
        lapl[-1,:] +=  rdy2*( p[-2,: ] -2*p[-1,:] + p[jt,:] ) 
    return lapl

def poisson_fft_prep(Nx,Ny,dx,dy,lapl='discrete'):
    # returns the coefficients to multiply the vorticity Fourier amplitudes
    L = dx*(Nx-1)
    W = dy*(Ny-1)
    
    Ka = np.arange(Nx-2) +1 # integer wavenumbers of the sine functions in the x-direction
    Ma = np.arange(Ny-2) +1 # integer wavenumbers of the sine functions in the y-direction
    ka = Ka*np.pi/L
    ma = Ma*np.pi/W
    
    lapl_op = np.zeros( (Ny-2,Nx-2) )
    if lapl == 'discrete':
        lapl_op[:] += (2*np.cos(ka*dx)-2)/dx**2 # add to every row
    else: # the calculus Laplacian
        lapl_op[:] += -ka**2 
    lapl_opT = lapl_op.T # reverse columns and rows
    if lapl == 'discrete':
        lapl_opT[:] += (2*np.cos(ma*dy)-2)/dy**2 # add to every row
    else: # the calculus Laplacian
        lapl_opT[:] += -ma**2
    lapl_op = lapl_opT.T # reverse columns and rows
    invlapl = 1./lapl_op #the coefficents for multiplying the vorticity Fourier amplitudes
    return invlapl 

def poisson_fft(vort, invlapl): 
    # solves for psi in del^2 psi = vort 
    cv = vort[1:-1,1:-1] # central vorticity
    
    #convert gridded vorticity to gridded Fourier coefficients A_k,m
    cvt = scipy.fft.dst( cv , axis=1 , type=1) 
    cvt = scipy.fft.dst( cvt , axis=0 , type=1)  
   
    cpsit = cvt*invlapl # Calculate B_k,m from A_k,m
    
    # convert array of Fourier coefficents for psi to gridded central psi
    cpsit = scipy.fft.idst(cpsit,axis=0,type=1) # inverse transform 
    cpsi = scipy.fft.idst(cpsit,axis=1,type=1) # inverse transform 

    sh = vort.shape
    psi = np.zeros(sh) # we need 0 on boundaries, next line fills the center
    psi[1:-1,1:-1] = cpsi/(4*(sh[0]-1)*(sh[1]-1)) # apply normalization convention of FFT
    return psi

def jacobian(Nx,Ny,dx,dy,psi,zeta):
    jacob = np.zeros(Ny,Nx)
    #compute Jacobian for psi, zeta/dx, dy using Arakawa fourth order scheme
    ip1 = 0
    ip2 = 0
    im1 = 0
    im2 = 0
    jp1 = 0
    jp2 = 0
    jm1 = 0
    jm2 = 0
    dh = dx*dy
    weight = 1.0/(12.0 * dh)
    weight2 = 1.0/(24.0 * dh)
    for j in range(Ny):
        for i in range(Nx):
            #calculate indices for doubly periodic domain
            if (i == 0) or (i == Nx-1):
                ip1 = 1
                im1 = Nx-2
                ip2 = 2
                im2 = Nx-3
            elif (i == 1):
                ip1 = 2
                im1 = Nx-1
                ip2 = 3
                im2 = Nx-2
            elif (i == Nx-2):
                ip1 = 0
                im1 = Nx-3
                ip2 = 1
                im2 = Nx-4
            else:
                ip1 = i+1
                im1 = i-1
                ip2 = i+2
                im2 = i-2
            if (j == 0) or (j == Ny-1):
                jp1 = 1
                jm1 = Ny-2
                jp2 = 2
                jm2 = Ny-3
            elif (j == 1):
                jp1 = 2
                jm1 = Ny-1
                jp2 = 3
                jm2 = Ny-2
            elif (j == Ny-2):
                jp1 = 0
                jm1 = Ny-3
                jp2 = 1
                jm2 = Ny-4
            else:
                jp1 = j+1
                jm1 = j-1
                jp2 = j+2
                jm2 = j-2
            a1 = (psi[i,jm1] + psi[ip1,jm1] - psi[i,jp1] - psi[ip1,jp1]) * \
                ( zeta[ip1,j] - zeta[i,j])
            a2 = (psi[im1,jm1] + psi[i,jm1] - psi[im1,jp1] - psi[i,jp1]) * \
                ( zeta[i,j] - zeta[im1,j])
            a3 = (psi[ip1,j] + psi[ip1,jp1] - psi[im1,j] - psi[im1,jp1]) * \
                ( zeta[i,jp1] - zeta[i,j])
            a4 = (psi[ip1,jm1] + psi[ip1,j] - psi[im1,jm1] - psi[im1,j]) * \
                ( zeta[i,j] - zeta[i,jm1])
            b1 = (psi[ip1,j] - psi[i,jp1]) * (zeta[ip1,jp1] - zeta[i,j])
            b2 = (psi[i,jm1] - psi[im1,j]) * (zeta[i,j] - zeta[im1,jm1])
            b3 = (psi[i,jp1] - psi[im1,j]) * (zeta[im1,jp1] - zeta[i,j])
            b4 = (psi[ip1,j] - psi[i,jm1]) * (zeta[i,j] - zeta[ip1,jm1])
            vd1f = weight*(a1 + a2 + a3 + a4 + b1 + b2 + b3 + b4)
            
            # recompute coefficients for 4th order
            a1 = (zeta[ip1,jp1] - zeta[im1,jm1]) * (psi[im1,jp1] - psi[ip1,jm1])
            a2 = -1.0*(zeta[im1,jp1] - zeta[ip1,jm1]) * (psi[ip1,jp1] - psi[im1,jm1])
            a3 = zeta[ip1,jp1] * (psi[i,jp2] - psi[ip2,j])
            a4 = -1.0*zeta[im1,jm1] * (psi[im2,j] - psi[i,jm2])
            a5 = -1.0*zeta[im1,jp1] * (psi[i,jp2] - psi[im2,j])
            a6 = zeta[ip1,jm1] * (psi[ip2,j] - psi[i,jm2])
            a7 = zeta[ip2,j] * (psi[ip1,jp1] - psi[ip1,jm1])
            a8 = -1.0*zeta[im2,j] * (psi[im1,jp1] - psi[im1,jm1])
            a9 = -1.0*zeta[i,jp2] * (psi[ip1,jp1] - psi[im1,jp1])
            a10 = zeta[i,jm2] * (psi[ip1,jm1] - psi[im1,jm1])

            vd2f = -1.0*weight2*( a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10)
            jacob[i,j] = 2.0*vd1f - vd2f          

    return jacob

# 
xmax = 1. # 0 <= x <= xmax
ymax = 1.
x1U = np.linspace(0,xmax,n)
y1U = np.linspace(0,ymax,m)
x1p = .5*(x1U[:-1]+x1U[1:])
y1p = .5*(y1U[:-1]+y1U[1:])
xU,yU = np.meshgrid(x1U,y1U)
xp,yp = np.meshgrid(x1p,y1p)

print("compute array of inverse Laplacian")
# An array of the inverse Laplacian, 
# to be applied to the Fourier components of the r.h.s. of the Poisson equation.
# This is calculated once, and used throughout the notebook.
invlapl = poisson_fft_prep(n,m,dx,dy)#,lapl='discrete') #lapl='calculus' or lapl='discrete'

print("solve for initial streamfunction")
psii = poisson_fft(zeta,invlapl) # solve for initial psi
lapl_psii = laplacian(psii,dx,dy) # should recover vorti
ui,vi = psi_to_uv(psii,dx,dy) # solve for initial u and initial v from initial psi

# # check if poisson solver is working
# print("vort, maximum value:", zeta.max() )
# print("root mean square vort:", np.sqrt( (zeta*zeta).mean() ) )
# print("calculated psi, min value:", psii.min() )
# print("laplacian of psi, max value:", lapl_psii.max() )
# sqdiff = ( lapl_psii - zeta )**2
# rms = np.sqrt(sqdiff.mean())
# print("root mean square (lapl_psi - vort):",rms)
# print("psi -> u,v -> vortU, max value:", vortU(ui,vi,dx,dy).max())
# print("max u:",ui.max())
# print("max v:",vi.max())

vd=4 # vector skip in arrow plots (vd=1 plots all arrows)

vortamp = zeta.max()
vortlevs = np.linspace(-0.1*vortamp,1.1*vortamp,7)
print('contour levels for vorticity:',vortlevs)

psilevs=np.linspace(3*psii.min(),psii.max(),30) # factor of 3 accomodates stengthening vortex
print('contour levels for psi:',psilevs)

blevs = np.linspace(-.05,1.05,12)
print('contour levels for b:',blevs)

plt.contourf(xU[100:200,100:200],yU[100:200,100:200],zeta[100:200,100:200],vortlevs)
plt.savefig("./vorticity_init.jpg")
plt.contourf(xU[100:200,100:200],yU[100:200,100:200],psii[100:200,100:200],psilevs)
plt.savefig("./psi_inti.jpg")
#if abs(psii).max()>1.e-12: 
#    simple.contour(xU[100:200,100:200],yU[100:200,100:200],
#                   psii[100:200,100:200],colors=('w',) )
#    simple.quiver(xU[100:200:vd,100:200:vd],yU[100:200:vd,100:200:vd],
#                  ui[100:200:vd,100:200:vd],vi[100:200:vd,100:200:vd])
#simple.set_title("Initial $\zeta$, $\psi$ and $(u,v)$")
#stop_clipping(simple)


# need to type in timestep routine and then create loop for total all iterations    
