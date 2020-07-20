# -*- coding: utf-8 -*-
#
# hmw4_q2.py
#
# purpose:  Plot the Eady Model solutions from question 2
#
# author:   Tiago Bilo

import numpy as np 
import matplotlib.pyplot as plt
from oceans_old.plotting import rstyle
import gsw

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True





## Parameters 
f = gsw.f(45)					# Coriolis parameter (s^-1)
H = 10.*1000					# Height of the domain (m)
u0 = 10. 						# Maximum velocity (m s^-1)
N = 0.01						# Brunt-Vaisala frequency (s^-1)
Rd = N*H/f 						# Rossby deformation radius (m)

## Domain axes (m) 
z = np.linspace(0,H,100)
x = np.linspace(-10/1.0e-6,10/1.0e-6,100)

## Zonal wave numbers (m^-1)
k = np.linspace(0,4.0e-6,100)

## Admensional wavenumber mu
mu = N*H*k/f


## Growth rate day^{-1}
gr = (k*u0/mu)*np.sqrt(((1./np.tanh(mu/2.))-(mu/2))*((mu/2.)-np.tanh(mu/2.)))
gr[np.isnan(gr)] = 0.

agr = gr.copy()*(Rd/u0)
gr =gr*60.*60.*24.


# fig1 = plt.figure(facecolor='w')
# ax1 = plt.gca()

# ax1.plot(mu,agr,color='darkkhaki',lw=3.5)
# ax1.plot(mu[agr == agr.max()],agr.max(),marker='o',mfc='k',
# 	markeredgecolor='darkkhaki',markersize=15)

# ax1.text(mu[agr == agr.max()]-0.6,agr.max()+0.02,
# 	ur'($\mu = $%1.2f, $\sigma_{i} = $ %1.2f)'%(mu[agr == agr.max()],agr.max()),
# 	fontsize=20,fontweight='demibold')

# ax1.text(mu[agr == agr.max()]-1.2,agr.max()+0.05,
# 	ur'($k =$ %1.2f $\times$ 10$^{6}$ m$^{-1}$, $\omega_{i} = $%1.2f day$^{-1}$)'%(k[agr == agr.max()]*1.0e6,gr.max()),
# 	fontsize=20,fontweight='demibold')


# ax1.set_xlim(0,4)
# ax1.set_ylim(0,0.4)

# rstyle(ax1)

# ax1.set_xlabel(ur"Dimensionless wavenumber $\mu = \frac{N H}{f_{0}} k$",
# 	fontsize=22,fontweight='demibold')

# ax1.set_ylabel(ur"Dimensionless growth rate $\sigma_{i} = \frac{N H}{f_{0} U_{0}} \omega_{i}$",
# 	fontsize=22,fontweight='demibold')

# ax1.tick_params(labelsize=20)


# fig1.savefig('q2_i.png',format='png',bbox_inches='tight',pad_inches=0.5)



## Streamfunction at t = 0 (m^2 s^-1)
z = z+ 1j*np.zeros(z.shape[0])

c = (u0/2.)+(u0/1.61)*np.sqrt((-((1.+0j)/np.tanh(1.61/2.))+(1.61/2))*((1.61/2.)-(1.+0j)*np.tanh(1.61/2.)))			

a = 1.
b = (-a)*u0/(1.61*c)

phi = (a*np.cosh(1.61*z/H))+(b*np.sinh(1.61*z/H))
phix = (np.cos(k[agr == agr.max()]*x) + 1j*np.sin(k[agr == agr.max()]*x))

phix,phi = np.meshgrid(phi,phix)

psi = (phi*phix).real

psip = psi.copy(); psip[psi < 0] = np.nan
psin = psi.copy(); psin[psi > 0] = np.nan

fig2 = plt.figure(facecolor='w',figsize=(12,4))
ax2 = plt.gca()

ax2.contour(x/1000.,z.real/1000,psip.transpose(),5,colors='k',linestyles='solid',linewidth=3.5)
ax2.contour(x/1000.,z.real/1000,psin.transpose(),5,colors='k',linestyles='dashed',linewidth=4.0)

ax2.set_title(ur"Vertical and horizontal structures of the $\psi$ ($k = $ %1.2f $\times$ 10$^{6}$ m$^{-1}$)"%(k[agr == agr.max()]*1.0e6),
	fontsize=20,fontweight='demibold')
ax2.set_ylabel('Height (km)',fontsize=25,fontweight='demibold')
ax2.set_xlabel('Zonal distance (km)',fontsize=25,fontweight='demibold')

ax2.tick_params(labelsize=20)


fig2.savefig('q2_ii.png',format='png',bbox_inches='tight',pad_inches=0.5)
