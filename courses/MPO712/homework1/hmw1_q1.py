# -*- coding: utf-8 -*-
# MPO 712 - Large Scale Ocean Circulation
# Homework 1 - Question 1
#
# Tiago Bilo - Spring 2016 

import numpy as np 
import matplotlib.pyplot as plt

plt.close('all')

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True

### Question 1
## (a) Wind stress field 

# Parameters
t0 = 0.1 				# N m-2 
L = 1000.0*1000.0 		# m 
f0 = 1.E-4				# s-1
b0 = 1.E-11				# m-1 s-1
rho0 = 1025 			# kg m-3

# Domain of L x 2L m 
dx = 50000.
x = np.arange(0,(1000.**2.)+dx,dx)
y = np.arange(0,(2.*(1000.**2.))+dx,dx)

x,y = np.meshgrid(x,y)


# Wind Stress field  
taux = -t0*np.cos(np.pi*y/L) 
tauy = t0*np.cos(np.pi*x/L)*np.sin(np.pi*y/L)


fig1 = plt.figure(1,figsize=(10,10),facecolor='w')
ax1 = plt.gca()

ax1.quiver(x/1000.,y/1000.,taux,tauy)
ax1.tick_params(labelsize=25)
ax1.grid('on')

ax1.set_xlabel('West-East distance x [km]',fontsize=30,fontweight='bold')
ax1.set_ylabel('South-North distance y [km]',fontsize=30,fontweight='bold')
ax1.set_title(ur'Wind Stress $\tau (x,y)$', fontsize=35,fontweight='bold')



# Sverdrup Transport My
My = (-t0/(L*b0))*np.sin(np.pi*y/L)*(1.+np.sin(np.pi*x/L))

fig2 = plt.figure(2,figsize=(10,10),facecolor='w')
ax2 = plt.gca()

cm2 = ax2.contourf(x/1000.,y/1000.,My,200,extend='both',cmap='seismic')
cb2 = fig2.colorbar(cm2,ax=ax2,fraction=0.045,extend='both')
cb2.ax.tick_params(labelsize=20)
ax2.tick_params(labelsize=25)
ax2.grid('on')

ax2.set_xlabel('West-East distance x [km]',fontsize=30,fontweight='bold')
ax2.set_ylabel('South-North distance y [km]',fontsize=30,fontweight='bold')
fig2.suptitle(ur'Sverdrup Transport $\times \rho _{0}$ [kg $s^{-1}$ $m^{-1}$]',fontsize=20,fontweight='bold')
ax2.set_title(ur'$\tau _{0} = 0$ Pa, $\beta _{0} = 10^{-11}$ m$^{-1}$ s$^{-1}$',
	fontsize=20,fontweight='bold')


# Streamfunction 
C = (-t0/b0)*(1.+np.pi)*np.sin(np.pi*y/L)
psi = ((-t0*np.pi/(b0*L))*np.sin(np.pi*y/L)*((np.cos(np.pi*x/L)*(L/np.pi))-x)) + C

psi_n = psi.copy(); psi_n[psi>=0.] = np.nan
psi_p = psi.copy(); psi_p[psi<0.] = np.nan

fig3 = plt.figure(3,figsize=(10,10),facecolor='w')
ax3 = plt.gca()

c3 = ax3.contour(x/1000.,y/1000.,psi_n/(rho0*1.E+6),10,linestyles='solid',colors='k',linewidths=2.5)
c31 = ax3.contour(x/1000.,y/1000.,psi_p/(rho0*1.E+6),10,linestyles='dashed',colors='k',linewidths=3.5)
ax3.clabel(c3,fmt='%i',colors='k',fontsize=20,fontweight='demibold',inline=1)
ax3.clabel(c31,fmt='%i',colors='k',fontsize=20,fontweight='demibold',inline=1)
ax3.tick_params(labelsize=25)
ax3.grid('on')

ax3.set_xlabel('West-East distance x [km]',fontsize=30,fontweight='bold')
ax3.set_ylabel('South-North distance y [km]',fontsize=30,fontweight='bold')
fig3.suptitle(ur'Sverdrup Transport Streamfunction [Sv]',fontsize=20,fontweight='bold')
ax3.set_title(ur'$\tau _{0} = 0$ Pa, $\beta _{0} = 10^{-11}$ m$^{-1}$ s$^{-1}$, $\rho _{0}$ = 1025 kg m$^{-3}$',
	fontsize=20,fontweight='bold')



# saving figures 
fig1.savefig('wind_field.png',format='png',bbox_inches='tight',pad_inches=0.5)
fig2.savefig('sverdrup_transport.png',format='png',bbox_inches='tight',pad_inches=0.5)
fig3.savefig('sverdrup_streamfuncion.png',format='png',bbox_inches='tight',pad_inches=0.5)
