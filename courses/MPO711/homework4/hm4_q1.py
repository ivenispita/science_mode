# -*- coding: utf-8 -*-
#
# hmw4_q1.py
#
# purpose:  Plot different profiles for the meridional potential vorticity gradient from question 1
#
# author:   Tiago Bilo

import numpy as np 
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True

## Parameters
r = 6.371e+6 														# Earth's radius (m)
omega = 7.292115e-5 												# Earth's rotation rate (s^-1)
beta = 2.*omega*np.cos(np.pi*45./180.)/r 							# m^-1 s^-1

L1 = 1.0e+6															# 2L = meridional width of the domain (m)
L2 = 1.0e+5



## "Amplitude" of the velocity u0
u01 = np.arange(-25,26,1)
u02 = np.arange(-2.5,2.6,0.1)

## Meridional axis y
y1 = np.arange(-L1,L1+5.0e+3,5.0e+3)
y2 = np.arange(-L2,L2+5.0e+2,5.0e+2)


## Potential gradient dqdy
dqdyA1 = np.ones((y1.shape[0],u01.shape[0]))*np.nan
dqdyA2 = np.ones((y2.shape[0],u02.shape[0]))*np.nan
dqdyB1 = dqdyA1.copy()
dqdyB2 = dqdyA2.copy()


for i in xrange(u01.shape[0]):
	dqdyA1[:,i] = beta + (6.0*u01[i]*y1/(L1*L1*L1))
	dqdyA2[:,i] = beta + (6.0*u02[i]*y2/(L2*L2*L2))	

	dqdyB1[:,i] = beta + (u01[i]*np.pi*np.pi/(4.0*L1*L1))*np.cos(np.pi*y1/(2.0*L1))
	dqdyB2[:,i] = beta + (u02[i]*np.pi*np.pi/(4.0*L2*L2))*np.cos(np.pi*y2/(2.0*L2))




## Item (a)
fig1, axes1 = plt.subplots(figsize=(4,8),nrows=2,ncols=1,facecolor='w')

axes1[0].contourf(u01,y1/1000.,dqdyA1,
	np.linspace(-np.abs(dqdyA1).max(),np.abs(dqdyA1).max(),num=200),cmap='seismic')

axes1[0].contour(u01,y1/1000,dqdyA1,levels=[0.],colors='k',linewidths=3.5)

axes1[0].set_yticks(range(-1000,1250,250))
axes1[0].set_xticks(range(-25,30,5))
axes1[0].tick_params(labelsize=13)

axes1[0].set_title(ur'$\frac{\partial \Pi}{\partial y} = \beta_{0} + \frac{6 U_{0}}{L^{3}} \times y$',fontsize=20)
axes1[0].set_ylabel(ur'y (km)',fontsize=20,fontweight='demibold')



axes1[1].contourf(u02,y2/1000.,dqdyA2,
	np.linspace(-np.abs(dqdyA2).max(),np.abs(dqdyA2).max(),num=200),cmap='seismic')

axes1[1].contour(u02,y2/1000,dqdyA2,levels=[0.],colors='k',linewidths=3.5)

axes1[1].set_yticks(range(-100,125,25))
axes1[1].set_xticks(np.arange(-2.5,3.0,0.5))
axes1[1].tick_params(labelsize=13)



axes1[1].set_xlabel(ur'U$_{0}$ (m s$^{-1}$)',fontsize=20,fontweight='demibold')
axes1[1].set_ylabel(ur'y (km)',fontsize=20,fontweight='demibold')


fig1.savefig('itemA.png',format='png',bbox_inches='tight',pad_inches=0.5)



plt.show(block=False)


## Item (b)
fig2, axes2 = plt.subplots(figsize=(4,8),nrows=2,ncols=1,facecolor='w')

axes2[0].contourf(u01,y1/1000.,dqdyB1,
	np.linspace(-np.abs(dqdyB1).max(),np.abs(dqdyB1).max(),num=200),cmap='seismic')

axes2[0].contour(u01,y1/1000,dqdyB1,levels=[0.],colors='k',linewidths=3.5)

axes2[0].set_yticks(range(-1000,1250,250))
axes2[0].set_xticks(range(-25,30,5))
axes2[0].tick_params(labelsize=13)

axes2[0].set_title(ur'$\frac{\partial \Pi}{\partial y} = \beta_{0} + \frac{U_{0} \pi^{2}}{4 L^{2}} \times cos(\frac{\pi}{2 L} y)$'
	,fontsize=20)
axes2[0].set_ylabel(ur'y (km)',fontsize=20,fontweight='demibold')



axes2[1].contourf(u02,y2/1000.,dqdyB2,
	np.linspace(-np.abs(dqdyB2).max(),np.abs(dqdyB2).max(),num=200),cmap='seismic')

axes2[1].contour(u02,y2/1000,dqdyB2,levels=[0.],colors='k',linewidths=3.5)

axes2[1].set_yticks(range(-100,125,25))
axes2[1].set_xticks(np.arange(-2.5,3.0,0.5))
axes2[1].tick_params(labelsize=13)



axes2[1].set_xlabel(ur'U$_{0}$ (m s$^{-1}$)',fontsize=20,fontweight='demibold')
axes2[1].set_ylabel(ur'y (km)',fontsize=20,fontweight='demibold')


fig2.savefig('itemB.png',format='png',bbox_inches='tight',pad_inches=0.5)


plt.show(block=False)

