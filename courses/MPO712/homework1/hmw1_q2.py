# -*- coding: utf-8 -*-
# MPO 712 - Large Scale Ocean Circulation
# Homework 1 - Question 2
#
# Tiago Bilo - Spring 2016 


import numpy as np 

## My Parameters
t0 = 0.1 				# N m-2 
L = 5.*1000.0*1000.0 	# m 
f0 = 1.E-4				# s-1
beta = 1.E-11			# m-1 s-1
rho0 = 1025.0 			# kg m-3
D = 200.0				# m 
b = 3.*1000.0*1000.0	# m 

Vmax = 2. 				# m s-1
dVmax = 0.05*2. 		# m s-1 (5% of Vmax)



# Stommel's friction term 
R = t0*np.pi*beta*L/((Vmax*b*rho0*D*beta)+(t0*np.pi))

Rplus = t0*np.pi*beta*L/(((Vmax+dVmax)*b*rho0*D*beta)+(t0*np.pi)) 
Rminus = t0*np.pi*beta*L/(((Vmax-dVmax)*b*rho0*D*beta)+(t0*np.pi))


# Munk's friction term 
C = Vmax*np.sqrt(3.0)*beta*b*D*rho0/(2.*t0*np.pi*np.exp(-np.pi/np.sqrt(3.)))
Cplus = (Vmax+dVmax)*np.sqrt(3.0)*beta*b*D*rho0/(2.*t0*np.pi*np.exp(-np.pi/np.sqrt(3.)))
Cminus = (Vmax-dVmax)*np.sqrt(3.0)*beta*b*D*rho0/(2.*t0*np.pi*np.exp(-np.pi/np.sqrt(3.)))


A = beta*((L/(C+1.))**3.)
Aplus = beta*((L/(Cplus+1.))**3.)
Aminus = beta*((L/(Cminus+1.))**3.)