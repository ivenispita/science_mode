# -*- coding: utf-8 -*-
#
# philips_model_project.py
#
# purpose:  Compute and plot the linear and Non-linear terms of the QG-PV equation for the 2-layer system (i. e., the Phillips model)   
# author:   Tiago Bilo


# Scientific packages
import numpy as np 
import pylab as py
from scipy import interpolate
from scipy.optimize import curve_fit
import fortranfiles

# Ocean packages
import seawater as sw
import gsw

# Plotting packages
import matplotlib.pyplot as plt
from matplotlib import ticker, dates
from oceans_old.plotting import rstyle

# Paralelization packages
from joblib import Parallel, delayed  
import multiprocessing

# Warning off
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True


### Classes 
# Model Parameters 
class model_parameters(object):
	def __init__(self):


		self.U = np.array([0.06,0.])				# background velocities at each layer (m s-1)
		self.L = 1800000. 							# Meridional width of the channel (m)
		self.H = np.array([1000.,3000.])			# Thickness of each layer (m)
		self.f0 = 0.83e-4							# Coriolis parameter (s-1) - Channel centered at approximately 35N  
		self.beta = 2.0e-11							# Meridional variation of f0 (m-1 s-1)
		self.Rd = 25000. 							# First baroclinic Rossby deformation radius (m)
		
		self.delta_rho = ((25000.*0.83e-4)**2.)*1025./(9.8*(1000.*3000./(4000.))) # density jump (kg m-3)

		self.lateral_friction = 50. 				# Lateral Laplacian viscosity	(m2 s-1)				  
		self.bottom_friction = 1.0e-7 				# Bottom friction coeff. (s-1)

		self.g = 9.82 								# Acceleration due to gravity (m s-2)
		self.rho0 = 1025. 							# Average density of seawater (kg m-3)

# Scales of motion
class model_scales(object):
	def __init__(self):

		self.L0 = p.L/256. 							# Length scale (m)
		self.U0 = 0.01 								# Velocity scale (m s-1)
		self.T0 = (p.L/256.)/0.01 					# Advective time scale (s)


# Dimensionless energy
class load_energy(object):
	def __init__(self,filename):

		dtype = [("potential","double"),("kinetic_layer1","double"),("kinetic_layer2","double")]
		data = np.loadtxt(filename, dtype=dtype)

		self.potential = data['potential']		
		self.kinetic_layer1 = data['kinetic_layer1']
		self.kinetic_layer2 = data['kinetic_layer2']
		self.time = np.arange(10,1210,10)					# Time axis in days


# Dimensionless stream function at each layer
class load_psi(object):
	def __init__(self):
		num_cores = 6
		data = Parallel(n_jobs=num_cores)(delayed(load_fortranfile)(t) for t in xrange(10,1210,10))

		psi = np.ones((2,120,257,512))*np.nan

		for t in xrange(120):
			psi[0,t,:,:] = data[t][0] 
			psi[1,t,:,:] = data[t][1]

		self.psi = psi*s.U0*s.L0							# Perturbation stream function (m2 s-1)
		self.time = np.arange(10,1210,10)					# Time axis (days) 
		self.y = np.linspace(0.,p.L,257.)/1000. 			# meridional-axis (km)
		self.x = np.linspace(0.,s.L0*512,512.)/1000. 		# zonal-axis (km)


# getting a crest in around x = 1500 km 
class get_crest(object):
	def __init__(self,X,T,PSI):
		PSImax = 719.5

		X_I = X[np.abs(PSI[0,:]-PSImax).argmin()]
		x = np.ones(T.shape[0])*np.nan
		psi_crest = x.copy()

		x[0] = X_I
		psi_crest[0] = PSI[0,np.abs(PSI[0,:]-PSImax).argmin()]


		for i in xrange(1,x.shape[0],1):
			j = (X >= x[i-1]) & (X <= 1550)
			
			XX = X[j]
			PSII = PSI[i,j]

			jj = PSII == PSII.max()


			x[i] = XX[jj]
					
			PSImax = PSII[jj]
			psi_crest[i] = PSImax			


		self.x   = x
		self.psi = psi_crest


class compute_jacobian(object):
	def __init__(self,PSI):	

		global dx,dy,stf,q1

		R1 = p.delta_rho*p.g*p.H[0]/(p.rho0*p.f0*p.f0)
		Rterm = (PSI.psi[0,:,:,:]-PSI.psi[0,:,:,:])/(R1*R1)



		# Nice and precise way to compute derivatives (however it requires loops)
		dx = np.gradient(PSI.x)
		dy = np.gradient(PSI.y)


		stf = PSI.psi[0,:,:,:].copy()


		# Laplacian operator 
		num_cores = 6
		laplacian = Parallel(n_jobs=num_cores)(delayed(hlaplacian)(t) for t in xrange(PSI.time.shape[0]))
		laplacian = np.array(laplacian)

		q1 = laplacian-Rterm


		# Jacobian operator J(psi1,q1)
		num_cores = 6
		jacobian = Parallel(n_jobs=num_cores)(delayed(jacobian_layer1)(t) for t in xrange(PSI.time.shape[0]))	
		jacobian = np.array(jacobian)


		self.q1 = q1.copy()
		self.jacobian = jacobian



		del dx,dy,stf,q1



### Subroutines 
def load_fortranfile(T):

	if T < 100:
		file1 = '../outputs/psi1/psi1_0'+np.str(T)+'.bin'
		file2 = '../outputs/psi2/psi2_0'+np.str(T)+'.bin'

	else:
		file1 = '../outputs/psi1/psi1_'+np.str(T)+'.bin'
		file2 = '../outputs/psi2/psi2_'+np.str(T)+'.bin'		


	PSI1 = fortranfiles.FortranFile(file1) 
	PSI1 = PSI1.readReals()


	PSI2 = fortranfiles.FortranFile(file2) 
	PSI2 = PSI2.readReals()



	# Redimensionalizing files
	PSI1 = PSI1.reshape((257,512))
	PSI2 = PSI2.reshape((257,512))



	return PSI1,PSI2


def plt_energy(fign):

	fig = plt.figure(num=fign,figsize=(8,5),facecolor='w')
	ax = plt.gca()


	ax.plot(energy.time,energy.potential,'k',lw=3.5,label='EPE')
	ax.plot(energy.time,energy.kinetic_layer1,'r',lw=3.5,label='Layer 1 - EKE')
	ax.plot(energy.time,energy.kinetic_layer2,'b',lw=3.5,label='Layer 2 - EKE')


	ax.set_ylim(0,20)
	ax.set_xlim(0,1210)
	ax.set_yticks(np.arange(0,22.5,2.5))
	ax.tick_params(labelsize=22)
	rstyle(ax)

	ax.set_ylabel('Nondimensionalized energy',fontsize=25,fontweight='demibold')
	ax.set_xlabel(ur'Time (days)',fontsize=30,fontweight='demibold')

	ax.legend(numpoints=1,loc=0,prop={'size':18},frameon=False)

	# plt.draw()
	# plt.show(block=False)

	return fig 


def plt_psi(LAYER):

	if LAYER == 'Layer 1':
		levels = np.arange(-18.,18.2,0.2)
		ticks = np.arange(-18,20,2)		
		PSI = psi.psi[0,:,:,:].copy()
	else:
		levels = np.arange(-9.,9.1,0.1)
		ticks = np.arange(-9,10.,1.)		
		PSI = psi.psi[1,:,:,:].copy()


	# Snapshots to be plotted
	T = [10,50,100,150,250,350,550,800,1200]


	fig, axes = plt.subplots(figsize=(11,7.3),nrows=3,ncols=3,facecolor='w')


	t = 0

	for i in xrange(3):
		for j in xrange(3):

			k = psi.time == T[t]

			cm = axes[i,j].contourf(psi.x,psi.y,PSI[k,:,:].squeeze()/1000.,levels=levels,cmap='seismic',extend='both')


			axes[i,j].set_xlim(psi.x.min(),psi.x.max())
			axes[i,j].set_ylim(psi.y.min(),psi.y.max())
			axes[i,j].tick_params(labelsize=15)

			axes[i,j].set_title(np.str(T[t])+' days',fontsize=22,fontweight='demibold')

			if i < 2:
				axes[i,j].set_xticks([])
			else:
				axes[i,j].set_xticks(range(0,4800,1200))


			if j > 0:
				axes[i,j].set_yticks([])
			else:
				axes[i,j].set_yticks(range(0,2100,300))				


			t += 1


	axes[1,0].set_ylabel('y (km)',fontsize=30,fontweight='demibold')
	axes[2,1].set_xlabel('x (km)',fontsize=30,fontweight='demibold')
	
	fig.suptitle(LAYER,fontsize=32,fontweight='demibold')
				
	cbaxes = fig.add_axes([.92, 0.1, 0.03, 0.8])
	cb = fig.colorbar(cm,cax=cbaxes,extend='both',ticks=ticks) 
	cb.ax.tick_params(labelsize=20)
	cb.set_label(ur"$\psi \times 10^{3}$ (m$^{2}$ s$^{-1}$)",fontsize=25,fontweight='bold')

	# plt.draw()
	# plt.show(block=False)

	return fig 


def plt_jacobian(fign):

	levels = np.arange(-15.,15.25,0.25)
	ticks = np.arange(-15,17.5,2.5)		


	# Snapshots to be plotted
	T = [10,50,100,150,250,350,550,800,1200]


	fig, axes = plt.subplots(figsize=(11,7.3),nrows=3,ncols=3,facecolor='w')


	t = 0

	for i in xrange(3):
		for j in xrange(3):

			k = psi.time == T[t]

			cm = axes[i,j].contourf(psi.x,psi.y,pv.jacobian[k,:,:].squeeze(),levels=levels,cmap='seismic',extend='both')

			axes[i,j].set_xlim(psi.x.min(),psi.x.max())
			axes[i,j].set_ylim(psi.y.min(),psi.y.max())
			axes[i,j].tick_params(labelsize=15)

			axes[i,j].set_title(np.str(T[t])+' days',fontsize=22,fontweight='demibold')

			if i < 2:
				axes[i,j].set_xticks([])
			else:
				axes[i,j].set_xticks(range(0,4800,1200))


			if j > 0:
				axes[i,j].set_yticks([])
			else:
				axes[i,j].set_yticks(range(0,2100,300))				


			t += 1


	axes[1,0].set_ylabel('y (km)',fontsize=30,fontweight='demibold')
	axes[2,1].set_xlabel('x (km)',fontsize=30,fontweight='demibold')
	
	cbaxes = fig.add_axes([.92, 0.1, 0.03, 0.8])
	cb = fig.colorbar(cm,cax=cbaxes,extend='both',ticks=ticks) 
	cb.ax.tick_params(labelsize=20)
	cb.set_label(ur"Layer 1 - J(q,$\psi$) $= \frac{\partial \psi}{\partial x} \frac{\partial q}{\partial y} - \frac{\partial \psi}{\partial y} \frac{\partial q}{\partial x}$ (s$^{-2}$)"
		,fontsize=22,fontweight='bold')

	# plt.draw()
	# plt.show(block=False)

	return fig 


def plt_transition(fign):

	# Snapshots to be plotted
	T = range(150,270,20)

	fig, axes = plt.subplots(figsize=(11,5),nrows=2,ncols=3,facecolor='w')

	t = 0

	for i in xrange(2):
		for j in xrange(3):

			k = psi.time == T[t]

			cm = axes[i,j].contourf(psi.x,psi.y,psi.psi[0,k,:,:].squeeze()/1000.,
				levels=np.arange(-18.,18.2,0.2),cmap='seismic',extend='both')


			axes[i,j].set_xlim(psi.x.min(),psi.x.max())
			axes[i,j].set_ylim(psi.y.min(),psi.y.max())
			axes[i,j].tick_params(labelsize=17)

			axes[i,j].set_title(np.str(T[t])+' days',fontsize=22,fontweight='demibold')

			if i < 1:
				axes[i,j].set_xticks([])
			else:
				axes[i,j].set_xticks(range(0,4800,1200))


			if j > 0:
				axes[i,j].set_yticks([])
			else:
				axes[i,j].set_yticks(range(0,2100,300))				


			t += 1

	axes[1,0].set_ylabel('y (km)',fontsize=30,fontweight='demibold',y=1.02)
	axes[1,1].set_xlabel('x (km)',fontsize=30,fontweight='demibold')
	
				
	cbaxes = fig.add_axes([.92, 0.1, 0.03, 0.8])
	cb = fig.colorbar(cm,cax=cbaxes,extend='both',
		ticks=np.arange(-18,20,2)) 
	cb.ax.tick_params(labelsize=20)
	cb.set_label(ur"Layer 1 - $\psi \times 10^{3}$ (m$^{2}$ s$^{-1}$)",fontsize=25,fontweight='bold')

	# plt.draw()
	# plt.show(block=False)

	return fig



def plt_jacobian_transition(fign):

	# Snapshots to be plotted
	T = range(150,270,20)

	fig, axes = plt.subplots(figsize=(11,5),nrows=2,ncols=3,facecolor='w')

	t = 0

	for i in xrange(2):
		for j in xrange(3):

			k = psi.time == T[t]

			cm = axes[i,j].contourf(psi.x,psi.y,pv.jacobian[k,:,:].squeeze(),
				levels=np.arange(-15.,15.25,0.25),cmap='seismic',extend='both')


			axes[i,j].set_xlim(psi.x.min(),psi.x.max())
			axes[i,j].set_ylim(psi.y.min(),psi.y.max())
			axes[i,j].tick_params(labelsize=17)

			axes[i,j].set_title(np.str(T[t])+' days',fontsize=22,fontweight='demibold')

			if i < 1:
				axes[i,j].set_xticks([])
			else:
				axes[i,j].set_xticks(range(0,4800,1200))


			if j > 0:
				axes[i,j].set_yticks([])
			else:
				axes[i,j].set_yticks(range(0,2100,300))				


			t += 1

	axes[1,0].set_ylabel('y (km)',fontsize=30,fontweight='demibold',y=1.02)
	axes[1,1].set_xlabel('x (km)',fontsize=30,fontweight='demibold')
	
				
	cbaxes = fig.add_axes([.92, 0.1, 0.03, 0.8])
	cb = fig.colorbar(cm,cax=cbaxes,extend='both',
		ticks=np.arange(-15,17.5,2.5)) 
	cb.ax.tick_params(labelsize=20)
	cb.set_label(ur"Layer 1 - J(q,$\psi$) $= \frac{\partial \psi}{\partial x} \frac{\partial q}{\partial y} - \frac{\partial \psi}{\partial y} \frac{\partial q}{\partial x}$ (s$^{-2}$)"
		,fontsize=22,fontweight='bold')

	# plt.draw()
	# plt.show(block=False)

	return fig





def plt_hovmoller(fign):

	t = psi.time<=200
	j = py.find(psi.y == 900)[0]

	T = psi.time[t].copy()
	PSI = psi.psi[0,t,j,:].squeeze().copy()/1000.


	# Computing the phase speed
	crest = get_crest(psi.x,T,PSI*1000)

	c = ((crest.x[-1]-crest.x[0])*1000)/((T[-1]-T[0])*24.*60.*60.)

	fig = plt.figure(num=fign,figsize=(7,6),facecolor='w')
	ax = plt.gca()

	cm = ax.contourf(psi.x,T,PSI,
		levels=np.arange(-9.,9.1,0.1),cmap='seismic',extend='both')

	ax.plot(crest.x,T,'k--',lw=5.)
	ax.set_title(ur'Phase Speed $\sim$ %1.2f $\times 10^{-2}$ m s$^{-1}$'%(c*100.),fontsize=25,fontweight='demibold')


	ax.set_xlim(psi.x.min(),psi.x.max())
	ax.set_ylim(10.,200)
	ax.set_xticks(range(0,4200,600))
	ax.set_yticks(range(10,200,15))

	ax.tick_params(labelsize=22)


	ax.set_ylabel('Time (days)',fontsize=30,fontweight='demibold')
	ax.set_xlabel('x (km)',fontsize=30,fontweight='demibold')

	cbaxes = fig.add_axes([.93, 0.1, 0.04, 0.8])
	cb = fig.colorbar(cm,cax=cbaxes,extend='both',
		ticks=np.arange(-9,10,1)) 
	cb.ax.tick_params(labelsize=20)
	cb.set_label(ur"Layer 1 - $\psi \times 10^{3}$ (m$^{2}$ s$^{-1}$)",fontsize=25,fontweight='bold')

	# plt.draw()
	# plt.show(block=False)


	return fig, crest 


def plt_gr(fign):

	t = psi.time<=200

	T = psi.time[t].copy()	

	# Fitting the growth curve
	popt, pcov = curve_fit(exponential, T, crest_growth.psi, p0=(1, 1e-6, 1))
	curve = exponential(T, *popt)


	fig = plt.figure(num=fign,figsize=(8,5),facecolor='w')
	ax = plt.gca()


	ax.plot(T,crest_growth.psi/1000,'bo',mfc='b',markersize=10,label=ur'Crest $\psi$')
	ax.plot(T,curve/1000,'k--',lw=4.5,label=ur'Exponential Fit')

	ax.set_xlim(10.,200)
	ax.set_xticks(range(10,200,15))	
	ax.tick_params(labelsize=22)
	rstyle(ax)

	ax.set_ylabel(ur"$\psi(y=900) \times 10^{3}$ (m$^{2}$ s$^{-1}$)",fontsize=25,fontweight='demibold')
	ax.set_xlabel(ur'Time (days)',fontsize=30,fontweight='demibold')
	ax.set_title(ur'Growth Rate $\sim$ %1.2f $\times 10^{-2}$ day$^{-1}$'%(popt[1]*100.),fontsize=25,fontweight='demibold')

	ax.legend(numpoints=1,loc=0,prop={'size':25},frameon=False)

	plt.draw()
	plt.show(block=False)

	return fig 


def exponential(t, a, b, d):

	fit_psi = a*np.exp(b*t)+d

	return fit_psi


# Horizontal Laplacian Operator
def hlaplacian(T):
	STF = stf[T,:,:].copy()

	d2pdx2 = np.ones(STF.shape)*np.nan
	d2pdy2 = d2pdx2.copy()

	for i in xrange(STF.shape[0]):
		d2pdx2[i,:] = np.gradient(STF[i,:],dx)
		d2pdx2[i,:] = np.gradient(d2pdx2[i,:],dx)	

	for j in xrange(STF.shape[1]):
		d2pdy2[:,j] = np.gradient(STF[:,j],dy)
		d2pdy2[:,j] = np.gradient(d2pdy2[:,j],dy)


	LAPLACIAN = d2pdx2 + d2pdy2

	return LAPLACIAN


# Jacobian operator
def jacobian_layer1(T):
	STF = stf[T,:,:].copy()
	Q = q1[T,:,:].copy()


	dpdx = Q.copy()*np.nan
	dqdy = dpdx.copy()

	dpdy = dpdx.copy()
	dqdx = dpdx.copy()


	for i in xrange(STF.shape[0]):
		dpdx[i,:] = np.gradient(STF[i,:],dx)
		dqdx[i,:] = np.gradient(Q[i,:],dx)

	for j in xrange(STF.shape[1]):
		dpdy[:,j] = np.gradient(STF[:,j],dy)
		dqdy[:,j] = np.gradient(Q[:,j],dy)

	JACOBIAN = (dpdx*dqdy)-(dqdx*dpdy)

	return JACOBIAN 


### Main Program
## Parameters and scales 
p = model_parameters()
s = model_scales()
XXXXXX

## Loading the model outputs (perturbation fields)
psi = load_psi()
energy = load_energy('../outputs/eddy_energy.dat') 



## Computing the PV-equation non-linear jacobian
pv = compute_jacobian(psi)



# Plotting results: General stuff
# Eddy Energy development 
fig1 = plt_energy(1)

# Temporal development of the perturbations
fig2 = plt_psi('Layer 1')
fig3 = plt_psi('Layer 2')

# Non-linearity becoming important
fig4 = plt_transition(4)


# Plotting results: Terms and linear instability properties

# Hovmoller of the first 200 days
fig5,crest_growth = plt_hovmoller(5)

# Growth rate
fig6 = plt_gr(6)


# Jacobian
fig7 = plt_jacobian(7)
fig8 = plt_jacobian_transition(8)



# Saving figures
fig1.savefig('../figures/energy.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig2.savefig('../figures/psi_layer1.png',format='png',bbox_inches='tight',
	transparent=True,pad_inches=0.5)

fig3.savefig('../figures/psi_layer2.png',format='png',bbox_inches='tight',
	transparent=True,pad_inches=0.5)

fig4.savefig('../figures/psi_transition_layer1.png',format='png',bbox_inches='tight',
	transparent=True,pad_inches=0.5)

fig5.savefig('../figures/psi_layer1_y900_hovmoller.png',format='png',bbox_inches='tight',
	transparent=True,pad_inches=0.5)

fig6.savefig('../figures/psi_layer1_y900_growth.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig7.savefig('../figures/jacobian_layer1.png',format='png',bbox_inches='tight',
	transparent=True,pad_inches=0.5)

fig8.savefig('../figures/jacobian_transition_layer1.png',format='png',bbox_inches='tight',
	transparent=True,pad_inches=0.5)
