# -*- coding: utf-8 -*-
#
# compute_geostrophic_velocity.py
#
# purpose:  Compute geostrophic velocities from the SOEST ADT climatology (from ARGO data)
#
# author:   Tiago Bilo

import numpy as np
import pylab as py
import seawater as sw
import models 
import scipy.io

# Dates manipulation
from dateutil.parser import parse
from netCDF4 import num2date
import datetime 

# Paralelization of the code 
from joblib import Parallel, delayed  
import multiprocessing



import warnings
warnings.filterwarnings('ignore')




### Classes 
class constants(object):
	def __init__(self):
		# Beta = meridional rate of change of the Coriolis parameters (m^-1 s^-1)
		r = 6.371e+6 														# Earth's radius (m)
		omega = 7.292115e-5 												# Earth's rotation rate (s^-1)
		beta = 2.*omega*np.cos(np.pi*soest.LAT16_100/180.)/r 

		# Coriolis parameter (s^-1)
		f = sw.f(soest.LAT16_100)

		# Rho = typical density of the sea water above the pycnocline (kg m^-3)
		# Since we are considering a incompressible ocean (i. e., Boussinesq approximation) it will be 
		# constant
		rho = 1025.0

		self.beta = beta 
		self.rho  = rho 
		self.f = f



class compute_geos(object):
	def __init__(self):

		global longitude, v, psi

		# Months
		time = num2date(soest.TIME.copy(),'days since 0001-01-01 00:00:00')


		# Merional geostrophic velocities
		data = Parallel(n_jobs=num_cores)(delayed(dm)(t) for t in xrange(time.shape[0]))

		latitude = soest.LAT16_100.copy()
		longitude = data[0][0]

		v = np.ones((time.shape[0],soest.LEV.shape[0],latitude.shape[0],longitude.shape[0]))*np.nan

		months = []
		for t in xrange(time.shape[0]):
			time1 = parse(np.str(time[t]),dayfirst=True)
			months.append(datetime.datetime.strftime(time1,'%B'))

			v[t,:,:,:] = data[t][1]

		del data 


		# Geostrophic velocities stream function at each depth
		data = Parallel(n_jobs=num_cores)(delayed(cpsi)(i) for i in xrange(latitude.shape[0]))

		psi = np.array(data).transpose((1,2,0,3))

		del data


		# Vertically integrating psi  
		data = Parallel(n_jobs=num_cores)(delayed(vpsi)(t) for t in xrange(time.shape[0]))

		psi_int = np.array(data)

		del data


		self.months = months
		self.latitude = latitude
		self.longitude = longitude
		self.depth = soest.LEV.copy()
		self.v = v
		self.psi = psi
		self.psi_int = psi_int

		del longitude, v, psi			


### Subroutines 
## Classic Dynamic Method from ARGO 
## dynamic height 
def dm(I):	
	H = soest.ADDEP[I,:,:,:].copy()
	LAT = soest.LAT16_100.copy()
	LON = soest.LONN69_25.copy()
	DEPTH = soest.LEV.copy()


	# Replacing the mask for NaN
	mask = H.mask
	H[mask] = np.nan


	LAT,DEPTH1,LON1 = np.meshgrid(LAT,DEPTH,LON)
	f0,DEPTH,LON = np.meshgrid(c.f,DEPTH,LON); f0 = f0[:,:,1:]

	del DEPTH1,LON1

	DX = np.diff(LON,axis=2)*np.cos(LAT[:,:,1:]*np.pi/180.)*(60.*1852.)  
	DH = np.diff(H,axis=2)


	# Geostrophic velocity computation
	# Obs: No reference is required, because its Absolute dynamic topography
	V = 10.*DH/(DX*f0)


	# Longitude axis
	LON = 0.5*(LON[0,0,1:]+LON[0,0,:-1])


	# Excluding the tropical areas
	ibad = np.abs(LAT[:,:,1:]) <= 5
	V[ibad] = np.nan


	return LON,V


## Horizontal Integration at each depth (Stream Function at different depths)
def cpsi(I):
	MONTH = soest.TIME.copy()
	DEPTH = soest.LEV.copy()
	LON = longitude.copy()
	LAT = soest.LAT16_100[I]	
	V = v[:,:,I,:].copy()


	# Closing the boundaries (continents + ocean at latitudes higher than Cape Town)
	ibad = np.isnan(V)
	V[ibad] = 0.

	ibad1 = LON >= 18
	V[:,:,ibad1] = 0.


	# Spatial resolution
	DX = np.cos(LAT*np.pi/180.)*(60.*1852.)


	# Horizontal integration	
	PSI = V.copy()*np.nan
	
	for j in xrange(1,PSI.shape[2],1):
		PSI[:,:,-j-1] = np.sum(V[:,:,-j:]*DX,axis=2)/1.e+6 


	PSI[ibad] = np.nan
	PSI[:,:,ibad1] = np.nan

	return PSI



## Vertical Integration of PSI
def vpsi(T):
	DEPTH = soest.LEV.copy()
	LON = longitude.copy()
	LAT = soest.LAT16_100.copy()

	LAT,DEPTH,LON = np.meshgrid(LAT,DEPTH,LON)

	PSI = psi[T,:,:,:].copy()
	PSI_INT = PSI.copy()*np.nan; PSI_INT[0,:,:] = 0.

	DZ = np.diff(DEPTH,axis=0)
	PSIA = 0.5*(PSI[1:,:,:]+PSI[:-1,:,:]) 

	for k in xrange(1,PSI.shape[0],1):
		PSI_INT[k,:,:] = np.sum(PSIA[:k,:,:]*DZ[:k,:,:],axis=0)

	return PSI_INT




### Main Program 
num_cores = 6

print "Loading data"
soest = models.load_netcdf('../data/soest_hw_1x1ADT_climatology.nc')

print "Computing constants to be used"
c = constants()

print "Computing relative geostrophic velocities, stream function and vertically integrated stream funcion"
geostrophic  = compute_geos()


# Replacing mask for Nan to save it
mask = soest.ADDEP.mask
soest.ADDEP[mask] = np.nan


print "Saving data"
scipy.io.savemat('../data/SOEST_SouthAtlantic.mat',
	{'ADDEP':np.array(soest.ADDEP),'latitude':soest.LAT16_100,'depth':soest.LEV,
	'longitude':soest.LONN69_25,'months':geostrophic.months})


scipy.io.savemat('../data/SOESTGeostrophic_SouthAtlantic.mat',
	{'psi':geostrophic.psi,'psi_int':geostrophic.psi_int,'v':geostrophic.v,
	'latitude':geostrophic.latitude,'depth':geostrophic.depth,'longitude':geostrophic.longitude,
	'months':geostrophic.months})
 

