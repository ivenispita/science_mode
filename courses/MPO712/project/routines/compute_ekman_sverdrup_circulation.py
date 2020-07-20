# -*- coding: utf-8 -*-
#
# compute_ekman_sverdrup_circulation.py
#
# purpose:  Compute the ekman, sverdrup and vertically integrated geostrophic transport and Stream Function
#			from the SCOW climatology in the South Atlantic
#
# author:   Tiago Bilo

# Scientific packages

import numpy as np
import pylab as py
import seawater as sw
from netCDF4 import Dataset
import scipy.io


# Paralelization of the code 
from joblib import Parallel, delayed  
import multiprocessing



import warnings
warnings.filterwarnings('ignore')


### Classes 
class load_scow(object):
	def __init__(self,path):
		file1 = path+'scow_taux_climatology.nc'
		file2 = path+'scow_tauy_climatology.nc'
		file3 = path+'scow_wind_stress_curl_climatology.nc'

		variables = ['taux','tauy','stress_curl']

		data1 = Dataset(file1,'r')
		data2 = Dataset(file2,'r')
		data3 = Dataset(file3,'r')

		# Coordinates systems - South Atlantic Coordinates 
		latitude = data1.variables.pop('latitude')[:].squeeze()
		longitude = data1.variables.pop('longitude')[:].squeeze()

		jgood = (longitude>=(-69.5+360.)) | (longitude<=24.5)
		jgood = py.find(jgood)

		shift = len(longitude[longitude>=(-69.5+360.)])
		jgood = np.roll(jgood,shift)

		igood = (latitude>=-74.5) & (latitude<=9.5)
		igood = py.find(igood)


		latitude = latitude[igood]
		longitude = longitude[jgood]; longitude[longitude>180] = longitude[longitude>180]-360.


		igood,jgood = np.meshgrid(igood,jgood)


		# Organizing the data
		months = data1.variables.keys()

		data = np.ones((len(variables),len(months),latitude.shape[0],longitude.shape[0]))*np.nan

		for i in xrange(len(months)):
			taux = data1[months[i]][:][igood,jgood].transpose()
			tauy = data2[months[i]][:][igood,jgood].transpose()
			curl = data3[months[i]][:][igood,jgood].transpose()

			# Replacing the missing filling number -9999.0 by nan
			ibad = taux == -9999.0
			taux[ibad] = np.nan
			ibad = tauy == -9999.0			
			tauy[ibad] = np.nan
			ibad = curl == -9999.0			
			curl[ibad] = np.nan

			# Fixing units (everything should be in the SI)
			curl = curl*1.e-7

			month = np.array([taux,tauy,curl])
			data[:,i,:,:] = month

		self.data = data 
		self.months = months
		self.variables = variables
		self.latitude = latitude
		self.longitude = longitude

class constants(object):
	def __init__(self):
		# Beta = meridional rate of change of the Coriolis parameters (m^-1 s^-1)
		r = 6.371e+6 														# Earth's radius (m)
		omega = 7.292115e-5 												# Earth's rotation rate (s^-1)
		beta = 2.*omega*np.cos(np.pi*scow.latitude/180.)/r 

		# Coriolis parameter (s^-1)
		f = sw.f(scow.latitude)

		# Rho = typical density of the sea water above the pycnocline (kg m^-3)
		# Since we are considering a incompressible ocean (i. e., Boussinesq approximation) it will be 
		# constant
		rho = 1025.0


		self.beta = beta 
		self.rho  = rho 
		self.f = f

# Sverdrup Transport + Stream function (Geostrophic + Ekman)
class compute_sv(object):
	def __init__(self):
		global transport 

		variables = ['meridional_transport','psi']
		num_cores = 6

		data = np.ones((len(variables),len(scow.months),scow.latitude.shape[0],scow.longitude.shape[0]))*np.nan

		beta = c.beta.repeat(scow.longitude.shape[0]).reshape((c.beta.shape[0],scow.longitude.shape[0]))

		for i in xrange(scow.data.shape[1]):
			transport = scow.data[2,i,:,:]/beta

			psi = Parallel(n_jobs=num_cores)(delayed(integration)(lat) for lat in scow.latitude)
			psi = np.array(psi)

			D = np.array([transport.copy()/(c.rho*1.e+6),psi.copy()/(c.rho*1.e+6)])
			data[:,i,:,:] = D			

		del transport 

		# Here I could derive psi in y and get the zonal sverdrup transport (I think I won't need it)


		# "Isolating" the subtropical gyre 
		ibad = (np.abs(scow.latitude) <= 5) | (np.abs(scow.latitude) >= 50)
		data[:,:,ibad,:] = np.nan

		self.latitude = scow.latitude
		self.longitude = scow.longitude
		self.variables = variables
		self.data = data


class compute_ekman(object):
	def __init__(self):
		global transport 

		variables = ['zonal_transport','meridional_transport','psi']
		num_cores = 6

		data = np.ones((len(variables),len(scow.months),scow.latitude.shape[0],scow.longitude.shape[0]))*np.nan

		f = c.f.repeat(scow.longitude.shape[0]).reshape((c.f.shape[0],scow.longitude.shape[0]))

		for i in xrange(scow.data.shape[1]):
			zonal_transport = scow.data[1,i,:,:]/(f*c.rho)
			meridional_transport = -scow.data[0,i,:,:]/(f*c.rho)

			transport = meridional_transport.copy()

			psi = Parallel(n_jobs=num_cores)(delayed(integration)(lat) for lat in scow.latitude)
			psi = np.array(psi)

			D = np.array([zonal_transport.copy()/1.e+6,meridional_transport.copy()/1.e+6,psi.copy()/1.e+6])
			data[:,i,:,:] = D			

		del transport

		# "Isolating" the subtropical gyre 
		ibad = (np.abs(scow.latitude) <= 5) | (np.abs(scow.latitude) >= 50)
		data[:,:,ibad,:] = np.nan

		self.latitude = scow.latitude
		self.longitude = scow.longitude
		self.variables = variables
		self.data = data


class compute_geos(object):
	def __init__(self):

		data = sverdrup.data.copy()-ekman.data[1:,:,:,:].copy()

		self.data = data
		self.latitude = scow.latitude
		self.longitude = scow.longitude
		self.variables = ['meridional_transport','psi']


### Subroutines 
def integration(LAT):
	i = np.abs(scow.latitude-LAT).argmin() 
	Tr = transport[i,:].copy()

	ibad = np.isnan(Tr)
	Tr[ibad] = 0.

	# Closing the eastern boundary at latitudes higher than Cape Town's
	ibad1 = scow.longitude >= 18
	Tr[ibad1] = 0. 


	# Spatial Resolution 
	DX = sw.dist([scow.latitude[i],scow.latitude[i]],
		[scow.longitude[0],scow.longitude[1]])[0][0]*1000.


	# Horizontal integration	
	PSI = np.zeros(Tr.shape[0])

	for j in xrange(1,Tr.shape[0],1):
		PSI[-j-1] = np.sum(Tr[-j:]*DX) 

	PSI[ibad] = np.nan
	PSI[ibad1] = np.nan
	return PSI



### Main Program
print "Loading data"
scow = load_scow('../data/')

print "Computing constants to be used"
c = constants()

print "Computing Sverdrup Transport + Stream function"
sverdrup = compute_sv()

print "Computing Ekman transport + Stream function"
ekman = compute_ekman()

print "Geostrophic transport + Stream function"
geostrophic  = compute_geos()


print "Saving data"
scipy.io.savemat('../data/Scow_SouthAtlantic.mat',
	{'data':scow.data,'latitude':scow.latitude,'longitude':scow.longitude,'variables':scow.variables,'months':scow.months})

scipy.io.savemat('../data/ScowSverdrup_SouthAtlantic.mat',
	{'data':sverdrup.data,'latitude':sverdrup.latitude,'longitude':sverdrup.longitude,'variables':sverdrup.variables,'months':scow.months})

scipy.io.savemat('../data/ScowEkman_SouthAtlantic.mat',
	{'data':ekman.data,'latitude':ekman.latitude,'longitude':ekman.longitude,'variables':ekman.variables,'months':scow.months})

scipy.io.savemat('../data/ScowGeostrophic_SouthAtlantic.mat',
	{'data':geostrophic.data,'latitude':geostrophic.latitude,'longitude':geostrophic.longitude,'variables':geostrophic.variables,'months':scow.months})
