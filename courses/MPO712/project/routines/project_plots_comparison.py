# -*- coding: utf-8 -*-
#
# project_plots_comparison.py
#
# purpose: Plotting and comparing the Sverdrup theory and geostrophic computations
# for the class project
#
# author: Tiago Bilo

# Science packages
import numpy as np 
import seawater as sw
import pylab as py
from scipy.io import loadmat
from netCDF4 import Dataset
from dateutil.parser import parse


# Plotting packages
import matplotlib.pyplot as plt
from matplotlib import ticker, dates
from mpl_toolkits.basemap import Basemap
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
class load_data(object):
	def __init__(self,filename):
		data = loadmat(filename)

		varnames = data.keys()

		for i in xrange(len(varnames)):
			if varnames[i] != '__header__' and varnames[i] != '__globals__' and varnames[i] != '__version__':
				exec("self."+varnames[i]+" = data['"+varnames[i]+"'].squeeze()")


### Subroutines 
## Study region map 
def plt_area(fign):

	fig = plt.figure(fign,facecolor='w')
	ax = plt.gca()

	# Projection + Map properties
	m = Basemap(projection='ortho',lat_0=-30,lon_0=-15,resolution='l',ax=ax)

	# Meridians and parallels
	m.drawmeridians(np.arange(0,360,30))
	m.drawparallels(np.arange(-90,90,30))

	# Topography
	m.etopo()

	# Study Area Polygon 
	rlon = [-70.,-60,-50,-40,-30,-20,-10,0,10,20,
	20,20.,20,20,20,20,20,20,20,
	10,0,-10,-20,-30,-40,-50,-60,-70,
	-70.,-70.,-70,-70,-70,-70,-70,-70,-70]

	rlat = [-5.,-5.,-5,-5,-5,-5,-5,-5,-5,-5,
	-10,-15,-20,-25,-30,-35,-40,-45,-50.,
	-50,-50,-50,-50,-50,-50,-50,-50.,-50,
	-45,-40,-35,-30,-25,-20,-15,-10,-5]


	rlon,rlat = m(rlon,rlat)
	ax.plot(rlon,rlat,color='k',linestyle='dashed',lw=5.)

	tkw = dict(fontsize=22,fontweight='demibold')

	# Box info
	tlon,tlat = m(25,-8)
	ax.text(tlon,tlat,ur'5$^{\circ}$S',**tkw)
	tlon,tlat = m(25,-53)
	ax.text(tlon,tlat,ur'50$^{\circ}$S',**tkw)

	tlon,tlat = m(-80,8)
	ax.text(tlon,tlat,ur'70$^{\circ}$W',rotation=25,**tkw)
	tlon,tlat = m(20,-2)
	ax.text(tlon,tlat,ur'20$^{\circ}$E',**tkw)


	# plt.draw()
	# plt.show(block=False)
	return fig


## SCOW - Sverdrup transport 
def map_labelaxis(N):
	"""	
	return a list of zeros and ones containing 
	where the axis labels of a map will be drawn 
	or not on a figure panel.

	n is the number of panels of each figure 

	"""


	if N == 12:
		LABELS = [[1,0,1,0],[0,0,1,0],[0,1,1,0],
		[1,0,0,0],[0,0,0,0],[0,1,0,0],
		[1,0,0,0],[0,0,0,0],[0,1,0,0],
		[1,0,0,1],[0,0,0,1],[0,1,0,1]]
	elif N == 16:
		LABELS = [[1,0,1,0],[0,0,1,0],[0,0,1,0],[0,1,1,0],
		[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,1,0,0],
		[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,1,0,0],
		[1,0,0,1],[0,0,0,1],[0,0,0,1],[0,1,0,1]]






	return LABELS



def plt_SvPsi(fign):

	fig, axes = plt.subplots(figsize=(10,10),nrows=4,ncols=3,facecolor='w')

	labels = map_labelaxis(12)
	lon,lat = np.meshgrid(scow_sv.longitude,scow_sv.latitude)

	k = 0

	for i in xrange(4):
		for j in xrange(3): 
			ax = axes[i,j]

			m = Basemap(projection='merc',
				resolution='l',
				llcrnrlat=llcrnrlat,
				llcrnrlon=llcrnrlon,
				urcrnrlat=urcrnrlat,
				urcrnrlon=urcrnrlon,
				ax=ax) 		

			m.drawcoastlines()
			m.fillcontinents(**kwc)			
			m.drawcountries()

			m.drawparallels(parallels,labels=labels[k],**kwl)
			m.drawmeridians(meridians,labels=labels[k],**kwl)


			if k == 0:
				kk = -1
				lon,lat = m(lon,lat)				
			else:
				kk = k-1



			cm = ax.contourf(lon,lat,scow_sv.data[1,kk,:,:],**kwpsi)
			ax.contour(lon,lat,scow_sv.data[1,kk,:,:],**kwpsic)

			# Info
			if k == 0:
				ax.set_ylabel('Summer',**kwseasons12)
			elif k == 3:
				ax.set_ylabel('Autumn',**kwseasons12)
			elif k == 6:
				ax.set_ylabel('Winter',**kwseasons12)
			elif k == 9:
				ax.set_ylabel('Spring',**kwseasons12)

			ax.yaxis.labelpad = 40

			lont,latt = m(-68,-15)			
			ax.text(lont,latt,soest_gs.months[kk][:3]+'.',**kwmonths12)


			# Colorbar
			if k == 10:
				#l,b,w,h
				cbaxes = fig.add_axes([0.15, 0.01, 0.75, 0.03])
				cb = fig.colorbar(cm,cax=cbaxes,extend='max',
					ticks=np.arange(0.,40.,5.),orientation='horizontal') 
				cb.ax.tick_params(labelsize=20)
				cb.set_label(ur'$\psi _{Sverdrup}$ (Sv)',fontsize=30,fontweight='bold')			


			k += 1

	fig.tight_layout()

	return fig


def plt_WGeosPsi(fign):

	fig, axes = plt.subplots(figsize=(10,10),nrows=4,ncols=3,facecolor='w')

	labels = map_labelaxis(12)
	lon,lat = np.meshgrid(scow_gs.longitude,scow_gs.latitude)

	k = 0

	for i in xrange(4):
		for j in xrange(3): 
			ax = axes[i,j]

			m = Basemap(projection='merc',
				resolution='l',
				llcrnrlat=llcrnrlat,
				llcrnrlon=llcrnrlon,
				urcrnrlat=urcrnrlat,
				urcrnrlon=urcrnrlon,
				ax=ax) 		

			m.drawcoastlines()
			m.fillcontinents(**kwc)			
			m.drawcountries()

			m.drawparallels(parallels,labels=labels[k],**kwl)
			m.drawmeridians(meridians,labels=labels[k],**kwl)


			if k == 0:
				kk = -1
				lon,lat = m(lon,lat)				
			else:
				kk = k-1



			cm = ax.contourf(lon,lat,scow_gs.data[1,kk,:,:],**kwpsi)
			ax.contour(lon,lat,scow_gs.data[1,kk,:,:],**kwpsic)

			# Info
			if k == 0:
				ax.set_ylabel('Summer',**kwseasons12)
			elif k == 3:
				ax.set_ylabel('Autumn',**kwseasons12)
			elif k == 6:
				ax.set_ylabel('Winter',**kwseasons12)
			elif k == 9:
				ax.set_ylabel('Spring',**kwseasons12)

			ax.yaxis.labelpad = 40

			lont,latt = m(-68,-15)			
			ax.text(lont,latt,soest_gs.months[kk][:3]+'.',**kwmonths12)


			# Colorbar
			if k == 10:
				#l,b,w,h
				cbaxes = fig.add_axes([0.15, 0.01, 0.75, 0.03])
				cb = fig.colorbar(cm,cax=cbaxes,extend='max',
					ticks=np.arange(0.,40.,5.),orientation='horizontal') 
				cb.ax.tick_params(labelsize=20)
				cb.set_label(ur'$\psi _{Sverdrup} - \psi _{Ekman}$ (Sv)'
					,fontsize=30,fontweight='bold')			


			k += 1

	fig.tight_layout()

	return fig


def plt_GeosPsi(fign):

	fig, axes = plt.subplots(figsize=(10,10),nrows=4,ncols=4,facecolor='w')

	labels = map_labelaxis(16)
	lon,lat = np.meshgrid(soest_gs.longitude,soest_gs.latitude)

	k = 0

	for i in xrange(4):
		for j in xrange(4): 
			ax = axes[i,j]

			m = Basemap(projection='merc',
				resolution='l',
				llcrnrlat=llcrnrlat,
				llcrnrlon=llcrnrlon,
				urcrnrlat=urcrnrlat,
				urcrnrlon=urcrnrlon,
				ax=ax) 		

			m.drawcoastlines()
			m.fillcontinents(**kwc)			
			m.drawcountries()

			m.drawparallels(parallels,labels=labels[k],**kwl)
			m.drawmeridians(meridians,labels=labels[k],**kwl)


			if k == 0:
				lon,lat = m(lon,lat)				


			# Selection of months and depths to be plotted
			if k < 4:
				M = 0
			elif k >= 4 and k < 8:
				M = 3
			elif k >= 8 and k < 12:
				M = 6
			elif k > 12:
				M = 9


			if k == 0 or k == 4 or k == 8 or k == 12:
				d = py.find(soest_gs.depth == 50)[0]
				kwpsi1 = dict(levels=np.arange(-0.02,0.1002,0.002)*100,extend='both',zorder=3)				
			elif k == 1 or k == 5 or k == 9 or k == 13:
				d = py.find(soest_gs.depth == 250)[0]
				kwpsi1 = dict(levels=np.arange(-0.01,0.051,0.001)*100,extend='both',zorder=3)									
			elif k == 2 or k == 6 or k == 10 or k == 14:
				d = py.find(soest_gs.depth == 800)[0]
				kwpsi1 = dict(levels=np.arange(-0.01,0.031,0.001)*100,extend='both',zorder=3)
			elif k == 3 or k == 7 or k == 11 or k == 15:
				d = py.find(soest_gs.depth == 2000)[0]						
				kwpsi1 = dict(levels=np.arange(-0.01,0.0155,0.0005)*100,extend='both',zorder=3)



			# Contour maps 
			cm = ax.contourf(lon,lat,soest_gs.psi[M,d,:,:]*100,**kwpsi1)
			ax.contour(lon,lat,soest_gs.psi[M,d,:,:],levels=[0.],
				colors='k',linewidths=1.5,zorder=4)

			# Info
			if k == 0:
				ax.set_ylabel('Summer (Jan.)',**kwseasons12)
			elif k == 4:
				ax.set_ylabel('Autumn (Apr.)',**kwseasons12)
			elif k == 8:
				ax.set_ylabel('Winter (Jul.)',**kwseasons12)
			elif k == 12:
				ax.set_ylabel('Spring (Oct.)',**kwseasons12)

			ax.yaxis.labelpad = 40

			lont,latt = m(-68,-15)			
			ax.text(lont,latt,"%i"%(soest_gs.depth[d]),**kwmonths12)


			# Colorbars
			if k >=12:
				if k < 14:
					cbn = 15-k+0.08
				else:
					cbn=15-k
				cbaxes = fig.add_axes([0.79-0.22*cbn, 0.04, 0.17, 0.01])				
				cb = fig.colorbar(cm,cax=cbaxes,extend='both',format="%1.1f"
					,orientation='horizontal') 
				cb.ax.tick_params(labelsize=8)

				# Setting number of ticks of the colorbar
				tick_locator = ticker.MaxNLocator(nbins=5)
				cb.locator = tick_locator
				cb.update_ticks()

				cb.set_label(ur'$\times$ 10$^{4}$ (m$^{2}$ s$^{-1}$)',
					fontsize=14,fontweight='bold')			


			k += 1

	fig.suptitle(ur"ARGO-derived $\psi _{geostrophic}$",fontsize=30,fontweight='demibold')
	fig.tight_layout()

	plt.draw()
	plt.show(block=False)
	return fig


def plt_psi_profile(fign):

	# Line properties
	colors = ['b','r','g','k','orange','m']
	linestyles = ['solid','dashed','dashdot',
	'dotted','dashed','solid']
	linewidths = [3,3.5,3.5,3.5,3.5,3.5,3]

	lat = np.arange(-15.625,-45.625,-5)


	fig = plt.figure(num=fign,figsize=(5,8),facecolor='w')
	ax = plt.gca()


	for i in xrange(lat.shape[0]):
		isoest = np.abs(soest_gs.latitude - lat[i]).argmin()

		psi = soestgs_annual_psi[:,isoest,:]
		j = py.find(~np.isnan(psi[0,:]))[0]
		psi = psi[:,j]

		ax.plot(psi*100,soest_gs.depth,color=colors[i],linestyle=linestyles[i],
			lw = linewidths[i],label=ur'%1.1f$^{\circ}$ S'%(np.abs(soest_gs.latitude[isoest])))

		ax.plot([0,0],[0,2000],'k',lw=5,zorder=1)

		del psi


	ax.set_ylim(2000,0)
	ax.set_xlim(-4,10)

	ax.set_yticks(range(0,2200,200))
	ax.set_xticks(np.arange(-4,12,2))
	ax.tick_params(labelsize=22)
	rstyle(ax)

	ax.set_ylabel('Depth (m)',fontsize=30,fontweight='demibold')
	ax.set_xlabel(ur'ARGO-derived $\bar{\psi}_{geostrophic} \times 10^{4}$ (m$^{2}$ s$^{-1}$)',fontsize=18,fontweight='demibold')

	ax.legend(numpoints=1,loc=0,prop={'size':18},frameon=False)

	# plt.draw()
	# plt.show(block=False)

	return fig 




def plt_IGeosPsi(fign):

	fig, axes = plt.subplots(figsize=(10,10),nrows=4,ncols=3,facecolor='w')

	labels = map_labelaxis(12)
	lon,lat = np.meshgrid(soest_gs.longitude,soest_gs.latitude)

	k = 0

	for i in xrange(4):
		for j in xrange(3): 
			ax = axes[i,j]

			m = Basemap(projection='merc',
				resolution='l',
				llcrnrlat=llcrnrlat,
				llcrnrlon=llcrnrlon,
				urcrnrlat=urcrnrlat,
				urcrnrlon=urcrnrlon,
				ax=ax) 		

			m.drawcoastlines()
			m.fillcontinents(**kwc)			
			m.drawcountries()

			m.drawparallels(parallels,labels=labels[k],**kwl)
			m.drawmeridians(meridians,labels=labels[k],**kwl)


			if k == 0:
				kk = -1
				lon,lat = m(lon,lat)				
			else:
				kk = k-1


			d = py.find(soest_gs.depth == 1200)[0]
			cm = ax.contourf(lon,lat,soest_gs.psi_int[kk,d,:,:],**kwpsi)
			ax.contour(lon,lat,soest_gs.psi_int[kk,d,:,:],**kwpsic)

			# Info
			if k == 0:
				ax.set_ylabel('Summer',**kwseasons12)
			elif k == 3:
				ax.set_ylabel('Autumn',**kwseasons12)
			elif k == 6:
				ax.set_ylabel('Winter',**kwseasons12)
			elif k == 9:
				ax.set_ylabel('Spring',**kwseasons12)

			ax.yaxis.labelpad = 40

			lont,latt = m(-68,-15)			
			ax.text(lont,latt,soest_gs.months[kk][:3]+'.',**kwmonths12)


			# Colorbar
			if k == 10:
				#l,b,w,h
				cbaxes = fig.add_axes([0.15, 0.01, 0.75, 0.03])
				cb = fig.colorbar(cm,cax=cbaxes,extend='max',
					ticks=np.arange(0.,40.,5.),orientation='horizontal') 
				cb.ax.tick_params(labelsize=20)
				cb.set_label(ur'ARGO-derived $\int_{0 m}^{1200 m} \psi_{geostrophic} dz$ (Sv)',fontsize=25,fontweight='bold')			


			k += 1

	fig.tight_layout()

	# plt.draw()
	# plt.show(block=False)
	return fig


def plt_WBCt(fign):

	lat = np.arange(-15.625,-45.625,-5)


	fig, axes = plt.subplots(figsize=(5,10),nrows=lat.shape[0],ncols=1,facecolor='w')

	time = []
	for t in xrange(soest.months.shape[0]):
		time.append(parse(soest.months[t]+'-01-01'))
	time = np.array(time)

	for i in xrange(lat.shape[0]):
		ax = axes[i]

		iscow = py.find(scow_sv.latitude == lat[i])[0]
		isoest = np.abs(soest_gs.latitude - lat[i]).argmin()

		argo = np.ones(time.shape[0])*np.nan
		wind_sv = argo.copy()
		wind_ek = argo.copy()
		wind_gs = argo.copy()		

		for t in xrange(time.shape[0]):
			d = py.find(soest_gs.depth == 1200)[0]
			gs = soest_gs.psi_int[t,d,isoest,:]
			argo[t] = gs[~np.isnan(gs)][0]

			gs = scow_sv.data[1,t,iscow,:]
			wind_sv[t] = gs[~np.isnan(gs)][0] 

			gs = scow_ek.data[2,t,iscow,:]
			wind_ek[t] = gs[~np.isnan(gs)][0]

			gs = scow_gs.data[1,t,iscow,:]
			wind_gs[t] = gs[~np.isnan(gs)][0]


		ax.plot(time,argo,'b',lw=2.5,label="ARGO")
		ax.plot(time,wind_sv,'k',lw=2.5,label="Sverdrup")
		ax.plot(time,wind_ek,'g--',lw=2,label="Ekman")			
		ax.plot(time,wind_gs,'r--',lw=2,label="Sverdrup - Ekman")

		ax.set_xlim(time[0], time[-1])
		ax.set_title(ur"%1.1f$^{\circ}$S"%(np.abs(soest_gs.latitude[isoest])),
			fontsize=22,fontweight='demibold')



		if i == lat.shape[0]-1:
			adjust_spines(ax, ['left', 'bottom'])

			months = dates.MonthLocator()
			dfmt = dates.DateFormatter('%b')

			ax.xaxis.set_major_locator(months)
			ax.xaxis.set_major_formatter(dfmt)

			ax.legend(numpoints=1,loc=3,prop={'size':13},frameon=False,mode="expand",
				ncol=2,bbox_to_anchor=(0., -1., 1., .102))	


		else:
			adjust_spines(ax, ['left'])

	fig.suptitle(ur"Transport at X$_{w}$ (Sv)"
		,fontsize=25,fontweight='demibold')

	# plt.draw()
	# plt.show(block=False)
	return fig


# Nice Spine axis adjustment
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])




def plt_psi_annual_mean(fign):

	fig, axes = plt.subplots(figsize=(10,5),nrows=1,ncols=2,facecolor='w')

	ax1 = axes[0]
	ax2 = axes[1]


	# Maps background
	m1 = Basemap(projection='merc',
		resolution='l',
		llcrnrlat=llcrnrlat,
		llcrnrlon=llcrnrlon,
		urcrnrlat=urcrnrlat,
		urcrnrlon=urcrnrlon,
		ax=ax1) 		

	m2 = Basemap(projection='merc',
		resolution='l',
		llcrnrlat=llcrnrlat,
		llcrnrlon=llcrnrlon,
		urcrnrlat=urcrnrlat,
		urcrnrlon=urcrnrlon,
		ax=ax2)

	m1.drawcoastlines()
	m1.fillcontinents(**kwc)			
	m1.drawcountries()

	m1.drawparallels(parallels,labels=[1,0,1,1],**kwl1)
	m1.drawmeridians(meridians,labels=[1,0,1,1],**kwl1)


	m2.drawcoastlines()
	m2.fillcontinents(**kwc)			
	m2.drawcountries()

	m2.drawparallels(parallels,labels=[0,1,1,1],**kwl1)
	m2.drawmeridians(meridians,labels=[0,1,1,1],**kwl1)


	# Contours
	lon1,lat1 = np.meshgrid(scow_gs.longitude,scow_gs.latitude)
	lon2,lat2 = np.meshgrid(soest_gs.longitude,soest_gs.latitude)

	lon1,lat1 = m1(lon1,lat1)
	cm1 = ax1.contourf(lon1,lat1,scowgs_annual_psi,**kwpsi)
	c1 = ax1.contour(lon1,lat1,scowgs_annual_psi,**kwpsic1)

	d = py.find(soest_gs.depth == 1200)[0]
	lon2,lat2 = m2(lon2,lat2)	
	cm2 = ax2.contourf(lon2,lat2,soestgs_annual_Ipsi[d,:,:],**kwpsi)
	c2 = ax2.contour(lon2,lat2,soestgs_annual_Ipsi[d,:,:],**kwpsic1)



	# Title of each panel
	ax1.set_title(ur'$\bar{\psi}_{Sverdrup} - \bar{\psi}_{Ekman}$',
		fontsize=25,fontweight='demibold',y=1.12)

	ax2.set_title(ur'$\int_{0 m}^{1200 m} \bar{\psi}_{geostrophic} dz$',
		fontsize=25,fontweight='demibold',y=1.12)	

	# Colorbar
	#l,b,w,h
	cbaxes = fig.add_axes([0.15, 0.16, 0.75, 0.05])
	cb = fig.colorbar(cm1,cax=cbaxes,extend='max',
		ticks=np.arange(0.,37.5,2.5),orientation='horizontal')
	cb.ax.tick_params(labelsize=18)
	cb.set_label(ur'$\bar{\psi}$ (Sv)',fontsize=32,fontweight='bold')	



	# Contour labels
	pos1 = [m1(5,-15),m1(-20,-13),m1(-7,-22),m1(-20,-25),m1(-40,-27)]
	ax1.clabel(c1,fmt='%i',colors='k',fontsize=18,
		fontweight='demibold',inline=1,manual=pos1)

	pos2 = [m2(5,-15),m2(-10,-20),m2(-20,-25)]
	ax2.clabel(c2,fmt='%i',colors='k',fontsize=18,
		fontweight='demibold',inline=1,manual=pos2)


	plt.draw()
	plt.show(block=False)
	return fig



def plt_annual_WBCt(fign):


	fig = plt.figure(num=fign,figsize=(8,5),facecolor='w')
	ax = plt.gca()


	argo = np.ones(soest_gs.latitude.shape[0])*np.nan
	wind = np.ones(scow_gs.latitude.shape[0])*np.nan


	d = py.find(soest_gs.depth == 1200)[0]
	for i in xrange(argo.shape[0]):
		gs = soestgs_annual_Ipsi[d,i,:]
		good  = py.find(~np.isnan(gs))

		if len(good) > 0:
			argo[i] = gs[good[0]]

		del gs, good


	for i in xrange(wind.shape[0]):
		gs = scowgs_annual_psi[i,:]
		good  = py.find(~np.isnan(gs))

		if len(good) > 0:		
			wind[i] = gs[good[0]]

		del gs, good

	ax.plot(soest_gs.latitude,argo,'b',lw=3.5,label='ARGO')
	ax.plot(scow_gs.latitude,wind,'r',lw=3.5,label='Sverdrup - Ekman')

	ax.legend(numpoints=1,loc=0,prop={'size':20},frameon=False)
	ax.plot([-60,0],[0,0],'k--',lw=3.5)

	ax.set_ylim(-40,40)
	ax.set_xlim(-5,-50)

	ax.set_yticks(range(-40,50,10))
	ax.set_xticks(range(-50,0,5))
	ax.tick_params(labelsize=18)
	rstyle(ax)

	labels = [ur'5$^{\circ}$S',ur'10$^{\circ}$S',ur'15$^{\circ}$S',ur'20$^{\circ}$S',
	ur'25$^{\circ}$S',ur'30$^{\circ}$S',ur'35$^{\circ}$S',ur'40$^{\circ}$S',ur'45$^{\circ}$S',
	ur'50$^{\circ}$S']

	labels.reverse()

	ax.set_xticklabels(labels)

	ax.set_ylabel('Transport at X$_{w}$ (Sv)',fontsize=30,fontweight='demibold')
	ax.set_xlabel(ur'Latitude',fontsize=30,fontweight='demibold')

	plt.draw()
	plt.show(block=False)

	return fig







### Main program 
print "Loading Scow-derived Products"
scow = load_data('../data/Scow_SouthAtlantic.mat')
scow_sv = load_data('../data/ScowSverdrup_SouthAtlantic.mat')
scow_ek = load_data('../data/ScowEkman_SouthAtlantic.mat')
scow_gs = load_data('../data/ScowGeostrophic_SouthAtlantic.mat')

print "Loading the ADDEP and ADDEP-derived Products"
soest = load_data('../data/SOEST_SouthAtlantic.mat')
soest_gs = load_data('../data/SOESTGeostrophic_SouthAtlantic.mat')


# Annual means
scowsv_annual_psi = np.nanmean(scow_sv.data[1,:,:,:],axis=0)
scowgs_annual_psi = np.nanmean(scow_gs.data[1,:,:,:],axis=0)
soestgs_annual_psi = np.nanmean(soest_gs.psi,axis=0)
soestgs_annual_Ipsi = np.nanmean(soest_gs.psi_int,axis=0)


print "Plotting Introduction figures"

## Defining Maps limits, masks and etc
llcrnrlat = -50. 
urcrnrlat = -5.
llcrnrlon = -70.
urcrnrlon = 20.

# Continents
kwc = dict(color='wheat',lake_color='w')

# Parallels and meridians
parallels = range(-10,-60,-10)
meridians = range(-60,40,20)

kwl = dict(linewidth=0.1,fontsize=14,
	fontweight='demibold')	

kwl1 = dict(linewidth=0.1,fontsize=20,
	fontweight='demibold')


# Psi
kwpsi=dict(levels=np.arange(0.,35.25,0.25),extend='max',zorder=3)
kwpsic=dict(levels=[0.,2.,5.,10.,20.,30.,40.],
				colors='k',linewidths=1.5,zorder=4)

kwpsic1=dict(levels=[0.,2.,5.,10.,20.,30.,40.],
				colors='k',linewidths=3.5,zorder=4)

# Axis labels and info
kwseasons12 = dict(fontsize=18,fontweight='demibold')
kwmonths12 = dict(color='k',fontweight='demibold',fontsize=16,zorder=8)




# Plots 
# fig1 = plt_area(1)
fig2 = plt_SvPsi(2)
fig3 = plt_WGeosPsi(3)
fig4 = plt_GeosPsi(4)
fig5 = plt_psi_profile(5)
fig6 = plt_IGeosPsi(6)
fig7 = plt_WBCt(7)
fig8 = plt_psi_annual_mean(8)
fig9 = plt_annual_WBCt(9)

print "Saving Figures"
# fig1.savefig('../figures/area.png',format='png',bbox_inches='tight',
# 	transparent=True,pad_inches=0.5)

fig2.savefig('../figures/sverdrup_psi.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig3.savefig('../figures/wind_geostrophic_psi.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig4.savefig('../figures/argo_geostrophic_psi.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig5.savefig('../figures/argo_geostrophic_psi_prof.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig6.savefig('../figures/argo_zintegrated_psi.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig7.savefig('../figures/wbl_transp.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig8.savefig('../figures/psi_annual_mean.png',format='png',bbox_inches='tight',pad_inches=0.5)

fig9.savefig('../figures/wbl_annual_transp.png',format='png',bbox_inches='tight',pad_inches=0.5)

