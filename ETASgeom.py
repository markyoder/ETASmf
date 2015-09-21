#
# a study of aftershock geometry, namely large-aftershock distributions.
#
#

import BASScast as bcp
import ANSStools as atp


import numpy
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from geographiclib.geodesic import Geodesic as ggp
#
import datetime as dtm
import matplotlib.dates as mpd
import pytz
#
import os
import math
import glob
#
import multiprocessing as mpp
#
lon2km = 111.1
deg2rad = 2.*math.pi/360.

# tohoku:
# a=egp.etas_auto(to_dt=None, lat_center=38.32, lon_center=142.37, Lr_map_factor=2., mc=5.0, catlen=5.*365., dt_0=None)
#tohoku_params = {}

# gorkha:



def etas_auto(lon_center=None, lat_center=None, d_lat_0=.25, d_lon_0=.5, dt_0=10, Lr_map_factor=5.0, mc=2.5, mc_0=None, dm_cat=2.0, gridsize=.1, to_dt=None, fnameroot='etas_auto', catlen=5.0*365.0, d_lambda=1.76, doplot=False):
	'''
	# , Lr_as_factor=.5
	# a version of this exists also in ETASscripts.py, but for now let's just develop separately.
	#	
	# a starter script to auto-select some parameters for an ETAS run. in practice, give this script a center location (probably a mainshock epicenter).
	# the script will find the largest earthquake in that region and scale up an ETAS parameter set accordingly.
		# d_lat_0, d_lon_0, dt_0 are the starting catalog parameters (largest earthquake in d_lat_0 x d_lon_0 x dt_0 cube).
	'''
	#
	lw=2.5
	# catch some default value exceptions:
	if dt_0==None: dt_0=10
	#
	_colors =  mpl.rcParams['axes.color_cycle']
	Lr_as_factor=.5
	#
	#if to_dt == None: to_dt = dtm.datetime.now(pytz.timezone('UTC'))
	to_dt = (to_dt or dtm.datetime.now(pytz.timezone('UTC')))
	mc_0  = (mc_0 or mc)
	#
	if lon_center==None and lat_center==None:
		# let's look for any large earthquake in the world. assume for this, mc
		mc_0=6.0
		lat_center = 0.
		lon_center = 0.
		d_lat_0 = 88.
		d_lon_0 = -180.
	#
	# get a preliminary catalog:
	cat_0 = atp.catfromANSS(lon=[lon_center-d_lon_0, lon_center+d_lon_0], lat=[lat_center - d_lat_0, lat_center+d_lat_0], minMag=mc_0, dates0=[to_dt-dtm.timedelta(days=dt_0), to_dt], fout=None, rec_array=True)
	# diagnostic: output catalog length
	print "catalog length: %d" % len(cat_0)
	#
	# if there are no events, we probably are looking for a long range ETAS, so let's look for ALL events in the full catlen interval:
	if len(cat_0)==1:
		print "empty catalog. search over the whole catalog length..."
		# note: when we return as a recarray, an empty array has length 1 (pull up an empty array and figure out the details some time).
		cat_0 = atp.catfromANSS(lon=[lon_center-d_lon_0, lon_center+d_lon_0], lat=[lat_center - d_lat_0, lat_center+d_lat_0], minMag=mc_0, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)
	#
	#biggest_earthquake = filter(lambda x: x['mag']==max(cat_0['mag']), cat_0)[0]
	print "fetch preliminary catalog; find *mainshock*"
	mainshock = {cat_0.dtype.names[j]:x for j,x in enumerate(filter(lambda x: x['mag']==max(cat_0['mag']), cat_0)[0])}
	#print "biggest event(s): ", mainshock
	#
	# now, get new map domain based on rupture length, etc.
	#L_r = .5*10.0**(.5*mainshock['mag'] - 1.76)
	L_r = Lr(mainshock['mag'], fact=Lr_as_factor, d_lambda=d_lambda)
	d_lat = Lr_map_factor*L_r/lon2km
	d_lon = Lr_map_factor*L_r/(lon2km*math.cos(deg2rad*mainshock['lat']))
	lats = [mainshock['lat']-d_lat, mainshock['lat']+d_lat]
	lons = [mainshock['lon']-d_lon, mainshock['lon']+d_lon]
	print "lats, lons: ", lats, lons, L_r
	print "mainshock found: ", mainshock, ".\nNow, find aftershocks..."
	#
	#working_cat = atp.catfromANSS(lon=[mainshock['lon']-d_lon, mainshock['lon']+d_lon], lat=[mainshock['lat']-d_lat, mainshock['lat']+d_lat], minMag=mc, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)
	aftershock_cat = atp.catfromANSS(lon=lons, lat=lats, minMag=mainshock['mag']-dm_cat, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)
	#
	#cat_map = get_catalog_map(lats=lats, lons=lons, eq_cat = aftershock_cat)
	cat_map = get_catalog_map(lats=lats, lons=lons)
	#
	# now, let's draw circles around (some of) these events. this will be a bit sloppy, since we need to convert coordinates from dist -> lat/lon -> map.
	#
	for k,ev in enumerate(aftershock_cat):
		print "event: %d: " % k, ev
		for j,lr_fact in enumerate([.5, 1.0]):
			lr = Lr(ev['mag'], fact=lr_fact, d_lambda=d_lambda)
			#as_circle = circle_geo(lon0=ev['lon'], lat0=ev['lat'], R=lr*1000., d_theta=None, N=250, units_theta='deg', R_units='m')
			#Xa, Ya = cat_map(*zip(*as_circle))
			#
			if j==0:
				line_style = '--'
			else:
				line_style='-'
			#
			Xa, Ya = cat_map(*zip(*circle_geo(lon0=ev['lon'], lat0=ev['lat'], R=lr*1000., d_theta=None, N=250, units_theta='deg', R_units='m')))
			plt.plot(Xa, Ya, '%s' % line_style, zorder=7, color=_colors[k%len(_colors)], lw=lw)
		x,y = cat_map(ev['lon'], ev['lat'])
		ev_dtm = ev['event_date'].tolist()
		lbl_str = 'm=%.2f, %d-%d-%d %d:%d:%d' % (ev['mag'], ev_dtm.year, ev_dtm.month, ev_dtm.day, ev_dtm.hour,ev_dtm.minute,ev_dtm.second)
		#
		# plot pre-mainshock events as squares, post-mainshock events as circles:
		if ev['event_date']<mainshock['event_date']: marker_str = 's'
		if ev['event_date']==mainshock['event_date']: marker_str = '*'
		if ev['event_date']>mainshock['event_date']: marker_str = 'o'
		#
		plt.plot([x], [y], marker_str, ms=15.*ev['mag']/8.0, color=_colors[k%len(_colors)], label=lbl_str)
	plt.legend(loc=0, numpoints=1)
	#
	# and for now, let's plot time-links (lines betwen sequential events):
	for k, ev in enumerate(aftershock_cat[1:]):
		X,Y = cat_map([aftershock_cat[k]['lon'], ev['lon']], [aftershock_cat[k]['lat'], ev['lat']])
		plt.plot(X,Y, '.-', ms=8, color=_colors[k%len(_colors)])
	#
	f=plt.figure(1)
	plt.clf()
	ax3d = f.add_axes([.1, .1, .8, .8], projection='3d')
	ax3d.plot(aftershock_cat['lon'], aftershock_cat['lat'], aftershock_cat['event_date_float'], 'o-')
	#z = [mpd.date2num(dt.tolist()) for dt in aftershock_cat['event_date']]
	#ax3d.plot(aftershock_cat['lat'][1:], aftershock_cat['lon'][1:], [x for x in (aftershock_cat['event_date_float'][1:]-mainshock['event_date_float'])], 'o-')
	#ax3d.plot(aftershock_cat['lat'][1:], aftershock_cat['lon'][1:], [math.log10(x-z[0]) for x in z[1:]], 'o-')
	ax3d.set_ylabel('latitude')
	ax3d.set_xlabel('longitude')
	ax3d.set_zlabel('time')
	#print "times: ", [x-z[0] for x in z[1:]]
	#
	f=plt.figure(2)
	plt.clf()
	ax3d2 = f.add_axes([.1, .1, .8, .8], projection='3d')
	ax3d2.plot(aftershock_cat['lon'][1:], aftershock_cat['lat'][1:], numpy.log10(aftershock_cat['event_date_float'][1:]-aftershock_cat['event_date_float'][0]), 'o-')
	
	ax3d2.set_ylabel('latitude')
	ax3d2.set_xlabel('longitude')
	ax3d2.set_zlabel('$log(\\Delta t)$')
		
	print "log(dt): ", numpy.log10(aftershock_cat['event_date_float'][1:]-aftershock_cat['event_date_float'][0])
	
	f=plt.figure(3)
	plt.clf()
	ax3d3 = f.add_axes([.1, .1, .8, .8], projection='3d')
	ax3d3.plot(aftershock_cat['lon'][0:], aftershock_cat['lat'][0:], aftershock_cat['depth'][0:], 'o-')
	ax3d3.set_ylabel('latitude')
	ax3d3.set_xlabel('longitude')
	ax3d3.set_zlabel('depth $z$')
	#
	#
	#print "mainshock: ", mainshock, mainshock['event_date']
	cat = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[mpd.num2date(mainshock['event_date_float']-1.0), to_dt], fout=None, rec_array=True)
	f=plt.figure(4)
	plt.clf()
	ax3d = f.add_axes([.1, .1, .8, .8], projection='3d')
	ax3d.plot(cat['lon'], cat['lat'], cat['event_date_float'], 'o')
	#
	#
	ax3d.set_title('All quakes, lat, lon, time')
	ax3d.set_ylabel('latitude')
	ax3d.set_xlabel('longitude')
	ax3d.set_zlabel('time')
	
	return aftershock_cat
	#return atp.catfromANSS(lon=lons, lat=lats, minMag=mainshock['mag']-dm_cat, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=False)
	
	#return biggest_earthquake
#
def Lr(m, fact=.5, d_lambda=1.76):
	return fact*10.**(.5*m - d_lambda)
#
def get_catalog_map(lats=None, lons=None, eq_cat=[], map_res='i', map_projection='cyl', fignum=0, ax=None, do_clf=True):
	#
	if lats==None: lats = [31., 42.]
	if lons==None: lons = [-125., -114.]
	#
	if fignum!=None: plt.figure(fignum)
	if do_clf: plt.clf()
	if ax==None: ax=plt.gca()
	#
	cm = Basemap(llcrnrlon=lons[0], llcrnrlat=lats[0], urcrnrlon=lons[1], urcrnrlat=lats[1], resolution=map_res, projection=map_projection, lon_0=numpy.mean(lons), lat_0=numpy.mean(lats), ax=ax)
	cm.drawcoastlines(color='gray', zorder=1)
	cm.drawcountries(color='gray', zorder=1)
	cm.drawstates(color='gray', zorder=1)
	cm.drawrivers(color='gray', zorder=1)
	cm.fillcontinents(color='beige', zorder=0)
	# drawlsmask(land_color='0.8', ocean_color='w', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution='l', grid=5, **kwargs)
	#cm.drawlsmask(land_color='0.8', ocean_color='c', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution=mapres, grid=5)

	print "lat, lon ranges: ", lats, lons
	cm.drawmeridians(range(int(lons[0]), int(lons[1])), color='k', labels=[0,0,1,1])
	cm.drawparallels(range(int(lats[0]), int(lats[1])), color='k', labels=[1, 1, 0, 0])
	#
	if eq_cat!=None:
		# can we also assign sizes dynamically, like colors?
		if hasattr(eq_cat, 'dtype'):
			# it's a recrray:
			X,Y = cm(eq_cat['lon'], eq_cat['lat'])
		else:
			# it's probably a list... though we should check for dicts, etc.
			X,Y = cm([rw[2] for rw in eq_cat], [rw[1] for rw in eq_cat])
		#
		cm.plot(X,Y, '.')
				
	#
	return cm

def circle_xy(x0=0., y0=0., R=1.0, d_theta=1., N=None, units_theta='deg'):
	if d_theta==None or N!=None:
		d_theta = 2.0*math.pi/float(max(1,N-1))
		units_theta = 'rad'
	if units_theta in ('deg', 'degrees', 'd'):
		d_theta*=deg2rad
	#
	return [[x0+math.cos(theta), y0+math.sin(theta)] for theta in numpy.arange(0., 2.0*math.pi+d_theta, d_theta)]
#
def circle_geo(lon0=0., lat0=0., R=1.0, d_theta=1., N=None, units_theta='deg', R_units='km'):
	# geodesic uses degrees...
	#
	if d_theta==None or N!=None:
		d_theta = 360./float(max(1,N-1))
	#
	if units_theta.lower() in ('rad', 'radians'):
		d_theta/=deg2rad
	#
	#print "geo d_theta: ", d_theta
	R_factor = 1.0
	if R_units.lower()=='km': R_factor = 1000.
	#return [[x0+math.cos(theta), y0+math.sin(theta)] for theta in numpy.arange(0., 2.0*math.pi+d_theta, d_theta)]
	#
	XY = []
	for theta in numpy.arange(0., 360.+d_theta, d_theta):
		z = ggp.WGS84.Direct(lat0, lon0, theta%360., R*R_factor)
		XY += [[z['lon2'], z['lat2']]]
	return XY
		
	#return [[z['lon2'], z['lat2']] for z in [ggp.WGS84.Direct(lat0, lon0, theta, R*R_factor) for theta in numpy.arange(0., 2.0*math.pi+d_theta, d_theta)]]
	

if __name__=='__main__':
	print "do __main__ stuff..."
	vc_parser.mpl.use('Agg')
else:
	plt.ion()
	
	
	
	
	
	
