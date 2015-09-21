#
# a study of aftershock geometry, namely large-aftershock distributions.
#
#

import BASScast as bcp
import ANSStools as atp
import ETASscripts as esp

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

# mainshock data+1
sischuan_prams = {'to_dt':dtm.datetime(2008,6,12, tzinfo=pytz.timezone('UTC')), 'mainshock_dt':dtm.datetime(2008,5,13, tzinfo=pytz.timezone('UTC')), 'lat_center':31.021, 'lon_center':103.367, 'Lr_map_factor':4.0, 'mc':4.0, 'mc_0':None, 'dm_cat':2.0, 'gridsize':.1, 'fnameroot':'etas_auto_sichuan', 'catlen':10.0*365., 'd_lambda':1.76, 'doplot':True}
#sischuan_prams['to_dt'] = dtm.datetime(2014,4,20, tzinfo=pytz.timezone('UTC'))
sischuan_prams['todt'] = dtm.datetime.now(pytz.timezone('UTC'))

def mod_kwargs(prams_dict, **kwargs):
	r_dict = prams_dict.copy()
	r_dict.update(kwargs)
	return r_dict

def chengdus():
	A=chengdu_etas(prams=sischuan_prams, catlen=5.*365.)
	plt.figure(1)
	plt.title('Chengdu ETAS: %s\n' % str(sischuan_prams['todt']))
	plt.savefig(os.path.join('/home/myoder/Dropbox/Research/ACES/China2015/talks/nepal/images', 'chengdu_etas_5yr_b.png'))
	f3d = contour_3d_etas(A)
	
	A=chengdu_etas(prams=sischuan_prams, catlen=10.*365.)
	plt.figure(1)
	plt.title('Chengdu ETAS: %s\n' % str(sischuan_prams['todt']))
	plt.savefig(os.path.join('/home/myoder/Dropbox/Research/ACES/China2015/talks/nepal/images', 'chengdu_etas_10yr_b.png'))
	
	return A


def chengdu_etas(prams=sischuan_prams, dlat=3., dlon=3., **kwargs):
	# ,to_dt=dtm.datetime.now(pytz.timezone('UTC'))
	#A=esp.makeETASFCfiles(**mod_kwargs(sichuan_prams, catlen=3650.))
	#
	prams = mod_kwargs(prams,**kwargs)
	#prams.update({'lats':[31.021-dlat, 31.021+dlat], 'lons':[103.367-dlon, 103.367+dlon]})
	print "prams: ", prams['todt']
	#
	A=esp.make_etas_fcfiles(prams, lats=[31.021-dlat, 31.021+dlat], lons=[103.367-dlon, 103.367+dlon])
	#A=esp.makeETASFCfiles(todt=None, gridsize=.1, contres=3, mc=4.0, kmldir='kml', catdir='kml', fnameroot='chengdu_etas', catlen=prams['catlen'], doplot=False, lons=prams['lons'], lats=prams['lats'], bigquakes=[], bigmag=6.5, addquakes=[], eqeps=None, eqtheta=None, fitfactor=5.0, cmfnum=0, fignum=1, colorbar_fontcolor='k', contour_intervals=None, rtype='ssim', contour_top=1.0, contour_bottom=0.0, p_quakes=None, p_map=None, fnameroot_suffix='', maxNquakes=None)
	#
	chengdu_lon = 104.+4./60.
	chengdu_lat = 30. + 2./3.
	#
	A.BASScastContourMap(maxNquakes=5)	# could also use lats=[], lons=[]
	x,y = A.cm(chengdu_lon, chengdu_lat)
	A.cm.plot([x], [y], 'r*', ms=15, zorder=8, alpha=.7)
	A.cm.plot([x], [y], 'r*', ms=18, zorder=7, alpha=.7)
	#
	for rw in A.catalog:
		if rw[3]<6.5: continue
		if rw[0]<mpd.date2num(dtm.datetime(2007,1,1, tzinfo=pytz.timezone('UTC'))): continue
		#
		x,y = A.cm(rw[2], rw[1])
		dt = mpd.num2date(rw[0])
		plt.plot([x], [y], 'o', ms=9.*rw[3]/8., zorder=5, label='m=%.2f, %d/%d/%d' % (rw[3], dt.year, dt.month, dt.day))
	#
	plt.legend(loc=0, numpoints=1)
	#
	return A
>>>>>>> ab2fc4c2b1af22439cf699efdca7a18cdc3b8980

def etas_auto(lon_center=None, lat_center=None, d_lat_0=.25, d_lon_0=.5, dt_0=10, Lr_map_factor=5.0, mc=2.5, mc_0=None, dm_cat=2.0, gridsize=.1, to_dt=None, mainshock_dt=None, fnameroot='etas_auto', catlen=5.0*365.0, d_lambda=1.76, doplot=False, show_legend=True):
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
	mainshock_dt = (mainshock_dt or to_dt)
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
	#cat_0 = atp.catfromANSS(lon=[lon_center-d_lon_0, lon_center+d_lon_0], lat=[lat_center - d_lat_0, lat_center+d_lat_0], minMag=mc_0, dates0=[to_dt-dtm.timedelta(days=dt_0), to_dt], fout=None, rec_array=True)
	cat_0 = atp.catfromANSS(lon=[lon_center-d_lon_0, lon_center+d_lon_0], lat=[lat_center - d_lat_0, lat_center+d_lat_0], minMag=mc_0, dates0=[to_dt-dtm.timedelta(days=dt_0), mainshock_dt], fout=None, rec_array=True)
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
	primary_cat = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)
	#
	if not doplot: return primary_cat
	aftershock_cat = numpy.core.records.fromarrays(zip(*filter(lambda rw: rw['mag']>=(mainshock['mag']-dm_cat), primary_cat)),dtype=primary_cat.dtype)
	
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
		plt.plot([x], [y], marker_str, ms=15.*ev['mag']/8.0, color=_colors[k%len(_colors)], label=lbl_str, zorder=7)
	if show_legend: plt.legend(loc=0, numpoints=1)
	plt.plot(*cat_map(primary_cat['lon'], primary_cat['lat']), marker='.', ls='', ms=3., zorder=5, alpha=.6)
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
	
	#return aftershock_cat
	return primary_cat
	#return atp.catfromANSS(lon=lons, lat=lats, minMag=mainshock['mag']-dm_cat, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=False)
	
	#return biggest_earthquake
#
def moment_figure():
	c1 = etas_auto(**mod_kwargs(sichuan_prams, doplot=False))
	M = lambda x: 10.**(1.0*x)
	#
	M_cum = [M(c1['mag'][0])]
	for m in c1['mag'][1:]:
		M_cum+=[M_cum[-1]+M(m)]
	#
	#f_cat = [[mpd.date2num(rw[0].tolist()), mpd.date2num(rw['event_date'].tolist())-mpd.date2num(c1[j]['event_date'].tolist()), rw['mag']] for j,rw in c1]
	#
	plt.figure(17)
	plt.clf()
	ax=plt.gca()
	mcy = min(c1['mag'])
	print "mcy: ", mcy
	#return c1['event_date'], [mcy for x in c1['mag']]
	ax.vlines([x.tolist() for x in c1['event_date']], [mcy for x in c1['mag']], c1['mag'], lw=2., color='r')
	ax2=ax.twinx()
	ax2.set_yscale('linear')
	ax2.plot([x.tolist() for x in c1['event_date']], M_cum, '.-', lw=2.5)
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
	
	
	
	
	
	
