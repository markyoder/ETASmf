# china_2015_scripts.py
#
# Primary Author: Mark R. Yoder, Ph.D.
#                 mryoder@ucdavis.edu
#                 mark.yoder@gmail.com
#
# A bunch of scripts to develop figures, data, and materials for talks for talks given at the 2015 ACES meeting in Chengdu China.
# These scripts are pretty specific to their purpose, but they also provide some examples of how to use more generalized scripts
# and tools, including generating ETAS forecasts/figures.
#
# some scripts and notes for China_2015 talk(s):

import ETASscripts as esp
import pylab as plt
import os
import datetime as dtm
import pytz
import numpy

import BASScast as bcp
import ANSStools as atp

import matplotlib.dates as mpd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.dates as mpd
import matplotlib.mpl as mpl
#
#import shapely.geometry as sgp
#
from mpl_toolkits.basemap import Basemap as Basemap

# copy some parameters from ETASscripts:
nepal_epi_lon = 84.698
nepal_epi_lat = 28.175
nepal_dlon = 5.
nepal_dlat = 5.
kmldir = 'ACES_2015/kml'
catdir = kmldir
nepal_ETAS_prams = {'todt':None, 'gridsize':.1, 'contres':5, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'nepal', 'catlen':5.0*365.0, 'doplot':False, 'lons':[nepal_epi_lon-nepal_dlon, nepal_epi_lon+nepal_dlon], 'lats':[nepal_epi_lat-nepal_dlat, nepal_epi_lat+nepal_dlat], 'bigquakes':None, 'bigmag':7.00, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}

lats_nepal=[26., 30.]
lons_nepal=[83., 89.]

def my_now():
	return dtm.datetime.now(pytz.timezone('UTC'))

def nepal_etas_set(big_mag=6.28, cat_len = 5.*365., kmldir = 'ACES_2015/png_etas', rtype='ssim_inv_gamma'):
	# make a bunch of ETAS, basically for each earthquake in the Gorkha series. todt = {just before the next earthquake}, start with just-before the mainshock.
	# let's not make a full set of files; just .png... or anyway, don't use the canned function. script .png and maybe .kml here. also script decoration (aka, 
	# plotting earthquakes, etc.).
	# so, call ETAS, returns a BASScast object. then, call nepal_bm_decorator(basscast_object.cm), then save...
	#
	#kmldir = 'ACES_2015/kml_prime'
	#kmldir = 'ACES_2015/png_etas'
	if not os.path.isdir(kmldir): os.makedirs(kmldir)
	#
	mc = nepal_ETAS_prams['mc']
	lons = nepal_ETAS_prams['lons']
	lats = nepal_ETAS_prams['lats']
	#
	lats_map = lats_nepal
	lons_map = lons_nepal
	#
	dt1 = my_now()
	dt0 = dt1-dtm.timedelta(days=cat_len)
	#
	cat_bass =  bcp.getMFETAScatFromANSS(lons=lons, lats=lats, dates=[dt0, my_now()], mc=mc)
	big_quakes = [rw for rw in cat_bass if rw[3]>big_mag]
	#
	# find mainshock:
	mainshock_mag = cat_bass[0][3]
	for j,rw in enumerate(cat_bass):
		if rw[3]>mainshock_mag:
			mainshock=rw
			mainshock_mag=rw[3]
			mainshock_j = j
		#
	#
	print("mainshock info: (%d/%d:%f) %s" % (mainshock_j, len(cat_bass), mainshock_mag, str(mainshock)))
	#
	for j,rw in enumerate(cat_bass[mainshock_j:]):
		B = bcp.BASScast(incat = cat_bass[:j+mainshock_j], fcdate=numpy_date_to_datetime(rw[0]), gridsize=.1, contres=10,mc=mc,eqtheta=None, eqeps=None, fitfactor=5.0, lats=lats, lons=lons, doBASScast=True, rtype=rtype, map_projection='cyl')
		plt.figure(1)
		plt.close()		
		q = B.BASScastContourMap(maxNquakes=0., lats=lats_map, lons=lons_map, fignum=1, fig_size=(12,8))
		#
		# now decorate:
		# (i think this returns the cm object it injests, so we might be careful to not make a weird memory leak.
		#B.cm=nepal_bm_decorator(cm=B.cm, prams = nepal_ETAS_prams, fnum=0, map_res='i', big_mag=6.28, hours_after=None, **kwargs)
		print("hours after: %f - %f = %f :: %f" % (rw[0], mainshock[0],rw[0]-mainshock[0], (rw[0]-mainshock[0])*24.))
		B.cm=nepal_bm_decorator(cm=B.cm, todt=None, hours_after=(rw[0]-mainshock[0])*24.)
		#
		plt.title('Nepal ETAS: %s\n\n' % str(numpy_date_to_datetime(rw[0])))
		fig_name = os.path.join(kmldir, 'nepal_%s_etas_set_%s.png' % (rtype, str('000000%d' % j)[-4:]))
		plt.savefig(fig_name)
		
	
	return big_quakes

def nepal_etas():
	# first, use existing machinery to make a BASScast object complete to today. we can then use this object to make custom forecsts... or maybe we'll just run a
	# massive set of these things and output to a special dir...
	to_dt = dtm.datetime.now(pytz.timezone('UTC'))
	z=esp.make_etas_fcfiles(root_prams=nepal_ETAS_prams, todt=to_dt, rtype='ssim_inv_gamma', fnameroot='nepal_china2015')

	z.BASScastContourMap(maxNquakes=5., lats=lats_nepal, lons=lons_nepal)
	fig_fname = os.path.join(kmldir, 'nepal_inv_gamma_etas_%s.png' % to_dt.isoformat())
	plt.savefig(fig_fname)
	#
	# get mainshock:
	mainshock_mag = max([rw[3] for rw in z.catalog])
	for j,rw in enumerate(z.catalog):
		if rw[3]==mainshock_mag:
			mainshock = rw
			first_aftershock = z.catalog[j+1]
			break
		#
	#
	# so now, we'll run ETAS for just before the first aftershock, then maybe a couple hours, 6 hours, day, then 7 may, etc.
	# show largest events (mainshock, largest aftershocks) for each pane and also aftershocks that have occrred in that time window.
	
	#
	'''
	big_quakes = [rw for rw in z.catalog if rw[3]>=7.3]

	B1=esp.bcp.BASScast(incat=big_quakes, fcdate=esp.dtm.datetime(2015,5,25, tzinfo=esp.pytz.timezone('UTC')), contres=5, mc=5.0, eqtheta=None, eqeps=None, lats=[24., 32.], lons=[81., 90.], rtype='ssim_inv_gamma')

	B1.BASScastContourMap(lats=[26., 30.], lons=[83., 19.])

	B2=esp.bcp.BASScast(incat=big_quakes[0:1], fcdate=esp.dtm.datetime(2015,5,25, tzinfo=esp.pytz.timezone('UTC')), contres=5, mc=5.0, eqtheta=None, eqeps=None, lats=[24., 32.], lons=[81., 90.], rtype='ssim_inv_gamma')
	B2.BASScastContourMap()

	B3=esp.bcp.BASScast(incat=z.catalog, fcdate=esp.dtm.datetime(2015,5,7, tzinfo=esp.pytz.timezone('UTC')), contres=5, mc=5.0, eqtheta=None, eqeps=None, lats=[24., 32.], lons=[81., 90.], rtype='ssim_inv_gamma')
	B3.BASScastContourMap()
	'''
	#
	return z

def nepal_basemap(prams = nepal_ETAS_prams, fnum=0, map_res='i', **kwargs):
	# hours_after: upper time limit for aftershocks to plot.
	prams.update(kwargs)
	#
	lons_nepal = [83., 87.]
	lats_nepal = [26., 30.]
	
	todt=prams.get('todt', None)
	catlen=prams.get('catlen', 5.*365.)
	lons=prams['lons']
	lats=prams['lats']
	mc = prams['mc']
	#
	if todt==None: todt = dtm.datetime.now(pytz.timezone('UTC'))
	dt0 = todt - dtm.timedelta(days=catlen)
	#
	plt.figure(fnum)
	plt.clf()
	ax1=plt.gca()
	cntr = [.5*(lons[0]+lons[1]), .5*(lats[0]+lats[1])]
	#
	cm=Basemap(llcrnrlon=lons_nepal[0], llcrnrlat=lats_nepal[0], urcrnrlon=lons_nepal[1], urcrnrlat=lats_nepal[1], resolution=map_res, projection='cyl', lon_0=cntr[0], lat_0=cntr[1])
	cm.drawcoastlines(color='gray', zorder=1)
	cm.drawcountries(color='gray', zorder=1)
	cm.drawstates(color='gray', zorder=1)
	cm.drawrivers(color='gray', zorder=1)
	cm.fillcontinents(color='beige', zorder=0)
	#
	cm.drawmeridians(list(range(int(lons[0]), int(lons[1]))), color='k', labels=[0,0,1,1])
	cm.drawparallels(list(range(int(lats[0]), int(lats[1]))), color='k', labels=[1, 1, 0, 0])
	#
	return cm
#
def nepal_bm_decorator(cm=None, prams = nepal_ETAS_prams, fnum=0, map_res='i', big_mag=6.28, hours_after=None, **kwargs):
	cm = (cm or nepal_basemap())
	prams.update(kwargs)
	#
	lons_nepal = [83., 87.]
	
	todt=prams.get('todt', None)
	catlen=prams.get('catlen', 5.*365.)
	lons=prams['lons']
	lats=prams['lats']
	mc = prams['mc']
	#
	if todt==None: todt = dtm.datetime.now(pytz.timezone('UTC'))
	dt0 = todt - dtm.timedelta(days=catlen)
	
	# get a catalog:
	cat_0 = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[dt0, todt], fout=None, rec_array=True)
	#print cat_0[0:5]
	#print list(cm(cat_0[0][2], cat_0[0][1]))
	quakes_map = [[rw[0]] + list(cm(rw[2], rw[1])) + [rw[3]] for rw in cat_0]
	
	max_mag = max([m for m in zip(*cat_0)[3]])
	print("max magnitude: ", max_mag)
	for rw in cat_0:
		if rw[3]==max_mag:
			mainshock = rw
			break
	# make a datetime type for mainshock:
	mainshock_dtm = numpy_date_to_datetime(mainshock[0])
	if hours_after!=None:
		time_after = dtm.timedelta(hours=hours_after)
	else:
		time_after = my_now()-mainshock_dtm
	#
	cat_big = [x for x in cat_0 if (x[3]>big_mag and x[0]>mainshock[0] and x!=mainshock)]		# just big aftershocks...
	#
	plt.plot([mainshock[2]], [mainshock[1]], 'r*', ms=17, zorder=6)
	plt.plot([mainshock[2]], [mainshock[1]], 'k*', ms=20, zorder=5)
	plt.plot(*list(zip(*[[rw[1], rw[2]] for rw in quakes_map if rw[0]<mainshock[0]])), marker='o', ms=5, ls='', color='r', zorder=5, label='befores')
	#
	plt.plot(*list(zip(*[[rw[1], rw[2]] for rw in quakes_map if (rw[0]>mainshock[0] and rw[0]<mainshock[0].tolist()+time_after)])), marker='o', ls='', ms=5, color='b', zorder=5, label='afters')
	#
	# ... and the big ones:
	for rw in cat_big:
		print(rw)
		dt=rw[0].tolist()
		#
		plt.plot([rw[2]], [rw[1]], 'o', ms=15.*(rw[3]/8.), label='m=%.2f, %d/%d' % (rw[3], dt.month, dt.day), zorder=6, alpha=.6)
	
	#
	plt.legend(loc=0, numpoints=1)
	
	return cm
#
def numpy_date_to_datetime(numpy_date, tz='UTC'):
	#
	if isinstance(numpy_date, dtm.datetime): return numpy.date
	#
	if isinstance(numpy_date,float): return mpd.num2date(numpy_date)
	#
	return dtm.datetime(*list(numpy_date.tolist().timetuple())[:6] + [numpy_date.tolist().microsecond], tzinfo=pytz.timezone(tz))
	
def nepal_intervals_figure(root_prams=nepal_ETAS_prams, make_figs=False, fig_dir = '/home/myoder/Dropbox/Research/ACES/China2015/talks/nepal', **kwargs):
	# 
	my_prams = root_prams.copy()
	my_prams.update(kwargs)
	#
	todt=root_prams.get('todt', None)
	catlen=root_prams.get('catlen', 5.*365.)
	lats = root_prams['lats']
	lons = root_prams['lons']
	mc = root_prams['mc']
	
	if todt==None: todt = dtm.datetime.now(pytz.timezone('UTC'))
	dt0 = todt - dtm.timedelta(days=catlen)
	#
	cat = bcp.getMFETAScatFromANSS(lons=lons, lats=lats, dates=[dt0, todt], mc=mc, rec_array=True)
	#
	T, dt, M = list(zip(*[[rw[0], rw[0]-cat[j][0], rw[3]] for j,rw in enumerate(cat[1:])]))
	big_mags = [[rw[0], mc, rw[3]] for rw in cat if rw[3]>=6.5]
	#
	mainshock_mag = max(cat['mag'])
	mainshock_index, mainshock = [(j,rw) for j,rw in enumerate(cat) if rw['mag']==mainshock_mag][0]
	print("mainshock: ", mainshock)
	
	f0=plt.figure(0)
	plt.clf()
	#ax1=plt.gca()
	#ax2=ax1.twinx()
	ax1 = f0.add_axes([.1, .5, .85, .45])
	ax2 = f0.add_axes([.1, .05, .85, .45], sharex=ax1)
	#
	ax1.set_ylabel('Intervals $\\Delta t_j$', size=20)
	ax2.set_ylabel('Magnitudes $m$', size=20)
	ax2.set_xlabel('Time $t$', size=20)
	#
	ax1.set_yscale('linear')
	ax1.plot(T,dt, '-', lw=2, label='intervals')
	ax2.vlines(T, [mc for x in M], M, color='g', linestyle='solid', lw=2.5)
	ax2.vlines(*list(zip(*big_mags)), color='r', linestyle='solid', lw=2.5)
	#
	if make_figs:
		plt.savefig(os.path.join(fig_dir, 'nepal_intervals_full.png'))
		#
		ax1.set_xlim(735539.)
		plt.savefig(os.path.join(fig_dir, 'nepal_intervals_close.png'))
		ax1.set_xlim(735711., 735715)
		plt.savefig(os.path.join(fig_dir, 'nepal_intervals_super_close.png'))
		ax1.set_xlim(7.35711e5+2.2, 7.35711e5+2.5)
		plt.savefig(os.path.join(fig_dir, 'nepal_intervals_super_duper_close.png'))
	
	#
	# let's look at an x,y,z t plot...
	#cat_3d = numpy.rec.array(filter(lambda rw: (rw['lon']>83. and rw['lon']<88.) and (rw['lat']>26.5 and rw['lat']<29.5) and rw['event_date']>=(mainshock[0]-1.), cat), dtype=cat.dtype)
	cat_3d = [rw for rw in cat if ((rw['lon']>83. and rw['lon']<88.) and (rw['lat']>26.5 and rw['lat']<29.5) and (rw['event_date']>=(mainshock[0]-1.) and rw['event_date']<mpd.datestr2num('2015-05-14')))]
	#return cat_3d, cat.dtype
	cat_3d = numpy.core.records.fromarrays(list(zip(*cat_3d)), dtype=cat.dtype)
	#
	f3d = plt.figure(1)
	plt.clf()
	ax3d = f3d.add_axes([.1, .1, .8, .8], projection='3d')
	ax3d.plot(cat_3d['lon'], cat_3d['lat'], cat_3d['event_date'], '.')
	ax3d.plot([mainshock['lon']], [mainshock['lat']], [mainshock['event_date']], 'rD', ms=12)
	ax3d.plot(*list(zip(*[[rw['lon'], rw['lat'], rw['event_date']] for rw in cat_3d if rw['mag']>6.28 and rw!=mainshock])), marker='o', color='m', ms=7)
	#
	ax3d.set_xlabel('longitude', size=18)
	ax3d.set_ylabel('latitude', size=18)
	ax3d.set_zlabel('time $t$ $years$\n', size=18)
	#
	#return cat

if __name__=='__main__':
	# running as main (called form terminal) script.
	z=doit()
else:
	plt.ion()
	

