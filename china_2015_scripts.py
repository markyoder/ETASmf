# some scripts and notes for China_2015 talk(s):

import ETASscripts as esp
import pylab as plt
import os
import datetime as dtm
import pytz
import BASScast as bcp
import ANSStools as atp

import matplotlib.dates as mpd
import matplotlib
import matplotlib.pyplot as plt
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
	cm.drawmeridians(range(int(lons[0]), int(lons[1])), color='k', labels=[0,0,1,1])
	cm.drawparallels(range(int(lats[0]), int(lats[1])), color='k', labels=[1, 1, 0, 0])
	#
	# get a catalog:
	cat_0 = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[dt0, todt], fout=None, rec_array=True)
	#print cat_0[0:5]
	#print list(cm(cat_0[0][2], cat_0[0][1]))
	quakes_map = [[rw[0]] + list(cm(rw[2], rw[1])) + [rw[3]] for rw in cat_0]
	
	max_mag = max([m for m in zip(*cat_0)[3]])
	print "max magnitude: ", max_mag
	for rw in cat_0:
		if rw[3]==max_mag:
			mainshock = rw
			break
	#
	plt.plot([mainshock[2]], [mainshock[1]], 'r*', ms=17, zorder=6)
	plt.plot([mainshock[2]], [mainshock[1]], 'k*', ms=20, zorder=5)
	plt.plot(*zip(*[[rw[1], rw[2]] for rw in quakes_map if rw[0]<mainshock[0]]), marker='o', ms=5, ls='', color='r', zorder=5, label='befores')
	time_after = dtm.timedelta(hours=36)
	plt.plot(*zip(*[[rw[1], rw[2]] for rw in quakes_map if (rw[0]>mainshock[0] and rw[0]<mainshock[0].tolist()+time_after)]), marker='o', ms=5, ls='', color='b', zorder=5, label='afters')
	#
	plt.legend(loc=0, numpoints=1)
	
	return quakes_map

def nepal_intervals_figure(root_prams=nepal_ETAS_prams, **kwargs):
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
	cat = bcp.getMFETAScatFromANSS(lons=lons, lats=lats, dates=[dt0, todt], mc=mc)
	#
	T, dt, M = zip(*[[rw[0], rw[0]-cat[j][0], rw[3]] for j,rw in enumerate(cat[1:])])
	big_mags = [[rw[0], mc, rw[3]] for rw in cat if rw[3]>=6.5]
	
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
	ax2.vlines(*zip(*big_mags), color='r', linestyle='solid', lw=2.5)
	#
	fig_dir = '/home/myoder/Dropbox/Research/ACES/China2015/talks/nepal'
	plt.savefig(os.path.join(fig_dir, 'nepal_intervals_ful.png'))
	ax1.set_xlim(735539.)
	plt.savefig(os.path.join(fig_dir, 'nepal_intervals_close.png'))
	ax1.set_xlim(735711., 735715)
	plt.savefig(os.path.join(fig_dir, 'nepal_intervals_super_close.png'))
	ax1.set_xlim(7.35711e5+2.2, 7.35711e5+2.5)
	plt.savefig(os.path.join(fig_dir, 'nepal_intervals_super_duper_close.png'))
	
	#
	#return cat

if __name__=='__main__':
	# running as main (called form terminal) script.
	z=doit()
else:
	plt.ion()
	

