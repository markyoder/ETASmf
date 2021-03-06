#!/usr/bin/python
#
# ETASscripts.py:
#
# Primary Author: Mark R. Yoder, Ph.D.
#                 mryoder@ucdavis.edu
#                 mark.yoder@gmail.com
#
# This is the main, go-to script for ETAS map generation. The core ETAS module is BASScast.py (aka, BASS-forecast); in this module, we have scripts
# to gather specific catalogs, produce ETAS forecasts, plot the forecast on maps, decorate maps, generate output files (figures, .kml, .csv, etc.).
# this includes some general tools; see epecially:
# - etas_auto()
# - makeETASFCfiles()
# - make_etas_fcfiles()  (which is basically a wrapper for makeETASFCfiles() with different/generalized input formatting.
#
# then, see a bunch of specific scripts (sadly, many evolving generations of them) that call some of these functions. the idea, in general, is to develop
# scripts that are easy to use (aka, a single command-line call), but that don't end up being a nasty spaghetti code or requiring too much cut-n-pasting.
# unfortunately, there does not appear to be a magic formula for this, but we seem to be closing in on something.
# 
#
import matplotlib as mpl
mpl.use('Agg')
#
import BASScast as bcp
import ANSStools as atp     # a yoder et al. script in my "common" git repository. simple tools to fetch catalogs from ANSS.
#
import datetime as dtm
import pytz
import os
import math
import glob
import numpy
#
import matplotlib.pyplot as plt
import matplotlib.dates as mpd
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import axes3d
#
import multiprocessing as mpp
import pickle

kmldir='kml'
catdir='kml'

lon2km=111.1
deg2rad = math.pi*2.0/360.
#
# some parameter sets for important earthquakes:
#
tohoku_ETAS_prams = {'todt':dtm.datetime.now(pytz.timezone('UTC')), 'gridsize':.1, 'contres':3, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'tohoku', 'catlen':5.0*365.0, 'doplot':False, 'lons':[135., 146.], 'lats':[30., 41.5], 'bigquakes':None, 'bigmag':7.50, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}
#
chengdu_ETAS_prams = {'todt':dtm.datetime.now(pytz.timezone('UTC')), 'gridsize':.1, 'contres':5, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'chengdu', 'catlen':5.0*365.0, 'doplot':False, 'lons':[100.367, 106.367], 'lats':[31.06-3., 31.06+3.], 'bigquakes':None, 'bigmag':6.7, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}

#'lat_center':31.021, 'lon_center':103.367
#
nepal_epi_lon = 84.698
nepal_epi_lat = 28.175
nepal_dlon = 5.
nepal_dlat = 5.
nepal_ETAS_prams = {'todt':None, 'gridsize':.1, 'contres':5, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'nepal', 'catlen':5.0*365.0, 'doplot':False, 'lons':[nepal_epi_lon-nepal_dlon, nepal_epi_lon+nepal_dlon], 'lats':[nepal_epi_lat-nepal_dlat, nepal_epi_lat+nepal_dlat], 'bigquakes':None, 'bigmag':7.00, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}#
#
#
def japan_pacific_2015_m78():
	E=esp.makeETASFCfiles(todt=esp.dtm.datetime.now(esp.pytz.timezone('UTC')), contres=8, mc=4.0, bigmag=6.5, lons=[140.5-5., 150.5+5.], lats=[27.83-5., 27.83+5.])
	mymapy = E.BASScastContourMap(maxNquakes=3)
	
	return E
#	
# Some scripts to run some standard california ETAS sets.
def makeCaliFiles(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, fcmc=2.0, kmldir='kml', catdir='kml', catlen=5.0*365.0):
	norcalkml = kmldir + '/norcalconts.kml'
	socalkml = kmldir + '/socalconts.kml'
	calnevkml = kmldir+ '/calnevconts.kml'
	#
	norcalcat = catdir + '/norcal.cat'
	socalcat = catdir + '/socal.cat'
	calnevcat = catdir + '/calnev.cat'
	#
	print("make parkfield.")
	c=makeParkfieldETAS()
	print("make el mayor.")
	d=makeElMayorETAS()
	print("make socal")
	a=makeSocalETAS()
	print("make norcal")
	b=makeNorcalETAS()
	#
	print("finished.")
	#
	return None

def pfmovie(dtstart=dtm.datetime(2002,1,1, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC')), dtend=dtm.datetime(2004,12,31, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC')), dt=10, outdir='pfmovie'):
	'''
	# make a bunch of parkfield ETAS frames; put them together into a movie.
	'''
	mylons = [-120.75, -119.5]	# -120.75
	mylats = [35.5, 36.5]		# 35.5
	#
	gridsize=.1
	contres=3
	mc=2.0
	fout=open('%s/movielist.lst' % outdir, 'w')
	fout.close()
	#
	#cat0 = getMFETAScatFromANSS(lons=mylons, lats=mylats, dates=[dtstart-dtm.timedelta(days=365*5), dtend], mc=3.0)
	#cat=[cat0.pop(0)]	# a smarter way to do this is to just use an index and pass a subset of the array.
	deltat=dtm.timedelta(days=dt)
	catdt=dtm.timedelta(days=365*5)
	#
	fstart=bcp.mpd.date2num(dtstart)
	fstop=bcp.mpd.date2num(dtend)
	fdt=float(dt)	# in days...
	#	
	#thisdt=fstart
	thisdt = dtstart
	framecount=0
	#while thisdt<=fstop:
	while thisdt<dtend:
		#while cat0[0][0]<mpd.date2num(thisdt):
		#	cat+=[cat0.pop(0)]
		cat = bcp.getMFETAScatFromANSS(lons=mylons, lats=mylats, dates=[thisdt-catdt, thisdt], mc=3.0)
		obj=bcp.BASScast(incat=cat, fcdate=thisdt, gridsize=gridsize, contres=contres, mc=mc)
		# or, we could instantiate the object once, then add to its catalog and run calcBASS()...
		#
		#thisdt=thisdt+fdt
		thisdt+=deltat
		framecount+=1
		#
		plt.figure(1)
		plt.clf()
		obj.BASScastContourMap()
		# parkfield, san sim, coalinga
		x,y=obj.cm([-120.37, -121.10, -120.32], [35.82, 35.71, 36.23] )		
		plt.plot([x], [y], '*', ms=15, alpha=.7)
		plt.title('date: %s' % thisdt)
		#
		strIndex='00000000%d' % framecount
		strIndex=strIndex[-6:]
		framename='%s/pf-%s.png' % (outdir, strIndex)
		framejpg ='%s/pf-%s.jpg' % (outdir, strIndex)
		plt.savefig(framename)
		os.system('convert %s %s' % (framename, framejpg))
		#
		fout=open('%s/movielist.lst' % outdir, 'a')
		fout.write('%s\n' % framejpg)
		fout.close()
	#
	#os.system('pfmovie$ mencoder mf://*.jpg -mf w=800:h=600:fps=5:type=jpg -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o pfmovie.avi')

#
def etas_auto(lon_center=None, lat_center=None, d_lat_0=.25, d_lon_0=.5, dt_0=10, Lr_factor=4.0, mc=2.5, mc_0=4.5, gridsize=.1, to_dt=None, fnameroot='etas_auto', catlen=5.0*365.0, doplot=False, kmldir='kml_auto', **kwargs):
	'''
	# this is probably a major staple of future work in ETAS; when in doubt, start here and point it at a big earthquake.
	#
	# a starter script to auto-select some parameters for an ETAS run. in practice, give this script a center location (probably a mainshock epicenter).
	# the script will find the largest earthquake in that region and scale up an ETAS parameter set accordingly.
		# d_lat_0, d_lon_0, dt_0 are the starting catalog parameters (largest earthquake in d_lat_0 x d_lon_0 x dt_0 cube).
	'''
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
	#
	#biggest_earthquake = filter(lambda x: x['mag']==max(cat_0['mag']), cat_0)[0]
	mainshock = {cat_0.dtype.names[j]:x for j,x in enumerate(filter(lambda x: x['mag']==max(cat_0['mag']), cat_0)[0])}
	#
	# now, get new map domain based on rupture length, etc.
	L_r = .5*10.0**(.5*mainshock['mag'] - 1.76)
	d_lat = Lr_factor*L_r/lon2km
	d_lon = Lr_factor*L_r/(lon2km*math.cos(deg2rad*mainshock['lat']))
	print("mainshock data: ", mainshock, L_r, d_lat, d_lon)
	#
	working_cat = atp.catfromANSS(lon=[mainshock['lon']-d_lon, mainshock['lon']+d_lon], lat=[mainshock['lat']-d_lat, mainshock['lat']+d_lat], minMag=mc, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)
	#
	print("biggest event(s): ", [rw for rw in working_cat if rw['mag']==max(working_cat['mag'])])
	#
	# now, do some ETAS:
	# skip working_cat above, but parse lon, lat, etc. parameters similarly. pass those (and other) params to make_etas_fcfiles()
	# looks like this: make_etas_fcfiles(root_prams=nepal_ETAS_prams, **kwargs), and:
	# nepal_ETAS_prams = {'todt':None, 'gridsize':.1, 'contres':5, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'nepal', 'catlen':5.0*365.0, 'doplot':False, 'lons':[nepal_epi_lon-nepal_dlon, nepal_epi_lon+nepal_dlon], 'lats':[nepal_epi_lat-nepal_dlat, nepal_epi_lat+nepal_dlat], 'bigquakes':None, 'bigmag':7.00, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}
	root_prams = {'todt':None, 'gridsize':gridsize, 'contres':10, 'mc':mc, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':fnameroot, 'catlen':catlen, 'doplot':False, 'lons':[mainshock['lon']-d_lon, mainshock['lon']+d_lon], 'lats':[mainshock['lat']-d_lat, mainshock['lat']+d_lat], 'bigquakes':None, 'bigmag':mainshock['mag']-1.5, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}
	#root_prams = {}
	print("now execute with root_prams: ", root_prams)
	#my_kwargs = {}
	etas = make_etas_fcfiles(root_prams=root_prams, **kwargs)
	plt.figure(1)
	plt.title('ETAS to: %s' % str(to_dt if to_dt!=None else dtm.datetime.now(pytz.timezone('UTC'))))
	plt.savefig(os.path.join(kmldir, fnameroot+'.png'))
	#
	#return biggest_earthquake
	return etas
#
def contour_3d_etas(etas_object, fignum=0):
    # make a 3d ETAS figure (aka, a 3d surface showing ETAS rates).
	# take a BASScast object as an input (for now).
	#
	X0,Y0,Z = etas_object.X_i, etas_object.Y_i, etas_object.Z2d
	#
	min_z = min([min(rw) for rw in Z])
	max_z = max([max(rw) for rw in Z])
	#
	# get 2d mesh-arrays:
	X,Y = numpy.meshgrid(X0,Y0)
	#
	# find and plot mainshock:
	ms_mag = max([rw[3] for rw in etas_object.catalog])
	mainshock = [rw for rw in etas_object.catalog if rw[3]==ms_mag]
	chengdu = {'lon':104.+4./60., 'lat':30. + 2./3.}
	# 
	# find z-value:
	#delta_xs = sorted([[j,(mainshock[2]-x)**2.] for j,x in enumerate(X0)], key=lambda rw:rw[1])
	#delta_ys = sorted([[j,(mainshock[1]-y)*2.] for j,y in enumerate(Y0)], key=lambda rw:rw[1])
	delta_xs = sorted([[j,(chengdu['lon']-x)**2.] for j,x in enumerate(X0)], key=lambda rw:rw[1])
	delta_ys = sorted([[j,(chengdu['lat']-y)**2.] for j,y in enumerate(Y0)], key=lambda rw:rw[1])
	j_ms = delta_xs[0][0]
	k_ms = delta_ys[0][0]
	#z_ms = delta_ys[k_ms][j_ms]
	z_chengdu = Z[k_ms][j_ms] + .5
	#print "xy: ", X[j_ms][k_ms], Y[j_ms][k_ms], " :: ", X[k_ms][j_ms], Y[k_ms][j_ms]
	#
	fg=plt.figure(fignum)
	plt.clf()
	ax=fg.gca(projection='3d')
	mycm = cmx.spectral
	ax.plot_surface(X,Y,Z, rstride=1, cstride=1, alpha=.6, cmap=mycm)
	#ax.plot([chengdu['lon']], [chengdu['lat']], [z_chengdu], marker='o', color='r', ms=12)
	#ax.plot([chengdu['lon']], [chengdu['lat']], [-12.], marker='s', color='r', ms=12)
	
	cset = ax.contourf(X, Y, Z, zdir='z', offset=min_z-0.*.25*(max_z-min_z), cmap=mycm)
	#cset = ax.contourf(X, Y, Z, zdir='x', offset=min(X0), cmap=mycm)
	#cset = ax.contourf(X, Y, Z, zdir='y', offset=max(Y0), cmap=mycm)
	
	#cset = ax.contourf(X, Y, Z, zdir='z', cmap=cmx.coolwarm)
	#cset = ax.contourf(X, Y, Z, zdir='x', cmap=cmx.coolwarm)
	#cset = ax.contourf(X, Y, Z, zdir='y', cmap=cmx.coolwarm)


	ax.set_xlabel('lon')
	#ax.set_xlim(-40, 40)
	ax.set_ylabel('lat')
	#ax.set_ylim(-40, 40)
	ax.set_zlabel('$\\log(rate)$')
	#ax.set_zlim(-100, 100)
#
def make_etas_fcfiles(root_prams=nepal_ETAS_prams, **kwargs):
	'''
	# a wrapper for fast etas code development.
	# wraps makeETASFCfiles() but takes generalized input format.
	'''
	if not isinstance(root_prams, dict): root_prams = {}
	#
	my_prams = root_prams.copy()
	my_prams.update(kwargs)
	#
	return makeETASFCfiles(**my_prams)
#
def circle_poly(x=0., y=0., R=1.0, n_points=100):
	'''
	# make a circle-polygon witn n_points points. return array like [[x0,y1], [x1,y1], ...] in proper order. vectors can be constructedy by the consumer.
	'''
	#
	d_theta = math.pi*2.0/float(n_points)
	poly = []
	#
	for n in range(n_points):
		theta = n*d_theta
		poly += [[x+R*math.cos(theta), y+R*math.sin(theta)]]
	#
	return poly + [poly[0]]
#
def makeETASFCfiles(todt=dtm.datetime(2004, 9, 15, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC')), gridsize=.1, contres=3, mc=1.5, kmldir='kml', catdir='kml', fnameroot='parkfield', catlen=5.0*365.0, doplot=False, lons=[-120.75, -119.5], lats=[35.75, 36.5], bigquakes=[], bigmag=3.5, addquakes=[], eqeps=None, eqtheta=None, fitfactor=5.0, cmfnum=0, fignum=1, colorbar_fontcolor='k', contour_intervals=None, rtype='ssim', contour_top=1.0, contour_bottom=0.0, p_quakes=None, p_map=None, fnameroot_suffix='', maxNquakes=None, **kwargs):
	# general script to make ETAS forecast files (default values -> parkfield).
	# this is the real go-to script. this takes in a bunch of parameters (where, what range, when, how many contours, etc.), generates a BASScast forecast,
	# then generates some standard output files.
	#
	if todt==None: todt=dtm.datetime.now(pytz.timezone('UTC'))
	maxNquakes == int((1 or maxNquakes))
	#
	if kmldir==None and catdir==None:
		kmldir, fnameroot = os.path.split(fnameroot)
		catdir=kmldir
	if kmldir==None: kmldir=catdir
	if catdir==None: catdir=kmldir
		
	if not os.path.isdir(kmldir): os.makedirs(kmldir)
	if not os.path.isdir(catdir): os.makedirs(catdir)
	
	dt1=todt
	dt0=dt1-dtm.timedelta(days=catlen)
	fnameroot=fnameroot+str(fnameroot_suffix)
	
	print("fnameroot: ", fnameroot)
	#return None
	
	#bigquakes=[]
	#print dt0, dt1
	if bigquakes==None: bigquakes=[]
	#
	figname = kmldir + '/' + fnameroot + '_fig.png'
	kmlfile = kmldir + '/' + fnameroot + 'conts.kml'
	kmzfile = kmldir + '/' + fnameroot + '.kmz'
	catfile = kmldir + '/' + fnameroot + '-BASS' + '.cat'
	xyzfile = kmldir + '/' + fnameroot + '-ETASgrid' + '.xyz'
	cont_shapes_file = kmldir + '/' + fnameroot + '-contshapes' + '.txt'
	cbarname = '%sscale.png' % fnameroot
	cbarhname = '%sscale-h.png' % fnameroot
	
	cat=[]
	cat = bcp.getMFETAScatFromANSS(lons=lons, lats=lats, dates=[dt0, dt1], mc=mc)
	# note: this is likely now returning a recarray and the dates will be numpy.datetime64 types, not datetime.datetime.
	#
	# bandaid:
	# this query pulls [float-date, lat, lon, mag, depth?].
	# BUT, there's a bit that checks for epsilon, theta orientation in the ETAS script, and of course
	# it expects these to be in rw[4], rw[5]. this query puts depth in rw[4]... so we get a mess.
	# for now, truncate depth. we'll fix this later by converting the catalog rows to dictionary objects... or by finding the scripts that apply
	# epsilon, theta and correcting them to use at least [5],[6]... or we just add a dict. object as an item in the list?
	# ... you know what, screw it. let's just call the correct cols here. we'll fix the scenario later if need be.
	#cat = map(operator.itemgetter(0,1,2,3), cat)
	# note: a newer version of getMFETAScatFromANSS() returns a recarray, so we can be more explicit about which columns we use.
	#
	if len(cat)==0:
		print("no available data.")
		return None
	# returns "catalog" like: catalog+=[[mpd.date2num(rw[0]), rw[1], rw[2], rw[3], rw[4]]] (aka, a list of [float-dtm, lat, lon, mag, dpth]
	#
	# addquakes: do we want to add any other earthquakes into the catalog before we calculate ETAS?
	print("addquakes: %d" % len(addquakes))
	if len(addquakes)>0:
		for qk in addquakes:
			# using "hasattr()" would be a better way of determining data type (date or float).
			try:
				dtyr = qk[0].year
			except:
				# not a datetime type object?
				# good. it's probably a float.
				# (and be more robust with this later)
				fdt=float(qk[0])
			else:
				fdt=mpd.date2num(qk[0])
			#
			#
			cat+=[[fdt] + qk[1:]]
		#
		# and sort it by time (for now):
		cat.sort(key=lambda x: x[0])
	#
	mfb=None
	mfb=bcp.BASScast(incat=cat, fcdate=dt1, gridsize=gridsize, contres=contres, mc=mc, eqeps=eqeps, eqtheta=eqtheta, fitfactor=fitfactor, contour_intervals=contour_intervals, lats=lats, lons=lons, rtype=rtype, p_quakes=p_quakes, p_map=p_map)
	print("quakeLen: %d, %d" % (len(mfb.quakes),len(cat)))
	#
	#cbs=mfb.makeColorbar(cset=mfb.conts, colorbarname='%s/%s' % (kmldir, cbarname))
	mfb.writeCat(fout=catfile)
	#
	for rw in cat:
		#print "bigmag check: ", rw[0], mpd.date2num(dt1)
		if rw[0]>mpd.date2num(dt1): continue
		if rw[3]>=bigmag:
			bigquakes +=[[rw[0], rw[1], rw[2], rw[3]]]
			#print rw
	print("mfb.conts: ", mfb.conts)
	mfb.writeKMLfile(cset=mfb.conts, fout=kmlfile, colorbarname=cbarname, equakes=bigquakes, top=contour_top, bottom=contour_bottom)
	#str_cont_shapes = mfb.contsToStr(top=contour_top, bottom=contour_bottom)
	#str_cont_shapes = mfb.contsToStr(top=1.0, bottom=0.0)
	#f_shapes = open(cont_shapes_file, 'w')
	#f_shapes.write(str_cont_shapes)
	#f_shapes.close()
	str_cont_shapes = mfb.fixedContsToFile(cset=None, file_out=cont_shapes_file, top=1.0, bottom=0.0)	# and users can easily select conts.
	mfb.xyztofile(outfile=xyzfile)
	
	print("len(bigquakes): %d" % len(bigquakes))
	#
	os.system('zip -q %s %s %s' % (kmzfile, kmlfile, cbarname))
	#cbs=mfb.makeColorbar(cset=mfb.conts, colorbarname='%s/%s' % (kmldir, cbarname), reffMag=bigmag)
	cbs=mfb.makeColorbar(colorbarname='%s/%s' % (kmldir, cbarname), reffMag=bigmag, fnum=cmfnum, fontcolor=colorbar_fontcolor)
	#
	if doplot:
		plt.figure(fignum)
		plt.clf()
		mfb.BASScastContourMap(maxNquakes=maxNquakes)
		x,y=mfb.cm(bigquakes[0][2], bigquakes[0][1])
		plt.plot([x], [y], '*', ms=15, alpha=.7)
		plt.savefig(figname)
		#
	#
	return mfb
#
def make_capstone_2015_scenario(fnum=0, fnameroot = 'capstone/capstone', gridsize=.1, contres=3, mc=3.0, catlen=3.0*365.0, doplot=False, lons=[-119.5, -114.5], lats=[32.0, 35.5], bigquakes=[], bigmag=6.0, scalemag=5.0, addquakes=[], eqeps=None, eqtheta=None, rtype='ssim', maxNquakes=0):
	'''
	# 2015 capstone ETAS script:
	# scenario: mainshock: m=7.8, lat=33.35, lon=-115.71, 5:32 am PDT, 11 May 2015.
 	# rupture_SE = [-115.71, 33.35] (epicenter), rupture_NW = [-118.51, 34.72] (endpoint linear approximation).
 	# get ETAS for 1 hour, 6 hours, 12 hours after mainshock. since there are no aftershocks in the scenario, this will nominally be a silly exercise,
 	# and we might otherwise pick a delta_t to normalize the aftershock halo to something reasonable.
	'''
	#
	# check the output directory:
	# if not os.path.isdir(output_dir):
	#	os.makedirs(output_dir)
	out_dir, out_froot = os.path.split(fnameroot)
	if out_dir == '': out_dir='.'
	if not os.path.isdir(out_dir): os.makedirs(out_dir)	
	kmldir = out_dir
	catdir = out_dir
	fnameroot = out_froot
	#
	mainshock = {'lat': 33.35, 'lon':-115.71, 'depth':7.6, 'mag':7.8, 'event_time':dtm.datetime(2015, 5, 11, 5, 32, tzinfo=pytz.timezone('US/Pacific')), 'rupture_lat_start': 33.35, 'rupture_lon_start':-115.71, 'rupture_lat_end':34.72, 'rupture_lon_end':-118.51}
	#
	# do a sloppy job of estimating the mainshock rupture angle (be careful of the phase. i think our non-standard convention is CCW from E).
	dy = mainshock['rupture_lat_end']-mainshock['rupture_lat_start']
	dx = math.cos(math.pi*2.0*mainshock['rupture_lat_start']/360.)*(mainshock['rupture_lon_end']-mainshock['rupture_lon_start'])
	L_r = math.sqrt(dx**2. + dy**2.)
	theta_mainshock =  math.atan(dy/dx)*360./(2.0*math.pi)
	#
	epsilon_mainshock = 2.0
	#theta_mainshock = 15.
	#
	t0 = mainshock['event_time']
	time_points = [t0 - dtm.timedelta(days=1), t0 + dtm.timedelta(hours=1), t0+dtm.timedelta(days=1), t0+dtm.timedelta(days=2), t0 + dtm.timedelta(days=3), t0 + dtm.timedelta(days=4), t0 + dtm.timedelta(days=5)]
	#
	addquakes=[]		# add just the mainshock; no aftershocks in this scenario. eventually, this will be redone in a smarter way. for now, [dt, lat, lon, mag, depth]
						#    ... and assume that we identified the rupture, and the aftershock hazard is centered on the rupture center, not the epicenter:
	addquakes+=[[mainshock['event_time'], mainshock['lat'] + .5*dy, mainshock['lon'] + .6*dx, mainshock['mag'], mainshock['depth'], epsilon_mainshock, theta_mainshock]]
	#
	print("begin ETAS scenario.")
	print("mainshock: ", mainshock)
	print(" theta: %f, epsilon: %f" % (theta_mainshock, epsilon_mainshock))
	print("rupture: ", dx, dy)
	print("L_r: ", L_r * 111.)
	print("time_points: ")
	for tp in time_points: print(tp)
	
	#
	fn=fnum
	#
	for dfn, fcdt in enumerate(time_points):
		for j in range(2):
			plt.figure(j)
			plt.clf()
		#
		print("forecastdate: %s" % fcdt)
		#
		thisfnameroot = fnameroot + '_%d' % dfn
		bc1=None
		bc1 = makeETASFCfiles(todt=fcdt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=thisfnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes[:], bigmag=scalemag, addquakes=addquakes, eqeps=eqeps, eqtheta=eqtheta, cmfnum=0, fignum=1, colorbar_fontcolor='w', rtype=rtype)
		#
		mfb=bc1.BASScastContourMap(fignum=1, maxNquakes=maxNquakes)
		# and we'll want to draw the fault here as well...
		#Fault endpoints: 36.82 -121.51, 40.25 -124.41
		x1,y1=bc1.cm(mainshock['rupture_lon_start'], mainshock['rupture_lat_start'])
		x2, y2=bc1.cm(mainshock['rupture_lon_end'], mainshock['rupture_lat_end'])
		#
		xms, yms=bc1.cm(mainshock['lon'], mainshock['lat'])
		#xas1, yas1=bc1.cm(addquakes[1][2], addquakes[1][1])
		#xas2, yas2=bc1.cm(as2lon, as2lat)
		plt.plot([x1, x2], [y1, y2], 'm-s', lw=3, alpha=.7, label='Mainshock rupture (approximate)')
		plt.plot([xms], [yms], 'k*', alpha=.8, zorder=7, ms=20)
		plt.plot([xms], [yms], 'r*', alpha=.9, zorder=9, ms=18)
		#
		
		plt.title('fcdate: %s, m>%.2f\n\n\n' % (fcdt, mc))
		plt.savefig('%s/%s_etas_contours-%d.png' % (kmldir,fnameroot, dfn))
		
	
	return bc1
	
	return None
#
def makeGG2013(todt=dtm.datetime(2013, 5, 13, tzinfo=pytz.timezone('US/Pacific-New')), gridsize=.1, contres=3, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='goldenguardian2013', catlen=3.0*365.0, doplot=False, lons=[-126.5, -118.5], lats=[34.0, 43.5], bigquakes=[], bigmag=6.0, scalemag=5.0, addquakes=[], eqeps=1.414, eqtheta=None, rtype='ssim'):
	#
	# eqeps normally
	#dt2=dtm.datetime.now(pytz.timezone('UTC'))
	dtms  = dtm.datetime(2013, 5, 13, 11, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
	dtas1 = dtm.datetime(2013, 5, 15, 12, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
	dtas2 = dtm.datetime(2013, 5, 20, 12, 0,0,0,tzinfo=pytz.timezone('US/Pacific-New'))
	as2lon, as2lat = -122.08, 38.04
	#
	if addquakes==None:
		addquakes=[]
	if addquakes==[]:
		#dt1=dtm.datetime(2013, 5, 13, 11, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
		#dt2=dtm.datetime(2013, 5, 15, 12, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
		#
		addquakes+=[[mpd.date2num(dtms), 37.7, -122.5, 7.9, 15]]			# epicenter (rupture initiation)
		#addquakes+=[[mpd.date2num(dt1), 38.265, -122.96, 7.9, 15]]		# rupture center
		addquakes+=[[mpd.date2num(dtas1), 37.99, -122.59, 6.5, 10.0]]
		# green valley fault scenario:
		addquakes+=[[mpd.date2num(dtas2), as2lat, as2lon, 6.9, 10.0]]
		#
	print("addquakes:", addquakes)
	if bigquakes==None:
		bigquakes=[]
	if dtas2>todt: todt=dtas2
	#
	faultSouth=[-121.51, 36.28]
	faultNorth=[-124.41, 40.25]
	dy=faultNorth[1]-faultSouth[1]
	dx=(faultNorth[0]-faultSouth[0])*math.cos(2.0*math.pi*faultSouth[1]/360.)
	#
	faultTheta = math.atan(dy/dx)*360./6.282	# though this might be negative or should be the other angle (it should be obvious).
	print("fault theta: %f degrees" % faultTheta)
	#faultTheta = (90.+60.)
	#
	fcdate0=dtm.datetime(2013, 5, 13, 00, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))	# nominally, a "standard" time.
	#
	fcdates = [fcdate0] 
	#fcdates0=[dtm.datetime(2013, 5, 13, 00, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))]
	while (fcdates[-1]<=dtm.datetime(2013, 5, 22, 00, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))):
		fcdates+=[fcdates[-1]+dtm.timedelta(hours=12.0)]
		#print "adding fcdate: %s" % str(fcdates[-1])
	# and special event driven fcs:
	fcdates+=[dtms+dtm.timedelta(hours=2.0)]
	fcdates+=[dtas1+dtm.timedelta(hours=1.0)]
	fcdates.sort()
	for rw in fcdates:
		print("adding fcdate: %s" % str(rw))
	fcindex=0
	#
	dome=1
	plt.close(0)
	plt.close(1)
	for fcdt in fcdates:
	#while fcdt<=dtm.datetime(2013, 5, 17, 0, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New')):
	#while dome==1:
		#
		#dome=0
		#fcdt=dtas1 - dtm.timedelta(hours=2.0)
		#fcdt=dtm.datetime(2013, 5, 15, 11, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
		#
		#plt.figure(num=fcindex, figsize=(7,10))
		#plt.close(0)
		#plt.close(1)
		plt.figure(num=1, figsize=(8,10))
		plt.clf()
		#
		print("forecastdate: %s" % fcdt)
		#
		activeEps=eqeps
		activeTheta=eqtheta
		if fcdt>dtms: 
			activeTheta=faultTheta
			activeEps=2.0
		#
		thisfnameroot = fnameroot + '-%d' % fcindex
		bc1=None
		bc1 = makeETASFCfiles(todt=fcdt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=thisfnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes[:], bigmag=scalemag, addquakes=addquakes, eqeps=activeEps, eqtheta=activeTheta, cmfnum=0, fignum=1, colorbar_fontcolor='w', rtype=rtype)
		#
		mfb=bc1.BASScastContourMap(fignum=1)
		# and we'll want to draw the fault here as well...
		#Fault endpoints: 36.82 -121.51, 40.25 -124.41
		x1,y1=bc1.cm(-121.51, 36.82)
		x2, y2=bc1.cm(-124.41, 40.25)
		#xms, yms=bc1.cm(addquakes[0][2], addquakes[0][1])
		xms, yms=bc1.cm(-122.5, 37.7)
		xas1, yas1=bc1.cm(addquakes[1][2], addquakes[1][1])
		xas2, yas2=bc1.cm(as2lon, as2lat)
		plt.plot([x1, x2], [y1, y2], 'm-s', lw=3, alpha=.7)
		#
		if fcdt>=dtms: plt.plot([xms], [yms], 'b*', alpha=.8, ms=2.5*addquakes[0][3])
		if fcdt>=dtas1: plt.plot([xas1], [yas1], 'g*', alpha=.8, ms=2.5*addquakes[1][3])
		if fcdt>=dtas2: plt.plot([xas1], [yas1], 'c*', alpha=.8, ms=2.5*addquakes[1][3])
		plt.title('fcdate: %s, m>%.2f\n\n\n' % (fcdt, mc))
		plt.savefig('%s/%sconts-%d.png' % (kmldir,fnameroot, fcindex))
		fcindex+=1
	
	return bc1

#
def makeShakeout1(todt=dtm.datetime(2012, 10, 19, tzinfo=pytz.timezone('UTC')), gridsize=.1, contres=3, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='shakeout1', catlen=3.0*365.0, doplot=False, lons=[-122.5, -114.0], lats=[31.0, 37.25], bigquakes=None, bigmag=6.0, scalemag=4.0, addquakes=[], eqeps=4.0, eqtheta=-40.0):
	#
	#dt2=dtm.datetime.now(pytz.timezone('UTC'))
	if addquakes==None:
		addquakes=[]
	if addquakes==[]:
		dt1=dtm.datetime(2012, 10, 18, 11, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
		dt2=dtm.datetime(2012, 10, 18, 12, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
		addquakes+=[[mpd.date2num(dt1), 33.35, -115.71, 7.8, 7.6]]
		addquakes+=[[mpd.date2num(dt2), 33.55, -115.89, 7.0, 6.2]]
		#
	print("addquakes:", addquakes)
	if bigquakes==None:
		bigquakes=[]
	if dt2>todt: todt=dt2
	#
	# issue maps for these date-times:
	#fcdates = [dtm.datetime(2012, 10, 18, 8, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New')), dtm.datetime(2012, 10, 18, 11, 10, 0, 0, tzinfo=pytz.timezone('US/Pacific-New')), dtm.datetime(2012, 10, 18, 11, 40, 0, 0, tzinfo=pytz.timezone('US/Pacific-New')), dtm.datetime(2012, 10, 18, 12, 10, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))]
	fcdates = [dtm.datetime(2012, 10, 18, 11, 10, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))]
	#
	fcindex=0
	for fcdt in fcdates:
		plt.figure(fcindex)
		plt.clf()
		#
		thisfnameroot = fnameroot + '-%d' % fcindex
		bc1 = makeETASFCfiles(todt=fcdt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=thisfnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, addquakes=addquakes, eqeps=eqeps, eqtheta=eqtheta)
		#
		bc1.BASScastContourMap(fignum=fcindex)
		fcindex+=1
	
	return bc1
#
def makeParkfieldETAS(todt=dtm.datetime(2004, 9, 15, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC')), gridsize=.1, contres=3, mc=1.5, kmldir='kml', catdir='kml', fnameroot='parkfield', catlen=5.0*365.0, doplot=False, lons=[-120.95, -119.5], lats=[35.75, 36.5], bigquakes=None, bigmag=4.0, eqtheta=None, eqeps=None, fitfactor=5.0, rtype='ssim'):
	if bigquakes==None:
		# but i don't think bigquakes does anything yet...
		bigquakes=[[mpd.date2num(dtm.datetime(2004, 9, 28, 17, 15, 24,0, tzinfo=pytz.timezone('UTC'))), 35.82, -120.37, 6.0]]
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor, rtype=rtype)
#
class getter(object):
	# a scilly container class to semi-emulate an MPP results 
	def __init__(self, obj):
		self.obj=obj
	def get(self):
		return self.obj
	#
#
def Napa_ApplGeo_sequence(n_cpus=None, gridsize=.1, mc=2.0, lats = [35.3667, 39.7400], lons = [-124.1636, -119.0167], dates=None, figsize=[8,8]):
	# Bob Anderson didn't like the EMC sequence because it was not "California". so we'll do napa. let's pull back a bit to show
	# the potentially related Clear Lake, North-o-Clear Lake, and Parkfieldish events as well (apparenly along the same
	# not fault?)
	#
	# local catalog:
	'''
	2014/01/02 09:32:27.53  38.7867 -122.7408   2.20  3.25   Mw   93  26    1 0.05  NC      72134726
	2014/04/04 04:04:54.97  38.4493 -122.2540   8.22  3.63   Mw  213  67    1 0.20  NC      72193350
	2014/08/05 12:40:01.22  38.2557 -122.3232   8.21  3.03   Mw  174  23    5 0.15  NC      72272361
	2014/08/24 10:20:44.06  38.2155 -122.3117  11.25  6.02   Mw  397  28    4 0.18  NC      72282711
	2014/08/24 10:21:10.84  38.7602 -122.7257   0.93  4.38   Md   26  57    2 0.13  NC      71086369
	2014/08/24 10:21:45.44  38.2350 -122.3198   9.09  3.81   ML  122  50    4 0.17  NC      72282716
	2014/08/24 10:24:44.25  38.2588 -122.3433  10.50  3.47   ML   96 102    5 0.11  NC      72282751
	2014/08/24 12:47:12.56  38.2380 -122.3483   8.47  3.60   Mw  137  30    2 0.17  NC      72283201
	2014/08/26 12:33:16.86  38.1782 -122.3015  12.33  3.90   Mw  282  41    7 0.17  NC      72284586
	2014/08/31 08:56:20.84  38.2352 -122.3293   9.69  3.24   Mw  300  22    3 0.16  NC      72288561
	'''
	napa_eq = {'lat':38.2155, 'lon':-122.3117, 'mag':6.02, 'date_time':dtm.datetime(2014, 8, 24, 10, 20,44, int(.06*10**6.), tzinfo=pytz.timezone('UTC')), 'depth':11.25}
	#
	#lats = [35.3667, 39.7400]		# lat of Chico, lon of Eureka
	#lons = [-124.1636, -119.0167]	# Bakersfield
	#
	#mainshock_datetime = dtm.datetime(2014, 8, 24, 3, 20, 44, tzinfo=pytz.timezone('US/Pacific'))
	mainshock_datetime = napa_eq['date_time']
	#
	my_dt=dtm.timedelta
	#
	if dates == None:
		dates = [mainshock_datetime - my_dt(days=1), mainshock_datetime + my_dt(hours=1)]
		#
		while dates[-1]<mainshock_datetime+my_dt(days=3):
			dates+=[dates[-1]+my_dt(hours=4)]
		while dates[-1]<mainshock_datetime+my_dt(days=5):
			dates+=[dates[-1]+my_dt(hours=6)]
		while dates[-1]<mainshock_datetime+my_dt(days=10):
			dates+=[dates[-1]+my_dt(days=1)]
		while dates[-1]<dtm.datetime(2014,10,14, tzinfo=pytz.timezone('UTC')):
			dates+=[dates[-1]+my_dt(days=3)]
		#
		# any "response" events?
		# all the biggest events happen within minutes of the mainshock, but we can trigger a special run for
		# the 26-aug m=3.9 event:
		dates += [dtm.datetime(2014, 10, 26, 12, 33, 17, tzinfo=pytz.timezone('UTC')) + my_dt(hours=1)]	# and set the ETAS for 1 hour after the event, so we don't
																# get full saturation.
	#
	dates.sort()
	#
	prams_dict = {'doplot':False, 'kmldir':'napa_ag', 'catdir':'napa_ag', 'lats':lats, 'lons':lons, 'gridsize':gridsize, 'contres':10, 'mc':mc, 'bigquakes':None, 'bigmag':5., 'eqtheta':None, 'eqeps':None, 'fitfactor':5.}
	#
	# now, get to work.
	if n_cpus==None: n_cpus=mpp.cpu_count()
	#my_cpu_count = max(1, n_procs-1)	# always leave me at least one cpu ?
	#mypool = mpp.Pool(n_cpus)
	return_etases = []
	#
	for i, this_date in enumerate(dates):
		prams_dict.update({'todt':this_date, 'fignum':i+2, 'fnameroot':'Napa-AG' + str(this_date), 'kmldir':'napa_ag', 'catdir':'napa_ag'})
		#
		# this almost works, but makeETASFCfiles() does some contouring, which invokes pyplot, which breaks...
		#pool_results += [mypool.apply_async(makeETASFCfiles, (), prams_dict)]
		# makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor)
		#
		#SPP it:
		this_etas = makeETASFCfiles(**prams_dict)
		#print "this_etas: ", this_etas
		#
#return_etases += [getter(this_etas)]
		#
		fnum = 2
		#
		plt.figure(fnum, figsize=figsize)
		plt.clf()
		bcm = this_etas.BASScastContourMap(maxNquakes=0, fignum=fnum)
		x,y=this_etas.cm(napa_eq['lon'], napa_eq['lat'])
		plt.figure(fnum)
		#
		plt.plot([x], [y], 'r*', ms=15, alpha=.7, zorder=11)
		plt.title('Napa ETAS: %s\n\n' % str(this_etas.fcdate))
		#
		plt.savefig('%s/napa_etas_%s.png' % (prams_dict['kmldir'], str(this_etas.fcdate)))
		try:
			with open('%s/BASS_napa_%s.pkl' % (prams_dict['kmldir'], str(this_etas.fcdate)), 'w') as f:
				pickle.dump(this_etas, f)
		except Exception as excep:
			print("failed to pickle: ", str(excep))
		#
	#
	#return return_etases
	return this_etas, bcm
#
def EMC_ApplGeo_sequence():
	# a sequence of EMC ETAS (potentially) for the Applied Geology chapter.
	dates = [dtm.datetime(2010,4,1,  tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,4,5, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,4,7, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,4,10, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,4,20, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,4,30, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,5,30, tzinfo=pytz.timezone('UTC')) ]
	#
	#dates=[]		# temporarily..
	#
	# event specific ETAS:
	# following the mainshock (and nearly simultaneous m=5.7 SE of the mainshock)
	dates += [dtm.datetime(2010,4,4, 0, 40,  tzinfo=pytz.timezone('UTC'))] # this is the wrong date, but interestingly enough, it
																			# appears to show a foreshock.
	dates += [dtm.datetime(2010,4,5, 0, 40,  tzinfo=pytz.timezone('UTC'))]
	dates += [dtm.datetime(2010,4,5, 1, 25, tzinfo=pytz.timezone('UTC'))]
	dates += [dtm.datetime(2010,4,8, 18, 44,  tzinfo=pytz.timezone('UTC'))]
	dates += [dtm.datetime(2010,6,15, 6, 26, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC'))]
	#
	dates.sort()
	#
	this_kml_dir = kmldir+'-EMC-AG'
	# set up a pool for multi-processing.
	n_procs = mpp.cpu_count()
	my_cpu_count = max(1, n_procs-1)	# always leave me at least one cpu...
	print("(hopefully) using %d cpus" % my_cpu_count)
	mypool = mpp.Pool(my_cpu_count)
	#
	for i, this_date in enumerate(dates):
		#A = makeElMayorETAS(todt=this_date, doplot=True, fnameroot='EMC-AG')
		prams_dict = {'todt':this_date, 'doplot':False, 'fignum':i+2, 'fnameroot':'EMC-AG' + str(this_date), 'kmldir':this_kml_dir, 'catdir':this_kml_dir}
		#
		# pass positional arguments in a tuple (), then key_work args in a dict.
		# this is the correct syntax, but pyplot requires more careful handling (loop not in main thread or something--
		# pyplot is modal). for now, let's just run them one at a time.
		#mypool.apply_async(makeElMayorETAS, (), prams_dict)		# this falls to pieces during the plots. it is possible to
																 	# pool() these and just pull them back for plotting.
		A = makeElMayorETAS(**prams_dict)
		#
		# sloppy, but we'll get away with it:
		plt.savefig(this_kml_dir + '/EMC-ApplGeo-%s.png' % str(this_date))
		#
		print("queued ETAS: %s" % str(this_date))
	mypool.close()
	mypool.join()
	#
	print("finished.")
	#
	return None
#
def makeElMayorETAS(todt=dtm.datetime(2010,4,1, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC')), gridsize=.1, contres=5, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='elmayor', catlen=5.0*365.0, doplot=False, fignum=1, lons=[-121.0, -114.0], lats=[30.0, 35.25], bigquakes=None, bigmag=5.0, eqtheta=None, eqeps=None, fitfactor=5.0, rtype='ssim', p_quakes=None, p_map=None):
	#
	if bigquakes==None:
		bigquakes=[[mpd.date2num(dtm.datetime(2010, 4, 4, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC'))), 32.286200, -115.295300, 7.2]]
	z= makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor, p_quakes=p_quakes, p_map=p_map)
	#
	z.BASScastContourMap(maxNquakes=10)
	x,y=z.cm(-115.303, 32.128)
	plt.figure(1)
	plt.plot([x], [y], 'r*', ms=15, alpha=.7, zorder=11)
	return z



def makeSocalETAS(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='socal', catlen=5.0*365.0, doplot=False, lons=[-122.5, -114.0], lats=[30.5, 37.25], bigquakes=None, bigmag=5.0, scalemag=4.0, eqtheta=None, eqeps=None, fitfactor=5.0,rtype='ssim'):
	#
	if bigquakes==None:
		bigquakes=[]
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor, rtype=rtype)

def makeNorcalETAS(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='norcal', catlen=5.0*365.0, doplot=False, lons=[-127.5, -117.5], lats=[36.75, 42.0], bigquakes=None, bigmag=5.5, eqtheta=None, eqeps=None, fitfactor=5.0, addquakes=[]):
	#
	if bigquakes==None:
		bigquakes=[]
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor, addquakes=addquakes)
#
def makeCalnevFiles(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='calnev', catlen=5.0*365.0, doplot=False, lons=[-127.5, -114.0], lats=[30.5, 42.5], bigquakes=None, bigmag=5.0):
	#
	if bigquakes==None:
		bigquakes=[]
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag)
#
# , lats=[46.0, 55.0], lons=[-134.0, -124.0]
def makeCascadiaETAS(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, mc=4.25, kmldir=kmldir, catdir=kmldir, fnameroot='cascadia', catlen=5.0*365.0, doplot=False, lats=[47.0, 56.0], lons=[-136.0, -126.0], bigquakes=None, bigmag=6.7, eqtheta=None, eqeps=None, fitfactor=5.0, contour_bottom=.3, contour_top=1.0):
	#
	if bigquakes==None:
		bigquakes=[]
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor, contour_top=contour_top, contour_bottom=contour_bottom)
#
def makeChileETAS(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, mc=4.5, kmldir=kmldir, catdir=kmldir, fnameroot='cascadia', catlen=5.0*365.0, doplot=False, lats=[-30, -10.], lons=[-80.0, -65.0], bigquakes=None, bigmag=6.7, eqtheta=None, eqeps=None, fitfactor=5.0):
	#
	if bigquakes==None:
		bigquakes=[]
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor)
#
def makeTohokuETAS(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=5, mc=4.5, kmldir=kmldir, catdir=kmldir, fnameroot='tohoku', catlen=5.0*365.0, doplot=False, lons=[135., 146.], lats=[30., 41.5], bigquakes=None, bigmag=7.50, eqtheta=None, eqeps=None, fitfactor=5.0, cmfnum=0, fignum=1, contour_intervals=None):
	#
	if bigquakes==None:
		bigquakes=[]
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes[:], bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor, cmfnum=cmfnum, fignum=fignum, contour_intervals=contour_intervals)
#
def makeNZ2013ETAS(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='nz2013', catlen=5.0*365.0, doplot=False, lons=[171.7, 173.7], lats=[42.58, 44.58], bigquakes=None, bigmag=5.5, eqtheta=None, eqeps=None, fitfactor=5.0):
	#
	if bigquakes==None:
		bigquakes=[]
	#
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor)
#
def makeChina2013a(todt=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=3, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='china2013a', catlen=5.0*365.0, doplot=False, lons=[103.5, 105.0], lats=[33.75, 35.25], bigquakes=None, bigmag=5.5, eqtheta=None, eqeps=None, fitfactor=5.0):
	#	 2013-07-21 16:45:56 UTC-07:00
	#		Location 34.499 N 104.243 E
	if bigquakes==None:
		bigquakes=[]
	#
	return makeETASFCfiles(todt=todt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=fnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes, bigmag=bigmag, eqtheta=eqtheta, eqeps=eqeps, fitfactor=fitfactor)
#.
def makeCascadia2014(todt=dtm.datetime(2013, 5, 13, tzinfo=pytz.timezone('US/Pacific-New')), gridsize=.1, contres=5, mc=3.0, kmldir=kmldir, catdir=kmldir, fnameroot='CascadiaScenario2014', catlen=3.0*365.0, doplot=False, lons=[-132.5, -115.0], lats=[35.0, 53.], bigquakes=[], bigmag=7.0, scalemag=5.0, addquakes=[], eqeps=1.414, eqtheta=None, contour_bottom=.5, contour_top=1.0, start_index=0, stop_index=None):
	#
	# 2014 Cascadia scenario for JPL, E-DECIDER.
	if glob.glob(kmldir)==[]: os.system('mkdir %s' % kmldir)	# take a cheap shot at it.
	if glob.glob(catdir)==[]: os.system('mkdir %s' % catdir)
	#
	# eqeps normally
	#dt2=dtm.datetime.now(pytz.timezone('UTC'))
	#dtms  = dtm.datetime(2013, 5, 13, 11, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
	#dtas1 = dtm.datetime(2013, 5, 15, 12, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
	#dtas2 = dtm.datetime(2013, 5, 20, 12, 0,0,0,tzinfo=pytz.timezone('US/Pacific-New'))
	#as2lon, as2lat = -122.08, 38.04
	#
	# i expect this to change, but actaully it's not a bad way to go. this way we have good pre
	# and post event seismicity (though it won't reflect the event's actual aftershock sequence).
	# 45.73 N, 125.12 W)
	#
	# USGS does appear to give us a rupture polygon, but let's just guess a linear approximation.
	# (so far, we don't have any real data, just the shake-map picture).
	faultSouth=[-124.0, 40.]
	faultNorth=[-125.5, 49.25]
	dy=faultNorth[1]-faultSouth[1]
	dx=(faultNorth[0]-faultSouth[0])*math.cos(2.0*math.pi*faultSouth[1]/360.)
	#
	theta_mainshock = math.atan(dy/dx)*360./(2.0*math.pi)	# though this might be negative or should be the other angle (it should be obvious).
	print("fault theta: %f degrees" % theta_mainshock)
	#
	as_faultSouth = [-120.820685, 35.622952]
	as_faultNorth = [-121.102000, 35.706000]
	dy_as=as_faultNorth[1]-as_faultSouth[1]
	dx_as=(as_faultNorth[0]-as_faultSouth[0])*math.cos(2.0*math.pi*as_faultSouth[1]/360.)
	#
	theta_as1 = math.atan(dy_as/dx_as)*360./(2.0*math.pi)	# though this might be negative or should be the other angle (it should be obvious).
	
	#
	event_timezone = pytz.timezone('US/Pacific-New')
	#eq_mainshock = [dtm.datetime(2011, 9, 20, 12, 0,0, tzinfo=event_timezone), 45.1, -125., 9.0]
	eq_mainshock = [dtm.datetime(2014, 5, 12, 9, 41, 0,0, tzinfo=event_timezone), 45.73, -125.12, 9.03, eqeps, theta_mainshock]
	#
	# aftershock: 40.802914 and -124.161832, depth 5 km, Mw = 7.2  # (time does not appear to be provided. make it 12 hours?)
	#
	if addquakes==None:
		addquakes=[]
	if addquakes==[]:
		addquakes += [eq_mainshock[:]]
		addquakes[-1][0]=mpd.date2num(addquakes[-1][0])
		#
		# date-numbers are in days (keep it simple with 12 hour increments for the simulation):
		addquakes += [[mpd.date2num(dtm.datetime(2014,5,14,10, tzinfo=event_timezone)), 40.802914, -124.161832, 7.2, eqeps, theta_as1]]
		#addquakes += [[addquakes[-1][0]+3.5, 40.802914, -124.161832, 7.2, eqeps, theta_as1]]   # "official(?)" aftershock?
																			# falut will be: [35.706000 -121.102000], [35.622952, -120.820685] ~ 19.91 degrees S of E
		# my proposed aftershocks:
		#addquakes += [[addquakes[-1][0]+.5, 46.5, -125.5, 7.2]]
		#addquakes += [[addquakes[-1][0]+.5, 40., -125., 7.5]]
		#addquakes += [[addquakes[-1][0]+.5, 46.0, -124.5, 7.0]]
		#addquakes += [[addquakes[-1][0]+.5, 42., -124.0, 7.7]]
		#
		#
	print("addquakes:", addquakes)
	if bigquakes==None:
		bigquakes=[]
	#
	#faultTheta = (90.+60.)
	#
	delta_t = dtm.timedelta(hours=12.0)
	fcdate0 = eq_mainshock[0]-delta_t
	#fcdate0=dtm.datetime(2013, 5, 13, 00, 0, 0, 0, tzinfo=event_timezone)	# nominally, a "standard" time.
	#
	# date-times at which ETAS will be run (assuming we have fast, parallel machines so that the
	# run-time reflects the current -- rapidly changing catalog).
	#
	# so far, this scenario is pretty simple. let's do 1 before, 1 immediately after,
	# maybe a few after.
	#while (fcdates[-1]<=dtm.datetime(2013, 5, 22, 00, 0, 0, 0, tzinfo=event_timezone)):
	#	fcdates+=[fcdates[-1]+dtm.timedelta(hours=12.0)]
		#print "adding fcdate: %s" % str(fcdates[-1])
	fcdates = [fcdate0] 
	fcdates += [eq_mainshock[0]+dtm.timedelta(seconds=15*60)]
	for i in range(2):
		fcdates += [fcdates[-1] + dtm.timedelta(hours=1)]
	for i in range(10):
		fcdates += [fcdates[-1] + dtm.timedelta(hours=8)]
	#
	# now, add "response" forecasts, run after a large aftershock.
	for i in range(1, len(addquakes)):
		# about 14.4 minutes after each aftershock.
		fcdates += [mpd.num2date(addquakes[i][0]+.01)]
	#
	fcdates.sort()
	for rw in fcdates:
		print("adding fcdate: %s" % str(rw))
	fcindex=start_index
	#
	plt.close(0)
	plt.close(1)
	if stop_index==None or stop_index==0: stop_index=len(fcdates)
	for fcdt in fcdates[start_index:stop_index]:
		plt.figure(num=1, figsize=(8,10))
		plt.clf()
		#
		print("forecastdate: %s" % fcdt)
		#
		#activeEps=eqeps
		#activeTheta=eqtheta
		activeEps=None
		activeTheta=None
		#
		#if fcdt>eq_mainshock[0]: 
		#	activeTheta=faultTheta
		#	#activeEps=2.0	# anyway, for cascadia, we expect a broader rupture than CA... but check out the AK Shield exercise to specify orientation for specific earthquakes...
		#
		#fcdt_string = str(fcdt).replace(' ', '_').replace(':', '_')
		fcdt_string = '%d-%d-%d_%d_%d_%d' % (fcdt.year, fcdt.month, fcdt.day, fcdt.hour, fcdt.minute, fcdt.second)
		#thisfnameroot = fnameroot + '-%d' % fcindex
		thisfnameroot = '%s-%s' % (fcdt_string, fnameroot)
		bc1=None
		bc1 = makeETASFCfiles(todt=fcdt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=thisfnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes[:], bigmag=bigmag, addquakes=addquakes, eqeps=activeEps, eqtheta=activeTheta, cmfnum=0, fignum=1, colorbar_fontcolor='w', contour_bottom=contour_bottom, contour_top=contour_top)
		#
		# for batched jobs on remote machines, it may be necessary to forego the plots (sometimes this generates a "can't render graphics"
		# type error).
		#
		mfb=bc1.BASScastContourMap(fignum=1)
		# and we'll want to draw the fault here as well...
		#Fault endpoints: 36.82 -121.51, 40.25 -124.41
		x1,y1=bc1.cm(faultSouth[0], faultSouth[1])
		x2, y2=bc1.cm(faultNorth[0], faultNorth[1])
		#xms, yms=bc1.cm(addquakes[0][2], addquakes[0][1])
		xms, yms=bc1.cm(eq_mainshock[2], eq_mainshock[1])
		#xas1, yas1=bc1.cm(addquakes[1][2], addquakes[1][1])
		#xas2, yas2=bc1.cm(as2lon, as2lat)
		plt.plot([x1, x2], [y1, y2], 'm-s', lw=3, alpha=.7)
		#
		#if fcdt>=dtms: 
		plt.plot([xms], [yms], 'b*', alpha=.8, ms=2.5*addquakes[0][3])
		#if fcdt>=dtas1: plt.plot([xas1], [yas1], 'g*', alpha=.8, ms=2.5*addquakes[1][3])
		#if fcdt>=dtas2: plt.plot([xas1], [yas1], 'c*', alpha=.8, ms=2.5*addquakes[1][3])
		plt.title('fcdate: %s, m>%.2f\n\n\n' % (fcdt, mc))
		plt.savefig('%s/%sconts-%d.png' % (kmldir,fnameroot, fcindex))
		#
		fcindex+=1
	
	return bc1
#
def makeAKshield2014(todt=dtm.datetime(2013, 5, 13, tzinfo=pytz.timezone('US/Pacific-New')), gridsize=.1, contres=5, mc=3.0, kmldir=kmldir+'AK', catdir=kmldir+'AK', fnameroot='AKScenario2014', catlen=3.0*365.0, doplot=False, lons=[-168., -138.0], lats=[50.0, 64.], bigquakes=[], bigmag=7.0, scalemag=5.0, addquakes=[], eqeps=1.414, eqtheta=30., starting_i=0):
	#
	# 2014 Cascadia scenario for JPL, E-DECIDER.
	if glob.glob(kmldir)==[]: os.system('mkdir %s' % kmldir)	# take a cheap shot at it.
	if glob.glob(catdir)==[]: os.system('mkdir %s' % catdir)
	#
	# eqeps normally
	#dt2=dtm.datetime.now(pytz.timezone('UTC'))
	#dtms  = dtm.datetime(2013, 5, 13, 11, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
	#dtas1 = dtm.datetime(2013, 5, 15, 12, 0, 0, 0, tzinfo=pytz.timezone('US/Pacific-New'))
	#dtas2 = dtm.datetime(2013, 5, 20, 12, 0,0,0,tzinfo=pytz.timezone('US/Pacific-New'))
	#as2lon, as2lat = -122.08, 38.04
	#
	# i expect this to change, but actaully it's not a bad way to go. this way we have good pre
	# and post event seismicity (though it won't reflect the event's actual aftershock sequence).
	event_timezone = pytz.timezone('US/Pacific-New')
	eq_mainshock = [dtm.datetime(2014, 0o3, 27, 10, 54,0, tzinfo=event_timezone), 55.2, -156.7, 8.2, 1.0, 0.]
	#
	# mainshock magnitude updates.
	mainshocks = []
	mainshocks += [[dtm.datetime(2014, 0o3, 26, 10, 31, 0, tzinfo=event_timezone), 55.2, -156.7, 4.0, eqeps, eqtheta]]
	mainshocks += [eq_mainshock]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 8.6, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	mainshocks += [[dtm.datetime(2014, 0o3, 27, 10, 54, 0, tzinfo=event_timezone), 55.2, -156.7, 9.1, eqeps, eqtheta]]
	#
	update_times = []
	update_times += [dtm.datetime(2014, 0o3, 26, 11, 31, 0, tzinfo=event_timezone)]
	update_times += [eq_mainshock[0]]
	update_times += [dtm.datetime(2014, 0o3, 27, 11, 31, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 12,  3, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 13,  5, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 14,  1, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 15,  1, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 16,  2, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 17,  2, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 18,  2, 0, tzinfo=event_timezone)]
	update_times += [dtm.datetime(2014, 0o3, 27, 19,  2, 0, tzinfo=event_timezone)]
	#
	# aftershocks:
	aftershocks=[]
	aftershocks += [[dtm.datetime(2014, 3, 27, 10, 40, 0, tzinfo=event_timezone), 61.5839, -150.1473, 6.1, eqeps, 258.]]
	aftershocks += [[dtm.datetime(2014, 3, 27, 13, 55, 0, tzinfo=event_timezone), 61.5612, -150.1170, 5.1, None, None]]
	aftershocks += [[dtm.datetime(2014, 3, 28,  5, 30, 0, tzinfo=event_timezone),   52.07, -167.7291, 7.8, eqeps, 248.]]
	aftershocks += [[dtm.datetime(2014, 3, 28,  8,  0, 0, tzinfo=event_timezone),  57.824,  -152.650, 6.5, eqeps, 240.]]
	aftershocks += [[dtm.datetime(2014, 3, 28,  9,  0, 0, tzinfo=event_timezone),  60.7203, -145.5957, 6.4, eqeps, 258.]]
	aftershocks += [[dtm.datetime(2014, 3, 28,  9, 18, 0, tzinfo=event_timezone), 53.6187, -133.0468, 6.8, eqeps, 333.]]
	aftershocks += [[dtm.datetime(2014, 3, 28,  9, 45, 0, tzinfo=event_timezone), 63.6128, -147.9054, 7.8, eqeps, 298.]]
	aftershocks += [[dtm.datetime(2014, 3, 28, 10,  0, 0, tzinfo=event_timezone), 61.7404,  -151.006, 5.6, None, None]]
	aftershocks += [[dtm.datetime(2014, 3, 28, 13, 00, 0, tzinfo=event_timezone), 60.3494, -145.8620, 6.4, None, None]]
	aftershocks += [[dtm.datetime(2014, 3, 28, 13, 30, 0, tzinfo=event_timezone), 60.5363, -145.491, 5.2, None, None]]
	aftershocks += [[dtm.datetime(2014, 3, 28, 13, 30, 0, tzinfo=event_timezone), 61.7238, -150.1509, 6.8, eqeps, 290.]]
	aftershocks += [[dtm.datetime(2014, 3, 29,  9, 15, 0, tzinfo=event_timezone),  61.558, -149.6591, 6.8, eqeps, 250.]]	
	aftershocks += [[dtm.datetime(2014, 3, 29,  9, 56, 0, tzinfo=event_timezone), 60.9892, -147.0379, 7.5, eqeps, 230.]]
	#
	# let's add a bunch of new update times for the new aftershocks (this is a holy mess, but there's a deadline):
	# assume we update 15 mins after each large aftershock (let's just do this independent of the scheduled updates).
	for aftershock in aftershocks:
		#print "aftershock: %s" % str(aftershock)
		update_times += [aftershock[0] + dtm.timedelta(seconds=15*60)]
		#
	update_times.sort()
	#
	'''
	if addquakes==None:
		addquakes=[]
	if addquakes==[]:
		addquakes += [eq_mainshock[:]]
		addquakes[-1][0]=mpd.date2num(addquakes[-1][0])
		#
		# date-numbers are in days (keep it simple with 12 hour increments for the simulation):
		addquakes += [[addquakes[-1][0]+.5, 46.5, -125.5, 7.2]]
		addquakes += [[addquakes[-1][0]+.5, 40., -125., 7.5]]
		addquakes += [[addquakes[-1][0]+.5, 46.0, -124.5, 7.0]]
		addquakes += [[addquakes[-1][0]+.5, 42., -124.0, 7.7]]
		#
		#
	3print "addquakes:", addquakes
	if bigquakes==None:
		bigquakes=[]
	#
	'''
	# Guess 30 deg N of E...
	# (and i don't remember which way the angle goes...)
	faultTheta = -30.0
	'''
	faultSouth=[-124.0, 40.]
	faultNorth=[-125.5, 49.25]
	dy=faultNorth[1]-faultSouth[1]
	dx=(faultNorth[0]-faultSouth[0])*math.cos(2.0*math.pi*faultSouth[1]/360.)
	#
	faultTheta =math.atan(dy/dx)*360./(2.0*math.pi)	# though this might be negative or should be the other angle (it should be obvious).
	print "fault theta: %f degrees" % faultTheta
	#faultTheta = (90.+60.)
	'''
	#
	delta_t = dtm.timedelta(hours=12.0)
	fcdate0 = eq_mainshock[0]-delta_t
	#
	plt.close(0)
	plt.close(1)
	#fcindex=0
	fcindex = starting_i
	# this is sloppy, but we'll get away with it (clean up before cloning this)
	#for i in xrange(starting_i, len(mainshocks)):
	for i in range(starting_i, len(update_times)):
	#for rw in mainshocks:
		mainshock_index = min(i, len(mainshocks)-1)
		this_mainshock=mainshocks[mainshock_index]
		fcdt = update_times[i]
		
		#fcdt = rw[0]
		plt.figure(num=1, figsize=(8,10))
		plt.clf()
		#
		print("forecastdate: %s" % fcdt)
		#addquakes=[rw[0:4]]
		addquakes = [this_mainshock]
		j=0
		while j<len(aftershocks) and aftershocks[j][0]<= fcdt:
			addquakes+=[aftershocks[j]]
			j+=1
		#
		# how hard will it be to specify epsilon, theta for each earthquake?
		#activeEps=this_mainshock[4]
		#activeTheta=this_mainshock[5]
		activeEps=None
		activeTheta=None
		#
		#if fcdt>eq_mainshock[0]: 
		#	activeTheta=faultTheta
		#	#activeEps=2.0	# anyway, for cascadia, we expect a broader rupture than CA.
		#
		thisfnameroot = fnameroot + '-%d' % fcindex
		bc1=None
		bc1 = makeETASFCfiles(todt=fcdt, gridsize=gridsize, contres=contres, mc=mc, kmldir=kmldir, catdir=catdir, fnameroot=thisfnameroot, catlen=catlen, doplot=doplot, lons=lons, lats=lats, bigquakes=bigquakes[:], bigmag=scalemag, addquakes=addquakes, eqeps=activeEps, eqtheta=activeTheta, cmfnum=0, fignum=1, colorbar_fontcolor='w')
		#
		mfb=bc1.BASScastContourMap(fignum=1)
		# and we'll want to draw the fault here as well...
		'''
		#Fault endpoints: 36.82 -121.51, 40.25 -124.41
		x1,y1=bc1.cm(faultSouth[0], faultSouth[1])
		x2, y2=bc1.cm(faultNorth[0], faultNorth[1])
		'''
		#xms, yms=bc1.cm(addquakes[0][2], addquakes[0][1])
		xms, yms=bc1.cm(this_mainshock[2], this_mainshock[1])
		#xas1, yas1=bc1.cm(addquakes[1][2], addquakes[1][1])
		#xas2, yas2=bc1.cm(as2lon, as2lat)
		#plt.plot([x1, x2], [y1, y2], 'm-s', lw=3, alpha=.7)
		#
		#if fcdt>=dtms: 
		plt.plot([xms], [yms], 'b*', alpha=.8, ms=2.5*addquakes[0][3])
		#if fcdt>=dtas1: plt.plot([xas1], [yas1], 'g*', alpha=.8, ms=2.5*addquakes[1][3])
		#if fcdt>=dtas2: plt.plot([xas1], [yas1], 'c*', alpha=.8, ms=2.5*addquakes[1][3])
		plt.title('fcdate: %s, m>%.2f\n\n\n' % (fcdt, mc))
		plt.savefig('%s/%sconts-%d.png' % (kmldir,fnameroot, fcindex))
		fcindex+=1
	
	return bc1
#
def tohokuETASmovie(gridsize=.2, maxNquakes=10, moviedir='tohokumovie'):
	#deltat=dtm.timedelta(hours=3)
	# which r-distribution are we using?
	# ssim
	deltat=dtm.timedelta(days=1)
	#date0=dtm.datetime(2009,4,8, 0, tzinfo=pytz.timezone('UTC'))
	date0=dtm.datetime(2012,1,7, 0, tzinfo=pytz.timezone('UTC'))
	datef=dtm.datetime(2012,6,14, 0, tzinfo=pytz.timezone('UTC'))
	#datef=dtm.datetime(2009,4,19, 0, tzinfo=pytz.timezone('UTC'))
	#
	# define contours:
	contints=[-13.0]
	while contints[-1]<-7.5:
		contints+=[contints[-1]+.5]
	#
	i=1003
	nprocs=2
	plt.close(1)
	while date0<datef:
		indstr='0000000%d' % i
		indstr=indstr[-4:]
		fnameroot='tohoku%s' % indstr
		#
		#plt.close(1)
		plt.figure(num=1, figsize=(13,13))
		plt.clf()
		a=None
		a=makeTohokuETAS(todt=date0, gridsize=gridsize, kmldir=moviedir, catdir=kmldir, fnameroot=fnameroot, cmfnum=0, fignum=1, contour_intervals=contints)
		a.BASScastContourMap(maxNquakes=maxNquakes)
		plt.figure(1)
		plt.title('%s\n\n\n' % str(date0))
		#
		
		plt.savefig('%s/tohokuetas%s.png' % (moviedir,indstr))
		date0=date0+nprocs*deltat
		i+=1*nprocs
	# avconv -i tohokumovie/tohokuetas%04d.png -r 25 tohokuetasmovie.mp4
	
def makeColormaps(colormap='spectral', nColors=1024):
	# a colormap for javascript...
	cm = plt.get_cmap(colormap)
	cNorm = colors.Normalize(0, nColors-1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	#
	colorList=["#000000"]
	i=1
	#while colorList[-1]!='#ffffff':
	for i in range(nColors):
		colorVal = scalarMap.to_rgba(float(i))
		# now, convert to hex:
		colorValHex = '#'
		for z in colorVal:
			zbit=hex(255*z)[2:4]
			if zbit[-1]=='L': zbit='0' + zbit[0]
			colorValHex+=zbit
		colorList+=[colorValHex[0:-2]]
		i+=1
	#
	return colorList

	
def etasFigs6(todt=dtm.datetime(2013,4,17, tzinfo=pytz.timezone('UTC'))):
	# socal ETAS
	if todt=='now': todt=dtm.datetime.now(pytz.timezone('UTC'))
	
	b=makeSocalETAS(todt=todt, fnameroot='socalPAGEOPH', rtype='ssim')
	b.contres=7		# or i think we can use the contour_intervals= keyword and use an integer (or a set of specific contour values).
	b.BASScastContourMap(maxNquakes=10)
	dtstr='%d%s%s' % (todt.year, ('00'+str(todt.month))[-2:], ('00'+str(todt.day))[-2:])
	plt.savefig('figs/socalETAS-ssim-%s.png' % dtstr)
	#
	b=makeSocalETAS(todt=todt, fnameroot='socalPAGEOPH', rtype='omorisat')
	b.contres=7		# or i think we can use the contour_intervals= keyword and use an integer (or a set of specific contour values).
	b.BASScastContourMap(maxNquakes=10)
	dtstr='%d%s%s' % (todt.year, ('00'+str(todt.month))[-2:], ('00'+str(todt.day))[-2:])
	plt.savefig('figs/socalETAS-threshSatPL-%s.png' % dtstr)
	
	
def etasFigs7():
	# El-Mayor ETAS
	c=makeElMayorETAS(contres=7, rtype='ssim')
	c.BASScastContourMap(maxNquakes=10)
	x,y=c.cm(-115.303, 32.128)
	plt.figure(1)
	plt.plot([x], [y], 'r*', ms=15, alpha=.7, zorder=11)
	plt.plot([x], [y], 'k*', ms=18, alpha=.9, zorder=10)
	#
	plt.savefig('figs/elmayorETAS-ssim-20100401.png')
	#
	c=makeElMayorETAS(contres=7, rtype='omorisat')
	c.BASScastContourMap(maxNquakes=10)
	x,y=c.cm(-115.303, 32.128)
	plt.figure(1)
	plt.plot([x], [y], 'r*', ms=15, alpha=.7, zorder=11)
	plt.plot([x], [y], 'k*', ms=18, alpha=.9, zorder=10)
	#
	plt.savefig('figs/elmayorETAS-threshSatPL-20100401.png')
	
if __name__=='__main__':
	pass
	# but eventually, do the mpl.use('Agg') or whatever puts mpl into the background...
else:
	plt.ion()

