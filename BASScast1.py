import datetime as dtm
import pytz
import operator
import math
import random
import numpy
import scipy
import os
from PIL import Image as ipp
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mpd
import matplotlib.mpl as mpl
#
import shapely.geometry as sgp
import polytest as ptp

#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import ANSStools as atp
#
days2secs = 60.*60.*24.
year2secs = 60.*60.*24.*365.
deg2km=111.12
deg2rad = 2.0*math.pi/360.
alpha_0 = 2.28

class BASScast(object):
	# BASS/ETAS type forecast container.
	sites=[]
	quakes=[]
	conts=[]			# contour (local) objects (from sites x,y,z data).
	latrange=[]
	lonrange=[]		# lat/lon range (rectangular) of catalog
	gridsize=.1		# degrees
	midlat=None		# average latitude of the catalog, though i think we'll employ a more accurate calc. that uses both lat coords (dx = cos(l2)*x2 - cos(l1)*x1)
	mc=2.0			# minimum magnitude (at least for BASS/ETAS calculations)
	contres=2
	contfact=5
	fcdate=None
	fcdatef=None	# forecast date, float(forecast date)
	deg2km=111.12	# degrees (latitute) to km conversion
	zresolution=9	# decimal resolution of Z_ij array from which contours are calculated.
	#
	reffMag=5.0	# equivalent number of this magnitude/year earthquakes.
							# so, the contours can be read as " n m>5/km^2. nominally
							# we can integrate the contours as well if we like...
							# but from this we should be able to talk about critical rates et al.
							# n -> n0*10^(mc-rateMag.)
	mapres='i'		# basemap map resolution.
	#
	def __init__(self, incat=[], fcdate=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=2, mc=2.0):
		return self.initialize(incat=incat, fcdate=fcdate, gridsize=gridsize, contres=contres, mc=mc)
		#
	#
	def initialize(self, incat=[], fcdate=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=2, mc=2.0):
		#
		# incat: [ [dtm, lat, lon, mag, (depth?)], [], ... ] (consistent with ypp.eqcatalog() format).
		# note: catalog dates/times must be in seconds (or convert all the scaling bits to days-- yuck).
		self.mc=mc
		self.contres=contres
		self.gridsize=gridsize
		self.fcdate=fcdate
		self.fcdatef = mpd.date2num(fcdate)*days2secs
		#
		# read a catalog
		if type(incat)==type('astring'):
			# it's a filename. assume it's formatted for a catalog object...
			# ... and we're not sure what to do with this just yet. for now, require a list
			# [ [dtm, lat, lon, mag, (depth?)], [], ... ] (consistent with ypp.eqcatalog() format.
			#
			# load catalog into [quakes] and determine X, Y range from catalog.
			print "we need a list-catalog. files are not supported at this time."
			#
			return None
		#
		self.catalog=incat	# dt, lat, lon, mag
		#
		# permit an empty catalog:
		if len(self.catalog)==0:
			# nothing to do, but we can still use the member functions.
			return None
			
		# initialize the lat/lon range arrays:
		self.latrange = [incat[0][1], incat[0][1]]
		self.lonrange = [incat[0][2], incat[0][2]]
		#
		# add earthquakes and other catalog related bits.
		for rw in incat:
			self.quakes+=[earthquake(mag=rw[3], loc=[rw[2], rw[1]], evtime=rw[0], mc=mc)]
			#
			# rupture length, to modify catalog dimensions (we set the grid to the lat/lon from the catalog + rupture-len)
			rl=self.quakes[-1].getruptlen()
			dlon=rl/(self.deg2km*math.cos(rw[1]*2.0*math.pi/360.0))		# noting the deg -> rad and the km -> deg. conversions.
			dlat=rl/self.deg2km
			# adjust the forecast extents according to catalog events and their rupture lens.
			if self.latrange[0]>(rw[1]-dlat): self.latrange[0] = (rw[1]-dlat)		# expand min-lat?
			if self.latrange[1]<(rw[1]+dlat): self.latrange[1] = (rw[1]+dlat)		# expand max-lat?
			if self.lonrange[0]>(rw[2]-dlon): self.lonrange[0] = (rw[2]-dlon)
			if self.lonrange[1]<(rw[2]+dlon): self.lonrange[1] = (rw[2]+dlon)
			#
		self.midlat = .5*(sum(self.latrange))
		#
		# for now, let's use gridsize x gridsize (in degrees) size elements (rather than project into a proper rectangular coord-system,
		# though this is potentially not terribly difficult for this application since the grid is not strictly enforced (aka, we
		# have a set of center-points for which we calc. z-values based on distance between [(x1,y1), (x2, y2)]. the potential difficulty arises
		# from the fact that the array from which contours are calculated must be square... so maybe it is difficult.
		#
		# set up forecast sites:
		x0 = self.lonrange[0]-self.lonrange[0]%gridsize
		xmax = gridsize + self.lonrange[1]-self.lonrange[1]%gridsize
		self.nLon = int(math.ceil((xmax-x0)/gridsize))
		#
		y0 = self.latrange[0]-self.latrange[0]%gridsize
		ymax = gridsize + self.latrange[1]-self.latrange[1]%gridsize
		self.nLat = int(math.ceil((ymax-y0)/gridsize))
		#
		#print self.lonrange, self.latrange
		y=y0
		#
		self.sites=[]
		for i in xrange(self.nLat*self.nLon):
			self.sites+=[forecastsite(loc=[x0+gridsize*(i%self.nLon), y0+gridsize*(i/int(self.nLon))], dxdy=[gridsize, gridsize], evtime=self.fcdatef, mc=mc)]
		#		
		#print "lengths: %d, %d, %d" % (self.nLon, self.nLat, len(self.sites))
		#
		# and now, this object should br ready to calculate forecasts...
		#bc=self.calcBASScast()
		#bc=self.plotContours()
		bcast=self.calcBASScast()
		#self.conts = self.getContourSet(X_i=self.X_i, Y_i=self.Y_i, Z_ij=self.Z2d.round(self.zresolution), contres=self.contfact*self.contres)
		#self.conts = self.getContourSet(X_i=bcast[0], Y_i=bcast[1], Z_ij=bcast[2].round(4), contres=self.contfact*self.contres)
		self.conts = self.getContourSet(X_i=bcast[0], Y_i=bcast[1], Z_ij=bcast[2], contres=self.contfact*self.contres)
		
		return None
	#
	def resetzvals(self):
		# set all z to 0 or None.
		for i in xrange(len(self.sites)):
			self.sites[i].z=None
	#
	def calcBASScast(self, grid=None, quakes=None, nx=None, ny=None, gridsize=None):
		if grid==None: grid = self.sites
		if quakes==None: quakes = self.quakes
		if nx==None: nx = self.nLon
		if ny==None: ny = self.nLat
		if gridsize==None: gridsize=self.gridsize
		#
		X=[]
		Y=[]
		Z=[]	
		for site in grid:
			site.setz(equakes=quakes)
			# note, z values (eq counts) have been summed at this point
			#
			X+=[site.loc[0]]
			Y+=[site.loc[1]]
			if site.z!=0: Z+=[math.log10(site.z)]
			if site.z==0: Z+=[None]
			#
		#	
		#Z=map(operator.itemgetter(2), cdata)
		Z2d=scipy.array(Z)
		#
		#Z2d.shape=(len(Y), len(X))
		Z2d.shape=(ny, nx)		# aka, (ny rows, of nx elements) such that ny*nx=N
		#X1, Y1=numpy.meshgrid(X, Y)
		#
		X_i=numpy.array(map(float, range(nx)))
		X_i*=gridsize
		X_i+=grid[0].loc[0]
		#
		Y_i=numpy.array(map(float, range(ny)))
		Y_i*=gridsize
		Y_i+=grid[0].loc[1]
		#
		self.X_i = X_i
		self.Y_i = Y_i
		self.Z2d=Z2d
		#
		return [X_i, Y_i, Z2d]
	
	#self.fillHexString=fillHexString
	
	def gridplot(self, ncontours=25):
		plt.figure(0)
		plt.clf()
		minval=self.Z2d.min()
		maxval=self.Z2d.max()
		zrange=maxval-minval
		#
		hexColorMin='0x000000'
		hexColorMax='0xFFFFFF'
		colorRange=int('0xFFFFFF', 0)
		#
		#print "gridlens: %d, %d" % (len(self.Z2d), len(self.Z2d[0]))
		#		
		#dz=(maxval-minval)/float(ncontours)
		
		for y in xrange(len(self.Z2d)):
			for x in xrange(len(self.Z2d[0])):
				colorint=int(colorRange*(self.Z2d[y,x]-minval)/zrange)
				#colorhex=self.fillHexString(hex(colorint), 6)
				colorhex=fillHexString(hex(colorint), 6)
				colorhexstr='#' + colorhex.split('x')[1]
				plt.plot([x], [y], '.', color=colorhexstr, ms=10)
				
	
	'''
	def plotContours(self, grid=None, quakes=None, nx=None, ny=None, gridsize=None):
		XYZ=self.calcBASScast(grid=grid, quakes=quakes, nx=nx, ny=ny, gridsize=gridsize)
		return self.BASScastContours(XYZ[0], XYZ[1], XYZ[2])
	
	def mapContours(self, grid=None, quakes=None, nx=None, ny=None, gridsize=None):
		XYZ=self.calcBASScast(grid=grid, quakes=quakes, nx=nx, ny=ny, gridsize=gridsize)
		return self.BASScastContourMap(XYZ[0], XYZ[1], XYZ[2])
		
	def BASScastContours(self, X_i, Y_i, Z2d, fignum=0):
		#
		if fignum!=None: plt.figure(fignum)
		#plt.clf()
		cnts = self.getContourSet(X_i, Y_i, Z2d, self.contfact*self.contres)
		#
		self.conts = cnts
		#
		return None
	'''
		
	def BASScastContourMap(self, X_i=None, Y_i=None, Z2d=None, fignum=1, maxNquakes=250.0):
		#
		if X_i==None: X_i = self.X_i
		if Y_i==None: Y_i = self.Y_i
		if Z2d==None: Z2d = self.Z2d
		#
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		#
		cntr = [self.lonrange[0] + .5*(self.lonrange[1]-self.lonrange[0]), self.latrange[0] + .5*(self.latrange[1]-self.latrange[0])]
		cm = Basemap(llcrnrlon=self.lonrange[0], llcrnrlat=self.latrange[0], urcrnrlon=self.lonrange[1], urcrnrlat=self.latrange[1], resolution=self.mapres, projection='tmerc', lon_0=cntr[0], lat_0=cntr[1])
		self.cm=cm
		cm.drawcoastlines(color='gray', zorder=1)
		cm.drawcountries(color='gray', zorder=1)
		cm.drawstates(color='gray', zorder=1)
		cm.drawrivers(color='gray', zorder=1)
		cm.fillcontinents(color='beige', zorder=0)
		
		cm.drawmeridians(range(int(self.lonrange[0]), int(self.lonrange[1])), color='k', labels=[1,1,1,1])
		cm.drawparallels(range(int(self.latrange[0]), int(self.latrange[1])), color='k', labels=[1, 1, 1, 1])
		#
		# get X,Y for contours:
		#
		X,Y=cm(*numpy.meshgrid(X_i, Y_i))
		#X,Y=cm(X_i, Y_i)
		#		
		cnts = self.getContourSet(X, Y, Z2d, self.contfact*self.contres)
		plt.colorbar()
		plt.spectral()
		#
		#maxNquakes=250.0	# max num. equakes to plot (approximately)
		Nquakes = float(len(self.quakes))
		#
		mthresh=math.log10(Nquakes/maxNquakes) + self.mc	# we can take  out the "if" bc for N<Nmax, mthresh<mc
		# now, plot the earthquakes:
		for q in self.quakes:
			if q.mag<mthresh: continue
			qx, qy, qm = q.loc[0], q.loc[1], q.mag
			qxp, qyp = cm(qx, qy)
			cm.plot([qxp], [qyp], 'ro', ms=qm, alpha=.7)
		#
		#self.conts = cnts
		#
		return None
	
	def getContourSet(self, X_i=None, Y_i=None, Z_ij=None, contres=None):
		#
		if X_i==None: X_i=self.X_i
		if Y_i==None: Y_i=self.Y_i
		if Z_ij==None: Z_ij=self.Z2d
		if contres==None: contres=self.contfact*self.contres
		#
		# X_i, Y_i are x,y coordinates (aka, like range(100), range(200) ), Z_ij is a 2D array that is actually contoured.
		# contres is the resolution of the contoursl; numConts -> 5*contres
		#
		resolution = contres
		LevelsNumber = self.contfact * resolution
		warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#
		# cnts=plt.contourf(gridsize*numpy.array(range(nelements)), gridsize*numpy.array(range(nelements)), Z2d,10)
		# retrieve the collections() object from the contourf() function which returns a matplotlib.contour.ContourSet
		#
		#cs = plt.contourf(X_i, Y_i, Z_ij, LevelsNumber, cm=plt.spectral()).collections
		cs = plt.contourf(X_i, Y_i, Z_ij, LevelsNumber, cm=plt.spectral(), alpha=.3)
	
		return cs

	def arrayFromSet(self, cs=None):
		# get an array (list) from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
		# format:
		# [ [z0, [[x,y], [x,y],...], [z1, [paths1]], 
		#
		if cs==None: cs=self.conts
		#
		levels=cs.levels	# lower bound of contours
		layers=cs.layers	# upper bound of contours (at least for contours <0)
		dz= layers[0] - levels[0]
		collects=cs.collections
		carray = []
		for i in xrange(0, len(collects)):
			bgra_array = 255*collects[i].get_facecolor()[0]
			strclr = '7d%02x%02x%02x' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
			#
			for trace in collects[i].get_paths():
				#carray+=[[levels[i], layers[i], '%s' % strclr, []]]
				# objContour(low=None, high=None, alpha=None, R=None, G=None, B=None, verts=None, RGB=None)
				carray += [objContour(low=cs.levels[i], high=cs.layers[i], RGB=strclr, verts=[])]
			
				for lng, lat in trace.vertices:
					#carray[-1][-1]+=[[lng, lat]]
					carray[-1].verts+=[[lng, lat]]
					#
				#
			#
		#
		#
		return carray
	#
	def contsToStr(self, cset=None):
		if cset==None: cset=self.conts
		cs=cset.collections
		outstr='# contour coordinates\n'
		#
		for i in xrange(len(cs)):
			outstr+='#contour\t%d\n' % i
			itrace=0
			for trace in cs[i].get_paths():
				#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				outstr+='#trace\t%d\n' % itrace
				for lng,lat in trace.vertices:
					#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
					outstr+='%f,%f,0\n' % (lng, lat)
					#tmpLL+=[[lng, lat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])
		#
		return outstr

	#self.plotcolor=plotcolor
	
	def Z2str(self):
		outstr=''
		for rw in self.Z2d:
			for elem in rw:
				outstr+='%f\t' % elem
			outstr=outstr[:-1] + '\n'
		#
		return outstr
	#
	def contsToPlotLists(self, cset=None, getlevels=False):
		# this corrects a mpl/kml mismatch. mpl likes to plot inner polygons by drawing a line from the outside
		# ring to the inner poly, draw the inenr poly, then back to the outer ring over the first line. mpl
		# interprets the line as infinitely thin and basically ignores it. KML will not ignore that line, and so
		# you get a big fat mess. This script separates the outer/inner groups into outer and inner poly rings.
		#
		# BUT, (2012 09 05) this needs to be revisited. inner-inner polys are not plotting correctly. aka,
		# the z-profile like: /-\_/-\_/-\ (like an s0 within an s1) need to be considered. in such a case,
		# we get: <outer><inner><inner></inner></inner></outer>, but there should be 1 outer/inner poly then a
		# separate poly entirely. the second "inner" needs to be pulled out as a separate poly. (see similar note below)
		#
		# so, polys that are inside an odd number of polys are "inner" polys; polys inside an even number are "outside" 
		# any poly (as per knot theory), or otherwise constitute an "outer" poly. poly1 is an "inner" poly of poly2 if:
		#   1: it is inside poly2
		#   2: its "inner index n2 = n1 - 1
		if cset==None: cset=self.conts
		cs=cset.collections
		outlist=[]	# the length of this will be the number of levels: (or contours?)
						# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
		levels=[]
		#contlevel=0
		#
		for i in xrange(len(cs)):
			# level-level:
			outlist+=[[]]
			#contcount=0
			# each level will have multiple polygons:
			for trace in cs[i].get_paths():
				# each "trace" will be a polygon (there might be multiple polygons per level), which might have "internal" polygons.
				tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				#outlist[-1]+=[[[],[]]]
				newpoly=[ [[], []] ]	# newpoly[poly-index][x=0 or y=1][x_i or y_i]. note, outer poly is first element.
				# draw the polygons (first the outside one, then the inside ones -- presumably)
				for lng,lat in trace.vertices:
					#
					# fix mangled contours:
					#
					newpoly[-1][0]+=[lng]
					newpoly[-1][1]+=[lat]
					#
					# have we closed the polygon (in the middle?):
					# what is happening is that the inner polys are being drawn screwy ways.
					# the errors are not a total mess; basically, the link between then inner and outer poly is drawn
					# from the last/first point rather than clostest two points.
					# the first poly is probably the outer-most poly. subsequent closed sequences are inside polys.
					# however, see comments above. the inner-inner->outer polys must be separated. do the following procedure
					# first, then pass that list of polys for further parsing.
					#
					# SO, make a list of all closed objects in a single poly object, then insert them into the outer sequence
					# (noting the inner poly has to close and then return to the starting point on the outer ring.
					# is newpoly[-1] == newpoly[0] ?
					if newpoly[-1][0][-1] == newpoly[-1][0][0] and newpoly[-1][1][-1] == newpoly[-1][1][0] and len(newpoly[-1][0])>1:
						# the polygon has closed. pinch it off and start a new poly.
						newpoly+=[ [[],[]] ]
						tmpLL=[]
						#contcount+=1
				#
				# now, clean up a bit:
				# if it's a bogus poly -- probably a line back to a prev. poly (each poly is like: poly0, line (segment) to poly1,
				#		poly1, either: {line back to poly0, line to poly2, nothing?}
				# i think the problem we're having is multi-inner polys aka, profile like 1 - 0 - 1. (ring=1, dip-ring=0, bump in middle=1)
				# ___|-\_/-\_/-|___ . the inner-inner bit is causing a problem. to fix this, i think, we have to draw all the polys
				# and explicitly determine which polys they are inside. (<outer><inner><inner></innter></inner></outer> ??)
				# note: a poly-ring inside an even  number of polys is an outer ring; a poly inside an odd number of rings is an inner ring.
				if len(newpoly[-1][0])<2: newpoly.pop()
				if newpoly[-1][0][-1] != newpoly[-1][0][0] and newpoly[-1][1][-1] != newpoly[-1][1][0]:
					newpoly.pop()
				#
				# at this point, newpoly is like: [[likely-outermost],[maybe-inner], [maybe-inner], ... ]
				# where we infer, or have in the past, that the first row represents the outer-most poly; subsequent polys are inners.
				# this is not the cast. find the 2xinner (aka, outer) polys within as per comments above.
				#
				# break up newpoly into inner-outer groups
				newpolylist=self.innerouterpolys(newpoly)
				#if len(newpoly)>2:
				#	print "len(newpoly), len(newpolylist): %d, %d" % (len(newpoly), len(newpolylist))
					#print newpoly
					#print newpolylist
				#newpolylist = [newpoly]
				# return the full list and interpret i=0 as the outer ring; i>0 as the inner boundaries.
				#
				#outlist[-1]+=[newpoly]	# add this trace(set) to the last level of outlist.
				for ply in newpolylist:
					#outlist[-1]+=[newpoly]
					outlist[-1]+=[ply]
				#
				#contcount+=1
			#contlevel+=1
		if getlevels==False:
			return outlist
		if getlevels==True:
			return [levels, outlist]
	
	def innerouterpolys(self, polylist):
		#print "do nothing yet."
		# but, eventually: separate polygons that are "inside" two two polys (aka, and "outside" poly), and
		# associate their inner polys with them.
		# use shapely.geometry (sgp) module?? maybe not, as per compatiblity.
		#
		# note, pollys are coming like [ [[Xs], [Ys]], [] ], 
		# so the jth vertex of the ith poly is (x,y) = (polylist[i][0][j], polylist[i][1][j])
		polylistplus=[]		# indexed entries: [polyindex, [list of inners], [verts] ]
		#outerpolys=[]			# these will be lists. each entry is like: [[outer],[true-inner],[true-inner],..]
		#
		#for ply in polylist:
		for i in xrange(len(polylist)):
			polylistplus += [[i, [], polylist[i]]]
			#
			# shorthand:
			# which poly each polygon it is inside, and how many polys it is inside (aka, len(that list)).
			#x0,y0=poly[i][0][0], poly[i][1][0]	# we only need to test one point since we're starting with contours (don't cross).
			x0,y0=polylistplus[-1][2][0][0], polylistplus[-1][2][1][0]	# we only need to test one point since we're starting with contours (don't cross).
			#print "x0, y0: ", x0, y0
			# in general, we'd need to test all points to be inside.
			#
			# for each polygon in this level:
			# is the ith ("top") polygon inside the j'th poly?
			for j in xrange(len(polylist)):
				if j==i: continue # (and also use exclusive inclusion (yi>y_test, not >=) to exclude self-insidedness).
				X,Y = polylist[j][0][:], polylist[j][1][:]
				#if x0>=max(X) or x0<=min(X) or y0>max(Y) or y0<min(Y): 
				#	print "outside max/min..."
				#	continue
				#
				if X[0]!=X[-1]:
					X+=[X[0]]
					Y+=[Y[0]]	# complete the poly...
				#
				N=len(X)
				ncrossings = 0
				# how many poly boundaries do we cross if we draw a line out of the poly in one direction.
				# equivalently (and in computer language), how many segments at y1 < y <y2 (or upside down)
			# are to the right of the point (or to the left, or up/down -- pick one).
				for k in xrange(1,N):
					k1 = k-1
					#k2 = (k+1)%N	# note the k%N: the poly does not have to be closed
					k2 = k	# but it should be, or you can count a crossing twice and get a bogus answer.
					x1,y1 = polylist[j][0][k1], polylist[j][1][k1]
					x2,y2 = polylist[j][0][k2], polylist[j][1][k2]
					#if y0>=min(y1, y2) and y0<=max(y1, y2) and x0<max(x1, x2):
					#
					if x0>=min(x1, x2) and x0<max(x1, x2):	# note, one must be <= and the other < or, if we're on a grid -- and
						fx = (y2-y1)/(x2-x1)*(x0-x1)			# we're always on a grid, we'll count each crossing (left seg then right).
						if fx>=(y0-y1):							# that's why the test script worked (x1,y1 != x2, y2) and the production
							ncrossings += 1						# failed...
				#
				#print i,j,j,ncrossings, ncrossings%2
				if ncrossings%2==1:
					# i'th poly is inside j'th poly...
					polylistplus[-1][1] += [j]	
				#
			#
			# so now, we have a list of polygons and the polys they are inside.
			# for outerPolys, len(innersList)%2==0. if polyA is inside polyB and nA-nB=1, polyA is an inner-poly to polyB
			#
			#for rw in polylistplus:
			#	print rw
			outerpolys=[]
			for ply1 in polylistplus:
				#print "inner-len: ", len(ply1[1]), len(ply1[1])%2
				if len(ply1[1])%2==0:
					#print "***", len(ply1[1])
					# it's an outer poly...
					outerpolys+=[[ply1[2]]]
					for ply2 in polylistplus:
						# find its inners:
						if ply2==ply1: continue
						if len(ply2[1])%2==0: continue	# skip outers...
						#
						#print len(ply2[1]), (len(ply1[1])+1)
						if ply1[0] in ply2[1] and len(ply2[1])==(len(ply1[1])+1):
							# the outer poly's index is in ply2's "inside-list", then ply2 is inside ply...
							# AND, ply2 is one "deeper" (inside exactly one more poly) thatn ply
							outerpolys[-1]+=[ply2[2]]
		#
		#return polylist
		return outerpolys
	#
	def contsToPlotListsOuter(self, cset=None):
		# simple polygons, not distinguishing inner/outer.
		if cset==None: cset=self.conts
		cs=cset.collections
		outlist=[]	# the length of this will be the number of levels: (or contours?)
						# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
		contlevel=0
		#
		for i in xrange(len(cs)):
			# level-level:
			outlist+=[[]]
			contcount=0
			for trace in cs[i].get_paths():
				tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				outlist[-1]+=[[[],[]]]
				for lng,lat in trace.vertices:
					#
					# fix mangled contours:
					#
					outlist[-1][-1][0]+=[lng]
					outlist[-1][-1][1]+=[lat]
					tmpLL+=[[lng, lat]]
					#
					# have we closed the polygon (in the middle?):
					# what is happening is that the inner polys are being drawn screwy ways.
					# the errors are not a total mess; basically, the link between then inner and outer poly is drawn
					# from the last/first point rather than clostest two points.
					# the first poly (i think) is the outer poly. subsequent closed sequences are inner polys (excluded regions).
					# SO, make a list of all closed objects in a single poly object, then insert them into the outer sequence
					# (noting the inner poly has to close and then return to the starting point on the outer ring.
					if outlist[-1][-1][0][-1] == outlist[-1][-1][0][0] and outlist[-1][-1][1][-1] == outlist[-1][-1][1][0] and len(outlist[-1][-1][0])>1:
						# the polygon has closed. pinch it off and start a new poly.
						outlist[-1]+=[[[],[]]]
						tmpLL=[]
						contcount+=1

				if len(outlist[-1][-1][0])<2: outlist[-1].pop()
				if outlist[-1][-1][0][-1] != outlist[-1][-1][0][0] and outlist[-1][-1][1][-1] != outlist[-1][-1][1][0]:
					outlist[-1].pop()
				# (outlist[contour level][polynum][x or y][x,y index] )
				
				
				'''
				#if tmpLL[0]!=tmpLL[-1]:
				while tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
					#fixedpolys+=1
					#outlist[-1][-1][0] += [tmpLL[0][0]]
					#outlist[-1][-1][1] += [tmpLL[0][1]]
					outlist[-1][-1][0].pop()
					outlist[-1][-1][1].pop()
					tmpLL.pop()
				'''
					#
				contcount+=1
			contlevel+=1
		return outlist

	def plotPolyList(self, polys=None, markerstr='.-', fignum=3, pllevels=None):
		if polys==None: polys=self.contsToPlotLists()
		if pllevels==None: pllevels = range(len(polys))
		if type(pllevels).__name__ in ('float', 'int'): pllevels=[pllevels]
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		nlevels=len(polys)
		
		#
		ilevel=0
		for level in polys:
			#lvlclr=self.plotcolor(ilevel, nlevels)
			if ilevel not in pllevels:
				ilevel+=1
				continue
				
				 
			lvlclr=plotcolor(ilevel, nlevels*2)
			print ilevel, lvlclr
			for xy in level:
				for ring in xy:
				#fixedpolys = checkPoly(inpoly=xy, fignum=None)
				#for fxy in fixedpolys:
					plt.plot(ring[0], ring[1], markerstr, color=lvlclr)
				#plt.plot(xy[0], xy[1], markerstr, color=lvlclr)
			ilevel+=1
		return polys			
	
	def plotPolyListOuter(self, polys=None, markerstr='.-', fignum=3, pllevels=None):
		# plot the simpler "just outers" version of the polygon list.
		if polys==None: polys=self.contsToPlotListsOuter()
		if pllevels==None: pllevels = range(len(polys))
		if type(pllevels).__name__ in ('float', 'int'): pllevels=[pllevels]
		plt.figure(fignum)
		plt.clf()
		plt.ion()
		nlevels=len(polys)
		
		#
		ilevel=0
		for level in polys:
			#lvlclr=self.plotcolor(ilevel, nlevels)
			if ilevel not in pllevels:
				ilevel+=1
				continue
				
				 
			lvlclr=plotcolor(ilevel, nlevels*2)
			print ilevel, lvlclr
			for xy in level:
				#fixedpolys = checkPoly(inpoly=xy, fignum=None)
				#for fxy in fixedpolys:
				#	plt.plot(fxy[0], fxy[1], markerstr, color=lvlclr)
				plt.plot(xy[0], xy[1], markerstr, color=lvlclr)
			ilevel+=1
		return polys
	
	def getmarkerKMLstr(self, markers):
		#kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
		kmlstr=''
		for ev in markers:
			thisdt=str(mpd.num2date(ev[0]))
			datestr=thisdt.split('.')[0]
			kmlstr+='<Placemark><name>%s, %.2f</name><Point><coordinates>%f,%f,0</coordinates></Point></Placemark>' % (datestr, ev[3], ev[2], ev[1])
		
		return kmlstr
		
	def writeKMLfile(self, cset=None, fout='conts.kml', colorbarname='scale.png', equakes=None):
		#kmlstr=self.KMLstrFromConts(cset)
		kmlstr=self.KMLstrFromConts1(cset, openfile=True, closefile=False)
		if equakes !=None:
			# we've been given some earthquakes to plot as well.
			kmlstr+='\n'
			kmlstr+=self.getmarkerKMLstr(equakes)
		kmlstr+='</Document>\n'
		kmlstr+='</kml>'
		#
		f=open(fout, 'w')
		f.write(kmlstr)
		f.close()
		return None
	
	def KMLstrFromConts(self, cset=None, colorbarname='scale.png', openfile=True, closefile=True, warnings=None):
		# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
		# (will it be necessary to write directly to file? is there a string length limit problem?)
		#
		if cset==None:
			cset=self.conts
			#cset=plt.contourf(self.X_i, self.Y_i, self.Z2d, self.contres*self.contfact, cm=plt.spectral())
		cs=cset.collections
		polys = self.contsToPlotLists(cset=cset)	# this function fixes multi-closed polys. use for actual kml contours.
																# still has same polys[cont. level][contour index][x=0,y=1][x,y index] format
																# the polys array will differ from cset.collections in the number of polys
																# per level; the number of levels should be the same.
		#
		#resolution = 5
		#LevelsNumber = 5 * resolution
		if warnings==None: warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#contoursStart=int(.2*len(cs))
		#resolution = int(len(cs)/len(warnings))
		resolution = int(len(polys)/len(warnings))
		#startindex=resolution
		startindex=0
		kmlstr=''
		#
		if openfile==True:
			kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
		#
		# styles come from the contours collection (cs):
		for i in range(0,len(cs)):
			bgra_array = 255*cs[i].get_facecolor()[0]
			kmlstr+='<Style id="l%d">\n' % i
			kmlstr+='<LineStyle><color>00000000</color></LineStyle>\n'
			kmlstr+='<PolyStyle><color>7d%02x%02x%02x</color></PolyStyle>\n' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
			kmlstr+='</Style>\n'
		#
		#
		kmlstr+='<ScreenOverlay id="scale">\n'
		kmlstr+='<name>Color Scale</name>\n'
		#kmlstr+='<Icon><href>scale.png</href></Icon>\n'
		kmlstr+='<Icon><href>%s</href></Icon>\n' % colorbarname
		kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
		kmlstr+='</ScreenOverlay>\n'
		#
		#fixedpolys=0
		# now, stake out the contours from the polys array.
		# note: len(cs) = len(polys) (same number of levels)
		#       len(cs[i]) != len(polys[i]) (in fact, len(polys[i]>=len(cs[i]) ), which is to say
		#       polys may have more distinct polygons per contour level. of course, the individual polys
		#       will have differend length as well.
		#for i in xrange(startindex, len(cs)):
		for i in xrange(startindex, len(polys)):
			# each i is a contour level.
			kmlstr+='<Placemark>\n'
			#print i, resolution, len(warnings), len(polys)
			warningindex=int(float(i)/float(resolution))
			if warningindex>(len(warnings)-1): warningindex=len(warnings)-1	# in the event that we have off-integer numbers.
			kmlstr+='<name>%s Risk</name>\n' % warnings[warningindex]
			kmlstr+='<styleUrl>#l%d</styleUrl>\n' % i
			kmlstr+='<MultiGeometry>\n'
			 
			#for trace in cs[i].get_paths():
			for ii in xrange(len(polys[i])):
				# each ii is a polygon (set).
				kmlstr+='<Polygon>\n'
				kmlstr+='<extrude>0</extrude>\n'
				kmlstr+='<altitudeMode>clampToGround</altitudeMode>\n'
				kmlstr+='<outerBoundaryIs>\n'
				kmlstr+='<LinearRing>\n'
				kmlstr+='<coordinates>\n'
				
				#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				#for lng,lat in trace.vertices:
				for ill in xrange(len(polys[i][ii][0][0])):
					# first set is the outerBoundary
					# noting that each polygon is stored as [[x], [y]], so len(polys[i][ii])=2 always.
					# len(polys[i][ii][0])=len(polys[i][ii][1]) is the length or size (num. verts.) of the polygon.
					lng=polys[i][ii][0][0][ill]
					lat=polys[i][ii][0][1][ill]
					#
					kmlstr+='%f,%f,0\n' % (lng, lat)
					#tmpLL+=[[lng, lat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])

				kmlstr+='</coordinates>\n</LinearRing>\n</outerBoundaryIs>\n'
				#
				# inner polys?
				for iInner in xrange(1,len(polys[i][ii])):
					thispoly=polys[i][ii][iInner]
					# if any, these will be inner polys, each like [[x], [y]]
					kmlstr+='<innerBoundaryIs>\n<LinearRing>\n<coordinates>\n'
					# coords like 'lon,lat,alt\n'
					# -77.05668055019126,38.87154239798456,100
					for ill in xrange(len(thispoly[0])):
						kmlstr+='%f,%f,0\n' % (thispoly[0][ill], thispoly[1][ill])
					kmlstr+='</coordinates>\n</LinearRing>\n</innerBoundaryIs>\n'
				#
				kmlstr+='</Polygon>\n'
			#
			kmlstr+='</MultiGeometry>\n'
			kmlstr+='</Placemark>\n'
		#
		if closefile==True:
			kmlstr+='</Document>\n'
			kmlstr+='</kml>'
		#
		return kmlstr

	def KMLstrFromConts1(self, cset=None, colorbarname='scale.png', openfile=True, closefile=True, warnings=None):
		# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
		# (will it be necessary to write directly to file? is there a string length limit problem?)
		#
		if cset==None:
			cset=self.conts
			#cset=plt.contourf(self.X_i, self.Y_i, self.Z2d, self.contres*self.contfact, cm=plt.spectral())
		cs=cset.collections
		polys = self.contsToPlotLists(cset=cset)	# this function fixes multi-closed polys. use for actual kml contours.
																# still has same polys[cont. level][contour index][x=0,y=1][x,y index] format
																# the polys array will differ from cset.collections in the number of polys
																# per level; the number of levels should be the same.
		#
		#resolution = 5
		#LevelsNumber = 5 * resolution
		if warnings==None: warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#contoursStart=int(.2*len(cs))
		#resolution = int(len(cs)/len(warnings))
		resolution = int(len(polys)/len(warnings))
		#startindex=resolution
		startindex=0
		kmlstr=''
		#
		if openfile==True:
			kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
		#
		# styles come from the contours collection (cs):
		for i in range(0,len(cs)):
			bgra_array = 255*cs[i].get_facecolor()[0]
			kmlstr+='<Style id="l%d">\n' % i
			kmlstr+='<LineStyle><color>00000000</color></LineStyle>\n'
			kmlstr+='<PolyStyle><color>7d%02x%02x%02x</color></PolyStyle>\n' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
			kmlstr+='</Style>\n'
		#
		#
		kmlstr+='<ScreenOverlay id="scale">\n'
		kmlstr+='<name>Color Scale</name>\n'
		#kmlstr+='<Icon><href>scale.png</href></Icon>\n'
		kmlstr+='<Icon><href>%s</href></Icon>\n' % colorbarname
		kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
		kmlstr+='</ScreenOverlay>\n'
		#
		#fixedpolys=0
		# now, stake out the contours from the polys array.
		# note: len(cs) = len(polys) (same number of levels)
		#       len(cs[i]) != len(polys[i]) (in fact, len(polys[i]>=len(cs[i]) ), which is to say
		#       polys may have more distinct polygons per contour level. of course, the individual polys
		#       will have differend length as well.
		#for i in xrange(startindex, len(cs)):
		for i in xrange(startindex, len(polys)):
			# each i is a contour level.
			kmlstr+='<Placemark>\n'
			#print i, resolution, len(warnings), len(polys)
			warningindex=int(float(i)/float(resolution))
			if warningindex>(len(warnings)-1): warningindex=len(warnings)-1	# in the event that we have off-integer numbers.
			kmlstr+='<name>%s Risk</name>\n' % warnings[warningindex]
			kmlstr+='<styleUrl>#l%d</styleUrl>\n' % i
			kmlstr+='<MultiGeometry>\n'
			 
			#for trace in cs[i].get_paths():
			for ii in xrange(len(polys[i])):
				# each ii is a polygon (set).
				kmlstr+='<Polygon>\n'
				kmlstr+='<extrude>0</extrude>\n'
				kmlstr+='<altitudeMode>clampToGround</altitudeMode>\n'
				kmlstr+='<outerBoundaryIs>\n'
				kmlstr+='<LinearRing>\n'
				kmlstr+='<coordinates>\n'
				
				#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				#for lng,lat in trace.vertices:
				for ill in xrange(len(polys[i][ii][0][0])):		#[level][trace][outer(top)-poly][Xcoords] (last one could be 1)
					# first set is the outerBoundary
					# noting that each polygon is stored as [[x], [y]], so len(polys[i][ii])=2 always.
					# len(polys[i][ii][0])=len(polys[i][ii][1]) is the length or size (num. verts.) of the polygon.
					lng=polys[i][ii][0][0][ill]
					lat=polys[i][ii][0][1][ill]
					#
					kmlstr+='%f,%f,0\n' % (lng, lat)
					#tmpLL+=[[lng, lat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])

				kmlstr+='</coordinates>\n</LinearRing>\n</outerBoundaryIs>\n'
				#
				# inner polys?
				for iInner in xrange(1,len(polys[i][ii])):
					thispoly=polys[i][ii][iInner]
					# if any, these will be inner polys, each like [[x], [y]]
					kmlstr+='<innerBoundaryIs>\n<LinearRing>\n<coordinates>\n'
					# coords like 'lon,lat,alt\n'
					# -77.05668055019126,38.87154239798456,100
					for ill in xrange(len(thispoly[0])):
						kmlstr+='%f,%f,0\n' % (thispoly[0][ill], thispoly[1][ill])
					kmlstr+='</coordinates>\n</LinearRing>\n</innerBoundaryIs>\n'
				#
				kmlstr+='</Polygon>\n'
			#
			kmlstr+='</MultiGeometry>\n'
			kmlstr+='</Placemark>\n'
		#
		if closefile==True:
			kmlstr+='</Document>\n'
			kmlstr+='</kml>'
		#
		return kmlstr
		
	#
	def KMLstrFromConts_raw(self, cset=None, colorbarname='scale.png'):
		# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
		# (will it be necessary to write directly to file? is there a string length limit problem?)
		#
		if cset==None:
			cset=self.conts
			#cset=plt.contourf(self.X_i, self.Y_i, self.Z2d, self.contres*self.contfact, cm=plt.spectral())
		cs=cset.collections
		#
		#resolution = 5
		#LevelsNumber = 5 * resolution
		warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#contoursStart=int(.2*len(cs))
		resolution = int(len(cs)/self.contfact)
		#startindex=resolution
		startindex=0
		#
		#
		kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
		#
		for i in range(0,len(cs)):
			bgra_array = 255*cs[i].get_facecolor()[0]
			kmlstr+='<Style id="l%d">\n' % i
			kmlstr+='<LineStyle><color>00000000</color></LineStyle>\n'
			kmlstr+='<PolyStyle><color>7d%02x%02x%02x</color></PolyStyle>\n' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
			kmlstr+='</Style>\n'
		#
		#
		kmlstr+='<ScreenOverlay id="scale">\n'
		kmlstr+='<name>Color Scale</name>\n'
		#kmlstr+='<Icon><href>scale.png</href></Icon>\n'
		kmlstr+='<Icon><href>%s</href></Icon>\n' % colorbarname
		kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
		kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
		kmlstr+='</ScreenOverlay>\n'
		#
		fixedpolys=0
		for i in xrange(startindex, len(cs)):
			kmlstr+='<Placemark>\n'
			kmlstr+='<name>%s Risk</name>\n' % warnings[i/resolution]
			kmlstr+='<styleUrl>#l%d</styleUrl>\n' % i
			kmlstr+='<MultiGeometry>\n'
			 
			for trace in cs[i].get_paths():
				kmlstr+='<Polygon>\n'
				kmlstr+='<extrude>0</extrude>\n'
				kmlstr+='<altitudeMode>clampToGround</altitudeMode>\n'
				kmlstr+='<outerBoundaryIs>\n'
				kmlstr+='<LinearRing>\n'
				kmlstr+='<coordinates>\n'
				
				tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				for lng,lat in trace.vertices:
					#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
					
					# fixing contours mangled by pyplot:
					# an older attempted fixing routine that didn't quite work out.
					# keep this version as the "raw" form; see other KML function that 
					# uses a repaired set of polygons.
					#
					thislng=round(lng, self.zresolution)
					thislat=round(lat, self.zresolution)
					kmlstr+='%f,%f,0\n' % (thislng, thislat)
					tmpLL+=[[thislng, thislat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])

				kmlstr+='</coordinates>\n'
				kmlstr+='</LinearRing>\n'
				kmlstr+='</outerBoundaryIs>\n'
				kmlstr+='</Polygon>\n'
			#
			kmlstr+='</MultiGeometry>\n'
			kmlstr+='</Placemark>\n'
		#
		kmlstr+='</Document>\n'
		kmlstr+='</kml>'
		#
		return kmlstr
	
	def makeColorbar(self, cset=None, colorbarname=None, reffMag=5.0):
		#
		# reffMag:
		if reffMag==None: reffMag=self.reffMag
		#
		if colorbarname==None: colorbarname='scale.png'
		if cset==None: cset=self.conts
		cs=cset.collections
		#resolution = 5
		#LevelsNumber = 5 * resolution
		warnings = ['No','Low','Guarded','Elevated','High','Severe']
		#contoursStart=int(.2*len(cs))
		resolution = int(len(cs)/self.contfact)
		#startindex=resolution
		startindex=0
		#
		fig = plt.figure(figsize=(1,3))
		axe = fig.add_axes([0.05, 0.05, 0.2, 0.9])
		#
		fig.figurePatch.set_alpha(0)
		#
		matplotlib.rc('font', family='monospace', weight='black')
		#
		cmap = mpl.cm.spectral
		#norm = mpl.colors.Normalize(vmin=Z_i.min(), vmax=Z_i.max())
		#ratefactor=10**(self.mc-self.reffMag)
		ratefactorExp=self.mc-reffMag	# so the raw rate probably should be with mc_catalog = mc_etas.
													# but, we'll want to report some rate of m>reffMag
		norm = mpl.colors.Normalize(vmin=self.Z2d.min(), vmax=self.Z2d.max())
		#
		#print self.Z2d.min(), self.Z2d.max(), norm
		#tics = [0,10,20,40,80]
		#tics = [round(self.Z2d.min(),3), round(self.Z2d.max(),3)]
		timefactExp = math.log10(year2secs)
		t1='%.2f' % (self.Z2d.min() + timefactExp + ratefactorExp)
		t2='%.2f' % (self.Z2d.max() + timefactExp + ratefactorExp)
		print "max set: ", self.Z2d.max(), timefactExp, ratefactorExp
		tics = [self.Z2d.min(), self.Z2d.max()]
		#print t1, t2
		#
		#cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%g%%", orientation='vertical')
		cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%.2f", orientation='vertical')
		cb1.set_ticklabels([t1, t2], update_ticks=True)
		cb1.set_label('ETAS rate: $\\log({N(m \geq %.1f)}/{yr \cdot km^2})$' % reffMag)
		#cb1.set_label('ETAS rate')
		plt.savefig(colorbarname)	# vertical
		#
		#suffix=''
		i=len(colorbarname)-1
		#while colorbarname[i]!='.':
		#	#suffix.insert(0, colorbarcopy[-i])
		#	# this will get us the last '.'
		#	i-=1
		i=colorbarname.rfind('.')
		hname=colorbarname[:i] + '-h' + colorbarname[i:]
		#
		im = ipp.open(colorbarname, 'r')
		im.load()							# just to be sure. open() only creates a pointer; data are not loaded until analyzed.
		im=im.rotate(-90)
		#im.show()
		im.save(hname)
		#
		# depricated method:
		#os.system('convert -rotate 90 %s %s' % (colorbarname, hname))
		#os.system('cp %s %s' % (colorbarname, hname))
		#
		return [colorbarname, hname]
	#
	def writeCat(self, fout='cat.cat'):
		# dt, lat, lon, mag
		f=open(fout, 'w')
		f.write('#BASScast catalog\n#dtm, lat, lon, mag\n')
		f.write('#!writedate\t%s\n' % str(dtm.datetime.now(pytz.timezone('UTC'))))
		f.write('#!lats=%s\tlons=%s\tfcdate=%s\t%f\n' % (str(self.latrange), str(self.lonrange), self.fcdate, self.fcdatef))
		for rw in self.catalog:
			f.write('%f\t%f\t%f\t%f\n' % (rw[0], rw[1], rw[2], rw[3]))
		f.close()

class locbase(object):
	loc = [0,0]		# [lon, lat]
	eventDtm = None
	eventftime = None # float version of time
	mc=2.0
	#
	# location related functions:
	def timesince(self, intime):
		# return time since the in-event (in seconds):
		#print type(self)
		# this will be a guess until "type()" can be properly implemented, but...
		if type(intime)==type(self):
			# we've been asked to compare with an earthquake object.
			# but the way we SHOULD do this is to check for a member function/variable...
			intime=intime.eventftime
		#
		# otherwise, assume float ops will work.
		tsince = intime - self.eventftime
		#
		return tsince
		
	def dispFrom(self, inloc, coordtype='cart'):
		# displacement from inloc...
		# inloc is a vector [lon, lat]
		# return a vector [dLon, dLat] or [r, theta]
		# return distances in km.
		#
		#meanlat=.5*(inloc[1]+self.loc[1])*2.0*math.pi/360.0
		latfactSelf=math.cos(self.loc[1]*2.0*math.pi/360.0)
		latfactin = math.cos(inloc[1]*2.0*math.pi/360.0)
		#
		#dx = (inloc[0] - self.loc[0])*111.12*latfact
		dx = (inloc[0]*latfactin - self.loc[0]*latfactSelf)*111.12
		dy = (inloc[1] - self.loc[1])*111.12
		rvect=[dx, dy]
		#
		if coordtype=='radial' or coordtype=='rad':
			R=math.sqrt(dx*dx + dy*dy)
			if R==0:
				# should not happen often except in demos. just return 0.
				theta=0.0
			else:
				theta = math.acos(dx/R)
			rvect=[R,theta]
		#
		return rvect
	
	def distFrom(self, inloc, coordtype='cart'):
		return self.dispFrom(inloc, coordtype)

	def setEventTime(self, evtime):
			#print type(evtime).__name__
		if type(evtime).__name__ == 'datetime': # == type(dtm.datetime):
			self.eventDtm = evtime
			self.eventftime = mpd.date2num(evtime) * days2secs	# catalog must be in seconds for scaling to work.
			#self.eventftime=self.datetimeToFloat(evtime)
		if type(evtime) == type(1.0) or type(evtime) == type(1):
			self.eventftime = float(evtime)
			#
			if self.eventftime < (mpd.date2num(dtm.datetime.now(pytz.timezone('UTC')))*100.0):
				# date is in days, not seconds
				self.eventftime*=days2secs
				#
			self.eventDtm = mpd.num2date(evtime/days2secs)	# assuming it is properly in seconds, as it should be.
			# self.eventDtm = self.datetimeFromFloat(evtime)
		if type(evtime) == type('abc'):
			# convert the string, etc.
			#self.eventDtm = self.datetimeFromString(strDtin=evtime)
			self.eventftime = mpd.datestr2num(evtime)*days2secs
			self.eventDtm = mpd.num2date(self.eventftime/days2secs)			
			# mpd.num2date
			#self.eventftime = self.datetimeToFloat(self.eventDtm)
			#self.eventftime = mpd.date2num(self.eventDtm)
		return None

	# utils we probably don't need any longer (use matplotlib.date.{various date conversions}

class objContour(object):
	lwr=None
	upr=None
	#clr='7d000000'	# white? (alpha - R - G - B)
	alpha = '7d'
	R='00'
	G='00'
	B='00'
	
	verts=[]
	
	def __init__(self, low=None, high=None, alpha=None, R=None, G=None, B=None, verts=None, RGB=None):
		return self.initialize(low=low, high=high, alpha=alpha, R=R, G=G, B=B, verts=verts, RGB=RGB)
		#
	#
	def initialize(self, low=None, high=None, alpha=None, R=None, G=None, B=None, verts=None, RGB=None):
		#print "RGBinit: %s" % RGB
		if RGB!=None:
			self.setRGB2(RGB=RGB)
		else:
			self.setRGB(R=R, G=G, B=B, alpha=alpha)
		self.lwr=low
		self.upr=high
		#
		self.verts = verts
		#
		return None
		#
	#
	def getRGB(self):
		return '%s%s%s%s' % (self.alpha, self.R, self.G, self.B)
	
	def setRGB(self, R='00', G='00', B='00', alpha='7d'):
		R='00' + str(R)
		G='00' + str(G)
		B='00' + str(B)
		alpha='00' + str(alpha)
		#
		R=R[-2:]
		G=G[-2:]
		B=B[-2:]
		alpha=alpha[-2:]
		#
		self.R=R
		self.G=G
		self.B=B
		self.alpha=alpha
		#
		return None
	
	def setRGB2(self, RGB='7d000000'):
		#if type(RGB).__name__ != 'string': return None
		#print 'RGB: %s' % RGB
		#
		self.alpha=RGB[0:2]
		self.R=RGB[2:4]
		self.G=RGB[4:6]
		self.B=RGB[6:8]
		#
		return None
	
	def aryRow(self):
		return [self.lwr, self.upr, self.getRGB(), self.verts]

class forecastsite(locbase):
	# site, or grid-cell, for a forecast.
	# basically, it needs to hava location (position and time), size, and z-value (forecast metric).
	# nominally, it should also have some catalog prams, like mc to determine an earthquake count.
	# it does beg the question as to whether scaling prams like rho, alpha, etc. are properties of the
	# earthquake or the local earth (aka, earthquake object or forecastsite object).
	#
	#loc = [0,0]		# [lon, lat]
	#eventDtm = None
	#eventftime = None # float version of time
	#mc=2.0
	#z=1.0
	z=None
	#
	#
	def __init__(self, loc=[0.0,0.0], dxdy=[.1, .1], evtime=1.0, mc=2.0):
		return self.initialize(loc=loc, dxdy=dxdy, evtime=evtime, mc=mc)
	
	def initialize(self, loc=[0.0,0.0], dxdy=[.1, .1], evtime=1.0, mc=2.0):
		# note scaling exponents are pre-declared.
		self.loc=loc
		self.latfactor = math.cos(2.0*math.pi*loc[1]/360.0)
		#
		for i in xrange(len(self.loc)):
			self.loc[i]=float(self.loc[i])
		# element size:
		if type(dxdy).__name__==type(.1).__name__:
			dxdy=[dxdy, dxdy]
		for i in xrange(len(dxdy)):
			dxdy[i]=float(dxdy[i])
		self.dxdy=dxdy	# note: these values are in lat/lon at this point.
		#
		self.mc=float(mc)
		#
		self.setEventTime(evtime)
		#
		#
		return None
	
	def getIntensity(self, equake):
		# dN/dt for this site. localIntensity() returns a value like n/[sec][km^2]. here, we itegrate that over area.
		dz=equake.localIntensity(self.loc, self.eventftime)
		#
		# note this is an approximation in which dz/dA is calculated from circular coords. and now we interpolate using rect. coords.
		dZ = dz*self.dxdy[0]*self.dxdy[1]*self.latfactor # though a more accurate calc. can be made by doing the integral (where the
																			# the approximation dA = drdr,
																			# Ndot = (orate/2*pi*r0^(-(1+rho))) * (1/(rho*(rho-1)) * [r2^{1-rho} - r1^{1-rho}]
		return dZ
	
	def getIntensities(self, equakes):
		# for this individual site (presumably, for all sites but righ here for this site), calculate
		# contribution from all earthquakes.
		# note that this is not integrated -- z ~ n/[sec][km^2]
		dzs=[]
		for equake in equakes:
			dzs+=[equake.localIntensity(self.loc, self.eventftime)]
		return dzs
	
	def setz(self, equakes):
		#print 'setting z.'
		z = sum(self.getIntensities(equakes))
		self.z=z
		return z
	
class earthquake(locbase):
	# (probably should have made forecastsite() first and then derived earthquake from it).
	#
	# earthquake parameters (purely observable):
	loc = [0.0,0.0]		# [lon, lat]
	mag=None
	eventDtm = None
	eventftime = None # float version of time
	#mc=2.0	# obviously, this is not an earthquake propertie, but we will need it to calculate rates, etc.
				# another approach might be to return self-rates, aka, the rate at which an earthquake produces
				# itself. this could be converted to catalog rates, of course, by subtracting mc.
	#
	# scaling exponents:
	rho=1.37				# spatial scaling exponent
	#p = 1.1  			# Omori scaling exponent
	sig = 1.68			# aftershock surface density exponent.
	#alpha = 2.28		# rupture duration scaling exponent, nominally 
	b=1.0
	#dmstar=1.2
	dl0=1.0
	#drho=5.25
	
	# calculated earthquake properties:
	ltau = None  # Omori pram (log[tau]) in seconds
	lrupttime=None	# duration of rupture (log[rupture duration]) in seconds
	lt0 = None	# Omori pram derrived from rupt. duration (typically log(t0) = lrupttime + 1.0)
	lruptlen = None	# log rupt. len (in km).
	#
	def __init__(self, mag=5.0, loc=[0.0,0.0], evtime=1.0, mc=2.0, drho=5.25, p=1.1, dmstar=1.2, dm=1.0, alpha = 2.28):
		return self.initialize(mag=mag, loc=loc, evtime=evtime, mc=mc, drho=drho, p=p, dmstar=dmstar, alpha=alpha)
	
	def initialize(self, mag=5.0, loc=[0.0,0.0], evtime=1.0, mc=2.0, drho=5.25, p=1.1, dmstar=1.2, dm=1.0, alpha = 2.28):
		# note scaling exponents are pre-declared.
		self.drho=float(drho)
		self.p=float(p)
		self.dmstar=dmstar	# mag difference for a primary sequence (a mainshock and its aftershocks)
		self.dm=dm				# mag difference for a full (compound) aftershock sequence
		self.mag=float(mag)
		self.loc=loc
		self.alpha=alpha
		for i in xrange(len(self.loc)):
			self.loc[i]=float(self.loc[i])
		#
		self.mc=float(mc)
		#
		nothin=self.setEventTime(evtime)
		'''
		if type(evtime).__name__ == 'datetime': # == type(dtm.datetime):
			self.eventDtm = evtime
			self.eventftime = nbd.date2num(evtime)
			#self.eventftime=self.datetimeToFloat(evtime)
		if type(evtime) == type(1.0):
			self.eventftime = evtime
			self.eventDtm = mpd.num2date(evtime)
			# self.eventDtm = self.datetimeFromFloat(evtime)
		if type(evtime) == type('abc'):
			# convert the string, etc.
			#self.eventDtm = self.datetimeFromString(strDtin=evtime)
			self.eventftime = mpd.datestr2num(evtime)
			self.eventDtm = mpd.num2date(self.eventftime)			
			# mpd.num2date
			#self.eventftime = self.datetimeToFloat(self.eventDtm)
			#self.eventftime = mpd.date2num(self.eventDtm)
		'''
		#
		self.dmp=1.0-math.log10(self.p-1.0)/self.b	# scaling exponent for initial rate 1/tau?
		#
		# and of course, here is the debate. how do we calc. tau/t0? this method might work actually...
		# the main problem from before was that the catalog was not being converted from days -> seconds.
		# we do, however, introduce a new method of calculating t0/tau. rather than use the theoretical maximum
		# (GR) rupture rate, we use the event-resolution time (the time when the omori-interval is greater than
		# the rupture duration (or a factor thereof).
		#
		#self.ltau=self.ltauYoder(m=mag, alpha=self.alpha, p=self.p, b=self.b, dmp=self.dmp, mc=mc)
		#self.lt0 = self.getlt0(m=mag, alpha=self.alpha)
		self.lt0 = self.getlt02 (mstar=mag, alpha=self.alpha, dmstar=self.dmstar, drho=self.drho)
		self.ltau = self.getltau2(mstar=mag, mc=self.mc, alpha=self.alpha, p=self.p, drho=self.drho, thislt0=self.lt0)
		self.lrupttime=self.lDeltat(m=mag, alpha=self.alpha)
		self.lruptlen=self.getlruptLen(m=mag, dl0=None)
		#
		return None
	#
	# earthquake utility functions:
	def omoriRate(self, t, t0=None, tau=None, p=None):
		if t0==None: t0=10.0**self.lt0
		if tau==None: tau=10.0**self.ltau
		if p==None: p=self.p
		#
		rinv=tau*(t0+t)**p	# "interval" form of Omori
		tfact=1.0				# don't remember what this is for... some factor of the rate.
		#tfact=5.0
		return tfact/rinv
	#
	def omoriN(self, t1=0., t2=None, t0=None, tau=None, p=None):
		# integral of omoriRate (N occurred since t0) for p!=1.0 (p>1.0).
		if t0==None: t0=10.0**self.lt0
		if tau==None: tau=10.0**self.ltau
		if p==None: p=self.p
		if t2==None: t2=t1+t0	# or should we let this error off?
		#
		if p==1.0:
			return omoriNln(t=t, t0=t0, tau=tau)
		#
		N=( (t0 + t2)**(1.0-p) - (t0 + t1)**(1.0-p) )/(tau*(1.0-p))
		#
		return N
	#			
	def omoriNln(self, t1=0., t2=None, t0=None, tau=None):
		# integral of omoriRate (N occurred since t0) with p=1.0
		if t0==None: t0=10.0**self.lt0
		if tau==None: tau=10.0**self.ltau
		if t2==None: t2=t1+t0	# or should we let this error off?
		p=1.0
		#
		#N=(math.log((1.0 + t/t0))/tau	# note: log(x) = ln(x)
		N=(math.log((t2 + t0)/(t1 + t0)))/tau
		#
		return N
	#	
	def radialDensity(self, m=None, r=0.0, lrupt=None, rho=None, dm=None):
		# lrupt -> l_rupture (rupture length). l_rupture = 10**lruptlen
		if lrupt==None: 
			lrupt=10.0**self.lruptlen
			if dm==None: dm=self.dm
		#
		if m==None: m=self.mag
		if rho==None: rho=self.rho
		if dm==None: dm=self.dm
		Nom=10.0**(m-dm-self.mc)
		#
		#note: is this r=0 at the edge of the rupt. or r=0 at the center?
		# the (r+r0) formulation implies the "scaling" distance. we want just plain distance, and assume
		# the rate is constant for 0 < r < lrupt
		#
		#
		# omori-like formulation. r is the "scaling" distance; r=0 indicates the edge of the rupture.
		# the "scaling" distance r' = r - lrupt
		#rate = (rh0 - 1.0)*(lrupt**(rho-1.0))*Nom*(lrupt + r)**(-rho)
		#
		# equivalently, if r is the dist. from the epicenter,
		if r<lrupt: r=lrupt		# note: this will produce an upper-truncated distribution (which is good).
		rate = (rho-1.0)*(lrupt**(rho-1.0))*Nom*r**(-rho)
		# noting that scaling will start at r=lrupt, so effectively r = lrupt + r'
		#
		return rate	
	#	
	def dmprime(self, p=None, b=None):
		# scaling decay exponent for initial omori rate (?), 1/tau:
		if p==None: p=self.p
		if b==None: b=self.b
		dmp = 1.0-math.log10(p-1.0)/b	# usually dmp~2.0
		#
		return dmp
	
	def getruptlen(self, m=None, dl0=None):
		if m==None: m=self.mag
		if dl0==None: dl0=self.dl0
		#
		return 10.0**self.getlruptLen(m,dl0)
	
	def getlruptLen(self, m=None, dl0=None):
		if dl0==None: dl0=self.dl0
		if m==None:
			m=self.mag
		return (m/2.0) - 1.755 - dl0

	#	method1: assuming some maximum rate at the end of rupture rate = N_GR/Delta-t_R * (some factor)
	def gett0(self, lt0=None):
		if lt0==None:
			lt0=self.lt0
		return 10**lt0
	
	def gettau(self, ltau=None):
		if ltau==None:
			ltau=self.ltau
		return 10**ltau
	
	def ltauYoder(self, m=None, alpha=None, p=None, b=1.0, dmp=None, mc=None):
		if mc==None: mc=self.mc
		if p==None: p=self.p
		if dmp==None: dmp=self.dmp
		if m==None: m=self.mag
		if alpha==None: alpha=self.alpha
		if p==None: p=self.p
		if b==None: b=self.b
		if dmp==None: dmp=self.dmp
		if mc==None: mc=self.mc
		#
		# actualy log(tau)
		ltau=alpha*(p-1.0) + b*(dmp+mc) + m*((1.0-p)/2.0 - b)
		#
		return ltau

	def lDeltat(self, m=None, alpha=None):
		if m==None: m=self.mag
		if alpha==None: alpha=self.alpha
		
		return (m/2.0) - alpha

	def getlt0(self, m=None, alpha=None):
		# experimentally, 10*deltaT seem to make nice t0 values. let's explicitly make the distinction between
		# these two parameters and see where it gets us...
		if m==None: m=self.mag
		if alpha==None: alpha=self.alpha
		#
		return self.lDeltat(m=m, alpha=alpha) + 1.0
	###
	##
	# method 2: upper rate limit is where events become resolvable (dt_omori > dt_rupt.)
	# ultimately, i think these two approaches both work (except it is better to estimate one
	# parameter from intrinsic/intensive physics (max rate, etc.) and then constrain the second
	# by integrating Omori.
	#
	# ... and of course, now all of this needs to be cleaned up a bit.
	def getlt02 (self, mstar=None, alpha=alpha_0, dmstar=None, drho=None):
		# the value drho=5.25 appears to be good within about +/- 0.5 for medium earthquakes.
		if mstar==None: mstar = self.mag
		if alpha==None: alpha=self.alpha
		if dmstar==None: dmstar=self.dmstar
		if drho==None: drho=self.drho
		#
		x=(.5*mstar) - alpha - dmstar - 1.0 + drho
		return x

	def getltau2(self, mstar=None, mc=None, alpha=alpha_0, p=None, drho=None, thislt0=None):
		if mstar==None: mstar = self.mag
		if alpha==None: alpha=self.alpha
		if drho==None: drho=self.drho
		if p==None: p=self.p
		if thislt0==None:
			thislt0=getlt02(mstar=mstar, alpha=alpha, dmstar=1.2, drho=drho)
		x = mc-.5*mstar-alpha - p*thislt0 + drho
		return x
	
	def lruptLen(self, m, dl0=1.0):
		return (m/2.0) - 1.755 - dl0
	
	def localIntensity(self, inloc, intime, rho=None):
		#
		if rho==None: rho=self.rho	# radial density scaling exponent
		#
		# this is the real thing. return estimated current rate of seismicity based on
		# Omori and spatial decay.
		#
		# specifically, report dN/(dt * dA), or d(omori rate)/dA.
		#
		# this is a little bit tricky with the modified Omoir type formats.
		# the easiest way is probably to get the first (Omori) intensity directly.
		# then, calculate the relative spatial intensity by comparing N'(R) to N'(R->l_r).
		# otherwise, we have to solve for the spatial (or Omori) parameters based on the 
		# decayed rate (aka, new t0, tau parameters).
		# also, assume the decay rate of the full (compound) sequence is p=1.0 (not 1.1),
		# and something similar for spatial decay (though the 1.37 numbers from Frolich probably
		# represent the full sequence.
		# do we integrate over the area here, or is that done by the calling function? (which
		# probably assumes linearity across the bin).
		t=self.timesince(intime)
		if t<0.0: return 0.0	# earthquake has not yet happened.
		#
		#if t<10.0**self.lt0:
		#
		# t is the actual time since the (begining of?) the earthquake. convert to t', the "scaling time".
		# 
		tcrit = 10.0**(1.0+self.lrupttime)	# but i think i just made up this number. let's try something more
														# scientific. how'bout dt_omori(t=0) = tau * t0**p
														# that said, this does not seem to be terribly important, except as 
														# an mechanism to handle the singularity.
		#tcrit = 10**(self.drho+self.mc-self.alpha-self.mag/2.0)
		#print "tcrit: %f, %f" % (tcrit/days2secs, (10.0**(1.0+self.lrupttime))/days2secs)
		if t<tcrit:
			#t=10.0**self.lt0	# flat bit of the Omori distribution.
			t = 0.0
		else:
			t=t-tcrit
		#
		# print "t: %f" % t
		#
		lrupture=10**self.lruptlen
		R=self.dispFrom(inloc=inloc, coordtype='rad')	# displacement...
		r=R[0]
		#if r<10.0**(self.lruptlen):
		#	r=10.0**(self.lruptlen)
		if r<lrupture:
			r=lrupture		# this implies that the spatial-rate distribution saturates at (some factor of) lrupture
		#
		# the basic idea is to
		#1) calculate the absolute Omori or spatial rate/density
		#2) then, calculate the expected decay in the other.
		# aka, total omori rate at t=t', then the relative decay in the spatial dimension
		# 
	
		#orate=self.omoriRate(t=t, t0=self.gett0(), tau=self.gettau(), p=1.0)
		orate=self.omoriRate(t=t, t0=10**self.lt0, tau=10**self.ltau, p=self.p)
		#
		# note: the omori rate is being estimated directly from the magnitude
		# (by estimating the parameters tau, t0).
		#Nremaining=(10**(self.mag - 1.0 - self.mc)) - self.omoriN(t=t, t0=self.t0, tau=self.tau, p=self.p)
		#Nomori = 10.0**(self.mag - 1.0 - self.mc)
		#Nremaining = Nomori - self.omoriNln(t=t, t0=10.0**self.lt0, tau=self.tau)
		# and we use that as N in the spatial density formula.
		# the relative intensity factor will be
		#intensityFactor=Nremaining/Nomori
		#lrupture=10**self.lruptlen
		#spatialDecayFactor = (1.0 + r/lrupture)**(1.0-self.rho)
		#spatialDecayFactor = (1.0 + r/lrupture)**(-self.rho)		# this is the relative linear density (decay) at r. "1+" comes from using the CDF version.
																					# but it is probably bettter to just calc. the density at r, then integrate around the area.
																					# in fact, it is arguably the job of forecastlocation() to do that. earthquake just reports
																					# the local intensity at r
		#start with linear density (Felzer):
		# lambda = N/l_r = c * 1/(r0+r)**rho
		# integrate and set equal to Omori:
		# N_omori = integral [0 to inf] [ lambda dr ] = (c/(rho-1)) * r0**(1-rho)
		# solve for c:
		# c = N_omori * (rho-1) * r0**(rho-1)
		# r0 = l_rupt.
		# now, substitute omori-rate for N_omori to give a linear rate. to find the local rate at some radial distance
		# and some radial location theta, divide the total rate by 2pi*r
		#
		#spatialdensity = (orate/(2.0*math.pi*(lrupture**(-self.rho)))) * r**(-(1.0 + self.rho))
		#
		#spatialdensity = ((orate * (rho-1.0)*lrupture**(rho-1.0)))/( 2.0*math.pi*r**(1.0+rho))
		radialDens = self.radialDensity(r=r, lrupt=None, rho=None, dm=None)	# all events at distance r
		spatialdensity = radialDens/(2.0*math.pi*r)
		# noting that the 1+rho exponent comes from 1/r**rho * 1/2pi*r
		# also note that we are using r instead of (r0+r'). we correct for this 
		#
		#localIntensity = orate*spatialDecayFactor
		localIntensity = spatialdensity*orate
		
		#
		# so this is the total rate of aftershock production.
		# now, calculate the spatial distribuiton. the total number of remaining (expected) aftershocks is:
		#
		#lruptTime = self.lDeltat(m, alpha)
		
		#return 1.0+lruptTime			
		
		return localIntensity


def plotcolor(z, Z):
	# z is value, Z is max value
	#
	colorRange=int('0xFFFFFF', 0)
	#
	colorint=int(colorRange*z/Z)
	#colorhex=self.fillHexString(hex(colorint), 6)
	colorhex=fillHexString(hex(colorint), 6)
	colorhexstr='#' + colorhex.split('x')[1]
	#
	return colorhexstr

def fillHexString(hexstring='', strlen=6):
	strs=hexstring.split('x')
	while len(strs[1])<strlen:
		strs[1]='0'+strs[1]
	return strs[0] + 'x' + strs[1]

def getMFETAScatFromANSS(lons=[-121.0, -114.0], lats=[31.0, 37.0], dates=[None, None], mc=4.0):
	if dates==None: dates=[None, None]
	if dates[1]==None: dates[1]=dtm.datetime.now(pytz.timezone('UTC'))
	if dates[0]==None or dates[0]>=dates[1]: dates[0]=dates[1]-dtm.timedelta(days=840)
	#
	print dates
	clist1=atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=dates, Nmax=999999, fout=None)
	catalog=[]
	#X,Y,M = [], [], []
	#
	for rw in clist1:
		#catalog+=[[mpd.datestr2num(rw[0]), rw[1], rw[2], rw[4]]]
		catalog+=[[mpd.datestr2num(rw[0]), rw[1], rw[2], rw[3], rw[4]]]
		#X+=[rw[2]]	# lon
		#Y+=[rw[1]]	# lat
		#M+=[rw[4]]
	return catalog

def checkPoly(inpoly, fignum=None):
	# polygons are getting mangled by the "find closest neighbor" algarithm.
	# try looking for dangling elements after a polygon has closed. we'll then either
	# discard them or insert them next to their closest neighbor.
	#
	# given local syntax, inpoly -> [[x], [y]]
	#
	x0=inpoly[0][0]
	y0=inpoly[1][0]
	outpoly=[[[inpoly[0][0]], [inpoly[1][0] ]]]
	#outpoly=[[[], []]]
	for i in xrange(1,len(inpoly[0])):
		thisx=inpoly[0][i]
		thisy=inpoly[1][i]
		#
		outpoly[-1][0]+=[thisx]
		outpoly[-1][1]+=[thisy]
		#if thisx==x0 and thisy==y0 and len(outpoly[-1][0])>1:
		if thisx==outpoly[-1][0][0] and thisy==outpoly[-1][1][0] and len(outpoly[-1][0])>1:
			print "(improperly) closed poly, %d/%d" % (i, len(inpoly[0])-1)
			outpoly+=[[[], []]]
	
	if len(outpoly[-1][0])<2: outpoly.pop()
	if outpoly[-1][0][0]!=outpoly[-1][0][-1] and outpoly[-1][1][0]!=outpoly[-1][1][-1]: outpoly.pop()
	#
	#plt.plot(inpoly[0], inpoly[1], 'o-', lw=3, alpha=.4)
	if fignum!=None:
		plt.figure(fignum)
		plt.clf()
	
		for ply in outpoly:
			plt.plot(ply[0], ply[1], '.--', alpha=.8)

	return outpoly

