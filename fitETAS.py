# fitETAS.py
# Mark R. Yoder, Ph.D.
# use these scripts to generate figures for PAGEOPH paper:
# Near-field ETAS constraints and applications to seismic hazard assessment
#
import yodapy as ypp
import BASScast as bcp
import ANSStools as atp
import eqcatalog as eqp
import linefit as lft
import yodabass as ybp
import ETASscripts as esp

import datetime as dtm
import pytz
import operator
import math
import random
import numpy
import scipy
import scipy.optimize as spo
import os
import multiprocessing as mcp
#
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.cm as cmx
#
import matplotlib.dates as mpd
import matplotlib.mpl as mpl
import scipy.optimize as spo
#
from geographiclib.geodesic import Geodesic as ggp
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

tzUTC=pytz.timezone('UTC')
days2secs=3600.0*24.0
plt.ion()

def fitdoodle(nits=10000, arange=[-1.5, .5], brange=[-1.0, 1.0], crange=[-10.0, 10.0]):
	nits=int(nits)
	# a quick fit doodle...
	#alphas = scipy.ones(4)
	#betas = scipy.ones(4)
	#
	equakes=['parkfield', 'sansim', 'coalinga', 'elmayor30', 'elmayor25', 'hmine']
	ms=  scipy.array([6.0, 6.5,   6.7,  7.2,  7.2,  7.3])
	mcs=  scipy.array([1.5, 3.5,   3.0,  3.0,  2.5,  3.0 ])
	#ndots=scipy.array([3.07, 2.36, 2.37, 1.8,  1.76, 2.08])
	ndots=scipy.array([2.19, 2.36, 2.37, 1.8,  1.76, 2.08])
	#ndots=scipy.array([1.63, 2.36, 2.37, 1.8,  1.76, 2.08])
	sigs= scipy.array([.60, .43,    .46,  3.8,  .36, .52])
	Ws=1.0/sigs**2.0
	#
	# these data, from Enescu2009, are great, but the fitting method seems to be in conflict with
	# Shcherbakof 2006,etc. They get much smaller t0 values, on par with rupture duration.
	msjapan =  [6.6, 6.6, 6.7, 6.6]
	mcsjapan = [3.1, 3.3, 3.3, 3.3]
	ld2s=math.log10(24.*3600.)
	ndotsjap = [-4.5+ld2s, -3.9+ld2s, +4.0-ld2s, +4.0-ld2s]	
	#equakes += ['Tottori', 'Fukuoka', 'Noto', 'Chuetsu-oki']
	
	#ms = scipy.array(ms.tolist() + msjapan)
	#mcs = scipy.array(mcs.tolist() + mcsjapan)
	#ndots = scipy.array(ndots.tolist() + ndotsjap)
	#Ws=scipy.ones(len(ndots))
	#
	'''
	ms=scipy.array(msjapan)
	mcs=scipy.array(mcsjapan)
	ndots = scipy.array(ndotsjap)
	equakes=['Tottori', 'Fukuoka', 'Noto', 'Chuetsu-oki']
	sigs=scipy.ones(len(equakes))
	Ws=sigs
	'''
	#
	err=270000.0
	amin=0.
	bmin=0.
	cmin=0.
	#
	errw=270000.0
	aminw=0.
	bminw=0.
	cminw=0.
	#
	alpha=1.0
	beta=1.0
	Ra=random.Random()
	Rb=random.Random()
	Rc=random.Random()
	#
	print "arange: [%f, %f], brange: [%f, %f], crange: [%f, %f]" % (arange[0], arange[1], brange[0], brange[1], crange[0], crange[1])
	#
	# converging fit:
	# plsq=spo.leastsq(fitres, p, args=(scipy.array(Y), scipy.array(X), scipy.array(W)), full_output=1)
	p0=scipy.array([arange[0] + .5*(arange[1]-arange[0]), brange[0] + .5*(brange[1]-brange[0]), crange[0] + .5*(crange[1]-crange[0])])
	#print "p0: ", p0
	#
	plsq = spo.leastsq(doodleerror2, p0, args=(scipy.array(ms), scipy.array(mcs), scipy.array(ndots), Ws), full_output=1 )
	#
	# now, calc the mean-error, etc.
	print "lsq fit: " , plsq[0]
	errplsq=doodleerror2(plsq[0], ms, mcs, ndots)
	#print "error array: ", errplsq
	dlen=float(len(ms))
	print "mean errors: ", errplsq.mean(), sum(errplsq)/(dlen-3.), errplsq.std()
	#
	# omoriError(prams, t, dts, m, mc, p0=1.2)
	#
	for i in xrange(nits):
		# so we guess values for a*m + b*mc = ndot (logs anyway) and minimize the total error.
		thisa = arange[0] + (arange[1]-arange[0])*Ra.random()
		thisb = brange[0] + (brange[1]-brange[0])*Rb.random()
		thisc = crange[0] + (crange[1]-crange[0])*Rc.random()
		#
		#thisa=-1.5 + 2.0*Ra.random()
		#thisb= -1.0 + 2.0*Rb.random()
		#A=10.
		#thisc=A*Rc.random()-A/2.
		#thisc=-10.0 + 20.0**Rc.random()
		#
		thiserr=0.0
		thiserrw=0.0
		for j in xrange(len(ms)):
			derr=doodleerror(ms[j], mcs[j], ndots[j], thisa, thisb, thisc)
			thiserr+=derr
			thiserrw+=derr*Ws[j]
		if thiserr<err:
			#print "new UNweighted best fit: %.4f, %.4f, %.4f, %.4f" % (thisa, thisb, thisc, thiserr)
			err=thiserr
			amin=thisa
			bmin=thisb
			cmin=thisc
		if thiserrw<errw:
			#print "new WEIGHTED best fit: %.4f, %.4f, %.4f, %.4f" %(thisa, thisb, thisc, thiserr)
			errw=thiserrw
			aminw=thisa
			bminw=thisb
			cminw=thisc
		#
	#
	print "finished."
	print "new UNweighted best fit: %.4f, %.4f, %.4f, %.4f" % (amin, bmin, cmin, err)
	print "new WEIGHTED best fit: %.4f, %.4f, %.4f, %.4f" %(aminw, bminw, cminw, errw)
	for j in xrange(len(ms)):
		m=ms[j]
		mc=mcs[j]
		ndot=ndots[j]
		#print "(%s): m=%f, mc=%f, ndot=%f :: %f/%f" % (equakes[j], m, mc, ndot, (amin*m + bmin*mc + cmin), (amin*m + bmin*mc + cmin)-ndots[j])
		print "(%s): m=%f, mc=%f, ndot=%f :: %f/%f" % (equakes[j], m, mc, ndot, (plsq[0][0]*m + plsq[0][1]*mc + plsq[0][2]), (plsq[0][0]*m + plsq[0][1]*mc + plsq[0][2])-ndots[j])


def getETAScat(catname='parkfield', dorefresh=True, mc=None):
	if catname=='coalinga':
		filename='cats/coalinga-mfetas.cat'
		if mc==None: mc=3.5
		lats=[36.22-.7, 36.22+.7]
		lons=[-120.32-.7, -120.32+.7]
		daterange=[dtm.datetime(1983, 5, 2, 23, 42, 38, 60000, tzUTC), dtm.datetime(2003,12,31, 0, 0, 0, 0, tzinfo=tzUTC)]
		mainshock=[dtm.datetime(1983, 5, 2, 23, 42, 38, 60000, tzinfo=tzUTC), 36.2317, -120.312, 6.7, 10.18]
	#
	if catname=='sansim':
		filename='cats/sansim-mfetas.cat'
		if mc==None: mc=3.0
		lats=[35.0, 36.4]
		lons=[-121.8, -120.4]
		daterange=[dtm.datetime(2003, 12, 22, 19, 15, 56, 240000, tzinfo=pytz.timezone('UTC')), dtm.datetime(2012,12,31, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC'))]
		mainshock=[dtm.datetime(2003, 12, 22, 19, 15, 56, 240000, tzinfo=tzUTC), 35.7005, -121.1005, 6.5, 8.7]
	#	
	if catname=='parkfield':
		filename='cats/parkfield-mfetas.cat'
		#filename='cats/parkfield.cat'
		if mc==None: mc=1.5
		lats=[35.45, 36.4]
		lons=[-121.0, -120.0]
		daterange=[dtm.datetime(2004, 9, 28, 17, 15, 24, 250000, tzinfo=tzUTC), dtm.datetime(2009,12,31, 0, 0, 0, 0, tzinfo=tzUTC)]
		mainshock=[dtm.datetime(2004, 9, 28, 17, 15, 24, 250000, tzinfo=tzUTC), 35.8182, -120.366, 5.97, 8.58]
	#
	if catname=='hmine':
		filename='cats/hmine-mfetas.cat'
		if mc==None: mc=3.0
		lons=[-117.25, -115.5]
		lats=[33.5, 35.5]
		daterange=[dtm.datetime(1999, 10, 16, 9, 46, 44, 130000, tzinfo=tzUTC), dtm.datetime(2008,12,31, 0, 0, 0, 0, tzinfo=tzUTC)]
		mainshock=[dtm.datetime(1999, 10, 16, 9, 46, 44, 130000, tzinfo=tzUTC), 34.594, -116.271, 7.1, 0.02]
		# note landers: [datetime.datetime(1992, 6, 28, 11, 57, 34, 130000, tzinfo=<UTC>), 34.2, -116.437, 7.3, 0.97]
	#				
	if catname=='tohoku':
		filename='cats/tohoku-mfetas.cat'
		if mc==None: mc=4.5
		lons=[135.0, 150.0]
		lats=[30.0, 41.5]
		daterange=[dtm.datetime(2011, 3, 11, 5, 46, 24, 120000, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(tzUTC)]
		mainshock=[dtm.datetime(2011, 3, 11, 5, 46, 24, 120000, tzinfo=tzUTC), 38.297, 142.373, 9.1, 29.0]
	#
	if catname=='elmayor':
		filename='cats/elmayor-mfetas.cat'
		if mc==None: mc=2.5
		lons=[-118.0, -113.5]
		lats=[31.0, 35.0]
		daterange=[dtm.datetime(2010, 4, 4, 22, 40, 42, 360000, tzinfo=pytz.timezone('UTC')), dtm.datetime(2012,12,31, tzinfo=pytz.timezone('UTC'))]
		mainshock=[dtm.datetime(2010, 4, 4, 22, 40, 42, 360000, tzinfo=tzUTC), 32.2862, -115.2953, 7.2, 10.0]
	#
	# now, get the catalog:
	if dorefresh==True:
		hmcat = bcp.atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=daterange, Nmax=999999, fout=filename)
		cpf=ypp.eqcatalog(hmcat)
	else:
		cpf=ypp.eqcatalog()
		cpf.loadCatFromFile(filename, minmag=mc)
	cpf.mc=mc
	cpf.mainshock=mainshock
	
	#
	# special cases:
	if catname=='parkfield':
		cpf.addEllipCat('PFshock (.4 x .15)', cpf.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
	#
	return cpf
	
def getEqInts(objCat, mc=None, mainshock=None, catnum=0, fnum=None):
	eqints=objCat.getIntervals(catList=objCat.getcat(catnum), interval_length=1).tolist()	# note: numpy arrays have weird sorting.
	eqints.sort(key = lambda x: x[0])
	while eqints[0][0]<=(objCat.mainshock[0]): eqints.pop(0)
	dt0=mpd.date2num(objCat.mainshock[0])
	#fT=(scipy.array(map(mpd.date2num, map(operator.itemgetter(0), eqints))-dt0)*days2secs
	#
	T=(scipy.array(map(mpd.date2num, map(operator.itemgetter(0), eqints)))-mpd.date2num(objCat.mainshock[0]))*days2secs
	dt=scipy.array(map(operator.itemgetter(1), eqints))*days2secs
	#
	if fnum!=None:
		f1=plt.figure(fnum)
		plt.clf()
		plt.plot(map(operator.itemgetter(0), eqints), map(operator.itemgetter(1), eqints), '.-')
	
		meandt=scipy.array(map(operator.itemgetter(1), eqints[0:10])).mean()
		print meandt
	'''
	f1=plt.figure(0)
	ax=plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.plot(T,dt, '.-.')
	'''
	#	
	#return eqints
	return [T,dt]

def getNN(eqset='parkfield', catnum=0, mc=None, fnum=0, directional=False):
	# directional=True: the second element must be farther (r_i > r_(i-1)) away.
	# nominally, this would be a great opportunity to use this as an index, but we
	# only run this now and again...
	cpf=getETAScat(catname=eqset, mc=mc)
	if catnum>len(cpf.subcats): catnum=len(cpf.subcats)
	#
	# assume catalogs are mc selected...
	thismc=min(map(operator.itemgetter(3), cpf.getcat(catnum)))
	#
	mev = cpf.getMainEvent(cpf.getcat(catnum))
	phis0 = mev[2]					
	lambs0 = mev[1]				#phis0, lambs0 still in degrees...
	#
	#mstar = cpf.getMainEvent()[3]
	mstar=mev[3]
	while cpf.getcat(catnum)[0][0]<mev[0]: cpf.getcat(catnum).pop(0)
	bathlen=10**(mstar-1.0-cpf.mc)
	while len(cpf.getcat(catnum))>bathlen: cpf.getcat(catnum).pop(-1)
	#
	# getDistances(self, catList=None, center=None):
	# will return like: dists += [[rw[0], thisx, thisy, R0, R1, R2, R3, R4, rw[4]]]
	#D0=cpf.getDistances(catList=cpf.getcat(catnum), center=None)	# we'll need r_NN(r) -- aka, \delta_r(r)
	#
	# now, get the distance of each event from each event.
	NNs=[]	# [ r_mev, \delta_r ]
	distkey=7
	for ev in cpf.getcat(catnum):
		#print ev
		# get distance from mainshock:
		g1=ggp.WGS84.Inverse(lambs0, phis0, ev[1], ev[2])
		r=g1['s12']/1000.0	# and converting from m to km.
		#if r<=0.0: continue
		#
		#print "center: ", [ev[2], ev[1]]
		dists=[]
		for i in xrange(0, len(cpf.getcat(catnum))):
		#for ev2 in cpf.gatcat(catnum):
			ev2=cpf.getcat(catnum)[i]
			if ev2==ev: continue
			if directional==True:
				g20=ggp.WGS84.Inverse(lambs0, phis0, ev2[1], ev2[2])	# this event's distance from center.
				r20=g20['s12']/1000.0	# and converting from m to km.
				if r20<r: continue
			
			g2=ggp.WGS84.Inverse(ev2[1], ev2[2], ev[1], ev[2])
			dr=g2['s12']/1000.0
			#
			if len(dists)==0 or dr<dists[-1]: dists+=[dr]
		#	
		#thisdists=cpf.getDistances(center=[ev[2], ev[1]])
		# let's use the geographiclib distance (R4 = rw[7]). the best spherical solution is R3, rw[6]
		#thisdists.sort(key=lambda x: x[distkey])
		#NNs+=[r, thisdists[-1][distkey]]
		#
		if len(dists)>0: NNs+=[[r, dists[-1]]]
		#
	foutname = 'data/%snns.dat' % eqset
	if directional==True: foutname = 'data/%snnsDir.dat' % eqset
	fout=open(foutname, 'w')
	fout.write('#%s\n' % eqset)
	fout.write('#r\tdr\n')
	for rw in NNs:
		fout.write('%f\t%f\n' % (rw[0], rw[1]))	
	fout.close()
	#	
	return NNs

def plotNNfile(nnfile, mc=None, mstar=None):
	nns=[]
	f=open(nnfile)
	for rw in f:
		if rw[0]=='#': continue
		nns+=[map(float, rw.split())]
	#
	return plotNNs(nns, mc, mstar)

def plotNNs(nns, mc=None, mstar=None):
	X=map(operator.itemgetter(0), nns)
	Y=map(operator.itemgetter(1), nns)
	Yns=scipy.array(Y).copy()
	Yns.sort()
	Ns=range(1, len(X)+1)
	Ns.reverse()
	#
	dens=[]
	xdens=[]
	for i in xrange(len(X)):
		if Y[i]==0.0: continue
		dens+=[1.0/Y[i]]
		xdens+=[X[i]]
	#
	plt.figure(0)
	plt.clf()
	plt.plot(X,Y, '.')
	if mc!=None:
		dx= max(xdens) - min(xdens)
		plt.plot([min(xdens) + dx*.01, max(xdens)-dx*.1], [10**(mc/2.0 - 1.76), 10**(mc/2.0 - 1.76)], 'd--', lw=2, label='$L_{c}$')
	#	
	plt.title('NN distances of aftershocks')
	plt.xlabel('distance from mainshock $r$')
	plt.ylabel('NN distance $\\Delta r$')
	#
	plt.figure(1)
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(X,Y, '.')
	if mc!=None:
		dx= max(xdens) - min(xdens)
		y0=10**(mc/2.0 - 1.76)
		plt.plot([min(xdens) + dx*.01, max(xdens)-dx*.1], [y0, y0], 'd-', lw=2, label='$L_{c}$')
		if mstar!=None:
			# tohoku Lr correction: (378./670.)
			Lr=10**(mstar/2.0 - 1.76)
			Lc=10**(mc/2.0 - 1.76)
			#plt.plot([Lr, Lr], [10.*min(Y), .9*max(Y)], 's--', lw=2)
			plt.plot([Lr, Lr], [.01, 100], 's--', lw=2)
			plt.plot([Lr/10., Lr/10.], [.01, 100], 's--', lw=2)
			plt.plot([Lr, Lr], [.01, 100], 's--', lw=2)
			#plt.plot([Lc, Lc], [.01, 100], 's--', lw=2)
			newX=Lr+(scipy.arange(500.0)/500.)*(max(X)-Lr)
			newY=map(pow, (newX/Lr), scipy.ones(len(newX))*1.35)
			newY=y0*scipy.array(newY)
			plt.plot(newX, newY, 'c--', lw=2)
			#
			# sloped parts?
			newXa=Lr/10.+(scipy.arange(500.0)/500.)*(max(X)-Lr/10.)
			newYa=map(pow, (newXa/(Lr/10.)), scipy.ones(len(newX))*1.35)
			newYa=y0*scipy.array(newYa)
			plt.plot(newXa, newYa, 'r--', lw=2)
			#
			'''
			newXb=Lc+(scipy.arange(500.0)/500.)*(max(X)-Lc)
			newYb=map(pow, (newXb/Lc), scipy.ones(len(newX))*1.35)
			newYb=y0*scipy.array(newYb)			
			plt.plot(newXb, newYb, 'g--', lw=2)
			'''
			#
			# let's get a fit for the r<Lr elements:
			fXY=[]
			for i in xrange(len(X)):
				if X[i]==0. or Y[i]==0.: continue
				if X[i]<=(Lr/2.): fXY+=[[math.log10(X[i]), math.log10(Y[i])]]
			fXY.sort(key = lambda x: x[0])
			f1=lft.linefit(fXY)
			plsq=f1.doFit()[0]
			yfit=plsq[0]
			print "r<Lr fit: " , plsq, numpy.average(map(operator.itemgetter(1), fXY)), numpy.average(map(operator.itemgetter(1), fXY))/.85
			mcalpha=((1.76+yfit)/mc)
			print "estimated mc=%f, or estimated mc-fact: %f/%f/%f" % ((2.*(yfit + 1.76)), mcalpha, mcalpha/.9, mcalpha/.8)
			plt.plot([min(xdens) + dx*.01, max(xdens)-dx*.1], [10**yfit, 10**yfit], 'd--', lw=2, label='$fit: L_{c}$, $\\alpha_{mc}=%.2f$' % (mcalpha))
	#
	plt.title('NN distances of aftershocks')
	plt.xlabel('distance from mainshock $r$')
	plt.ylabel('NN distance $\\Delta r$')
	plt.legend(loc='best', numpoints=2)
	#
	plt.figure(2)
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(xdens, dens, '.')
	if mc!=None:
		dx= max(xdens) - min(xdens)
		plt.plot([min(xdens) + dx*.01, max(xdens)-dx*.1], [10**-(mc/2.0 - 1.76), 10**-(mc/2.0 - 1.76)], 'd--', lw=2, label='$L_{c}$')
	plt.title('NN aftershock densities')
	plt.xlabel('distance from mainshock $r$')
	plt.ylabel('NN density, $1/\\Delta r$')
	#
	plt.figure(3)
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(Yns, Ns, '.')
	plt.xlabel('NN distances $\\Delta r$')
	plt.ylabel('Cumulative number $N$')
	plt.title('NN distances CDF: $N(\\Delta r)$')
	#
	return nns

def maxDensityExp(mstar, mc, dm=1.0, D=1.8):
	# i guess not "density" per se, but actually N_as, the number of aftershocks
	# crammed into a rupture.
	return (2.0/(2.0+D))*math.log10(1.0 + D/2.0) + (D/(2.0 + D))*(mstar - dm - mc)

def q(mstar, mc, D=1.5, dm=1.0, dlambd=1.76, b=1.0):
	Nom=10.0**(b*(mstar-dm-mc))
	N0p = 10.0**maxDensityExp(mstar=mstar, mc=mc, dm=dm, D=D)
	r0 = 10**(mstar/2.0 - dlambd)
	rc = 10**(mc/2.0 - dlambd)
	#
	q=Nom/(Nom-(N0p*r0*.5))
	r0prime=(r0/rc)**(2.0/(2.0+D))
	q=Nom/(Nom-(N0p*r0prime*.5))
	#
	return q
	

def mcDsolver(N, mstar, mc, dm=1.0, Drange=[0.0, 2.0], nits=10**6):
	##
	# returns the estimated dimension D.
	# note: N is really log(N)... 
	#
	drange=Drange[1]-Drange[0]
	R1=random.Random()
	minerr=N*10.0**10.0
	for i in xrange(nits):
		thisD = Drange[0] + drange*R1.random()
		err=abs(N-maxDensityExp(mstar, mc, dm, D=thisD))	# a little faster than x**2.0 i think.
		#
		if err<minerr:
			minD = thisD
			minerr = err
		#
	#
	return [minD, minerr]

	
def fitEqDists(eqset='parkfield', catnum=None, q=1.35, mc=None, dalpha=2.0, nbins=100, fnum=0, rlenfac=.5, dmrl=1.0, D=1.5, hidelegend=True):
	if catnum==None:
		catnum=0
		if eqset=='parkfield': catnum=1
	#
	cpf=getETAScat(catname=eqset, mc=mc)
	b=1.0
	if catnum>len(cpf.subcats): catnum=len(cpf.subcats)
	#
	# assume catalogs are mc selected...
	thismc=min(map(operator.itemgetter(3), cpf.getcat(catnum)))
	#
	mstar = cpf.getMainEvent(cpf.getcat(catnum))[3]
	#if eqset=='coalinga': mstar=6.3
	
	print "eqset: %s, mstar: %f" % (eqset, mstar)
	#
	while cpf.getcat(catnum)[0][0]<cpf.getMainEvent()[0]: cpf.getcat(catnum).pop(0)
	Nom     = 10**(b*(mstar-1.0-cpf.mc))
	bathlen = Nom
	while len(cpf.getcat(catnum))>bathlen: cpf.getcat(catnum).pop(-1)
	
	
	lrlen=mstar/2. - 2.76			# the "rupture core" length
	fullruptlen=10**(mstar/2. - 1.76)
	ruplenFull = fullruptlen		# for consistency with ruplenGraph
	#ruplen=rlenfac*10.0**(mstar/2. - 1.76)
	ruplenOmori=rlenfac*10.0**(mstar/2. - 1.76)	# rupture length for Omori-like. basically, for center-ruptures,
																# rlenfac=.5. For functions with circular symm., this is probably
																# the right value to use. for functions that are specifically about
																# cartesian coords/area, productivity, etc.... well, it depends.
																# note that the omori-like functs. are parameterized as dN/dr.
	#
	print "mstar=%f, mc=%f, lrlen=%f, rlen=%f/%f" % (mstar, thismc, lrlen, 10**lrlen, 10**(lrlen+1.))
	#
	dists=cpf.getDistances(cpf.getcat(catnum))
	#plt.title('distances for %s' % eqset)
	# R0: rect-lin distance, R1: approx1, R2: approx2 (better than 1), R3: bomb-proof spherical, R4: from geographicLib,
	#dists += [[rw[0], thisx, thisy, R0, R1, R2, R3, R4, rw[4]]]	# and R3 is the only distance we really need; rw[4]->depth.
	#
	r0s=map(operator.itemgetter(3), dists)
	r3s=map(operator.itemgetter(6), dists)
	r4s=map(operator.itemgetter(7), dists)
	#
	# 3d-dist:
	X,Y,Z = scipy.array(map(operator.itemgetter(1), dists)), scipy.array(map(operator.itemgetter(2), dists)), scipy.array(map(operator.itemgetter(8), dists))
	r3d=map(math.sqrt, (X**2. + Y**2. + Z**2.))
	#
	Ninsides=[0.,0.,0.,0.]
	for j in xrange(len(r0s)):
		if r0s[j]<(ruplenOmori): Ninsides[0]+=.5*(1.0/rlenfac)*1./ruplenOmori		# note: N/{r from epicenter}
		if r3s[j]<(ruplenOmori): Ninsides[1]+=.5*(1.0/rlenfac)*1./ruplenOmori		# so, if rlenfac=.5, each dr -> 2*dr.
		if r4s[j]<(ruplenOmori): Ninsides[2]+=.5*(1.0/rlenfac)*1./ruplenOmori
		if r3d[j]<(ruplenOmori): Ninsides[3]+=.5*(1.0/rlenfac)*1./ruplenOmori
	#
	#print b, mstar, mc, ruplen
	#chi = gettau(N=10.0**(b*(mstar-1.0-thismc)), t0=ruplen, p=q)
	N=10**(b*(mstar-1.0-thismc))
	chi = ruplenOmori**(1.0-q)/(N*(q-1.0))
	#
	ruplen2=rlenfac*10.0**((mstar-1.0)/2. - 1.76)
	#ruplen2=1.0*10.0**((mstar-1.0)/2. - 1.76)
	#chi2=ruplen2**(1.0-q)/(N*(q-1.0))
	#
	ruplencore=rlenfac*10**(mstar/2. -2.76)
	#ruplencore=1.0*10**(mstar/2. - 2.76)
	chicore=ruplencore**(1.0-q)/(N*(q-1.0))
	#
	# t, tau, t0, p
	nelements=len(r4s)
	# omorirate(t, tau, t0, p)
	Yth=[]  #map(omorirate, r4s, scipy.one(nelements)*chi, 10**(ruplen), q)
	Yth2=[]
	Ythcore=[]
	YthPL=[]
	YthSSim=[]
	#rc=10.0**(thismc/2. - 1.76)
	#xth=scipy.array(r4s).copy()
	# : xth=min(X) + (scipy.arange(float(len(X)+1))/float(len(X)))*(10.-1.)
	xthrange=(max(r4s)-min(r4s))*2.0
	npoints=1000
	xth=min(r4s) + (scipy.arange(float(npoints+1))/float(npoints))*xthrange
	xth.sort()
	#
	# (estimated rate exponent)
	#estRateExp=(mstar-dmrl-thismc + math.log10(2.))*.5		# this is the right idea, but not precise with D
	#																			# and there's a constant part that gets multiplied...
	estRateExp = maxDensityExp(mstar=mstar, mc=thismc, dm=1.0, D=D)
	mt=7.6
	if mstar>mt:
		#meandm = ((1.5*(mstar-mt) + 1.0*(mt-thismc))/(mstar-thismc) + (1.0))/2.0
		meandm = ((1.5*(mstar-mt) + 1.0*(mt-thismc))/(mstar-thismc)*(mstar-mt) + (1.0)*(mt-thismc))/(mstar-thismc)
		print "m>7.6 meandm: %f" % meandm
		#
		# following the same area integration methd, this overestimates:
		#estRateExp = (.75*mstar - .5*thismc - .25*mt - .75*meandm + math.log10(2.))
		#
		# this seems to be about right (just use m_eff -> m_t):
		#estRateExp=(mt-dmrl-thismc + math.log10(2.))*.5
		#
		estRateExp = maxDensityExp(mstar=mt, mc=thismc, dm=1.0, D=D)
		#
		#for k in xrange(len(Ninsides)):
		#	Ninsides[k]*=(10**((mt-mstar)/2.))
	#
	rate0 = 10**(1.0*estRateExp)
	r0ssim = Nom*(q-1.0)/rate0		# r0 in the ssim config.
	chi2 = 1.0/(rate0*ruplen2**q)	# scaling initiates at r_thresh ~ L_R, but then decays according to L_{r-m1}
	for x in xth:
		# this is supposed to be callibrated for rupture-centered earthquakes, so for rupture-tip earthquakes,
		# the density/intensity should be about 1/2 that of a centered earthquake.
		# this should probably be, then, .5*(1.0/rlenfac)*...
		Yth+=[.5*(1.0/rlenfac)*omorirate(x, chi, ruplenOmori, q)]
		#Ythcore+=[.5*(1.0/rlenfac)*omorirate(x, chicore, ruplencore, q)]
		# for now, hijack Ythcore for Ythssim:
		Ythcore+=[.5*(1.0/rlenfac)*omorirate(x, chi, r0ssim, q)]
		#
		#if x<=(ruplen*rlenfac):
		if x<=ruplenOmori:  
			#Yth2+=[2.*omorirate(0.0, chi, ruplen, q)]
			#Yth2+=[(1.0/rlenfac)*omorirate((0.), chi2, ruplen2, q)]
			Yth2 += [.5*(1.0/rlenfac)*rate0]
			YthPL += [.5*(1.0/rlenfac)*rate0]		# check this fore-factor. should it be .5*(1/rlenfac)?
		#
		#if x>(ruplen*rlenfac):
		if x>ruplenOmori:
			q2=q
			#Yth2+=[10.*omorirate(ruplen, chi, ruplenOmori, q)*(1.0*(rc+x-ruplen)**(-q2))]
			Yth2+=[.5*(1.0/rlenfac)*omorirate((x-ruplenOmori), chi2, ruplen2, q)]
			#Yth2+=[omorirate((x-ruplen), chi2, ruplen2, q)]
			YthPL+=[.5*(1.0/rlenfac)*rate0*(ruplenOmori**q)*(x**(-q)) ]
	
	Yth=scipy.array(Yth)
	Yth2=scipy.array(Yth2)	
	Ythcore=scipy.array(Ythcore)
	YthPL=scipy.array(YthPL)
	#print "some lens: ", len(xth), len(Yth), len(Ythcore), len(Yth2)
	#
	plt.figure(fnum, figsize=(7,5))
	plt.clf()
	plt.figure(fnum+6, figsize=(7,5))
	ax=plt.gca()
	ax.set_yscale('log')
	plt.clf()
	i=-1
	ax=plt.gca()
	#ax.set_xscale('log')
	#ax.set_yscale('log')
	#disttypes=['0','3','4','3d']
	#dtypeindex=0
	#Yth=[]
	#dtypes={0:'rect', 1:'spher.', 2:'lib', 3:'3D'}
	dtypes={0:'rect', 1:'spher.', 2:'2D', 3:'3D'}
	mysymbs=['o', 's', '^', 'h', '+']
	for r in [r0s, r3s, r4s, r3d]:
		i+=1
		#if r==r0s: continue
		if i in (0, 1): continue
		#
		#plt.figure()
		plt.figure(fnum, size=(7,10))
		#h1=plt.hist(r, bins=nbins, histtype='step', label='dist-type: %s' % dtypes[i])
		ax0=plt.gca()
		ax0.set_yscale('log')
		ax0.set_xscale('linear')
		h1=plt.hist(r, bins=nbins, histtype='bar', alpha=.25, label='dist-type: %s' % dtypes[i])
		#
		plt.figure(fnum+6, size=(7,10))
		ax=plt.gca()
		ax.set_yscale('log')
		ax.set_xscale('log')
		X=scipy.array(h1[1][0:-1])+(h1[1][1]-h1[1][0])*.5
		binwidth=h1[1][1]-h1[1][0]
		plt.plot(X, h1[0]/binwidth, '%s-' % mysymbs[i%len(mysymbs)], lw=2, label='dist-type: %s' % dtypes[i], alpha=.45, zorder=1)
		#plt.fill_between(X, h1[0]/binwidth, min(h1[0])/binwidth, alpha=.2)
		
		#plt.plot(xth, Yth, '--', lw=3, alpha=.7, zorder=42, label='SPL Theory')
		#plt.plot(xth, Yth2, '--', lw=3, alpha=.7, zorder=42, label='Sat. Theory')
		'''
		plt.figure(fnum+12)
		ax=plt.gca()
		ax.set_yscale('log')
		ax.set_xscale('log')
		#X=scipy.array(h1[1][0:-1])+(h1[1][1]-h1[1][0])*.5
		X2=X/(10.0**(lrlen+1))
		plt.plot(X2, h1[0], 'o-', label='dist-type: %s' % dtypes[i])
		'''
		#plt.clf()
		#
		# this plot is a nice idea, but because aftershock locations are 2D, we don't really get much from r_i-r_{i-1},
		# because the two events are not necessarily co-located (as are 1D temporal distributions). what i think i'm 
		# properly looking for is a nearest-neighbors plot, which is a bit more involved to produce and can arguably be
		# substituted with binned data. 
		'''
		r.sort()
		z=plt.semilogy(r[1:],scipy.array(r[1:])-scipy.array(r[0:-1]), '.', label='dist-type: %d' % i)
		i+=1
		#dtypeindex+=1
		plt.ylabel('$\\Delta r = r_{i+1} - r_i$')
		plt.xlabel('$r$')
		plt.title('%s' % eqset)
		'''
		#i+=1
		plt.figure(fnum)
	vz=max(h1[0])*2.0
	zcopy=scipy.array(h1[0]).copy()
	zcopy.sort()
	for icpy in xrange(len(zcopy)):
		if zcopy[icpy]>0:
			vzmin=zcopy[icpy+0]
			break
		#
	#
	vzmin*=.5
	zcopy=None
	#
	for fn in [fnum, fnum+6]:
		plt.figure(fn)
		#plt.plot(xth, Yth, '--', lw=2, alpha=.7, zorder=42)
		#plt.plot(xth, Yth2, '--', lw=2, alpha=.7, zorder=42)
		if fn==fnum:
			binwidthfact=binwidth
			'''
			plt.plot([.5*10**lrlen, .5*10**lrlen], [vzmin, vz], '--d', lw=3, label='$0.05 \cdot L_r$')
			plt.plot([5.*10**lrlen, 5.*10**lrlen], [vzmin, vz], '--d', lw=3, label='$L_r/2$')
			plt.plot([10**(1.0+lrlen), 10**(1.0+lrlen)], [vzmin, vz], '--d', lw=3, label='$L_r$')
			plt.plot(xth, Yth*binwidth, 'm--', lw=3, alpha=.9, zorder=9, label='SPL Theory $L_r$')
			plt.plot(xth, Ythcore*binwidth, 'c--', lw=3, alpha=.9, zorder=9, label='SPL Theory ssim')
			plt.plot(xth, Yth2*binwidth, 'k--', lw=3, alpha=.9, zorder=9, label='Sat. Theory')
			plt.plot(xth, YthPL*binwidth, 'b-', lw=2, alpha=.9, zorder=9, label='PL-thresh Theory')
			'''
		else:
			binwidthfact=1.0
			'''
			plt.plot([.5*10**lrlen, .5*10**lrlen], scipy.array([vzmin, vz])/binwidth, '--d', lw=3, label='$0.05 \cdot L_r$')
			plt.plot([5.*10**lrlen, 5.*10**lrlen], scipy.array([vzmin, vz])/binwidth, '--d', lw=3, label='$L_r/2$')
			plt.plot([10**(1.0+lrlen), 10**(1.0+lrlen)], scipy.array([vzmin, vz])/binwidth, '--d', lw=3, label='$L_r$')
			plt.plot(xth, Yth, 'r-', lw=4, alpha=1., zorder=9, label='SPL Theory $L_r$')
			plt.plot(xth, Ythcore, 'm-', lw=4, alpha=.9, zorder=9, label='SPL Theory ssim')
			plt.plot(xth, Yth2, 'k-', lw=4, alpha=.9, zorder=9, label='Sat. Theory')
			plt.plot(xth, YthPL, 'b-', lw=2, alpha=.9, zorder=9, label='PL Theory')
			'''
		plt.plot([.5*10**lrlen, .5*10**lrlen], scipy.array([vzmin, vz])/binwidth, '--d', lw=3, label='$0.05 \cdot L_r$')
		plt.plot([5.*10**lrlen, 5.*10**lrlen], scipy.array([vzmin, vz])/binwidth, '--d', lw=3, label='$L_r/2$')
		plt.plot([10**(1.0+lrlen), 10**(1.0+lrlen)], scipy.array([vzmin, vz])/binwidth, '--d', lw=3, label='$L_r$')
		plt.plot(xth, Yth*binwidthfact, 'r-', lw=4, alpha=1., zorder=9, label='SPL Theory $L_r$')
		plt.plot(xth, Ythcore*binwidthfact, 'm-', lw=4, alpha=.9, zorder=9, label='SPL Theory ssim')
		plt.plot(xth, Yth2*binwidthfact, 'k-', lw=4, alpha=.9, zorder=9, label='Sat. Theory')
		plt.plot(xth, YthPL*binwidthfact, 'b-', lw=2, alpha=.9, zorder=9, label='PL Theory')
		#
		# and  horizontal lines (estimated initial rates):
		dx=(max(X)-min(X))
		hfact=1.0
		if fn==fnum: hfact=binwidth
		#
		hXs=[min(X)*.5, max(X)*2.]
		#
		#plt.plot(hXs, hfact*scipy.array([10**estRateExp, 10**estRateExp]), 's-.', lw=3, label='$\\log(N\') = m/2$')
		#plt.plot(hXs, hfact*scipy.array([10**(.8*estRateExp), 10**(.8*estRateExp)]), 's-.', lw=3, label='$\\log(N\') = 0.4 \\cdot m$')
		#
		#x1box=10**(-2.0+lrlen)
		#x2box=10**(1.0+lrlen)
		x1box, x2box = .01*ruplenOmori, ruplenOmori
		#plt.fill_between([.05*ruplen, ruplen], hfact*scipy.array([10**estRateExp, 10**estRateExp]), [10**(.8*estRateExp), 10**(.8*estRateExp)], alpha=.25, color='m')
		#
		#
		fillmstar = mstar
		if mstar>mt: fillmstar=mt
		#
		topfillY = 10.**maxDensityExp(mstar=fillmstar, mc=thismc, dm=1.0, D=2.0)
		bottomfillY = 10.**maxDensityExp(mstar=fillmstar, mc=thismc, dm=1.0, D=1.0)
		plt.fill_between([x1box, x2box], .5*(1.0/rlenfac)*hfact*scipy.array([topfillY, topfillY]), .5*(1.0/rlenfac)*hfact*scipy.array([bottomfillY, bottomfillY]), alpha=.25, color='m')
		plt.plot([x1box, x2box, x2box, x1box, x1box], .5*(1.0/rlenfac)*hfact*scipy.array([topfillY, topfillY, bottomfillY, bottomfillY, topfillY]), 'm-', label='$\\log(\\lambda(r)_{max} ) = %.2f$' % estRateExp )
		#
		#plt.fill_between([x1box, x2box], .5*(1.0/rlenfac)*hfact*scipy.array([10**estRateExp, 10**estRateExp]), .5*(1.0/rlenfac)*hfact*scipy.array([10**(.8*estRateExp), 10**(.8*estRateExp)]), alpha=.25, color='m')
		#plt.plot([x1box, x2box, x2box, x1box, x1box], .5*(1.0/rlenfac)*hfact*scipy.array([10**estRateExp, 10**estRateExp, 10**(.8*estRateExp), 10**(.8*estRateExp), 10**estRateExp]), 'm-', label='$\\log(\\lambda(r)_{max} ) = %.2f$' % estRateExp )
		#plt.fill_between([x1box, x2box], scipy.array([10**estRateExp, 10**estRateExp]), scipy.array([10**(.8*estRateExp), 10**(.8*estRateExp)]), alpha=.25, color='m')
		# , label='$\\log(\\lambda(r)_{max} ) = %.2f$' % estRateExp 
		#plt.plot([x1box, x2box, x2box, x1box, x1box], scipy.array([10**estRateExp, 10**estRateExp, 10**(.8*estRateExp), 10**(.8*estRateExp), 10**estRateExp]), 'm-', label='$\\log(\\lambda(r)_{max} ) = %.2f$' % estRateExp )
		#
		# Ninsides
		#for k in [0,2,3]:
		#for k in [1,2,3]:
		for k in [2,3]:
			#plt.plot([.02*ruplenOmori, ruplenOmori*1.25], hfact*scipy.array([Ninsides[k], Ninsides[k]]), '--', lw=3, label='mean,type: %s' % dtypes[k])
			plt.plot([.02*ruplenOmori, ruplenOmori*10.], hfact*scipy.array([Ninsides[k], Ninsides[k]]), '--', lw=3, label='mean (%s)' % dtypes[k])
		#
		if not hidelegend: plt.legend(loc='best')
		#
	#
	#plt.figure(fnum+6)
	# plt.legend(loc='best')
	
	#return None
	return [mstar, thismc, rlenfac, q, estRateExp, Ninsides]

def fitEqInts(eqset='parkfield', nrange=[1,250], catnum=None, p=1.01, mc=None, dalpha=2.0, dmstar=1.0, forcep=False, doMCfit=False, mcnits=10**5, mcNprocs=None):
	plt.ion()
	if catnum==None:
		catnum=0
		if eqset=='parkfield': catnum=1
	#
	namedict={'parkfield': 'Parkfield', 'hmine': "Hector Mine", 'elmayor': 'El Mayor-Cucapah', 'tohoku': 'Tohoku-oki', 'sansim': 'San Simeon', 'coalinga': 'Coalinga'}
	#
	cpf=getETAScat(catname=eqset, mc=mc)
	eqints = getEqInts(cpf, catnum=catnum)		# [t, dt]	where t is seconds since ms.
	nskip=int(float(len(eqints[0]))/100.)
	#
	X=scipy.array([])
	Y=scipy.array([])
	for i in xrange(len(eqints[0])):
		thisdt=(eqints[1][i])
		if thisdt==0.0: continue
		Y+=[math.log10(thisdt)]
		X+=[math.log10(eqints[0][i])]
	#
	mstar=cpf.mainshock[3]
	prams0 = [-.5, 0., 5.]
	prams0 = scipy.array([2.0, 0.001])
	if forcep==True: prams0 = scipy.array([2.0])
	#prams0 = scipy.array([2.0])
	# omoriError(prams, t, dts, m, mc, p0=1.2):
	myones=scipy.ones(len(X))
	fp=p
	plsqInts=spo.leastsq(omoriError2, prams0, args=(eqints[0][nskip:], eqints[1][nskip:], cpf.mainshock[3], cpf.mc, fp))[0]
	#plsqInts=spo.leastsq(omoriError2, prams0, args=(eqints[0][0:], eqints[1][0:], cpf.mainshock[3], cpf.mc, fp))[0]
	#plsqInts=spo.leastsq(omoriError, prams0, args=(eqints[0], eqints[1], cpf.mainshock[3], cpf.mc, 1.2))[0]
	#
	ft0=10**(plsqInts[0])
	#
	print "plsqInts: ", plsqInts
	#
	yth=[]
	yth2=[]
	ythraw=[]
	ythssim=[]
	ythssim2=[]
	ythssimcum=[]
	xth=[]
	mprime=[]
	mfactor=.5	# scaling factor for rupture length/duration
	alpha=2.28
	#
	# Tohoku seems to length (and therefore delta-t probably) scale in 2D, so let's skip this and see how it goes.
	'''
	if mstar>7.6:
		mfactor=1.0
		alpha=6.28
		# nominally, m>7.6 should use dtau~2.5 as well (maybe)
		# and what's up? the Tohoku-oki earthquake seems to length-scale more like a 2D event. what abour delta t?
	'''
	#
	Ngr = 10**(1.0*(mstar-dmstar-cpf.mc))
	yodat0=10**(mstar*mfactor - alpha + dalpha)
	yodatau = (yodat0**(1.0-p))/((p-1.0)*Ngr)
	lt0raw = mstar*mfactor - alpha
	t0raw=10**lt0raw
	ltauraw = (1.0-p)*lt0raw - math.log10(p-1.0) - math.log10(Ngr)
	# t0 from self-sim arguments...
	Nom = 10.0**(1.0*(mstar- dmstar - mc))
	dtau = 2.28
	#
	# Self Similar formulations:
	# use this funct. for uniformity...
	dmtau=0.0
	t0tau = t0taussim(mstar, mc, dmtau=dmtau, dmstar=dmstar, b=1.0, p=p)
	t0ssim = t0tau[0]
	taussim = t0tau[1]
	#
	dmtau=6.0	# eventually, let's use the mean <dmtau> value.
	#dmtau = 6.84	# (elmayor)
	t0tau = t0taussim(mstar, mc, dmtau=dmtau, dmstar=dmstar, b=1.0, p=p)
	t0ssim2 = t0tau[0]
	taussim2 = t0tau[1]
	#
	print "m=%f: t0=%f, p=%f, deltat= %f - %f, t0ss=%f" % (cpf.mainshock[3], ft0, fp, 10**(mstar/2.-2.28), 10**(mstar/2.-2.58), t0ssim)
	
	print "yodat0=%f, yodatau=%f, dt0=%f, rupture=%f" % (yodat0, yodatau, yodatau*yodat0**p, 10**(mstar/2.0 - alpha))
	#for i in xrange(len(eqints[0])):
	lxth=math.log10(eqints[0][0])-4.0
	dlxth=.1
	while lxth<math.log10(eqints[0][-1]):
		#x=eqints[0][i]
		x=10**lxth
		xth+=[x]
		yth+=[((ft0**(1.0-fp))*(ft0+x)**fp)/((fp-1.)*Ngr)]		# note: intervals, not rate...
		ythssim+=[((t0ssim**(1.0-fp))*(t0ssim+x)**fp)/((fp-1.)*Ngr)]
		ythssim2+=[((t0ssim2**(1.0-fp))*(t0ssim2+x)**fp)/((fp-1.)*Ngr)]
		#
		ythssimcum+=[((t0ssim+x)**(1.-p) - t0ssim**(1.0-p))/((1.0-p)*taussim)]
		#
		#myx=-ft0+eqints[0][i]
		#if myx<0: continue
		#yth2+=[((ft0**(1.0-fp))*(myx)**fp)/((fp-1.)*Ngr)]
		yth2+=[(yodatau*(yodat0 + x)**p)]
		ythraw+=[(10**(ltauraw))*(t0raw+x)**p]
		mprime += [(1./mfactor)* (alpha + ltauraw + p*math.log10(t0raw+x))]
		lxth+=dlxth
	#
	yth=scipy.array(yth)
	#plt.figure(2)
	#plt.loglog(eqints[0], yth, '-')
	#plt.loglog(eqints[0], yth2, '--', lw=3, alpha=.7)
	#
	# Interval plot:
	#
	plt.figure(1)
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	#
	plt.plot(eqints[0], eqints[1], 'b.', label='$intervals$ $\\Delta t$', alpha=.7, lw=2.5)
	# skipints=scipy.array(eqints[0][nskip:]) - eqints[0][nskip-1]
	#plt.loglog(eqints[0][nskip:], eqints[1][nskip:], 'y.', alpha=.7)
	plt.plot(xth, yth, 'g-', label='$t_{0 fit}=%.3f$' % ft0, lw=2.5)
	plt.plot(xth, ythraw, 'r--', label='$t_{0} = \Delta t_r=%.3f$' % t0raw, lw=2.5)
	#plt.loglog(xth, yth2, 'm--', label='$t_{0,\\Delta \\alpha=%.2f}=%.3f$' % (dalpha, yodat0), lw=2.5, alpha=.7)
	plt.plot(xth, yth2, 'm--', label='$%.2f \cdot \Delta t_{r}=%.3f$' % (10**dalpha, yodat0), lw=2.5, alpha=.7)
	plt.plot(xth, ythssim, 'c--', label='$t_{0 ssim} = %.2f$' % (t0ssim), lw=2.5, alpha=.7)
	plt.loglog(xth, ythssim2, 'y--', label='$t_{0 ssim2} = %.2f$' % (t0ssim2), lw=2.5, alpha=.7)
	thism=cpf.mc
	tc = 10**((thism*mfactor - alpha - ltauraw)/p)
	dmstar=1.
	thism=mstar-dmstar
	tbath = 10**((thism*mfactor - alpha - ltauraw)/p)
	#ttrans = 10**((7.6 - alpha - ltauraw)/p)
	#
	#plt.plot(xth, mprime, '-', label='mprime')
	xstar = [tc, tbath]
	ystar = [10**(cpf.mc*mfactor - alpha), 10**((mstar-dmstar)*mfactor - alpha)]
	#
	plt.scatter([tc], [10**(cpf.mc*mfactor - alpha)], marker='*', s=175, c='m', label='$m_{max} = %.2f$' % cpf.mc, zorder=10)
	plt.scatter([tbath], [10**(thism*mfactor - alpha)], marker='*', s=175, c='c', label='$m_{max} = %.2f$' % thism, zorder=10)
	
	print "tc-tbath bit: ", math.log10(tc), math.log10(tbath), math.log10(tbath)/math.log10(tc)

	plt.legend(loc='best', title='%s, $p=%.2f$' % (namedict[eqset], p), scatterpoints=1, numpoints=1)
	plt.xlabel('time $t$ (seconds) since mainshock', size=20)
	plt.ylabel('interval $\\Delta t$ (seconds)', size=20)
	#
	# ############################
	# Cumulative Dists.:
	plt.figure(2)
	plt.clf()
	ax2=plt.gca()
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	plt.xlabel('time $t$ since mainshock', size=16)
	plt.ylabel('Number of earthquakes $N$', size=16)
	Ns=range(1, len(ythssimcum)+1)
	Nsints = range(1, len(eqints[0])+1)
	#
	plt.plot(eqints[0], Nsints, 'bs', label='raw data', zorder=2, ms=8, alpha=.6)
	#nskip=int(float(len(eqints[0]))/100.)
	skipints=scipy.array(eqints[0][nskip:]) - eqints[0][nskip-1]
	Nskips= range(1, len(skipints)+1)
	plt.plot(skipints,Nskips, 'go', ms=8, alpha=.4, label='skip-data', zorder=2)
	plt.plot([min(skipints), max(skipints)], [Nom, Nom], 'k--', lw=3, alpha=.7, label='$N_{Omori}=%d$' % int(Nom))
	plt.plot(scipy.array(xth), scipy.array(ythssimcum), 'r--', label='Raw theory: $p=%.2f$, $\Delta m_\\tau = 0.0$' % fp, lw=2, zorder=5)
	#
	# fitting:
	# plsq=spo.leastsq(resAlphaOmori, P, args=(scipy.array(Y), scipy.array(X), mstar, mc, p, scipy.ones(len(X)), alpha, beta, gamma))
	#plsq=spo.leastsq(resOmori1, [0.0], args=(scipy.array(Nsints), scipy.array(eqints[0]), scipy.ones(len(Nsints))*mstar, scipy.ones(len(Nsints))*mc))
	
	if doMCfit:
		miniq=mcp.Queue()
		fitnits=10**4
		nprocs=mcNprocs
		if nprocs==None: nprocs=mcp.cpu_count()
		Procs = []
		for pcount in xrange(nprocs):
			Procs+=[mcp.Process(target=mcomorifit1MCPwrapper, args=(scipy.array(skipints), scipy.array(Nskips), mstar, mc, fitnits, 1.08, -2.0, 2.28, 1.0, 1.0, miniq))]
			Procs[-1].start()
		#
		for pcount in xrange(nprocs):
			Procs[pcount].join()
		#
		#P1=mcp.Process(target=mcomorifit1MCPwrapper, args=(scipy.array(skipints), scipy.array(Nskips), mstar, mc, fitnits, 1.08, -2.0, 2.28, 1.0, 1.0, miniq))
		#P2=mcp.Process(target=mcomorifit1MCPwrapper, args=(scipy.array(skipints), scipy.array(Nskips), mstar, mc, fitnits, 1.08, -2.0, 2.28, 1.0, 1.0, miniq))
		#P1.start()
		#P2.start()
		#P1.join()
		#P2.join()
		#
		#	
		#print plsq1, plsq2, miniq
		#plsq1=miniq.get()
		#plsq2=miniq.get()
		#print plsq1, plsq2
		#
		plsqCum=miniq.get()
		for i in xrange(nprocs-1):
			thisPLSQ=miniq.get()
			#print thisPLSQ
			if thisPLSQ[4]<plsqCum[4]: plsqCum=thisPLSQ
		
		#if plsq2[4]<plsq1[4]: plsqCum=plsq2
		#print plsq1, plsq2
		print "mc-plsq: ", plsqCum
	else:
		# and now we've done the fit. this is a stupid way to do this, but for now, just record the fit-prams here in hard-code.
		# [p, dmtau, t0, tau, err]
		# these fits from a 2*10**6 nits MC fit (see above)
		# [ [p, dmtau, t0, tau, meanErr], []...]
		# and these are all wrong now... we'll write them to file and pull from that.
		eqprams = {'hmine':[1.0800709883478472, -1.8629651135639544, 10749.57000220929, 0.004717667640749667, 0.00053782841592927568], 'sansim':[1.219638848698415, -2.2719093509266317, 11021.851443270456, 0.0018640377582269854, 0.0070431527275239166], 'coalinga':[1.3297329878916768, -2.7923573636720778, 29219.523838853547, 0.0006447002764860759, 0.0073651599763037241], 'parkfield':[1.080034455491961, -3.0986557494550753, 34407.51608162181, 0.001834993064142981, 0.033581325763624943],'tohoku':[1.255982454583056, -0.69629903172771, 123504.2912057097, 4.8798731442385375e-05, 0.001196232406482376], 'elmayor':[1.2133037703624776, 6.835948557592083, 87818.74207724657, 8.251172686901148e-05, 0.0017337615307428645]}
		plsqCum=eqprams[eqset]
	#
	#plsq = mcomorifit1(scipy.array(eqints[0]), scipy.array(Nsints), mstar, mc, nits=10000)
	#plsq = mcomorifit1(scipy.array(skipints), scipy.array(Nskips), mstar, mc, nits=10000)
	#print "plsq: ",  plsq
	#
	# Plot the MC fit:
	Yfit = ybp.omoriRateInt(0., scipy.array(xth), t0=plsqCum[2], tau=plsqCum[3], p=plsqCum[0])
	plt.plot(xth, Yfit, 'c-', label='MC fit: $p=%.2f, \Delta m_\\tau = %.2f$' % (plsqCum[0], plsqCum[1]), lw=2, zorder=5 )
	#
	plt.legend(loc='best', numpoints=1, title = namedict[eqset])
	
	#vars: interval fits: {ft0, t0ssim}, cumfits: {p=plsqCum[0], dmtau=plsqCum[1]}
	#			# but note we get the cum. fits from a MC method, so we just catalogued them above.
	return [cpf.mainshock[3], cpf.mc, p, ft0, t0ssim, plsqCum[0], plsqCum[1]]

def mcomorifit1MCPwrapper(X, Y, m, mc, nits=10**5, p=1.08, dmtau=-2.0, dtau=2.28, b=1.0, dmstar=1.0, miniq=None):
	z=mcomorifit1(X, Y, m, mc, nits, p, dmtau, dtau, b, dmstar)
	miniq.put(z)
	
def mcomorifit1(X, Y, m, mc, nits=10**5, p=1.08, dmtau=-2.0, dtau=2.28, b=1.0, dmstar=1.0):
	#
	R1=random.Random()
	R2=random.Random()
	#
	minerr=10**12.
	minprams=[]
	for i in xrange(nits):
		# random sample the fitting parameters {dmtau, p}
		dmtau = 4.5 + 4.0*R1.random()
		p = 1.01 + .30*R2.random()
		#
		thist0tau = t0taussim(m, mc, dmtau=dmtau, dmstar=1.0, b=1.0, p=p)
		t0=thist0tau[0]
		tau = thist0tau[1]
		#
		# 
		#y = ((t0 + X)**pm - t0**pm)/(pm*tau)
		y = ybp.omoriRateInt(0., scipy.array(X), t0=t0, tau=tau, p=p)
		#
		errs = scipy.array(map(math.log10, Y)) - scipy.array(map(math.log10, y))
		thiserr = sum(errs**2.0)/float(len(y)-3.)
		#
		if thiserr<minerr:
			print "new min: %f, (%f, %f)" % (thiserr, p, dmtau)
			minerr=thiserr
			minprams=[p, dmtau, t0, tau, minerr]
		#
	#print "miniq: ", miniq
	#if miniq!=None: miniq.put(minprams)
	#print "miniq: ", miniq
	return minprams

def t0taussim(mstar, mc, dmtau=0.0, dmstar=1.0, b=1.0, p=1.08):
	dtau = 2.28
	lNgr = b*(mstar - dmstar - mc)
	Ngr=10.0**lNgr
	#dmtau=-2.0
	#lt0 = 7.0*mstar/6.0 - (2.0*mc/3.0) + (dmtau/3.0) - (2./3.)*math.log10(1.5) - dtau + math.log10(p-1.)
	lt0 = 7.0*mstar/6.0 - (2.0*mc/3.0) - dtau - dmstar - (2./3.)*math.log10(1.5) + (dmtau/3.0) + math.log10(p-1.0)
	t0 = 10.**lt0
	tau = (t0**(1.0-p))/(Ngr*(p-1.0))
	#
	return [t0, tau]

def resOmori1(prams, Y, X, mstar, mc):
	# and this just isn't working. we'll probably need a MC fit.
	# prams are [dmtau, p]
	# assuming self-sim:
	dtau=2.28
	Ngr = 1.0*(mstar - 1.0 - mc)
	print 10.**prams[0]
	dmtau=10**prams[0]
	#p=prams[1]
	p=1.08
	pm = 1.0-p
	#
	t0tau = t0taussim(mstar, mc, dmtau, 1.0, 1.0, p)
	t0 = t0tau[0]
	tau = t0tau[1]
	#lt0 = 7.0*mstar/6.0 - (2.0/3.0)*(dmtau+mc) - (2./3.)*math.log10(1.5) - dtau + math.log10(p-1.)
	#t0 = 10.**lt0
	#tau = (t0**pm)/(-Ngr*pm)
	#
	# this is a cumulative residual...
	y = ((t0 + X)**pm - t0**pm)/(pm*tau)
	#
	errs=scipy.array(map(math.log10, y))-scipy.array(map(math.log10, Y))
	return errs

#######
def fitscript(prange=[1.01, 1.26], dp=.05):
	eqs=['parkfield', 'coalinga', 'sansim', 'elmayor', 'hmine']
	#prange=[1.01, 1.26]
	#dp=.05
	outses=[]	# [ [name, m, mc, p, t0], ..., ]
	thisp=prange[0]
	while thisp<prange[1]:
		print "p=%f" % thisp
		ms=[]
		mcs=[]
		t0s=[]
		for eq in eqs:
			cnum=0
			if eq=='parkfield': cnum=1
			thisprams=fitEqInts(eq, catnum=cnum, p=thisp, mc=3.5)
			outses+=[thisprams]
			#
			ms+=[thisprams[0]]
			mcs+=[thisprams[1]]
			t0s+=[math.log10(thisprams[3])]
		#
		thisplsq = spo.leastsq(doodleerror2, [0., 0., 0.], args=(scipy.array(ms), scipy.array(mcs), scipy.array(t0s), scipy.ones(len(t0s))), full_output=0 )
		print "p, plsq:", thisp, thisplsq
		#
		thisp+=dp
	#
	Ms=map(operator.itemgetter(0), outses)
	Mcs=map(operator.itemgetter(1), outses)
	ps=map(operator.itemgetter(2), outses)
	t0s=map(operator.itemgetter(3), outses)
	aryt0s=scipy.array(t0s)
	#
	'''
	MMs, pps1=numpy.meshgrid(set(Ms), set(ps))
	MMcs, pps2=numpy.meshgrid(set(Mcs), set(ps))
	Zs1 = scipy.array(t0s).copy()
	Zs1.shape=len(ps), len(Ms)
	Zs1=
	'''
	#
	#plsq = spo.leastsq(doodleerror2, p0, args=(scipy.array(ms), scipy.array(mcs), scipy.array(ndots), Ws), full_output=1 )
	#
	plt.figure(0)
	plt.clf()
	#plt.contourf(Ms, ps, t0s, contours=10)
	#
	plt.figure(1)
	plt.clf()
	#plt.contourf(Mcs, ps, t0s, contours=10)
	#
	plt.figure(2)
	

	return outses
	#
def fitEqIntsXY(X,Y, nrange=[1,250]):
	# monte-carlo fit
	bestprams=[]	#[a,b,n,a,b]
	minerr=27000.0
	n=nrange[0]
	Ndof = float(len(X)-4.0)
	#
	# for customizing, as necessary:
	#while X[-1]>6.25:
	#	X.pop()
	#	Y.pop()
	#
	while n<nrange[-1]:
		n+=1
		f1=lft.linefit([X[0:n], Y[0:n]])
		f2=lft.linefit([X[n:], Y[n:]])
		f1.doFit()
		f2.doFit()
		#print "errs: %f/%f, %f/%f" % (f1.totalVar, f1.totalVar/n, f2.totalVar, f2.totalVar/(float(len(X)-n)))
		#
		err = (f1.totalVar + f2.totalVar)/Ndof
		#err = f1.meanVar() + f2.meanVar()
		if err<minerr:
			minerr=err
			bestprams=[f1.a, f1.b, n, f2.a, f2.b]
			print "(%d) newBestprams (%f): %s" % (n, minerr, bestprams)
		#
	#
	print "meanDX0 = %f, %f" % (scipy.array(Y[0:bestprams[2]]).mean(), scipy.array(Y[0:bestprams[2]]).std())
	#
	fig1=plt.figure(0)
	plt.clf()
	ax=plt.gca()
	#ax.set_yscale('log')
	#ax.set_xscale('log')
	#
	#plt.plot(eqints[0][bestprams[2]], eqints[1][bestprams[2]], 'r*', ms=15, zorder=11)
	#plt.plot(eqints[0], eqints[1], 'b.-.')
	#x=[eqints[0][0], eqints[0][-1]]
	#
	plt.plot(X[bestprams[2]], Y[bestprams[2]], '*', ms=15, zorder=11)
	plt.plot(X, Y, '.-.')
	x=[X[0], X[-1]]
	#
	#x=[math.log10(eqints[0][0]), math.log10(eqints[0][-1])]
	#plt.figure(1)
	#y=[(10**f1.a)*x[0]**f1.b, (10**f1.a)*x[-1]**f1.b]
	y=[f1.a + x[0]*f1.b, f1.a + x[-1]*f1.b]
	plt.plot(x,y, '-')
	y=[f2.a + x[0]*f2.b, f2.a + x[-1]*f2.b]
	plt.plot(x,y, '-')
	#plt.plot(x, [(10**f1.a)*x[0]**f1.b, (10**f1.a)*x[-1]**f1.b], '-')
	#plt.plot(x, [(10**f2.a)*x[0]**f2.b, (10**f2.a)*x[-1]**f2.b], '-')
	#plt.plot(x, [(f1.a) + x[0]*f1.b, f1.a + x[-1]*f1.b], '-')
	#plt.plot(x, [(f2.a) + x[0]*f2.b, (f2.a) + x[-1]*f2.b], '-')
	
	return [X,Y]

def omoriError(prams, t, dts, m, mc, p0=1.2):
	alpha=prams[0]
	beta=prams[1]
	gamma=prams[2]
	#
	dm=1.0
	b=1.0
	#
	lt0=alpha*m + beta*mc + gamma
	t0=10**lt0
	#print "prams: ", alpha, beta, gamma, m, mc, p0, lt0
	#print "type: ", type(lt0).__name__, len(lt0), " :: ", lt0
	lNGR = b*(m-dm-mc)
	
	lts=scipy.array(map(math.log10, (t+t0)))
	ldts=scipy.array(map(math.log10, dts))
	#lts=t
	#print "dologp0"
	#logp0=map(math.log10, (p0[0]-1.))
	logp0 = math.log10(p0-1.)
	#
	#thiserr=0.0
	ltheory=logp0*lt0 + p0*lts - logp0 - lNGR
	#ldeltats=(1.0-p0)*lt0 + p0*lts - logp0 - lNGR
	#
	return (ldts-ltheory)**2.0

def omoriError2(prams, t, dts, m, mc, p):
	lt0=prams[0]
	t0=10**lt0
	if len(prams)>1: 
		p=1.0 + prams[1]		# replace the p parameter with the dynamic(fitting) pram.
	#p=1.2
	dm=1.0
	b=1.0
	#	
	#lt0=t0
	lNGR = b*(m-dm-mc)
	#
	lts, ldts=[],[]
	for i in xrange(len(t)):
		if dts[i]==0: continue
		lts+=[math.log10(t[i]+t0)]
		#lts+=[math.log10(t[i])]
		ldts+=[math.log10(dts[i])]
	lts=scipy.array(lts)
	ldts=scipy.array(ldts)
	#
	#lts=scipy.array(map(math.log10, (t+t0)))
	#ldts=scipy.array(map(math.log10, dts))
	#
	#print "p: ", p
	#logp0 = math.log10(p-1.)
	#
	#thiserr=0.0
	ltheory = (1.-p)*lt0 + p*(lts) - math.log10(p-1.) - lNGR
	#
	errs=(ldts-ltheory)
	
	#return (ldts-ltheory)**2.0
	return errs

def omorirate(t, tau, t0, p):
	return (1.0/tau) * (t0+t)**(-p)

def omoriN(t, tau, t0, p):
	N=(1.0/((1.0-p)*tau)) * ((t0+t)**(1.0-p) - t0**(1.0-p))
	return N

def gettau(N, t0, p):
	return (t0**(1.0-p))/((1.0-p)*N)
	
def doodleerror(m, mc, ndot, a, b, c):
	return (a*m+b*mc+c - ndot)**2.0
def doodleerror2(p, m, mc, ndot, ws=None):
	a=p[0]
	b=p[1]
	c=p[2]
	if ws==None: ws=scipy.ones(len(ndot))
	#print "internal: ", a,b,c
	return ws*(a*m+b*mc+c - ndot)**2.0

def makeDistTable(fname='data/productionfits/distfits-rlf50-q135.dat'):
	f=open(fname, 'r')
	filebits=fname.split('-')
	outfile='data/productionfits/disttable-temp-rl%s-q%s.txt' % (filebits[1], filebits[2])
	fout = open(outfile, 'w')
	#
	#fout.write('\% latex for distance distribution table, q=%s, rlr=%s\n\n' % (filebits[1], filebits[2]))
	fout.write('\\begin{table}\\begin{tabular}{|r||c|c|c|c|c||c|c|c|} \\hline textbf{earthquake} & \\textbf{$m$}  & \\textbf{$m_c$} & \\textbf{$\log ($ \\Delta L_r) $} & \\textbf{\\alpha_{Lr}} & \\textbf{$q$} & \\textbf{Eq. ()} \\textbf{$\\log \\left( \\left \\langle \\frac{N}{L_r} \\right \\rangle \\right )$} & textbf{$D$} \\\\\n')
	for rw in f:
		if rw[0]=='#': continue
		#
		rws=rw.split()
		eqname=rws[0]
		mstar, mc, lrlen1, lrlen2, rlenfac, q, lestrate = map(float, rws[1:8])
		Ns=map(float, rws[8:])
		D=float(rws[12])
		Derr=float(rws[13])
		#
		#thislN = math.log10(.5*(Ns[1]+Ns[2]))	# algebraic mean of spher + lib methods
		thislN = .5*(math.log10(Ns[1]*Ns[2]))	# geometric mean of spher+lib methods
		print "Nvals (%s): %f, %f, <%f>" % (rws[0], Ns[1], Ns[2], thislN)
		foutstr = '\\textbf{%s} & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f  & %.2f $\pm$ %.2f\\\\ \n' % (eqname, mstar, mc, lrlen1, rlenfac, q, lestrate, thislN, D, Derr)
		if mstar>7.6:
			foutstr = '\\textbf{%s} & %.2f & %.2f & %.2f (%.2f) & %.2f & %.2f & %.2f & %.2f & %.2f $\pm$ %.2f \\\\\n' % (eqname, mstar, mc, lrlen1, lrlen2, rlenfac, q, lestrate, thislN, D, Derr)
		fout.write(foutstr)
	#
	fout.write('\\hline \\end{tabular} \\caption{\\label{tbl:intOmoriFits} Best fit variables, cumulative interval distributions. } \\end{table}')
#
# production plots for figs. 4. do we also have a script for figs. 5? (PF: r0-> Lr, not Lr/2, and for Tohoku-oki: q=2.0(b), for c)D=1.0).
def makeDistFigs(qs=[1.01, 1.05, .01], dalpha=2.0, nbins=100, q=1.35, rlenfac=.5, dosave=True, D=1.5):
	eqDict = {'hmine': ['Hector Mine', 7.1, 3.0, 0], 'parkfield':['Parkfield', 5.97, 1.5, 1], 'sansim':['San Simeon', 6.5, 3.0, 0], 'coalinga':['Coalinga', 6.7, 3.5, 0], 'tohoku':['Tohoku-oki', 9.0, 4.5, 0], 'elmayor':['El Mayor-Cucapah', 7.2, 2.5, 0]}
	#
	#a=fep.fitEqInts('hmine', p=p, dalpha=dalpha)
	dfile='data/productionfits/distfits-rlf%d-q%d.dat' % (int(100*rlenfac), int(100*q))
	fout=open(dfile, 'w')
	fout.write('#ETAS distance fits, q=%f, rlenfac=%f\n' % (q, rlenfac))
	fout.write('#eq, mstar, mc, lrlen1, lrlen2, rlenfac, q, estRateExp, Ninsides (N/L_r), <D>\n')
	fout.close()
	#
	fnum=0
	for eq in eqDict:
		# eq will be the key...
		mstar = eqDict[eq][1]
		mc    = eqDict[eq][2]
		print "Making %s plot, fignums (%d, %d)" % (eq, fnum, fnum+6)
		#plt.figure()
		#plt.clf()
		z=fitEqDists(eqset=eq, nbins=nbins, fnum=fnum, catnum=eqDict[eq][3], q=q, rlenfac=rlenfac, D=D)
		# returns: return [mstar, thismc, rlenfac, q, estRateExp, [Ninsides]]
		#
		# write some data to files:
		# z like: return [mstar, mc, rlenfac, q, estRateExp, Ninsides]
		lruplen1=z[0]/2.0 - 1.76
		lruplen2=lruplen1
		mt=7.6
		if z[0]>mt:
			lruplen1 = 1.5*z[0] - 9.36
			lruplen2 = mt/2.0 - 1.76		# = 1.5*mt - 9.36
		#
		# (m, mc, lr, lr, l-factor, q, log(N') estimated)
		print (z[0], z[1], lruplen1, lruplen2, z[2], z[3], z[4])
		outstr='%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t' % (eq, z[0], z[1], lruplen1, lruplen2, z[2], z[3], z[4])
		for n in z[5]:
			outstr+='%f\t' % n
		#
		# now, estimate the fractal dimension D:
		#meanN = math.log10(.5*(z[5][1]+z[5][2]))	# Algebraic mean N from the spherical and geographicLib methods. and of course, this is mood.
		#meanN = .5*(math.log10(z[5][1]*z[5][2]))	# quasi-Geometric mean						# because these yield almost identical results now.
		meanN = math.log10(z[5][2])	# either is fine.
		#print "solve D for meanN = ", meanN
		Ds = mcDsolver(meanN, min(mt,mstar), mc, dm=1.0, Drange=[0.0, 2.0], nits=10**5)
		print "Approx. Fract D: %f/%f" % (Ds[0], Ds[1])
		print "meanN_estimated: %f, meanN_obs: %f" % (z[4], meanN)
		outstr+='%f\t%f\t' % (Ds[0], Ds[1])
		#
		outstr=outstr[0:-1]+'\n'
		fout = open(dfile, 'a')
		fout.write(outstr)
		fout.close()
		#
		for fn in [fnum+6, fnum]:
			plt.figure(fn)
			if eq=='sansim':
				#plt.legend(title='%s' % eqDict[eq][0], loc='best', scatterpoints=1)
				#plt.legend(loc='best', scatterpoints=1)
				plt.legend(loc='upper left', scatterpoints=1)
			#else:
			plt.title('%s' % eqDict[eq][0])
			#
			if dosave: 
				if fn==fnum: plt.savefig('figs/distfits/%s-OmoriDistsBinned-rlf%d-q%d.png' % (eq, int(100*rlenfac), int(100*q)))
				if fn!=fnum: plt.savefig('figs/distfits/%s-OmoriDistsBinnedLL-rlf%d-q%d.png' % (eq, int(100*rlenfac), int(100*q)))
		fnum+=1

def makeDistFigsSpecial():
	print "making tohoku, q=2.0, D=1.0"
	z=fep.fitEqDists(eqset='tohoku', nbins=100, fnum=15, catnum=0, q=2.0, rlenfac=.5, D=1.0)
	lN=math.log10(z[5][2])
	Ds = mcDsolver(lN, min(7.6, 9.0), mc, dm=1.0, Drange=[0.0, 2.0], nits=10**5)
	print (z[0], z[1], z[2], z[3], z[4], lN)
	#
	print "making tohoku, q=2.0, D=1.5"
	A=fep.fitEqDists(eqset='tohoku', nbins=100, fnum=16, catnum=0, q=2.0, rlenfac=.5, D=1.5)
	lN=math.log10(z[5][2])
	Ds = mcDsolver(lN, min(7.6, 9.0), mc, dm=1.0, Drange=[0.0, 2.0], nits=10**5)
	print (z[0], z[1], z[2], z[3], z[4], lN)
	#
	for i in xrange(15, 16,21, 22):
		plt.figure(i)
		plt.title("Tohoku-oki")
	#
	print "making parkfield, rlenfac=1.0"
	A=fep.fitEqDists(eqset='parkfield', nbins=100, fnum=15, catnum=0, rlenfac=1.0, D=1.5)
	lN=math.log10(z[5][2])
	Ds = mcDsolver(lN, min(7.6, 5.96), mc, dm=1.0, Drange=[0.0, 2.0], nits=10**5)
	print (z[0], z[1], z[2], z[3], z[4], lN)
	for i in xrange(17,23):
		plt.figure(i)
		plt.title("Parkfield")
	#
	
# production Omori figures (Figs. 1,2,3 (i think))
def makeOmoriFigses(ps=[1.01, 1.05, 1.08, 1.1, 1.2, 1.25], dalpha=2.0, forcep=False, doMCfit=False, mcnits=10**5):
	#while ps[0]<ps[1]:
	tratios=[]
	f1=plt.figure(num=0, figsize=(6,5))
	f2=plt.figure(num=3, figsize=(6,5))
	for p in ps:
		fitfile = 'data/productionfits/intervalfits-p%d.dat' % (int(100*p))
		tr=makeOmoriFigs(p=p, dalpha=dalpha, forcep=True, doMCfit=doMCfit, datafile=fitfile)
		#ps[0]+=ps[2]
		tratios+=[tr]
	plt.figure(3)
	plt.clf()
	plt.figure(0)
	plt.clf()
	
	ax=plt.gca()
	ax.set_xscale('linear')
	ax.set_yscale('log')
	for rw in tratios:
		print rw
		plt.figure(0)
		plt.plot(map(operator.itemgetter(0), rw), map(operator.itemgetter(1), rw), 'o-')
		plt.figure(3)
		plt.plot(map(operator.itemgetter(2), rw), map(operator.itemgetter(1), rw), 'o')

def makeOmoriFigs(p=1.05, dalpha=2.0, fnum=0, forcep=False, datafile='data/productionfits/intervalfits.dat', doMCfit=False, mcnits=10**5):	
	eqDict = {'hmine': ['Hector Mine', 7.1, 3.0], 'parkfield':['Parkfield', 5.97, 1.5], 'sansim':['San Simeon', 6.5, 3.0], 'coalinga':['Coalinga', 6.7, 3.5], 'tohoku':['Tohoku-oki', 9.0, 4.5], 'elmayor':['El Mayor-Cucapah', 7.2, 2.5]}
	#
	#a=fep.fitEqInts('hmine', p=p, dalpha=dalpha)
	# for each eq. collect t0_fit/t0_rupt
	#datafiles=datafile.split('.')
	#datafile = datafile[0]+('p%d' % (int(100*p))) + '.' + datafile[1]
	fout=open(datafile, 'w')
	fout.write('# production interval fits. First set from intervals, second from CDF\n')
	fout.write('#mstar, mc, p_fixed, t0(p_fixed), t0ssim(p=p_fixed, dmtau=0), p(CDF fit), dmtau(CDF fit)\n')
	fout.close()
	#
	tratios=[]
	for eq in eqDict:
		# eq will be the key...
		print "Making %s plot" % eq
		a=fitEqInts(eq, p=p, dalpha=dalpha, mc=eqDict[eq][2], forcep=forcep, doMCfit=doMCfit, mcnits=mcnits)
		fout = open(datafile, 'a')
		outstr=''
		for itm in a:
			outstr += '%f\t' % itm
		outstr=outstr[0:-1]+'\n'
		fout.write(outstr)
		fout.close()
		#
		# write a file. fitEqInts reutns: [cpf.mainshock[3], cpf.mc, p, ft0, t0ssim, plsqCum[0], plsqCum[1]]
		# where p is a fixed p, ft0 is the pdf fit val for t0 given p, and plsqCum[0,1] are the p, dmtau from the cum fits.
		#
		plt.figure(1)	# we know the invervals plot is fig 1...
		info = eqDict[eq]
		plt.legend(title='%s: $m=%.2f, m_c=%.2f$' % (info[0], info[1], info[2]), loc='upper left', scatterpoints=1, numpoints=1)
		#print "savename: " + 'figs/omorifits/%s-OmoriInts-p%s.png' % (eq, str(p).replace('.', ''))
		plt.savefig('figs/omorifits/%s-OmoriInts-p%s.png' % (eq, str(p).replace('.', '')))
		#
		plt.figure(2)
		plt.legend(title='%s:\n $m=%.2f, m_c=%.2f, p=%.2f$' % (info[0], info[1], info[2], p), loc='best', scatterpoints=1, numpoints=1)
		plt.savefig('figs/omoriCumFits/%s-OmoriCum-p%s.png' % (eq, str(p).replace('.', '')))
		#
		#
		tratios+=[[p, a[3]/(10**(info[1]/2. - 2.28)), a[3]]]	# so this is t0/(dt-rupt)
																		# note the plots will probably be in p-sequences
																		#
		#fnum+=1
		
	return tratios

def tahirETAS(rtype='ssimbump', gridsize=.05, eqtheta=40., eqeps=1.0/2.0, deltat0=5, deltat1=1, dmstar=1.0, usecatlen=None):
	# demo of Tahir2012 and Zalohar2013 type distributions.
	#
	# deltat: number of days between the mainshock and aftershocks. use this to spoof weighting the hazard fields of
	# 			 these events (aka, large deltat effectively weights the aftershocks over the mainshock.
	# let's create an El Mayor-Cucapah scenario...
	lat0 = 32.128
	lon0 = -115.303
	m0=7.2
	m1=m0-dmstar
	Lr=10.0**(m0/2.0 - 1.76)
	#
	lats=[lat0-1.5, lat0+1.5]
	lons=[lon0-1.5, lon0+1.5]
	#
	deltaTime = dtm.timedelta(days=deltat0)
	#
	# make proposed aftershock positions:
	rad2deg=math.pi*2.0/360.0
	dx1 = Lr*math.cos(rad2deg*eqtheta)/(111.1*math.cos(lat0*rad2deg))
	dy1 = Lr*math.sin(rad2deg*eqtheta)/111.1
	#
	lat1a=lat0+dy1*1.3
	lat1b=lat0-dy1*1.3
	lon1a=lon0-dy1*1.3
	lon1b=lon0+dy1*1.3
	#
	t0=dtm.datetime(2010,4,4, 15, 40, tzinfo=pytz.timezone('UTC'))	# mainshock time	(change this notation of we introduce omori prams).
	t1=t0 + deltaTime
	tfc = t1+dtm.timedelta(days=deltat1)	# a guess for proper asthetics...
	#
	mycat = [[t0, lat0, lon0, m0]]
	mycat+= [[t1, lat1a, lon1a, m1]]
	mycat+= [[t1, lat1b, lon1b, m1]]
	
	b=bcp.BASScast(incat=mycat, fcdate=tfc, gridsize=gridsize, contres=7, mc=2.5, eqtheta=eqtheta, eqeps=eqeps, fitfactor=5.0, contour_intervals=None, lats=lats, lons=lons, doBASScast=False, rtype=rtype)
	# incat=[], fcdate=dtm.datetime.now(pytz.timezone('UTC')), gridsize=.1, contres=2, mc=2.0, eqtheta=0, eqeps=1.0, fitfactor=5.0, contour_intervals=None, lats=None, lons=None, doBASScast=True, rtype='ssim')
	#
	#myrtype='ssimbump'
	myrtype=rtype
	
	#b.rtype='omorisatbump'
	b.quakes[0].rtype=myrtype
	
	if len(b.quakes)>1:
		b.quakes[1].eqeps=1.0
		b.quakes[1].rtype='ssim'
	if len(b.quakes)>2:
		b.quakes[2].eqeps=1.0
		b.quakes[2].rtype='ssim'
	#
	while (usecatlen!=None and len(b.quakes)<=usecatlen):
		b.quakes.pop()
		#b.quakes.pop()
	#
	b.calcBASScast()
	#b.BASScastContourMap()
	#
	b.BASScastContourMap(maxNquakes=0)
	#
	#
	# now, let's do a 3d plot for fun:
	fig3d=plt.figure(2) 
	plt.clf()
	ax1 = fig3d.gca(projection='3d')
	
	Z=[]
	X=[]
	Y=[]
	ix=0
	iy=0
	for rw in b.Z2d:
		ix=0
		for elem in rw:
			Z+=[elem]
			X+=[b.X_i[ix]]
			Y+=[b.Y_i[iy]]
			ix+=1
		iy+=1
	#
	ax1.plot(X,Y,Z, 'b.-')
	#
	Z2=[]
	X2=[]
	Y2=[]
	ix=0
	iy=0
	Z2d2=b.Z2d.copy()
	Z2d2.transpose()
	for rw in b.Z2d:
		iy=0
		for elem in rw:
			Z2+=[elem]
			X2+=[b.X_i[ix]]
			Y2+=[b.Y_i[iy]]
			iy+=1
		ix+=1
	#	
	ax1.plot(X2,Y2,Z2, 'g-')
	#
	fig3=plt.figure(3)
	plt.clf()
	ax3 = fig3.gca(projection='3d')
	
	mycm = plt.get_cmap('spectral')
	cNorm  = colors.Normalize(vmin=min(Z), vmax=max(Z))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=mycm)
	print "len(Z): ", len(Z)
	
	zlen=len(Z)
	for i in xrange(zlen):
		colorVal = scalarMap.to_rgba(Z[i])
		#retLine, = ax.plot(line, color=colorVal)
		ax3.plot([X[i]], [Y[i]], [Z[i]], '.', color=colorVal)
		#print "arow: ", i, X[i], Y[i], Z[i], colorVal
		if i>0:
			z2val=scalarMap.to_rgba(.5*(Z[i]+Z[i-1]))
			ax3.plot([X[i-1], X[i]], [Y[i-1], Y[i]], [Z[i-1], Z[i]], '-', color=z2val, alpha=.6)
	for i in xrange(1,zlen):
		#colorVal = scalarMap.to_rgba(Z[i])
		#retLine, = ax.plot(line, color=colorVal)
		#ax3.plot([X[i]], [Y[i]], [Z[i]], '.', color=colorVal)
		#print "arow: ", i, X[i], Y[i], Z[i], colorVal
		z2val=scalarMap.to_rgba(.5*(Z2[i]+Z2[i-1]))
		ax3.plot([X2[i-1], X2[i]], [Y2[i-1], Y2[i]], [Z2[i-1], Z2[i]], '-', color=z2val, alpha=.6)
 	
 	#ax3.plot(X,Y,Z, 'k--', alpha=.4)
 	#ax3.plot(X2,Y2,Z2, 'k--', alpha=.4)
	#
	# and annotate a bit:
	#x0, y0 = 
	
	return b
	
def lognorm(x, mu, sigma):
	denom1 = x*math.sqrt(2.0*math.pi)*sigma
	expon = ((math.log(x)-mu)**2.0)/(2.0*sigma*sigma)
	z= (1.0/denom1)*math.exp(-expon)
	#z = (1.0/(x*sigma*math.sqrt(2.0*math.pi)))*math.exp(-((math.log(x)-mu)**2.0)/2.0*sigma**2.0)
	return z

def lognorm2(x, m, v):
	mu = math.log(m*m/math.sqrt(v+m*m))
	sigma = math.sqrt(math.log(1.0+v/(m*m)))
	return lognorm(x, mu, sigma)

#def lognormstd(x, r0)

def lognorm1d(X=10, mu=0.0, sigma=1.0, N=1000, fnum=1, doclf=True):
	#z=(1.0/(x*math.sqrt(2.0*math.pi)*sigma))*math.exp((-(math.log(x)-mu)**2.0)/(2.0*sigma*sigma))
	dx=float(X)/float(N)
	Xs=[0.0]
	Ys=[0.0]
	for i in xrange(N):
		Xs+=[Xs[-1]+dx]
		Ys+=[lognorm(Xs[-1], mu, sigma)]
	#
	plt.figure(fnum)
	if doclf: plt.clf()
	plt.ion()
	ax=plt.gca()
	ax.set_xscale('linear')
	ax.set_yscale('linear')
	plt.plot(Xs, Ys, '-')

#def lognorm1d2(X=10, mag=6.0, fnum=1, doclf=True):
#	r0=10.0**(mag/2.0 - 1.76)
#	mu1=.5*( 1.0 + 1.0-4.0*math.log(10.)*(mag/2.0 - 1.76))
#	mu2=.5*( 1.0 - 1.0-4.0*math.log(10.)*(mag/2.0 - 1.76))
	
	
	

def lognormtest():
		# make a lattice and contour-plot it.
		mu=0.0
		sigma=1.0
		N=100.
		Z=[]
		x0, y0=0.0, 0.0
		rangefact=.05
		#
		for i in scipy.arange(N+1):
			Z+=[[]]
			y=(i-N/2.0)*rangefact
			for j in scipy.arange(N+1):
				x=(j-N/2.0)*rangefact
				r=math.sqrt((y-y0)**2.0 + (x-x0)**2.0)
				#
				if r==0.0:
					Z[-1]+=[None]
					continue
				#	
				z=lognorm(r, mu, sigma)
				Z[-1]+=[z]
			#
		#
		plt.figure(0)
		plt.ion()
		plt.clf()
		plt.contourf(Z, 15)
		plt.colorbar()
	
	
	
	
	
	
	
		
