import yodapy as ypp
import BASScast as bcp
import ANSStools as atp
import eqcatalog as eqp
import linefit as lft

import datetime as dtm
import pytz
import operator
import math
import random
import numpy
import scipy
import os
import glob
import multiprocessing as mcp
#
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mpd
import scipy.optimize as spo
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

#alpha0=1.19
# (note, we might have been using alpha_0=1.28 for some development bits.)
alpha_0 = 2.28	#log(dt) ~ (m/2) - alpha0; 2.28 < alpha0 < 2.57 ? (from Vidale and Houston [1993]) # where does 1.28 come from?
#alpha_0 = 2.28 - 1.0	# where this notation is to show that we're doing an alpha0 - DeltaAlpha0 type scaling...
							# this is similar to the notion that Omori scaling starts approximately 10*rupture-duration after mainshock ends.
							# or maybe we make the distinction between rupture time and t0...
							# and alpha_0 -> 2.6 might be closer to correct (as per SCEC poster, v_rupt = 1.5 km/sec, log(v).
							# because: dt = L/v; m/2 - alpha = m/2 - l0 - nu, (v=log(nu)). so, alpha = nu - l0 and v~1.5 km/sec.
globrand=random.Random()
mrand=random.Random()
days2secs=3600.0*24.0
plt.ion()

# interval and distance testing scripts:
def eqintshell(alpha0=None, omoriWaitExp=None, mc=1.5, drho=5.0, p=1.1, eqname='parkfield', dorefresh=False):
	# shell for earthquake intervals.
	# defaults:
	catnum=0
	p=1.1
	alpha0=None
	omoriWaitExp=None
	if eqname=='parkfield':
		fname='cats/parkfield.cat'
		#alpha0=None
		#omoriWaitExp=None
		mc=1.5
		drho=5.5
		#p=1.1
		catnum=1	# deviating from default
		cutDtm = dtm.datetime(2004,9,27, tzinfo=pytz.timezone('UTC'))
		mstar=6.0
	#
	# other earthquakes...
	#
	# main bits:
	cpf=ypp.eqcatalog()
	cpf.loadCatFromFile(fname=fname, minmag=mc)
	if alpha0==None:
		alpha0=alpha_0
	
	mev=cpf.getMainEvent(cpf.getcat(catnum))
	
	#while cpf.cat[0][0]<dtm.datetime(2004,9,27, 0, 0, 0, 0, pytz.timezone('UTC')): cpf.cat.pop(0)
	while cpf.cat[0][0]<mev[0]: cpf.cat.pop(0)
	
	# eq specific catalogs, etc.
	if eqname=='parkfield':
		cpf.addEllipCat('PFshock (.4 x .15)', cpf.cat, 40.0, 35.9, -120.5, 0.8, 0.3)
	#
	if mstar==None: mstar=mev[3] #?
	
	return eqintervals(cat=cpf, mstar=mstar, mc=mc, mev=mev, catnum=catnum, alpha0=alpha0, omoriWaitExp=omoriWaitExp, drho=drho)

def eqdistances(c1=None, mc=2.0, catnum=0, mev=None, maxlen=None, cantame=None):
	# get earthquake distance distribution. if mev is none, getMainEvent() from the catalog; dump everything before MEV
	#
	if c1==None:
		c1=eqp.eqcatalog([])
		c1.loadCatFromFile(catname)
	#
	while len(c1.subcats)<catnum:
		catnum-=1
	#
	if maxlen==None: maxlen=len(c1.getcat(0))
	#
	if mev==None:
		mev=c1.getMainEvent()
		#print mev
	mycat=[]
	# keep only earthquakes after the mainshock
	for rw in c1.getcat(catnum):
		if rw[0]<=mev[0]: continue
		if rw[3]>=mc: mycat+=[rw]
		if len(mycat)>=maxlen: break
	#print "mycat: ", mycat[0:4]
	#
	d1=c1.getDistances(catList=mycat, center=[mev[2], mev[1]])
	#
	return [d1, mev]
#
def pfdisttest(catname='cats/parkfield.cat', mc=1.5, mev=None, len0=None, q=1.34, lenfact=1.5, Nfact=0.8):
	# a special script because the parkfield cat probably includes foreign earthquakes.
	c1=eqp.eqcatalog([])
	c1.loadCatFromFile(catname)
	c1.addEllipCat('PFshock (.4 x .15)', c1.cat, 40.0, 35.9, -120.5, 0.8, 0.3)
	thismev=c1.getMainEvent(thiscat=c1.getcat(1))
	#print "pfmev: ", thismev
	while c1.cat[0][0]<thismev[0]: c1.cat.pop(0)
	#
	return eqdistTest(c1=c1, mc=mc, mev=thismev, catnum=1, len0=len0, q=q, lenfact=lenfact, Nfact=Nfact)
#
#def fomorilike(a, lx, lx0, p):
def fomorilike(X, tau, x0, p):
	#lx = scipy.array(lx)
	X=scipy.array(X)
	#thisNth = Nfact*((r0 + r2)**(1.0-q) - r0**(1.0-q)) / (chi * (1.0-q))
	
	Yout = ((X+x0)**(1.0-p) - x0**(1.0-p))/(tau*(1.0-p))
	#Yout = (A/(1.0-p))*((X+X0)**(1.0-p) - X0**(1.0-p))
	#Yout = ((10**a)/(1.0-p)) * ((10**lx + 10**lx0)**(1.0-p) - 10**(lx0*(1.0-p)))
	return Yout
#
def MComoriFit(X, Y, A0=None, x0=None, p0=None, nits=10*5):
	# brute-foce fit to omori type function:
	# X, Y should be scipy-arrays. it may be necessary to loop through the lists semi-manually.
	# assume X arrives in linear form.
	#a0=math.log10(A0)
	#lx00 = math.log10(x0)
	#nits=int(10**5)
	p0=1.001
	R1 = random.Random()
	R2 = random.Random()
	R3 = random.Random()
	#
	#minprams=[a, lx0, p]
	minprams=None
	#
	lX = []
	lY = []
	for x in X:
		lX+=[math.log10(x)]
	for y in Y:
		lY += [math.log10(y)]
	#
	lxreff = math.log10(.75*max(X))
	la0reff = math.log10(max(Y))
	#
	#Ytest = fomorilike(a0, lX, lx0, p)
	#minerr = sum((Ytest-Y)**2.0)
	minerr=None
	#
	plt.figure(7)
	plt.clf()
	plt.ion()
	for i in range(nits):
		# do a mc fit...
		#dummyvar=0.0
		p=1.25 + .75*R1.random()
		#lx0 = lxreff*R2.random()
		lx0 = -2.0 + 3.0*R2.random() 
		#a = la0reff*R3.random()
		ltau = -6.0 + 6.0*R3.random()
		#a = 0.0 + 4.0*R3.random()
		#
		Xprime = scipy.array(X, copy=True)-10**lx0
		Ytest = fomorilike(Xprime, 10**ltau, 10**lx0, p)
		thiserr = 0.
		NerrCount=0
		for i in range(len(Y)):
			#print Y[i], Ytest[i]
			if Y[i]<=0 or Ytest[i]<=0:
				thiserr += 0
				continue
			thiserr += (math.log10(Y[i]) - math.log10(Ytest[i]))**2.0 
			NerrCount+=1
		thiserr=thiserr/float(NerrCount)
		#
		if minerr==None or thiserr<minerr:
			minerr = thiserr
			minprams=[ltau, lx0, p]
			print(minerr, minprams)
			#plt.loglog(X, fomorilike(X, 10.0**minprams[0], 10.0**minprams[1], minprams[2]), '-') 
		#
	plt.figure(7)
	plt.clf()
	plt.loglog(X, Y, '.')
	#	
	#plt.loglog(Xprime + 10**minprams[1], fomorilike(X, 10.0**minprams[0], 10.0**minprams[1], minprams[2]), '-', label='fomoriX2')
	#plt.loglog(X, fomorilike(X, 10.0**minprams[0], 10.0**minprams[1], minprams[2]), '-', label='fomoriX2')
	plt.loglog(X, fomorilike(X-10**minprams[1], 10.0**minprams[0], 10.0**minprams[1], minprams[2]), '-', label='fomoriX3')
	#plt.loglog(Xprime, fomorilike(X, 10.0**minprams[0], 10.0**minprams[1], minprams[2]), '-', label='fomoriX')
	#plt.loglog(X, fomorilike(X, 10.0**minprams[0], 10.0**minprams[1], minprams[2]), '-', label='fomoriX')
	plt.legend(loc='best')
	#Y2 = fomorilike(10**2.16, X, .99, 1.1)
	#plt.loglog(X, Y2, '.-.')
	#
	print(minprams)
	print("lt0: %f, ltau: %f, p: %f" % (minprams[1], minprams[0], minprams[2]))
	#
	return minprams
#
def calcEqDists(cat=None, mstar=None, mc=None, b=1.0, C1=2.77, C2=3.77, q=1.5, skipexp=3.0, lw=2, mev=None, Nfact=1.0, fnum=1, titlename='distances', mt=7.6, len0=None, catnum=0):
	# ralpha=.5, rbeta=0.0, K1=3., K2=.1, K3=1.0
	
	#
	# this is an experimental eqDists version, consistent with present write-up (Dec. 2012).
	# default cat for dev:
	#if cat==None:
	#	cat=eqp.eqcatalog([])
	#	cat.loadCatFromFile('cats/sansim-mfetas.cat', minmag=3.0)
	#	mc=3.0
	#
	# note: dmstar=1.0 gives full sequence type behavior
	#       dmstar=1.2 gives "one earthquake" behavior (aka, sum(all events with dmstar=1.2) -> dmstar=1.0
	#if alpha0==None: alpha0=alpha_0
	#alpha=alpha0
	#b=1.0
	#
	cpf=cat # catalog object
	if type(cpf).__name__ == 'list':
		#cpf=ypp.eqcatalog(cpf)
		eqp=ypp.eqcatalog(cpf)
		cpf=eqp
	if len(cpf.subcats)==0: catnum=0
	if catnum>len(cpf.subcats): catnum=len(cpf.subcats)
	if mev==None: mev=cat.getMainEvent(cat.getcat(catnum))
	#
	if mstar==None:
		# assume largest mag in catalog:
		mstar=cpf.getMainEvent(cpf.getcat(catnum))[3]
		print("mstar chosen: ", mstar)
	#
	lNom =  b*(mstar - 1.0 - mc)
	Nom = 10**lNom
	#	
	cpf.plotGRdist(fignum=0, doShow=True, fname=None)
	#eqints=cpf.getIntervals(catList=cpf.getcat(catnum), winLen=1)	# intervals are in days, right?
	#eqints.reverse()
	#D1 = eqdistances(c1=cpf, mc=mc, mev=mev, catnum=catnum, maxlen=round(10**lNom))
	D1 = eqdistances(c1=cpf, mc=mc, mev=mev, catnum=catnum)
	#
	# eqdistances returns: [dists, mev]
	# dists (from eqcatalog.getDistances()): [ [rw[0], thisx, thisy, R0, R1, R2, R3, R4, depth], .. [] ]
	distances=D1[0]		#
	mev=D1[1]
	#
	# be sure distances is date-sorted:
	distances.sort(key=lambda x:x[0])
	while distances[0][0]<=(mev[0]):
		distances.pop(0)
	# cat. length is important. we want only aftershocks, so select by length N
	if len0==None:
		len0=int(10**(mev[3]-mc-1.0))
	if len0==-1:
		len0=len(distances)
	if len0>len(distances): len0=len(distances)
	#
	#skipIntervals=int(10**(mstar-mc-skipexp))
	skipIntervals=int(10**(mstar-mc-skipexp))
	#
	Rs=list(map(operator.itemgetter(6), distances))
	Rskip = Rs[skipIntervals:]
	Rs.sort()	# now, in distance-order
	Rskip.sort()
	#Rskip = Rs[skipIntervals:len0]	# skip in time or in r? how 'bout skip the whole thing for now...
	#
	R3D=[]
	for rw in distances:
		#print rw
		#R3D += [math.sqrt(rw[6]*rw[6] + rw[8]*rw[8])]	# rw[3]-rw[7] are distances; rw[8] is depth.
		R3D += [math.sqrt(rw[3]*rw[3] + rw[8]*rw[8])]	# rw[3]-rw[7] are distances; rw[8] is depth.
	R3D.sort()
	#RskipR = Rs[skipIntervals:]
	RskipR = R3D[skipIntervals:]
	#r02d = 10**(mstar/2.0 - C2)	# define second r0 below...
	#
	RdR =  []
	Nrdr = []
	#
	Ns=list(range(1, len0+1))
	
	Nactual=[]
	Ntheory=[]
	NexpTh=[]
	NTPL=[]
	NPLexp=[]
	#
	lr0 = mstar/2.0 - C1
	lr2 = mstar/2.0 - C2
	if mstar>mt:
		lr0 = mstar - C1
		lr2 = mstar - C2
	# though 3*mstar/2 might actually be correct, in the sense that nominally dN/dr ~ L_r.
	r0=10**lr0
	r02 = 10**lr2
	rtrunc = 10**(lr2 + 1.5)
	#	
	NactSkip=[0]
	rskip=[0]
		#
	# this is a genaral solution for chi (from GR normalization) given an lr0:
	if mstar<=mt:
		#ltau = mc*(p+1.)/2.0 - mstar -math.log10(p-1.) + C2
		lchi = (1.-q)*lr0 - math.log10(q-1.) - lNom	# and note that this gives the (q+1) behavior above.
		lchi2 = (1.-q)*lr2 - math.log10(q-1.) - lNom
	else:
		lr0 = mstar - C1
		lr2 = mstar - C2
		lchi = (1.-q)*lr0 - 1.0*(mt-mc) - 1.5*(mstar-mt) - math.log(q-1.0)	# where we have assumed b1=1.0, b2=1.5
		lchi2 = (1.-q)*lr2 - 1.0*(mt-mc) - 1.5*(mstar-mt) - math.log(q-1.0)		# it might also be necessary to average over
																										# the large/small mag. domain (to get a more
																										# accurate result).
	chi = 10**lchi
	chi2 = 10**lchi2
	
	#Atpl = (Nom/(r02**(-q) - rtrunc**(-q)))
	lr0tpl = q*(mstar/2.0 - C2 + .5)
	Atpl = (Nom*(10**lr0tpl)) - rtrunc**(-q)       # truncated power law
	Aexp = Nom*(10**lr0tpl)*math.exp(r02/rtrunc)   # exponential funct.
	rs=[0.]
	print("mstar=%.2f, Lr=%f, r0=%f, r2=%f, chi=%f, chi2=%f, Ngr=%f:" % (mstar, 10**(mstar/2.0-2.77), r0, r02, chi, chi2, 10**lNom))
	#
	#
	Ntheory=[]
	Ntheory2=[]
	NTPL=[]
	#Ntheory3=[0]	# "local rate" method.
	#Ntheory4=[0]
	RsTheory=[0.]
	#
	logmaxR=math.log10(Rs[-1])
	#
	# note: mean rate for san-sim catalog approx 15.52 eq>3/yr, 4.92*10**-7 eq>3/sec
	lRtheory=-2.0	# starting log-time
	#lttheory = lt0+6.
	dRtheory=(logmaxR - lRtheory)/10000	# Dist interval.
	thisindex=0
	#
	#rmc=10**(mc/2.0 - C2)
	#
	bgdensity = 0.0*.001
	while lRtheory<(logmaxR + .25):
		r1=10**(lRtheory)
		r2=10**(lRtheory+dRtheory)
		#
		#Ntheory  += [Nfact*((r0 + r2)**(1.0-q) - r0**(1.0-q)) / (chi * (1.0-q)) + bgdensity*math.pi*(r0+r2)**2.]
		#Ntheory2 += [Nfact*((r02 + r2)**(1.0-q) - r02**(1.0-q)) / (chi2 * (1.0-q))  + bgdensity*math.pi*(r0+r2)**2.]
		Ntheory  += [Nfact*(Nom - (((r0 + r2)**(1.0-q) - r0**(1.0-q)) / (chi * (1.0-q)) + bgdensity*math.pi*(r0+r2)**2.))]
		Ntheory2 += [Nfact*(Nom - (((r02 + r2)**(1.0-q) - r02**(1.0-q)) / (chi2 * (1.0-q))  + bgdensity*math.pi*(r0+r2)**2.))]
		#
		#		
		#thisNexp = Ntheory2[-1]*(math.exp(-(r2-r02)/(10.*q)))
		#if thisNexp>.5: NexpTh += [thisNexp]
		#
		NTPL+=[Atpl*((r2 + 2.*r02)**(-q) - rtrunc**(-q))]
		#NTPL+=[Atpl*((r2 + rtrunc)**(-q) - rtrunc**(-q))]
		#NPLexp+=[Nom*((r02/r2)**q)*math.exp((r02-r2)/rtrunc)]
		NPLexp+=[Aexp*((r2+2.0*r02)**(-q))*math.exp(-r2/rtrunc)]
		#NTPL += [(Nom/ (chi2 * (1.0-q)))*((r02+r2)**(1.0-q) - (r2 + 30.*r02)**(1.0-q))]
		#if (r2-r02)>0: 
		#NTPL += [Nom*((r2)**(1.0-q) - (100.)**(1.0-q))]
		#
		'''
		if r2<r0:
			RdR +=  [r2]
			Nrdr += [1.]
		if r2>=r0:
			RdR += [r2]
			Nrdr += [Nrdr[-1] + Nfact*((r0 + RdR[-1])**(1.0-q) - (r0+ RdR[-2])**(1.0-q))/(chi*(1.0-q))]
		'''
		#
		RsTheory += [RsTheory[-1] + (r2-r1)]
		thisindex+=1
		lRtheory+=dRtheory
	RsTheory.pop(0)
	#
	#dfits=MComoriFit(X=scipy.array(Rs), Y=scipy.array(range(1, len(Rs)+1)), nits=10**5)
	dfits=[0., 0., 0.]
	#print "MC fits imply: lr0=%f, lchi=%f, p=%f, C1=(dunno), C2=%f" % (dfits[1], dfits[0], dfits[2], (mstar/2.0-dfits[1]))
	#
	#plfunct=plt.plot
	plfunct=plt.loglog
	#
	plt.figure(fnum)
	plt.clf()
	#
	# the following plots are in "scaling" space, rprime = r-r0 (aka, dist. from rupture-horizon).
	#plfunct(scipy.array(R3D)-r02, range(1, len(R3D)+1), 'b.-.', label='nAct3D rprime',lw=lw, zorder=4)
	#plfunct(scipy.array(RsTheory)+r02*0., Ntheory2, 'b-.', label='theory2', alpha=.7, lw=2*lw, zorder=3)
	#
	# here's the 3D-dist plot in "actual" space, r = rprime + r0
	Rs=scipy.array(Rs)
	RsTheory = scipy.array(RsTheory)
	R3D=scipy.array(R3D)
	Ntheory=scipy.array(Ntheory)
	Ntheory2 = scipy.array(Ntheory2)
	#
	# note sorting syntax: X[::-1] --> " by negative 1's", so sort backwards. every-other backwards is X[::-2].
	plfunct(Rs, list(range(1, len(Rs)+1))[::-1], 'r.', label='nActual r', lw=lw, zorder=4)
	plfunct(RsTheory + r0, Ntheory, 'r-.', label='theory', alpha=.7, lw=2*lw, zorder=3)

	plfunct(R3D, list(range(1, len(R3D)+1))[::-1], 'g.', label='nAct3D, r',lw=lw, zorder=4)
	plfunct(RsTheory + r02, Ntheory2, 'g-.', label='theory2, r', alpha=.7, lw=2*lw, zorder=3)
	
	#plfunct((Rs + R3D)/2.0,range(1, len(R3D)+1)[::-1], 'c--')
	
	plfunct(RsTheory[:len(NexpTh)], NexpTh, 'y--', label='exp')
	plfunct(RsTheory[:len(NTPL)], NTPL, 'c--', label='TPL')
	plfunct(RsTheory[:len(NPLexp)], NPLexp, 'k--', label='PLexp')
	#
	#print "lens: %d, %d" % (len(RdR), len(Nrdr))
	#plfunct(RdR, Nrdr, '.-.', label='theory2', alpha=1.0, lw=2)
	#plfunct(RdR, range(1, len(RdR)+1), '-.', label='theory-RdR', lw=2.0*lw, zorder=4)
	#plt.plot(scipy.array(tsTheory), scipy.array(Ntheory2), '--', label='theory2', lw=2*lw, zorder=3)
	plt.xlabel('distance $r$')
	plt.ylabel('Number $N(<r)$')
	#strtitle = 
	plt.title('%s Distances, $m_c=%.2f$, $C_1 = %.2f$, $C_2 = %.2f$ \n $lr0=%f$, $lr2d=%f$, $l\\chi=%f$, $q=%.2f$, $N_{fact}=%.2f$' % (titlename, mc, C1, C2, lr0, math.log10(r02), math.log10(chi), q, Nfact))
	plt.legend(loc='best')
	#
	# transposed:
	plt.figure(2)
	plt.clf()
	print("nths:", Ntheory[0:3], Ntheory2[0:3])
	plfunct(Rs, list(range(1, len(Rs)+1)), 'r.-.', label='nActual r', lw=lw, zorder=4)
	plfunct(RsTheory + r0, Nom -Ntheory/Nfact, 'r-.', label='theory', alpha=.7, lw=2*lw, zorder=3)

	plfunct(R3D, list(range(1, len(R3D)+1)), 'g.-.', label='nAct3D, r',lw=lw, zorder=4)
	plfunct(RsTheory + r02, -((Ntheory2/Nfact)-Nfact*Nom), 'g-.', label='theory2, r', alpha=.7, lw=2*lw, zorder=3)
	#
	return distances
#
def eqdistTestF(catname='cats/parkfield.cat', mc=2.0, mev=None, len0=None, catnum=0, q=1.34, lenfact=3.07, Nfact=.8, maxdt=None, lw=1, C1=0.0, C2=0.0, K1=1.0, K2=0.0):
	# note:
	# we originally parameterize lrupt at m/2 - 2.77, so these 3.0 numbers are pretty interesting.
	# we might just fix at 2.77 (or so) and see what we can do with the other prams as well.
	# a=ybp.eqdistTestF(catname='cats/sansim-mfetas.cat', mc=2.5, q=1.5, lenfact=3.0, Nfact=.61)
	# a=ybp.eqdistTestF(catname='cats/sansim-mfetas.cat', mc=3.0, q=1.65, lenfact=3.0, Nfact=1.95)
	#  a=ybp.eqdistTestF(catname='cats/pfshock-20120109.cat', mc=1.5, q=1.25, lenfact=2.59, Nfact=.9)
	#
	c1=eqp.eqcatalog([])
	c1.loadCatFromFile(catname, minmag=mc)
	#
	return eqdistTest(c1=c1, mc=mc, mev=mev, len0=len0, catnum=catnum, q=q, lenfact=lenfact, Nfact=Nfact, maxdt=maxdt, lw=lw, C1=C1, C2=C2, K1=K1, K2=K2)

def eqdistTest(c1=None, mc=2.0, mev=None, len0=None, catnum=0, q=1.34, lenfact=3.07, Nfact=.8, maxdt=None, lw=1, C1=0.0, C2=0.0, K1=1.0, K2=0.0):
	# len0 is a catalog length (aka, for skipping -- skip m-mc-2.0 early events)
	# c1 is a catalog...
	if lenfact==None:
		# see: lruptLen(m, dl0=1.0)
		# this gives r=L_r/2, so aftershocks and scaling initiate outside the rupture zone.
		# note, this is from kagan 2002 (1.77), dL=1.0 (yoder et al. basic SS scaling theory), then r=D/2 = L_r/2
		#lenfact = 1.77 + 1.0 + math.log10(2.0)
		lenfact = 3.07
	# Nfact: this is just a numerical efficiency factor. 1) we might not have seen or be representing all
	# the aftershocks (aka, recent earthquake). 2) there is statistical variation in the actual numbers
	# of detectable aftershocks., 3) it may be necessary to mask out some aftershocks, so N is affected.
	if c1==None:
		c1=eqp.eqcatalog([])
		c1.loadCatFromFile('cats/sansim-mfetas.cat', minmag=mc)
	#
	# max deta t:
	if maxdt!=None:
		maxdt=int(maxdt)
		#
		maxdate=c1.getcat(catnum)[0][0]+dtm.timedelta(days=maxdt)
		#print maxdate, c1.cat[0][0]
		while c1.getcat(catnum)[-1][0]>(maxdate): c1.getcat(catnum).pop()
		print("end-date: ", c1.getcat(catnum)[-1][0])
		#
	#
	D1 = eqdistances(c1=c1, mc=mc, mev=mev, catnum=catnum)
	# eqdistances returns: [dists, mev]
	# dists comes from eqcatalog.getDistances()
	# rw[0] is the datetime, x,y, are the position prams, R_i are differend distance calculations.
	# dists += [[rw[0], thisx, thisy, R0, R1, R2, R3, R4]]
	d1=D1[0]		# eqdistances returns [distances, main-event]
	mev=D1[1]
	#
	if len0==None:
		len0=int(10**(mev[3]-mc-1.0))
	if len0==-1:
		len0=len(d1)
	if len0>len(d1): len0=len(d1)
	#
	# first, plot the actual distances. then, we'll compare to ETAS models.
	#
	# d1 is like [T, x, y, R0, R1, R2, R3], where R_i are in order of reliability.
	T=list(map(operator.itemgetter(0), d1))
	Tf = mpd.date2num(scipy.array(T))	# float-times, to evaluate scaling time.
	Tf=Tf-Tf[0]
	# full spherical expression:
	R=list(map(operator.itemgetter(6), d1))	# dist-type R3
	Ns=list(range(1, len0+1))
	Ns.reverse()
	thisR=R[:]
	thisR.sort()
	# get 3D distance (probably a more compact way to do this...)::
	R3D=[]
	for rw in d1:
		#R3D += [math.sqrt(rw[6]*rw[6] + rw[8]*rw[8])]	# rw[3]-rw[7] are distances; rw[8] is depth.
		R3D += [math.sqrt(rw[3]*rw[3] + rw[8]*rw[8])]	# rw[3]-rw[7] are distances; rw[8] is depth.
	R3D.sort()
	#
	R0 = list(map(operator.itemgetter(3), d1))
	R0.sort()
	# geographic-lib distances:
	Rgl = list(map(operator.itemgetter(7), d1))
	Rgl.sort()
	#
	Nth = []
	Nth1 = []
	Nth2 = []
	Nth3 = []
	#
	mstar=mev[3]
	b=1.0
	dmstar=1.0
	lNom=(b*(mstar-dmstar-mc))
	Nom = 10**lNom
	print("C1: ", C1)
	lr0 = mstar/2.0 - C1
	lchi = (1.0-q)*lr0 - lNom - math.log10(q-1.0) - C2
	r0 = 10**lr0
	chi = 10**lchi
	#print "lr0, lchi, lNom, q, r0+r-r0", lr0, lchi, math.log10(Nom), q
	for r in thisR:
		Nth+=[RDensity(mstar=mstar, r=r, mc=mc, q=q, lenfact=lenfact)]
		Nth2+=[Nfact*RDensity(mstar=mstar, r=r, mc=mc, q=q, lenfact=lenfact)]
		
		Nth3 += [Nom - (K2 +K1*((r0 + r)**(1.-q) - r0**(1.-q))/(chi*(1.0-q)))]
		#Nth2 += [Nfact * Nth[-1]]
		##Nth2 += [Nfact*(Nom - (r0 + r)**(1.-q) - r0**(1.-q))/(chi*(1.0-q))]
	#
	# and the other formulations check out pretty well
	#R2 = map(operator.itemgetter(5), d1)
	#R2.sort()	
	#
	plt.ion()
	#plt.figure()
	#plt.clf()
	#
	c1.plotCatMap(catalog=c1.cat[-len0:], fignum=0)
	X=list(map(operator.itemgetter(1), d1))
	Y=list(map(operator.itemgetter(2), d1))
	plt.figure(1)
	plt.clf()
	plt.plot(X,Y, '.')
	plt.plot([0.], [0.], '*r', ms=15)
	
	plt.figure(2)
	plt.clf()
	#plfunct=plt.plot
	plfunct=plt.loglog
	#plfunct=plt.semilogx
	#plfunct=plt.semilogy
	plfunct(thisR[0:len0],Ns, '-', lw=lw, label='dist1-propSph.', zorder=1)	# dist-type 1
	# -10.0**((mstar/2.)-lenfact)
	plfunct(scipy.array(R3D[0:len0])/1., Ns, 'c.-', lw=lw, label='dist3D')
	plfunct(R0[0:len0],Ns, '-', lw=lw, label='dist-RL')		# dist-type 2
	plfunct(Rgl[0:len0],Ns, '-', lw=lw, label='dist3-GeoLib')		# dist-type 3
	plt.xlabel('distance R')
	plt.ylabel('Number N(>R)')
	#plfunct(R2[0:len0],Ns, '--')
	#print len(thisR), len0, len(Nth)
	plfunct(thisR, Nth, 'b--', lw=2*lw, label='theory-raw', zorder=2)
	plfunct(thisR, Nth2, 'r--', lw=2*lw, label='theory-nfact', zorder=2)
	plfunct(scipy.array(thisR)+5., Nth3, 'g--', lw=2*lw, label='theory-omtype', zorder=2)
	plt.legend(loc='best', numpoints=1)
	#
	plt.figure(3)
	plt.clf()
	plt.plot(T, R, '.')
	plt.xlabel('float-time (days)')
	plt.ylabel('Distance $R$')
	#
	plt.figure(4)
	plt.clf()
	plt.loglog(Tf, R, '.')
	plt.xlabel('log(float-time[days])')
	plt.ylabel('Distance $R$')
	#
	#
	return d1

def rdensity(mstar, r, mc=2.0, dmstar=1.0, q=1.35, lenfact=1.5):
	r0 = 10.0**(.5*mstar - lenfact)
	Nom = 10**(mstar - dmstar - mc)
	#
	# handle the scaling break point.
	rprime=r-r0
	if rprime<0.0: rprime=0.0
	#
	lamb = ((1.0-q) * Nom * r0**(1.0-q))/((r0 + rprime)**q)
	
	return lamb
	
def RDensity(mstar, r, mc=2.0, dmstar=1.0, q=1.35, lenfact=1.5):
	r0 = 10.0**(.5*mstar - lenfact)
	Nom = 10**(mstar - dmstar - mc)
	#
	# handle the scaling break point.
	rprime=r-r0
	if rprime<0.0: rprime=0.0
	#
	#lamb = Nom * (r0**(1.0-q)) * (r0 + r)**(1.0-q)
	lamb = Nom * (r0**(1.0-q)) * (r0 + rprime)**(1.0-q)
	
	return lamb
	
#def pfintervals(alpha0=None, omoriWaitExp=None, mc=1.5, drho=5.0, p=1.1):
def pfintervals(catname='cats/parkfield-mfetas.cat', dorefresh=True, mc=1.5, mstar=None, C1=5.95, C2=0.0, p=1.1, skipexp=3.0, lw=2., mev=None, Nfact=1.0, fnum=1, titlename='Parkfield Intervals', useCatnum=1, ralpha=-5., rbeta=-.5, fitFract=1.0):
	cpf=ypp.eqcatalog()
	#cpf.loadCatFromFile(fname='cats/parkfield.cat', minmag=mc)
	cpf.loadCatFromFile(fname=catname, minmag=mc)
	#if alpha0==None:
	#	alpha0=alpha_0
	#
	#while cpf.cat[0][0]<dtm.datetime(2004,9,27, 0, 0, 0, 0, pytz.timezone('UTC')): cpf.cat.pop(0)
	#print "alpha0: %f" % alpha0
	# get the pf aftershocks:	
	cpf.addEllipCat('PFshock (.4 x .15)', cpf.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
	#cpf.addEllipCat('PFshock (.4 x .15)', cpf.cat, 40.0, 35.9, -120.5, 0.8, 0.3)
	#useCatnum=1
	mev=cpf.getMainEvent(cpf.getcat(useCatnum))
	print(mev)
	while cpf.getcat(useCatnum)[0][0]<=mev[0]: cpf.getcat(useCatnum).pop(0)
	mstar=6.0
	
	#return cpf
	
	#
	pfwaittime=None
	
	#return eqintervals(cat=cpf, mstar=mstar, mc=mc, mev=mev, catnum=useCatnum, alpha0=alpha0, omoriWaitExp=omoriWaitExp, drho=drho)
	return calcEqints(cat=cpf, mstar=mstar, mc=mc, C1=C1, C2=C2, p=p, skipexp=skipexp,lw=lw, catnum=useCatnum, ralpha=ralpha, rbeta=rbeta, fitFract=fitFract)

#def hmineintervals(alpha0=None, omoriWaitExp=None, catname='cats/hmine-hmetas.cat', dorefresh=False, useCatnum=0, mc=3.0, drho=5.0, p=1.1):
def hmineintervals(catname='cats/hmine-mfetas.cat', dorefresh=True, mc=3.0, mstar=None, C1=5.95, C2=0.0, p=1.1, skipexp=3.0, lw=2., mev=None, Nfact=1.0, fnum=1, titlename='Hmine Intervals', ralpha=-5., rbeta=-.5, fitFract=1.0):

	# try p=1.01
	# #lon=[-117.25, -115.5]	lat=[33.5, 35.5]	m0=3.000000	dates=[datetime.datetime(1980, 1, 1, 0, 0), datetime.date(2011, 5, 12)]
	
	if dorefresh==True:
		print("getting catalog")
		hmcat = bcp.atp.catfromANSS(lon=[-117.25, -115.5], lat=[33.5, 35.5], minMag=mc, dates0=[dtm.datetime(1980, 1,1, tzinfo=pytz.timezone('UCT')), dtm.datetime(2009,12,31, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC'))], Nmax=999999, fout=catname)
		cpf=ypp.eqcatalog(hmcat)
	else:
		cpf=ypp.eqcatalog()
		cpf.loadCatFromFile(fname=catname, minmag=mc)
	
	useCatnum=0
	
	#if alpha0==None:
	#	alpha0=alpha_0
	#
	#print "alpha0: %f" % alpha0
	# get the hmine aftershocks:
	# 67.4, 34.594, -116.271, .65, .2, datetime.datetime(1999, 10, 16, 02, 46, 45)
	cpf.addEllipCat('shock (.65 x .2)', cpf.cat, 67.4, 34.594, -116.271, 0.65, 0.2)
	mev=cpf.getMainEvent(cpf.getcat(useCatnum))
	#	
	print(mev)
	#mstar=7.1	# hmine is 7.1 according to ANSS; we might check with USGS source.
	#
	#myt0=10.0**(lDeltat(m=mstar, alpha=alpha0))
	#myt0=10.0**(lt0(m=mstar, alpha=alpha0))   
	
	#return eqintervals(cat=cpf, mstar=mstar, mc=mc, mev=mev, catnum=useCatnum, alpha0=alpha0, omoriWaitExp=omoriWaitExp, drho=drho, p=p)
	return calcEqints(cat=cpf, mstar=mstar, mc=mc, C1=C1, C2=C2, p=p, skipexp=skipexp,lw=lw, ralpha=ralpha, rbeta=rbeta, fitFract=fitFract)

#def elmayorineintervals(alpha0=None, omoriWaitExp=None, catname='cats/elmayor-mfetas.cat', dorefresh=True, mc=3.0, drho=5.0, p=1.1):
def elmayorintervals(catname='cats/elmayor-mfetas.cat', dorefresh=True, mc=3.0, mstar=None, C1=5.95, C2=0.0, p=1.15, skipexp=3.0, lw=2., mev=None, Nfact=1.0, fnum=1, titlename='Elmayor Intervals', ralpha=-5., rbeta=-.5, fitFract=1.0):
	# good fits for 3D dist. using: a=ybp.elmayorDists(C1=3.7, C2=2.59, q=1.5, mc=2.5)
	# good fits for 2D dist using: a=ybp.elmayorDists(C1=3.03, C2=2.59, q=1.33, mc=2.5)
	# this suggests that the intrinsic, 3D, slope ia about q=1.5
	#[31.0, 35.0]	[-118.0, -113.5]
	#
	#print "doing el-mayor..."
	cpf = elmayorCat(catname=catname, dorefresh=dorefresh, mc=mc)
	return calcEqints(cat=cpf, mstar=mstar, mc=mc, C1=C1, C2=C2, p=p, skipexp=skipexp,lw=lw, ralpha=ralpha, rbeta=rbeta, fitFract=fitFract)

def elmayorDists(catname='cats/elmayor-mfetas.cat', mstar=None, mc=2.5, C1=2.8, C2=3.77, q=1.5, skipexp=3.0, lw=2, mev=None, Nfact=1.0, fnum=1, titlename='El Mayor distances', mt=7.6, len0=None, catnum=0, dorefresh=True, ralpha=.5, rbeta=0.0):
	cpf = elmayorCat(catname=catname, dorefresh=dorefresh, mc=mc)
	
	return calcEqDists(cat=cpf, mstar=mstar, mc=mc, C1=C1, C2=C2, q=q, skipexp=skipexp, lw=lw, mev=mev, Nfact=Nfact, fnum=fnum, titlename=titlename, mt=mt, len0=len0, catnum=catnum, ralpha=ralpha, rbeta=rbeta)


def elmayorCat(catname='cats/elmayor-mfetas.cat', dorefresh=True, mc=None):
	if dorefresh==True:
		hmcat = bcp.atp.catfromANSS(lon=[-118.0, -113.5], lat=[31.0, 35.0], minMag=mc, dates0=[dtm.datetime(2010, 4,1, tzinfo=pytz.timezone('UTC')), dtm.datetime(2012,12,31, tzinfo=pytz.timezone('UTC'))], Nmax=999999, fout=catname)
		cpf=ypp.eqcatalog(hmcat)
	else:
		cpf=ypp.eqcatalog()
		cpf.loadCatFromFile(catname, minmag=mc)
	#	
	thiscatnum=0
	mev=cpf.getMainEvent(cpf.getcat(thiscatnum))
	#	
	print(mev)
	mstar=mev[3] #7.1
	#
	return cpf


def sansimintervals(catname='cats/sansim-mfetas.cat', dorefresh=True, mc=3.0, mstar=None, C1=5.95, C2=0.0, p=1.1, skipexp=3.0, lw=2., mev=None, Nfact=1.0, fnum=1, titlename='Sansim Intervals', ralpha=-5., rbeta=-.5, fitFract=1.0):
	#
	cpf = sansimCat(catname=catname, dorefresh=dorefresh, mc=mc)
	#
	thiscatnum=0
	mev=cpf.getMainEvent(cpf.getcat(thiscatnum))
	#	
	print(mev)
	if mstar==None: mstar=mev[3] # 6.5
	#
	return calcEqints(cat=cpf, mstar=mstar, mc=mc, C1=C1, C2=C2, p=p, skipexp=skipexp,lw=lw, ralpha=ralpha, rbeta=rbeta, fitFract=fitFract)

def sansimDists(catname='cats/sansim-mfetas.cat', mstar=None, mc=3.5, b=1.0, C1=2.8, C2=0., q=1.5, skipexp=2.0, lw=2, mev=None, Nfact=1.0, fnum=1, titlename='SanSim distances', mt=7.6, len0=None, catnum=0, dorefresh=True, ralpha=.5, rbeta=1.0):
	cpf = sansimCat(catname=catname, dorefresh=dorefresh, mc=mc)
	#
	return calcEqDists(cat=cpf, mstar=mstar, mc=mc, b=b, C1=C1, C2=C2, q=q, skipexp=skipexp, lw=lw, mev=mev, Nfact=Nfact, fnum=fnum, titlename=titlename, mt=mt, len0=len0, catnum=catnum, ralpha=ralpha, rbeta=rbeta)
	
def sansimCat(catname='cats/sansim-mfetas.cat', dorefresh=True, mc=3.0):
	dLL=.7
	x0=-121.10
	y0=35.71
	if dorefresh==True:
		hmcat = bcp.atp.catfromANSS(lon=[x0-dLL-.15, x0+dLL-.15], lat=[y0-dLL-.2, y0+dLL-.2], minMag=mc, dates0=[dtm.datetime(2003,12,12, tzinfo=pytz.timezone('UTC')), dtm.datetime(2012,12,31, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC'))], Nmax=999999, fout=catname)
		cpf=ypp.eqcatalog(hmcat)
	else:
		cpf=ypp.eqcatalog()
		cpf.loadCatFromFile(catname, minmag=mc)
	
	return cpf
	

#def coalingaintervals(alpha0=None, catname='cats/coalinga-mfetas.cat', dorefresh=True, mc=3.0, drho=5.0, p=1.1, lw=1):
def coalingaintervals(catname='cats/coalinga-mfetas.cat', dorefresh=True, mc=3.5, mstar=None, C1=5.95, C2=0.0, p=1.15, skipexp=3.0, lw=2., mev=None, Nfact=1.0, fnum=1, titlename='Coalinga Intervals', ralpha=-5., rbeta=-.5, fitFract=1.0):
	#[31.0, 35.0]	[-118.0, -113.5]
	#
	cpf = coalingaCat(catname=catname, dorefresh=dorefresh, mc=mc)
	thiscatnum=0
	mev=cpf.getMainEvent(cpf.getcat(thiscatnum))
	#	
	print(mev)
	if mstar==None: mstar=mev[3] # 6.5
	#
	return calcEqints(cat=cpf, mstar=mstar, mc=mc, C1=C1, C2=C2, p=p, skipexp=skipexp,lw=lw, ralpha=ralpha, rbeta=rbeta, fitFract=fitFract)

def coalingaDists(catname='cats/coalinga-mfetas.cat', mstar=None, mc=3.0, b=1.0, C1=2.8, C2=0., q=1.5, skipexp=2.0, lw=2, mev=None, Nfact=1.0, fnum=1, titlename='Coalinga distances', mt=7.6, len0=None, catnum=0, dorefresh=True, ralpha=.5, rbeta=1.0):
	cpf = coalingaCat(catname=catname, dorefresh=dorefresh, mc=mc)
	
	#dfits=omoriFit(X, Y)
	
	return calcEqDists(cat=cpf, mstar=mstar, mc=mc, b=b, C1=C1, C2=C2, q=q, skipexp=skipexp, lw=lw, mev=mev, Nfact=Nfact, fnum=fnum, titlename=titlename, mt=mt, len0=len0, catnum=catnum, ralpha=ralpha, rbeta=rbeta)

def coalingaCat(catname='cats/coalinga-mfetas.cat', dorefresh=True, mc=3.0):
	dLL=.7
	x0=-120.32
	y0=36.22
	#
	xshift=0.0
	yshift=0.0
	if dorefresh==True:
		hmcat = bcp.atp.catfromANSS(lon=[x0-dLL-xshift, x0+dLL-xshift], lat=[y0-dLL-yshift, y0+dLL-yshift], minMag=mc, dates0=[dtm.datetime(1983, 5, 1, tzinfo=pytz.timezone('UTC')), dtm.datetime(2003,12,31, 0, 0, 0, 0, tzinfo=pytz.timezone('UTC'))], Nmax=999999, fout=catname)
		cpf=ypp.eqcatalog(hmcat)
	else:
		cpf=ypp.eqcatalog()
		cpf.loadCatFromFile(catname, minmag=mc)
	#	
	#cpf.loadCatFromFile(fname='cats/hmine.cat', minmag=3.0)
	#if alpha0==None:
	#	alpha0=alpha_0
	#
	return cpf
	

#def tohokuintervals(alpha0=None, omoriWaitExp=None, catname='cats/tohoku-mfetas.cat', dorefresh=False, mc=4.5, drho=5.0, p=1.1):
def tohokuintervals(catname='cats/tohoku-mfetas.cat', dorefresh=True, mc=4.5, mstar=None, C1=9.5, C2=0.0, p=1.1, skipexp=3.0, lw=2., mev=None, Nfact=1.0, fnum=1, titlename='Tohoku', ralpha=-5., rbeta=-.5, fitFract=1.0):
	#
	# we'll need, at least, a large earthquake renormalization for Tohoku (aka, log(N) = b1(mt-mc) + b2(m-mt)
	# #lon=[-117.25, -115.5]	lat=[33.5, 35.5]	m0=3.000000	dates=[datetime.datetime(1980, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), datetime.date(2011, 5, 12, 0, 0, 0, 0, pytz.timezone('UTC'))]
	# note: tohoku also needs a modified BASS formulation to reflect the very large magnitude (b=1.5 for m>7.6).
	#
	if dorefresh==True:
		#print "script catalog fetch."
		# lon=[135.0, 150.0]	lat=[30.0, 41.5]	m0=4.000000	dates=[datetime.date(2005, 1, 1, 0, 0, 0, 0, pytz.timezone('UTC')), datetime.date(2011, 7, 12, 0, 0, 0, 0, pytz.timezone('UTC'))]
		hmcat = bcp.atp.catfromANSS(lon=[135.0, 150.0], lat=[30.0, 41.5], minMag=mc, dates0=[dtm.datetime(2011, 3,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], Nmax=9999999, fout=catname)
		cpf=ypp.eqcatalog(hmcat)
	else:
		cpf=ypp.eqcatalog()
		cpf.loadCatFromFile(fname=catname, minmag=mc)
	
	useCatnum=0
	#
	mev=cpf.getMainEvent(cpf.getcat(useCatnum))
	#	
	print(("mainshock: ", mev, len(cpf.getcat(useCatnum))))
	mstar=9.0	# hmine is 7.1 according to ANSS; we might check with USGS source.
	#
	effmstar=None
	#
	#return eqintervals(cat=cpf, mstar=mstar, mc=mc, mev=mev, catnum=useCatnum, alpha0=alpha0, omoriWaitExp=omoriWaitExp, pfwaittime=thispfwait)
	#return eqintervals(cat=cpf, mstar=mstar, mc=mc, mev=mev, catnum=useCatnum, alpha0=alpha0, omoriWaitExp=omoriWaitExp, effmstar=effmstar, drho=drho, p=p)
	return calcEqints(cat=cpf, mstar=mstar, mc=mc, C1=C1, C2=C2, p=p, skipexp=skipexp,lw=lw, ralpha=ralpha, rbeta=rbeta, fitFract=fitFract)
	#
def resAlphaOmori(p, y, x, mstar, mc, pq=1.2, w=1.0, alpha=None, beta=None, gamma=None):
	#p[0]=-.5
	#p[2]=10.2
	if alpha!=None: p[0] = alpha
	if beta!=None: p[1] = beta
	if gamma!=None: p[2] = gamma
	lt0 = p[0]*mstar + p[1]*mc + p[2]
	t0=10**lt0
	ltau = (1.0-pq)*lt0 - 1.0*(mstar-1.0-mc) - math.log10(pq-1.0)
	tau = 10**ltau
	#
	# fomorilike(X, tau, x0, p)
	#errs = (y-fomorilike(x, tau, t0, pq))*w
	errs=[]
	for i in range(len(x)):
		ymodel = fomorilike(x[i], tau, t0, pq)
		if ymodel==0.0 or y[i]==0.0: continue
		#
		errs+=[math.log10(float(y[i])/ymodel)]	
	errs=scipy.array(errs)**2.0
	#
	#print("ressing...", len(errs), sum(errs**2.0), sum(errs))
	return errs

#
def fitInts(X,Y, prams0=[None, None, None], mstar=None, mc=None, p=1.2, alpha=None, beta=None, gamma=None):
	# this'll be a bit difference. prams are [alpha, beta, gamma] where
	# t0 = alpha*m + beta*mc + gamma. we'll assume p=1.2 for now (not a big difference) and tau from GR constraints.
	P=scipy.array(prams0, copy=True)
	if P[0]==None: P[0] = -.5
	if P[1]==None: P[1] = -.5
	if P[2]==None: P[2] = 3.0
	#
	# resAlphaOmori(p, y, x, mstar, mc, pq=1.2, w=1.0)
	plsq=spo.leastsq(resAlphaOmori, P, args=(scipy.array(Y), scipy.array(X), mstar, mc, p, scipy.ones(len(X)), alpha, beta, gamma))
	#
	return plsq

def calcEqints(cat=None, mstar=None, mc=None, C1=5.9, C2=0., p=1.1, skipexp=3.0, lw=2, mev=None, Nfact=1.0, fnum=1, titlename='intervals', mt=7.6, catnum=0, ralpha=0.0, rbeta=-.5,fitFract=1.0):
	#
	# this probably needs to be, more or less, rewritten. recode the parameter estimates based on the new write-up (is this consistent with 27 aug 2012 below?).
	# the free prams will be [p, q, lambda, alpha, dm, lambda' alpha'] -- where one of those can be eliminated from:
	# log((1-p)(1-q)) = ( lambda + alpha + dm + - lambda' - alpha' )
	#
	# fitFract: fraction of data to fit (for now, front-fraction).
	# this is an experimental eqintervals version, consistent with present write-up (27 aug 2012).
	# default cat for dev:
	#if cat==None:
	#	cat=eqp.eqcatalog([])
	#	cat.loadCatFromFile('cats/sansim-mfetas.cat', minmag=3.0)
	#	mc=3.0
	#
	# note: dmstar=1.0 gives full sequence type behavior
	#       dmstar=1.2 gives "one earthquake" behaviro (aka, sum(all events with dmstar=1.2) -> dmstar=1.0
	#if alpha0==None: alpha0=alpha_0
	#alpha=alpha0
	#b=1.0
	cpf=cat # catalog object
	if type(cpf).__name__ == 'list':
		#cpf=ypp.eqcatalog(cpf)
		eqp=ypp.eqcatalog(cpf)
	if len(cpf.subcats)==0: catnum=0
	if catnum>len(cpf.subcats): catnum=len(cpf.subcats)
	if mev==None: mev=cat.getMainEvent(cat.getcat(catnum))
	#
	if mstar==None:
		# assume largest mag in catalog:
		mstar=cpf.getMainEvent(cpf.getcat(catnum))[3]
		print(("mstar chosen: ", mstar))
	#
	print(("catrange: ", cpf.cat[0][0], cpf.cat[-1][0]))
	#
	cpf.plotGRdist(fignum=0, doShow=True, fname=None)
	eqints=cpf.getIntervals(catList=cpf.getcat(catnum), winLen=1)	# intervals are in days, right?
	#eqints.reverse()
	eqints.sort(key=lambda x: x[0])
	#
	# we can also do this; it's probably faster (should probably optimize eqcatalog.getIntervals() )
	#eqdts=map(mpd.date2num, map(operator.itemgetter(0), cpf.getcat(0)))
	#eqdts.sort()
	#eqints=eqdts[1:] - eqdts[:-1]
	#
	while eqints[0][0]<=(mev[0]):
		eqints.pop(0)
	# clean up dt=0 (screws up log-transforms). this might screw up counting, but dt=0 might also just be bad data.
	i=1
	while i<len(eqints):
		while eqints[i][1]==0.0: eqints.pop(i)
		i+=1
	#
	#skipIntervals=int(10**(mstar-mc-skipexp))
	skipIntervals=int(10**(mstar-mc-skipexp))
	#
	ts=[0]
	Nactual=[0]
	Ntheory=[1]
	
	NactSkip=[0]
	tskip=[0]
	#
	#ralpha = -.5
	#rbeta = math.log10(2.0)
	#rbeta = -.5
	lt0 = ralpha*mstar + rbeta*mc + C1
	#lt0 = mc + C1
	#ltau = mc*(3.0-p)/2.0 - mstar -math.log10(p-1.) # + C2
	#
	if mstar<=mt:
		#ltau = mc*(p+1.)/2.0 - mstar -math.log10(p-1.) + C2
		ltau = (1.-p)*lt0 - math.log10(p-1.) - 1.0*(mstar - 1.0 - mc)
		
	else:
		#lt0 = -mc + C1
		ltau = (1.-p)*lt0 - 1.0*(mt-mc) - 1.5*(mstar-mt) - math.log(p-1.0) + C2	# where we have assumed b1=1.0, b2=1.5
																										# it might also be necessary to average over
	t0=10**lt0																						# the large/small mag. domain (to get a more
	tau = 10**ltau																					# accurate result).
	#
	for dt in eqints:
		#thisint = dt[1]
		thisint=dt[1]*days2secs		# convert intervals (and occurrence times) to seconds.
		#ts+=[(ts[-1]+dt[1])]		# getIntervals() returns intervals in days; BASS caclucates rates in seconds.]
		ts+=[(ts[-1] + thisint)]	
		Nactual+=[Nactual[-1]+1]
		
		#thisNth = ((myt0 + ts[-1])**(1.0-p) - myt0**(1.0-p)) / (mytau * (1.0-p))
		#Ntheory+=[thisNth]
		#t=ts[-2]+dt[1]/2.0
		if len(ts)>=skipIntervals:
			tskip+=[(tskip[-1] + thisint)]
			NactSkip+=[NactSkip[-1]+1]
	#
	# fitInts(X,Y, prams0=[None, None, None], mstar=None, mc=None, p=1.2)
	#plsq = fitInts(X=scipy.array(ts), Y=scipy.array(Nactual), prams0=[ralpha, rbeta, C1], mstar=mstar, mc=mc, p=p)
	if fitFract==0.0: fitFact=1.0
	fitlen=int(len(ts)*fitFract)
	#	
	if fitlen>0:
		plsq = fitInts(X=scipy.array(ts)[0:fitlen], Y=scipy.array(Nactual[0:fitlen]), prams0=[ralpha, rbeta, C1], mstar=mstar, mc=mc, p=p, alpha=ralpha, beta=rbeta, gamma=None)
	if fitlen<0:
		plsq = fitInts(X=scipy.array(ts)[fitlen:], Y=scipy.array(Nactual[fitlen:]), prams0=[ralpha, rbeta, C1], mstar=mstar, mc=mc, p=p, alpha=ralpha, beta=rbeta, gamma=None)
	#
	print("plsq: ", plsq, mc, mstar)
	lt0fit = plsq[0][0]*mstar + plsq[0][1]*mc + plsq[0][2]
	ltaufit = (1.-p)*lt0fit - math.log10(p-1.) - 1.0*(mstar - 1.0 - mc)
	t0fit = 10**lt0fit
	taufit = 10**ltaufit
	#
	Ntheory=[0]
	Ntheory2=[0]
	#Ntheory3=[0]	# "local rate" method.
	Ntheoryfit=[0]
	tsTheory=[0.]
	#
	logmaxtime=math.log10(ts[-1])
	#print "times: %f, %f, %f" % (ts[0], ts[1], ts[-1])
	#print "logmaxtime: %f" % logmaxtime
	#
	# note: mean rate for san-sim catalog approx 15.52 eq>3/yr, 4.92*10**-7 eq>3/sec
	meanBGrate = 4.98*10**-7
	lttheory=-4.0	# starting log-time
	#lttheory = lt0+6.
	dttheory=(logmaxtime - lttheory)/10000	# time interval.
	thisindex=0
	#
	while lttheory<logmaxtime:
		t1=10**(lttheory)
		t2=10**(lttheory+dttheory)
		#
		thisNth = Nfact*((t0 + t2)**(1.0-p) - t0**(1.0-p)) / (tau * (1.0-p))
		#if thisNth<1.0: thisNth=1.0
		#thisdNth = ((myt0 + t2)**(1.0-p) - (myt0 + t1)**(1.0-p)) / (mytau * (1.0-p)) + meanBGrate*(t2-t0)
		#thisNth2 = Nfact*((myt0 + t2)**(1.0-p) - myt0**(1.0-p)) / (mytau * (1.0-p)) + meanBGrate*(t2-myt0)
		#
		#thisNth2 = thisNth + meanBGrate*(t2-t0)
		thisNth2 = thisNth + meanBGrate*(t2)
		#
		#thisNth2 =((myt0 + t2)**(1.0-p) - myt02**(1.0-p)) / (mytau2 * (1.0-p)) 
		#thisNth2 = thisNth  + meanBGrate*(t2-t1)
		#
		Ntheory  += [thisNth]
		#Ntheory2 += [Ntheory[-1] + thisdNth]
		Ntheory2  += [thisNth2]
		Ntheoryfit += [Nfact*((t0fit + t2)**(1.0-p) - t0fit**(1.0-p)) / (taufit * (1.0-p))]
		#
		#Ntheory += [Ntheory[-1] + omoriRateInt(t1=t1, t2=t2,  t0=t0fact*myt0, tau=mytau, p=p)]
		#
		#Ntheory2 += [Ntheory2[-1] + omoriRateInt(t1=t1, t2=t2,  t0=t02, tau=tau2, p=p)]
		#
		#Ntheory3 += [Ntheory3[-1] + eqobj.omoriN(t1=t1, t2=t2,  t0=None, tau=None, p=None)]
		#Ntheory4 += [eqobj.omoriN(t1=0.0, t2=t2)]
		#
		tsTheory += [tsTheory[-1] + (t2-t1)]
		thisindex+=1
		lttheory+=dttheory
	#
	
	#
	Nom=10**(1.0*(mstar-1.0-mc))
	print("Ns: Ngr=%f, Nth=%f, Nmax=%d, Rth=%f, Rmax=%f" % (Nom, Ntheory[-1], Nactual[-1], Ntheory[-1]/Nom, Nactual[-1]/Nom))
	plt.figure(fnum)
	plt.clf()
	plt.loglog(ts, Nactual, '-.', label='nActual', lw=lw, zorder=4)
	plt.loglog(tskip, NactSkip, '--', label='nSkipping', lw=lw, zorder=5)
	plt.loglog(scipy.array(tsTheory), scipy.array(Ntheory), '-.', label='theory', alpha=.7, lw=2*lw, zorder=3)
	plt.loglog(scipy.array(tsTheory), scipy.array(Ntheoryfit), '-.', label='fit', alpha=.7, lw=lw, zorder=3)
	
	#plt.loglog(scipy.array(tsTheory), scipy.array(Ntheory2), '--', label='theory2', lw=2*lw, zorder=3)
	plt.xlabel('time $t$')
	plt.ylabel('$N(<t)$')
	#strtitle = 
	plt.title('%s Intervals, $m_c=%.2f$, $C_1 = %.2f$, $C_2 = %.2f$ \n $p=%.2f$, $N_{fact}=%.2f$' % (titlename, mc, C1, C2, p, Nfact))
	plt.legend(loc='best')
	#
	plt.figure(2)
	plt.clf()
	plt.loglog(ts, Nom - scipy.array(Nactual), '-.', label='nActual', lw=lw, zorder=4)
	plt.loglog(tskip, Nom - scipy.array(NactSkip), '--', label='nSkipping', lw=lw, zorder=5)
	plt.loglog(scipy.array(tsTheory), Nom - scipy.array(Ntheory), '-.', label='theory', alpha=.7, lw=2*lw, zorder=3)
	#
	############################
	# rates and intervals:
	nbins=100
	donorm=False
	dts=[]
	rts=[]
	#for t in ts:
	avlen=1
	#print "ts:", ts[0:10], ts[-10:]
	for i in range(avlen,len(ts)):
		#if t<=0: continue
		#lts+=[math.log10(t)]
		#dts+=[(abs(ts[i]-ts[i-avlen]))/float(avlen)]
		thisdt=(abs(ts[i]-ts[i-avlen]))
		if thisdt<=0:
			#thisdt=None
			dts+=[None]
			rts+=[None]
			continue
		dts+=[thisdt]
		rts+=[1.0/dts[-1]]
	#dts=scipy.array(ts)
	#rts = 1.0/dts
	#
	t0len=50
	plt.figure(3)
	plt.clf()
	plt.title("Interval")
	plt.ylabel('$\\Delta t$')
	plt.xlabel('time $t$')
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(ts[avlen:], dts, '.-.')
	plt.plot(ts[avlen:avlen+t0len], dts[:t0len], '.')
	
	#print "mean Initial Interval: %f/%f" % (sum(dts[:t0len])/float(t0len), math.log10(sum(dts[:t0len])/float(t0len)))
	lmean=sum(map(math.log10, dts[:t0len]))/float(t0len)
	lsdev= scipy.array(list(map(math.log10, dts[:t0len]))).std()
	print("mean log-Initial Interval: %f/%f, %f (%.2f, %.2f)" % (10**(lmean),lmean, lsdev, mstar, mc ))
	
	#
	plt.figure(4)
	plt.clf()
	plt.title("Rate")
	#plt.ylabel('$\\frac{dN}{dt}$')
	plt.xlabel('time $t$')
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(ts[avlen:], rts, '.-.')
	plt.plot(ts[avlen:avlen+t0len], rts[:t0len], '.')
	#T=map(operator.itemgetter(0), eqints)
	#T=mpd.date2num(scipy.array(T))*days2secs
	#plt.plot(T, map(operator.itemgetter(1), eqints))
	#
	return [eqints, plsq]
	#
#
def intervalFitFig(flist=None, fdir='data'):
	if flist==None or 1==1:
		# we're going to have to hard-code some of these bits.
		flist = ['coalingafits.txt', 'elmayorfits.txt', 'hminefits.txt', 'sansimfits.txt', 'pfbetafits.txt', 'tohokufits.txt']
		#
	
	ms=[6.7,7.2, 7.1, 6.5, 6.0, 9.0]
	
	fnum=27
	
	plt.figure(fnum)
	plt.clf()

	#
	bints=[]
	bslopes=[]
	beta0s=[]
	#
	for fl in flist:
		fname = fdir + '/' + fl
		mcfitmin=None
		mcfitmax=None
		if 'tohoku' in fl:
			mcfitmax=6.1
		
		A=plotcalcfits(fname=fname, mcfitmin=mcfitmin, mcfitmax=mcfitmax)
		plt.figure(fnum)
		flabel=fl.split('fits')[0]
		plt.plot(A[0], A[2], '.-', label=flabel)
		#
		# now, where is the b=0 intercept?
		lf=lft.linefit([A[0], A[2]])
		plsq=lf.doFit()[0]
		print(plsq)
		bints += [plsq[0]]
		bslopes += [plsq[1]]
		#b0intercept = -plsq[0]/plsq[1]
		beta0s+=[-plsq[0]/plsq[1]]
	#
	plt.figure(fnum)
	plt.legend(loc='best', numpoints=1)
	plt.xlabel("$\\beta$")
	plt.ylabel("$b$ : $C = C_0 + b \\cdot m_c$")
	#
	plt.figure(1)
	plt.clf()
	#plt.plot(ms, beta0s, 'o', label='$\\beta_0$')
	lf2=lft.linefit([ms, beta0s])
	plsq2=lf2.doFit()[0]
	plt.plot([min(ms), max(ms)], [plsq2[0] + min(ms)*plsq2[1], plsq2[0] + max(ms)*plsq2[1]], lw=3)
	for i in range(len(ms)):
		plt.plot([ms[i]], [beta0s[i]], 'd', ms=10, label='$$m=%.2f$, \\beta_0 = %.2f$' % (ms[i], beta0s[i]))
	meanm=numpy.mean(ms[:-2])
	meanbeta0 = numpy.mean(beta0s[:-2])
	#
	plt.plot([meanm], [meanbeta0], '*', ms=18, label='$<\\beta_0> = %.2f$' % meanbeta0)
	plt.legend(loc='best', numpoints=1)
	plt.xlabel('m_c')
	plt.ylabel('$\\beta_0$')

		
#
def plotcalcfits(fname='data/pfbetafits.txt', fnum=0, mcfitmin=None, mcfitmax=None):
	# plot a calcfit file...
	## El Mayor-Cucapah
	#2.500000	-1.250000	13.993566
	#2.600000	-1.250000	13.955066
	plt.ion()
	X=[]
	Y=[]
	bs=[]
	avals=[]
	means=[]
	betas=[]
	thisbeta=None
	#fnum=0
	plt.figure(fnum)
	plt.clf()
	f=open(fname, 'r')
	for rw in f:
		if '#!dataname' in rw: plotname=rw.split()[1][:-1]
		if rw[0]=='#': continue	# though eventually we'll want meta-data lines with #!
		#
		rws = list(map(float, rw.split()))
		mc, newbeta, C = rws[0], rws[1], rws[2]
		if newbeta!=thisbeta:
			# new plot-series.
			
			if len(X)!=0:
				#
				print("betas: %f, %f, %d" % (newbeta, thisbeta, len(X)))
				lf=lft.linefit([X,Y])
				fprams = lf.doFit(xmin=mcfitmin, xmax=mcfitmax)[0]
				a=fprams[0]
				b=fprams[1]
				mean = sum(Y)/float(len(Y))		# mean value of C or gamma variable.
				avals+=[a]
				bs+=[b]
				means+=[mean]
				betas+=[thisbeta]
				#
				# plot data:
				plt.plot(X,Y, '.-', label='$\\beta=%f: a=%f, b=%f, \\mu=%.2f$' % (thisbeta, a, b, mean))
				plt.grid(True, which='both')
			#
			# reset"
			X=[]
			Y=[]
			#bs=[]
			#means=[]
			thisbeta = newbeta
		X+=[mc]
		Y+=[C]
	f.close()
	#
	plt.legend(loc='best')
	#
	plt.figure(fnum+1)
	plt.clf()
	plt.plot(betas, bs, '.-')
	#plt.legend(loc='best')
	plt.title('b-values')
	plt.xlabel('$\\beta$')
	plt.ylabel('slope $b$ of $C(\\beta)$')
	#
	plt.figure(fnum+2)
	plt.clf()
	plt.plot(betas, means, '.-')
	#plt.legend(loc='best')
	plt.xlabel('$\\beta$')
	plt.ylabel('mean $C(\\beta)$')
	
	return [betas, avals, bs, means]
#

def plotAlphaBetaData(fitset='parkfield', p=1.08, b=1.0, dm=1.0, alpha0=None, beta0=None, gamma0=None, debug=False, fnum=0):
	#
	if beta0==None: beta0=.5
	if alpha0==None: alpha0=.5
	if gamma0==None: gamma0=0.
	#
	cpf=ypp.eqcatalog()
	useCatnum=0
	if fitset=='parkfield':
		mc=1.5
		mstar=5.97
		cpf.loadCatFromFile(fname='cats/parkfield-mfetas.cat', minmag=mc)
		cpf.addEllipCat('PFshock (.4 x .15)', cpf.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		useCatnum=1
	#
	if fitset=='elmayor':
		mc=2.5
		mstar=7.2
		cpf.loadCatFromFile(fname='cats/elmayor-mfetas.cat', minmag=mc)
	#
	if fitset=='tohoku':
		mc=4.5
		mstar=9.0
		cpf.loadCatFromFile(fname='cats/tohoku-mfetas.cat', minmag=mc)
	#
	if fitset=='sansim':
		mc=3.0
		mstar=6.5
		cpf.loadCatFromFile(fname='cats/sansim-mfetas.cat', minmag=mc)
	#
	if fitset=='coalinga':
		mc=3.5
		mstar=6.7
		cpf.loadCatFromFile(fname='cats/coalinga-mfetas.cat', minmag=mc)
	#	
	if fitset=='hmine':
		mc=3.0
		mstar=7.1
		cpf.loadCatFromFile(fname='cats/hmine-mfetas.cat', minmag=mc)
	#	
	mev=cpf.getMainEvent(cpf.getcat(useCatnum))
	cpf.getcat(useCatnum).sort(key=lambda x: x[0])
	#
	while(cpf.getcat(useCatnum)[0][0]<mev[0]): 
		cpf.getcat(useCatnum).pop(0)
	#
	fmev=mpd.date2num(mev[0])
	X = (scipy.array(list(map(mpd.date2num, list(map(operator.itemgetter(0), cpf.getcat(useCatnum)))))) - fmev)*days2secs
	Y = numpy.arange(1., len(X)+1)
	#lY = map(math.log10, Y)
	ty = []
	tx=[]
	for t in X:
		ly=logNOmoriScaling(mstar, mc, t, p=p, alpha=alpha0, beta=beta0, gamma=gamma0, dm=dm, b=b)
		if ly==None: continue
		ty+=[10**(ly)]
		tx+=[t]
	#
	plt.figure(fnum)
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(X,Y, '.')
	plt.plot(tx, ty, '-')
	#
	return [X,Y,ty]

def showbrutefits():
	fls=glob.glob('data/brute/*.dat')
	fls.sort()
	alphas1, alphas2, alphas3=[],[],[]
	betas1, betas2, betas3=[],[],[]
	gammas=[]
	sets=[]
	for fl in fls:
		setname=fl.split('/')[-1].split('.')[0]
		sets+=[setname]
		#
		f=open(fl)
		for rw in f:
			thisrw=rw
		print(fl + ": " + rw)
		f.close()
		#
		rws=list(map(float, rw.split()))
		if 'alpha' in setname:
			alphas2+=[rws[0]]
			betas2 +=[rws[1]]
		elif 'beta' in setname:
			alphas3+=[rws[0]]
			betas3 +=[rws[1]]
		else:
			alphas1+=[rws[0]]
			betas1 +=[rws[1]]
		gammas +=[rws[2]]
	#
	plt.figure(0)
	plt.clf()
	plt.ion()
	plt.plot(alphas1, betas1, 'o')
	plt.plot(alphas2, betas2, 'd')
	plt.plot(alphas3, betas3, 's')

def logNOmoriScaling(m, mc, t, p=1.08, alpha=-.5, beta=-.5, gamma=9.0, dm=1.0, b=1.0):
	t0=10**(alpha*m + beta*mc + gamma)
	#print "val lNomori: ", (1.0 + t/t0)**(1.0-p)
	if (1.0 + t/t0)**(1.0-p)>=1.0: return None
	#N=b*(m-dm-mc) + math.log10(p-1.0) + math.log10( (1.0 + t/t0)**(1.0-p) - 1.0)
	N=b*(m-dm-mc) + math.log10( 1.0 - (1.0 + t/t0)**(1.0-p))
	return N

def brutefitAlphaBetaGamma(fitset='parkfield', nits=10**6, p=1.08, b=1.0, dm=1.0, alpha0=None, beta0=None, gamma0=None, debug=False, writefile=True):
	# just MC the hell out of this and see what we get.
	# fit cumulative distributions to:	logNOmoriScaling()
	#
	#debug=True
	if nits!=None: nits=int(nits)
	#	
	alphas=[-1.5, 1.5]
	betas = [-1.5, 1.5]
	gammas=[-10., 10.]
	R1 = random.Random()
	R2 = random.Random()
	R3 = random.Random()
	#
	fnamestr=''
	if alpha0!=None: fnamestr=fnamestr+'alpha0'
	if beta0!=None: fnamestr=fnamestr+'beta0'
	if gamma0!=None: fnamestr=fnamestr+'gamma0'
	foutname='data/brute/%s%s.dat' % (fitset, fnamestr)
	if debug: print("fout: %s" % foutname)
	#print "select data set."
	#
	if writefile==True:
		fout=open(foutname, 'w')
		fout.write('#MC alpha, beta, gamma fit: a0=%s, b0=%s, g0=%s\n' % (alpha0, beta0, gamma0))
		fout.close()
	#
	XYlY=getXYlY(fitset)
		#
	if debug:
		plt.figure(0)
		plt.clf()
		ax=plt.gca()
		ax.set_xscale('log')
		ax.set_yscale('log')
		plt.plot(X,Y, '.')
	#
	minerr=10**12.
	minprams=[0., 0., 0.]
	i=0
	while nits==None or i<(nits):
	#for i in range(nits):
		if nits!=None: i+=1
		#
		thisalpha=alpha0
		thisbeta=beta0
		thisgamma=gamma0
		#
		if alpha0==None: thisalpha = alphas[0] + (alphas[1]-alphas[0])*R1.random()
		if beta0==None:  thisbeta  = betas[0] + (betas[1]-betas[0])*R2.random()
		if gamma0==None: thisgamma = gammas[0] + (gammas[1]-gammas[0])*R3.random()
		#
		#print "nits: ", i
		err=0
		jcount=0
		for j in range(len(Y)):
			t=X[j]
			thisN=logNOmoriScaling(mstar, mc, t, p=p, alpha=thisalpha, beta=thisbeta, gamma=thisgamma, dm=dm, b=b)
			if thisN==None: continue
			#
			#thiserr = math.log10(Y[j]) - thisN
			thiserr = lY[j] - thisN
			err+=thiserr*thiserr
			jcount+=1
			#
		#print "jcount: ", jcount
		if jcount==0: continue
		err/=float(jcount)
		if err<minerr:
			minprams=[thisalpha, thisbeta, thisgamma]
			if debug: print("new min err: %f, %f, %s" % (minerr, err, str(minprams)))
			if writefile==True:
				fout=open(foutname, 'a')
				fout.write('%f\t%f\t%f\t%f\n' % (minprams[0], minprams[1], minprams[2], err) )
				fout.close()
			#
			minerr=err
	print("minerr: ", minerr, minprams, fitset)
	#
	newNs = []
	for x in X:
		print(newNs[-1], logNOmoriScaling(mstar, mc, x, p=p, alpha=minprams[0], beta=minprams[1], gamma=minprams[2], dm=dm, b=b))	
		newNs += [10**logNOmoriScaling(mstar, mc, x, p=p, alpha=minprams[0], beta=minprams[1], gamma=minprams[2], dm=dm, b=b)]
	if debug: plt.plot(X, newNs, '-', lw=2)
	print("return X,Y,newNs.", len(X), len(Y), len(newNs))
	return [X, Y, newNs]
#
#
def mpbrutes():
	print("execute a set of MC fits as MP processes.")
	#
	mcps=[]
	mcsets=['parkfield', 'elmayor', 'sansim', 'coalinga', 'tohoku', 'hmine']
	#
	for mc in mcsets:
		mcps+=[mcp.Process(target=brutefitAlphaBetaGamma, args=(mc, None, 1.08, 1.0, 1.0, None, None, None, False))]
		mcps+=[mcp.Process(target=brutefitAlphaBetaGamma, args=(mc, None, 1.08, 1.0, 1.0, .5, None, None, False))]
		mcps+=[mcp.Process(target=brutefitAlphaBetaGamma, args=(mc, None, 1.08, 1.0, 1.0, None, .5, None, False))]
		#
	#
	for pr in mcps:
		print("starting process")
		pr.start()
		pr.join
	#
	#return len(mcps)
	
	
#
def calcfitAlphas(rootdir='data', fitfract=1.0):
	#
	# streamlining:
	#fitses=[['parkfield', '%s/pfbetafits.txt' % rootdir, [-1.25, .25],  [1.5, 3.5], pfintervals], ['SanSim', '%s/sansimfits.txt' % rootdir, [-1.25, .25], [3.0, 4.5], sansimintervals], ['coalinga', '%s/coalingafits.txt' % rootdir, [-1.25, .25],  [3.5, 4.5], coalingaintervals], ['Tohoku', '%s/tohokufits.txt' % rootdir, [-1.25, .25],  [3.5, 4.5], tohokuintervals]]
	# [ [name, fname, betarange, mcrange, funct], ...]
	fitses=[['parkfield', '%s/pfbetafitsAlpha.txt' % rootdir, [-1.25, .25],  [1.5, 3.5], pfintervals], ['ElMayor','%s/elmayorfitsAlpha.txt' % rootdir, [-1.25, .25],  [2.5, 4.5], elmayorintervals], ['hmine','%s/hminefitsAlpha.txt' % rootdir, [-1.25, .25],  [3.0, 4.5], hmineintervals], ['SanSim', '%s/sansimfitsAlpha.txt' % rootdir, [-1.25, .25], [3.0, 4.5], sansimintervals], ['coalinga', '%s/coalingafitsAlpha.txt' % rootdir, [-1.25, .25],  [3.5, 4.5], coalingaintervals], ['Tohoku', '%s/tohokufitsAlpha.txt' % rootdir, [-1.75, -.5],  [4.5, 6.0], tohokuintervals]]
	#fitses=[ ['Tohoku', '%s/tohokufitsAlpha.txt' % rootdir, [-1.75, -.5],  [4.5, 6.0], tohokuintervals]]
	fnum=4
	#thisralpha = -.5
	thisrbeta = -.5
	for fprams in fitses:
		thisalpha = fprams[2][0]
		print(fprams)
		Cs=[]
		#
		fname=fprams[1]
		fout=open(fname, 'w')
		fout.write('#beta=%f\n' % thisrbeta)
		fout.write('#mc\talpha\tC1(gamma)\talpha\tbeta\n')
		fout.write('#!dataname\t%s\n' % fitses[0])
		fout.close()
		plt.figure(fnum)
		plt.clf()
		plt.ion()
		#while thisbeta<=fprams[2][1]:
		while thisalpha<=fprams[2][1]:
			print("thisalpha: %f " % thisalpha)
			plotX=[]
			plotY=[]
			thismc=fprams[3][0]
			while thismc<(fprams[3][1]-.5):
				print("thismc: %f" % thismc)
				calcfunct=fprams[4]
				#a=pfintervals(p=1.2, mc=thismc, C1=9.94, ralpha=thisralpha, rbeta=thisbeta)
				a=calcfunct(p=1.2, mc=thismc, C1=9.94, ralpha=thisalpha, rbeta=thisrbeta, fitFract=fitfract)
				plsq=a[1][0]
				Cs+=[[thismc, thisalpha, plsq[2]]]
				fout=open(fname, 'a')
				fout.write('%f\t%f\t%f\t%f\t%f\n' % (thismc, thisalpha, plsq[2], plsq[0], plsq[1]))
				fout.close()
				plotX+=[thismc]
				plotY+=[plsq[2]]
				#
				thismc+=.1
			plt.figure(fnum)
			plt.plot(plotX, plotY, '.-', label='$\\alpha=%f$' % thisalpha)
			plt.title('%s' % fprams[0])
			thisalpha+=.1
		plt.legend(loc='best')
		#
#
def calcfit2(rootdir='data', fitfract=1.0):
	#
	# streamlining:
	#fitses=[['parkfield', '%s/pfbetafits.txt' % rootdir, [-1.25, .25],  [1.5, 3.5], pfintervals], ['SanSim', '%s/sansimfits.txt' % rootdir, [-1.25, .25], [3.0, 4.5], sansimintervals], ['coalinga', '%s/coalingafits.txt' % rootdir, [-1.25, .25],  [3.5, 4.5], coalingaintervals], ['Tohoku', '%s/tohokufits.txt' % rootdir, [-1.25, .25],  [3.5, 4.5], tohokuintervals]]
	# [ [name, fname, betarange, mcrange, funct], ...]
	fitses=[['parkfield', '%s/pfbetafits.txt' % rootdir, [-1.25, .25],  [1.5, 3.5], pfintervals], ['ElMayor','%s/elmayorfits.txt' % rootdir, [-1.25, .25],  [2.5, 4.5], elmayorintervals], ['hmine','%s/hminefits.txt' % rootdir, [-1.25, .25],  [3.0, 4.5], hmineintervals], ['SanSim', '%s/sansimfits.txt' % rootdir, [-1.25, .25], [3.0, 4.5], sansimintervals], ['coalinga', '%s/coalingafits.txt' % rootdir, [-1.25, .25],  [3.5, 4.5], coalingaintervals], ['Tohoku', '%s/tohokufits.txt' % rootdir, [-1.75, -.5],  [4.5, 6.0], tohokuintervals]]
	#fitses=[ ['Tohoku', '%s/tohokufits.txt' % rootdir, [-1.75, -.5],  [4.5, 6.0], tohokuintervals]]
	fnum=4
	thisralpha = -.5
	for fprams in fitses:
		thisbeta = fprams[2][0]
		print(fprams)
		Cs=[]
		#
		fname=fprams[1]
		fout=open(fname, 'w')
		fout.write('#mc\tbeta\tC1\tralpha=%f\n' % (thisralpha))
		fout.write('#!dataname\t%s\n' % fitses[0])
		fout.close()
		plt.figure(fnum)
		plt.clf()
		plt.ion()
		while thisbeta<=fprams[2][1]:
			print("thisbeta: %f " % thisbeta)
			plotX=[]
			plotY=[]
			thismc=fprams[3][0]
			while thismc<(fprams[3][1]-.5):
				print("thismc: %f" % thismc)
				calcfunct=fprams[4]
				#a=pfintervals(p=1.2, mc=thismc, C1=9.94, ralpha=thisralpha, rbeta=thisbeta)
				a=calcfunct(p=1.2, mc=thismc, C1=9.94, ralpha=thisralpha, rbeta=thisbeta, fitFract=fitfract)
				plsq=a[1][0]
				Cs+=[[thismc, thisbeta, plsq[2]]]
				fout=open(fname, 'a')
				fout.write('%f\t%f\t%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2], plsq[0], plsq[1]))
				fout.close()
				plotX+=[thismc]
				plotY+=[plsq[2]]
				#
				thismc+=.1
			plt.figure(fnum)
			plt.plot(plotX, plotY, '.-', label='$\\beta=%f$' % thisbeta)
			plt.title('%s' % fprams[0])
			thisbeta+=.1
		plt.legend(loc='best')
		#
def calcfit():
	thisbeta = -1.25
	Cs=[]
	fnum=3
	
	fname='data/pfbetafits.txt'
	fout=open(fname, 'w')
	fout.write('#mc\tbeta\tC1\n')
	fout.write('# Parkfield\n')
	fout.close()
	plt.figure(fnum)
	plt.clf()
	plt.ion()
	while thisbeta<=.25:
		print("thisbeta: %f " % thisbeta)
		plotX=[]
		plotY=[]
		thismc=1.5
		while thismc<3.5:
			print("thismc: %f" % thismc)
			a=pfintervals(p=1.2, mc=thismc, C1=9.94, ralpha=-.5, rbeta=thisbeta)
			plsq=a[1][0]
			Cs+=[[thismc, thisbeta, plsq[2]]]
			fout=open(fname, 'a')
			#fout.write('%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2]))
			fout.write('%f\t%f\t%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2], plsq[0], plsq[1]))
			fout.close()
			plotX+=[thismc]
			plotY+=[plsq[2]]
			#
			thismc+=.1
		plt.figure(fnum)
		plt.plot(plotX, plotY, '.-', label='$\\beta=%f$' % thisbeta)
		plt.title('Parkfield')
		thisbeta+=.1
	plt.legend(loc='best')
	
	#
	fnum+=1
	plt.figure(fnum)
	fname='data/elmayorbetafits.txt'
	fout=open(fname, 'w')
	fout.write('# El Mayor-Cucapah\n')
	plt.clf()
	plt.ion()
	thisbeta = -1.25
	Cs=[]
	while thisbeta<=.25:
		print("thisbeta: %f " % thisbeta)
		plotX=[]
		plotY=[]
		thismc=2.5
		while thismc<5.0:
			print("thismc: %f" % thismc)
			a=elmayorintervals(p=1.2, mc=thismc, C1=9.94, ralpha=-.5, rbeta=thisbeta, dorefresh=False)
			plsq=a[1][0]
			Cs+=[[thismc, thisbeta, plsq[2]]]
			fout=open(fname, 'a')
			#fout.write('%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2]))
			fout.write('%f\t%f\t%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2], plsq[0], plsq[1]))
			fout.close()
			plotX+=[thismc]
			plotY+=[plsq[2]]
			#
			thismc+=.1
		plt.figure(fnum)
		plt.plot(plotX, plotY, '.-', label='$\\beta=%f$' % thisbeta)
		thisbeta+=.1
	plt.title('El Mayor-Cucapah')
	plt.legend(loc='best')
	plt.show()
	
	#
	fnum+=1
	plt.figure(fnum)
	fname='data/coalingabetafits.txt'
	fout=open(fname, 'w')
	fout.write('# Coalinga\n')
	plt.clf()
	plt.ion()
	thisbeta = -1.25
	Cs=[]
	while thisbeta<=.25:
		print("thisbeta: %f " % thisbeta)
		plotX=[]
		plotY=[]
		thismc=3.0
		while thismc<4.5:
			print("thismc: %f" % thismc)
			a=elmayorintervals(p=1.2, mc=thismc, C1=9.94, ralpha=-.5, rbeta=thisbeta, dorefresh=False)
			plsq=a[1][0]
			Cs+=[[thismc, thisbeta, plsq[2]]]
			fout=open(fname, 'a')
			#fout.write('%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2]))
			fout.write('%f\t%f\t%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2], plsq[0], plsq[1]))
			fout.close()
			plotX+=[thismc]
			plotY+=[plsq[2]]
			#
			thismc+=.1
		plt.figure(fnum)
		plt.plot(plotX, plotY, '.-', label='$\\beta=%f$' % thisbeta)
		thisbeta+=.1
	plt.legend(loc='best')
	plt.title('Coalinga')
	plt.show()

	#
	fnum+=1
	plt.figure(fnum)
	fname='data/sansimbetafits.txt'
	fout=open(fname, 'w')
	fout.write('# San Sim.\n')
	plt.clf()
	plt.ion()
	thisbeta = -1.25
	Cs=[]
	while thisbeta<=.25:
		print("thisbeta: %f " % thisbeta)
		plotX=[]
		plotY=[]
		thismc=3.0
		while thismc<4.5:
			print("thismc: %f" % thismc)
			a=elmayorintervals(p=1.2, mc=thismc, C1=9.94, ralpha=-.5, rbeta=thisbeta, dorefresh=False)
			plsq=a[1][0]
			Cs+=[[thismc, thisbeta, plsq[2]]]
			fout=open(fname, 'a')
			#fout.write('%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2]))
			fout.write('%f\t%f\t%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2], plsq[0], plsq[1]))
			fout.close()
			plotX+=[thismc]
			plotY+=[plsq[2]]
			#
			thismc+=.1
		plt.figure(fnum)
		plt.plot(plotX, plotY, '.-', label='$\\beta=%f$' % thisbeta)
		plt.title('San Simion')
		thisbeta+=.1
	plt.legend(loc='best')
	plt.show()
	#
	fnum+=1
	plt.figure(fnum)
	fname='data/tohokubetafits.txt'
	fout=open(fname, 'w')
	fout.write('# San Sim.\n')
	plt.clf()
	plt.ion()
	thisbeta = -1.25
	Cs=[]
	while thisbeta<=0.:
		print("thisbeta: %f " % thisbeta)
		plotX=[]
		plotY=[]
		thismc=4.5
		while thismc<6.5:
			print("thismc: %f" % thismc)
			a=tohokuintervals(p=1.2, mc=thismc, C1=9.94, ralpha=-.5, rbeta=thisbeta, dorefresh=False)
			plsq=a[1][0]
			Cs+=[[thismc, thisbeta, plsq[2]]]
			fout=open(fname, 'a')
			#fout.write('%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2]))
			fout.write('%f\t%f\t%f\t%f\t%f\n' % (thismc, thisbeta, plsq[2], plsq[0], plsq[1]))
			fout.close()
			plotX+=[thismc]
			plotY+=[plsq[2]]
			#
			thismc+=.1
		plt.figure(fnum)
		plt.plot(plotX, plotY, '.-', label='$\\beta=%f$' % thisbeta)
		thisbeta+=.1
	plt.legend(loc='best')
	plt.title('Tohoku')
	plt.show()
	
def omoriRate(t, t0=None, tau=None, p=1.1):
	rinv=tau*(t0+t)**p
	tfact=1.0
	#tfact=5.0
	return tfact/rinv

def omoriRateInt(t1, t2, t0=None, tau=None, p=1.1):
	deltaN = 1.0/(tau*(1.0-p))*((t0 + t2)**(1.0-p) - (t0+t1)**(1.0-p) )
	return deltaN

def ltauYoder(m=5.0, alpha=alpha_0, p=1.1, b=1.0, dmp=2.0, mc=3.0):
	# actualy log(tau)
	ltau=alpha*(p-1.0) + b*(dmp+mc) + m*((1.0-p)/2.0 - b)
	# to avoid confusion when messing wiht t0, just use t0:
	#thisLogt0 = lt0(m=m)
	#thisLogDr = lDeltat(m=m)
	#ltau = m*(.5-b) + b*(dmp+mc) - alpha - p*thisLogt0		# note, tau will change with t0.
	#ltau = thisLogDr - p*thisLogt0 - b*(m-dmp-mc)
	#
	return ltau

def lDeltat(m, alpha=alpha_0):
	return (m/2.0) - alpha

def lt0(m, alpha=alpha_0):
	# experimentally, 10*deltaT seem to make nice t0 values. let's explicitly make the distinction between
	# these two parameters and see where it gets us...
	lruptTime = lDeltat(m, alpha)
	return 0.0 + lruptTime

#lt0b = (.5*mstar) - alpha0 - dmstar - 1.0 + drho
	#ltaub = mc-.5*mstar-alpha0 - p*lt0b + drho

def lt02 (mstar, alpha0=alpha_0, dmstar=1.2, drho=0.0):
	x=(.5*mstar) - alpha0 - dmstar - 1.0 + drho
	return x

def ltau2(mstar=None, mc=None, alpha0=alpha_0, p=1.1, drho=0.0, thislt0=None):
	if thislt0==None:
		thislt0=lt02(mstar=mstar, alpha0=alpha0, dmstar=1.2, drho=drho)
	x = mc-.5*mstar-alpha0 - p*thislt0 + drho
	return x


# t0 and tau from "local" intensive formulation:
def lt0local(mstar=None, mc=None, l0=2.76, alpha=2.3, drho=0.0, dmstar=1.2, p=1.1):
	#return mc-l0-alpha-1.0 + drho
	x=.5*(mstar + mc) - dmstar + math.log10(p-1.) - (l0+alpha) + drho
	return x

def ltaulocal(mc=None, l0=2.76, alpha=2.3, p=1.1, lt0=None, drho=0.0):
	if lt0==None:
		lt0=lt0local(mc-mc, l0=l0, alpha=alpha, drho=drho)
	return mc-l0-alpha - p*lt0 + drho

	
def getDetBranchingMags(BASScat=[[5.0, 1]], mc=3.0, b=1.0, dmstar=1.0, B=9.0):
	# create a deterministic BASS set. mstar might be a fractional magnitude, so the final magnitude m_f might be <m_c.
	# B -> branching parameter (for realy earthquakes, B=9; in summation, this gives GR with base 10.
	#
	# revise thsis structure:
	# pass an initial BASScat = [mstar]
	# from there, use a recursive loop so that the first iteration is identical to the final iteration.
	# the recursive stack should take care of everything from there (aka, the second loop will walk through the first generated catalog).
	#
	#
	# in this case, go back and correct:
	# m_f -> m_c
	# N_f -> N_f * 10^{b(m_f - m_c)}
	#mc=1.5
	#
	#BASScat=[]
	#thism=mstar
	for mg in BASScat:
		#mstar=BASScat[0][0]
		#thism=BASScat[0][0]
		#thisNeq=BASScat[0][1]
		mstar=mg[0]
		thism=mstar
		thisNeq=mg[1]
		newEvents=[]
		while thism>=(mc+dmstar):
			bpower=(mstar-thism)/dmstar
			thism-=dmstar
			#
			N=thisNeq*(B**bpower)
			newEvents+=[[thism, N]]
		#
		BASScat+=newEvents	
	
	return BASScat

def plotBASScats(BASScat0=[[0, [0,0], 0, 0, 6.0]], mc=3.0, b=1.0, dmstar=1.2, B=9.0, alpha0=alpha_0, nits=100, limMags=True):
	plt.figure(1)
	plt.clf()
	#
	for i in range(nits):
		BASScat=[]
		for rw in BASScat0:
			BASScat+=[rw]
		#
		print("set %d at mc=%f" % (i, mc))
		#a=getDetBranchingCat(BASScat=[[0, [0,0], 0, 0, 6.0]], mc=mc)
		a=getBASS1(BASScat=BASScat, mc=mc, b=b, dmstar=dmstar, alpha0=alpha0, limMags=limMags)
		b=BASScatPlots(BASScat=a, doCLF=(i==0))
	return None
'''
class fcset(object):
	# thoug, in retrospect, i'm not sure i need this class for anything. it can probably be removed.
	lats=[None, None]
	lons=[None, None]
	mc=3.0
	dts=[None, None]
	kmlfile=None	# write to this...
	catalog=[]
	#
	__init__(self, lats=[None, None], lons=[None, None], mc=3.0, dts=[None, None], kmlfile='kml/socal.kml', catfile='kml/socal.cat', catalog=None):
		self.lats=lats
		self.lons=lons
		self.mc=mc
			if dts[1]==None: dts[1]=dtm.datetime.utcnow()
			if dts[0]==None: dts[0]=dtm.datetime.utcnow()-dtm.timedelta(years=5)
		self.dts=dts
		self.kmlfile=kmlfile
		self.catfile=catfile
		self.catalog=catalog
'''

def BASS3demo():
	cat0=[[0, [0,0], 0, 0, 0, 7.0]]
	bcat=getBASS4(BASScat=cat0, mc=1.0, maxDeltaM=3.0, q=1.25)
	mgs=list(map(operator.itemgetter(-1), bcat))
	mgs.sort()
	
	bcat.sort(key = lambda x: x[0])
	
	intervals=getints(bcat)
	justints=intervals[1][:]
	justints.sort()
	
	#print intervals[0][0:5], intervals[0][-5:], intervals[1][0:5], intervals[1][-5:]
	
	Rs=list(map(operator.itemgetter(2), bcat))
	Rs2=list(map(operator.itemgetter(3), bcat))
	Rs.sort()
	Rs2.sort()
	rs=[]
	rs2=[]
	for R in Rs:
		rs+=[R**.5]
	for R in Rs2:
		rs2+=[R**.5]
	#	
	ns=list(range(1, len(Rs)+1))
	ns.reverse()
	#
	plt.figure(0)
	plt.clf()
	plt.ion()
	plt.loglog(rs, ns)
	plt.loglog(rs2, ns)
	plt.xlabel('expected distance $<r>$')
	plt.ylabel('Number $N(>r)$')
	
	plt.figure(1)
	plt.clf()
	plt.semilogy(mgs, ns)
	plt.xlabel('magnitude $m$')
	plt.ylabel('Number $N(>m)$')
	#
	plt.figure(2)
	plt.clf()
	plt.loglog(justints, ns[0:-1], '-o', alpha=.4)
	plt.xlabel('intervals $\Delta t$')
	plt.ylabel('Number $N$')
	#
	plt.figure(3)
	plt.clf()
	a=plotBASSLocs2(BASScat=bcat, fignum=3)
	#
	plt.figure(4)
	plt.clf()
	plt.ion()
	intNs=list(range(1, len(intervals[1])+1))
	plt.semilogy(intNs, intervals[1], 'o', ms=2)
	plt.xlabel('Natural time (event cout) $n$')
	plt.ylabel('interval $\\Delta t$')
	#
	plt.figure(5)
	plt.clf()
	plt.ion()
	#intNs=range(1, len(intervals[1])+1)
	plt.loglog(intervals[0], intervals[1], 'o', ms=2)
	plt.xlabel('time $t$')
	plt.ylabel('interval $\\%Delta t$')
	#
	return bcat

def getBASS4(BASScat=[[0, [0,0], 0, 0, 0, 5.0]], mc=3.0, b=1.0, dmstar=1.2, alpha0=alpha_0, doplots=False, maxDeltaM=3.0, q=1.35):
	# catalog format: [ [time, center, sumR^2, radius, seqNum, mag], ]
	# note, this will not give an [x,y] map. we get a center (defined by the input data) and a distance from that center.
	# the expected distance will be propagated assuming random walk (aka, <r> -> sqrt(sum(r_i^2)).
	# from this, we can draw hazard profiles.
	#
	# same as BASS3 (basicall, retaining only radial component of spatial distribution, for purposes of forecasting), except
	# t, R are calculated deterministically. only magnutide is random, but it is constrained (basically, we shuffle the order).
	# even more optimal would probably be to determine R,t,m all deterministically (using i (0,N) sequentially instead of N*rand())
	# and then shuffle the order, but this should be sufficiently similar to that process, particularly where foreshocks are suppressed
	# or significantly limited.
	#
	# maxDeltaM: use the limited magnitude distribution; this parameter defines the (relative) maximum permitable magnitude for
	# aftershocks (aka, m_max -> (m - dmstar + maxDeltaM)
	#
	if type(BASScat).__name__ in ('float', 'int'):
		BASScat=[[0, [0,0], 0, 0, 0, 0, BASScat]]
	mrand=random.Random()
	#
	for ev in BASScat:
		mstar=ev[-1]
		#thism=mstar
		
		Nbass = 10**(b*(mstar - dmstar - mc))	# Nbass new "child" earthquakes
		if Nbass<=1: continue
		
		parentTime=ev[0]
		parentLoc=ev[1]
		parentDistSq=ev[2]	# distance from original "parent" event (from data).
		
		# Omori Parameters:
		#myt0=10.0**lDeltat(m=mstar, alpha=alpha0)
		myt0=10.0**lt0(m=mstar, alpha=alpha0)
		p=1.1
		b=1.0
		mydmp=1.0-math.log10(p-1.0)/b		# Delta m prime (from Yoder omori fomulation)
		#print "deltaM' %f" % mydmp
		mytau = 10**(ltauYoder(m=mstar, mc=mc, p=p, dmp=mydmp))
		#
		newEvents=[]
		newmags=[]
		# generate events from a branching model
		#NthisEarthquake = int(10**((b/1.1)*(thism - dmstar*1.1 - mc)))	# number of aftershocks for this earthquake.
																							# i think this is approximately right for branching.
		Nsequence = 0	# this will be the number for a single A.S. sequence, excluding sub-sequences.
		#
		# make a list of magnitudes; remove randomness except for order
		maglist=[]
		#print "mstar, Nbass, log(Nbass): %f, %f, %f " % (mstar, Nbass, math.log10(Nbass))
		for i in range(int(Nbass)):
			#
			#newm = randmagLim(mstar=(mstar-dmstar + maxDeltaM), mc=mc, b=b, localrand=mrand)
			newm = invGR(N=i+1.0, mc=mc, dmstar=0.0, b=b)
			maglist+=[newm]

			#maglist+=[[invGR(N=i+1.0, mc=mc, dmstar=dmstar, b=b), mrand.random()]]
			#maglist+=[invGR(N=i+1.0, mc=mc, dmstar=0.0, b=b)]
			#maglist+=[randmagLim(mstar=(mstar-dmstar + maxDeltaM), mc=mc, b=b, localrand=mrand)]
		#print "lens: %d, %d" % (len(maglist), Nbass)
		#maglist.sort(key = lambda x: x[0])
		#random.shuffle(maglist)
		#
		#print deltaN
		for i in range(int(Nbass)):
			Nsequence+=1
			#thist = invOmori(i, mytau, p, myt0)
			#newEvents+=[[parentTime + thist, [0,0], 0, i, thism]]
			#
			# random magnutide between mc and mstar; does not permit foreshocks and mitigates foreshock
			# and BASS singularity related problems. excluding foreshocks does not appear to significantly alter
			# the model, but it sure does cause problems.
			# use the limited magnitude distribution, but add a foreshock facilitator maxDeltaM
			thism = randmagLim(mstar=(mstar-dmstar + maxDeltaM), mc=mc, b=b, localrand=mrand) 	# note we can add a magitude delta to permit foreshocks
			#thism=maglist[i]
			#thist = invOmori((Nbass)*mrand.random(), mytau, p, myt0)
			#thisR=invR(N=(Nbass)*(.95*mrand.random()), m=mstar, mc=mc, dmstar=dmstar, b=b, q=q)
			thist = invOmori(i*.95, mytau, p, myt0)
			thisR=invR(N=i*.95, m=mstar, mc=mc, dmstar=dmstar, b=b, q=q)	# note factor to mitigate infinte distance/time
			#					
			#newmags+=[thism]
			BASScat += [[parentTime + thist, parentLoc, parentDistSq + thisR**2.0, thisR**2.0, i, thism]]
		#
	
	return BASScat	

def getBASS3(BASScat=[[0, [0,0], 0, 0, 0, 5.0]], mc=3.0, b=1.0, dmstar=1.2, alpha0=alpha_0, doplots=False, maxDeltaM=3.0, q=1.35, drho=5.5):
	# catalog format: [ [time, center, sumR^2, radius, seqNum, mag], ]
	# note, this will not give an [x,y] map. we get a center (defined by the input data) and a distance from that center.
	# the expected distance will be propagated assuming random walk (aka, <r> -> sqrt(sum(r_i^2)).
	# from this, we can draw hazard profiles.
	#
	# maxDeltaM: use the limited magnitude distribution; this parameter defines the (relative) maximum permitable magnitude for
	# aftershocks (aka, m_max -> (m - dmstar + maxDeltaM)
	#
	if type(BASScat).__name__ in ('float', 'int'):
		BASScat=[[0, [0,0], 0, 0, 0, 0, BASScat]]
	mrand=random.Random()
	#
	for ev in BASScat:
		mstar=ev[-1]
		#thism=mstar
		parentTime=ev[0]
		parentLoc=ev[1]
		parentDistSq=ev[2]	# distance from original "parent" event (from data).
		
		# Omori Parameters:
		#myt0=10.0**lDeltat(m=mstar, alpha=alpha0)
		#myt0=10.0**lt0(m=mstar, alpha=alpha0)
		p=1.1
		b=1.0
		mydmp=1.0-math.log10(p-1.0)/b		# Delta m prime (from Yoder omori fomulation)
		#print "deltaM' %f" % mydmp
		#mytau = 10**(ltauYoder(m=mstar, mc=mc, p=p, dmp=mydmp))
		mylt0=lt02(mstar=mstar, alpha0=alpha0, dmstar=dmstar, drho=drho)
		myt0 = 10.**(mylt0)
		mytau = 10.**ltau2(mstar=mstar, mc=mc, alpha0=alpha0, p=p, drho=drho, thislt0=mylt0)
		#
		newEvents=[]
		newmags=[]
		# generate events from a branching model
		#NthisEarthquake = int(10**((b/1.1)*(thism - dmstar*1.1 - mc)))	# number of aftershocks for this earthquake.
																							# i think this is approximately right for branching.
		Nsequence = 0	# this will be the number for a single A.S. sequence, excluding sub-sequences.
		#
		Nbass = 10**(b*(mstar - dmstar - mc))	# Nbass new "child" earthquakes
		#print deltaN
		for i in range(int(Nbass)):
			Nsequence+=1
			#thist = invOmori(i, mytau, p, myt0)
			#newEvents+=[[parentTime + thist, [0,0], 0, i, thism]]
			#
			# random magnutide between mc and mstar; does not permit foreshocks and mitigates foreshock
			# and BASS singularity related problems. excluding foreshocks does not appear to significantly alter
			# the model, but it sure does cause problems.
			# use the limited magnitude distribution, but add a foreshock facilitator maxDeltaM
			thism = randmagLim(mstar=(mstar-dmstar + maxDeltaM), mc=mc, b=b, localrand=mrand) 	# note we can add a magitude delta to permit foreshocks
			thist = invOmori((Nbass)*mrand.random(), mytau, p, myt0)
			thisR=invR(N=(Nbass)*(.95*mrand.random()), m=mstar, mc=mc, dmstar=dmstar, b=b, q=q)
			#					
			#newmags+=[thism]
			BASScat += [[parentTime + thist, parentLoc, parentDistSq + thisR**2.0, thisR**2.0, i, thism]]
		#
	
	return BASScat	

def getBASS1(BASScat=[[0, [0,0], 0, 0, 5.0]], mc=3.0, b=1.0, dmstar=1.2, alpha0=alpha_0, doplots=False, limMags=True, timefact=2000.0, drho=0.0, p=1.1, taumethod=0):
	# limMags==True: use limiting-magnitude formula (aka, mc < m <m_star); limmags==False: use unlimited mags, m in (m_c, infty)
	if type(BASScat).__name__ in ('float', 'int'):
		BASScat=[[0, [0,0], 0, 0, BASScat]]
	mrand=random.Random()
	#
	for ev in BASScat:
		mstar=ev[-1]
		thism=mstar
		parentTime=ev[0]
		parentLoc=ev[1]
		#
		# Omori Parameters:
		#myt0=10.0**lDeltat(m=mstar, alpha=alpha0)
		#
		b=1.0
		mydmp=1.0-math.log10(p-1.0)/b		# Delta m prime (from Yoder rupture-duration type omori fomulation)
		#if taumethod>1: taumethod=1	# catch exceptions...
		#
		if taumethod==0:
			# "original" rupture duration method (dt -> Ngr/dt_rupt.)
			myt0=10.0**(lt0(m=mstar, alpha=alpha0))
			mytau = 10**(ltauYoder(m=mstar, mc=mc, p=p, dmp=mydmp))
		if taumethod==1:
			# temporal separation method -- aftershocks begin when omori-interval > rupture-duration
			# (using full sequence parameterization)
			#print "using taumethd 1: %f, %f -- %f, %f" % (lt02(mstar=mstar, alpha0=alpha0, dmstar=1.1, drho=drho), ltau2(mstar=mstar, mc=mc, alpha0=alpha0, p=p, drho=drho, thislt0=None), lt0(m=mstar), ltauYoder(m=mstar, mc=mc, p=p, dmp=mydmp))
			mylt0=lt02(mstar=mstar, alpha0=alpha0, dmstar=dmstar, drho=drho)
			myt0 = 10.**(mylt0)
			mytau = 10.**ltau2(mstar=mstar, mc=mc, alpha0=alpha0, p=p, drho=drho, thislt0=mylt0)
		if taumethod==2:
			# temporal separation method -- aftershocks begin when omori-interval > rupture-duration
			# (using local rate parameterization)
			atmp=0
		#	
		newEvents=[]
		newmags=[]
		# generate events from a branching model
		#NthisEarthquake = int(10**((b/1.1)*(thism - dmstar*1.1 - mc)))	# number of aftershocks for this earthquake.
		#																				# i think this is approximately right for branching.
		Nsequence = 0	# this will be the number for a single A.S. sequence, excluding sub-sequences.
		#
		Nbass = 10**(b*(thism - dmstar - mc))	# Nbass new "child" earthquakes
		deltaN=Nbass
		#print "Nbass: ", Nbass
		for i in range(int(Nbass)):
			Nsequence+=1	# number of events (index of event) in a sequence.
			#
			# random magnutide between mc and mstar; does not permit foreshocks and mitigates foreshock
			# and BASS singularity related problems. excluding foreshocks does not appear to significantly alter
			# the model, but it sure does cause problems.
			#
			if limMags==True:
				thism = randmagLim(mstar=(mstar-dmstar), mc=mc, b=1.0, localrand=mrand) 	# note we can add a magitude delta to permit foreshocks
			# just ranom m. might have foreshocks and foreshock related problems that must be handled (see accompanying code).
			if limMags==False:
				thism = randmag(mc=mc, b=1.0, localrand=mrand)
			#
			# and in the case of unlimited magnitudes, in case we get something too big that's going to make a big mess for us,
			correctionindex=0
			while thism>=3.0*mstar:
				print("correcting for massive eq.,%d: %f/%f" % (correctionindex, thism, mstar)) 
				thism = randmag(mc=mc, b=1.0, localrand=mrand)
				correctionindex+=1
				
			newmags+=[thism]
		#
		#
		magindex=0
		seqLen=len(newmags)
		for nm in newmags:
			#irand=(Nsequence*mrand.random())
			rindex=magindex
			rindex=float(seqLen)*mrand.random()
			#print "Nsequence, irand, mytau, p, myt0",  Nsequence, irand, mytau, p, myt0, mytau/myt0
			#thist = invOmori(irand, mytau, p, myt0)
			thist=None
			# when we permit foreshocks (i think this is a consequence of foreshocks), we sometimes get nonsensical
			# results -- when N*tau(p-1)>t0**(1-p). not sure yet what this means; for now, we either dump those events
			# or reduce irand (effectively the "N" value in the inverted Omori's law) until we get something valid.
			# so what does this mean? some parameters do not permit long sequences? overlapping earthquakes?
			# note: this will never happen if we limit m_rand < mstar.
			while thist==None:
				#print "Neff: %f" % irand
				#irand-=1.0
				#if irand<0: continue
				irand=(Nsequence*mrand.random())
				#thist = timefact*(invOmori(irand, mytau, p, myt0) - myt0)	# imvOmori gives real or "catalog" time; this gives "scaling" time.
				thist = invOmori(irand, mytau, p, myt0) - myt0
				#thist = invOmori(irand, mytau, p, myt0)	# imvOmori gives real or "catalog" time; this gives "scaling" time.
				#print "Nsequence, irand, mytau, p, myt0",  Nsequence, irand, mytau, p, myt0, mytau/myt0
				#continue
			#locVect=vectR(Nsequence, nm, mc=mc, dmstar=dmstar, b=b, q=1.35, xrand=mrand)
			#print magindex
			locVect=vectR2(rindex, mstar, mc=mc, dmstar=dmstar, b=b, q=1.35, xrand=mrand)
			#locVect=vectR(Nsequence*mrand.random(), nm, mc=mc, dmstar=dmstar, b=b, q=1.35, xrand=mrand)
			thisloc=[parentLoc[0] + locVect[0], parentLoc[1]+locVect[1]] 
			newEvents += [[parentTime + thist, thisloc, Nsequence, irand, nm]]
			magindex+=1
		#
		BASScat+=newEvents
	BASScat.pop(0)		# get rid of teh seed event (we want just aftershocks).
	return BASScat

def getDetBranchingCat(BASScat=[[0, [0,0], 0, 0, 5.0]], mc=3.0, b=1.0, dmstar=1.0, B=9.0, alpha0=alpha_0, doplots=False):
	if type(BASScat).__name__ in ('float', 'int'):
		BASScat=[[0, [0,0], 0, 0, BASScat]]
	mrand=random.Random()
	print("basscat len: %d" % (len(BASScat)))
	#
	# for this version, randomize only the magnitude order. there is a way to generalize the whole process from mean values
	# but it appears to bn un-simple.
	#
	# events are [ [time, center, r, n, mag] ] where r is the distance from the "center" location ??, n is the local sequence number
	# 		( aka, the 1st, 2nd, etc. event in its local (parent) sequence).
	# create a deterministic BASS set. mstar might be a fractional magnitude, so the final magnitude m_f might be <m_c.
	# B -> branching parameter (for realy earthquakes, B=9; in summation, this gives GR with base 10.
	#
	# revise thsis structure:
	# pass an initial BASScat = [mstar]
	# from there, use a recursive loop so that the first iteration is identical to the final iteration.
	# the recursive stack should take care of everything from there (aka, the second loop will walk through the first generated catalog).
	#
	#
	# in this case, go back and correct:
	# m_f -> m_c
	# N_f -> N_f * 10^{b(m_f - m_c)}
	#mc=1.5
	#
	#BASScat=[]
	#thism=mstar
	for ev in BASScat:
		#mstar=BASScat[0][0]
		#thism=BASScat[0][0]
		#thisNeq=BASScat[0][1]
		mstar=ev[4]
		thism=mstar
		parentTime=ev[0]
		parentLoc=ev[1]
		
		# Omori Parameters:
		#myt0=10.0**lDeltat(m=mstar, alpha=alpha0)
		myt0=10.0**lt0(m=mstar, alpha=alpha0)
		p=1.1
		b=1.0
		mydmp=1.0-math.log10(p-1.0)/b		# Delta m prime (from Yoder omori fomulation)
		#print "deltaM' %f" % mydmp
		mytau = 10**(ltauYoder(m=mstar, mc=mc, p=p, dmp=mydmp))
		
		newEvents=[]
		newmags=[]
		# generate events from a branching model
		#NthisEarthquake = int(10**((b/1.1)*(thism - dmstar*1.1 - mc)))	# number of aftershocks for this earthquake.
																							# i think this is approximately right for branching.
		Nsequence = 0	# this will be the number for a single A.S. sequence, excluding sub-sequences.

		while thism>=(mc+dmstar):
			bpower=(mstar-thism)/dmstar
			thism-=dmstar
			Nfact=1.0
			if thism<mc:
				Nfact=10**(thism-mc)
				thism=mc
			#
			N=int(Nfact*1.0*(B**bpower))
			for i in range(N):
				Nsequence+=1
				#thist = invOmori(i, mytau, p, myt0)
				#newEvents+=[[parentTime + thist, [0,0], 0, i, thism]]
				newmags+=[thism]
			#
		#
		for nm in newmags:
			irand=(Nsequence*mrand.random())
			thist = invOmori(irand, mytau, p, myt0) - myt0	#invOmori gives actual time; this gives "scaling time"
			locVect=vectR(Nsequence, nm, mc=mc, dmstar=dmstar, b=b, q=1.35, xrand=mrand)
			thisloc=[parentLoc[0] + locVect[0], parentLoc[1]+locVect[1]] 
			newEvents += [[parentTime + thist, thisloc, Nsequence, irand, nm]]
		#
		BASScat+=newEvents
		#while len(newEvents)>0:
		#	irand=int(len(newEvents)*mrand.random())
		#	BASScat+=[newEvents.pop(irand)]
	if doplots:
		a=BASScatPlots(BASScat)
	
	return BASScat

def BASScatPlots(BASScat=None, doCLF=True):
	# some standard plots for a BASS catalog:
	#
	Ts = list(map(operator.itemgetter(0), BASScat))
	Ms = list(map(operator.itemgetter(4), BASScat))
	
	# return [Ts, Ms]
	#
	#if doplots:
	plt.ion()
	plt.figure(0)
	if doCLF: plt.clf()
	plt.semilogx(Ts, Ms, 'o')
	plt.ylabel('magnitude $m$')
	plt.xlabel('time $t$')
	#
	
	locs=list(map(operator.itemgetter(1), BASScat))
	X=list(map(operator.itemgetter(0), locs))
	Y=list(map(operator.itemgetter(1), locs))
	#
	#if doplots:
	plt.figure(1)
	if doCLF: plt.clf()
	plt.plot(X,Y, ',')
	plt.xlabel('longitude ($km$)')
	plt.ylabel('latitude ($km$)')
	#
	BASScat.sort(key=operator.itemgetter(0))
	#
	myints=getints(BASScat)
	print(len(myints[0]), len(BASScat))
	# just get our own intervals:
	#intsT=[BASScat[0][0]]
	#intsInt=[0]
	#for i in range(1, len(BASScat)):
	#	thisint=BASScat[i][0]-BASScat[i-1][0]
	#	intsT+=[BASScat[i][0]]
	#	print thisint, BASScat[i][0]
	#	intsInt+=[thisint]	
	
	#if doplots:
	plt.figure(2)
	if doCLF: plt.clf()
	plt.loglog(myints[0], myints[1], ',')
	#plt.plot(intsT, intsInt, '-')
	plt.ylabel('interval duration $t_0$ (seconds)')
	plt.xlabel('time $t$ (seconds)')
	
	#
	#if doplots:
	# and a cumulative time distribution:
	cdfT=list(map(operator.itemgetter(0), BASScat))
	cdfN=list(range(1, len(cdfT)+1))
	#
	plt.figure(3)
	if doCLF: plt.clf()
	plt.loglog(cdfT, cdfN, '-o')
	plt.xlabel('time $t$ (seconds)')
	plt.ylabel('Total Number of Earthquakes $N$')
	
	#return BASScat
	return None

def xyfromR(R):
	#x=2.0*mrand.random() - 1.0
	#y=(1.0-x*x)**.5
	#if mrand.random()>=.5:
	#	y=-y
	theta = 2.0*math.pi*mrand.random()
	x=math.cos(theta)
	y=math.sin(theta)
	#	
	newx=x*R
	newy=y*R
	#
	return [newx, newy]

def plotBASSLocs2(BASScat=None, fignum=3, maxfact=100.0):
	plt.ion()
	f=plt.figure(fignum)
	plt.clf()
	#
	#if maxfact!=None:
	#	maxr=maxfact*10.0**lruptLen(m, dl0=1.0)
	
	centers=list(map(operator.itemgetter(1), BASScat))
	locs=[]
	Rs = list(map(operator.itemgetter(3), BASScat))
	Rs0 = list(map(operator.itemgetter(2), BASScat))
	#
	for i in range(len(Rs)):
		r2=Rs[i]**.5
		r1=Rs0[i]**.5
		
		thisxy1=xyfromR(r1)
		thisxy2 = xyfromR(r2)
		
		#
		locs+=[[centers[i][0] + thisxy1[0] + thisxy2[0], centers[i][1] + thisxy1[1] + thisxy2[1]]]
		#
	#
	X=list(map(operator.itemgetter(0), locs))
	Y=list(map(operator.itemgetter(1), locs))
	#
	#plt.clf()
	plt.plot(X,Y, ',')
	#
	return None

def plotBASSLocs(BASScat=None, fignum=0, maxfact=100.0):
	plt.ion()
	f=plt.figure(fignum)
	plt.clf()
	#
	#if maxfact!=None:
	#	maxr=maxfact*10.0**lruptLen(m, dl0=1.0)
	
	locs=list(map(operator.itemgetter(1), BASScat))
	X=list(map(operator.itemgetter(0), locs))
	Y=list(map(operator.itemgetter(1), locs))
	#
	#plt.clf()
	plt.plot(X,Y, ',')
	#
	return None

def getints(incat=None):
	if incat==None:
		incat=getDetBranchingCat()
	theseints=[[],[]]
	for i in range(1, len(incat)):
		theseints[0]+=[incat[i][0]]
		theseints[1]+=[incat[i][0]-incat[i-1][0]]
	#
	return theseints

def invOmori(N, tau, p, t0):
	#print t0, tau, t0**(1-p)
	#if (N*tau*(1.0-p) + t0**(1.0-p))<0:
	#	return None
	# 
	# "scaling time", aka, 0 is at the end of the earthquake where scaling starts (lop off eq. duration)
	#t=(N*tau*(1.0-p) + t0**(1.0-p))**(1.0/(1.0-p)) - t0
	#
	# just plain time: t=0 indicates begining or middle of earthquake (aka, we don't lop off the eq. duration).
	q=1.0-p
	#t=(N*tau*(1.0-p) + t0**(1.0-p))**(1.0/(1.0-p))
	t=(N*tau*q + t0**q)**(1.0/q)
	#print t0, tau, t0**(1-p), t
	return t

def invGR(N, mc=3.0, dmstar=1.0, b=1.0):
	# this equation only makes sense when N is the total number of earthquakes and the returned
	# value equals the the mag. of the mainshock (mstar) of the sequence.
	# aka, given a full sequence length N and mc, dmstar, what is the magnitude of the parent event?
	# this can also be used to generate a set of magnitudes, but understand that it is sort of inverted.
	# aka, "if we have N events, the largest magnitude is probably...". note then that for N=1, the 
	# expected largest magnitude is mc=dmstar + log(N)/b. so, if using this function for this purpose,
	# use dmstar=0 to get just a GR distribution from the bottom up. an alternative formulation is given below.
	#
	# this will give all the small earthquakes first, so in a sense, it is more forgiving thatn invGR2().
	mstar = dmstar + mc + (math.log10(N)/b)
	return mstar

def invGR2(N, mstar=5.0, mc=3.0, b=1.0):
	# This gives the rank-order magnitudes from top to bottom. aka, for N=1, m=mstar. for n=N*, 
	# this can be used to generate a sequence of earthquakes.
	# note that by design, N {1, 10**(mstar-mc)
	#
	m = mstar - (math.log10(N))/b
	#
	return mstar

def vectR2(N, m, mc=3.0, dmstar=1.0, b=1.0, q=1.35, xrand=None):
	# vector from radius R:
	if xrand==None:
		xrand=random.Random()
	R=invR(N=N, m=m, mc=mc, dmstar=dmstar, b=b, q=q)
	#print R
	x=2.0*xrand.random() - 1.0 #random.random()
	theta=math.acos(x)
	yparity=1.0
	if xrand.random()>=.5:
		yparity=-1.0
	y=yparity*math.sin(theta)
	#
	rvect=[R*x, R*y]
	#
	return rvect
	

def vectR(N, m, mc=3.0, dmstar=1.0, b=1.0, q=1.35, xrand=None):
	# vector from radius R:
	if xrand==None:
		xrand=random.Random()
	R=invR(N=N, m=m, mc=mc, dmstar=dmstar, b=b, q=q)
	#print R
	x=R*xrand.random() #random.random()
	#y=R*(1.0-x*x)**.5
	y=(R*R - x*x)**.5
	#y=R*xrand.random()
	#x=R*x
	yparity=1.0
	if random.random()>.5:
		yparity=-1.0
	xparity=1.0
	if random.random()>.5:
		xparity=-1.0
	#
	rvect=[x*xparity, y*yparity]
	#print rvect, (rvect[0]**2 + rvect[1]**2)/(R*R)
	
	return rvect

	
def checkrandmag(N=10000):
	mstar=6.0
	mc=2.0
	b=1.0
	plt.figure(6)
	plt.clf()
	plt.ion()
	
	localrand=mrand
	#
	for j in range(10):
		Ms=[]
		for i in range(N):
			Ms+=[randmagLim(mstar, mc, b, localrand)]
		Ns=list(range(1, len(Ms)+1))
		Ms.sort()
		Ns.reverse()
		#
		plt.semilogy(Ms, Ns, ',')
		plt.semilogy(Ms, Ns, '-', alpha=.25)

def randmag(mc=3.0, b=1.0, localrand=mrand):
	m = mc - (1.0/b)*math.log10(mrand.random())
	return m
	
def randmagLim(mstar=None, mc=3.0, b=1.0, localrand=mrand):
	# an upper limited GR distribution, to avoid the singularity problem.
	# all magnitudes are m<mstar	
	nmax=10**(b*(mstar-mc))
	randN=math.ceil(nmax*mrand.random())
	#print nmax, randN
	rmag = mstar - ((1.0/b)*math.log10(randN))
	#rmag = mc - (1.0/b)*math.log10(mrand.random())
	
	return rmag
	
def invR(N, m, mc=3.0, dmstar=1.0, b=1.0, q=1.35, ruptfact=0.0):
	# (note that we use a random number generator to fix  a weird problem, so this is only good for ETAS sims.
	# and related applications).
	#
	# derrived from Felzer and Broadsky (2006), this gives us an estimated (max?) distance from the epicenter.
	# the rupture length is from Kagan 2002 and Yoder et al. 2012a. 
	# note that if we provide N from a random distribution, we get a random (radial) position generator.
	lrupt=10.0**(lruptLen(m, dl0=1.0) + ruptfact)
	#r=lrupt*(N*10.0**(b*(mc + dmstar - m)) +1.0)**(-1.0/(1.0-q))
	# minor variation (see notes) in linear density.
	Nomori = 10.0**(b*(m-dmstar-mc))
	qexp=(1.0/(1.0-q))
	#print (1.0-float(N)/Nomori)
	#print math.log10(Nomori), math.log10(N)
	#print Nomori, N
	
	if (1.0 - (float(N)/Nomori) )<0:
		#print (1.0 - (float(N)/Nomori) ), (float(N)/Nomori), N, Nomori
		# and just fix it:
		N=Nomori*mrand.random()
	
	#r = lrupt*((1.0-float(N)/Nomori)**qexp - 1.0)
	r = lrupt*(1.0 - (float(N)/Nomori) )**qexp

	return r

def lruptLen(m, dl0=1.0):
	return (m/2.0) - 1.755 - dl0

def forecastdemo2(gridsize=.1, contres=5, thisdt=None, eqtheta=0, eqeps=1.0):
	# and let's develop a more versitile script that outsources the contours.
	#resolution = 5
	# looks like this is defunct after we creaded the BASScast module...
	#
	plt.figure(1)
	plt.clf()
	plt.ion()
	if thisdt==None: thisdt=dtm.datetime.now(pytz.timezone('UTC'))
	ql=getTestQuakes2(floatnow=mpd.date2num(thisdt)*days2secs - 0.0)
	bc2=1.0
	bc2=None
	bc2=bcp.BASScast(incat=ql, mc=2.0, gridsize=gridsize, contres=contres, fcdate=thisdt, eqtheta=eqtheta, eqeps=eqeps)

	return bc2

def getMFETAScatFromANSS(lons=[-121.0, -114.0], lats=[31.0, 37.0], dates=[None, None], mc=4.0):
	if dates==None: dates=[None, None]
	if dates[1]==None: dates[1]=dtm.datetime.now(pytz.timezone('UTC'))
	if dates[0]==None or dates[0]>=dates[1]: dates[1]-dtm.timedelta(days=840)
	#
	clist1=atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=dates, Nmax=9999999, fout=None)
	catalog=[]
	#X,Y,M, Z = [], [], [], []
	#
	for rw in clist1:
		catalog+=[[mpd.date2num(rw[0]), rw[1], rw[2], rw[3], rw[4]]]
		#X+=[rw[2]]	# lon
		#Y+=[rw[1]]	# lat
		#M+=[rw[4]]
	return catalog
	
def forecastdemo3(gridsize=.1, contres=5, lons=[-121.0, -114.0], lats=[31.0, 37.0], dates=[None, None], mc=4.0, fromfname=None):
	# fromfname='cats/socal.cat' should work...
	# get a catalog from ANSS and apply a forecast.
	# plot plain and on a map.
	# dtm.datetime.now(pytz.timezone('UTC'))-dtm.timedelta(days=60), dtm.datetime.now(pytz.timezone('UTC')
	dt1=dtm.datetime.now(pytz.timezone('UTC'))-dtm.timedelta(days=int(5*360))
	dt2=dtm.datetime.now(pytz.timezone('UTC')) + dtm.timedelta(days=1)
	#dt1=dtm.date(dt1.year, dt1.month, dt1.day)
	#dt2=dtm.date(dt2.year, dt2.month, dt2.day)
	#print dt1, dt2
	#
	if fromfname!=None:
		thiscat=ypp.eqcatalog()
		thiscat.loadCatFromFile(fromfname)
		clist1=[]
		for i in range(len(thiscat.cat)):
			if thiscat.cat[i][0]<dt1 or thiscat.cat[i][3]<mc: continue
			clist1+=[thiscat.cat[i]]
			clist1[-1].insert(3,None)	# so some later stuff will work.
	else:
		clist1=ypp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[dt1, dt2], Nmax=999999, fout='tmpcat.cat')
		#
	#
	print(" catlen: %d" % len(clist1))
	# [ [dtm, lat, lon, depth(?), mag, ...], ...]
	catalog=[]
	X,Y,M = [], [], []
	plt.figure(0)
	plt.clf()
	
	maxNquakes=250.0	# max num. equakes to plot (approximately)
	Nquakes = float(len(clist1))
	#
	mthresh=math.log10(Nquakes/maxNquakes) + mc	# we can take  out the "if" bc for N<Nmax, mthresh<mc
																# (plot some max. number of earthquakes)
	print("spinning clist1.")
	for rw in clist1:	
		thisdt=rw[0]
		if type(thisdt) == type('abc'): thisdt = mpd.datestr2num(thisdt)
		#catalog+=[[mpd.datestr2num(rw[0]), rw[1], rw[2], rw[4]]]
		catalog+=[[thisdt, rw[1], rw[2], rw[4]]]
		X+=[rw[2]]	# lon
		Y+=[rw[1]]	# lat
		M+=[rw[4]]
		#
		if rw[4]<mthresh: continue
		plt.plot([rw[2]], [rw[1]], 'o', ms=1*rw[4], alpha=.7, zorder=10)
		#
	#
	print("getting basscast...")
	#return catalog
	bc3=bcp.BASScast(incat=catalog, gridsize=gridsize, contres=contres, mc=mc-2.0)
	plt.figure(1)
	plt.clf()
	print("doing map...")
	bc3.BASScastContourMap(fignum=1, X_i=bc3.X_i, Y_i=bc3.Y_i, Z2d=bc3.Z2d)
	#plt.plot(X,Y, 'o', ms=5 )
	kmstr1=bc3.KMLstrFromConts()
	fout=open('tmp/kml1.kml', 'w')
	fout.write(kmstr1)
	fout.close()
	#
	km2 = kmlConts(X_i=bc3.X_i, Y_i=bc3.Y_i, Z_i=bc3.Z2d)
	
	#
	#cnts = bc3.getContourSet(X_i=bc3.X_i, Y_i=bc3.Y_i, Z2d=bc3.Z2d, 25)
	#
	return bc3

def kmlConts(X_i, Y_i, Z_i):
#def kmlConts(cset=None, colorbarname=None):
	#-----------------------------------------------------------------------------#

	resolution = 5

	LevelsNumber = 5 * resolution
	warnings = ['No','Low','Guarded','Elevated','High','Severe']
	
	#-----------------------------------------------------------------------------#

	# Create your data here
	#
	# I assume there are three arrays: X_i, Y_i, and Z_i

	#-----------------------------------------------------------------------------#

	# Create the KML file
	cs = plt.contourf(X_i, Y_i, Z_i, LevelsNumber,
		               cm=plt.spectral()).collections
	
	#cs=cset
	
	file = open('tmp/kml2.kml','w')
	file.write('<?xml version="1.0" encoding="UTF-8"?>\n')
	file.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
	file.write('<Document>\n')

	for i in range(0,len(cs)):
		 bgra_array = 255*cs[i].get_facecolor()[0]
		 file.write('<Style id="l%d">\n' % i)
		 file.write('<LineStyle><color>00000000</color></LineStyle>\n')
		 file.write('<PolyStyle><color>7d%02x%02x%02x</color></PolyStyle>\n' %\
		                ( bgra_array[2] , bgra_array[1] , bgra_array[0] ) )
		 file.write('</Style>\n')

	file.write('<ScreenOverlay id="scale">\n')
	file.write('<name>Color Scale</name>\n')
	file.write('<Icon><href>scale.png</href></Icon>\n')
	file.write('<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n')
	file.write('<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n')
	file.write('<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n')
	file.write('</ScreenOverlay>\n')

	for i in range(resolution,len(cs)):
		 file.write('<Placemark>\n')
		 file.write('<name>%s Risk</name>\n' % warnings[i/resolution])
		 file.write('<styleUrl>#l%d</styleUrl>\n' % i)
		 file.write('<MultiGeometry>\n')
		 
		 for trace in cs[i].get_paths():
		     file.write('<Polygon>\n')
		     file.write('<extrude>0</extrude>\n')
		     file.write('<altitudeMode>clampToGround</altitudeMode>\n')
		     file.write('<outerBoundaryIs>\n')
		     file.write('<LinearRing>\n')
		     file.write('<coordinates>\n')

		     for lng,lat in trace.vertices:
		         file.write('%s,%s,0\n' % (round(lng,3),round(lat,3)))

		     file.write('</coordinates>\n')
		     file.write('</LinearRing>\n')
		     file.write('</outerBoundaryIs>\n')
		     file.write('</Polygon>\n')

		 file.write('</MultiGeometry>\n')
		 file.write('</Placemark>\n')

	file.write('</Document>\n')
	file.write('</kml>')
	file.close()

	fig = plt.figure(figsize=(1,3))
	axe = fig.add_axes([0.05, 0.05, 0.2, 0.9])

	fig.figurePatch.set_alpha(0)

	matplotlib.rc('font', family='monospace', weight='black')

	cmap = mpl.cm.spectral
	norm = mpl.colors.Normalize(vmin=Z_i.min(), vmax=Z_i.max())
	tics = [0,10,20,40,80]

	#cb1 = mpl.colorbar.ColorbarBase(axe, normorm, ticks=tics, format="%g%%")
	#?
	#cb1 = mpl.colorbar.ColorbarBase(axe, norm=mpl.colors.NoNorm, ticks=tics, format="%g%%")
	cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%g%%")
	cb1.set_label('Probability')

	plt.savefig('scale.png')

	#-----------------------------------------------------------------------------#

	os.system('zip -q forecast.kmz doc.kml scale.png')
	os.system('rm -rf doc.kml scale.png')

###########

def getContourSet(X_i, Y_i, Z_ij, contres):
	# X_i, Y_i are x,y coordinates (aka, like range(100), range(200) ), Z_ij is a 2D array that is actually contoured.
	# contres is the resolution of the contoursl; numConts -> 5*contres
	#
	resolution = contres
	LevelsNumber = 5 * resolution
	warnings = ['No','Low','Guarded','Elevated','High','Severe']
	#
	# cnts=plt.contourf(gridsize*numpy.array(range(nelements)), gridsize*numpy.array(range(nelements)), Z2d,10)
	# retrieve the collections() object from the contourf() function which returns a matplotlib.contour.ContourSet
	#
	#cs = plt.contourf(X_i, Y_i, Z_ij, LevelsNumber, cm=plt.spectral()).collections
	cs = plt.contourf(X_i, Y_i, Z_ij, LevelsNumber, cm=plt.spectral())
	
	return cs

def arrayFromSet(cs):
	# get an array (list) from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
	# format:
	# [ [z0, [[x,y], [x,y],...], [z1, [paths1]], 
	#
	levels=cs.levels	# lower bound of contours
	layers=cs.layers	# upper bound of contours (at least for contours <0)
	dz= layers[0] - levels[0]
	collects=cs.collections
	carray = []
	for i in range(0, len(collects)):
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

def KMLstrFromConts(cset):
	# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
	# (will it be necessary to write directly to file? is there a string length limit problem?)
	#
	cs=cset.collections
	#resolution = 5
	#LevelsNumber = 5 * resolution
	warnings = ['No','Low','Guarded','Elevated','High','Severe']
	#contoursStart=int(.2*len(cs))
	resolution = int(len(cs)/5)
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
	kmlstr+='<Icon><href>scale.png</href></Icon>\n'
	kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
	kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
	kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
	kmlstr+='</ScreenOverlay>\n'
	#
	for i in range(resolution, len(cs)):
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

			for lng,lat in trace.vertices:
				kmlstr+='%s,%s,0\n' % (round(lng,3),round(lat,3))

			kmlstr+='</coordinates>\n'
			kmlstr+='</LinearRing>\n'
			kmlstr+='</outerBoundaryIs>\n'
			kmlstr+='</Polygon>\n'

		kmlstr+='</MultiGeometry>\n'
		kmlstr+='</Placemark>\n'
	#
	kmlstr+='</Document>\n'
	kmlstr+='</kml>'
	#
	return kmlstr	

def getTestQuakes(mc=2.0, floatnow=None):
	quakes=[]
	if floatnow==None: floatnow=mpd.date2num(dtm.datetime.now(pytz.timezone('UTC')))*days2secs
	#for i in range(Neq):
	#	thismag=mc + 2.0 + 4.0*mrand.random()
	#	quakes += [earthquake(mag=thismag, loc=[X*mrand.random(), Y*mrand.random()], evtime=1.0, mc=mc)]
	#	#plt.plot([quakes[-1].loc[0]], [quakes[-1].loc[1]], 'o', ms=2*quakes[-1].mag)
	quakes+=[bcp.earthquake(mag=4.0, loc=[6.5, 9.], evtime=floatnow, mc=mc)]
	
	quakes+=[bcp.earthquake(mag=4.0, loc=[.5, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=4.5, loc=[.7, 8.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=5.0, loc=[1.5, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=5.5, loc=[2.0, 6.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=6.0, loc=[2.5, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=6.5, loc=[3.5, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=7.0, loc=[5.0, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=8.0, loc=[7.0, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=8.5, loc=[8.0, 7.], evtime=floatnow-1.0, mc=mc)]
	#
	testmag=7.0
	testquake = bcp.earthquake(mag=testmag, loc=[0.,0.],evtime=floatnow, mc=mc)
	thist=10**(testquake.lrupttime)
	
	print(floatnow, thist)
		 
	#print thist
	quakes+=[bcp.earthquake(mag=testmag, loc=[2.0, 3.], evtime=floatnow+1., mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[3.0, 3.], evtime=floatnow-0.1*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[3.0, 4.], evtime=floatnow-0.5*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[4.0, 3.], evtime=floatnow-1.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[5.0, 3.], evtime=floatnow-5.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[6.0, 3.], evtime=floatnow-10.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[7.0, 3.], evtime=floatnow-100.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[8.0, 3.], evtime=floatnow-5000.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[9.0, 3.], evtime=floatnow-10000.0*thist, mc=mc)]	
	
	quakes+=[bcp.earthquake(mag=testmag, loc=[2.0, 2.], evtime=floatnow-100.0*thist, mc=mc)]
	#
	qlist=[]
	for q in quakes:
		qlist += [[q.eventftime, q.loc[1], q.loc[0], q.mag]]
		#
	#
	#return quakes
	return qlist

def getTestQuakes2(mc=2.0, floatnow=None):
	quakes=[]
	if floatnow==None: floatnow=mpd.date2num(dtm.datetime.now(pytz.timezone('UTC')))*days2secs
	#for i in range(Neq):
	#	thismag=mc + 2.0 + 4.0*mrand.random()
	#	quakes += [earthquake(mag=thismag, loc=[X*mrand.random(), Y*mrand.random()], evtime=1.0, mc=mc)]
	#	#plt.plot([quakes[-1].loc[0]], [quakes[-1].loc[1]], 'o', ms=2*quakes[-1].mag)
	#quakes+=[bcp.earthquake(mag=4.0, loc=[6.5, 9.], evtime=floatnow, mc=mc)]
	
	quakes+=[bcp.earthquake(mag=2.0, loc=[.5, .5], evtime=floatnow-1.0, mc=mc)]
	'''
	quakes+=[bcp.earthquake(mag=4.5, loc=[.7, 8.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=5.0, loc=[1.5, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=5.5, loc=[2.0, 6.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=6.0, loc=[2.5, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=6.5, loc=[3.5, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=7.0, loc=[5.0, 7.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=8.0, loc=[7.0, 7.], evtime=floatnow-1.0, mc=mc)]
	'''
	quakes+=[bcp.earthquake(mag=2.0, loc=[8.0, 8.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=8.0, loc=[4.0, 4.], evtime=floatnow-1.0, mc=mc)]
	
	quakes+=[bcp.earthquake(mag=8.0, loc=[2.0, 2.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=8.0, loc=[3.0, 3.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=7.0, loc=[5.0, 5.], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=7.0, loc=[7.0, 7.], evtime=floatnow-1.0, mc=mc)]
	
	m0=3.0
	quakes+=[bcp.earthquake(mag=m0, loc=[2.125, 2.125], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=m0, loc=[2.25, 2.25], evtime=floatnow-1.0, mc=mc)]
	quakes+=[bcp.earthquake(mag=m0, loc=[2.375, 2.375], evtime=floatnow-1.0, mc=mc)]	
	quakes+=[bcp.earthquake(mag=m0, loc=[2.5, 2.5], evtime=floatnow-1.0, mc=mc)]	
	quakes+=[bcp.earthquake(mag=m0, loc=[2.625, 2.625], evtime=floatnow-1.0, mc=mc)]	
	quakes+=[bcp.earthquake(mag=m0, loc=[2.75, 2.75], evtime=floatnow-1.0, mc=mc)]	

	
	#
	testmag=7.0
	testquake = bcp.earthquake(mag=testmag, loc=[-1.,-1.],evtime=floatnow, mc=mc)
	testquake = bcp.earthquake(mag=testmag, loc=[-.5,-.5],evtime=floatnow, mc=mc)
	testquake = bcp.earthquake(mag=testmag, loc=[0.,0.],evtime=floatnow, mc=mc)
	testquake = bcp.earthquake(mag=testmag, loc=[1.,1.],evtime=floatnow, mc=mc)
	testquake = bcp.earthquake(mag=testmag, loc=[.5,.5],evtime=floatnow, mc=mc)
	thist=10**(testquake.lrupttime)
	
	'''
	print floatnow, thist
		 
	#print thist
	quakes+=[bcp.earthquake(mag=testmag, loc=[2.0, 3.], evtime=floatnow+1., mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[3.0, 3.], evtime=floatnow-0.1*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[3.0, 4.], evtime=floatnow-0.5*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[4.0, 3.], evtime=floatnow-1.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[5.0, 3.], evtime=floatnow-5.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[6.0, 3.], evtime=floatnow-10.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[7.0, 3.], evtime=floatnow-100.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[8.0, 3.], evtime=floatnow-5000.0*thist, mc=mc)]
	quakes+=[bcp.earthquake(mag=testmag, loc=[9.0, 3.], evtime=floatnow-10000.0*thist, mc=mc)]	
	
	quakes+=[bcp.earthquake(mag=testmag, loc=[2.0, 2.], evtime=floatnow-100.0*thist, mc=mc)]
	#
	'''
	qlist=[]
	for q in quakes:
		qlist += [[q.eventftime, q.loc[1], q.loc[0], q.mag]]
		#
	#
	#return quakes
	return qlist

def forecastdemo1(gridsize=.1, Neq=10, contres=2):
	# and this will get more versitile as we go...
	#
	mc=2.0
	plt.ion()
	plt.figure(0)
	plt.clf()
	#
	xmax=10
	# set up a bunch of forecastsite() objects.
	nelements = int(xmax/gridsize)	# not actually the number of elements; the the number of rows/cols
	#X=100.0*gridsize
	#Y=100.0*gridsize
	X=float(xmax)
	Y=float(xmax)
	grid=[]
	floatnow=mpd.date2num(dtm.datetime.now(pytz.timezone('UTC')))
	for i in range(nelements*nelements):
		grid+=[forecastsite(loc=[gridsize*(i%nelements), gridsize*(i/int(nelements))], dxdy=[gridsize, gridsize], evtime=floatnow, mc=2.0)]
	#
	# and make some earthquakes:
	quakes=getTestQuakes(floatnow=floatnow)
	#
	# now, get intensities and make a contour plot.
	# for now, plot intensities (and related calcs) here. eventually, these should be forecastsite() member functions.
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
	#Z=map(operator.itemgetter(2), cdata)
	Z2d=scipy.array(Z)
	#Z2d.shape=(len(Y), len(X))
	Z2d.shape=(nelements, nelements)
	#X1, Y1=numpy.meshgrid(X, Y)
	
	thisXs=numpy.array(list(range(nelements)))
	thisXs*=gridsize
	#thisYs=thisXs
	
	X_i = gridsize*numpy.array(list(range(nelements)))
	Y_i = gridsize*numpy.array(list(range(nelements)))
	# cnts=plt.contourf(X_i, Y_i, Z2d,10)
	cnts = getContourSet(X_i, Y_i, Z2d, 5*contres)
	#
	# get contours without plotting?
	#import matplotlib._cntr as cntr
	#c2= cntr.Cntr(thisXs, thisXs, numpy.array(Z2d))
	#plt.plot(X,Y, 's', alpha=.2, ms=10)
	
	#return [grid, quakes]
	#return c2
	return cnts

def mag5test(catname='cali5.cat', dorefresh=False, mc=2.5, mtest=5.0, drho=5.0, p=1.1):
	if catname==None:
		dorefresh=True
		catname='cali5.cat'
	if dorefresh==True:
		#print "script catalog fetch."
		# lon=[135.0, 150.0]	lat=[30.0, 41.5]	m0=4.000000	dates=[datetime.date(2005, 1, 1), datetime.date(2011, 7, 12)]
		hmcat = bcp.atp.catfromANSS(lon=[-110.0, -126.0], lat=[30.0, 42.5], minMag=mc, dates0=[dtm.datetime(2005, 1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], Nmax=9999999, fout=catname)
		cpf2=ypp.eqcatalog(hmcat)
	else:
		cpf2=ypp.eqcatalog()
		cpf2.loadCatFromFile(fname=catname, minmag=mc)
	#
	# just in case:
	cpf2.cat.sort(key = lambda x: x[0])
	cpf=ypp.eqcatalog(cpf2.getMagSubcat(minmag=mtest, fullcat=cpf2.getcat(0)))
	cpf.cat.sort(key = lambda x: x[3])
	cpf.cat.reverse()
	
	#
	plt.figure(1)
	plt.clf()
	plt.ion()
	cpf.plotCatMap(fignum=1)
	cpf.catmap.drawstates()
	cpf.catmap.drawcountries()
	#
	maincat=cpf.getcat(0)
	simplecat=[]
	for ev in maincat:
		x,y=cpf.catmap(ev[2], ev[1])
		t=mpd.date2num(ev[0])*60.*60.*24.
		#t=mpd.date2num(ev[0])
		simplecat+=[[t,y,x,ev[3]]]
		plt.plot([x], [y], '*', ms=5+5*(ev[3]-5))
		#
		# now, draw the sub-cat square (and eventually fetch the subcat).
		Lkm = 10**(ev[3]/2.0 - 1.655)
		latrads=ev[1]*2.0*math.pi/360.
		Lm=3.*Lkm*1000.
		#
		#Ldlon = Lkm*111./math.cos(latrads)
		#ldlat = Lkm*111.
		X=[x-Lm/2., x+Lm/2., x+Lm/2., x-Lm/2., x-Lm/2.]
		Y=[y-Lm/2., y-Lm/2., y+Lm/2., y+Lm/2., y-Lm/2.]
		plt.plot(X,Y, 'b.--')
	
	# decluster:
	thiscat = simplecat
	#keepquakes=[thiscat[0] + [0]]	# earthquakes to keep for analysis.
	# for now, let's keep:
	# 1) any lone earthquakes (not inside the bounds of any other eq.)
	# 2) the largest of any overlapping eqs. ? this is actually a bit harder. we need
	#     to make the initial list and then iterate over that list until overlappers are gone.
	keepquakes=[]
	for i in range(0, len(thiscat)):
		useeq=True
		#for ii in range(i):
		for ii in range(0, len(thiscat)):
			# is this earthquake in the catalog of a previous earthquake? aka, is it an aftershock (is the catalog contaminated)?
			dx=thiscat[ii][2] - thiscat[i][2]
			dy=thiscat[ii][1] - thiscat[i][1]
			mii, mi = thiscat[ii][3], thiscat[i][3]
			dtii, dti = thiscat[ii][0], thiscat[i][0]
			rsq=dx*dx + dy*dy
			#
			Lm = 3.*10**(3.0 + thiscat[ii][3]/2.0 - 1.655)	# approximate "range" of an earthquake (this could be calucaled
																				# more precisely using an omori-like distribution.
			#
			# ths i event occurred after the ii events. here we would nominally quantify how long after.
			#
			#
			# is it out of spatial range?
			#if Lm*Lm>rsq and mi<=mii:
			#if rsq<Lm*Lm and dti>dtii and mi<mii:
			if rsq<Lm*Lm and mi<mii:
				# keeping the biggest, independent earthquakes.
				useeq=False
				continue
			#
		if useeq==True:
			keepquakes += [thiscat[i] + [i]]
			plt.plot([keepquakes[-1][2]], [keepquakes[-1][1]], 'o', ms=15, alpha=.25)
			#print keepquakes[-1]
		#
	#
	print("{ %d	} distinct events." % len(keepquakes))
	#
	#catalogs=[]
	for rw in keepquakes:
		Lm = 1.5*10**(3.0 + rw[3]/2.0 - 1.655)
		x1, x2 = rw[2]-Lm, rw[2]+Lm
		y1, y2 = rw[1]-Lm, rw[1]+Lm
		lon1, lat1 = cpf.catmap(x1, y1, inverse=True)
		lon2, lat2 = cpf.catmap(x2, y2, inverse=True)
		#print lon1, lat1, lon2, lat2, maincat[rw[4]][0]
		#
		# getxytSubcat(self, fullcat=None, dts=[], lats=[], lons=[], llcols=[1,2]
		# addxytSubcat(self, subcatname='xytsubcat', fullcat=None, dts=[], lats=[], lons=[], llcols=[1,2])
		#cpf2.addxytSubcat(subcatname='i%d-m%.2f' % (rw[4], rw[3]), fullcat=cpf2.getcat(0), dts=[maincat[rw[4]][0], dtm.datetime.utcnow()], lats=[lat1, lat2], lons=[lon1, lon2])
		cpf2.addxytSubcat(subcatname='%f' % (rw[3]), fullcat=cpf2.getcat(0), dts=[maincat[rw[4]][0], dtm.datetime.now(pytz.timezone('UTC'))], lats=[lat1, lat2], lons=[lon1, lon2])
		#cpf2.subcats+=['i%d-m%.2f' % (rw[4], rw[3]), cpf2.getxytSubcat(fullcat=cpf2.getcat(0), dts=[maincat[rw[4]][0], dtm.datetime.now(pytz.timezone('UTC'))], lats=[y1, y2], lons=[x1, x2])
	#
	# lt02 (mstar, alpha0=alpha_0, dmstar=1.2, drho=0.0)
	#  ltau2(mstar=None, mc=None, alpha0=alpha_0, p=1.1, drho=0.0, thislt0=None)
	# t0=lt02(mstar=
	# start with absolutely no overlapping...
	#
	for sc in cpf2.subcats:
		#print sc[0], len(sc[1])
		ct=sc[1]
		Xp=scipy.array(list(map(operator.itemgetter(2), ct)))
		Yp=scipy.array(list(map(operator.itemgetter(1), ct)))
		X,Y = cpf.catmap(Xp, Yp)
		mainmag=cpf2.getMainEvent(ct)[3]
		plt.plot(X,Y, '.', ms=2, alpha=.5, zorder=round(mainmag))
	# eqintervals(cat=None, mstar=None, mc=None, dmstar=1.0, mev=None, catnum=1, alpha0=None, omoriWaitExp=None, effmstar=None, drho=5.0, p=1.1)
	return [cpf2, cpf]
	

def plotFitFile(fname, fnum=0):
	X=[]
	Y=[]
	mc, beta, y = None, None, None
	f=open(fname)
	plt.figure(fnum)
	plt.clf()
	plt.ion()
	#
	for rw in f:
		if rw[0]=='#': continue
		rws=list(map(float, rw.split()))
		#print rws
		#print beta
		#
		if beta==None or beta!=rws[1]:
			if len(X)>0: 
				plt.plot(X,Y, label='$\\beta=%f$' % beta)
				plt.legend(loc='best')
			#print beta, rws[1]
			X, Y = [], []
			#if beta==None: continue
			
		X+=[rws[0]]
		Y+=[rws[2]]
		beta = rws[1]
		#
	f.close()
		
		
		
		

		
		
		
		
		
		
	
