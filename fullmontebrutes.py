from yodabass import *
import copy

def logNOmoriScaling(m, mc, t, p=1.08, alpha=-.5, beta=-.5, gamma=9.0, dm=1.0, b=1.0):
	# note: we can use this function exactly for r-scaling with r0=t0, chi=tau, etc.
	t0=10**(alpha*m + beta*mc + gamma)
	#print "val lNomori: ", (1.0 + t/t0)**(1.0-p)
	if (1.0 + t/t0)**(1.0-p)>=1.0: return None
	#N=b*(m-dm-mc) + math.log10(p-1.0) + math.log10( (1.0 + t/t0)**(1.0-p) - 1.0)
	N=b*(m-dm-mc) + math.log10( 1.0 - (1.0 + t/t0)**(1.0-p))
	return N

def brutefitres(fitset='parkfield',mstar=6.0, mc=2.5, p=1.08, b=1.0, dm=1.0, alpha0=None, beta0=None, gamma0=None, X=None, Y=None):
	# "results queue", "process queue"
	# (for omori t0)
	nits=1
	# provide "fitset", alpha, beta, gamma (and other prams as desried); get back an error.
	# get residual from one brutefit. alpha, beta, gamma should be fixed values (though they may be randomly generated
	# externally).
	#print "starting set... % s" % fitset
	#XYY=brutefitAlphaBetaGamma(fitset=fitset, nits=nits, p=p, b=b, dm=dm, alpha0=alpha0, beta0=beta0, gamma0=gamma0, debug=False, writefile=False)
	lyth=[]
	for x in X:
		#print 'about to calc...', alpha0, beta0, gamma0
		#print logNOmoriScaling(mstar, mc, x, alpha=alpha0, beta=beta0, gamma=gamma0)
		lyth+=[logNOmoriScaling(mstar, mc, x, alpha=alpha0, beta=beta0, gamma=gamma0)]
	
	#lyth = logNOmoriScaling(mstar, mc, scipy.array(X), alpha=alpha0, beta=beta0, gamma=gamma0)
	#print "len xyy: ", len(lyth)
	#
	ndof=0
	err=0.0
	#errs=(scipy(XYY[1])-scipy(XYY[2]))**2.0
	for i in xrange(len(lyth)):
		if Y[i]!=None and lyth[i]!=None:
			ndof+=1
			#print ndof, i, float(XYY[1][i]), float(XYY[2][i])
			err+=(float(Y[i])-float(lyth[i]))**2.0
	ndof-=3
	err/=float(ndof)
	#print "res err: ", err
	#
	#if rq!=None: rq.put(err)
	#
	return err

def getXY(fitset='parkfield'):
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
		mstar=7.2
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
	#print "data set selected (%s) %s." % (fitset, str(mev))
	#
	while(cpf.getcat(useCatnum)[0][0]<mev[0]): 
		cpf.getcat(useCatnum).pop(0)
	#
	#print "new catlen: %d" % len(cpf.getcat(useCatnum))
	fmev=mpd.date2num(mev[0])
	X = (scipy.array(map(mpd.date2num, map(operator.itemgetter(0), cpf.getcat(useCatnum)))) - fmev)*days2secs
	Y = numpy.arange(1., len(X)+1)
	lY = map(math.log10, Y)
	#
	return [X, Y, lY]
#
def myrand(R, rng):
	return rng[0] + (rng[1]-rng[0])*R.random()
	
def fullmontebruteRfit(nits=10**6, alpha0=None, beta0=None, gamma0=None, reduceErr=True):
	alphas=[-1.5, 1.5]
	betas = [-1.5, 1.5]
	gammas=[-10., 10.]
	R1 = random.Random()
	R2 = random.Random()
	R3 = random.Random()
	#
	nmproc=mcp.cpu_count()
	#
	# forget about queue()s and pool()s. let's just create N processes, initialize with alpha,beta,gamma,
	# start them -- then join() to wait for them to finish.
	#work_queue = multiprocessing.Queue()
	#work_queue = None
	result_queue = mcp.Queue()
	# P = mcp.Pool(processes=6)
	totallen=0
	workers = []	# will be an array of rfitWorker() class objects
	xysets={}
	#
	eqname, thism, thismc, thisq, thisb, thisdm='parkfield', 5.97, 1.5, 1.37, 1.0, 1.0
	XY=getXY(eqname)
	xysets[eqname]=XY
	totallen+=(len(XY[0]))
	workers +=[rfitWorker(result_queue=result_queue, name=eqname, m=thism, mc=thismc, q=thisq, b=thisb, dm=thisdm, X=scipy.array(XY[0]).copy(), Y=scipy.array(XY[2]).copy())]
	#
	eqname, thism, thismc, thisq, thisb, thisdm='hmine', 7.1, 3.0, 1.37, 1.0, 1.0
	XY=getXY(eqname)
	xysets[eqname]=XY
	totallen+=(len(XY[0]))
	workers +=[rfitWorker(result_queue=result_queue, name=eqname, m=thism, mc=thismc, q=thisq, b=thisb, dm=thisdm, X=scipy.array(XY[0]).copy(), Y=scipy.array(XY[2]).copy())]
	#	
	eqname, thism, thismc, thisq, thisb, thisdm='sansim', 6.5, 3.0, 1.37, 1.0, 1.0
	XY=getXY(eqname)
	xysets[eqname]=XY
	totallen+=(len(XY[0]))
	workers +=[rfitWorker(result_queue=result_queue, name=eqname, m=thism, mc=thismc, q=thisq, b=thisb, dm=thisdm, X=scipy.array(XY[0]).copy(), Y=scipy.array(XY[2]).copy())]
	#
	eqname, thism, thismc, thisq, thisb, thisdm='coalinga',6.7, 3.5, 1.37, 1.0, 1.0
	XY=getXY(eqname)
	xysets[eqname]=XY
	totallen+=(len(XY[0]))
	workers +=[rfitWorker(result_queue=result_queue, name=eqname, m=thism, mc=thismc, q=thisq, b=thisb, dm=thisdm, X=scipy.array(XY[0]).copy(), Y=scipy.array(XY[2]).copy())]
	#
	eqname, thism, thismc, thisq, thisb, thisdm='tohoku', 9.0, 4.5, 1.37, 1.0, 1.0
	XY=getXY(eqname)
	xysets[eqname]=XY
	totallen+=(len(XY[0]))
	workers +=[rfitWorker(result_queue=result_queue, name=eqname, m=thism, mc=thismc, q=thisq, b=thisb, dm=thisdm, X=scipy.array(XY[0]).copy(), Y=scipy.array(XY[2]).copy())]
	#	
	eqname, thism, thismc, thisq, thisb, thisdm='elmayor', 7.2, 2.5, 1.37, 1.0, 1.0
	XY=getXY(eqname)
	xysets[eqname]=XY
	totallen+=(len(XY[0]))
	workers +=[rfitWorker(result_queue=result_queue, name=eqname, m=thism, mc=thismc, q=thisq, b=thisb, dm=thisdm, X=scipy.array(XY[0]).copy(), Y=scipy.array(XY[2]).copy())]
	#
	totalNdof=float(totallen-18)
	#
	minerr=10.0**12.0
	minprams=[]	# does not need to be initialized.
	#
	i=0
	while nits==None or i<nits:
		i+=1
		totalerr=0.0
		myworkers=[]
		
		for w in workers:
			myworkers +=[rfitWorker(result_queue=result_queue, name=w.name, m=w.m, mc=w.mc, q=w.q, b=w.b, dm=w.dm, X=w.X.copy(), Y=w.Y.copy())]
			myworkers[-1].reduceErr=reduceErr
		#
		#resultset=[]
		alpha=alpha0
		beta=beta0
		gamma=gamma0
		#
		if alpha0==None: alpha=myrand(R1, alphas)
		if beta0==None: beta=myrand(R2, betas)
		if gamma0==None: gamma=myrand(R3, gammas)
		#
		#print "starting workers."
		for j in xrange(len(myworkers)):
			myworkers[j].setalphabetagamma(alpha, beta, gamma)
			myworkers[j].start()
			#myworkers[j].join()
			#print "myworker %d started..." % j
		#
		#print "workers should be started. now join."
		# i'm not sure this is fully necessary, but it's probably the better way to execute the code.
		for j in xrange(len(myworkers)):
			myworkers[j].join()	
		#print "joined and waiting?"
		#
		# results should be in results_que
		# for now, assume errors are reduced.
		#print "q-size: " , result_queue.qsize(), result_queue.empty()
		#for r in result_queue:
		while result_queue.empty()==False:
			thisval=result_queue.get()
			#print "rval: ", thisval
			#totalerr+=result_queue.get()
			totalerr+=thisval
		#
		if reduceErr==False:
			#"aftershock" weighted (vs. mainshock weighted)
			totalerr/=totalNdof
		#
		if totalerr<minerr:
			print "new error: %f, (%f, %f, %f)" % (totalerr, alpha, beta, gamma)
			#
			minerr=totalerr
			minprams=[alpha, beta, gamma]
		#
	print "finished."
	
#
def fullmontebrutefit(nits=10**6, alpha0=None, beta0=None, gamma0=None):
	# do a brute force MC fit on all earthquakes simultaneously (aka, minimize the error from all earthquakes).
	# pool() might be the best way to do this, but i think it can be done with just processes.
	#
	
	#pq = mcp.Queue(mcp.cpu_count())		# job/process queue
	#rq = mcp.Queue()
	#
	# set up random numbers...
	alphas=[-1.5, 1.5]
	betas = [-1.5, 1.5]
	gammas=[-10., 10.]
	R1 = random.Random()
	R2 = random.Random()
	R3 = random.Random()
	#
	totallen=0
	#P=mcp.Pool(processes=mcp.cpu_count())
	pramsets=[]
	XY=getXY('parkfield')
	totallen+=len(XY[0])
	pramsets+=[['parkfield', 5.97, 1.5, 1.08, 1.0, 1.0, alpha0, beta0, gamma0, XY[0][:], XY[2][:] ]]
	XY=getXY('hmine')
	totallen+=len(XY[0])
	pramsets+=[['hmine', 7.1, 3.0, 1.08, 1.0, 1.0, alpha0, beta0, gamma0, XY[0][:], XY[2][:]]]
	XY=getXY('sansim')
	totallen+=len(XY[0])
	pramsets+=[['sansim', 6.5, 3.0, 1.08, 1.0, 1.0, alpha0, beta0, gamma0, XY[0][:], XY[2]]]
	XY=getXY('coalinga')
	totallen+=len(XY[0])
	pramsets+=[['coalinga',6.7, 3.5, 1.08, 1.0, 1.0, alpha0, beta0, gamma0, XY[0][:], XY[2][:]]]
	XY=getXY('tohoku')
	totallen+=len(XY[0])	
	pramsets+=[['tohoku', 9.0, 4.5, 1.08, 1.0, 1.0, alpha0, beta0, gamma0, XY[0][:], XY[2][:]]]
	XY=getXY('elmayor')
	totallen+=len(XY[0])
	pramsets+=[['elmayor', 7.2, 2.5, 1.08, 1.0, 1.0, alpha0, beta0, gamma0, XY[0][:], XY[2][:]]]
	# brutefitres(fitset='parkfield', p=1.08, b=1.0, dm=1.0, alpha0=None, beta0=None, gamma0=None, debug=False, writefile=False)
	#
	totalndof=float(totallen-18)
	
	#getXYlY()
	
	minerr=10.0**12.
	minprams=[0., 0., 0.]
	i=0
	while nits==None or i<(nits+1):
		#print "starting nit %d" % i
		P=mcp.Pool(processes=mcp.cpu_count())
		i+=1
		#
		totalerr=0.0
		resultses=[]	# results objects...
		#
		if alpha0==None: thisalpha = alphas[0] + (alphas[1]-alphas[0])*R1.random()
		if beta0==None:  thisbeta  = betas[0] + (betas[1]-betas[0])*R2.random()
		if gamma0==None: thisgamma = gammas[0] + (gammas[1]-gammas[0])*R3.random()
		#
		for i in xrange(len(pramsets)):
		#for ps in pramsets:
			ps=pramsets[i]
			ps[6]=thisalpha
			ps[7]=thisbeta
			ps[8]=thisgamma
			resultses+=[P.apply_async(brutefitres, ps)]
		P.close()
		P.join()
		
		#
		for R in resultses:
			#print "ready? ", R.ready(), R.successful()
			if R.ready() and R.successful():
				totalerr+=float(R.get())
			else:
				print "problem.", R
		P.terminate()
		#
		totalndof=1.0
		totalerr/=totalndof
		#
		#print "err, minerr, prams: ", totalerr, minerr, thisalpha, thisbeta, thisgamma
		if totalerr<minerr:
			print "new min: %f, [%f, %f, %f]" % (totalerr, thisalpha, thisbeta, thisgamma)
			minprams=[thisalpha, thisbeta, thisgamma]
			minerr=totalerr
	#
	#
	return minerr
####
####
class rfitWorker(mcp.Process):
	def __init__(self, work_queue=None, result_queue=None, name='rfitworker', m=None,mc=None,q=None,b=None,dm=None, X=None, Y=None):
		mcp.Process.__init__(self)
		#
		self.name=name
		#
		# job management:
		self.work_queue = work_queue
		self.result_queue = result_queue
		self.kill_received = False
		#
		# permanent valued members:
		# these will all be required, but can be set after initialization.
		self.X=scipy.array(X).copy()
		self.Y=scipy.array(Y).copy()
		#print Y[0:10]
		#self.lY = map(math.log10, Y)	# we've already taken the log (for expediency)
		self.m=m
		self.mc=mc
		self.q=q
		self.dm=dm
		self.b=b
		#
		self.alpha=0.
		self.beta=0.
		self.gamma=0.
		#
		self.ndof=3.0
		self.reduceErr=True
		#
	def run (self):
		#
		while not self.kill_received:
			# get a task
			if self.work_queue!=None:
				# then we're getting prams from a queue, rather than having them assigned.
				try:
					pramset = self.work_queue.get_nowait()
				except Queue.Empty:
					break
				#
				# at this point, it looks like the strategy is to just provide alpha,beta,gamma in the pramset.
				self.setalphabetagamma(pramset[0], pramset[1], pramset[2])	
			#
			thiserr=self.brutefitres()
			if self.result_queue!=None: self.result_queue.put(thiserr)
			return thiserr
	#
	def setalphabetagamma(self, alpha=None, beta=None, gamma=None):
		self.alpha=alpha
		self.beta=beta
		self.gamma=gamma
	
	def logNOmoriScaling(self, m, mc, r, alpha=-.5, beta=-.5, gamma=9.0, q=1.08, dm=1.0, b=1.0):
		r0=10**(alpha*m + beta*mc + gamma)
		r=scipy.array(r)
		#print "val lNomori: ", (1.0 + t/t0)**(1.0-p)
		if (1.0 + r/r0)**(1.0-q)>=1.0: return None
		#N=b*(m-dm-mc) + math.log10(p-1.0) + math.log10( (1.0 + t/t0)**(1.0-p) - 1.0)
		N=b*(m-dm-mc) + math.log10( 1.0 - (1.0 + r/r0)**(1.0-q))
		return N
	
	def brutefitres(self, alpha0=None, beta0=None, gamma0=None, updateprams=False, reduceErr=None):
		#fitset='parkfield',mstar=6.0, mc=2.5, p=1.08, b=1.0, dm=1.0, alpha0=None, beta0=None, gamma0=None, X=None, Y=None):
		# "results queue", "process queue"
		# (for omori t0)
		# for now, assume we'll always use class-local values. we reserve the option to pass alpha,beta,gamma or set the class vals.
		# note, however, that when we run this as a parallel process, we'll need to populate the class alpha, beta, gamma prams.
		if alpha0==None: alpha0 = self.alpha
		if beta0==None: beta0 = self.beta
		if gamma0==None: gamma0 = self.gamma
		if reduceErr==None: reduceErr=self.reduceErr
		#
		if updateprams==True:
			self.setalphabetagamma(alpha0, beta0, gamma0)
		#
		# first, let's try to do this the fast way...
		lyth=[]
		for x in self.X:
		#	#print 'about to calc...', alpha0, beta0, gamma0
		#	#print logNOmoriScaling(mstar, mc, x, alpha=alpha0, beta=beta0, gamma=gamma0)
			lyth+=[self.logNOmoriScaling(self.m, self.mc, x, alpha=self.alpha, beta=self.beta, gamma=self.gamma, q=self.q, dm=self.dm, b=self.b)]
		#lyth = self.logNOmoriScaling(self.m, self.mc, scipy.array(self.X), alpha=self.alpha, beta=self.beta, gamma=self.gamma, q=self.q, dm=self.dm, b=self.b)
		#
		#print "len xyy: ", len(lyth)
		#
		#
		# can we get away with the fast way?
		#
		# we get NONE types that screw up the fast operation:
		'''
		lyth=scipy.array(lyth)
		print self.Y[0:10], lyth[0:10]
		errs=(self.Y-lyth)**2.0
		Ndof=len(errs)
		'''
		#
		Ndof=0
		err=0.0
		for i in xrange(len(lyth)):
			if self.Y[i]!=None and lyth[i]!=None:
				Ndof+=1
				#print ndof, i, float(XYY[1][i]), float(XYY[2][i])
				err+=(float(self.Y[i])-float(lyth[i]))**2.0
		
		Ndof-=self.ndof
		if Ndof<=0: Ndof=1
		
		#err=float(sum(errs))
		if reduceErr==True:
			#Ndof=float(len(err))-self.ndof
			err/=float(Ndof)
		#print "res err: ", err
		#
		#if rq!=None: rq.put(err)
		#
		return err
		
		
		
		
		
