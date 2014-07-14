import matplotlib.pyplot as plt
import scipy
import math

defaultmc=2.5

def lrlen(m, dlamb=1.76):
	return m/2.0 - dlamb

def rlen(m, dlamb=1.76):
	return 10.0**(lrlen(m, dlamb))

def ldtr(m, dtau=2.4):
	return m/2.0 - dtau
	
def dtr(m, dtau=2.4):
	return 10.0**ldtr(m, dtau)

def lNom(m, mc, b=1.0, dm=1.0):
	return float(b*(m-dm-mc))

def Nom(m, mc=None, b=1.0, dm=1.0):
	if mc==None: mc=defaultmc
	#
	return 10.0**lNom(m,mc,b,dm)

def mytau(t0=None, p=1.08, N=None, m=None, mc=None):
	if t0==None: t0=dtr(m)
	if N==None: N=Nom(m, mc)
	#
	tau = (t0**(1.0-p))/(N*(p-1.0))
	#tau = (t0**(1.0-p))/(N*(p-1.0))
	return tau

def dndt(t, tau, t0, p):
	return (1.0/tau) * (t0+t)**(-p)

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

def tautestses(m, mc, ps=[1.01, 2.51, .1], tends=[1., 2., 5., 10., 25.]):
	yrsec=float(3600*24*365.24)
	#
	plt.figure(3)
	plt.ion()
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	#
	for lt in tends:
		tend = lt*yrsec
		PN=tautests(m=m, mc=mc, ps=ps, ltend=(math.log10(lt*yrsec)))
		#print "PN: ", PN, (math.log10(lt*yrsec))
		#
		plt.figure(3)
		plt.plot(PN[0], PN[1], '-o', label='$\\Delta t = %.3f$ years' % lt)
	plt.legend(loc='best')
	
	
def tautests(m, mc, ps=[1.01, 2.51, .1], ltend=35):
	p=ps[0]
	plt.figure(0)
	plt.clf()
	psout, Ns=[],[]
	while p<=ps[1]:
		tt=tautest(m, mc, p, doclf=False, ltend=ltend)
		p+=ps[2]
		#
		psout+=[tt[0]]
		Ns+=[tt[1]]
		
	tautest(m,mc,p=1.08, doclf=False, linestyle='--', lw=3)
	plt.legend(loc='best', title='$m=%f, m_c=%f$ \n $lrl=%f$' % (m, mc, lrlen(m)))
	#
	return [psout, Ns]
def tautest(m, mc, p=1.08, doclf=True, linestyle='-', lw=1, ltend=35):
	# the rlen and rupt.dt tests should be functionally equivalent.
	ldt=.001	# even log-spaced.
	dldt=.005
	#
	myNom=Nom(m, mc)
	#
	# t0,tau calculations:
	#t0=dtr(m)
	#tau = mytau(t0, p, myNom)
	#tau = (t0**(1.0-p))/(myNom*(p-1.0))
	#
	t0tau = t0taussim(m, mc, dmtau=0.0, dmstar=1.0, b=1.0, p=p)
	t0, tau = t0tau[0], t0tau[1]
	#
	thisldt=ldtr(m)-5.0
	dN0 = dndt(0., tau, t0, p)
	thisdN=dN0
	dNs=[dN0]
	dNs2=[dN0]
	ts=[10.0**ldt]
	#while dNs[-1]>(dN0/10000000000000.0):
	#print "ltend: %d" % ltend
	while ldt<ltend:
		ldt+=dldt
		t=10**ldt
		#
		if t>=t0:
			thisdn=dndt((t-t0), tau, t0, p)
			#thisdn=(myNom*(p-1.0)*(t0**(p-1.0))*((t)**-p))/p
		else:
			thisdn=dndt(0.0, tau, t0, p)
		
		dNs+=[thisdn]
		ts+=[t]
		dNs2+=[dndt((t), tau, t0, p)]
	#
	total=0.0
	prats=[]
	pset=[]
	for i in xrange(1, len(ts)):
		total+=(ts[i]-ts[i-1])*dNs2[i]
		#total+=( ((t0+ts[i])**(1.0-p)) - ((t0+ts[i-1])**(1.0-p)))/(tau*(1.0-p))
		
	print "total %f: %f/%f: %f [%f, %f]" % (p, total, Nom(m,mc), float(total)/Nom(m,mc), t0, tau)
	#
	#
	plt.figure(1)
	#plt.clf()
	#plt.ion()
	ax1=plt.gca()
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	plt.plot([p], [float(total)/Nom(m,mc)], 'o')
	plt.title('Aftershock Recovery Efficiency')
	plt.xlabel('Temporal Scaling Exponent $p$')
	plt.ylabel('Recovery rate $N/N_{Omori}$')
	
	plt.figure(0)
	if doclf: plt.clf()
	plt.ion()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(ts, dNs, linestyle, label='$(t=t-t0): p=%.2f' % p, lw=lw)
	plt.plot(ts, dNs2, '--', label='$(t=t): p2=%.2f' % p, lw=lw)
	plt.plot([t0, t0], [min(dNs), max(dNs)], '--', lw=2, alpha=.6)
	#plt.plot([ts[0], t0], [1.0/(tau*t0**p), 1.0/(tau*t0**p)], 'd--', lw=3, alpha=.5)
	plt.xlabel('Time $t$', size=18)
	plt.ylabel('Rate $dN/dt$', size=18)
	#
	return [p, float(total)/Nom(m,mc)]

def fracrecoveries(fracs=[.19, .99, .1], m=6.0, mc=1.5, ps=[1.01, 2.5, .01], b=1.0, dmstar=1.0, dmtau=2.0, dtau=2.28):
	frac=fracs[0]
	while frac<=fracs[1]:
		doclf=False
		if frac==fracs[0]: doclf=True
		fractRecovery1(frac=frac, m=m, mc=mc, ps=ps, b=b, dmstar=dmstar, dmtau=dmtau, dtau=dtau, doclf=doclf)
		frac+=fracs[2]
	plt.legend(loc='best')
		
def fractRecovery1(frac=.5, m=6.0, mc=1.5, ps=[1.01, 2.5, .01], b=1.0, dmstar=1.0, dmtau=2.0, dtau=2.28, doclf=True):
	#
	ts=[]
	plotps=[]
	p=ps[0]
	while p<=ps[1]:
		lt0 = math.log10(p-1.0) + m*7.0/6.0 - (2.0*mc/3.0) + (dmtau/3.0) - dmstar - dtau - ((2./3.)*math.log10(1.5))	
		t0=10.0**lt0
		thisexp=1.0/(1.0-p)
		tfrac=t0*(1.0-frac)**thisexp
		#
		ts+=[tfrac]
		plotps+=[p]
		#
		p+=ps[2]
	#
	mint=ts[0]
	minindex=0
	for j in xrange(len(ts)):
		if ts[j]<mint:
			mint=ts[j]
			minindex=j
	#
	plt.ion()
	plt.figure(0)
	if doclf: plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	#
	plt.plot(plotps, ts, label='$f = %.2f, p_{min}=%.3f$' % (frac, plotps[minindex]))
	plt.plot([plotps[minindex]], [ts[minindex]], '*')
	
 	
	
	
	
	
	
	
