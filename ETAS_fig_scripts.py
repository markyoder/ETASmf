import BASScast as bcp
import matplotlib.pyplot as plt
import math
import scipy
import numpy

def omori_tahir_dist(alpha=.5, beta=.5):
	#
	# plot an omori_x + poisson_tahir type figure
	#
	X=numpy.arange(0., 30., .1)
	norm_fact = alpha+beta
	alpha /= norm_fact
	beta  /= norm_fact
	#
	Y_omori_x   = [omori_x(x) for x in X]		#omori_x(X)
	Y_poisson_x = [poisson_x(x) for x in X]	#	poisson_x(X)
	#
	YY = alpha*numpy.array(Y_omori_x) + beta*numpy.array(Y_poisson_x)
	#
	myplot=plt.loglog
	#
	plt.ion()
	plt.figure(0)
	plt.clf()
	myplot(X, Y_omori_x, label='Aftershock Seismicity', lw=3)
	myplot(X, Y_poisson_x, label='Largest Aftershock', lw=3)
	myplot(X, YY, label = 'combined, $\\alpha = %.2f$, $\\beta = %.2f$' % (alpha, beta), lw=3)
	#
	plt.legend(loc=0, numpoints=1)
	plt.xlabel('\'Normalized\' distance, $r/r_0$', size=16)
	plt.ylabel('Seismic Hazard', size=16)


def omori_x(r, r0=2., chi=1., q=1.5):
	#
	return 1.0/(chi*(r0+r)**q)

def poisson_x(r, r0=1.0, k=1):
	#
	return (math.exp(-r/r0)*(r/r0)**k)/math.factorial(k)

def nOmori(m, dm=1.0, mc=1.0, base=10):
	return base**(m-dm-mc)

def make_little_bass(fout='minibass.txt', base=3, m=6, mc=3):
	# do we have time to make a mini-bass and graph-viz it?
	# just the branching...
	#m=6.
	#mc=3.
	quakes = [[None, 0]]		# like [[src, target]]
	#
	
	for quake in quakes:
		N=nOmori(m=m-float(quake[1]), base=base)
		parent = quake[1]
		child = parent + 1
		for i in xrange(int(N)):
			quakes += [[parent, child]]
		if len(quakes)>500: 
			print "breaking"
			break
	#
	f=open(fout, 'w')
	
	for rw in quakes:
		f.write('%s -- %s\n' % (str(rw[0]), str(rw[1])))
	f.close()
	#
	return quakes
	
		
