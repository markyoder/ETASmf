from geographiclib.geodesic import Geodesic as ggp
import math


#def sphericalDist(phis, lambs, phif, lambf):
def sphericalDist(lambs, phis, lambf, phif):
	# displacement from inloc...
	# inloc is a vector [lon, lat]
	# return a vector [dLon, dLat] or [r, theta]
	# return distances in km.
	#
	# also, we need to get the proper spherical angular displacement (from the parallel)
	#
	Rearth = 6378.1	# km
	deg2rad=2.0*math.pi/360.
	#
	phif*=deg2rad
	phis*=deg2rad
	lambf*=deg2rad
	lambs*=deg2rad
	#
	#phif  = inloc[0]*deg2rad
	#lambf = inloc[1]*deg2rad
	#phis  = self.loc[0]*deg2rad
	#lambs = self.loc[1]*deg2rad
	#
	#print 'phif: ', phif
	#print 'lambf: ', lambf
	#
	dphi = (phif - phis)
	dlambda = (lambf - lambs)
	#this one is supposed to be bulletproof:
	sighat3 = math.atan( math.sqrt((math.cos(phif)*math.sin(dlambda))**2.0 + (math.cos(phis)*math.sin(phif) - math.sin(phis)*math.cos(phif)*math.cos(dlambda))**2.0 ) / (math.sin(phis)*math.sin(phif) + math.cos(phis)*math.cos(phif)*math.cos(dlambda))  )
	R3 = Rearth * sighat3
	#
	return R3
	
def distTests():
	# [[ x1,y1,x2,y2], ...]
	#
	x0=-121.0
	y0=39.0		# cali-ish...
	#
	xys=[[x0,y0,x0+1.,y0+0.],[x0,y0,x0-1.,y0+0.],[x0,y0,x0+0.,y0+1.],[x0,y0,x0+0.,y0-1.], [x0,y0,x0+1.,y0+1.], [x0,y0,x0+1.,y0-1.], [x0,y0,x0-1.,y0-1.], [x0,y0,x0-1.,y0+1.] ]
	#
	for xy in xys:
		sphdist=sphericalDist(xy[0], xy[1], xy[2], xy[3])
		g1=ggp.WGS84.Inverse(xy[1], xy[0], xy[3], xy[2])		# lat, lon, lat, lon (i think
		Rgeolib=g1['s12']/1000.0	# and converting from m to km.
		#Rgeolib=g1['s12']
		#
		print("spherical: %f, lib: %f" % (sphdist, Rgeolib))
	
	
	#g1=ggp.WGS84.Inverse(self.loc[1], self.loc[0], inloc[1], inloc[0])
	#Rgeolib=g1['s12']/1000.0	# and converting from m to km.
