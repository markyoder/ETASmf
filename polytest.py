import yodapoly as ypp
import shapely.geometry as sgp

import matplotlib.pyplot as plt
import operator


def polytest1():
	p=[]
	p += [sgp.Polygon([[0,0], [10,0], [10,10], [0,10], [0,0]])]
	p += [sgp.Polygon([[1,1], [7,1], [7,7], [1,7], [1,1] ])]
	p += [sgp.Polygon([[2,2], [5,2], [5,5], [2,5], [2,2] ])]
	p += [sgp.Polygon([[5.5,5.5], [6.5, 5.5], [6.5,6.5], [5.5, 6.5], [5.5,5.5]])]
	p += [sgp.Polygon([[2.5,2.5], [3.5, 2.5], [3.4,3.4], [2.5, 3.5], [2.5,2.5]])]
	p += [sgp.Polygon([[3.6,3.6], [4.6, 3.6], [4.6,4.6], [3.6, 4.6], [3.6,3.6]])]
	p += [sgp.Polygon([[13.6,13.6], [14.6, 13.6], [14.6,14.6], [13.6, 14.6], [13.6,13.6]])]
	#
	#polylist=[p1, p2, p3, p4]
	#
	print "is p2 in p1? " , p[0].contains(p[1])
	print "is p1 in p2? " , p[1].contains(p[2])
	
	print "is p3 in p1,p2?", p[0].contains(p[2]), p[0].contains(p[1]), p[1].contains(p[2])
	p2 = []
	for pp in p:
		p2+=[ [pp.exterior.xy[0].tolist(), pp.exterior.xy[1].tolist()] ]
	#
	c=plotpolys(p)
	
	return [p, p2]

def innerouterpolys(self, polylist):
	#
	polylistplus=[]		# indexed entries: [polyindex, [list of inners], [verts] ]
	outerpolys=[]			# these will be lists. each entry is like: [[outer],[true-inner],[true-inner],..]

def testpolys1():
	polys=[]
	#
	X=[-115.0, -114.75, -114.5, -114.27428996711201, -114.35813215408139, -114.5, -114.62156320467432, -114.75, -115.0, -115.02000478717538, -115.25, -115.5, -115.65040828540928, -115.5, -115.31821871871739, -115.5, -115.7085735116233, -115.75, -116.0, -116.25, -116.30824837902244, -116.29358780688531, -116.25, -116.0584245131186, -116.0, -115.75, -115.69392386639223, -115.5, -115.25, -115.15489394301572, -115.0]
	Y=[31.685018937456629, 31.62597355151933, 31.626091175355523, 31.75, 32.0, 32.131746497417744, 32.25, 32.343751947686862, 32.489055129171554, 32.5, 32.651504509265358, 32.718733806587899, 32.75, 32.76654905381168, 33.0, 33.15839050065258, 33.0, 32.797787829158125, 32.84353128477926, 32.783891953729167, 32.75, 32.5, 32.457735514850903, 32.25, 32.205954277573312, 32.031597849209717, 32.0, 31.881232067737159, 31.77493355694515, 31.75, 31.685018937456629]
	polys+=[[X,Y]]
	X,Y=[-115.5, -115.69610748055838, -115.5, -115.32908343054898, -115.5], [32.780501960509511, 33.0, 33.148923810054349, 33.0, 32.780501960509511]
	polys+=[[X,Y]]
	X,Y=[-115.0, -115.0643533209925, -115.25, -115.5, -115.63923114903437, -115.75, -115.99861212704302, -116.0, -116.23031709342109, -116.18084840781762, -116.0, -115.7866100373639, -115.75, -115.5, -115.25, -115.06720028983338, -115.0, -114.75, -114.6773251468034, -114.5, -114.41654315664718, -114.37180455215244, -114.5, -114.75, -115.0],   [31.723002513236668, 31.75, 31.798670212912693, 31.914728413861418, 32.0, 32.062416169317707, 32.25, 32.252025331111355, 32.5, 32.75, 32.801983458313551, 32.75, 32.743642036490606, 32.69237252427898, 32.620415464489533, 32.5, 32.463233875701256, 32.303048731228152, 32.25, 32.077502740145796, 32.0, 31.75, 31.679624097501033, 31.66455172039289, 31.723002513236668]
	polys+=[[X,Y]]
	
	plt.figure(0)
	plt.clf()
	plt.ion()
	for poly in polys:
		plt.plot(poly[0], poly[1], '--')
	
	iop=innerouterpolys1(polys)
	print len(iop)
	
	
	return [polys, iop]
	
def innerouterpolys1(polylist):
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
	print polylist[1]
	for i in xrange(len(polylist)):
		polylistplus += [[i, [], polylist[i]]]
		#
		# shorthand:
		# which poly each polygon it is inside, and how many polys it is inside (aka, len(that list)).
		#x0,y0=poly[i][0][0], poly[i][1][0]	# we only need to test one point since we're starting with contours (don't cross).
		#print polylistplus[-1][2]
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
				'''
				if y0>=min(y1, y2) and y0<=max(y1, y2):
					# this segment is in the y-range and its xmax is to the right of the point.
					fx = (y2-y1)/(x2-x1)*(x0-x1)	# maybe do this in transpose space...
					fy = (x2-x1)/(y2-y1)*(y0-y1)
					#if x0<=min(x1,x2) or (x0>=min(x1,x2) and x0<=max(x1,x2) and fx>(y0-y1)):
					if fy>x0:
						# the point is either left of the left-most point of the segment, or it's under the line...
						ncrossings+=1
					#
				'''
				if x0>=min(x1, x2) and x0<max(x1, x2):
					fx = (y2-y1)/(x2-x1)*(x0-x1)
					if fx>=(y0-y1):
						ncrossings += 1
			#
			#print i,j,j,ncrossings, ncrossings%2
			if ncrossings%2==1:
				# i'th poly is inside j'th poly...
				print "poly %d is inside poly %d" % (i, j)
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

def plotpolys(polylist, fnum=0):
	# assume polylist is like:
	# [ [poly1], [poly2], ...] where poly_i like [ [x0,y0], [x1,y1], ...]
	plt.figure(fnum)
	plt.ion()
	plt.clf()
	#
	for poly in polylist:
		#X,Y = map(operator.itemgetter(0), poly), map(operator.itemgetter(1), poly)
		X=poly.exterior.xy[0]
		Y=poly.exterior.xy[1]
		print X
		plt.plot(X,Y, 'o-')

def inoutplots(outers, fnum=1):
	# outers is a list of list (of lists...)
	plt.ion()
	plt.figure(fnum)
	plt.clf()
	# quick set of colors...
	clrs=['b', 'g', 'r', 'c', 'm', 'y', 'k']	#... i think
	#
	for i in xrange(len(outers)):
		outer=outers[i]
		clr=clrs[i%len(clrs)]
		# first poly is the outer boundary; subsequent polys are inners.
		X=outer[0][0]
		Y=outer[0][1]
		#
		plt.plot(X,Y, '-', color=clr)
		if len(outer)==1: continue
		for j in xrange(1,len(outer)):
			inner=outer[j]
			X=inner[0]
			Y=inner[1]
			#
			plt.plot(X,Y, '--', color=clr)	
		
		




