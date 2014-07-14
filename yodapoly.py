import shapely.geometry as sgp

class yodapoly(sgp.Polygon):
	
	innermembers=[]	# "self" is inside these polys. this list can be used to store anyting -- polys, or indeces of polys.
							# note that n=len(innermembers) implies if this is an "inner" our "outer" polygon.
							
	def __init__(self, verts):
		if len(verts)==2:
			# probably an [[X], [Y]] array:
			newverts=[]
			for i in xrange(len(verts[0])):
				newverts+=[[verts[0][i], verts[1][i]]]
			
