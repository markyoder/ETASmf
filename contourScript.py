#!/usr/bin/python
#
# contourScript.py:
# Primary Author: Mark R. Yoder, Ph.D.
#                 mryoder@ucdavis.edu
#                 mark.yoder@gmail.com
#
# this is an example script for transforming pyplot.contourf() contours to KML contours. see also my contours2kml repository:
# https://github.com/markyoder/contours2kml
#
#
import numpy

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot  as plt
import matplotlib.mpl     as mpl

import sys
import os

#-----------------------------------------------------------------------------#

resolution = 5

LevelsNumber = 5 * resolution
warnings = ['No','Low','Guarded','Elevated','High','Severe']

#-----------------------------------------------------------------------------#

# Create your data here
#
# I assume there are three arrays: X_i, Y_i, and Z_i
Z_i=numpy.loadtxt('tmp/Zarray.dat')
#X_i = numpy.array([-127.5, -127.4, -127.3, -127.2, -127.1, -127. , -126.9, -126.8, -126.7, -126.6, -126.5, -126.4, -126.3, -126.2, -126.1, -126. , -125.9, -125.8, -125.7, -125.6, -125.5, -125.4, -125.3, -125.2, -125.1, -125. , -124.9, -124.8, -124.7, -124.6, -124.5, -124.4, -124.3, -124.2, -124.1, -124. , -123.9, -123.8, -123.7, -123.6, -123.5, -123.4, -123.3, -123.2, -123.1, -123. , -122.9, -122.8, -122.7, -122.6, -122.5, -122.4, -122.3, -122.2, -122.1, -122. ,-121.9, -121.8, -121.7, -121.6, -121.5, -121.4, -121.3, -121.2, -121.1, -121. , -120.9, -120.8, -120.7, -120.6, -120.5, -120.4, -120.3, -120.2, -120.1, -120. , -119.9, -119.8, -119.7, -119.6, -119.5, -119.4, -119.3])
#Y_i = numpy.array([ 37.1,  37.2,  37.3,  37.4,  37.5,  37.6,  37.7,  37.8,  37.9, 38. ,  38.1,  38.2,  38.3,  38.4,  38.5,  38.6,  38.7,  38.8, 38.9,  39. ,  39.1,  39.2,  39.3,  39.4,  39.5,  39.6,  39.7, 39.8,  39.9,  40. ,  40.1,  40.2,  40.3,  40.4,  40.5,  40.6, 40.7,  40.8,  40.9,  41. ,  41.1,  41.2,  41.3,  41.4,  41.5, 41.6,  41.7,  41.8,  41.9,  42. ])

Y_i=range(len(Z_i))
X_i=range(len(Z_i[0]))

#-----------------------------------------------------------------------------#

# Create the KML file
cs = plt.contourf(X_i, Y_i, Z_i, LevelsNumber,
                  cm=plt.spectral()).collections

file = open('tmp/doc.kml','w')
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

