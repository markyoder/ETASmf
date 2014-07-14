#!/usr/bin/python

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

Z_i = numpy.loadtxt('tmp/Zarray.dat')

Y_i =   37.0 + .1*numpy.arange(Z_i.shape[0])
X_i = -127.5 + .1*numpy.arange(Z_i.shape[1])

print Z_i.shape, len(Y_i), len(X_i)

#-----------------------------------------------------------------------------#

# Create the KML file
cs = plt.contourf(X_i, Y_i, Z_i, LevelsNumber,
                  cm=plt.spectral()).collections

file = open('doc.kml','w')
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

for i in range(len(cs)):
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

cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%g%%")
cb1.set_label('Probability')

plt.savefig('scale.png')

#-----------------------------------------------------------------------------#

os.system('zip -q tmp/forecast.kmz doc.kml scale.png')
os.system('rm -rf doc.kml scale.png')

