#!/usr/bin/env python

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import os, sys, inspect
import warnings
import matplotlib.ticker as ticker

#show_rN = True
show_rN = False

#makePDF = True
makePDF = False
for arg in sys.argv:
   if arg.lower()=='pdf':
      makePDF = True

if makePDF:
   matplotlib.use('PDF')
else:
   matplotlib.use('qt5agg')

import matplotlib.pyplot as plt

print("This is "+ inspect.getfile(inspect.currentframe()))

#########
##INPUT##
#########

quantityToPlot = "BHat"

filename = 'sfincsOutput.h5'

FigSize = (12,10)

font = {'size':25}
matplotlib.rc('font', **font)
matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
matplotlib.rc('axes',linewidth=1.5)

#matplotlib.rcParams['mathtext.default'] = 'it'
#matplotlib.rcParams['text.usetex'] = True

zFactor = 1 ##T
###W7-X##
#xAxisTicks = [r'$0$', r'$\pi/10$', r'$2\pi/10$', r'$3\pi/10$', r'$4\pi/10$']
###LHD
##xAxisTicks = [r'$0$', r'$\pi/20$', r'$2\pi/20$', r'$3\pi/20$', r'$4\pi/20$']
#
#yAxisTicks = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']


fig = plt.figure(figsize=FigSize)
fig.patch.set_facecolor('white')
numRows = 1
numCols = 1
#iteration = 0
numContours = 100
#ContourLevels = [2.7, 2.8, 2.9, 3.0, 3.1, 3.2]
numLevels = 5

ColorMap = 'rainbow'

#############
##END INPUT##
#############

#def fmt_cbar(x, pos):
#   if x == 0.0:
#      #return r'${}$'.format(x)
#      return r'{}'.format(x)
#   a, b = '{:.1e}'.format(x).split('e')
#   b = int(b)
#   #return r'${} \cdot 10^{{{}}}$'.format(a, b)
#   return r'${} \cdot 10^{{{}}}$'.format(a, b)

def fmt_xy_axis(x, pos):
   #return r'${}$'.format(x)
   #return r'${}$'.format('{:1.2f}'.format(x))
   return r'{}'.format('{:1.2f}'.format(x))

#for i in range(6):
print ("Processing file ",filename)
f = h5py.File(filename,'r')
theta = f["theta"][()]
zeta = f["zeta"][()]
BHat = f[quantityToPlot][()]
rN = f["rN"][()]
f.close()

print ("theta max: " + str(numpy.amax(theta)))
print ("zeta max: " + str(numpy.amax(zeta)))

zMinData = zFactor*numpy.amin(BHat[:,:])
zMaxData = zFactor*numpy.amax(BHat[:,:])
print ("zMin = " + str(zMinData))
print ("zMax = " + str(zMaxData))


delta = (numpy.amax(BHat) - numpy.amin(BHat)) / numLevels
ContourLevels = numpy.arange(numpy.amin(BHat), numpy.amax(BHat) + delta/2.0, delta)
ContourLevels = zFactor*ContourLevels
    
ax = plt.subplot(numRows,numCols,1)
    #plt.contourf(zeta,theta,1000*numpy.fliplr(BHat[:,:,iteration].transpose()),numContours)
BPlot = plt.contourf(zeta,theta,zFactor*BHat.transpose(),numContours, cmap=plt.get_cmap(ColorMap))
BPlot2 = plt.contour(BPlot,levels=ContourLevels, colors='k', hold='on')
#BPlot2 = plt.contour(BPlot,levels=BPlot.levels[::2], colors='k', hold='on')
#plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
#plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')
plt.xlabel(r'zeta' + " " + r'[rad]')
plt.ylabel(r'theta'+ " " + r'[rad]')
#plt.zlabel(r'$B$'+ ' [T]')
plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)])
plt.yticks([0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)])
#plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
#plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

#plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
#plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)
plt.gca().axes.xaxis.set_label_coords(0.5,-0.05)
plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

if show_rN:
    plt.title('rN = '+str(rN))

#cbar = plt.colorbar(BPlot, label=r'$B$'+ ' [T]', ticks=ContourLevels)
#cbar = plt.colorbar(BPlot, label=r'$\Phi_1$'+ ' [V]', ticks=BPlot.levels[::2])
#cbar.add_lines(BPlot2)
#cbar = plt.colorbar(BPlot, format=ticker.FuncFormatter(fmt_xy_axis), ticks=ContourLevels)
cbar = plt.colorbar(BPlot, ticks=ContourLevels)
#cbar.ax.set_ylabel(r'$B$'+ " " + r'$\mathrm{[T]}$', rotation=0, labelpad=10)
cbar.ax.set_ylabel(r'B'+ " " + r'[T]', rotation=0, labelpad=10)

#with warnings.catch_warnings():
#    warnings.simplefilter("always")
#plt.clabel(BPlot2, fmt='%2.1f', colors='k', fontsize=14)
#plt.clabel(BPlot2, fmt=ticker.FuncFormatter(fmt_xy_axis), colors='k', fontsize=18, inline=False)
plt.clabel(BPlot2, colors='k', fontsize=18, inline=False)

#plt.subplots_adjust(wspace=0.27)

print (BHat.shape)

if makePDF:
    print ("Saving PDF")

    if len(sys.argv)>2 : #Use the substituted name as file name
       print ("Writing plot to " + os.getcwd() + "/" + sys.argv[2] + ".pdf.")
       plt.savefig(sys.argv[2] + ".pdf", orientation = 'landscape', papertype='letter')
    else :
       head, tail = os.path.split(inspect.getfile(inspect.currentframe()))
       print ("Writing plot to " + os.getcwd() + "/" + tail + ".pdf.")
       plt.savefig(tail+'.pdf', orientation = 'landscape', papertype='letter')
else:
    plt.show()
