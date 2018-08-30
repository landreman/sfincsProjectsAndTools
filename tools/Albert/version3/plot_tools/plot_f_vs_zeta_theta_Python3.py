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
   #import matplotlib.pyplot as plt
   matplotlib.use('qt5agg')

import matplotlib.pyplot as plt

print ("This is "+ inspect.getfile(inspect.currentframe()))

#########
##INPUT##
#########

#quantityToPlot = "full_f"
quantityToPlot = "delta_f"

xv = 2

xi = 2

species = 1

xLabel = r'$\zeta$' + " " + r'$\mathrm{[rad]}$'
yLabel = r'$\theta$'+ " " + r'$\mathrm{[rad]}$'
#zLabel = r'$f_{C}$'
#zLabel = r'$f_{1C}$'
#zLabel = r'$f_{i}$'
zLabel = r'$f_{1i}$'
zUnits = r'$[10^{-3} \bar{n} \bar{v}^{-3}]$'

filename = 'sfincsOutput.h5'

FigSize = (12,10)

font = {'size':25}
matplotlib.rc('font', **font)
matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
matplotlib.rc('axes',linewidth=1.5)

matplotlib.rcParams['mathtext.default'] = 'it'
matplotlib.rcParams['text.usetex'] = True

zFactor = 1000 ##kV -> V
#zFactor = 1 ##
##W7-X##
#xAxisTicks = [r'$0$', r'$\pi/10$', r'$2\pi/10$', r'$3\pi/10$', r'$4\pi/10$']
##LHD
xAxisTicks = [r'$0$', r'$\pi/20$', r'$2\pi/20$', r'$3\pi/20$', r'$4\pi/20$']

yAxisTicks = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

fig = plt.figure(figsize=FigSize)
fig.patch.set_facecolor('white')
numRows = 1
numCols = 1
#iteration = 0
numContours = 100
#ContourLevels = [-3.0, -1.5, 0.0, 1.5, 3.0, 4.5, 6.0]
numLevels = 5


#############
##END INPUT##
#############



def fmt_cbar(x, pos):
   if x == 0.0:
      return r'${}$'.format(x)
   a, b = '{:.1e}'.format(x).split('e')
   b = int(b)
   return r'${} \cdot 10^{{{}}}$'.format(a, b)

def fmt_xy_axis(x, pos):
   return r'${}$'.format(x)

#for i in range(6):
print ("Processing file ",filename)
f = h5py.File(filename,'r')
#theta = f["theta"][()]
#zeta = f["zeta"][()]
export_f_x = f["export_f_x"][()]
export_f_xi = f["export_f_xi"][()]
export_f_zeta = f["export_f_zeta"][()]
export_f_theta = f["export_f_theta"][()]
Zs = f["Zs"][()]

OutputQuantity = f[quantityToPlot][()]
iteration = f["NIterations"][()] - 1 #Results from last iteration
rN = f["rN"][()]
f.close()

print ("x max: " + str(numpy.amax(export_f_x)))
print ("x min: " + str(numpy.amin(export_f_x)))
print ("xs: " + str(export_f_x))
print ("x: " + str(export_f_x[xv]))

print ("xi max: " + str(numpy.amax(export_f_xi)))
print ("xi min: " + str(numpy.amin(export_f_xi)))
print ("xis: " + str(export_f_xi))
print ("xi: " + str(export_f_xi[xi]))

print ("theta max: " + str(numpy.amax(export_f_theta)))
print ("theta min: " + str(numpy.amin(export_f_theta)))
print ("thetas: " + str(export_f_theta))

print ("zeta max: " + str(numpy.amax(export_f_zeta)))
print ("zeta min: " + str(numpy.amin(export_f_zeta)))
print ("zetas: " + str(export_f_zeta))

zMinData = zFactor*numpy.amin(OutputQuantity[xv, xi, :, :, species,iteration])
zMaxData = zFactor*numpy.amax(OutputQuantity[xv, xi, :, :, species,iteration])
print ("zMin = " + str(zMinData))
print ("zMax = " + str(zMaxData))

print ("Zs: " + str(Zs))
print ("Z: " + str(Zs[species]))

zLabel = zLabel + r'$(x = ' + str(export_f_x[xv]) + r', \xi = ' + str(export_f_xi[xi]) + r', \zeta, \theta)$'

zLabel = zLabel + " " + zUnits

delta = (numpy.amax(OutputQuantity[xv, xi, :, :, species,iteration]) - numpy.amin(OutputQuantity[xv, xi, :, :, species,iteration])) / numLevels
ContourLevels = numpy.arange(numpy.amin(OutputQuantity[xv, xi, :, :, species,iteration]), numpy.amax(OutputQuantity[xv, xi, :, :, species,iteration]) + delta/2.0, delta)
ContourLevels = zFactor*ContourLevels
    
ax = plt.subplot(numRows,numCols,1)
    #plt.contourf(zeta,theta,1000*numpy.fliplr(OutputQuantity[xs, xis, zetas, thetas, species,iteration].transpose()),numContours)
QuantityPlot = plt.contourf(export_f_zeta,export_f_theta,zFactor*OutputQuantity[xv, xi, :, :, species,iteration].transpose(),numContours, cmap=plt.get_cmap('gist_rainbow'))
#QuantityPlot2 = plt.contour(QuantityPlot,levels=ContourLevels, colors='k', hold='on')
QuantityPlot2 = plt.contour(QuantityPlot,levels=ContourLevels, colors='k')
#QuantityPlot2 = plt.contour(QuantityPlot,levels=QuantityPlot.levels[::2], colors='k', hold='on')
plt.xlabel(xLabel)
plt.ylabel(yLabel)
#plt.zlabel(r'$\Phi_1$'+ ' [V]')

plt.xticks([0,max(export_f_zeta)/4,max(export_f_zeta)/2,3*max(export_f_zeta)/4,max(export_f_zeta)])
plt.yticks([0.0,max(export_f_theta)/4,max(export_f_theta)/2,3*max(export_f_theta)/4,max(export_f_theta)])
plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

#plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
#plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)
plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

#ax.xaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
#ax.xaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
#ax.yaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis)) 
#ax.yaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))

#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

if show_rN:
    plt.title('rN = '+str(rN))

#cbar = plt.colorbar(QuantityPlot, label=r'$\Phi_1$'+ ' [V]', ticks=ContourLevels)
#cbar = plt.colorbar(QuantityPlot, label=r'$\Phi_1$'+ ' [V]', ticks=QuantityPlot.levels[::2])
#cbar.add_lines(QuantityPlot2)
cbar = plt.colorbar(QuantityPlot, format=ticker.FuncFormatter(fmt_cbar), ticks=ContourLevels)
cbar.ax.set_ylabel(zLabel, rotation=90, labelpad=10, fontsize=30)

#with warnings.catch_warnings():
#    warnings.simplefilter("always")
#plt.clabel(QuantityPlot2, fmt='%2.1f', colors='k', fontsize=14)
plt.clabel(QuantityPlot2, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=18, inline=False)

#plt.subplots_adjust(wspace=0.27)

print (OutputQuantity.shape)

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
