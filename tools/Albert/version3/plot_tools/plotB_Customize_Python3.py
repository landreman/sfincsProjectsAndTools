#!/usr/bin/env python

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import os, sys, inspect
import warnings
import matplotlib.ticker as ticker
import subprocess

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

sfincsHome = os.environ.get('SFINCS_HOME')
sfincsProjectsAndToolsHome = os.environ.get('SFINCS_PROJECTS_AND_TOOLS_HOME')

exec(open(sfincsProjectsAndToolsHome + "/tools/Albert/version3/plot_tools"  + "/RadialScanPlotOptions.py").read())

#########
##INPUT##
#########

quantityToPlot = "BHat"

filename = 'sfincsOutput.h5'

#FigSize = (12,10)

#font = {'size':25}
#matplotlib.rc('font', **font)
#matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
#matplotlib.rc('axes',linewidth=1.5)

#matplotlib.rcParams['mathtext.default'] = 'it'
#matplotlib.rcParams['text.usetex'] = True

zFactor = 1 ##T
##W7-X##
xAxisTicks = [r'$0$', r'$\pi/10$', r'$2\pi/10$', r'$3\pi/10$', r'$4\pi/10$']
##LHD
#xAxisTicks = [r'$0$', r'$\pi/20$', r'$2\pi/20$', r'$3\pi/20$', r'$4\pi/20$']

yAxisTicks = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

fig = plt.figure(figsize=FigSize)
fig.patch.set_facecolor('white')
numRows = 1
numCols = 1
#iteration = 0
numContours = 100
#ContourLevels = [2.7, 2.8, 2.9, 3.0, 3.1, 3.2]
#numLevels = 5
numShowLevels = 6

ShowColorbar = False

#cbarTicks = [2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75] # LHD discharge 113208 at t = 4.64 s
cbarTicks = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0] # W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#cbarTicks =None # None for automatic

#zMin = 2.2 # LHD discharge 113208 at t = 4.64 s
#zMax = 3.8 # LHD discharge 113208 at t = 4.64 s
zMin = 2.0 # W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
zMax = 3.0 # W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#zMin = None #None to get the default value 
#zMax = None #None to get the default value

##LOGARITHMIC Z-AXIS NOT IMPLEMENTED##

TickFormat = r'$%2.2f$'

zLabelPad = 30
ExtendRectangle = True

ContourLabelSize = 30
ContourLabelRemoveLine = True
ContourThickness = 5

ColorMap = 'gist_rainbow'
#ColorMap = 'hsv'
ColorMap = 'rainbow'

AddMaxMinBox = True
#MaxMinBoxXcoord = 0.415 # LHD discharge 113208 at t = 4.64 s
MaxMinBoxXcoord = 2.016*0.415 # W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
#MaxMinBoxYcoord = -1.185 # LHD discharge 113208 at t = 4.64 s
MaxMinBoxYcoord = -1.165 # W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
MaxMinBoxLabelSize = 40
MaxMinBoxFormat = '{:1.2f}'

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
   #return r'${}$'.format(x)
   return r'${}$'.format('{:1.2f}'.format(x))

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


if zMin == None:
   zMin = zFactor*numpy.amin(BHat[:,:])
if zMax == None:
   zMax = zFactor*numpy.amax(BHat[:,:])
print ("zMin = " + str(zMin))
print ("zMax = " + str(zMax))

zMinData = zFactor*numpy.amin(BHat[:,:])
zMaxData = zFactor*numpy.amax(BHat[:,:])
print ("zMinData = " + str(zMinData))
print ("zMaxData = " + str(zMaxData))

if zMax < zMin:
   print ("Error, zMax < zMin")
   sys.exit(1)


#delta = (numpy.amax(BHat) - numpy.amin(BHat)) / numLevels
#ContourLevels = numpy.arange(numpy.amin(BHat), numpy.amax(BHat) + delta/2.0, delta)
#ContourLevels = zFactor*ContourLevels

ContourLevels = numpy.linspace(zMin, zMax, num=numContours)
ShowLevels = numpy.linspace(zMin, zMax, num=numShowLevels)

print ("#############")
print ("ContourLevels = " + str( ContourLevels))
print ("")
print ("ShowLevels = " + str( ShowLevels))
print ("#############")
    
ax = plt.subplot(numRows,numCols,1)
    #plt.contourf(zeta,theta,1000*numpy.fliplr(BHat[:,:,iteration].transpose()),numContours)
    
#BPlot = plt.contourf(zeta,theta,zFactor*BHat.transpose(),numContours, cmap=plt.get_cmap(ColorMap))
BPlot = plt.contourf(zeta,theta,zFactor*BHat.transpose(),numContours, levels=ContourLevels, cmap=plt.get_cmap(ColorMap), vmin=zMin, vmax=zMax, extend='both')

#BPlot2 = plt.contour(BPlot,levels=ContourLevels, colors='k', hold='on')
BPlot2 = plt.contour(BPlot,levels=ShowLevels, colors='k', linewidths=ContourThickness)

plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$', fontsize=AxesLabelSize)
plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$', fontsize=AxesLabelSize)
#plt.zlabel(r'$B$'+ ' [T]')
plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)], fontsize=TickSize)
plt.yticks([0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)], fontsize=TickSize)
plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

#plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
#plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)
#plt.gca().axes.xaxis.set_label_coords(0.5,-0.05)
#plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

if show_rN:
    plt.title('rN = '+str(rN))

if ShowColorbar:

   #cbar = plt.colorbar(BPlot, label=r'$B$'+ ' [T]', ticks=ContourLevels)
   #cbar = plt.colorbar(BPlot, label=r'$\Phi_1$'+ ' [V]', ticks=BPlot.levels[::2])
   #cbar.add_lines(BPlot2)

   #cbar = plt.colorbar(BPlot, format=ticker.FuncFormatter(fmt_xy_axis), ticks=ContourLevels)
   cbar = plt.colorbar(BPlot, ticks=cbarTicks, format=TickFormat, extend='both', extendrect=ExtendRectangle)

   cbar.ax.set_ylabel(r'$B$'+ " " + r'$\mathrm{[T]}$', rotation=0, labelpad=zLabelPad, fontsize=AxesLabelSize)

#with warnings.catch_warnings():
#    warnings.simplefilter("always")
#plt.clabel(BPlot2, fmt='%2.1f', colors='k', fontsize=14)
#plt.clabel(BPlot2, fmt=ticker.FuncFormatter(fmt_xy_axis), colors='k', fontsize=18, inline=False)
plt.clabel(BPlot2, fmt=TickFormat, colors='k', fontsize=ContourLabelSize, inline=ContourLabelRemoveLine)


if ShowLegend:
    plt.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)#, fontsize=LegendFontSize)

plt.gca().axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
plt.gca().axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])

plt.tight_layout()

plt.subplots_adjust(left=LeftMargin, right=RightMargin, top=TopMargin, bottom=BottomMargin)

if ShowSubPlotLabel:
    plt.text(SubPlotLabelXcoord, SubPlotLabelYcoord, SubPlotLabel, fontsize=SubPlotLabelSize)

#plt.subplots_adjust(wspace=0.27)

if AddMaxMinBox:
   plt.text(MaxMinBoxXcoord, MaxMinBoxYcoord, r'$B_{\mathrm{min}} = $ ' + str(r'${}$'.format(MaxMinBoxFormat.format(zMinData))) + r' $\mathrm{T}$' + '\n' + r'$B_{\mathrm{max}} = $ ' + str(r'${}$'.format(MaxMinBoxFormat.format(zMaxData))) + r' $\mathrm{T}$', fontsize=MaxMinBoxLabelSize)

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
