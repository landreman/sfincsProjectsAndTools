#!/usr/bin/env python

import matplotlib
#import matplotlib.pyplot as plt
import h5py
import numpy
import os, sys, inspect
import warnings
import matplotlib.ticker as ticker

from VMECtoPEST_functions import interp2_cyclic, griddatacyclic
from transformVMECtoPEST_Python3 import transformVMECtoPEST

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

sfincsHome = os.environ.get('SFINCS_HOME')
sfincsProjectsAndToolsHome = os.environ.get('SFINCS_PROJECTS_AND_TOOLS_HOME')

exec(open(sfincsProjectsAndToolsHome + "/tools/Albert/version3/plot_tools"  + "/RadialScanPlotOptions.py").read())

#########
##INPUT##
#########

quantityToPlot = "Phi1Hat"

#filename = 'sfincsOutput.h5'

#W7X_180919.055
#filename = 'SFINCS_kinetic-e_Fokker-Planck_r0p1148.h5'
filename = 'SFINCS_kinetic-e_Fokker-Planck_r0p220343.h5'

#ncFilename = "/draco/u/almo/Phi1/LHD/lhd2_A_III/Input/wout_lhd2.nc"
#ncFilename = "C:/Users/almo/Desktop/svn/sfincs/Impurities/Phi1/Results/LHD_Velasco_PPCF18/input/wout_lhd_r3.60_0.0.nc"
#ncFilename = "C:/Users/almo/Desktop/svn/sfincs/TJ-II/Input_TJII_case_Regana_NF17/wout_tj20.nc"
#ncFilename = "C:/Users/almo/Desktop/svn/sfincs/W7-X/OP1.1/XICS_data_Langenberg/Equilibria/wout_w7x.1000_1000_1000_1000_+0390_+0000.05.0000.nc"
#ncFilename = "C:/Users/almo/Desktop/svn/sfincs/Impurities/Phi1/Results/LHD_Velasco_PPCF18/input/wout_lhd_r3.60_0.0.nc"

#W7X_180919.055
ncFilename = "C:/Users/almo/Desktop/svn/sfincs/W7-X/OP1.2/W7X_180919.055/Equilibrium/wout_w7x.1000_1000_1000_1000_+0000_+0000.01.00jh_l+252.nc"

#FigSize = (12,10)
FigSize = (10,12)
LeftMargin = 0.01
RightMargin = 0.5
TopMargin = 0.99
BottomMargin = 0.01
zLabelPad = 45
PhiLabelPad = 60

#font = {'size':25}
#matplotlib.rc('font', **font)
#matplotlib.rc('lines',markeredgewidth=0,markersize=3,linewidth=2.5)
#matplotlib.rc('axes',linewidth=1.5)

zFactor = 1000 ##kV -> V
##W7-X##
xAxisTicks = [r'$0$', r'$\pi/10$', r'$2\pi/10$', r'$3\pi/10$', r'$4\pi/10$']
##LHD
#xAxisTicks = [r'$0$', r'$\pi/20$', r'$2\pi/20$', r'$3\pi/20$', r'$4\pi/20$']
##TJ-II
#xAxisTicks = [r'$0$', r'$\pi/8$', r'$2\pi/8$', r'$3\pi/8$', r'$4\pi/8$']

yAxisTicks = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

fig = plt.figure(figsize=FigSize)
fig.patch.set_facecolor('white')
numRows = 1
numCols = 1
#iteration = 0
numContours = 100
#ContourLevels = [-3.0, -1.5, 0.0, 1.5, 3.0, 4.5, 6.0]
#numLevels = 5
numShowLevels = 10

ShowColorbar = True

#cbarTicks = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0] # TJ-II TTF poster r/a=0.6
#cbarTicks = [-12.0, -8.0, -4.0, 0.0, 4.0, 8.0, 12.0] # W7-X OP1.1 TTF poster r/a=0.5
#cbarTicks = [-20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0] # LHD inward shifted TTF poster r/a=0.8
#cbarTicks = [-40.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0] # LHD inward shifted TTF poster r/a=0.8 with magnetic drifts
#cbarTicks = [-10.0, -7.5, -5.0, -2.5, 0.0, 2.5, 5.0, 7.5, 10.0] # W7X_180919.055 r/a=0.11
cbarTicks = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0] # W7X_180919.055 r/a=0.21

#zMin = -2.0 # TJ-II TTF poster r/a=0.6
#zMax = 2.0 # TJ-II TTF poster r/a=0.6
#zMin = -12.0 # W7-X OP1.1 TTF poster r/a=0.5
#zMax = 12.0 # W7-X OP1.1 TTF poster r/a=0.5
#zMin = -20.0 # LHD inward shifted TTF poster r/a=0.8
#zMax = 10.0 # LHD inward shifted TTF poster r/a=0.8
#zMin = -45.0 # LHD inward shifted TTF poster r/a=0.8 with magnetic drifts
#zMax = 25.0 # LHD inward shifted TTF poster r/a=0.8 with magnetic drifts
#zMin = -10.0 # W7X_180919.055 r/a=0.11
#zMax = 10.0 # W7X_180919.055 r/a=0.11
zMin = -3.0 # W7X_180919.055 r/a=0.21
zMax = 3.0 # W7X_180919.055 r/a=0.21

zLogAxis = False
LinearThreshold = 1.0 #In symlog plot
LinearScale = 1.0 #When linscale == 1.0 (the default), the space used for the positive and negative halves of the linear range will be equal to one decade in the logarithmic range.

TickFormat = r'$%2.1f$'

ExtendRectangle = True

ContourLabelSize = 30
ContourLabelRemoveLine = True
ContourThickness = 5

ContourLabelShowLinearBorder = False
#ContourBorderColor = '#cccccc'
ContourBorderColor = '#EEEE33'
ContourBorderLineWidth = 4.

#ColorMap = 'gist_rainbow'
#ColorMap = 'hsv'
ColorMap = 'rainbow'

AddMaxMinBox = False
MaxMinBoxXcoord = 0.4177 # LHD discharge 113208 at t = 4.64 s
#MaxMinBoxXcoord = 2.016*0.4177 # W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
MaxMinBoxYcoord = -1.25 # LHD discharge 113208 at t = 4.64 s
#MaxMinBoxYcoord = -1.224 # W7-X_NBI_case_Q34Q78_Z10_Zeff2p0
MaxMinBoxLabelSize = 40
MaxMinBoxFormat = '{:1.1f}'

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
theta = f["theta"][()]
zeta = f["zeta"][()]
iteration = f["NIterations"][()] - 1 #Results from last iteration
Phi1Hat = ((f[quantityToPlot][()])[:,:,iteration]).transpose()
rN = f["rN"][()]

############################################
##Transform to PEST grid
Ntheta = f["Ntheta"][()]
Nzeta = f["Nzeta"][()]
psiN = f["psiN"][()]
f.close()
Pest_vmecu, Pest_vmecw, Vmec_vmecu, Vmec_vmecw, Geom_Nperiods = transformVMECtoPEST(ncFilename, psiN, theta, zeta, -1)

Pest_Phi1 = interp2_cyclic(Vmec_vmecu, Vmec_vmecw, Phi1Hat, Pest_vmecu, Pest_vmecw, Geom_Nperiods)
############################################

print ("psiN: " + str(psiN))
print ("theta max: " + str(numpy.amax(theta)))
print ("theta min: " + str(numpy.amin(theta)))
print ("zeta max: " + str(numpy.amax(zeta)))
print ("zeta min: " + str(numpy.amin(zeta)))
print ("Pest_vmecu max: " + str(numpy.amax(Pest_vmecu)))
print ("Pest_vmecu min: " + str(numpy.amin(Pest_vmecu)))
print ("Pest_vmecw max: " + str(numpy.amax(Pest_vmecw)))
print ("Pest_vmecw min: " + str(numpy.amin(Pest_vmecw)))

print ("Phi1 max: " + str(numpy.amax(Phi1Hat)))
print ("Phi1 min: " + str(numpy.amin(Phi1Hat)))
print ("PEST_Phi1 max: " + str(numpy.amax(Pest_Phi1)))
print ("PEST_Phi1 min: " + str(numpy.amin(Pest_Phi1)))


#delta = (numpy.amax(Pest_Phi1) - numpy.amin(Pest_Phi1)) / numLevels
#ContourLevels = numpy.arange(numpy.amin(Pest_Phi1), numpy.amax(Pest_Phi1) + delta/2.0, delta)
#ContourLevels = zFactor*ContourLevels

if zMin == None:
   zMin = zFactor*numpy.amin(Pest_Phi1)
if zMax == None:
   zMax = zFactor*numpy.amax(Pest_Phi1)
print ("zMin = " + str(zMin))
print ("zMax = " + str(zMax))

zMinData = zFactor*numpy.amin(Pest_Phi1)
zMaxData = zFactor*numpy.amax(Pest_Phi1)
print ("zMinData = " + str(zMinData))
print ("zMaxData = " + str(zMaxData))
#sys.exit(0)

if zMax < zMin:
   print ("Error, zMax < zMin")
   sys.exit(1)

SymLogFlag = False

if zLogAxis:
   if zMin > 0:
      ContourLevels = numpy.logspace(numpy.log10(numpy.absolute(zMin)), numpy.log10(numpy.absolute(zMax)), num=numpy.ceil(numContours/2.0), endpoint=True, base=10.0, dtype=None)/zFactor
      ShowLevels = numpy.logspace(numpy.log10(numpy.absolute(zMin)), numpy.log10(numpy.absolute(zMax)), num=numpy.ceil(numShowLevels/2.0), endpoint=True, base=10.0, dtype=None)/zFactor
   elif zMax < 0:
      ContourLevels = - numpy.logspace(numpy.log10(numpy.absolute(zMax)), numpy.log10(numpy.absolute(zMin)), num=numpy.ceil(numContours/2.0), endpoint=True, base=10.0, dtype=None)/zFactor
      ShowLevels = - numpy.logspace(numpy.log10(numpy.absolute(zMax)), numpy.log10(numpy.absolute(zMin)), num=numpy.ceil(numShowLevels/2.0), endpoint=True, base=10.0, dtype=None)/zFactor
      ContourLevels = numpy.flip(ContourLevels, 0)
      ShowLevels = numpy.flip(ShowLevels, 0)
   else:
      SymLogFlag = True
      LinearThreshold = numpy.absolute(LinearThreshold)
      LinearScale = numpy.absolute(LinearScale)

      NumberOfPositiveDecades = numpy.amax([numpy.log10(numpy.absolute(zMax)) - numpy.log10(LinearThreshold), 0.0])
      NumberOfNegativeDecades = numpy.amax([numpy.log10(numpy.absolute(zMin)) - numpy.log10(LinearThreshold), 0.0])
      NumberOfLinearDecades = LinearScale

      FractionOfPositiveDecades = NumberOfPositiveDecades / (NumberOfPositiveDecades + NumberOfNegativeDecades +  NumberOfLinearDecades)
      FractionOfNegativeDecades = NumberOfNegativeDecades / (NumberOfPositiveDecades + NumberOfNegativeDecades +  NumberOfLinearDecades)
      FractionOfLinearDecades = NumberOfLinearDecades / (NumberOfPositiveDecades + NumberOfNegativeDecades +  NumberOfLinearDecades)
      
      print ("NumberOfPositiveDecades = " + str(NumberOfPositiveDecades))
      print ("NumberOfNegativeDecades = " + str(NumberOfNegativeDecades))
      print ("NumberOfLinearDecades = " + str(NumberOfLinearDecades))

      print ("FractionOfPositiveDecades = " + str(FractionOfPositiveDecades))
      print ("FractionOfNegativeDecades = " + str(FractionOfNegativeDecades))
      print ("FractionOfLinearDecades = " + str(FractionOfLinearDecades))

      numPositiveContours = numpy.rint(FractionOfPositiveDecades*numContours)
      numNegativeContours = numpy.rint(FractionOfNegativeDecades*numContours)
      numLinearContours = numContours - numPositiveContours - numNegativeContours

      numPositiveShowLevels = numpy.rint(FractionOfPositiveDecades*numShowLevels)
      numNegativeShowLevels = numpy.rint(FractionOfNegativeDecades*numShowLevels)
      numLinearShowLevels = numShowLevels - numPositiveShowLevels - numNegativeShowLevels

      print ("numPositiveContours = " + str(numPositiveContours))
      print ("numNegativeContours = " + str(numNegativeContours))
      print ("numLinearContours = " + str(numLinearContours))

      print ("numPositiveShowLevels = " + str(numPositiveShowLevels))
      print ("numNegativeShowLevels = " + str(numNegativeShowLevels))
      print ("numLinearShowLevels = " + str(numLinearShowLevels))
      
      ContourLevelsPos = numpy.logspace(numpy.log10(LinearThreshold), numpy.log10(numpy.absolute(zMax)), num=numPositiveContours + 1, endpoint=True, base=10.0, dtype=None)/zFactor
      if ContourLevelsPos.size > 0:
         ContourLevelsPos = numpy.delete(ContourLevelsPos, 0)
      ContourLevelsNeg = - numpy.logspace(numpy.log10(LinearThreshold), numpy.log10(numpy.absolute(zMin)), num=numNegativeContours + 1, endpoint=True, base=10.0, dtype=None)/zFactor
      if ContourLevelsNeg.size > 0:
         ContourLevelsNeg = numpy.delete(ContourLevelsNeg, 0)
      ContourLevelsLin = numpy.linspace(-LinearThreshold, LinearThreshold, num=numLinearContours)/zFactor
      ContourLevels = numpy.append(numpy.append(numpy.flip(ContourLevelsNeg, 0), ContourLevelsLin), ContourLevelsPos)

      ShowLevelsPos = numpy.logspace(numpy.log10(LinearThreshold), numpy.log10(numpy.absolute(zMax)), num=numPositiveShowLevels + 1, endpoint=True, base=10.0, dtype=None)/zFactor
      if ShowLevelsPos.size > 0:
         ShowLevelsPos = numpy.delete(ShowLevelsPos, 0)
      ShowLevelsNeg = - numpy.logspace(numpy.log10(LinearThreshold), numpy.log10(numpy.absolute(zMin)), num=numNegativeShowLevels + 1, endpoint=True, base=10.0, dtype=None)/zFactor
      if ShowLevelsNeg.size > 0:
         ShowLevelsNeg = numpy.delete(ShowLevelsNeg, 0)
      ShowLevelsLin = numpy.linspace(-LinearThreshold, LinearThreshold, num=numLinearShowLevels)/zFactor
      if ShowLevelsLin.size > 0:
         ShowLevelsLin = numpy.delete(ShowLevelsLin, 0)
      if ShowLevelsLin.size > 0:
         ShowLevelsLin = numpy.delete(ShowLevelsLin, -1)
      ShowLevels = numpy.append(numpy.append(numpy.flip(ShowLevelsNeg, 0), ShowLevelsLin), ShowLevelsPos)
      
else:
   ContourLevels = numpy.linspace(zMin, zMax, num=numContours)/zFactor
   ShowLevels = numpy.linspace(zMin, zMax, num=numShowLevels)/zFactor

ContourLevels = zFactor*ContourLevels
ShowLevels = zFactor*ShowLevels

print ("#############")
print ("ContourLevels = " + str( ContourLevels))
print ("")
print ("ShowLevels = " + str( ShowLevels))
print ("#############")

LogNormToUse = matplotlib.colors.SymLogNorm(linthresh=LinearThreshold, linscale=LinearScale, vmin=zMin, vmax=zMax)


    
ax = plt.subplot(numRows,numCols,1)
    #plt.contourf(zeta,theta,1000*numpy.fliplr(Phi1Hat[:,:,iteration].transpose()),numContours)

if zLogAxis:
   Phi1Plot = plt.contourf(Vmec_vmecw.transpose(), Vmec_vmecu.transpose(), zFactor*Pest_Phi1.transpose(), numContours, levels=ContourLevels, cmap=plt.get_cmap(ColorMap), norm=LogNormToUse, extend='both')
else:
   Phi1Plot = plt.contourf(Vmec_vmecw.transpose(), Vmec_vmecu.transpose(), zFactor*Pest_Phi1.transpose(), numContours, levels=ContourLevels, cmap=plt.get_cmap(ColorMap), vmin=zMin, vmax=zMax, extend='both')

ax.set_visible(False)
    
#Phi1Plot = plt.contourf(Vmec_vmecw.transpose(), Vmec_vmecu.transpose(), zFactor*Pest_Phi1.transpose(),numContours, cmap=plt.get_cmap(ColorMap))
#Phi1Plot2 = plt.contour(Phi1Plot,levels=ContourLevels, colors='k', hold='on')
#plt.xlabel(r'$\zeta$' + " " + r'$\mathrm{[rad]}$')
#plt.ylabel(r'$\theta$'+ " " + r'$\mathrm{[rad]}$')

##plt.xticks([0,max(zeta)/4,max(zeta)/2,3*max(zeta)/4,max(zeta)])
##plt.yticks([0.0,max(theta)/4,max(theta)/2,3*max(theta)/4,max(theta)])
#plt.xticks([0,numpy.amax(Vmec_vmecw)/4,numpy.amax(Vmec_vmecw)/2,3*numpy.amax(Vmec_vmecw)/4,numpy.amax(Vmec_vmecw)])
#plt.yticks([0,numpy.amax(Vmec_vmecu)/4,numpy.amax(Vmec_vmecu)/2,3*numpy.amax(Vmec_vmecu)/4,numpy.amax(Vmec_vmecu)])
#plt.gca().axes.xaxis.set_ticklabels(xAxisTicks)
#plt.gca().axes.yaxis.set_ticklabels(yAxisTicks)

##plt.gca().axes.xaxis.set_label_coords(0.5,-0.09)
##plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)
#plt.gca().axes.xaxis.set_label_coords(0.5,-0.05)
#plt.gca().axes.yaxis.set_label_coords(-0.09,0.5)

##ax.xaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis))
##ax.xaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))
##ax.yaxis.set_major_formatter( ticker.FuncFormatter(fmt_xy_axis)) 
##ax.yaxis.set_minor_formatter( ticker.FuncFormatter(fmt_xy_axis))

##plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

if show_rN:
    plt.title('rN = '+str(rN))

##cbar = plt.colorbar(Phi1Plot, label=r'$\Phi_1$'+ ' [V]', ticks=ContourLevels)
##cbar = plt.colorbar(Phi1Plot, label=r'$\Phi_1$'+ ' [V]', ticks=Phi1Plot.levels[::2])
##cbar.add_lines(Phi1Plot2)
#cbar = plt.colorbar(Phi1Plot, format=ticker.FuncFormatter(fmt_cbar), ticks=ContourLevels)
#cbar.ax.set_ylabel(r'$\Phi_1$'+ " " + r'$\mathrm{[V]}$', rotation=0, labelpad=zLabelPad)

if ShowColorbar:

   cbar = plt.colorbar(Phi1Plot, ticks=cbarTicks, format=TickFormat, extend='both', extendrect=ExtendRectangle)

   #cbar.add_lines(Phi1Plot2)
   #cbar.ax.set_yticklabels(cbarTicks)
   
   cbar.ax.set_ylabel(r'$\Phi_1$'+ " " + r'$\mathrm{[V]}$', rotation=0, labelpad=PhiLabelPad, fontsize=AxesLabelSize)

##with warnings.catch_warnings():
##    warnings.simplefilter("always")
##plt.clabel(Phi1Plot2, fmt='%2.1f', colors='k', fontsize=14)
#plt.clabel(Phi1Plot2, fmt=ticker.FuncFormatter(fmt_cbar), colors='k', fontsize=18, inline=False)

#plt.subplots_adjust(wspace=0.27)

if SymLogFlag:
   #if ContourLabelShowLinearBorder:
   #   plt.clabel(Phi1Plot3, fmt=TickFormat, colors=ContourBorderColor, fontsize=ContourLabelSize, inline=ContourLabelRemoveLine)

   if ShowColorbar:
      #cbar.add_lines(Phi1Plot3)
      #cbar.ax.axhline(y=FractionOfNegativeDecades + 0.0/numContours,linewidth=4.)
      cbar.ax.axhline(y=1.0*numNegativeContours/numContours + 0.5/numContours, linewidth=ContourBorderLineWidth, color=ContourBorderColor, linestyle='dashed')
      #numNegativeContours
      cbar.ax.axhline(y=1.0 - 1.0*numPositiveContours/numContours - 0.5/numContours, linewidth=ContourBorderLineWidth, color=ContourBorderColor)
      cbar.ax.text(-0.5, (numLinearContours/2.0 + 1.0*numNegativeContours + 0.5/numContours)/numContours, 'Linear', rotation='vertical', rotation_mode='anchor', color=ContourBorderColor, horizontalalignment='center')

if ShowLegend:
    plt.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)#, fontsize=LegendFontSize)

plt.gca().axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
plt.gca().axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])

plt.tight_layout()

plt.subplots_adjust(left=LeftMargin, right=RightMargin, top=TopMargin, bottom=BottomMargin)

#plt.subplots_adjust(wspace=0.27)

if ShowSubPlotLabel:
    plt.text(SubPlotLabelXcoord, SubPlotLabelYcoord, SubPlotLabel, fontsize=SubPlotLabelSize)

if AddMaxMinBox:
   plt.text(MaxMinBoxXcoord, MaxMinBoxYcoord, r'$\Phi_{1}^{\mathrm{min}} = $ ' + str(r'${}$'.format(MaxMinBoxFormat.format(zMinData))) + r' $\mathrm{V}$' + '\n' + r'$\Phi_{1}^{\mathrm{max}} = $ ' + str(r'${}$'.format(MaxMinBoxFormat.format(zMaxData))) + r' $\mathrm{V}$', fontsize=MaxMinBoxLabelSize)

print ("PEST_Phi1 shape: " + str(Pest_Phi1.shape))

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
