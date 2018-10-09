#!/usr/bin/env python 

import matplotlib 
import h5py
import numpy 
import numpy as np 
import os, sys, inspect, math
import warnings  
import matplotlib.ticker as ticker  
import subprocess

makePDF = False 
for arg in sys.argv:
    if arg.lower()=='pdf':
        makePDF = True

if makePDF:
    matplotlib.use('PDF')

import matplotlib.pyplot as plt

print ("This is "+ inspect.getfile(inspect.currentframe()))

sfincsHome = os.environ.get('SFINCS_HOME') 
sfincsProjectsAndToolsHome = os.environ.get('SFINCS_PROJECTS_AND_TOOLS_HOME')

##PLOT OPTIONS##

#execfile(sfincsProjectsAndToolsHome + "/tools/Albert/version3/plot_tools"  + "/RadialScanPlotOptions.py")
exec(open(sfincsProjectsAndToolsHome + "/tools/Albert/version3/plot_tools"  + "/RadialScanPlotOptions.py").read())

################

##INPUTS##

externalDataFileType = '.newdat'
radiusColumn = 0
PlotColumn = 1

##NORMALIZATION FACTORS FOR SI UNITS##
######################################
qe = 1.6021766208*10**(-19) #Electron charge
mbar = 1.672621777*10**(-27) #Proton mass 
nbar = 10.0**20 # in 10^20 
Tbar = 1000*qe 
Rbar = 1.0  
Bbar = 1.0  
vbar = np.sqrt(2.0 * Tbar / mbar)  
######################################

externalNormalization = 1.0

MinFloat = pow(10, -sys.float_info.dig) 

##############################
##########END INPUTS##########
##############################

##READ AND PLOT THE DATA##
originalDirectory = os.getcwd() 
print ("Starting to create a plot in " + originalDirectory)


fig = plt.figure(figsize=FigSize) 
fig.patch.set_facecolor('white') 

ax = plt.subplot(1, 1, 1)

linenumber = 0

##ADD EXTERNAL DATA TO PLOT (E.G. DKES)##
os.chdir(originalDirectory)
externalInputFiles = [];

for externalfile in os.listdir(originalDirectory):
    if externalfile.endswith(externalDataFileType):
        try:
            inputParams = np.genfromtxt(externalfile, dtype=None, comments="#", skip_header=1)
            #print(inputParams[:,radiusColumn])
            #print(inputParams[:,ErColumn])
            try:
                LegendLabel = PlotLegendLabels[linenumber]
            except:
                LegendLabel = externalfile
            plt.plot(inputParams[:,radiusColumn], inputParams[:,PlotColumn]/externalNormalization, PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
            linenumber += 1
            
            externalInputFiles.append(externalfile)
        except Exception as e:
            print (e.__class__.__name__, ": ", e.message, "while reading %s" % inputfile)
            print ("Continuing with next file!")
            continue

if len(externalInputFiles) > 0 :
    print("Read external data files: " + str(externalInputFiles))
else:
    print("Could not read any external data files.")

#########################################


plt.xscale(xAxisScale) 
plt.yscale(yAxisScale)

if ShowGrid:
    plt.grid(which='both')

plt.xlabel(xAxisLabel, fontsize=AxesLabelSize)
plt.ylabel(yAxisLabel, fontsize=AxesLabelSize) 

if AxisLimAuto: 
    ymin,ymax = plt.ylim()
    xmin,xmax = plt.xlim() 
else : 
    ymin,ymax = plt.ylim(yAxisLim) 
    xmin,xmax = plt.xlim(xAxisLim) 

if ShowLegend:
    plt.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)#, fontsize=LegendFontSize)

plt.gca().axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
plt.gca().axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])

plt.tight_layout()

plt.subplots_adjust(left=LeftMargin, right=RightMargin, top=TopMargin, bottom=BottomMargin)

if ShowSubPlotLabel:
    plt.text(SubPlotLabelXcoord, SubPlotLabelYcoord, SubPlotLabel)

if NoScientificAxes :
    ax.get_xaxis().get_major_formatter().set_scientific(False)
    ax.get_yaxis().get_major_formatter().set_scientific(False)

os.chdir(originalDirectory) 

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
