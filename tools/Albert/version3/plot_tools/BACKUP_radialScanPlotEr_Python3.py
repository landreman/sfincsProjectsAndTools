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

print("This is "+ inspect.getfile(inspect.currentframe()))

sfincsHome = os.environ.get('SFINCS_HOME') 
sfincsProjectsAndToolsHome = os.environ.get('SFINCS_PROJECTS_AND_TOOLS_HOME')

##INPUTS##

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

filename = 'sfincsOutput.h5' ##Name for SFINCS output HDF5 files.

radiusName = "rN" ##Radial coordinate to use on x-axis. Must be "psiHat", "psiN", "rHat" or "rN".

plotVariableName = "Er" ##Parameter to plot on y-axis. In this version it must be "Er", "dPhiHatdpsiHat", "dPhiHatdpsiN", "dPhiHatdrHat" or "dPhiHatdrN" .

MinFloat = pow(10, -sys.float_info.dig) 

##PLOT OPTIONS##

exec(open(sfincsProjectsAndToolsHome + "/tools/Albert/version3/plot_tools"  + "/RadialScanPlotOptions.py").read())

################



##############################
##########END INPUTS##########
##############################

if radiusName != "psiHat" and radiusName != "psiN" and radiusName != "rHat" and radiusName != "rN":
    print ("Error! Invalid radial coordinate.")
    sys.exit(1)

if plotVariableName != "Er" and plotVariableName != "dPhiHatdpsiHat" and plotVariableName != "dPhiHatdpsiN" and plotVariableName != "dPhiHatdrHat" and plotVariableName != "dPhiHatdrN":
    print ("Error! Invalid variable name.")
    sys.exit(1)




##READ AND PLOT THE DATA##
originalDirectory = os.getcwd() 
print ("Starting to create a plot from directories in " + originalDirectory)

# Get a list of the subdirectories:
PlotDirectories = sorted(filter(os.path.isdir, os.listdir("."))) 

if len(PlotDirectories) < 1:
    print ("Error! Could not find any directories in " + originalDirectory) 
    sys.exit(1)

fig = plt.figure(figsize=FigSize) 
fig.patch.set_facecolor('white') 

ax = plt.subplot(1, 1, 1)

linenumber = 0

for directory in PlotDirectories:
    try:
        print ("*************************************************") 
        print ("Processing directory "+directory)
        print ("*************************************************")

        fullDirectory = originalDirectory + "/" + directory
        os.chdir(fullDirectory)

        # Get a list of the subdirectories containing the SFINCS output
        SubDirectories = sorted(filter(os.path.isdir, os.listdir(".")))

        if len(SubDirectories) < 1:
            print ("Could not find any directories in " + fullDirectory)
            print ("Continuing with next directory.") 
            continue

        Nradii = 0
        radii = []
        ydata = []
        
        for SubDirectory in SubDirectories:
            fullSubDirectory = fullDirectory + "/" + SubDirectory
            os.chdir(fullSubDirectory)
            
            try:
                file = h5py.File(fullSubDirectory + "/" + filename,'r')
                radiusValue = file[radiusName][()]
                VariableValue = file[plotVariableName][()] 

                finished = file["finished"][()] 
                integerToRepresentTrue = file["integerToRepresentTrue"][()]
                includePhi1 = file["includePhi1"][()] 

                if includePhi1 == integerToRepresentTrue:
                    didNonlinearCalculationConverge = file["didNonlinearCalculationConverge"][()]

                file.close()
                
                if includePhi1 == integerToRepresentTrue:
                    if didNonlinearCalculationConverge != integerToRepresentTrue:
                        print ("The nonlinear solver did not converge in " + fullSubDirectory)
                        print ("Continuing with next sub directory.")
                        continue
            except:
                print ("Error when reading from " + fullSubDirectory + "/" + filename)
                print ("Maybe the SFINCS run did not finish")
                print ("Continuing with next sub directory.")
                continue

            Nradii += 1
            radii.append(radiusValue)
            ydata.append(VariableValue)

        if Nradii < 1:
            print ("Could not read any data in " + fullDirectory) 
            print ("Continuing with next directory.") 
            continue

        ##Sort data after radii
        radii_sorted = sorted(radii)
        ydata_sorted = []
        for radius in radii_sorted:
            ydata_sorted.append(ydata[radii.index(radius)])
        
        print ("radii: " + str(radii))
        print ("")
        print ("ydata: " + str(ydata))
        print ("")
        print ("radii_sorted: " + str(radii_sorted))
        print ("")
        print ("ydata_sorted: " + str(ydata_sorted))
        print ("")

        print (np.array(radii_sorted))
        print ("")
        print (np.array(ydata_sorted))

        try:
            LegendLabel = PlotLegendLabels[linenumber]
        except:
            LegendLabel = directory

        plt.plot(np.array(radii_sorted), np.array(ydata_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel)
        linenumber += 1

    except:
        os.chdir(originalDirectory)
        print ("Unexpected error when processing " + directory)
        print ("Continuing with next directory.")
        continue


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

plt.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., fontsize=LegendFontSize)

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