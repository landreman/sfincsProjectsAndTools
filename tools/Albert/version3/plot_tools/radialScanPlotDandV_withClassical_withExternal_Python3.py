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

species = 3

withExternal = False
withClassicalExternal = False
externalDataFileType = '.dat'
radiusColumn = 0
aNorm = 0.51092 ##This is used to normalize if the same radial coordinate is not used in the external data as for the SFINCS results, e.g. r -> r/a. Put to 1.0 if same coordinate.
Dcolumn = 10
Vcolumn = 11
ClassicalDcolumn = 13
ClassicalVcolumn = 13
#FluxNorm = 10**20 ##Use if normalization in external fluxes is different than for the SFINCS results
#DensityColumn = 4

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
filename2 = 'sfincsOutput (2).h5' ##Name for second SFINCS output HDF5 files.

radiusName = "rN" ##Radial coordinate to use on x-axis. Must be "psiHat", "psiN", "rHat" or "rN".

TransformPlotVariableToOutputUnitsFactorD = vbar
TransformPlotVariableToOutputUnitsFactorV = vbar

MinFloat = pow(10, -sys.float_info.dig)

InputTol = 10**(-3)

##############################
##########END INPUTS##########
##############################

if radiusName != "psiHat" and radiusName != "psiN" and radiusName != "rHat" and radiusName != "rN":
    print ("Error! Invalid radial coordinate.")
    sys.exit(1)

##READ AND PLOT THE DATA##
originalDirectory = os.getcwd() 
print ("Starting to create a plot from directories in " + originalDirectory)

# Get a list of the subdirectories:
PlotDirectories = sorted(filter(os.path.isdir, os.listdir("."))) 

if len(PlotDirectories) < 1:
    print ("Error! Could not find any directories in " + originalDirectory)
    sys.exit(1)

#fig = plt.figure(figsize=FigSize)
fig, (axD, axV) = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, figsize=FigSize)
fig.patch.set_facecolor('white') 

#ax = plt.subplot(1, 1, 1)

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

        DneoData = []
        DclassicalData = []
        VneoData = []
        VclassicalData = []
        
        for SubDirectory in SubDirectories:
            fullSubDirectory = fullDirectory + "/" + SubDirectory
            os.chdir(fullSubDirectory)
            
            try:
                file = h5py.File(fullSubDirectory + "/" + filename,'r')
                file2 = h5py.File(fullSubDirectory + "/" + filename2,'r')
                radiusValue = file[radiusName][()]
                radiusValue2 = file2[radiusName][()]

                finished = file["finished"][()]
                finished2 = file2["finished"][()]
                integerToRepresentTrue = file["integerToRepresentTrue"][()]
                integerToRepresentTrue2 = file2["integerToRepresentTrue"][()]
                includePhi1 = file["includePhi1"][()]
                includePhi12 = file2["includePhi1"][()]

                assert ((finished == integerToRepresentTrue) and (finished2 == integerToRepresentTrue2)),"Error, both runs must finish!"
                assert ((includePhi1 == integerToRepresentTrue) == (includePhi12 == integerToRepresentTrue2)),"Error, seems like Phi1 is included in one calculation but not the other!"

                nHats = file["nHats"][()]
                nHats2 = file2["nHats"][()]
                THats = file["THats"][()]
                THats2 = file2["THats"][()]
                mHats = file["mHats"][()]
                mHats2 = file2["mHats"][()]
                Zs = file["Zs"][()]
                Zs2 = file2["Zs"][()]

                ##assert ((radiusValue - radiusValue2)/radiusValue <= InputTol),"Error, radii don't seem to match!"
                assert (numpy.allclose(radiusValue, radiusValue2, rtol=InputTol)),"Error, radii don't seem to match!"
                assert (numpy.allclose(nHats, nHats2, rtol=InputTol)),"Error, densities don't seem to match!"
                assert (numpy.allclose(THats, THats2, rtol=InputTol)),"Error, temperatures don't seem to match!"
                assert (numpy.allclose(mHats, mHats2, rtol=InputTol)),"Error, masses don't seem to match!"
                assert (numpy.allclose(Zs, Zs2, rtol=InputTol)),"Error, charges don't seem to match!"            

                if includePhi1 == integerToRepresentTrue:
                    didNonlinearCalculationConverge = file["didNonlinearCalculationConverge"][()]
                    didNonlinearCalculationConverge2 = file2["didNonlinearCalculationConverge"][()]
                    #if plotVariableName == "particleFlux_vm_rHat":
                    #    VariableValue = file["particleFlux_vd_rHat"][()]

                    assert ((didNonlinearCalculationConverge == integerToRepresentTrue) and (didNonlinearCalculationConverge2 == integerToRepresentTrue2)),"Error, both runs must finish!"

                particleFluxName = "particleFlux_vd_rHat"
                classicalparticleFluxName = "classicalParticleFlux_rHat"
                DensityGradientName = "dnHatdrHat"                 
                
                if (particleFluxName.find('Flux_vd') != -1 or particleFluxName.find('Flux_vE') != -1) and (includePhi1 != integerToRepresentTrue):
                    #print (particleFluxName + " only exists in nonlinear runs, but this is a linear run.")
                    particleFluxName = particleFluxName.replace('Flux_vd', 'Flux_vm').replace('Flux_vE', 'Flux_vm')
                    #print ("Reading " + particleFluxName + " instead.")
                    #VariableValue = file[plotVariableName.replace('Flux_vd', 'Flux_vm').replace('Flux_vE', 'Flux_vm')][()]
                #else:
                #    VariableValue = file[plotVariableName][()]

                classicalParticleFlux = file[classicalparticleFluxName][()]
                classicalParticleFlux = classicalParticleFlux[:, - 1]
                classicalParticleFlux = classicalParticleFlux[species - 1]

                classicalParticleFlux2 = file2[classicalparticleFluxName][()]
                classicalParticleFlux2 = classicalParticleFlux2[:, - 1]
                classicalParticleFlux2 = classicalParticleFlux2[species - 1]

                neoParticleFlux = file[particleFluxName][()]
                neoParticleFlux = neoParticleFlux[:, - 1]
                neoParticleFlux = neoParticleFlux[species - 1]

                neoParticleFlux2 = file2[particleFluxName][()]
                neoParticleFlux2 = neoParticleFlux2[:, - 1]
                neoParticleFlux2 = neoParticleFlux2[species - 1]

                densityGradient = file[DensityGradientName][()]
                densityGradient = densityGradient[species - 1]

                densityGradient2 = file2[DensityGradientName][()]
                densityGradient2 = densityGradient2[species - 1]
                nHat = nHats[species - 1]
                
                #classicalParticleFlux = TransformPlotVariableToOutputUnitsFactor * classicalParticleFlux
                #classicalParticleFlux = classicalParticleFlux / nHats[species - 1]

                #classicalParticleFlux2 = TransformPlotVariableToOutputUnitsFactor * classicalParticleFlux2
                #classicalParticleFlux2 = classicalParticleFlux2 / nHats2[species - 1] 

                file.close()
                file2.close()
                
                
                #if plotVariableName.find('Flux_v') != -1:
                #    VariableValue = VariableValue[:, -1]
                #    VariableValue = VariableValue[species - 1] 

                #VariableValue = TransformPlotVariableToOutputUnitsFactor * VariableValue
                    
            except AssertionError as err:
                print ("Error when reading from " + fullSubDirectory + "/" + filename + " and " + filename2)
                print (err.args[0])
                print ("Continuing with next sub directory.")
                continue
            except:
                print ("Error when reading from " + fullSubDirectory + "/" + filename + " and " + filename2)
                print ("Maybe the SFINCS run did not finish")
                print ("Continuing with next sub directory.")
                continue

            #print("nHat: " + str(nHats[species -1]))

            #VariableValue = VariableValue / nHats[species -1]

            Dneo = - (neoParticleFlux2 - neoParticleFlux) / (densityGradient2 - densityGradient)
            Dclassical = - (classicalParticleFlux2 - classicalParticleFlux) / (densityGradient2 - densityGradient)
            
            Vneo = (neoParticleFlux + Dneo * densityGradient) / nHat
            Vclassical = (classicalParticleFlux + Dclassical * densityGradient) / nHat

            Dneo = Dneo * TransformPlotVariableToOutputUnitsFactorD
            Dclassical = Dclassical * TransformPlotVariableToOutputUnitsFactorD
            
            Vneo = Vneo * TransformPlotVariableToOutputUnitsFactorV
            Vclassical = Vclassical * TransformPlotVariableToOutputUnitsFactorV

            Nradii += 1
            radii.append(radiusValue)
            
            #ydata.append(VariableValue)
            #ydata2.append(classicalParticleFlux)
            
            DneoData.append(Dneo)
            DclassicalData.append(Dclassical)
            VneoData.append(Vneo)
            VclassicalData.append(Vclassical)

        if Nradii < 1:
            print ("Could not read any data in " + fullDirectory) 
            print ("Continuing with next directory.") 
            continue

        ##Sort data after radii
        radii_sorted = sorted(radii)
        DneoData_sorted = []
        DclassicalData_sorted = []
        VneoData_sorted = []
        VclassicalData_sorted = []
        
        for radius in radii_sorted:
            DneoData_sorted.append(DneoData[radii.index(radius)])
            DclassicalData_sorted.append(DclassicalData[radii.index(radius)])
            VneoData_sorted.append(VneoData[radii.index(radius)])
            VclassicalData_sorted.append(VclassicalData[radii.index(radius)])
        
        print ("radii: " + str(radii))
        print ("")
        print ("DneoData: " + str(DneoData))
        print ("")        
        print ("DclassicalData: " + str(DclassicalData))
        print ("")
        print ("VneoData: " + str(VneoData))
        print ("")
        print ("VclassicalData: " + str(VclassicalData))
        print ("")
        print ("radii_sorted: " + str(radii_sorted))
        print ("")
        print ("DneoData_sorted: " + str(DneoData_sorted))
        print ("")        
        print ("DclassicalData_sorted: " + str(DclassicalData_sorted))
        print ("")
        print ("VneoData_sorted: " + str(VneoData_sorted))
        print ("")
        print ("VclassicalData_sorted: " + str(VclassicalData_sorted))

        #print (np.array(radii_sorted))
        #print ("")
        #print (np.array(ydata_sorted))
        #print ("")
        #print (np.array(ydata2_sorted))
        
        try:
            LegendLabel = PlotLegendLabels[linenumber]
        except:
            LegendLabel = directory
        
        #D neoclassical
        axD.plot(np.array(radii_sorted), np.array(DneoData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

        #V neoclassical
        axV.plot(np.array(radii_sorted), np.array(VneoData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)                
        
        linenumber += 1

        try:
            LegendLabel = PlotLegendLabels[linenumber]
        except:
            LegendLabel = directory

        #D neoclassical + classical
        axD.plot(np.array(radii_sorted), np.array(DneoData_sorted) + np.array(DclassicalData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

        #V neoclassical + classical
        axV.plot(np.array(radii_sorted), np.array(VneoData_sorted) + np.array(VclassicalData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
        
        linenumber += 1

        #try:
        #    LegendLabel = PlotLegendLabels[linenumber]
        #except:
        #    LegendLabel = directory
        #
        ##Neoclassical + classical flux
        #plt.plot(np.array(radii_sorted), np.array(ydata_sorted) + np.array(ydata2_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
        #linenumber += 1

    except:
        os.chdir(originalDirectory)
        print ("Unexpected error when processing " + directory)
        print ("Continuing with next directory.")
        continue


##ADD EXTERNAL DATA TO PLOT (E.G. DKES)##
if withExternal :
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
                plt.plot(inputParams[:,radiusColumn]/aNorm, inputParams[:,FluxColumn]/inputParams[:,DensityColumn]/FluxNorm, PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
                linenumber += 1

                if withClassicalExternal :
                    try:
                        LegendLabel = PlotLegendLabels[linenumber]
                    except:
                        LegendLabel = externalfile
                        plt.plot(inputParams[:,radiusColumn]/aNorm, inputParams[:,ClassicalFluxColumn]/inputParams[:,DensityColumn]/FluxNorm, PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
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

#plt.xscale(xAxisScale) 
#plt.yscale(yAxisScale)

axD.set_xscale(xAxisScale)
axD.set_yscale(yAxisScale)
axV.set_xscale(xAxisScale)
axV.set_yscale(yAxisScale)

if ShowGrid:
    #plt.grid(which='both')
    axD.grid(which='both')
    axV.grid(which='both')

#plt.xlabel(xAxisLabel, fontsize=AxesLabelSize)
#plt.ylabel(yAxisLabel, fontsize=AxesLabelSize)

axD.set_xlabel(xAxisLabel, fontsize=AxesLabelSize)
axD.set_ylabel(yAxisLabel[0], fontsize=AxesLabelSize)
axV.set_xlabel(xAxisLabel, fontsize=AxesLabelSize)
axV.set_ylabel(yAxisLabel[1], fontsize=AxesLabelSize)

if AxisLimAuto: 
    #ymin,ymax = plt.ylim()
    #xmin,xmax = plt.xlim()
    yminD,ymaxD = axD.get_ylim()
    xminD,xmaxD = axD.get_xlim()
    yminV,ymaxV = axV.get_ylim()
    xminV,xmaxV = axV.get_xlim()
else : 
    #ymin,ymax = plt.ylim(yAxisLim) 
    #xmin,xmax = plt.xlim(xAxisLim)
    yminD,ymaxD = plt.set_ylim(yAxisLim[0]) 
    xminD,xmaxD = plt.set_xlim(xAxisLim[0])
    yminV,ymaxV = plt.set_ylim(yAxisLim[1]) 
    xminV,xmaxV = plt.set_xlim(xAxisLim[1]) 

if ShowLegend:
    #plt.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)#, fontsize=LegendFontSize)
    axD.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)#, fontsize=LegendFontSize)

#plt.gca().axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
#plt.gca().axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])

axD.axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
axD.axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])
axV.axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
axV.axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])

plt.tight_layout()
#axD.tight_layout()
#axV.tight_layout()

plt.subplots_adjust(left=LeftMargin, right=RightMargin, top=TopMargin, bottom=BottomMargin)
#axD.subplots_adjust(left=LeftMargin, right=RightMargin, top=TopMargin, bottom=BottomMargin)
#axV.subplots_adjust(left=LeftMargin, right=RightMargin, top=TopMargin, bottom=BottomMargin)

if ShowSubPlotLabel:
    #plt.text(SubPlotLabelXcoord, SubPlotLabelYcoord, SubPlotLabel)
    axD.text(SubPlotLabelXcoord[0], SubPlotLabelYcoord[0], SubPlotLabel[0])
    axV.text(SubPlotLabelXcoord[1], SubPlotLabelYcoord[1], SubPlotLabel[1])

if NoScientificAxes :
    try:
        #ax.get_xaxis().get_major_formatter().set_scientific(False)
        #ax.get_yaxis().get_major_formatter().set_scientific(False)
        axD.get_xaxis().get_major_formatter().set_scientific(False)
        axD.get_yaxis().get_major_formatter().set_scientific(False)
        axV.get_xaxis().get_major_formatter().set_scientific(False)
        axV.get_yaxis().get_major_formatter().set_scientific(False)
    except:
        pass

os.chdir(originalDirectory) 

if makePDF: 
    print ("Saving PDF")  

    if len(sys.argv)>2 : #Use the substituted name as file name 
       print ("Writing plot to " + os.getcwd() + "/" + sys.argv[2] + ".pdf.")
       #plt.savefig(sys.argv[2] + ".pdf", orientation = 'landscape', papertype='letter')
       fig.savefig(sys.argv[2] + ".pdf", orientation = 'landscape', papertype='letter')
    else : 
       head, tail = os.path.split(inspect.getfile(inspect.currentframe()))
       print ("Writing plot to " + os.getcwd() + "/" + tail + ".pdf.")
       #plt.savefig(tail+'.pdf', orientation = 'landscape', papertype='letter')
       fig.savefig(tail+'.pdf', orientation = 'landscape', papertype='letter')
else:   
    #plt.show()
    fig.show()
