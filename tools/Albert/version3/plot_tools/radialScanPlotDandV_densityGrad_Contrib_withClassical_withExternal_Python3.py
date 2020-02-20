#!/usr/bin/env python 

import matplotlib 
import h5py
import numpy 
import numpy as np 
import os, sys, inspect, math
import warnings  
import matplotlib.ticker as ticker  
import subprocess
from scipy import interpolate

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

SFINCSplotWithClassical = [True] #Whether to include classical fluxes in SFINCS results
SFINCSplotFactors = [1.0] #Factor to multiply SFINCS results

legendAtDplot = True #If True the legend will be in the D plot, otherwise in the V plot

withExternal = True
withClassicalExternal = [False]
externalDataFileType = '.dat'
radiusColumn = [0]
aNorm = 1.0 #0.51092 ##This is used to normalize if the same radial coordinate is not used in the external data as for the SFINCS results, e.g. r -> r/a. Put to 1.0 if same coordinate.
Dcolumn = [3]
Vcolumn = [5]
sigmaDColumn = [4] #Put -1 if no sigma
sigmaVColumn = [6] #Put -1 if no sigma
classicalDcolumn = [-1] #Put -1 if no classical
classicalVcolumn = [-1] #Put -1 if no classical
sigmaClassicalDcolumn = [-1] #Put -1 if no classical sigma
sigmaClassicalVcolumn = [-1] #Put -1 if no classical sigma
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

filename1 = 'sfincsOutput.h5' ##Name for SFINCS output HDF5 files.
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
fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, figsize=FigSize)
axD = axes[0, 0]
axV = axes[0, 1]
axDensityGradient = axes[1, 0]
axDgrad = axes[1, 1]
fig.patch.set_facecolor('white')

LogarithmicDensityGradientData  = []
LogarithmicDensityGradientData_sorted = []


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
                file1 = h5py.File(fullSubDirectory + "/" + filename1,'r')
                file2 = h5py.File(fullSubDirectory + "/" + filename2,'r')
                radiusValue = file1[radiusName][()]
                radiusValue2 = file2[radiusName][()]

                finished = file1["finished"][()]
                finished2 = file2["finished"][()]
                integerToRepresentTrue = file1["integerToRepresentTrue"][()]
                integerToRepresentTrue2 = file2["integerToRepresentTrue"][()]
                includePhi1 = file1["includePhi1"][()]
                includePhi12 = file2["includePhi1"][()]

                assert ((finished == integerToRepresentTrue) and (finished2 == integerToRepresentTrue2)),"Error, both runs must finish!"
                assert ((includePhi1 == integerToRepresentTrue) == (includePhi12 == integerToRepresentTrue2)),"Error, seems like Phi1 is included in one calculation but not the other!"

                nHats = file1["nHats"][()]
                nHats2 = file2["nHats"][()]
                THats = file1["THats"][()]
                THats2 = file2["THats"][()]
                mHats = file1["mHats"][()]
                mHats2 = file2["mHats"][()]
                Zs = file1["Zs"][()]
                Zs2 = file2["Zs"][()]

                ##assert ((radiusValue - radiusValue2)/radiusValue <= InputTol),"Error, radii don't seem to match!"
                assert (numpy.allclose(radiusValue, radiusValue2, rtol=InputTol)),"Error, radii don't seem to match!"
                ##Note: density of studied species can be different since d(ln n)/dr is the thermodynamic force!
                assert (numpy.allclose(numpy.delete(nHats, species - 1), numpy.delete(nHats2, species - 1), rtol=InputTol)),"Error, densities don't seem to match!"
                assert (numpy.allclose(THats, THats2, rtol=InputTol)),"Error, temperatures don't seem to match!"
                assert (numpy.allclose(mHats, mHats2, rtol=InputTol)),"Error, masses don't seem to match!"
                assert (numpy.allclose(Zs, Zs2, rtol=InputTol)),"Error, charges don't seem to match!"            

                if includePhi1 == integerToRepresentTrue:
                    didNonlinearCalculationConverge = file1["didNonlinearCalculationConverge"][()]
                    didNonlinearCalculationConverge2 = file2["didNonlinearCalculationConverge"][()]
                    #if plotVariableName == "particleFlux_vm_rHat":
                    #    VariableValue = file1["particleFlux_vd_rHat"][()]

                    assert ((didNonlinearCalculationConverge == integerToRepresentTrue) and (didNonlinearCalculationConverge2 == integerToRepresentTrue2)),"Error, both runs must finish!"

                particleFluxName = "particleFlux_vd_rHat"
                classicalparticleFluxName = "classicalParticleFlux_rHat"
                DensityGradientName = "dnHatdrHat"                 
                
                if (particleFluxName.find('Flux_vd') != -1 or particleFluxName.find('Flux_vE') != -1) and (includePhi1 != integerToRepresentTrue):
                    #print (particleFluxName + " only exists in nonlinear runs, but this is a linear run.")
                    particleFluxName = particleFluxName.replace('Flux_vd', 'Flux_vm').replace('Flux_vE', 'Flux_vm')
                    #print ("Reading " + particleFluxName + " instead.")
                    #VariableValue = file1[plotVariableName.replace('Flux_vd', 'Flux_vm').replace('Flux_vE', 'Flux_vm')][()]
                #else:
                #    VariableValue = file1[plotVariableName][()]

                classicalParticleFlux = file1[classicalparticleFluxName][()]
                classicalParticleFlux = classicalParticleFlux[:, - 1]
                classicalParticleFlux = classicalParticleFlux[species - 1]

                classicalParticleFlux2 = file2[classicalparticleFluxName][()]
                classicalParticleFlux2 = classicalParticleFlux2[:, - 1]
                classicalParticleFlux2 = classicalParticleFlux2[species - 1]

                neoParticleFlux = file1[particleFluxName][()]
                neoParticleFlux = neoParticleFlux[:, - 1]
                neoParticleFlux = neoParticleFlux[species - 1]

                neoParticleFlux2 = file2[particleFluxName][()]
                neoParticleFlux2 = neoParticleFlux2[:, - 1]
                neoParticleFlux2 = neoParticleFlux2[species - 1]

                densityGradient = file1[DensityGradientName][()]
                densityGradient = densityGradient[species - 1]

                densityGradient2 = file2[DensityGradientName][()]
                densityGradient2 = densityGradient2[species - 1]
                nHat = nHats[species - 1]
                nHat2 = nHats2[species - 1]
                
                #classicalParticleFlux = TransformPlotVariableToOutputUnitsFactor * classicalParticleFlux
                #classicalParticleFlux = classicalParticleFlux / nHats[species - 1]

                #classicalParticleFlux2 = TransformPlotVariableToOutputUnitsFactor * classicalParticleFlux2
                #classicalParticleFlux2 = classicalParticleFlux2 / nHats2[species - 1] 

                file1.close()
                file2.close()
                
                
                #if plotVariableName.find('Flux_v') != -1:
                #    VariableValue = VariableValue[:, -1]
                #    VariableValue = VariableValue[species - 1] 

                #VariableValue = TransformPlotVariableToOutputUnitsFactor * VariableValue
                    
            except AssertionError as err:
                print ("Error when reading from " + fullSubDirectory + "/" + filename1 + " and " + filename2)
                print (err.args[0])
                print ("Continuing with next sub directory.")
                continue
            except:
                print ("Error when reading from " + fullSubDirectory + "/" + filename1 + " and " + filename2)
                print ("Maybe the SFINCS run did not finish")
                print ("Continuing with next sub directory.")
                continue

            #print("nHat: " + str(nHats[species -1]))

            #VariableValue = VariableValue / nHats[species -1]

            Dneo = - (neoParticleFlux2/nHat2 - neoParticleFlux/nHat) / (densityGradient2/nHat2 - densityGradient/nHat)
            Dclassical = - (classicalParticleFlux2/nHat2 - classicalParticleFlux/nHat) / (densityGradient2/nHat2 - densityGradient/nHat)
            
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

            LogarithmicDensityGradientData.append(densityGradient/nHat)

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
            LogarithmicDensityGradientData_sorted.append(LogarithmicDensityGradientData[radii.index(radius)])
        
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

        print ("")
        print ("LogarithmicDensityGradientData: " + str(LogarithmicDensityGradientData))
        print ("")
        print ("LogarithmicDensityGradientData_sorted: " + str(LogarithmicDensityGradientData_sorted))

        #print (np.array(radii_sorted))
        #print ("")
        #print (np.array(ydata_sorted))
        #print ("")
        #print (np.array(ydata2_sorted))
        
        try:
            LegendLabel = PlotLegendLabels[linenumber]
        except:
            LegendLabel = directory

        if SFINCSplotWithClassical[linenumber] :
            #D neoclassical + classical
            axD.plot(np.array(radii_sorted), SFINCSplotFactors[linenumber]*(np.array(DneoData_sorted) + np.array(DclassicalData_sorted)), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

            #V neoclassical + classical
            axV.plot(np.array(radii_sorted), SFINCSplotFactors[linenumber]*(np.array(VneoData_sorted) + np.array(VclassicalData_sorted)), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

            axDgrad.plot(np.array(radii_sorted), -np.array(LogarithmicDensityGradientData_sorted)*SFINCSplotFactors[linenumber]*(np.array(DneoData_sorted) + np.array(DclassicalData_sorted)), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
        
        else :
        
            #D neoclassical
            axD.plot(np.array(radii_sorted), SFINCSplotFactors[linenumber]*np.array(DneoData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

            #V neoclassical
            axV.plot(np.array(radii_sorted), SFINCSplotFactors[linenumber]*np.array(VneoData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

            axDgrad.plot(np.array(radii_sorted), -np.array(LogarithmicDensityGradientData_sorted)*SFINCSplotFactors[linenumber]*np.array(DneoData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

        #Plot Logarithmic Density Gradient
        axDensityGradient.plot(np.array(radii_sorted), -np.array(LogarithmicDensityGradientData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
        
        linenumber += 1

        ##External data must also have density gradient
        interpFunc = interpolate.interp1d(np.array(radii_sorted), np.array(LogarithmicDensityGradientData_sorted), fill_value="extrapolate")

        #try:
        #    LegendLabel = PlotLegendLabels[linenumber]
        #except:
        #    LegendLabel = directory
        #
        #D neoclassical + classical
        #axD.plot(np.array(radii_sorted), np.array(DneoData_sorted) + np.array(DclassicalData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
        #
        #V neoclassical + classical
        #axV.plot(np.array(radii_sorted), np.array(VneoData_sorted) + np.array(VclassicalData_sorted), PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
        #
        #linenumber += 1

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

    externalCounter = -1
    for externalfile in os.listdir(originalDirectory):
        if externalfile.endswith(externalDataFileType):
            externalCounter += 1
            try:
                inputParams = np.genfromtxt(externalfile, dtype=None, comments="#", skip_header=1)
                #print(inputParams[:,radiusColumn])
                #print(inputParams[:,ErColumn])
                try:
                    LegendLabel = PlotLegendLabels[linenumber]
                except:
                    LegendLabel = externalfile

                externalDtoPlot = inputParams[:,Dcolumn[externalCounter]]
                if (withClassicalExternal[externalCounter] and classicalDcolumn[externalCounter] != -1):
                    externalDtoPlot = externalDtoPlot + inputParams[:,classicalDcolumn[externalCounter]]

                externalVtoPlot = inputParams[:,Vcolumn[externalCounter]]
                if (withClassicalExternal[externalCounter] and classicalVcolumn[externalCounter] != -1):
                    externalVtoPlot = externalVtoPlot + inputParams[:,classicalVcolumn[externalCounter]]

                externalSigmaDtoPlot = np.array([])
                if (ErrorBars[linenumber] and (sigmaDColumn[externalCounter] != -1 or sigmaClassicalDcolumn[externalCounter] != -1)):
                    if sigmaDColumn[externalCounter] == -1:
                        externalSigmaDtoPlot = inputParams[:,sigmaClassicalDcolumn[externalCounter]]
                    elif sigmaClassicalDcolumn[externalCounter] == -1:
                        externalSigmaDtoPlot = inputParams[:,sigmaDColumn[externalCounter]]
                    else :
                        externalSigmaDtoPlot = inputParams[:,sigmaDColumn[externalCounter]] + inputParams[:,sigmaClassicalDcolumn[externalCounter]]

                externalSigmaVtoPlot = np.array([])
                if (ErrorBars[linenumber] and (sigmaVColumn[externalCounter] != -1 or sigmaClassicalVcolumn[externalCounter] != -1)):
                    if sigmaVColumn[externalCounter] == -1:
                        externalSigmaVtoPlot  = inputParams[:,sigmaClassicalVcolumn[externalCounter]]
                    elif sigmaClassicalVcolumn[externalCounter] == -1:
                        externalSigmaVtoPlot = inputParams[:,sigmaVColumn[externalCounter]]
                    else :
                        externalSigmaVtoPlot = inputParams[:,sigmaVColumn[externalCounter]] + inputParams[:,sigmaClassicalVcolumn[externalCounter]]

                externalLogarithmicDensityGradient = np.array(interpFunc(inputParams[:,radiusColumn[externalCounter]]/aNorm))
                #print(str(externalLogarithmicDensityGradient))
                
                        
                if (FilledErrors or (not ErrorBars[linenumber]) or (sigmaDColumn[externalCounter] == -1 and sigmaClassicalDcolumn[externalCounter] == -1)) :

                    #D external
                    axD.plot(inputParams[:,radiusColumn[externalCounter]]/aNorm, externalDtoPlot, PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

                    axDgrad.plot(inputParams[:,radiusColumn[externalCounter]]/aNorm, -externalLogarithmicDensityGradient*externalDtoPlot, PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
                    
                    if (ErrorBars[linenumber] and externalSigmaDtoPlot.size != 0):
                        axD.fill_between(inputParams[:,radiusColumn[externalCounter]]/aNorm, externalDtoPlot + externalSigmaDtoPlot, externalDtoPlot - externalSigmaDtoPlot, color=PlotLineColors[linenumber], alpha=ErrorBarAlpha)
                        axDgrad.fill_between(inputParams[:,radiusColumn[externalCounter]]/aNorm, -externalLogarithmicDensityGradient*(externalDtoPlot + externalSigmaDtoPlot), -externalLogarithmicDensityGradient*(externalDtoPlot - externalSigmaDtoPlot), color=PlotLineColors[linenumber], alpha=ErrorBarAlpha)
                else :
                    axD.errorbar(inputParams[:,radiusColumn[externalCounter]]/aNorm, externalDtoPlot, fmt=PlotLinespecs[linenumber], yerr=externalSigmaDtoPlot, color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

                    axDgrad.errorbar(inputParams[:,radiusColumn[externalCounter]]/aNorm, -externalLogarithmicDensityGradient*externalDtoPlot, fmt=PlotLinespecs[linenumber], yerr=np.absolute(externalLogarithmicDensityGradient*externalSigmaDtoPlot), color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)


                if (FilledErrors or (not ErrorBars[linenumber]) or (sigmaVColumn[externalCounter] == -1 and sigmaClassicalVcolumn[externalCounter] == -1)) :
                    
                    #V external
                    axV.plot(inputParams[:,radiusColumn[externalCounter]]/aNorm, externalVtoPlot, PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)
                    
                    if (ErrorBars[linenumber] and externalSigmaVtoPlot.size != 0):
                        axV.fill_between(inputParams[:,radiusColumn[externalCounter]]/aNorm, externalVtoPlot + externalSigmaVtoPlot, externalVtoPlot - externalSigmaVtoPlot, color=PlotLineColors[linenumber], alpha=ErrorBarAlpha)

                else :
                    axV.errorbar(inputParams[:,radiusColumn[externalCounter]]/aNorm, externalVtoPlot, fmt=PlotLinespecs[linenumber], yerr=externalSigmaVtoPlot, color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

                axDensityGradient.plot(inputParams[:,radiusColumn[externalCounter]]/aNorm, - externalLogarithmicDensityGradient, PlotLinespecs[linenumber], color=PlotLineColors[linenumber], markersize=PlotMarkerSize, markeredgewidth=PlotMarkerEdgeWidth[linenumber], markeredgecolor=PlotLineColors[linenumber], label=LegendLabel, linewidth=PlotLineWidth)

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
axDensityGradient.set_xscale(xAxisScale)
axDensityGradient.set_yscale(yAxisScale)
axDgrad.set_xscale(xAxisScale)
axDgrad.set_yscale(yAxisScale)

if ShowGrid:
    #plt.grid(which='both')
    axD.grid(which='both')
    axV.grid(which='both')
    axDensityGradient.grid(which='both')
    axDgrad.grid(which='both')

#plt.xlabel(xAxisLabel, fontsize=AxesLabelSize)
#plt.ylabel(yAxisLabel, fontsize=AxesLabelSize)

axD.set_xlabel(xAxisLabel, fontsize=AxesLabelSize)
axD.set_ylabel(yAxisLabel[0], fontsize=AxesLabelSize)
axV.set_xlabel(xAxisLabel, fontsize=AxesLabelSize)
axV.set_ylabel(yAxisLabel[1], fontsize=AxesLabelSize)
axDensityGradient.set_xlabel(xAxisLabel, fontsize=AxesLabelSize)
axDensityGradient.set_ylabel(yAxisLabel[2], fontsize=AxesLabelSize)
axDgrad.set_xlabel(xAxisLabel, fontsize=AxesLabelSize)
axDgrad.set_ylabel(yAxisLabel[3], fontsize=AxesLabelSize)

if AxisLimAuto: 
    #ymin,ymax = plt.ylim()
    #xmin,xmax = plt.xlim()
    yminD,ymaxD = axD.get_ylim()
    xminD,xmaxD = axD.get_xlim()
    yminV,ymaxV = axV.get_ylim()
    xminV,xmaxV = axV.get_xlim()
    DensityGradientymin,DensityGradientymax = axDensityGradient.get_ylim()
    DensityGradientxmin,DensityGradientxmax = axDensityGradient.get_xlim()
    Dgradymin,Dgradymax = axDgrad.get_ylim()
    Dgradxmin,Dgradxmax = axDgrad.get_xlim()
else : 
    #ymin,ymax = plt.ylim(yAxisLim) 
    #xmin,xmax = plt.xlim(xAxisLim)
    yminD,ymaxD = axD.set_ylim(yAxisLim[0]) 
    xminD,xmaxD = axD.set_xlim(xAxisLim[0])
    yminV,ymaxV = axV.set_ylim(yAxisLim[1]) 
    xminV,xmaxV = axV.set_xlim(xAxisLim[1])
    DensityGradientymin,DensityGradientymax = axDensityGradient.set_ylim(yAxisLim[2]) 
    DensityGradientxmin,DensityGradientxmax = axDensityGradient.set_xlim(xAxisLim[2])
    Dgradymin,Dgradymax = axDgrad.set_ylim(yAxisLim[3])
    Dgradxmin,Dgradxmax = axDgrad.set_xlim(xAxisLim[3])

if ShowLegend:
    #plt.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)#, fontsize=LegendFontSize)
    if legendAtDplot:
        axD.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)#, fontsize=LegendFontSize)
    else :
        axV.legend(bbox_to_anchor = LegendBBoxToAnchor, loc=LegendPosition, ncol=LegendNumberColumns, mode=None, borderaxespad=0., prop=LegendProperties)

#plt.gca().axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
#plt.gca().axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])

axD.axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
axD.axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])
axV.axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
axV.axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])
axDensityGradient.axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
axDensityGradient.axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])
axDgrad.axes.xaxis.set_label_coords(xAxisLabelCoords[0], xAxisLabelCoords[1])
axDgrad.axes.yaxis.set_label_coords(yAxisLabelCoords[0], yAxisLabelCoords[1])

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
    axDensityGradient.text(SubPlotLabelXcoord[2], SubPlotLabelYcoord[2], SubPlotLabel[2])
    axDgrad.text(SubPlotLabelXcoord[3], SubPlotLabelYcoord[3], SubPlotLabel[3])

if NoScientificAxes :
    try:
        #ax.get_xaxis().get_major_formatter().set_scientific(False)
        #ax.get_yaxis().get_major_formatter().set_scientific(False)
        axD.get_xaxis().get_major_formatter().set_scientific(False)
        axD.get_yaxis().get_major_formatter().set_scientific(False)
        axV.get_xaxis().get_major_formatter().set_scientific(False)
        axV.get_yaxis().get_major_formatter().set_scientific(False)
        axDensityGradient.get_xaxis().get_major_formatter().set_scientific(False)
        axDensityGradient.get_yaxis().get_major_formatter().set_scientific(False)
        axDgrad.get_xaxis().get_major_formatter().set_scientific(False)
        axDgrad.get_yaxis().get_major_formatter().set_scientific(False)
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
