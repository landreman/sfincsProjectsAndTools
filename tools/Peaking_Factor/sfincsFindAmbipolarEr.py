#!/usr/bin/env python
# USAGE: 
#
#    python sfincsFindAmbipolarEr.py <InputFile>

import sys
import os
import shutil
import subprocess
import numpy as np
from tempfile import mkstemp
import h5py
import matplotlib as mpl
mpl.use('PDF') #Use backend to plot in pdf file
import matplotlib.pyplot as plt

#----------------------------------------
# Capture command-line arguments
#

directory = os.getcwd() #Current working directory
sfincsInput = "input.namelist"
sfincsOutput = "sfincsOutput.h5"

plotOutput = "Plot_sfincsFindAmbipolarEr.pdf"
arrayOutput = "sfincsFindAmbipolarEr_Array.out"
resultOutput = "sfincsFindAmbipolarEr_Result.out"

if len(sys.argv) > 2:
	print "Too many input arguments for %s" % sys.argv[0]
	sys.exit(0)
elif len(sys.argv) == 2: 
	inputfile = sys.argv[1]
else:
	inputfile = "input.sfincsFindAmbipolarEr"

if not (os.path.exists(directory + "/" + inputfile) or os.path.exists(inputfile)):
	print "Input file %s not found" % inputfile
	sys.exit(0)


#To solve python rounding issue, i.e. floats are not represented by exact values (important for comparisons).
#MinFloat = sys.float_info.epsilon
MinFloat = pow(10, -sys.float_info.dig)

#Go to directory
os.chdir(directory)

print "**********************************************************"
print "Starting sfincsFindAmbipolarEr.py in directory %s" % os.getcwd()
print "**********************************************************"

#-------------------------------------------------------
# Read input parameters from file
print "Reading input from %s" % inputfile

#Standard values
NumVariables = 0
ABS_ERROR = 1e-8
REL_ERROR = 1e0
MAX_ITER = 10
NPROC = 1
CONTINUE_SCAN = 0
CONTINUE_DIR = 1
dPhiHatdpsiN_START = 0
PATH_TO_SFINCS = "sfincs"
USE_APRUN = False


try: #Try to read from inputfile
	inputParams = np.genfromtxt(inputfile, dtype='str', comments="#")
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while reading %s" % inputfile
	print "\nUSING STANDARD INPUT VALUES"
else:
	VariableName = inputParams[:,0]
	#VariableValue = inputParams[:,1].astype(np.float)
	VariableValue = inputParams[:,1]
	NumVariables = len(VariableName)

i=0
while i < NumVariables:
	if VariableName[i] == 'ABS_ERROR':
		ABS_ERROR = VariableValue[i].astype(np.float)
	elif VariableName[i] == 'REL_ERROR':
		REL_ERROR = VariableValue[i].astype(np.float)
	elif VariableName[i] == 'MAX_ITER':
		MAX_ITER = int(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'NPROC':
		NPROC = int(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'CONTINUE_SCAN':
		CONTINUE_SCAN = bool(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'CONTINUE_DIR':
		CONTINUE_DIR = int(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'dPhiHatdpsiN_START':
		dPhiHatdpsiN_START = VariableValue[i].astype(np.float)
	elif VariableName[i] == 'SFINCS_INPUT':
		sfincsInput = VariableValue[i]
	elif VariableName[i] == 'SFINCS_OUTPUT':
		sfincsOutput = VariableValue[i]
	elif VariableName[i] == 'PLOT_OUTPUT':
		plotOutput = VariableValue[i]
	elif VariableName[i] == 'ARRAY_OUTPUT':
		arrayOutput = VariableValue[i]
	elif VariableName[i] == 'RESULT_OUTPUT':
		resultOutput = VariableValue[i]
	elif VariableName[i] == 'USE_APRUN':
		USE_APRUN = bool(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'PATH_TO_SFINCS':
		PATH_TO_SFINCS = VariableValue[i]
	else:
		print "Got false variable: %s" % VariableName[i]
	i += 1

print "Using parameters:"
print "ABS_ERROR: %E" % ABS_ERROR
print "REL_ERROR: %E" % REL_ERROR
print "MAX_ITER: %i" % MAX_ITER
print "NPROC: %i" % NPROC
print "CONTINUE_SCAN: %i" % CONTINUE_SCAN
print "CONTINUE_DIR: %i" % CONTINUE_DIR
print "PATH_TO_SFINCS: %s" % PATH_TO_SFINCS
print "USE_APRUN: %i" % USE_APRUN
print "dPhiHatdpsiN_START: %E" % dPhiHatdpsiN_START
print "SFINCS_INPUT: %s" % sfincsInput
print "SFINCS_OUTPUT: %s" % sfincsOutput
print "PLOT_OUTPUT: %s" % plotOutput
print "ARRAY_OUTPUT: %s" % arrayOutput
print "RESULT_OUTPUT: %s" % resultOutput
print "**********************************************************"

##SFINCS INPUT/OUTPUT SECTION

from sfincsInputOutput import *


###########################################################################
##Start execution of program
###########################################################################
dirIter = 1
dPhiHatdpsiN = dPhiHatdpsiN_START
#previousdPhiHatdpsiN = dPhiHatdpsiN_START + 1
previousdPhiHatdpsiN = 2*dPhiHatdpsiN + 10*MinFloat*(1 + 2*np.sign(dPhiHatdpsiN))
sumCurrent = max(2*abs(ABS_ERROR), MinFloat)

#Go to directory
os.chdir(directory)

if readInput("programMode", 0, sfincsInput) != 1:
	print "sfincsFindAmbipolarEr.py has to be run with programMode = 1 in SFINCS"	
	sys.exit(0)

if CONTINUE_SCAN:
	dirIter = CONTINUE_DIR
	MAX_ITER = MAX_ITER + CONTINUE_DIR - 1 #Need to increase MAX_ITER so the effective number of iterations are correct
	if not os.path.exists(str(dirIter - 1)):
		print "Error: Can not find directory %s in current directory %s" % (str(dirIter - 1), os.getcwd())
		sys.exit(0)
	try:
		with h5py.File(str(dirIter - 1) + "/" + sfincsOutput, 'r') as sfincsOutputFile:
			sumCurrent = (np.array(sfincsOutputFile['run  1/particleFlux'].value)).dot(np.array(sfincsOutputFile['run  1/Zs'].value))
			dPhiHatdpsiN = sfincsOutputFile['run  1/d(PhiHat)d(psi_N)'].value
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while trying to read %s in %s, check if file exists or file access rights" % ((str(dirIter - 1) + "/" + sfincsOutput), os.getcwd())
		sys.exit(0)
	
	if dirIter < 3: #Only one iteration done earlier
		previousSumCurrent = max(2*abs(ABS_ERROR), MinFloat) 
		#previousdPhiHatdpsiN = dPhiHatdpsiN_START + 1
		previousdPhiHatdpsiN = 2*dPhiHatdpsiN + 10*MinFloat*(1 + 2*np.sign(dPhiHatdpsiN))
	else:
		try:
			with h5py.File(str(dirIter - 2) + "/" + sfincsOutput, 'r') as sfincsOutputFile:
				previousSumCurrent = (np.array(sfincsOutputFile['run  1/particleFlux'].value)).dot(np.array(sfincsOutputFile['run  1/Zs'].value))
				previousdPhiHatdpsiN = sfincsOutputFile['run  1/d(PhiHat)d(psi_N)'].value
		except Exception,e:
			print e.__class__.__name__, ": ", e.message, "while trying to read %s in %s, check if file exists or file access rights" % ((str(dirIter - 2) + "/" + sfincsOutput), os.getcwd())
			sys.exit(0)

	if previousSumCurrent == sumCurrent or previousdPhiHatdpsiN == dPhiHatdpsiN:
		print "Error: Iteration scheme can not update dPhiHatdpsiN"
		sys.exit(0)
	#Update radial electric field
	#From the assumption of a linear dependence
	temporary_dPhiHatdpsiN = previousdPhiHatdpsiN - (dPhiHatdpsiN - previousdPhiHatdpsiN) / (sumCurrent - previousSumCurrent) * previousSumCurrent
	previousdPhiHatdpsiN = dPhiHatdpsiN
	dPhiHatdpsiN = temporary_dPhiHatdpsiN	
	

relative_error = abs(sumCurrent) + abs(REL_ERROR) #To keep track of relative error between iterations step

while dirIter <= MAX_ITER and (abs(ABS_ERROR) <= abs(sumCurrent) or abs(REL_ERROR) <= abs(relative_error)):
	print "#######################################################################################"
	print "Find Ambipolar Electric Field Iteration %i" %dirIter
	print "Radial electric field %f" %dPhiHatdpsiN

	if os.path.exists(str(dirIter)):
		try:
			shutil.rmtree(str(dirIter))
		except Exception,e:
			print e.__class__.__name__, ": ", e.message, "while trying to remove %s in %s, check access rights" % (str(dirIter), os.getcwd())
			sys.exit(0)
	os.mkdir(str(dirIter))

	
	try:
		shutil.copyfile(sfincsInput, str(dirIter) + "/" +  sfincsInput)
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while trying to copy %s into directory %s, check access rights" % (sfincsInput, os.getcwd() + str(dirIter))
		sys.exit(0)
	
	writeInput("dPhiHatdpsiN", dPhiHatdpsiN, str(dirIter) + "/" +  sfincsInput)
	#Run SFINCS
	try:
		os.chdir(str(dirIter))
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while trying to enter %s from %s" % (str(dirIter), os.getcwd() )
		sys.exit(0)

	print "Entering %s and running SFINCS" % os.getcwd() #+ "/" + str(dirIter)


	try:
		#subprocess.check_call(["mpirun", "-n", "%i" %NPROC, "sfincs"])
                if USE_APRUN:
                        procSFINCS = subprocess.Popen(["aprun", "-n", "%i" %NPROC, PATH_TO_SFINCS], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                else:
                        procSFINCS = subprocess.Popen(["mpirun", "-n", "%i" %NPROC, PATH_TO_SFINCS], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#print (procSFINCS.communicate())[0]
		procSFINCSOutput = procSFINCS.communicate()
	      	print procSFINCSOutput[0]
		print procSFINCSOutput[1]
	except subprocess.CalledProcessError,e:
		print e.__class__.__name__, ": ", e.message, "SFINCS stopped with exit code %i" % e.returncode
		sys.exit(0)
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "when trying to run SFINCS"
		sys.exit(0)
			
	#Return to base directory
	os.chdir(directory)

	#Read output and update error
	previousSumCurrent = sumCurrent #To keep track of relative error between iterations step

	try:
		sfincsOutputFile = h5py.File(str(dirIter) + "/" + sfincsOutput, 'r')
		sumCurrent = (np.array(sfincsOutputFile['run  1/particleFlux'].value)).dot(np.array(sfincsOutputFile['run  1/Zs'].value))
		sfincsOutputFile.close()
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "when trying to read SFINCS output from directory %s" % os.getcwd() + "/" + str(dirIter)
		sys.exit(0)


	relative_error = abs((sumCurrent - previousSumCurrent) / previousSumCurrent)
	print "Current absolute error: %e" % abs(sumCurrent)
	print "Current relative error: %e" % relative_error

	if abs(sumCurrent) == 0.0 or abs(relative_error) == 0.0: # Solution found or can not continue iteration
		dirIter += 1
		break

	#Update radial electric field
	#From the assumption of a linear dependence
	temporary_dPhiHatdpsiN = previousdPhiHatdpsiN - (dPhiHatdpsiN - previousdPhiHatdpsiN) / (sumCurrent - previousSumCurrent) * previousSumCurrent
	previousdPhiHatdpsiN = dPhiHatdpsiN
	dPhiHatdpsiN = temporary_dPhiHatdpsiN

	dirIter += 1
#End while loop

###########################################################################
##Write results of the scan to file
###########################################################################

print "#######################################################################################"
print "Scan has finished, writing output to files"
print "#######################################################################################"

NumDir = dirIter - 1

#Go to directory
os.chdir(directory)

List_sumCurrent = []
List_dPhiHatdpsiN = []
List_dirNo = []

Result_fluxes = []
Result_Zs = []
Result_dPhiHatdpsiN = 0

for dirIter in xrange(1, NumDir + 1):

	if not os.path.exists(str(dirIter)):
		print "WARNING: Could not find directory %s in %s" % (str(dirIter), os.getcwd())
		continue
	
	List_dirNo.append(dirIter)
	try:
		sfincsOutputFile = h5py.File(str(dirIter) + "/" + sfincsOutput, 'r')

		Result_fluxes = np.array(sfincsOutputFile['run  1/particleFlux'].value)
		Result_Zs = np.array(sfincsOutputFile['run  1/Zs'].value)
		Result_dPhiHatdpsiN = sfincsOutputFile['run  1/d(PhiHat)d(psi_N)'].value

		List_sumCurrent.append(Result_fluxes.dot(Result_Zs) )
		List_dPhiHatdpsiN.append( Result_dPhiHatdpsiN )		

		sfincsOutputFile.close()
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "when trying to read SFINCS output from directory %s" % os.getcwd() + "/" + str(dirIter)
		sys.exit(0)


#End for loop

try:
	with open(arrayOutput, 'w') as f:
		f.write("\nIteration:\t")
		f.write(str(List_dirNo))
		f.write("\nCurrent:\t")
		f.write(str(List_sumCurrent))
		f.write("\ndPhiHatdpsiN:\t")
		f.write(str(List_dPhiHatdpsiN))
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while trying to open %s, check access rights" % arrayOutput
	sys.exit(0)
try:
	with open(resultOutput, 'w') as f:
		f.write("\ndPhiHatdpsiN*:\t")
		f.write(str(Result_dPhiHatdpsiN))
		f.write("\nZs:\t")
		f.write(str(Result_Zs))
		f.write("\nFluxes:\t")
		f.write(str(Result_fluxes))
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while trying to open %s, check access rights" % resultOutput
	sys.exit(0)

###########################################################################
##Plot results of the scan
###########################################################################
if len(List_dirNo) < 1: #Nothing to plot
	sys.exit(0)

plt.ioff() #Turn interactive mode off in matplotlib

#Set up plotting options:
font = {'family' : 'Sans Serif',
        'weight' : 'normal',
        'size'   : 30}

mpl.rc('font', **font)

mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['lines.linewidth'] = 4

#End plotting options

Array_dirNo = np.array(List_dirNo)
Array_sumCurrent = np.array(List_sumCurrent)
Array_dPhiHatdpsiN = np.array(List_dPhiHatdpsiN)

print "Directories to plot: %s" % str(Array_dirNo)
print "Radial currents: %s" % str(Array_sumCurrent)
print "Radial electric field: %s" % str(Array_dPhiHatdpsiN)

fig = plt.figure(figsize=(25,40))

##########
##Plot 1
##########

ax1 = plt.subplot(3, 1, 1)
plt.plot(Array_dirNo, Array_sumCurrent, 'mx-', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Radial current per iteration step', weight='bold')
plt.ylabel(r'$\sum_{\alpha} e_{\alpha} \Gamma_{\alpha}$', fontsize=40)
plt.xlabel(r'Iteration', fontsize=30)
currMin = np.amin(np.absolute(Array_sumCurrent))
currMax = np.amax(np.absolute(Array_sumCurrent))
space_lin1 = 1
linlog_limit1 = MinFloat
if 0 < currMax and 0 < currMin:
	space_lin1 = np.maximum(np.ceil(np.log10(currMax / currMin)/4.0), 1.0)
	linlog_limit1 = np.power(10, np.ceil(np.log10(currMin)))
linlog_limit1_array = np.empty(Array_dirNo.size); linlog_limit1_array.fill(linlog_limit1)
	
plt.yscale('symlog', linthreshy = linlog_limit1, linscaley = space_lin1)
plt.fill_between(Array_dirNo, linlog_limit1_array, -linlog_limit1_array, facecolor='green', alpha=0.08)
plt.text(0.10, 0.5, r'Linear scale', color='green', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)

##########
##Plot 2
##########

ax2 = plt.subplot(3, 1, 2)
Array_dPhiHatdpsiN_diff = Array_dPhiHatdpsiN - Array_dPhiHatdpsiN[-1]
plt.plot(Array_dirNo, Array_dPhiHatdpsiN_diff , 'rx-', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Radial electric field per iteration step', weight='bold')
plt.xlabel(r'Iteration', fontsize=30)
plt.ylabel(r'$d \hat{\Phi} / d \psi_N - \left(d \hat{\Phi} / d \psi_N\right)^\ast$', fontsize=40)

potMin = 0; potMax = 0
if Array_dPhiHatdpsiN_diff.size > 1:
	potMin = np.amin(np.absolute(Array_dPhiHatdpsiN_diff[0:-1]))
	potMax = np.amax(np.absolute(Array_dPhiHatdpsiN_diff[0:-1]))
space_lin2 = 1
linlog_limit2 = MinFloat
if 0 < potMax and 0 < potMin:
	space_lin2 = np.maximum(np.ceil(np.log10(potMax / potMin)/4.0), 1.0)
	linlog_limit2 = np.power(10, np.ceil(np.log10(potMin)))
linlog_limit2_array = np.empty(Array_dirNo.size); linlog_limit2_array.fill(linlog_limit2)

plt.yscale('symlog', linthreshy = linlog_limit2, linscaley = space_lin2)
plt.fill_between(Array_dirNo, linlog_limit2_array, -linlog_limit2_array, facecolor='green', alpha=0.08)
plt.text(0.10, 0.50, r'Linear scale', color='green', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes)
plt.text(0.80, 0.90, r'$\left(d \hat{\Phi} / d \psi_N\right)^\ast$ = %e' %Array_dPhiHatdpsiN[-1], color='black', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes)

##########
##Plot 3
##########

ax3 = plt.subplot(3, 1, 3)
plt.plot(Array_dPhiHatdpsiN_diff, Array_sumCurrent, 'cx', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Radial current as function of radial electric field', weight='bold')
plt.xlabel(r'$d \hat{\Phi} / d \psi_N - \left(d \hat{\Phi} / d \psi_N\right)^\ast$', fontsize=40)
plt.ylabel(r'$\sum_{\alpha} e_{\alpha} \Gamma_{\alpha}$', fontsize=40)
p2, = plt.plot(Array_dPhiHatdpsiN_diff[-1], Array_sumCurrent[-1], color='red', linestyle='none', marker='o', markersize=15, markeredgewidth=0)
plt.legend([p2], [r'Ambipolar $d \hat{\Phi} / d \psi_N$'], loc=0, numpoints=1)

plt.yscale('symlog', linthreshy = linlog_limit1, linscaley = space_lin1)
x_array_lim = np.empty(Array_dPhiHatdpsiN_diff.size); x_array_lim[0] = - linlog_limit2; x_array_lim[-1] = linlog_limit2;

plt.fill_between(x_array_lim, np.amin(Array_sumCurrent), np.amax(Array_sumCurrent), facecolor='green', alpha=0.08)
plt.fill_between(np.sort(Array_dPhiHatdpsiN_diff), linlog_limit1_array, -linlog_limit1_array, facecolor='green', alpha=0.08)

plt.fill_between(x_array_lim, linlog_limit1_array, -linlog_limit1_array, facecolor='green', alpha=0.16)
plt.text(0.10, 0.5, r'Linear scale', color='green', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)
plt.xscale('symlog', linthreshx = linlog_limit2, linscalex = space_lin2)
plt.text(0.80, 0.90, r'$\left(d \hat{\Phi} / d \psi_N\right)^\ast$ = %e' %Array_dPhiHatdpsiN[-1], color='black', horizontalalignment='center',verticalalignment='center', transform=ax3.transAxes)

print "\nWriting results of scan to %s" % plotOutput
plt.savefig(plotOutput)


sys.exit(0)


