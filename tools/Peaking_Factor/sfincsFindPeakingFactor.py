#!/usr/bin/env python

# USAGE: 
#
#    python sfincsFindPeakingFactor.py <InputFile>

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
sfincsFindAmbipolarEr_resultOutput = "sfincsFindAmbipolarEr_Result.out"

plotOutput_peaking = "Plot_sfincsFindPeakingFactor.pdf"
arrayOutput_peaking = "sfincsFindPeakingFactor_Array.out"
resultOutput_peaking = "sfincsFindPeakingFactor_Result.out"

if len(sys.argv) > 2:
	print "Too many input arguments for %s" % sys.argv[0]
	sys.exit(0)
elif len(sys.argv) == 2: 
	inputfile = sys.argv[1]
else:
	inputfile = "input.sfincsFindPeakingFactor"

if not (os.path.exists(directory + "/" + inputfile) or os.path.exists(inputfile)):
	print "Input file %s not found" % inputfile
	sys.exit(0)

#To solve python rounding issue, i.e. floats are not represented by exact values (important for comparisons).
#MinFloat = sys.float_info.epsilon
MinFloat = pow(10, -sys.float_info.dig)

#Go to directory
os.chdir(directory)

print "**********************************************************"
print "Starting sfincsFindPeakingFactor.py in directory %s" % os.getcwd()
print "**********************************************************"
#print "\n"

#-------------------------------------------------------
# Read input parameters from file
print "Reading input from %s" % inputfile


#Standard values
NumVariables = 0
sfincsFindAmbipolarEr_INPUT = "input.sfincsFindAmbipolarEr"
ABS_ERROR_PEAKING = 1e-4
REL_ERROR_PEAKING = 1e0
MAX_ITER_PEAKING = 10
CONTINUE_SCAN_PEAKING = False
CONTINUE_DIR_PEAKING = 1
dNHatdpsiN_START_PEAKING = -1.0
CHECK_QUASINEUTRALITY_PEAKING = True

try: #Try to read from inputfile
	inputParams = np.genfromtxt(inputfile, dtype='str', comments="#")
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while reading %s" % inputfile
	print "\nUSING STANDARD INPUT VALUES"
else:
	VariableName = inputParams[:,0]
	VariableValue = inputParams[:,1]
	NumVariables = len(VariableName)

i=0
while i < NumVariables:
	if VariableName[i] == 'ABS_ERROR_PEAKING':
		ABS_ERROR_PEAKING = VariableValue[i].astype(np.float)
	elif VariableName[i] == 'REL_ERROR_PEAKING':
		REL_ERROR_PEAKING = VariableValue[i].astype(np.float)
	elif VariableName[i] == 'MAX_ITER_PEAKING':
		MAX_ITER_PEAKING = int(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'CONTINUE_SCAN_PEAKING':
		CONTINUE_SCAN_PEAKING = bool(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'CONTINUE_DIR_PEAKING':
		CONTINUE_DIR_PEAKING = int(VariableValue[i].astype(np.float))
	elif VariableName[i] == 'dNHatdpsiN_START_PEAKING':
		dNHatdpsiN_START_PEAKING = VariableValue[i].astype(np.float)
	elif VariableName[i] == 'PLOT_OUTPUT_PEAKING':
		plotOutput_peaking = VariableValue[i]
	elif VariableName[i] == 'ARRAY_OUTPUT_PEAKING':
		arrayOutput_peaking = VariableValue[i]
	elif VariableName[i] == 'RESULT_OUTPUT_PEAKING':
		resultOutput_peaking = VariableValue[i]
	elif VariableName[i] == 'sfincsFindAmbipolarEr_INPUT':
		sfincsFindAmbipolarEr_INPUT = VariableValue[i]
	elif VariableName[i] == 'CHECK_QUASINEUTRALITY_PEAKING':
		CHECK_QUASINEUTRALITY_PEAKING = bool(VariableValue[i].astype(np.float))
	else:
		print "Got false variable: %s" % VariableName[i]
	i += 1

#Check if input file to FindAmbipolarEr exists, and then find the SFINCS input file from it
if not (os.path.exists(directory + "/" + sfincsFindAmbipolarEr_INPUT) or os.path.exists(sfincsFindAmbipolarEr_INPUT)):
	print "Input file %s to sfincsFindAmbipolarEr.py not found" % sfincsFindAmbipolarEr_INPUT
	sys.exit(0)
else: 
	try: #Try to read from inputfile
		inputParamsEr = np.genfromtxt(sfincsFindAmbipolarEr_INPUT, dtype='str', comments="#")
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while reading %s" % sfincsFindAmbipolarEr_INPUT
		sys.exit(0)
	else:
		VariableNameEr = inputParamsEr[:,0]
		VariableValueEr = inputParamsEr[:,1]
		NumVariablesEr = len(VariableNameEr)
		j = 0
		while j < NumVariablesEr:
			if VariableNameEr[j] == 'SFINCS_INPUT':
				sfincsInput = VariableValueEr[j]
			elif VariableNameEr[j] == 'RESULT_OUTPUT':
				sfincsFindAmbipolarEr_resultOutput = VariableValueEr[j]
			j += 1

print "Using parameters:"
print "ABS_ERROR_PEAKING: %E" % ABS_ERROR_PEAKING
print "REL_ERROR_PEAKING: %E" % REL_ERROR_PEAKING
print "MAX_ITER_PEAKING: %i" % MAX_ITER_PEAKING
print "CONTINUE_SCAN_PEAKING: %i" % CONTINUE_SCAN_PEAKING
print "CONTINUE_DIR_PEAKING: %i" % CONTINUE_DIR_PEAKING
print "dNHatdpsiN_START_PEAKING: %E" % dNHatdpsiN_START_PEAKING
print "PLOT_OUTPUT_PEAKING: %s" % plotOutput_peaking
print "ARRAY_OUTPUT_PEAKING: %s" % arrayOutput_peaking
print "RESULT_OUTPUT_PEAKING: %s" % resultOutput_peaking
print "sfincsFindAmbipolarEr_INPUT: %s" % sfincsFindAmbipolarEr_INPUT
print "SFINCS_INPUT: %s" % sfincsInput
print "sfincsFindAmbipolarEr_resultOutput: %s" % sfincsFindAmbipolarEr_resultOutput
print "CHECK_QUASINEUTRALITY_PEAKING: %i" % CHECK_QUASINEUTRALITY_PEAKING

print "**********************************************************"

from sfincsInputOutput import *

###########################################################################
##Define functions
###########################################################################

def readsfincsFindAmbipolarEr_Result(variableName, ResultFile=sfincsFindAmbipolarEr_resultOutput):
	"Reads an array of floats from the result output file of sfincsFindAmbipolarEr.py"
	if not os.path.exists(ResultFile):
		print "Error: Can not find output file %s from sfincsFindAmbipolarEr.py in directory %s" % (ResultFile, os.getcwd())
		sys.exit(0)
	variableName = variableName.strip()
	if variableName[-1] == ':':
		variableName = variableName[0:-1]
		variableName = variableName.strip()

	if not (variableName == 'dPhiHatdpsiN*' or variableName == 'Zs' or variableName == 'Fluxes'):
		print "WARNING: %s is not a valid output parameter from the result output file of sfincsFindAmbipolarEr.py, returning False" % variableName
		return False
	try:
		with open(ResultFile, 'r') as f:
			lines = [line.strip() for line in f]
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while trying to open %s, check access rights" % ResultFile
		sys.exit(0)

	numMatches = 0
	value = 0
	for line in lines:
		if len(line) == 0:
			continue
		elif line.find(variableName) == -1:
			continue
		elif line.split()[0] != variableName + ":":
			continue
		else:
			try:
				value = ' '.join(line.split()[1:]).replace('[', ' ').replace(']', ' ')
				value = [float(i) for i in value.split()]
				value = np.array(value)
			except Exception,e:
				print e.__class__.__name__, ": ", e.message, "while trying to read %s from %s" % (variableName, ResultFile)
				continue
			numMatches += 1

	if numMatches != 1:
		print "WARNING: %s seems to be corrupt, returning False" % ResultFile
		return False

	print "Read %s: %s" % (variableName, str(value))
	return value

###########################################################################
##Start execution of program
###########################################################################
dirIter = 1
dNHatdpsiN = dNHatdpsiN_START_PEAKING
#previousdNHatdpsiN = dNHatdpsiN + 1
previousdNHatdpsiN = 2*dNHatdpsiN + 10*MinFloat*(1 + 2*np.sign(dNHatdpsiN))
flux_peaking = max(2*abs(ABS_ERROR_PEAKING), MinFloat)

#Go to directory
os.chdir(directory)

if readInput("programMode", 0, sfincsInput) != 1:
	print "sfincsFindPeakingFactor.py has to be run with programMode = 1 in SFINCS"	
	sys.exit(0)

if CONTINUE_SCAN_PEAKING:
	dirIter = CONTINUE_DIR_PEAKING
	MAX_ITER_PEAKING = MAX_ITER_PEAKING + CONTINUE_DIR_PEAKING - 1 #Need to increase MAX_ITER_PEAKING so the effective number of iterations are correct
	if not os.path.exists("peaking" + str(dirIter - 1)):
		print "Error: Can not find directory %s in current directory %s" % ("peaking" + str(dirIter - 1), os.getcwd())
		sys.exit(0)

	flux_peaking = (readsfincsFindAmbipolarEr_Result("Fluxes", "peaking" + str(dirIter - 1) + "/" + sfincsFindAmbipolarEr_resultOutput))[-1]
	try:
		dNHatdpsiN = float(readInput("dNHatdpsiNs", 2, "peaking" + str(dirIter - 1) + "/" + sfincsInput).split()[-1].replace('d', 'e'))
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while reading dNHatdpsiNs in %s, control that it is an array of floats" % (os.getcwd() + "/peaking" + str(dirIter - 1) + "/" + sfincsInput)
		sys.exit(0)

	if dirIter < 3: #Only one iteration done earlier
		previousflux_peaking = max(2*abs(ABS_ERROR_PEAKING), MinFloat)
		previousdNHatdpsiN = 2*dNHatdpsiN + 10*MinFloat*(1 + 2*np.sign(dNHatdpsiN))
	else:
		previousflux_peaking = (readsfincsFindAmbipolarEr_Result("Fluxes", "peaking" + str(dirIter - 2) + "/" + sfincsFindAmbipolarEr_resultOutput))[-1]
		try:
			previousdNHatdpsiN = float(readInput("dNHatdpsiNs", 2, "peaking" + str(dirIter - 2) + "/" + sfincsInput).split()[-1].replace('d', 'e'))
		except Exception,e:
			print e.__class__.__name__, ": ", e.message, "while reading dNHatdpsiNs in %s, control that it is an array of floats" % (os.getcwd() + "/peaking" + str(dirIter - 2) + "/" + sfincsInput)
			sys.exit(0)
	
	if previousflux_peaking == flux_peaking or previousdNHatdpsiN == dNHatdpsiN:
		print "Error: Iteration scheme can not update dNHatdpsiN"
		sys.exit(0)

	#Update dNHatdpsiN
	temporary_dNHatdpsiN = previousdNHatdpsiN - (dNHatdpsiN - previousdNHatdpsiN) / (flux_peaking - previousflux_peaking) * previousflux_peaking
	previousdNHatdpsiN = dNHatdpsiN
	dNHatdpsiN = temporary_dNHatdpsiN	

#Get original lists of dNHatdpsiNs and Zs from SFINCS input file
ORIGINAL_dNHatdpsiNs_List = readInput("dNHatdpsiNs", 2, sfincsInput).replace('d', 'e').split()
try:
	ORIGINAL_dNHatdpsiNs_List = [float(i) for i in ORIGINAL_dNHatdpsiNs_List]
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while reading dNHatdpsiNs in %s, control that it is an array of floats" % (os.getcwd() + "/" + sfincsInput)
	sys.exit(0)
ORIGINAL_Zs_List = readInput("Zs", 2, sfincsInput).replace('d', 'e').split()
try:
	ORIGINAL_Zs_List = [int(i) for i in ORIGINAL_Zs_List]
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while reading Zs in %s, control that it is an array of integers" % (os.getcwd() + "/"  + sfincsInput)
	sys.exit(0)

if len(ORIGINAL_Zs_List) != len(ORIGINAL_dNHatdpsiNs_List): 
	print "Error: In %s Zs and dNHatdpsiNs seem to have different dimensions" % (os.getcwd() + "/"  + sfincsInput)


relative_error_peaking = abs(flux_peaking) + abs(REL_ERROR_PEAKING) #To keep track of relative error between iterations step

while dirIter <= MAX_ITER_PEAKING and (abs(ABS_ERROR_PEAKING) <= abs(flux_peaking) or abs(REL_ERROR_PEAKING) <= abs(relative_error_peaking)):
	#print "#######################################################################################"
	print "======================================================================================="
	print "Find Peaking Factor Iteration %i" %dirIter
	print "Last species density gradient %f" % dNHatdpsiN

	#Make directory
	if os.path.exists("peaking" + str(dirIter)):
		try:
			shutil.rmtree("peaking" + str(dirIter))
		except Exception,e:
			print e.__class__.__name__, ": ", e.message, "while trying to remove %s in %s, check access rights" % ("peaking" + str(dirIter), os.getcwd())
			sys.exit(0)
	os.mkdir("peaking" + str(dirIter))
	#Copy SFINCS input file into directory
	try:
		shutil.copyfile(sfincsInput, "peaking" + str(dirIter) + "/" +  sfincsInput)
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while trying to copy %s into directory %s, check access rights" % (sfincsInput, os.getcwd() + "peaking" + str(dirIter))
		sys.exit(0)


	#Update dNHatdpsiNs in SFINCS input file.
	dNHatdpsiNs_List = ORIGINAL_dNHatdpsiNs_List
	
	#Modify dNHatdpsiN for the last two species
	if len(ORIGINAL_dNHatdpsiNs_List) < 2: #Only one species so can not maintain quasineutrality
		dNHatdpsiNs_List = dNHatdpsiN
	else:
		dNHatdpsiNs_List[-1] = dNHatdpsiN
		if CHECK_QUASINEUTRALITY_PEAKING: #To fulfill quasineutrality we modify dNHatdpsiN for the last two species
			dNHatdpsiNs_List[-2] = - 1/ORIGINAL_Zs_List[-2] * ( np.array(ORIGINAL_Zs_List).dot(np.array(dNHatdpsiNs_List)) - ORIGINAL_Zs_List[-2]*dNHatdpsiNs_List[-2])

	writeInput("dNHatdpsiNs", str(dNHatdpsiNs_List).replace('[', ' ').replace(']', '').replace(',', ' '), "peaking" + str(dirIter) + "/" +  sfincsInput)


	#Enter directory and run sfincsFindAmbipolarEr.py
	try:
		os.chdir("peaking" + str(dirIter))
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while trying to enter %s from %s" % ("peaking" + str(dirIter), os.getcwd() )
		sys.exit(0)
	print "Entering %s and running sfincsFindAmbipolarEr.py" % os.getcwd() #+ "/" + str(dirIter)

	try:
		#subprocess.check_call(["python", "../sfincsFindAmbipolarEr.py", "../" + sfincsFindAmbipolarEr_INPUT])
	       	procAmbEr = subprocess.Popen(["python", "../sfincsFindAmbipolarEr.py", "../" + sfincsFindAmbipolarEr_INPUT], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		procAmbErOutput = procAmbEr.communicate()
	      	print procAmbErOutput[0]
		print procAmbErOutput[1]
	except subprocess.CalledProcessError,e:
		print e.__class__.__name__, ": ", e.message, "sfincsFindAmbipolarEr.py stopped with exit code %i" % e.returncode
		sys.exit(0)
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "when trying to run sfincsFindAmbipolarEr.py"
		sys.exit(0)

	#Return to base directory
	os.chdir(directory)

	#Read output and update error
	previousflux_peaking = flux_peaking

	flux_peaking = (readsfincsFindAmbipolarEr_Result("Fluxes", "peaking" + str(dirIter) + "/" + sfincsFindAmbipolarEr_resultOutput))[-1]

	relative_error_peaking = abs((flux_peaking - previousflux_peaking) / previousflux_peaking)
	print "Current absolute error: %e" % abs(flux_peaking)
	print "Current relative error: %e" % relative_error_peaking

	if abs(flux_peaking) == 0.0 or abs(relative_error_peaking) == 0.0: # Solution found or can not continue iteration
		dirIter += 1
		break

	#Update radial electric field
	#From the assumption of a linear dependence
	temporary_dNHatdpsiN = previousdNHatdpsiN - (dNHatdpsiN - previousdNHatdpsiN) / (flux_peaking - previousflux_peaking) * previousflux_peaking
	previousdNHatdpsiN = dNHatdpsiN
	dNHatdpsiN = temporary_dNHatdpsiN

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

List_flux_peaking = []
List_dNHatdpsiN = []
List_dirNo = []
List_dPhiHatdpsiN_amb = []

Result_fluxes = []
Result_Zs = []
Result_dNHatdpsiN = -1.0
Result_dPhiHatdpsiN_amb = 0.0

for dirIter in xrange(1, NumDir + 1):

	if not os.path.exists("peaking" + str(dirIter)):
		print "WARNING: Could not find directory %s in %s" % ("peaking" + str(dirIter), os.getcwd())
		continue

	List_dirNo.append(dirIter)
	
	Result_fluxes = readsfincsFindAmbipolarEr_Result("Fluxes", "peaking" + str(dirIter) + "/" + sfincsFindAmbipolarEr_resultOutput)
	Result_Zs = readsfincsFindAmbipolarEr_Result("Zs", "peaking" + str(dirIter) + "/" + sfincsFindAmbipolarEr_resultOutput)
	Result_dPhiHatdpsiN_amb = readsfincsFindAmbipolarEr_Result("dPhiHatdpsiN*", "peaking" + str(dirIter) + "/" + sfincsFindAmbipolarEr_resultOutput)[0]
	try:
		Result_dNHatdpsiN = float(readInput("dNHatdpsiNs", 2, "peaking" + str(dirIter) + "/" + sfincsInput).split()[-1].replace('d', 'e'))
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while reading dNHatdpsiNs in %s, control that it is an array of floats" % (os.getcwd() + "/peaking" + str(dirIter) + "/" + sfincsInput)
		sys.exit(0)
	
	List_flux_peaking.append( Result_fluxes[-1] )
	List_dNHatdpsiN.append(Result_dNHatdpsiN)
	List_dPhiHatdpsiN_amb.append(Result_dPhiHatdpsiN_amb)

#End for loop

try:
	with open(arrayOutput_peaking, 'w') as f:
		f.write("\nIteration:\t")
		f.write(str(List_dirNo))
		f.write("\nFlux_last_species:\t")
		f.write(str(List_flux_peaking))
		f.write("\ndNHatdpsiN_last_species:\t")
		f.write(str(List_dNHatdpsiN))
		f.write("\nAmbipolar_dPhiHatdpsiN:\t")
		f.write(str(List_dPhiHatdpsiN_amb))
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while trying to open %s, check access rights" % arrayOutput_peaking
	sys.exit(0)
try:
	with open(resultOutput_peaking, 'w') as f:
		f.write("\ndNHatdpsiN*:\t")
		f.write(str(Result_dNHatdpsiN))
		f.write("\ndPhiHatdpsiN*:\t")
		f.write(str(Result_dPhiHatdpsiN_amb))
		f.write("\nZs:\t")
		f.write(str(Result_Zs))
		f.write("\nFluxes:\t")
		f.write(str(Result_fluxes))
except Exception,e:
	print e.__class__.__name__, ": ", e.message, "while trying to open %s, check access rights" % resultOutput_peaking
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

Array_flux_peaking = np.array(List_flux_peaking)
Array_dNHatdpsiN = np.array(List_dNHatdpsiN)
Array_dirNo = np.array(List_dirNo)
Array_dPhiHatdpsiN_amb = np.array(List_dPhiHatdpsiN_amb)

print "Directories to plot: %s" % str(Array_dirNo)
print "Radial fluxes of last species: %s" % str(Array_flux_peaking)
print "Ambipolar radial electric fields: %s" % str(Array_dPhiHatdpsiN_amb)
print "Last species density gradients: %s" % str(Array_dNHatdpsiN)

fig = plt.figure(figsize=(25,60))

##########
##Plot 1
##########
ax1 = plt.subplot(5, 1, 1)
plt.plot(Array_dirNo, Array_flux_peaking, 'mx-', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Radial flux of last species per iteration step', weight='bold')
plt.ylabel(r'$\Gamma_{\alpha}$', fontsize=40)
plt.xlabel(r'Iteration', fontsize=30)
fluxMin = np.amin(np.absolute(Array_flux_peaking))
fluxMax = np.amax(np.absolute(Array_flux_peaking))
space_lin1 = 1
linlog_limit1 = MinFloat
if 0 < fluxMax and 0 < fluxMin:
	space_lin1 = np.maximum(np.ceil(np.log10(fluxMax / fluxMin)/4.0), 1.0)
	linlog_limit1 = np.power(10, np.ceil(np.log10(fluxMin)))
linlog_limit1_array = np.empty(Array_dirNo.size); linlog_limit1_array.fill(linlog_limit1)

plt.yscale('symlog', linthreshy = linlog_limit1, linscaley = space_lin1)
plt.fill_between(Array_dirNo, linlog_limit1_array, -linlog_limit1_array, facecolor='green', alpha=0.08)
plt.text(0.10, 0.5, r'Linear scale', color='green', horizontalalignment='center',verticalalignment='center', transform=ax1.transAxes)

##########
##Plot 2
##########
ax2 = plt.subplot(5, 1, 2)
Array_dNHatdpsiN_diff = Array_dNHatdpsiN - Array_dNHatdpsiN[-1]
plt.plot(Array_dirNo, Array_dNHatdpsiN_diff, 'rx-', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Last species density gradient per iteration step', weight='bold')
plt.xlabel(r'Iteration', fontsize=30)
plt.ylabel(r'$d \hat{n}_{\alpha} / d \psi_N - \left(d \hat{n}_{\alpha} / d \psi_N\right)^0$', fontsize=40)

dnMin = 0; dnMax = 0
if Array_dNHatdpsiN_diff.size > 1:
	dnMin = np.amin(np.absolute(Array_dNHatdpsiN_diff[0:-1]))
	dnMax = np.amax(np.absolute(Array_dNHatdpsiN_diff[0:-1]))
space_lin2 = 1
linlog_limit2 = MinFloat
if 0 < dnMax and 0 < dnMin:
	space_lin2 = np.maximum(np.ceil(np.log10(dnMax / dnMin)/4.0), 1.0)
	linlog_limit2 = np.power(10, np.ceil(np.log10(dnMin)))
linlog_limit2_array = np.empty(Array_dirNo.size); linlog_limit2_array.fill(linlog_limit2)

plt.yscale('symlog', linthreshy = linlog_limit2, linscaley = space_lin2)
plt.fill_between(Array_dirNo, linlog_limit2_array, -linlog_limit2_array, facecolor='green', alpha=0.08)
plt.text(0.10, 0.50, r'Linear scale', color='green', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes)
plt.text(0.80, 0.90, r'$\left(d \hat{n}_{\alpha} / d \psi_N\right)^0$ = %e' %Array_dNHatdpsiN[-1], color='black', horizontalalignment='center',verticalalignment='center', transform=ax2.transAxes)

##########
##Plot 3
##########
ax3 = plt.subplot(5, 1, 3)
plt.plot(Array_dirNo, Array_dPhiHatdpsiN_amb, 'cx-', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Ambipolar radial electric field per iteration step', weight='bold')
plt.xlabel(r'Iteration', fontsize=30)
plt.ylabel(r'$\left(d \hat{\Phi} / d \psi_N\right)^\ast$', fontsize=40)

##########
##Plot 4
##########
ax4 = plt.subplot(5, 1, 4)
plt.plot(Array_dNHatdpsiN_diff, Array_flux_peaking, 'bx', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Radial flux of last species as function of last species density gradient', weight='bold')
plt.xlabel(r'$d \hat{n}_{\alpha} / d \psi_N - \left(d \hat{n}_{\alpha} / d \psi_N\right)^0$', fontsize=40)
plt.ylabel(r'$\Gamma_{\alpha}$', fontsize=40)
p2, = plt.plot(Array_dNHatdpsiN_diff[-1], Array_flux_peaking[-1], color='red', linestyle='none', marker='o', markersize=15, markeredgewidth=0)
plt.legend([p2], [r'Last species peaking factor'], loc=0, numpoints=1)

plt.yscale('symlog', linthreshy = linlog_limit1, linscaley = space_lin1)
x_array_lim = np.empty(Array_dNHatdpsiN_diff.size); x_array_lim[0] = - linlog_limit2; x_array_lim[-1] = linlog_limit2;

plt.fill_between(x_array_lim, np.amin(Array_flux_peaking), np.amax(Array_flux_peaking), facecolor='green', alpha=0.08)
plt.fill_between(np.sort(Array_dNHatdpsiN_diff), linlog_limit1_array, -linlog_limit1_array, facecolor='green', alpha=0.08)

plt.fill_between(x_array_lim, linlog_limit1_array, -linlog_limit1_array, facecolor='green', alpha=0.16)
plt.text(0.10, 0.5, r'Linear scale', color='green', horizontalalignment='center',verticalalignment='center', transform=ax4.transAxes)
plt.xscale('symlog', linthreshx = linlog_limit2, linscalex = space_lin2)
plt.text(0.80, 0.90, r'$\left(d \hat{n}_{\alpha} / d \psi_N\right)^0$ = %e' %Array_dNHatdpsiN[-1], color='black', horizontalalignment='center',verticalalignment='center', transform=ax4.transAxes)

##########
##Plot 5
##########
ax5 = plt.subplot(5, 1, 5)
plt.plot(Array_dNHatdpsiN_diff, Array_dPhiHatdpsiN_amb, 'yx', markersize=20, markeredgewidth=5)
plt.grid(True)
plt.title('Ambipolar radial electric field as function of last species density gradient', weight='bold')
plt.xlabel(r'$d \hat{n}_{\alpha} / d \psi_N - \left(d \hat{n}_{\alpha} / d \psi_N\right)^0$', fontsize=40)
plt.ylabel(r'$\left(d \hat{\Phi} / d \psi_N\right)^\ast$', fontsize=40)

plt.fill_between(x_array_lim, np.amin(Array_dPhiHatdpsiN_amb), np.amax(Array_dPhiHatdpsiN_amb), facecolor='green', alpha=0.08)
plt.text(0.45, 0.15, r'Linear scale', color='green', horizontalalignment='center',verticalalignment='center', transform=ax5.transAxes)
plt.xscale('symlog', linthreshx = linlog_limit2, linscalex = space_lin2)
plt.text(0.80, 0.90, r'$\left(d \hat{n}_{\alpha} / d \psi_N\right)^0$ = %e' %Array_dNHatdpsiN[-1], color='black', horizontalalignment='center',verticalalignment='center', transform=ax5.transAxes)



print "\nWriting results of scan to %s" % plotOutput_peaking
plt.savefig(plotOutput_peaking)

sys.exit(0)












