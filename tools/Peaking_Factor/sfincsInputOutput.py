
# sfincsInputOutput.py
#SFINCS Input/Output module, for reading from and writing to SFINCS input file

########################################################################
##Define functions for reading from and writing to SFINCS input.namelist
########################################################################


import sys
import os
import shutil
from tempfile import mkstemp

sfincsInput = "input.namelist" #Standard name of SFINCS input file


def readInput(variableName, intOrFloat=2, InputFile=sfincsInput): # set intOrFloat = 0 to read an integer, 1 to read a float or 2 to read a string.
	"Reads a variable (integer, float or string) from SFINCS input file"
	if not (os.path.exists(os.getcwd() + "/" + InputFile) or os.path.exists(InputFile) ):
		print "Error: %s not found in current directory %s" % (InputFile, os.getcwd())
		sys.exit(0)

	if (not isinstance(intOrFloat, int)) or intOrFloat > 2 or intOrFloat < 0:
		print "Error: intOrFloat must be 0, 1 or 2."
		sys.exit(0)


	#Read file and split each line into a list
	try:
		with open(InputFile, 'r') as f:
			lines = [line.strip().lower() for line in f]	
	except Exception,e:
		print e.__class__.__name__, ": ", e.message, "while trying to read %s, check access rights" % InputFile
		sys.exit(0)


	numMatches = 0
	value = 0
	for line in lines:
		if len(line) == 0:
			continue
		elif line[0] == "!" or line[0] == "&" or line[0] == "/": #Ignore comments
			continue
		elif intOrFloat == 2: #Read string
			if line.find("!") != -1: #Remove comment in end of line
				line = line[0:line.find("!")] 
			if len(line) <= len(variableName) or line[0:(len(variableName))] != variableName.lower() or line.find("=") == -1 or (''.join(line.split()))[len(variableName)] != '=':
				continue
			else:
				value = str(line.split("=", 1)[1])
				numMatches += 1
		else: #Read int or float
			line = ''.join(line.split()) #Remove all white spaces

			if line.find("!") != -1: #Remove comment in end of line
				line = line[0:line.find("!")] 

			if len(line) <= len(variableName) or line[0:(len(variableName))] != variableName.lower():
				continue
			elif line[len(variableName)] != '=':
				continue
			else:
				numMatches += 1
				value = line.split("=", 1)[1]
				value = value.replace("d", "e") #Change from Fortran format to Python
				if value == ".true.": #Change from Fortran format to Python
					value = True
				if value == ".false.": #Change from Fortran format to Python
					value = False

				try:
					if intOrFloat: #float				
						value = float(value)
					else: #int
						try:
							value = int(value)
						except ValueError: 
							value = int(float(value))
				except Exception,e:
					print e.__class__.__name__, ": ", e.message, "Problem when reading value %s" % value
					sys.exit(0)

	if numMatches < 1:
		print "WARNING: No occurrence for %s was found in %s, returning False" % (variableName, InputFile)
		return False
	elif numMatches > 1:
		print "WARNING: several occurrences of %s was found in %s, returning the last" % (variableName, InputFile)

	print "Read %s = %s" % (variableName, str(value))
	return value

def writeInput(variableName, value, InputFile=sfincsInput): # variableName is a string, value is a string, float or int.
	"Writes variableName = value into SFINCS input file"
	if not (os.path.exists(os.getcwd() + "/" + InputFile) or os.path.exists(InputFile) ):
		print "Error: %s not found in current directory %s" % (InputFile, os.getcwd())
		sys.exit(0)

	variableName = variableName.strip()
	value = str(value)

	source_file_path = os.path.abspath(InputFile)
	fh, target_file_path = mkstemp()
	with open(target_file_path, 'w') as target_file:
		try:
			with open(source_file_path, 'r') as source_file:
				numMatches = 0
				for line in source_file:
					if len(line.strip()) == 0:
						target_file.write(line)
					elif len(line.strip()) < len(variableName) or (line.strip())[0] == "!" or (line.strip())[0] == "&" or (line.strip())[0] == "/" or line.find(variableName) == -1 or (line.strip()).find(variableName) != 0  or line.find("=") == -1 or (line.find(variableName + " ") == -1 and line.find(variableName + "=") == -1) or (line.find("!") != -1 and line.find("!") < line.find("=")): #Ignore comments
						target_file.write(line)
					else: #variableName assignment found and should be replaced

						if line.find("!") == -1:
							newline = line.split("=", 1)[0] + "= " + value + "\n"
						else:
							newline = line.split("=", 1)[0] + "= " + value + " !" + line.split("!", 1)[1] #+ "\n"

						target_file.write(newline)
						numMatches += 1

				if numMatches < 1:
					print "No occurrence for %s was found in %s, adding %s = %s in end of file" % (variableName, InputFile, variableName, value)
					target_file.write(" " + variableName + " = " + value + "\n")
				else:
					print "The value of %s was replaced with %s at %i locations" % (variableName, value, numMatches)
						
		except Exception,e:
			print e.__class__.__name__, ": ", e.message, "while trying to read %s, check access rights" % source_file_path
			sys.exit(0)

	os.remove(source_file_path)
	shutil.move(target_file_path, source_file_path)


