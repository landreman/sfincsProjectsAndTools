import sys, os
from numpy import exp, swapaxes, log10, sqrt, pi
import string
import h5py

def readVariable(varName, inputFilename,intOrFloatOrString, required=True):
    # This function reads normal fortran variables from the input.namelist file.
    # It is assumed that the input.namelist file has been loaded into the variable "inputFile".

    inputFile = file(inputFilename,'r')
    if (intOrFloatOrString != "int") and (intOrFloatOrString != "float") and (intOrFloatOrString != "string"):
        print "intOrFloatOrString must be int, float, or string."
        exit(1)

    originalVarName = varName
    #varName = varName.lower()
    returnValue = None
    numValidLines = 0
    for line in inputFile:
        #line3 = line.strip().lower()
        line3 = line.strip()
        if len(line3)<1:
            continue

        if line3[0]=="!":
            continue

        if len(line3) < len(varName)+2:
            continue

        if not line3[:len(varName)].lower()==varName.lower():
            continue

        line4 = line3[len(varName):].strip()

        if not line4[0] =="=":
            continue

        line5 = line4[1:].strip();
        if intOrFloatOrString != "string":
            # python does not recognize fortran's 1d+0 scientific notation
            line5 = line5.replace('d','e').replace('D','e')

        # Remove any comments:
        if "!" in line5:
            line5 = line5[:string.find(line5,"!")]

        if intOrFloatOrString=="int":
            try:
                returnValue = int(line5)
                numValidLines += 1
            except:
##                print "Warning! I found a definition for the variable "+originalVarName+" in "+filename+" but I was unable to parse the line to get an integer." #Commented by AM 2015-12
                print "Warning! I found a definition for the variable "+originalVarName+" in "+inputFilename+" but I was unable to parse the line to get an integer." #Added by AM 2015-12
                print "Here is the line in question:"
                print line
        elif intOrFloatOrString=="float":
            try:
                returnValue = float(line5)
                numValidLines += 1
            except:
##                print "Warning! I found a definition for the variable "+originalVarName+" in "+filename+" but I was unable to parse the line to get a float." #Commented by AM 2015-12
                print "Warning! I found a definition for the variable "+originalVarName+" in "+inputFilename+" but I was unable to parse the line to get a float." #Added by AM 2015-12
                print "Here is the line in question:"
                print line
        elif intOrFloatOrString=="string":
            returnValue = line5
            numValidLines += 1

    ##if returnValue==None: #Commented by AM 2015-12
        ##if required: #Commented by AM 2015-12
    if required and (returnValue is None):
        print "Error! Unable to find a valid setting for the variable "+originalVarName+" in "+inputFilename+"."
##        raise #Commented by AM 2015-12
        exit(1) #Added by AM 2015-12

    if numValidLines > 1:
        print "Warning! More than 1 valid definition was found for the variable "+originalVarName+". The last one will be used."

    #print "Read "+originalVarName+" = "+str(returnValue)
    inputFile.close()
    return returnValue


def species_array_read(varName,inputFilename):
    raw = readVariable(varName,inputFilename, "string", required=True)
    raw = raw.replace('d','e')
    raws = raw.split()
    l = [float(e) for e in raws]
    return l


def read_namelist(inputFilename):
    #execfile(os.path.dirname(os.path.abspath(__file__))+'/../../../fortran/version3/utils//sfincsScan_common')

    Zs = species_array_read("Zs",inputFilename)
    mHats = species_array_read("mHats",inputFilename)
    nHats = species_array_read("nHats",inputFilename)
    THats = species_array_read("THats",inputFilename)
    inputRadialCoordinateForGradients = readVariable("inputRadialCoordinateForGradients",inputFilename, "int", required=False)
    if inputRadialCoordinateForGradients is None:
        inputRadialCoordinateForGradients = 4

    if inputRadialCoordinateForGradients == 0:
        dstring = "dpsiHats"
    elif inputRadialCoordinateForGradients == 1:
        dstring = "dpsiNs"
    elif inputRadialCoordinateForGradients == 2 or inputRadialCoordinateForGradients == 4:
        dstring = "drHats"
    elif inputRadialCoordinateForGradients == 3:
        dstring = "drNs"
    else:
        raise ValueError("inputRadialCoordinateForGradients must be 0,1,2,3 or 4")
    
    
    dNHatdY = species_array_read("dNHat"+dstring,inputFilename)
    dTHatdY = species_array_read("dTHat"+dstring ,inputFilename)
    if inputRadialCoordinateForGradients in [0,1,2,3]:
        dPhiHatdY = readVariable("dPhiHat"+dstring[:-1],inputFilename, "float", required=False)
    elif inputRadialCoordinateForGradients == 4:
        dPhiHatdY = -readVariable("Er",inputFilename, "float", required=True)
    else:
        raise ValueError("inputRadialCoordinateForGradients must be 0,1,2,3 or 4")
        
    Delta = readVariable("Delta",inputFilename, "float", required=True)
    alpha = readVariable("alpha",inputFilename, "float", required=True)
    nu_n = readVariable("nu_n",inputFilename, "float", required=True)
    
    if nu_n < 0:
        """ assume RBar = 1m, mBar = proton mass, nBar = 1e20/m^3, TBar = 1keV; electrons first species"""
        if (abs(Delta-4.5694e-3) /4.5694e-3) > 0.1:
            raise ValueError("nu_n will be calculated with standard bars, but Delta is wrong.")
        lnLambda = (25.3e+0) - (1.15e+0)*log10(nHats[0]*(1e14)) + (2.30e+0)*log10(THats[0]*1000)
        eC = 1.6022e-19
        epsilon0 = 8.8542e-12
        mproton = 1.6726e-27
        nu_n = sqrt(mproton/(2*1000*eC)) * 4*sqrt(2*pi)*(1e20)*(eC**4)*lnLambda / (3*((4*pi*epsilon0)**2)*sqrt(mproton)*((1000*eC)**(1.5e+0)))
    
            
    # Er probably incorrect

    inputRadialCoordinate = readVariable("inputRadialCoordinate",inputFilename, "int", required=True)
    rN_wish = readVariable("rN_wish",inputFilename, "float", required=True)
    equilibriumFile = readVariable("equilibriumFile",inputFilename,"string", required=True)[1:-1]
    if equilibriumFile[0] is '.':
        #relative path from the input file directory
        equilibriumFile = os.path.dirname(os.path.abspath(inputFilename)) +"/" + equilibriumFile
    
    min_Bmn_to_load = readVariable("min_Bmn_to_load",inputFilename, "float", required=False)
    if min_Bmn_to_load is None:
        min_Bmn_to_load = 0.0
    Ntheta= readVariable("Ntheta",inputFilename, "int", required=True)
    Nzeta= readVariable("Nzeta",inputFilename, "int", required=True)

    # to read Phi1
    includePhi1 = readVariable("includePhi1",inputFilename,"string",required=False)

    if includePhi1 is None:
        includePhi1 = False
    else:
        if includePhi1.lower() == ".false.":
            includePhi1 = False
        elif includePhi1.lower() == ".true.":
            includePhi1 = True
        else:
            raise ValueError("includePhi1 in input file must be either '.true.' or '.false.' (not case sensitive)")


    outputFilename = readVariable("outputFilename",inputFilename,"string", required=False)
    if outputFilename is None:
        outputFilename = "sfincsOutput.h5"
    else:
        if outputFilename[0] == '"' and outputFilename[-1] == '"':
            outputFilename = outputFilename[1:-1]

    outputFilename = os.path.dirname(os.path.abspath(inputFilename)) +"/" + outputFilename

    NHats = nHats[:]
    if includePhi1:
        Phi1Hat = h5py.File(outputFilename,'r')["/Phi1Hat"][:,:,-1]
        Phi1Hat = swapaxes(Phi1Hat,0,1)
        #Nspecies = len(nHats)
        #for i in range(Nspecies):
        #    nHats[i] = nHats[i] * exp(-Zs[i]*alpha*Phi1Hat/THats[i])
    else:
        Phi1Hat = None
    
    return (Zs,mHats,NHats,THats,dNHatdY, dTHatdY, dPhiHatdY,Delta,alpha,nu_n,inputRadialCoordinate,inputRadialCoordinateForGradients,rN_wish,equilibriumFile,min_Bmn_to_load,Ntheta,Nzeta,Phi1Hat)


if __name__=="__main__":
    input_filename = sys.argv[1]
    (Zs,mHats,NHats,THats,dNHatdrHats, dTHatdrHats, dPhiHatdrHats,Delta,alpha,nu_n,inputRadialCoordinate,inputRadialCoordinateForGradients,rN_wish,equilibriumFile,min_Bmn_to_load,Ntheta,Nzeta,Phi1Hat) = read_namelist(input_filename)
    print "Zs: " + str(Zs)
    print "mHats: " + str(mHats)
    print "NHats: " + str(NHats)
    print "THats: " + str(THats)
    print "dNHatdrHats: " + str(dNHatdrHats)
    print "dTHatdrHats: " + str(dTHatdrHats)
    print "dPhiHatdrHats: " + str(dPhiHatdrHats)
    print "Delta: " + str(Delta)
    print "alpha: " + str(alpha)
    print "nu_n: " + str(nu_n)
    print "inputRadialCoordinate: " + str(inputRadialCoordinate)
    print "inputRadialCoordinateForGradients: " + str(inputRadialCoordinateForGradients)
    print "rN_wish: " + str(rN_wish)
    print "equilibriumFile: " + equilibriumFile
    print "min_Bmn_to_load: " + str(min_Bmn_to_load)
    print "Ntheta: " + str(Ntheta)
    print "Nzeta: " + str(Nzeta)
    print "Phi1Hat: " + str(Phi1Hat)
