#!/usr/bin/env ruby

# This ruby script launches a bunch of single species fortran sfincs jobs
# in a batch system, with the capability to perform several types of scans.

# This version includes a few more scan types than the basic one 
# in the sfincs repository, and produces .sge files for launch 
# on the IPP cluster

# Scan types (i.e. values for programMode) implemented in this script:
# 2  = Convergence scan
# 6  = Scan of Estar, at fixed resolutions.
# 7  = Simultaneous convergence scan and scan of Estar
# 8  = Scan of dPhiHatdpsi, at fixed resolutions.
# 9  = Simultaneous convergence scan and scan of dPhiHatdpsi
# 10 = Simultaneous scan of dphiHatdpsi, of collisionOperator=0 and 1, and of DKES vs full trajectories, keeping resolutions fixed.
# 20 = Scan radius, taking profiles from profiles.dat
# 21 = make a run for each line in the file runspec.dat

$filename="input.namelist"

$sfincs_system=ENV["SFINCS_SYSTEM"]
case $sfincs_system
when "hgw"
  $jobSubmitCommand="qsub"
when "hydra"
  $jobSubmitCommand="llsubmit"
else
  puts "Error! SFINCS_SYSTEM=#{$sfincs_system} is not implemented in the ruby script!"
      exit
end
jobFilename="job.#{$sfincs_system}"

$profilesdatfile="profiles.dat"
$runspecfile="runspec.dat"

#The following line is only relevant on the hgw cluster:
clusters=["wsclus","wsclus","edge01","edge01"]
#clusters=["wsclus","wsclus","edge01","wsclus","edge02","wsclus","wsclus","edge03"]
sizeFactor=0.1e-3

require 'fileutils'
include FileUtils        
            
if !File.exists?($filename)
  puts "Error! #{$filename} not found."
  exit
end
puts "File #{$filename} exists."

if !File.exists?(jobFilename)
  puts "Error! #{jobFilename} not found."
  exit
end
puts "File #{jobFilename} exists."

# Make sure job file does not already contain "#PBS -N"
inFile = File.open(jobFilename,"r")
lines = inFile.readlines
lines.each {|line| 
  case $sfincs_system
  when "hgw"
    if (line.include? "-N")
      puts "Error! #{jobFilename} should not include a -N line with job name."
      exit
    end
  when "hydra"
    #if (line.include? "output")
    #  puts "Error! #{jobFilename} should not include a output line with job name."
    #end
    #if (line.include? "error")
    #  puts "Error! #{jobFilename} should not include an error line with job name."
    #end
    if (line.include? "job_name")
      puts "Error! #{jobFilename} should not include a job_name line with job name."
    end
  else
      puts "Error! SFINCS_SYSTEM=#{$sfincs_system} is not implemented in the ruby script!"
      exit
  end
}
inFile.close

def namelistLineContains(line, variableName)
  line = line.strip.downcase
  variableName = variableName.downcase
  
  if (!line.include? variableName)
    return false 
  end
  if (line[0].chr == "!") 
    return false 
  end
  
  nextChar = line[variableName.size].chr
  if (line[0..(variableName.size-1)] == variableName) & ((nextChar==" ") | (nextChar == "="))
    return true
  else
    return false
  end
end

def readInput(variableName, intOrFloat)
  # set intOrFloat = 0 to read an integer or 1 to read a float.
  
  if !variableName.is_a?(String)
    puts "Error: variableName must be a string."
    exit
  end
  
  if !intOrFloat.is_a?(Integer)
    puts "Error: intOrFloat must be an integer."
    exit
  end
  if (intOrFloat > 1) | (intOrFloat < 0)
    puts "Error: intOrFloat must be 0 or 1."
    exit
  end
  
  variableName = variableName.downcase
  s = `grep -i #{variableName} #{$filename}`
  numMatches = 0
  value = 0
  s.each_line do |line|
    line = line.strip.downcase
    # If not a comment:
    if line[0].chr != "!"
      nextChar = line[variableName.size].chr
      if (line[0..(variableName.size-1)] == variableName) & ((nextChar==" ") | (nextChar == "="))
        # Ruby doesn't understand the fortran scientific notation 1d-0 so replace it with 1e-0.
        substring = line[(line.index("=")+1)..(line.size-1)]
        substring.gsub!("d-","e-")
        #value = line.scan(/\d+/)[0].to_f
        if intOrFloat==0
          value = substring.to_i
        else
          value = substring.to_f
        end
        numMatches = numMatches + 1
      end
    end
  end
  if numMatches < 1
    puts "Error! No lines in #{$filename} match #{variableName}."
    exit
  end
  if numMatches > 1
    puts "Warning: more than 1 line in #{$filename} matches #{variableName}."
  end
  if intOrFloat == 0
    puts "Read #{variableName} = " + value.to_i.to_s
    return value.to_i
  else
    puts "Read #{variableName} = " + value.to_s
    return value
  end
end

def readProfiles(profilefilename, noElems)
  if !File.exists?(profilefilename)
    puts "Error! #{profilefilename} not found."
    exit
  end
  puts "File #{profilefilename} exists."
  
  inFile = File.open(profilefilename,"r")
  lines = inFile.readlines
  inFile.close

  radparams=Array.new(lines.size){Array.new(noElems)}
  lind = 0
  lines.each do |line|
    if line.length>1
      if line[0].chr != "!"
        #puts line
        substr=line
        varind=0
        while varind<noElems-1 do
          varstr=substr[0..substr.index(" ")-1]
          radparams[lind][varind]=varstr.to_f
          substr=substr[(substr.index(" ")+1)..(substr.size-1)]
          while substr.index(" ")==0
            substr=substr[1..(substr.size-1)]
          end
          #puts "Element #{lind},#{varind} is #{radparams[lind][varind]}"
          varind += 1
        end
        radparams[lind][noElems-1]=substr.to_f
        #puts "Element #{lind},#{noElems-1} is #{radparams[lind][noElems-1]}"
        lind += 1
      end
    end
  end
  return radparams[0..(lind-1)]
end

def readRunspec(runspecfilename)
  if !File.exists?(runspecfilename)
    puts "Error! #{runspecfilename} not found."
    exit
  end
  puts "File #{runspecfilename} exists."
  
  inFile = File.open(runspecfilename,"r")
  lines = inFile.readlines
  inFile.close

  #runparams=Array.new(lines.size){Array.new(noElems)}
  lind = 0
  while lines[lind][0].chr =="!"
    lind=lind+1
  end
  #now read the names of the variables to read
  nameline=lines[lind-1][1..-1]
  runparamnames=nameline.split
  noElems=runparamnames.size
  datalines=lines[lind..-1]
  runparams=Array.new(datalines.size){Array.new(noElems)}
  lind=0
  datalines.each do |line|
    if line.length>1
      if line[0].chr != "!"
        #puts line
        substr=line
        varind=0
        while varind<noElems-1 do
          varstr=substr[0..substr.index(" ")-1]
          runparams[lind][varind]=varstr.to_f
          substr=substr[(substr.index(" ")+1)..(substr.size-1)]
          while substr.index(" ")==0
            substr=substr[1..(substr.size-1)]
          end
          #puts "Element #{lind},#{varind} is #{runparams[lind][varind]}"
          varind += 1
        end
        runparams[lind][noElems-1]=substr.to_f
        #puts "Element #{lind},#{noElems-1} is #{runparams[lind][noElems-1]}"
        lind += 1
      end
    end
  end
  return runparamnames, runparams[0..(lind-1)]
end

def linspace(min, max, nn)
  return (0..(nn-1)).collect{|x| x*(max-min)/(nn-1.0)+min}
end

def logspace(min, max, nn)
  if (min <= 0)
    raise "Error! In logspace, min must be positive"
  end
  if (max <= 0)
    raise "Error! In logspace, max must be positive"
  end
  logs = linspace(Math.log(min), Math.log(max), nn)
  return logs.collect {|x| Math.exp(x)}
end

##############################################################
# Load the data from input.namelist and prepare order of runs
##############################################################
programMode = readInput("programMode",0)

case programMode
when 2,7,9
  puts "Beginning convergence scan"
  
  Ntheta = readInput("Ntheta",0)
  NthetaMinFactor = readInput("NthetaMinFactor",1)
  NthetaMaxFactor = readInput("NthetaMaxFactor",1)
  NthetaNumRuns = readInput("NthetaNumRuns",0)
  Nthetas_tmp = logspace(Ntheta*NthetaMinFactor, Ntheta*NthetaMaxFactor, NthetaNumRuns).collect{|x| x.round}
  # Force Ntheta to be odd:
  Nthetas = Nthetas_tmp.collect{|x| if (x % 2 == 1) then x else x+1 end}
  puts "Nthetas:"
  p Nthetas
  
  Nzeta = readInput("Nzeta",0)
  NzetaMinFactor = readInput("NzetaMinFactor",1)
  NzetaMaxFactor = readInput("NzetaMaxFactor",1)
  NzetaNumRuns = readInput("NzetaNumRuns",0)
  Nzetas_tmp = logspace(Nzeta*NzetaMinFactor, Nzeta*NzetaMaxFactor, NzetaNumRuns).collect{|x| x.round}
  # Force Nzeta to be odd:
  Nzetas = Nzetas_tmp.collect{|x| if (x % 2 == 1) then x else x+1 end}
  puts "Nzetas:"
  p Nzetas
  
  Nxi = readInput("Nxi",0)
  NxiMinFactor = readInput("NxiMinFactor",1)
  NxiMaxFactor = readInput("NxiMaxFactor",1)
  NxiNumRuns = readInput("NxiNumRuns",0)
  Nxis = logspace(Nxi*NxiMinFactor, Nxi*NxiMaxFactor, NxiNumRuns).collect{|x| x.round}
  puts "Nxis:"
  p Nxis
  
  Nx = readInput("Nx",0)
  NxMinFactor = readInput("NxMinFactor",1)
  NxMaxFactor = readInput("NxMaxFactor",1)
  NxNumRuns = readInput("NxNumRuns",0)
  Nxs = logspace(Nx*NxMinFactor, Nx*NxMaxFactor, NxNumRuns).collect{|x| x.round}
  puts "Nxs:"
  p Nxs
  
  NL = readInput("NL",0)
  NLMinFactor = readInput("NLMinFactor",1)
  NLMaxFactor = readInput("NLMaxFactor",1)
  NLNumRuns = readInput("NLNumRuns",0)
  NLs = logspace(NL*NLMinFactor, NL*NLMaxFactor, NLNumRuns).collect{|x| x.round}
  puts "NLs:"
  p NLs
  
  NxPotentialsPerVth = readInput("NxPotentialsPerVth",1)
  NxPotentialsPerVthMinFactor = readInput("NxPotentialsPerVthMinFactor",1)
  NxPotentialsPerVthMaxFactor = readInput("NxPotentialsPerVthMaxFactor",1)
  NxPotentialsPerVthNumRuns = readInput("NxPotentialsPerVthNumRuns",0)
  NxPotentialsPerVths = logspace(NxPotentialsPerVth*NxPotentialsPerVthMinFactor, NxPotentialsPerVth*NxPotentialsPerVthMaxFactor, NxPotentialsPerVthNumRuns)
  puts "NxPotentialsPerVths:"
  p NxPotentialsPerVths
  
  xMax = readInput("xMax",1)
  xMaxMinFactor = readInput("xMaxMinFactor",1)
  xMaxMaxFactor = readInput("xMaxMaxFactor",1)
  xMaxNumRuns = readInput("xMaxNumRuns",0)
  xMaxs = logspace(xMax*xMaxMinFactor, xMax*xMaxMaxFactor, xMaxNumRuns)
  puts "xMaxs:"
  p xMaxs
  
  solverTolerance = readInput("solverTolerance",1)
  solverToleranceMinFactor = readInput("solverToleranceMinFactor",1)
  solverToleranceMaxFactor = readInput("solverToleranceMaxFactor",1)
  solverToleranceNumRuns = readInput("solverToleranceNumRuns",0)
  solverTolerances = logspace(solverTolerance*solverToleranceMinFactor, solverTolerance*solverToleranceMaxFactor, solverToleranceNumRuns)
  puts "solverTolerances:"
  p solverTolerances
  
  
  numRunsInScan = 1 + NthetaNumRuns + NzetaNumRuns + NxiNumRuns \
  + NLNumRuns + NxNumRuns + NxPotentialsPerVthNumRuns + xMaxNumRuns + solverToleranceNumRuns
  
  baseCase = [Ntheta,Nzeta,Nxi,NL,Nx,NxPotentialsPerVth,xMax,solverTolerance];
  parametersForScan = Array.new
  descriptions = Array.new
  for i in 1..numRunsInScan
    # Need to use .dup since otherwise each element of the outer array will point to the same instance of the inner array.
    parametersForScan << baseCase.dup
  end
  
  currentIndex = 1
  descriptions[0]="baseCase"
  
  for i in 1..NthetaNumRuns
    parametersForScan[currentIndex][0] = Nthetas[i-1]
    descriptions[currentIndex] = "Ntheta" + Nthetas[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NzetaNumRuns
    parametersForScan[currentIndex][1] = Nzetas[i-1]
    descriptions[currentIndex] = "Nzeta" + Nzetas[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NxiNumRuns
    parametersForScan[currentIndex][2] = Nxis[i-1]
    descriptions[currentIndex] = "Nxi" + Nxis[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NLNumRuns
    parametersForScan[currentIndex][3] = NLs[i-1]
    descriptions[currentIndex] = "NL" + NLs[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NxNumRuns
    parametersForScan[currentIndex][4] = Nxs[i-1]
    descriptions[currentIndex] = "Nx" + Nxs[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NxPotentialsPerVthNumRuns
    parametersForScan[currentIndex][5] = NxPotentialsPerVths[i-1]
    descriptions[currentIndex] = "NxPotentialsPerVth" + NxPotentialsPerVths[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..xMaxNumRuns
    parametersForScan[currentIndex][6] = xMaxs[i-1]
    descriptions[currentIndex] = "xMax" + xMaxs[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..solverToleranceNumRuns
    parametersForScan[currentIndex][7] = solverTolerances[i-1]
    descriptions[currentIndex] = "solverTolerance" + solverTolerances[i-1].to_s
    currentIndex += 1
  end
  
  if (currentIndex != numRunsInScan)
    puts "Error! Something went wrong."
    exit
  end
  
  # Now remove any duplicates
  
  i = 0
  while i < numRunsInScan-1
    j=i+1
    while j < numRunsInScan
      isThisADuplicate = true
      for k in 0...(parametersForScan[0].size)
        # If values are within 0.1%, the difference is probably just due to roundoff error,
        # so treat the values as the same:
        if (parametersForScan[i][k] - parametersForScan[j][k]*1.0).abs/parametersForScan[i][k] > 0.001
          isThisADuplicate = false
        end
      end
      if isThisADuplicate
        # Item j is a duplicate, so remove it.
        parametersForScan.delete_at(j)
        descriptions.delete_at(j)
        numRunsInScan -= 1
        j -= 1
      end
      j += 1
    end
    i += 1
  end
  
  puts "Parameters for scan:"
  p parametersForScan
  
  puts "Description of runs:"
  p descriptions
  
when 6
  # Scan of E_r (i.e. of Estar)
  
  NEStar = readInput("NEStar",0)
  EStarMin = readInput("EStarMin",1)
  EStarMax = readInput("EStarMax",1)
  EStars = logspace(EStarMin, EStarMax, NEStar)
  
  numRunsInScan = NEStar
  
  puts "EStars:"
  p EStars

when 8
  # Scan of E_r (i.e. of dPhiHatdpsi)
  
  # The following is ugly, but I do it to avoid changing in sfincs.
  # I call dPhiHatdpsi parameters EStar in input.parameters
  # and rename them here.
  #NErs = readInput("NErs",0)
  #dPhiHatdpsi_min = readInput("dPhiHatdpsi_min",1)
  #dPhiHatdpsi_max = readInput("dPhiHatdpsi_max",1)
  NErs = readInput("NEStar",0)
  dPhiHatdpsi_min = readInput("EStarMin",1)
  dPhiHatdpsi_max = readInput("EStarMax",1)
  dPhiHatdpsis = linspace(dPhiHatdpsi_min, dPhiHatdpsi_max, NErs)
  
  numRunsInScan = NErs
  
  puts "dPhiHatdpsis:"
  p dPhiHatdpsis

when 10
  descriptions = ["FP_full", "FP_DKES", "PAS_full", "PAS_DKES"]
  numRunsInScan = 4
  
when 20
  parametersForScan=readProfiles($profilesdatfile, 9)
  numRunsInScan = parametersForScan.size

when 21
  parameternamesForScan, parametersForScan=readRunspec($runspecfile)
  numRunsInScan = parametersForScan.size

else
  puts "I do not know what to do with programMode = "+programMode.to_s
  exit
end


##############################################################
# Make subdirs and write input.namelist in them, then launch.
##############################################################
case programMode
when 7,9,10
  # A 2-level scan, with Er as the inner directory.
  
  case programMode
  when 7
    NErs = readInput("NEStar",0)
    EStarMin = readInput("EStarMin",1)
    EStarMax = readInput("EStarMax",1)
    EStars = logspace(EStarMin, EStarMax, NErs)
    
    puts "EStars:"
    p EStars
    
  else
    NErs = readInput("NErs",0)
    dPhiHatdpsi_min = readInput("dPhiHatdpsi_min",1)
    dPhiHatdpsi_max = readInput("dPhiHatdpsi_max",1)
    dPhiHatdpsis = linspace(dPhiHatdpsi_min, dPhiHatdpsi_max, NErs)
    
    puts "dPhiHatdpsis:"
    p dPhiHatdpsis
  end

  puts "Scan will consist of #{numRunsInScan} runs for each of #{NErs} Ers, for a total of #{numRunsInScan*NErs} runs."
  
  for i in 0...numRunsInScan
    outerDirName = descriptions[i]
    if !File.exists?(outerDirName)
      mkdir(outerDirName)
    end
    case programMode
    when 7,9
      # Calculate  discrSize=Ntheta*Nzeta*Nxi*Nx
      discrSize=parametersForScan[i][0]*parametersForScan[i][1]*parametersForScan[i][2]*parametersForScan[i][4]
    end

    puts "Beginning to submit jobs for "+descriptions[i]
    
    for k in 0...NErs
      innerDirNum = k-1
      # Starting with 0, look for the first unused whole number to use for a directory name:
      begin
        innerDirNum += 1
      	innerDirName = innerDirNum.to_s
        if innerDirNum < 10
          innerDirName = "0" + innerDirName
        end
        dirName = outerDirName + "/" + innerDirName
      end until !File.exists?(dirName)
      mkdir dirName
      
      # Copy pbs file
      outFilename = dirName + "/" + jobFilename
      inFile = File.open(jobFilename,"r")
      outFile = File.open(outFilename,"w")
      lines = inFile.readlines
      outFile.write(lines[0])
      case $sfincs_system 
      when "hgw"
        outFile.write("#\$ -N #{outerDirName}.#{innerDirName}\n")
        mem_phys=0
        for j in 1..(lines.size-1)
          if lines[j].include? "mem_phys"
            mem_phys=lines[j][lines[j].index('=')+1..lines[j].index('G')-1]
          end
          if lines[j].include? "CLUSTR"
            clind=lines[j].index('CLUSTR')
            iall=i*NErs+k
            lines[j][clind..clind+5]=clusters[iall.modulo(clusters.size)]
            puts lines[j]
          end 
          case programMode
          when 7,9
            if lines[j].include? "ompi_b"
              if mem_phys>0
                lines[j][lines[j].index('ompi_b')+6..-1]=4*((discrSize*sizeFactor/mem_phys)/4).ceil
              else
                puts "WARNING: mem_phys and processor line in the wrong order!"
              end
            end
          end
          outFile.write(lines[j])
        end
      when "hydra"
        #outFile.write("# @ error = #{outerDirName}.#{innerDirName}.e$(jobid)\n")
        #outFile.write("# @ output = #{outerDirName}.#{innerDirName}.o$(jobid)\n")
        outFile.write("# @ job_name = #{outerDirName}.#{innerDirName}.j$(jobid)\n")
	for j in 1..(lines.size-1)
          outFile.write(lines[j])
        end
      end
      inFile.close
      outFile.close
      
      # Copy input.namelist
      outFilename = dirName + "/" + $filename
      inFile = File.open($filename,"r")
      outFile = File.open(outFilename,"w")
      lines = inFile.readlines
      for j in 0..(lines.size-1)
        line = lines[j].strip
        
        if namelistLineContains(line,"programMode")
          line = "programMode = 1"
        end

        case programMode
        when 7,9
          if namelistLineContains(line,"Ntheta")
            line = "Ntheta = " + parametersForScan[i][0].to_s
          end
          
          if namelistLineContains(line,"Nzeta")
            line = "Nzeta = " + parametersForScan[i][1].to_s
          end
          
          if namelistLineContains(line,"Nxi")
            line = "Nxi = " + parametersForScan[i][2].to_s
          end
          
          if namelistLineContains(line,"NL")
            line = "NL = " + parametersForScan[i][3].to_s
          end
          
          if namelistLineContains(line,"Nx")
            line = "Nx = " + parametersForScan[i][4].to_s
          end
          
          if namelistLineContains(line,"NxPotentialsPerVth")
            line = "NxPotentialsPerVth = " + parametersForScan[i][5].to_s
          end
          
          if namelistLineContains(line,"xMax")
            line = "xMax = " + parametersForScan[i][6].to_s
          end
          
          if namelistLineContains(line,"solverTolerance")
            line = "solverTolerance = " + parametersForScan[i][7].to_s
          end

        when 10
          if i==0 or i==1
            # Use Fokker-Planck collisions:
            if namelistLineContains(line,"collisionOperator")
              line = "collisionOperator = 0"
            end
          else
            # Use pure pitch-angle scattering collisions:
            if namelistLineContains(line,"collisionOperator")
              line = "collisionOperator = 1"
            end
          end

          if i==0 or i==2
            # Use full trajectories:

            if namelistLineContains(line,"includeXDotTerm")
              line = "includeXDotTerm = .true."
            end
            if namelistLineContains(line,"includeElectricFieldTermInXiDot")
              line = "includeElectricFieldTermInXiDot = .true."
            end
            if namelistLineContains(line,"useDKESExBDrift")
              line = "useDKESExBDrift = .false."
            end
            if namelistLineContains(line,"include_fDivVE_term")
              line = "include_fDivVE_term = .false."
            end

          else
            # Use DKES trajectories:

            if namelistLineContains(line,"includeXDotTerm")
              line = "includeXDotTerm = .false."
            end
            if namelistLineContains(line,"includeElectricFieldTermInXiDot")
              line = "includeElectricFieldTermInXiDot = .false."
            end
            if namelistLineContains(line,"useDKESExBDrift")
              line = "useDKESExBDrift = .true."
            end
            if namelistLineContains(line,"include_fDivVE_term")
              line = "include_fDivVE_term = .false."
            end

          end

        else
          raise "Program should not get here."
        end

        case programMode
        when 7
          if namelistLineContains(line,"EStar")
            line = "EStar = " + EStars[k].to_s
          end          
        else #9 or 10
          if namelistLineContains(line,"dPhiHatdpsi")
            line = "dPhiHatdpsi = " + dPhiHatdpsis[k].to_s
          end
        end
        
        outFile.write(line + "\n")
      end
      inFile.close
      outFile.close
      
      # Submit job!
      puts "Submitting job #{innerDirName} in the Er scan."
      # Make sure base case jobs get submitted first:
      #if i==0 and programMode == 9
      if true
        # In this way of submitting a job, ruby waits for qsub to complete before moving on.
        puts `cd #{dirName}; #{$jobSubmitCommand} #{jobFilename} &`
      else
        # In this way of submitting a job, ruby does not wait for qsub to complete before moving on.
        job1 = fork do
          exec "cd #{dirName}; #{$jobSubmitCommand} #{jobFilename}"
        end
        Process.detach(job1)
      end
      puts "Done submitting job #{innerDirName}"
    end # of loop j over Ers
  end # of loop i over convergence scan
  
else
  # A standard scan, like programMode = 2, 6, 8, 20 or 21 but not 7, 9 or 10.
  
  puts "Scan will consist of #{numRunsInScan} runs."
  
  
  for i in 0..(numRunsInScan-1)
    dirNum = i-1
    # Starting with 0, look for the first unused whole number to use for a directory name:
    begin
      dirNum += 1
      dirName = dirNum.to_s
      if dirNum < 10
        dirName = "0" + dirName
      end
    end until !File.exists?(dirName)
    mkdir(dirName)
    
    case programMode
    when 2
      jobName="SfxConv.#{dirName}"
      # Calculate  discrSize=Ntheta*Nzeta*Nxi*Nx
      discrSize=parametersForScan[i][0]*parametersForScan[i][1]*parametersForScan[i][2]*parametersForScan[i][4]
    when 6
      jobName="SfxEstr.#{dirName}"
    when 8
      jobName="SfxEr.#{dirName}"
    when 20
      jobName="SfxRad.#{dirName}"
    when 21
      jobName="Sfx.#{dirName}"
    else
      jobName="Sfx.#{dirName}"
    end

    # Copy sge file (pbs file in Matt's case)
    outFilename = dirName + "/" + jobFilename
    inFile = File.open(jobFilename,"r")
    outFile = File.open(outFilename,"w")
    lines = inFile.readlines
    outFile.write(lines[0])
    case $sfincs_system 
    when "hgw"
      outFile.write("#\$ -N #{jobName}\n")
      mem_phys=0
      for j in 1..(lines.size-1)
        if lines[j].include? "mem_phys"
          mem_phys=lines[j][lines[j].index('=')+1..lines[j].index('G')-1]
        end
        if lines[j].include? "CLUSTR"
          clind=lines[j].index('CLUSTR')
          lines[j][clind..clind+5]=clusters[i.modulo(clusters.size)]
          puts lines[j]
        end 
        case programMode
        when 2
          if lines[j].include? "ompi_b"
            if mem_phys>0
              lines[j][lines[j].index('ompi_b')+6..-1]=4*((discrSize*sizeFactor/mem_phys)/4).ceil
            else
              puts "WARNING: mem_phys and processor line in the wrong order!"
            end
          end
        end
        outFile.write(lines[j])
      end
    when "hydra"
      outFile.write("# @ error = #{jobName}.e$(jobid)\n")
      outFile.write("# @ output = #{jobName}.o$(jobid)\n")
      for j in 1..(lines.size-1)
        outFile.write(lines[j])
      end
    end

    inFile.close
    outFile.close
    
    # Copy input.namelist
    outFilename = dirName + "/" + $filename
    inFile = File.open($filename,"r")
    outFile = File.open(outFilename,"w")
    lines = inFile.readlines
    for j in 0..(lines.size-1)
      line = lines[j].strip
      
      if namelistLineContains(line,"programMode")
        line = "programMode = 1"
      end
      
      case programMode
      when 2
        if namelistLineContains(line,"Ntheta")
          line = "Ntheta = " + parametersForScan[i][0].to_s
        end
        
        if namelistLineContains(line,"Nzeta")
          line = "Nzeta = " + parametersForScan[i][1].to_s
        end
        
        if namelistLineContains(line,"Nxi")
          line = "Nxi = " + parametersForScan[i][2].to_s
        end
        
        if namelistLineContains(line,"NL")
          line = "NL = " + parametersForScan[i][3].to_s
        end
        
        if namelistLineContains(line,"Nx")
          line = "Nx = " + parametersForScan[i][4].to_s
        end
        
        if namelistLineContains(line,"NxPotentialsPerVth")
          line = "NxPotentialsPerVth = " + parametersForScan[i][5].to_s
        end
        
        if namelistLineContains(line,"xMax")
          line = "xMax = " + parametersForScan[i][6].to_s
        end
        
        if namelistLineContains(line,"solverTolerance")
          line = "solverTolerance = " + parametersForScan[i][7].to_s
        end
        
      when 6
        if namelistLineContains(line,"EStar")
          line = "EStar = " + EStars[i].to_s
        end
        
      when 8
        if namelistLineContains(line,"dPhiHatdpsi")
          line = "dPhiHatdpsi = " + dPhiHatdpsis[i].to_s
        end
        
      when 20
        if namelistLineContains(line,"normradius_wish")
          line = "normradius_wish = " + parametersForScan[i][0].to_s
        end

        if namelistLineContains(line,"THat")
          line = "THat = " + parametersForScan[i][1].to_s
        end

        if namelistLineContains(line,"nHat")
          line = "nHat = " + parametersForScan[i][2].to_s
        end

        if namelistLineContains(line,"dTHatdpsi")
          line = "dTHatdpsi = " + parametersForScan[i][3].to_s
        end

        if namelistLineContains(line,"dnHatdpsi")
          line = "dnHatdpsi = " + parametersForScan[i][4].to_s
        end

        if namelistLineContains(line,"dPhiHatdpsi")
          line = "dPhiHatdpsi = " + parametersForScan[i][5].to_s
        end

        if namelistLineContains(line,"nuN")
          line = "nuN = " + parametersForScan[i][6].to_s
        end

        if namelistLineContains(line,"Estar")
          line = "Estar = " + parametersForScan[i][7].to_s
        end

        if namelistLineContains(line,"nuPrime")
          line = "nuPrime = " + parametersForScan[i][8].to_s
        end
      when 21
        for nameind in 0..(parameternamesForScan.size-1)
          if namelistLineContains(line,parameternamesForScan[nameind])
            line = parameternamesForScan[nameind] + " = " + parametersForScan[i][nameind].to_s
          end
        end
      end
      
      outFile.write(line + "\n")
    end
    inFile.close
    outFile.close
    
    # Submit job!
    puts "Submitting job #{dirName}"
    # Make sure job 0 gets submitted first:
    #if i==0
    if true
      # In this way of submitting a job, ruby waits for qsub to complete before moving on.
      puts `cd #{dirName}; #{$jobSubmitCommand} #{jobFilename} &` #the ` signs make it happen!
    else
      # In this way of submitting a job, ruby does not wait for qsub to complete before moving on.
      job1 = fork do
        
        exec "cd #{dirName}; #{$jobSubmitCommand} #{jobFilename}"
      end
      Process.detach(job1)
    end
    puts "Done submitting job #{dirName}"
  end
  
end

puts "Finished submitting jobs for scan."
