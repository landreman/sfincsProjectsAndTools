#!/usr/bin/env ruby

#By AM 2014-09#

#This ruby script launches a bunch of multiple sspecies fortran sfincs job to a batch system.
#Each step in the scan is specified in the file input.sfincsScan

# This ruby script launches a bunch of single species fortran sfincs jobs
# in a batch system, with the capability to perform several types of scans.

$sfincs_inputfile="input.namelist"

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

$scanInput="input.sfincsScan"

#The following line is only relevant on the hgw cluster:
clusters=["wsclus","wsclus","edge01","edge01"]
#clusters=["wsclus","wsclus","edge01","wsclus","edge02","wsclus","wsclus","edge03"]

require 'fileutils'
include FileUtils        
            
if !File.exists?($sfincs_inputfile)
  puts "Error! #{$sfincs_inputfile} not found."
  exit
end
puts "File #{$sfincs_inputfile} exists."

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
  s = `grep -i #{variableName} #{$sfincs_inputfile}`
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
    puts "Error! No lines in #{$sfincs_inputfile} match #{variableName}."
    exit
  end
  if numMatches > 1
    puts "Warning: more than 1 line in #{$sfincs_inputfile} matches #{variableName}."
  end
  if intOrFloat == 0
    puts "Read #{variableName} = " + value.to_i.to_s
    return value.to_i
  else
    puts "Read #{variableName} = " + value.to_s
    return value
  end
end

#def readScriptInput(ScriptInputfile, noElems)
#  if !File.exists?(ScriptInputfile)
#    puts "Error! #{ScriptInputfile} not found."
#    exit
#  end
#  puts "File #{ScriptInputfile} exists."
#  
#  inFile = File.open(ScriptInputfile,"r")
#  lines = inFile.readlines
#  inFile.close
#
#  radparams=Array.new(lines.size){Array.new(noElems)}
#  lind = 0
#  lines.each do |line|
#    if line.length>1
#      if line[0].chr != "!"
#        #puts line
#        substr=line
#        varind=0
#        while varind<noElems-1 do
#          varstr=substr[0..substr.index(" ")-1]
#          radparams[lind][varind]=varstr.to_f
#          substr=substr[(substr.index(" ")+1)..(substr.size-1)]
#          while substr.index(" ")==0
#            substr=substr[1..(substr.size-1)]
#          end
#          #puts "Element #{lind},#{varind} is #{radparams[lind][varind]}"
#          varind += 1
#        end
#        radparams[lind][noElems-1]=substr.to_f
#        #puts "Element #{lind},#{noElems-1} is #{radparams[lind][noElems-1]}"
#        lind += 1
#      end
#    end
#  end
#  return radparams[0..(lind-1)]
#end

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
##This script only runs for programMode=1
if (programMode != 1)
  puts "Error! Script only implemented for programMode=1."
  exit
end

#case programMode

  
#when 20
#  parametersForScan=readProfiles($scanInput, 9)
#  numRunsInScan = parametersForScan.size

#else
#  puts "I do not know what to do with programMode = "+programMode.to_s
#  exit
#end


##############################################################
# Make subdirs and write input.namelist in them, then launch.
##############################################################
#case programMode
#when 7,9,10

  
#else
  # A standard scan, like programMode = 2, 6, 8 or 20 but not 7, 9 or 10.
  
  # puts "Scan will consist of #{numRunsInScan} runs."
  
  
  #for i in 0..(numRunsInScan-1)
numRunsInScan = 0
dirNum = 0

inScriptFile = File.open($scanInput,"r")
scriptLines = inScriptFile.readlines
inScriptFile.close

#scriptLines.each_line do |scriptLine|
for lineNum in 0..(scriptLines.size-1)
  
  scriptLine = scriptLines[lineNum].strip
  
  if (scriptLine[0].chr == "!" or scriptLine.empty?)
    next #Jump to next iteration
  end

  parametersToSet = scriptLine.split(',')

  for k in 0..(parametersToSet.size-1)
    puts parametersToSet[k].strip
  end
  exit
  
  # Starting with 0, look for the first unused whole number to use for a directory name:
  begin
    dirNum += 1
    dirName = dirNum.to_s
    #if dirNum < 10
    #  dirName = "0" + dirName
    #end
  end until !File.exists?(dirName)
  mkdir(dirName)
  
  jobName="MultiSfincsScan.#{dirName}"
  
  
  # Copy sge file (pbs file in Matt's case)
  outFilename = dirName + "/" + jobFilename
  inFile = File.open(jobFilename,"r")
  outFile = File.open(outFilename,"w")
  lines = inFile.readlines
  outFile.write(lines[0])
  case $sfincs_system 
  when "hgw"
    outFile.write("#\$ -N #{jobName}\n")
    for j in 1..(lines.size-1)
      if lines[j].include? "CLUSTR"
        clind=lines[j].index('CLUSTR')
        lines[j][clind..clind+5]=clusters[i.modulo(clusters.size)]
        puts lines[j]
      end 
      outFile.write(lines[j])
    end
  when "hydra"
    #outFile.write("# @ error = #{jobName}.e$(jobid)\n")
    #outFile.write("# @ output = #{jobName}.o$(jobid)\n")
    outFile.write("# @ job_name = #{jobName}.$(jobid)\n")
    for j in 1..(lines.size-1)
      outFile.write(lines[j])
    end
  end

  inFile.close
  outFile.close
    
  # Copy input.namelist
  outFilename = dirName + "/" + $sfincs_inputfile
  inFile = File.open($sfincs_inputfile,"r")
  outFile = File.open(outFilename,"w")
  lines = inFile.readlines

  for j in 0..(lines.size-1)
    line = lines[j].strip
      
    #if namelistLineContains(line,"programMode")
    #  line = "programMode = 1"
    #end
      

    if namelistLineContains(line,"normradius_wish")
      line = "normradius_wish = " + parametersForScan[i][0].to_s
    end
    
    if namelistLineContains(line,"THat")
      line = "THat = " + parametersForScan[i][1].to_s
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

#end

puts "Finished submitting jobs for scan."
