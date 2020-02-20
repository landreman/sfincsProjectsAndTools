#!/usr/bin/env python
from __future__ import division
import sys, copy, os, argparse, inspect, math, subprocess
import numpy as np
import datetime

import mnFourierlib
import fluxcoorddiscr
import netCDF4
from netCDF4 import Dataset
#import timeit
import scipy.integrate as integrate

#
# Note that ./blackbox cannot read the output from this,
# because m0b and n0b do not reflect the content anymore
#

###############################################################
#def parse_args(argv):
#    parser = argparse.ArgumentParser(description='bcfile filtering')
#    parser.add_argument('nickname', type=str)
#    return parser.parse_known_args(argv[1:])
#
###############################################################
def sscan(strng,dtype=float,count=-1):
    return np.fromstring(strng,dtype,count,sep=' ')

def filterbigbcfile(infile,outfile,min_Bmn,max_m=np.inf,maxabs_n=np.inf):

  if infile==outfile:
      sys.exit('The input and output files are the same!')
        
  fin = open(infile, 'r')

  if sys.version_info[0] == 2:
    fin.seek(-1,2)     # go to the file end.
  else:
    fin.seek(0, os.SEEK_END)
    if fin.tell()==0:
      fin.close()
      sys.exit('Input file was empty!')    
    fin.seek(fin.tell() - 1, os.SEEK_SET)
  eof = fin.tell()   # get the end of file location
  fin.seek(0,0)      # go back to file beginning

  fout= open(outfile, 'w')
  
  tmp_str=fin.readline()
  if tmp_str[0:2]=='CC':
    while tmp_str[0:2]=='CC':
      fout.write(tmp_str)
      tmp_str=fin.readline()

  fout.write(tmp_str) #variable name line
  fout.write(fin.readline()) #variables
  headertext_surfvars=fin.readline() #Variable name line  
  YTstyle=0
  if '[A]' in headertext_surfvars:
    #This is a file from Yuriy Turkin. Correct this line to the JG standard
    YTstyle=1 #indicates Yuriy Turkin style
    headertext_surfvars=(
               '       s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)')
    headertext_surfvarunits=(
               '                            [A]            [A]   dp/ds,[Pa] (dV/ds)/nper')

  endoffile=False
  rind=-1

  while not endoffile:
    rind=rind+1
    if not(YTstyle): #isempty(strfind(self.headertext.surfvars,'[A]'))
      headertext_surfvarunits=fin.readline() #unit line only in JG files

    surfheader=fin.readline()  
    fout.write(headertext_surfvars)
    fout.write(headertext_surfvarunits)
    fout.write(surfheader)
    tmpstrunits=fin.readline() #units line
    fout.write(tmpstrunits)

    if rind==0:
      StelSym=True
      if 'rmnc' in tmpstrunits:
        StelSym=False 

    position=fin.tell()

    if StelSym:
      tmp=sscan(fin.readline()) #,'%d %d %f %f %f %f',6)
      while not(tmp[0]==0 and tmp[1]==0):
        tmp=sscan(fin.readline(),count=6) #sscanf(tmp_str1,'%d %d %f %f %f %f',6)    
      B00=tmp[5]
    else:
      tmp=sscan(fin.readline())#sscanf(tmp_str1,'%d %d %f %f %f %f %f %f %f %f',10)
      while not(tmp[0]==0 and tmp[1]==0):
        tmp=sscan(fin.readline()) #sscanf(tmp_str1,'%d %d %f %f %f %f %f %f %f %f',10)
      B00=tmp[8]

    fin.seek(position) #Rewind to beginning of surface data

    proceed=True
    while proceed:
        tmp_str=fin.readline()
        if fin.tell() == eof+1:
            #print('eof found '+tmp_str)
            proceed=False
            endoffile=True
            #print('found end of file')
        if ('s' in tmp_str): #Next flux surface has been reached
            proceed=False
            #print('found next surface')
        else:
            #print(tmp_str)
            if StelSym:
                tmp=sscan(tmp_str,count=6)
                if ((abs(tmp[5])/B00>min_Bmn) and
                    (tmp[0]<=max_m) and (abs(tmp[1])<=maxabs_n)):
                    fout.write(tmp_str)

            else:
                tmp=sscan(tmp_str,count=10)
                if (tmp[0]<=max_m) and (abs(tmp[1])<=maxabs_n):
                    if (abs(tmp[8])/B00>min_Bmn) or (abs(tmp[9])/B00>min_Bmn):
                        fout.write(tmp_str)
                    #end
                #end
            #end
        #end
    #end while proceed

    
  fin.close()
  fout.close()
