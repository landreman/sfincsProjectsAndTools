#!/usr/bin/env python
import numpy, sys
from bcgeom import bcgeom

class mnlist:

    def __init__(self,input,n=None,data=None,cosparity=1,Nperiods=None,
                 rind=None,quantity='B'):
        if isinstance(input,bcgeom):
            geometry=input
            if rind is None:
                sys.exit('The radius index rind is required when taking input from a bcgeom!')

            self.Nperiods=geometry.Nperiods
            self.m=geometry.m[rind]
            self.n=geometry.n[rind]
            data=getattr(geometry,quantity)
            self.data=data[rind]
            if geometry.StelSym:
                if quantity=='B' or quantity=='R':
                    self.cosparity=numpy.ones(len(self.m))
                elif quantity=='Z' or quantity=='Dphi':
                    self.cosparity=numpy.zeros(len(self.m))
                else:
                    sys.exit('Undefined quantity!')
            else:
                if quantity=='B' or quantity=='R':
                    self.cosparity=geometry.parity
                elif quantity=='Z' or quantity=='Dphi':
                    self.cosparity=1-geometry.parity
                else:
                    sys.exit('Undefined quantity!')
            
        else: #input is supposed to be the m array
            m=input
            self.Nperiods=Nperiods
            if len(m)!=len(n):
                sys.exit("m and n have different lengths")
            if len(n)!=len(data):
                sys.exit("n and data have different lengths")
            if isinstance(cosparity,int):
                cosparity=float(cosparity)
            if isinstance(cosparity,float):
                self.cosparity=numpy.ones(len(m))*cosparity
            elif len(cosparity)!=len(m):
                sys.exit("m and cosparity have different lengths")
            else:
                self.cosparity=numpy.array(cosparity)
            self.m=numpy.array(m)
            self.n=numpy.array(n)
            self.data=numpy.array(data)

                
    def disp(self):
        print('---------')
        if not(self.Nperiods is None):
            print 'Nperiods='+str(self.Nperiods)
        print 'cosparity=\n'+str(self.cosparity)
        print 'm=\n'+str(self.m)
        print 'n=\n'+str(self.n)
        print 'data=\n'+str(self.data)        
        print('---------')
