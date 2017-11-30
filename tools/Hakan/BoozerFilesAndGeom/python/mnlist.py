#!/usr/bin/env python
import numpy as np
import sys
from bcgeom import bcgeom

class mnlist:

    def __init__(self,input,n=None,data=None,cosparity=1,Nperiods=None,
                 rind=None,quantity='B',vmecgrid='half'):
        if isinstance(input,bcgeom):
            #quantity can be R,Z,B,Dphi
            geometry=input
            if rind is None:
                sys.exit('The radius index rind is required when taking input from a bcgeom!')

            if Nperiods != None:
                if geometry.Nperiods != Nperiods:
                    sys.exit('Input Nperiods does not match the value in the bcgeom!')
                    
            self.Nperiods=geometry.Nperiods
            self.m=geometry.m[rind]
            self.n=geometry.n[rind]
            data=getattr(geometry,quantity)
            self.data=data[rind]
            if geometry.StelSym:
                if quantity=='B' or quantity=='R':
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='Z' or quantity=='Dphi':
                    self.cosparity=np.zeros(len(self.m))
                else:
                    sys.exit('Undefined quantity!')
            else:
                if quantity=='B' or quantity=='R':
                    self.cosparity=geometry.parity
                elif quantity=='Z' or quantity=='Dphi':
                    self.cosparity=1-geometry.parity
                else:
                    sys.exit('Undefined quantity!')
            
        elif isinstance(input,vmecgeom):
            #quantity can be R,Z,B,lambda,B_u or B_w where (u,w,s) is RH and w=-v. (u,v,s) is the LH VMEC system.
            wout=input
            if rind is None:
                sys.exit('The radius index rind is required when taking input from a bcgeom!')

            if Nperiods != None:
                if wout.Nperiods != Nperiods:
                    sys.exit('Input Nperiods does not match the value in the vmecgeom!')
                    
            signchange=float(wout.signgs) #is -1, because vmec is left handed
            if rind is None:
                sys.exit('The radius index rind is required when taking input from a vmecgeom!')
            if vmecgrid=='full':
                sys.exit('Extracting from full vmec grid is not implemented yet!')
            skip=wout.skip #=1, this is the skip in the beginning of the half grid vmec arrays
            skrind=skip+rind
            rindf_R = skrind
            rindf_L = skrind-1
            
            self.Nperiods=wout.Nperiods

            if wout.StelSym:
                if quantity=='B':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange/wout.Nperiods
                    self.data=wout.bmnc[skrind]
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='R':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange/wout.Nperiods
                    self.data=(wout.rmnc[rindf_L]+wout.rmnc[rindf_R])/2
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='Z':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange/wout.Nperiods
                    self.data=(wout.zmns[rindf_L]+wout.zmns[rindf_R])/2
                    self.cosparity=np.zeros(len(self.m))
                elif quantity=='lambda':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange/wout.Nperiods
                    self.data=(wout.lmns[rindf_L]+wout.lmns[rindf_R])/2
                    self.cosparity=np.zeros(len(self.m))
                elif quantity=='B_u':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange/wout.Nperiods
                    self.data=wout.bsubumnc[skrind]
                    self.cosparity=np.ones(len(self.m))                    
                elif quantity=='B_w':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange/wout.Nperiods
                    self.data=wout.bsubvmnc[skrind]*signchange
                    self.cosparity=np.ones(len(self.m))                    
                else:
                    sys.exit('Undefined quantity!')
            else:
                if quantity=='B':
                    self.m   =np.array([]).append(wout.xm_nyq).append(wout.xm_nyq)
                    self.n   =np.array([]).append(wout.xn_nyq).append(wout.xn_nyq)*signchange/wout.Nperiods
                    self.data=np.array([]).append(wout.bmnc[skrind]).append(wout.bmns[skrind])
                    self.cosparity=np.ones(len(self.m)).append(np.zeros(len(self.m)))
                elif quantity=='R':
                    self.m   =np.array([]).append(wout.xm).append(wout.xm)
                    self.n   =np.array([]).append(wout.xn).append(wout.xn)*signchange/wout.Nperiods
                    self.data=np.array([]).append((wout.rmnc[rindf_L]+wout.rmnc[rindf_R])/2).append(
                        (wout.rmns[rindf_L]+wout.rmns[rindf_R])/2)
                    self.cosparity=np.ones(len(self.m)).append(np.zeros(len(self.m)))
                elif quantity=='Z':
                    self.m   =np.array([]).append(wout.xm).append(wout.xm)
                    self.n   =np.array([]).append(wout.xn).append(wout.xn)*signchange/wout.Nperiods
                    self.data=np.array([]).append((wout.zmns[rindf_L]+wout.zmns[rindf_R])/2).append(
                        (wout.zmnc[rindf_L]+wout.zmnc[rindf_R])/2)                    
                    self.cosparity=np.zeros(len(self.m)).append(np.ones(len(self.m)))
                elif quantity=='lambda':
                    self.m   =np.array([]).append(wout.xm).append(wout.xm)
                    self.n   =np.array([]).append(wout.xn).append(wout.xn)*signchange/wout.Nperiods
                    self.data=np.array([]).append((wout.lmns[rindf_L]+wout.lmns[rindf_R])/2).append(
                        (wout.lmnc[rindf_L]+wout.lmnc[rindf_R])/2)                    
                    self.cosparity=np.zeros(len(self.m)).append(np.ones(len(self.m)))
                elif quantity=='B_u':
                    self.m   =np.array([]).append(wout.xm_nyq).append(wout.xm_nyq)
                    self.n   =np.array([]).append(wout.xn_nyq).append(wout.xn_nyq)*signchange/wout.Nperiods
                    self.data=np.array([]).append(wout.bsubumnc[skrind]).append(wout.bsubumns[skrind])
                    self.cosparity=np.ones(len(self.m)).append(np.zeros(len(self.m)))
                elif quantity=='B_w':
                    self.m   =np.array([]).append(wout.xm_nyq).append(wout.xm_nyq)
                    self.n   =np.array([]).append(wout.xn_nyq).append(wout.xn_nyq)*signchange/wout.Nperiods
                    self.data=np.array([]).append(wout.bsubvmnc[skrind]).append(wout.bsubvmns[skrind])*signchange
                    self.cosparity=np.ones(len(self.m)).append(np.zeros(len(self.m)))
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
                self.cosparity=np.ones(len(m))*cosparity
            elif len(cosparity)!=len(m):
                sys.exit("m and cosparity have different lengths")
            else:
                self.cosparity=np.array(cosparity)
            self.m=np.array(m)
            self.n=np.array(n)
            self.data=np.array(data)

                
    def disp(self):
        print('---------')
        if not(self.Nperiods is None):
            print 'Nperiods='+str(self.Nperiods)
        print 'cosparity=\n'+str(self.cosparity)
        print 'm=\n'+str(self.m)
        print 'n=\n'+str(self.n)
        print 'data=\n'+str(self.data)        
        print('---------')
