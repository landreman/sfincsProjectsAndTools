#!/usr/bin/env python
import numpy as np
import sys
import bcgeom
import vmecgeom

class mnlist:

    def __init__(self,input,n=None,data=None,cosparity=1,Nperiods=None,
                 rind=None,quantity='B',vmecgrid='half'):
        #if isinstance(input,mnmat):
        # instead of this, use e.g. Bmn.mnlist()
            
        if isinstance(input,bcgeom.bcgeom):
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
            if quantity=='Dphi':
                self.data=data[rind]*2*np.pi/self.Nperiods
            else:
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
                    self.cosparity=geometry.parity[rind]
                elif quantity=='Z' or quantity=='Dphi':
                    self.cosparity=1-geometry.parity[rind]
                else:
                    sys.exit('Undefined quantity!')
            
        elif isinstance(input,vmecgeom.vmecgeom):
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
            
            self.Nperiods=wout.nfp

            if wout.StelSym:
                if quantity=='B':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange/self.Nperiods #signchange because of toroidal direction swap
                    self.data=wout.bmnc[skrind]
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='R':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange/self.Nperiods
                    self.data=(wout.rmnc[rindf_L]+wout.rmnc[rindf_R])/2
                    self.cosparity=np.ones(len(self.m))
                elif quantity=='Z':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange/self.Nperiods
                    self.data=(wout.zmns[rindf_L]+wout.zmns[rindf_R])/2
                    self.cosparity=np.zeros(len(self.m))
                elif quantity=='lambda':
                    self.m   =wout.xm
                    self.n   =wout.xn*signchange/self.Nperiods
                    self.data=(wout.lmns[rindf_L]+wout.lmns[rindf_R])/2
                    self.cosparity=np.zeros(len(self.m))
                elif quantity=='B_u':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange/self.Nperiods
                    self.data=wout.bsubumnc[skrind]
                    self.cosparity=np.ones(len(self.m))                    
                elif quantity=='B_w':
                    self.m   =wout.xm_nyq
                    self.n   =wout.xn_nyq*signchange/self.Nperiods
                    self.data=wout.bsubvmnc[skrind]*signchange
                    self.cosparity=np.ones(len(self.m))                    
                else:
                    sys.exit('Undefined quantity!')
            else:
                if quantity=='B':
                    self.m   =np.concatenate((wout.xm_nyq,wout.xm_nyq))
                    self.n   =np.concatenate((wout.xn_nyq,wout.xn_nyq))*signchange/self.Nperiods
                    self.data=np.concatenate((wout.bmnc[skrind],wout.bmns[skrind]))
                    self.cosparity=np.concatenate((np.ones(wout.mnmax_nyq),np.zeros(wout.mnmax_nyq)))
                elif quantity=='R':
                    self.m   =np.concatenate((wout.xm,wout.xm))
                    self.n   =np.concatenate((wout.xn,wout.xn))*signchange/self.Nperiods
                    self.data=np.concatenate(((wout.rmnc[rindf_L]+wout.rmnc[rindf_R])/2,
                                              (wout.rmns[rindf_L]+wout.rmns[rindf_R])/2))
                    self.cosparity=np.concatenate((np.ones(wout.mnmax),np.zeros(wout.mnmax)))
                elif quantity=='Z':
                    self.m   =np.concatenate((wout.xm,wout.xm))
                    self.n   =np.concatenate((wout.xn,wout.xn))*signchange/self.Nperiods
                    self.data=np.concatenate(((wout.zmns[rindf_L]+wout.zmns[rindf_R])/2,
                                              (wout.zmnc[rindf_L]+wout.zmnc[rindf_R])/2))                    
                    self.cosparity=np.concatenate((np.zeros(wout.mnmax),np.ones(wout.mnmax)))
                elif quantity=='lambda':
                    self.m   =np.concatenate((wout.xm,wout.xm))
                    self.n   =np.concatenate((wout.xn,wout.xn))*signchange/self.Nperiods
                    self.data=np.concatenate(((wout.lmns[rindf_L]+wout.lmns[rindf_R])/2,
                                              (wout.lmnc[rindf_L]+wout.lmnc[rindf_R])/2))
                    self.cosparity=np.concatenate((np.zeros(wout.mnmax),np.ones(wout.mnmax)))
                elif quantity=='B_u':
                    self.m   =np.concatenate((wout.xm_nyq,wout.xm_nyq))
                    self.n   =np.concatenate((wout.xn_nyq,wout.xn_nyq))*signchange/self.Nperiods
                    self.data=np.concatenate((wout.bsubumnc[skrind],wout.bsubumns[skrind]))
                    self.cosparity=np.concatenate((np.ones(wout.mnmax_nyq),np.zeros(wout.mnmax_nyq)))
                elif quantity=='B_w':
                    self.m   =np.concatenate((wout.xm_nyq,wout.xm_nyq))
                    self.n   =np.concatenate((wout.xn_nyq,wout.xn_nyq))*signchange/self.Nperiods
                    self.data=np.concatenate((wout.bsubvmnc[skrind],wout.bsubvmns[skrind]))*signchange
                    self.cosparity=np.concatenate((np.ones(wout.mnmax_nyq),np.zeros(wout.mnmax_nyq)))
                else:
                    sys.exit('Undefined quantity!')

        #if isinstance(input,mnmat.mnmat):
        #    return input.mnlist()
        
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
