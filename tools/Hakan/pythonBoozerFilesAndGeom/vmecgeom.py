#!/usr/bin/env python
import numpy as np
import sys
from netCDF4 import Dataset
import timeit
#from ncdump import ncdump
#from mnlist import mnlist
#from mnmat import mnmat

class vmecgeom:
    
  def __init__(self,filename,min_Bmn=0,max_m=np.inf,maxabs_n=np.inf,
               symmetry='unknown'):
    # This reads the VMEC file of the wout type
    # Note that the wout file is left-handed and the output struct self
    # is right-handed.
    #
    # Only Fourier components with abs(Bmn)>min_Bmn, m<=max_m, |n|<=maxabs_n
    # are read.
    #
    # The input symmetry can be 'StelSym' to double-check that the 
    # input is stellarator symmetric
    #


    def ncdump(nc_fid, verb=False):
        '''
        ncdump outputs dimensions, variables and their attribute information.
        The information is similar to that of NCAR's ncdump utility.
        ncdump requires a valid instance of Dataset.

        Parameters
        ----------
        nc_fid : netCDF4.Dataset
            A netCDF4 dateset object
        verb : Boolean
            whether or not nc_attrs, nc_dims, and nc_vars are printed

        Returns
        -------
        nc_attrs : list
            A Python list of the NetCDF file global attributes
        nc_dims : list
            A Python list of the NetCDF file dimensions
        nc_vars : list
            A Python list of the NetCDF file variables
        '''
        def print_ncattr(key):
            """
            Prints the NetCDF file attributes for a given key

            Parameters
            ----------
            key : unicode
                a valid netCDF4.Dataset.variables key
            """
            try:
                print "\t\ttype:", repr(nc_fid.variables[key].dtype)
                for ncattr in nc_fid.variables[key].ncattrs():
                    print '\t\t%s:' % ncattr,\
                          repr(nc_fid.variables[key].getncattr(ncattr))
            except KeyError:
                print "\t\tWARNING: %s does not contain variable attributes" % key

        # NetCDF global attributes
        nc_attrs = nc_fid.ncattrs()
        if verb:
            print "NetCDF Global Attributes:"
            for nc_attr in nc_attrs:
                print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
        nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
        # Dimension shape information.
        if verb:
            print "NetCDF dimension information:"
            for dim in nc_dims:
                print "\tName:", dim 
                print "\t\tsize:", len(nc_fid.dimensions[dim])
                print_ncattr(dim)
        # Variable information.
        nc_vars = [var for var in nc_fid.variables]  # list of nc variables
        if verb:
            print "NetCDF variable information:"
            for var in nc_vars:
                if var not in nc_dims:
                    print '\tName:', var
                    print "\t\tdimensions:", nc_fid.variables[var].dimensions
                    print "\t\tsize:", nc_fid.variables[var].size
                    print_ncattr(var)
        return nc_attrs, nc_dims, nc_vars


    
    dataset=Dataset(filename)
    nc_attrs, nc_dims, nc_vars = ncdump(dataset)
    #print nc_vars
    
    self.StelSym=not('bmns' in dataset.variables)
    self.skip = 1 #=1,this is how many elements are skipped at low radii when going to half grid
    self.filename         =filename
    self.version_         =float(dataset.variables['version_'][:])
    self.input_extension  =''.join(dataset.variables['input_extension'][:])
    self.mgrid_file       =''.join(dataset.variables['mgrid_file'][:])
    self.pcurr_type       =''.join(dataset.variables['pcurr_type'][:])
    self.pmass_type       =''.join(dataset.variables['pmass_type'][:])
    self.piota_type       =''.join(dataset.variables['piota_type'][:])
    self.wb               =float(dataset.variables['wb'][:])
    self.wp               =float(dataset.variables['wp'][:])
    self.gamma            =float(dataset.variables['gamma'][:])
    self.rmax_surf        =float(dataset.variables['rmax_surf'][:])
    self.rmin_surf        =float(dataset.variables['rmin_surf'][:])
    self.zmax_surf        =float(dataset.variables['zmax_surf'][:])
    self.nfp              =float(dataset.variables['nfp'][:])
    self.ns               =int(dataset.variables['ns'][:])
    self.mpol             =int(dataset.variables['mpol'][:])
    self.ntor             =int(dataset.variables['ntor'][:])
    self.mnmax            =int(dataset.variables['mnmax'][:])
    self.mnmax_nyq        =int(dataset.variables['mnmax_nyq'][:])
    self.niter            =int(dataset.variables['niter'][:])
    self.itfsq            =float(dataset.variables['itfsq'][:])
    self.lasym__logical__ =bool(dataset.variables['lasym__logical__'][:])
    self.lrecon__logical__=bool(dataset.variables['lrecon__logical__'][:])
    self.lfreeb__logical__=bool(dataset.variables['lfreeb__logical__'][:])
    self.lrfp__logical__  =bool(dataset.variables['lrfp__logical__'][:])
    self.ier_flag         =int(dataset.variables['ier_flag'][:])
    self.aspect           =float(dataset.variables['aspect'][:])
    self.betatotal        =float(dataset.variables['betatotal'][:])
    self.betapol          =float(dataset.variables['betapol'][:])
    self.betator          =float(dataset.variables['betator'][:])
    self.betaxis          =float(dataset.variables['betaxis'][:])
    self.b0               =float(dataset.variables['b0'][:])
    self.rbtor0           =float(dataset.variables['rbtor0'][:])
    self.rbtor            =float(dataset.variables['rbtor'][:])
    self.signgs           =int(dataset.variables['signgs'][:])
    self.IonLarmor        =float(dataset.variables['IonLarmor'][:])
    self.volavgB          =float(dataset.variables['volavgB'][:])
    self.ctor             =float(dataset.variables['ctor'][:])
    self.Aminor_p         =float(dataset.variables['Aminor_p'][:])
    self.Rmajor_p         =float(dataset.variables['Rmajor_p'][:])
    self.volume_p         =float(dataset.variables['volume_p'][:])
    self.ftolv            =float(dataset.variables['ftolv'][:])
    self.fsql             =float(dataset.variables['fsql'][:])
    self.fsqr             =float(dataset.variables['fsqr'][:])
    self.fsqz             =float(dataset.variables['fsqz'][:])
    self.nextcur          =float(dataset.variables['nextcur'][:])
    self.extcur           =np.array(dataset.variables['extcur'][:]).astype(float)
    self.mgrid_mode       =''.join(dataset.variables['mgrid_mode'][:])
    self.xm               =np.array(dataset.variables['xm'][:]).astype(int)
    self.xn               =np.array(dataset.variables['xn'][:]).astype(int)
    self.xm_nyq           =np.array(dataset.variables['xm_nyq'][:]).astype(int)
    self.xn_nyq           =np.array(dataset.variables['xn_nyq'][:]).astype(int)
    self.raxis_cc         =np.array(dataset.variables['raxis_cc'][:]).astype(float)
    self.zaxis_cs         =np.array(dataset.variables['zaxis_cs'][:]).astype(float)
    if 'raxis_cs' in dataset.variables:
        self.raxis_cs         =np.array(dataset.variables['raxis_cs'][:]).astype(float)
        self.zaxis_cc         =np.array(dataset.variables['zaxis_cc'][:]).astype(float)
    else:
        self.raxis_cs=np.nan
        self.zaxis_cc=np.nan
    self.am               =np.array(dataset.variables['am'][:]).astype(float)
    self.ac               =np.array(dataset.variables['ac'][:]).astype(float)
    self.ai               =np.array(dataset.variables['ai'][:]).astype(float)
    self.am_aux_s         =np.array(dataset.variables['am_aux_s'][:]).astype(float)
    self.am_aux_f         =np.array(dataset.variables['am_aux_f'][:]).astype(float)
    self.ai_aux_s         =np.array(dataset.variables['ai_aux_s'][:]).astype(float)
    self.ai_aux_f         =np.array(dataset.variables['ai_aux_f'][:]).astype(float)
    self.ac_aux_s         =np.array(dataset.variables['ac_aux_s'][:]).astype(float)
    self.ac_aux_f         =np.array(dataset.variables['ac_aux_f'][:]).astype(float)
    self.iotaf            =np.array(dataset.variables['iotaf'][:]).astype(float)
    if 'q_factor' in dataset.variables:
      self.q_factor         =np.array(dataset.variables['q_factor'][:]).astype(float)
    else:
      self.q_factor         =np.nan*self.iotaf
    self.presf            =np.array(dataset.variables['presf'][:]).astype(float)
    self.phi              =np.array(dataset.variables['phi'][:]).astype(float)
    self.phipf            =np.array(dataset.variables['phipf'][:]).astype(float)
    if 'chi' in dataset.variables:
      self.chi              =np.array(dataset.variables['chi'][:]).astype(float)
      self.chipf            =np.array(dataset.variables['chipf'][:]).astype(float)
    else:
      self.chi              =np.nan*self.phi
      self.chipf            =np.nan*self.phi
    self.jcuru            =np.array(dataset.variables['jcuru'][:]).astype(float)
    self.jcurv            =np.array(dataset.variables['jcurv'][:]).astype(float)
    self.iotas            =np.array(dataset.variables['iotas'][:]).astype(float)
    self.mass             =np.array(dataset.variables['mass'][:]).astype(float)
    self.pres             =np.array(dataset.variables['pres'][:]).astype(float)
    self.beta_vol         =np.array(dataset.variables['beta_vol'][:]).astype(float)
    self.buco             =np.array(dataset.variables['buco'][:]).astype(float)
    self.bvco             =np.array(dataset.variables['bvco'][:]).astype(float)
    self.vp               =np.array(dataset.variables['vp'][:]).astype(float)
    self.specw            =np.array(dataset.variables['specw'][:]).astype(float)
    self.phips            =np.array(dataset.variables['phips'][:]).astype(float)
    self.over_r           =np.array(dataset.variables['over_r'][:]).astype(float)
    self.jdotb            =np.array(dataset.variables['jdotb'][:]).astype(float)
    self.bdotgradv        =np.array(dataset.variables['bdotgradv'][:]).astype(float)
    self.DMerc            =np.array(dataset.variables['DMerc'][:]).astype(float)
    self.DShear           =np.array(dataset.variables['DShear'][:]).astype(float)
    self.DWell            =np.array(dataset.variables['DWell'][:]).astype(float)
    self.DCurr            =np.array(dataset.variables['DCurr'][:]).astype(float)
    self.DGeod            =np.array(dataset.variables['DGeod'][:]).astype(float)
    self.equif            =np.array(dataset.variables['equif'][:]).astype(float)
    self.fsqt             =np.array(dataset.variables['fsqt'][:]).astype(float)
    self.wdot             =np.array(dataset.variables['wdot'][:]).astype(float)

    self.rmnc             =np.array(dataset.variables['rmnc'][:,:]).astype(float)
    self.zmns             =np.array(dataset.variables['zmns'][:,:]).astype(float)
    self.lmns             =np.array(dataset.variables['lmns'][:,:]).astype(float)
    self.gmnc             =np.array(dataset.variables['gmnc'][:,:]).astype(float)
    self.bmnc             =np.array(dataset.variables['bmnc'][:,:]).astype(float)
    self.bsubumnc         =np.array(dataset.variables['bsubumnc'][:,:]).astype(float)
    self.bsubvmnc         =np.array(dataset.variables['bsubvmnc'][:,:]).astype(float)
    self.bsubsmns         =np.array(dataset.variables['bsubsmns'][:,:]).astype(float)
    self.bsupumnc         =np.array(dataset.variables['bsupumnc'][:,:]).astype(float)
    self.bsupvmnc         =np.array(dataset.variables['bsupvmnc'][:,:]).astype(float)
    if self.StelSym:
      self.rmns             =np.nan
      self.zmnc             =np.nan
      self.lmnc             =np.nan
      self.gmns             =np.nan
      self.bmns             =np.nan
      self.bsubumns         =np.nan
      self.bsubvmns         =np.nan
      self.bsubsmnc         =np.nan
      self.bsupumns         =np.nan
      self.bsupvmns         =np.nan
    else:
      self.rmns             =np.array(dataset.variables['rmns'][:,:]).astype(float)
      self.zmnc             =np.array(dataset.variables['zmnc'][:,:]).astype(float)
      self.lmnc             =np.array(dataset.variables['lmnc'][:,:]).astype(float)
      self.gmns             =np.array(dataset.variables['gmns'][:,:]).astype(float)
      self.bmns             =np.array(dataset.variables['bmns'][:,:]).astype(float)
      self.bsubumns         =np.array(dataset.variables['bsubumns'][:,:]).astype(float)
      self.bsubvmns         =np.array(dataset.variables['bsubvmns'][:,:]).astype(float)
      self.bsubsmnc         =np.array(dataset.variables['bsubsmnc'][:,:]).astype(float)
      self.bsupumns         =np.array(dataset.variables['bsupumns'][:,:]).astype(float)
      self.bsupvmns         =np.array(dataset.variables['bsupvmns'][:,:]).astype(float)

    #sys.exit('tuuut')
    #self.tut=5

  ####################################################################
  
  def disp(self,verbose=False):
    frm='{:8.4f}'
    frme='{:10.4e}'
    frmi='{:8d}'
    frs='{:11s}'
    if verbose:
        print '-----------------------------------------------'
        print 'File info'
        print '-----------------------------------------------'
        print 'filename       = '+self.filename
        print 'version_       = '+str(self.version_)
        print 'input_extension= '+self.input_extension
        print 'mgrid_file     = '+self.mgrid_file
        print 'pcurr_type     = '+self.pcurr_type
        print 'pmass_type     = '+self.pmass_type
        print 'piota_type     = '+self.piota_type
        print 'mgrid_mode     = '+self.mgrid_mode
        print ' ' 
    print '-----------------------------------------------'
    print 'Scalars:'
    print '-----------------------------------------------'
    if verbose:
        print 'wb                   = '+frm.format(self.wb)+' :'
        print 'wp                   = '+frm.format(self.wp)+' :'
        print 'gamma                = '+frm.format(self.gamma)+' :'
        print 'rmax_surf            = '+frm.format(self.rmax_surf)+' :'
        print 'rmin_surf            = '+frm.format(self.rmin_surf)+' :'
        print 'zmax_surf            = '+frm.format(self.zmax_surf)+' :'
    print 'StelSym              =    '+str(self.StelSym)
    print 'nfp                  = '+frm.format(self.nfp)+' : Number of field periods'
    print 'ns                   = '+frmi.format(self.ns)+' : Number of flux surfaces'
    print 'mpol                 = '+frmi.format(self.mpol)+' :'
    print 'ntor                 = '+frmi.format(self.ntor)+' :'
    print 'mnmax                = '+frmi.format(self.mnmax)+' :'
    print 'mnmax_nyq            = '+frmi.format(self.mnmax_nyq)+' :'
    if verbose:
        print 'niter                = '+frmi.format(self.niter)+' :'
        print 'itfsq                = '+frm.format(self.itfsq)+' :'
        print 'lasym__logical__     =     '+str(self.lasym__logical__)+' :'
        print 'lrecon__logical__    =    '+str(self.lrecon__logical__)+' :'
        print 'lfreeb__logical__    =    '+str(self.lfreeb__logical__)+' :'
        print 'lrfp__logical__      =    '+str(self.lrfp__logical__)+' :'
        print 'ier_flag             = '+frmi.format(self.ier_flag)+' :'
        print 'aspect               = '+frm.format(self.aspect)+' :'
    print 'betatotal            = '+frm.format(self.betatotal)+' :'
    print 'betapol              = '+frm.format(self.betapol)+' :'
    print 'betator              = '+frm.format(self.betator)+' :'
    print 'betaxis              = '+frm.format(self.betaxis)+' :'
    print 'b0                   = '+frm.format(self.b0)+' :'
    if verbose:
        print 'rbtor0               = '+frm.format(self.rbtor0)+' :'
        print 'rbtor                = '+frm.format(self.rbtor)+' :'
        print 'signgs               = '+frmi.format(self.signgs)+' :'
        print 'IonLarmor            = '+frm.format(self.IonLarmor)+' :'
        print 'volavgB              = '+frm.format(self.volavgB)+' :'
        print 'ctor                 = '+frm.format(self.ctor)+' :'
    print 'Aminor_p             = '+frm.format(self.Aminor_p)+' : Minor radius'
    print 'Rmajor_p             = '+frm.format(self.Rmajor_p)+' : Major radius'
    print 'volume_p             = '+frm.format(self.volume_p)+' :'
    if verbose:
        print 'ftolv                = '+frm.format(self.ftolv)+' :'
        print 'fsql                 = '+frm.format(self.fsql)+' :'
        print 'fsqr                 = '+frm.format(self.fsqr)+' :'
        print 'fsqz                 = '+frm.format(self.fsqz)+' :'
        print 'nextcur              = '+frm.format(self.nextcur)+' :'

    print ' '
    if verbose:
        print '-----------------------------------------------'
        print 'Arrays:'
        print '-----------------------------------------------'
        print 'raxis_cc        '+frs.format(self.xm_nyq.shape)+  ': raxis (cosnv)'
        print 'zaxis_cs        '+frs.format(self.zaxis_cs.shape)+': zaxis (sinnv)'
        print 'raxis_cs        '+frs.format(self.raxis_cs.shape)+': raxis (sinnv)'
        print 'zaxis_cc        '+frs.format(self.zaxis_cc.shape)+': zaxis (cosnv)'
        print 'am              '+frs.format(self.am.shape)+':'
        print 'ac              '+frs.format(self.ac.shape)+':'
        print 'ai              '+frs.format(self.ai.shape)+':'
        print 'am_aux_s        '+frs.format(self.am_aux_s.shape)+':'
        print 'am_aux_f        '+frs.format(self.am_aux_f.shape)+':'
        print 'ai_aux_s        '+frs.format(self.ai_aux_s.shape)+':'
        print 'ai_aux_f        '+frs.format(self.ai_aux_f.shape)+':'
        print 'ac_aux_s        '+frs.format(self.ac_aux_s.shape)+':'
        print 'ac_aux_f        '+frs.format(self.ac_aux_f.shape)+':'
        print 'fsqt            '+frs.format(self.fsqt.shape)+':'
        print 'wdot            '+frs.format(self.wdot.shape)+':'
        print 'extcur          '+frs.format(self.extcur.shape)+':'
        
    print '-----------------------------------------------'
    print 'Radial arrays:'
    print '-----------------------------------------------'
    print 'iotaf        '+frs.format(self.iotaf.shape)+   ': q-factor on full mesh'
    if not((np.isnan(self.q_factor)).all()):
      print 'q_factor     '+frs.format(self.q_factor.shape)+':'
    print 'presf        '+frs.format(self.presf.shape)+': pressure on full mesh [Pa]'
    print 'phi          '+frs.format(self.phi.shape)+  ': Toroidal flux on full mesh [Wb]'
    print 'phipf        '+frs.format(self.phipf.shape)+': d(phi)/ds: Toroidal flux deriv on full mesh'
    if not((np.isnan(self.chi)).all()):
      print 'chi          '+frs.format(self.chi.shape)+  ': Poloidal flux on full mesh [Wb]'
      print 'chipf        '+frs.format(self.chipf.shape)+': d(chi)/ds: Poroidal flux deriv on full mesh'
    print 'jcuru        '+frs.format(self.jcuru.shape)+':'
    print 'jcurv        '+frs.format(self.jcurv.shape)+':'
    print 'iotas        '+frs.format(self.iotas.shape)+': iota half'
    print 'mass         '+frs.format(self.mass.shape)+ ': mass half'
    print 'pres         '+frs.format(self.pres.shape)+ ': pressure half [Pa]'
    print 'beta_vol     '+frs.format(self.beta_vol.shape)+':'
    print 'buco         '+frs.format(self.buco.shape)+':'
    print 'bvco         '+frs.format(self.bvco.shape)+':'
    print 'vp           '+frs.format(self.vp.shape)+':'
    print 'specw        '+frs.format(self.specw.shape)+':'
    print 'phips        '+frs.format(self.phips.shape)+':'
    print 'over_r       '+frs.format(self.over_r.shape)+':'
    print 'jdotb        '+frs.format(self.jdotb.shape)+':'
    print 'bdotgradv    '+frs.format(self.bdotgradv.shape)+':'
    if verbose:
        print 'DMerc        '+frs.format(self.DMerc.shape)+':'
        print 'DShear       '+frs.format(self.DShear.shape)+':'
        print 'DWell        '+frs.format(self.DWell.shape)+':'
        print 'DCurr        '+frs.format(self.DCurr.shape)+':'
        print 'DGeod        '+frs.format(self.DGeod.shape)+':'
        print 'equif        '+frs.format(self.equif.shape)+':'
    
    print ' '
    print '-----------------------------------------------'
    print 'Fourier quantities:'
    print '-----------------------------------------------'
    
    print 'xm        '+frs.format(self.xm.shape)+      ': Poloidal mode numbers'
    print 'xn        '+frs.format(self.xm.shape)+      ': Toroidal mode numbers'
    print 'xm_nyq    '+frs.format(self.xm_nyq.shape)+  ': Poloidal mode numbers (Nyquist)'
    print 'xn_nyq    '+frs.format(self.xm_nyq.shape)+    ': Toroidal mode numbers (Nyquist)'
    print 'rmnc      '+frs.format(self.rmnc.shape)+': cosmn component of cylindrical R, full mesh [m]'
    print 'zmns      '+frs.format(self.zmns.shape)+': sinmn component of cylindrical Z, full mesh [m]'
    print 'lmns      '+frs.format(self.lmns.shape)+': sinmn component of lambda, half mesh'
    print 'gmnc      '+frs.format(self.gmnc.shape)+': cosmn component of jacobian, half mesh'
    print 'bmnc      '+frs.format(self.bmnc.shape)+': cosmn component of mod-B, half mesh'
    print 'bsubumnc  '+frs.format(self.bsubumnc.shape)+': cosmn covariant u-component of B, half mesh'
    print 'bsubvmnc  '+frs.format(self.bsubvmnc.shape)+': cosmn covariant v-component of B, half mesh'
    print 'bsubsmns  '+frs.format(self.bsubsmns.shape)+': sinmn covariant s-component of B, full mesh'
    print 'bsupumnc  '+frs.format(self.bsupumnc.shape)+':'
    print 'bsupvmnc  '+frs.format(self.bsupvmnc.shape)+':'
    if not(self.StelSym):
      print 'rmns      '+frs.format(self.rmns.shape)+': sinmn component of cylindrical R, full mesh [m]'
      print 'zmnc      '+frs.format(self.zmnc.shape)+': cosmn component of cylindrical Z, full mesh [m]'
      print 'lmnc      '+frs.format(self.lmnc.shape)+': cosmn component of lambda, half mesh'
      print 'gmns      '+frs.format(self.gmns.shape)+': sinmn component of jacobian, half mesh'
      print 'bmns      '+frs.format(self.bmns.shape)+': sinmn component of mod-B, half mesh'
      print 'bsubumns  '+frs.format(self.bsubumns.shape)+': sinmn covariant u-component of B, half mesh'
      print 'bsubvmns  '+frs.format(self.bsubvmns.shape)+': sinmn covariant v-component of B, half mesh'
      print 'bsubsmnc  '+frs.format(self.bsubsmnc.shape)+': cosmn covariant s-component of B, full mesh'
      print 'bsupumns  '+frs.format(self.bsupumns.shape)+':'
      print 'bsupvmns  '+frs.format(self.bsupvmns.shape)+':'

    print '-----------------------------------------------'




    

    
 







      
   
