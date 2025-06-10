
source ~/sfincs/sfincsProjectsAndTools/tools/Hakan/raven/loadmodulesforsfincsonraven.common

setenv PATH ${PATH}:${HDF5_HOME}/bin  
setenv SFINCS_SYSTEM raven 
setenv LD_LIBRARY_PATH ${HDF5_HOME}/lib:${NETCDF_HOME}/lib
setenv MUMPS_OOC_TMPDIR /ptmp/smithh/mumps_temp

