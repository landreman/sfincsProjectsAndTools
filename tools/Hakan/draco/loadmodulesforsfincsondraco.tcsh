module load netcdf-serial
module load hdf5-serial
module load petsc-real
setenv PATH ${PATH}:${HDF5_HOME}/bin  
setenv SFINCS_SYSTEM draco 
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HDF5_HOME}/lib:${NETCDF_HOME}/lib
setenv MUMPS_OOC_TMPDIR /ptmp/smithh/mumps_temp
