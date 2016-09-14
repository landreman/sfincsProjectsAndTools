module load netcdf-serial
module load hdf5-serial
module load petsc-real
export PATH=${PATH}:${HDF5_HOME}/bin  #Or "export" if using bash 
export SFINCS_SYSTEM=draco 
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_HOME}/lib:${NETCDF_HOME}/lib
export MUMPS_OOC_TMPDIR=/ptmp/smithh/mumps_temp


