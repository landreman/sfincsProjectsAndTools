module switch intel/14.0
module unload mpi.ibm
module load mpi.intel/4.1.3
module switch mkl/11.1

module load petsc-real/3.5.2

# We must load the hdf5-serial module to set some environment variables:
module load hdf5-serial netcdf-serial

#The next line is not necessary according to Albert
#export HDF5_HOME=/hydra/u/system/SLES11/soft/hdf5/1.8.11/intel13.1/mpi.intel-4.1.0

export PATH=${PATH}:${HDF5_HOME}/bin

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_HOME}/lib:${NETCDF_HOME}/lib

export MUMPS_OOC_TMPDIR=/ptmp/smithh/mumps_temp