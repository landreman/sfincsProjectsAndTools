module load netcdf-serial
module load hdf5-serial
#module load petsc-real
#module load petsc-openmp-real/3.9.3

export PETSC_ARCH=20170728-10-removingMPI_shouldWork
export PETSC_DIR=/draco/u/mlan/petsc-3.7.6
export PETSC_HOME=/draco/u/mlan/petsc-3.7.6

export PATH=${PATH}:${HDF5_HOME}/bin  #Or "export" if using bash 
export SFINCS_SYSTEM=draco 
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_HOME}/lib:${NETCDF_HOME}/lib
export MUMPS_OOC_TMPDIR=/ptmp/smithh/mumps_temp


