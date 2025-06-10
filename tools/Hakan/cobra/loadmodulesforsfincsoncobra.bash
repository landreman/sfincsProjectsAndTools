#module load intel/19.1.2 openmpi/4 #this is necessary to be allowed to load the following modules
#module load hdf5-serial
#module load netcdf-serial
#module load petsc-real
#module load petsc-openmp-real/3.9.3
#module load intel/19.1.2 impi/2019.8 hdf5-mpi/1.8.21 netcdf-mpi/4.4.1 petsc-real/3.13.5

source ~/sfincs/sfincsProjectsAndTools/tools/Hakan/cobra/loadmodulesforsfincsoncobra.common

#export PETSC_ARCH=20170728-10-removingMPI_shouldWork
#export PETSC_DIR=/draco/u/mlan/petsc-3.7.6
#export PETSC_HOME=/draco/u/mlan/petsc-3.7.6

export PATH=${PATH}:${HDF5_HOME}/bin  #Or "export" if using bash 
export SFINCS_SYSTEM=cobra
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_HOME}/lib:${NETCDF_HOME}/lib
export MUMPS_OOC_TMPDIR=/ptmp/smithh/mumps_temp


