Python scripts to evaluate classical particle flux given a SFINCS input file.

Dependencies:
• scipy (for special functions, numerical integration)
• numpy (for data sctructures)
• Håkan's stellarator geometry script
• (netcdf for Håkan's geometry scripts)
• h5py for reading SFINCS outputs
• os,sys for glue code

Note: Add sfincsProjectsAndTools/tools/Hakan/pythonBoozerFilesAndGeom to the PYTHONPATH environment variable. For example, in .bashrc, add
export PYTHONPATH=$PYTHONPATH:/path/to/sfincsProjectsAndTools/tools/Hakan/pythonBoozerFilesAndGeom

Comments on benchmarkning:
I've benchmarked against the classical impurity particle flux from an earlier MATLAB script, which assumed m_z/m_i >> 1 and Z_z>> 1. I was able to come withing 1% of the asymptotic values by increasing m_z and Z_z, but with high mass-ratios I have some numerical problems with the Chandrasekhar-integrals that prevents the results from improving much beyond m_z/m_i ~ 300. This should not strongly affect the results I have obtained, as the errors in these integrals is only in the 3rd significant digit even at mass-ratios comparable to m_C/m_e, but it is enough to not yield a better match to the asymptotic results. At yet higher mass-ratios, the problems become severe; improvements to the numerics are likely needed if one wants to look at tungsten transport, for example.

USAGE:
Calculate the classical flux from a SFINCS input.namelist. Prints to STDOUT:
./classical_flux.py <path to a SFINCS input.namelist> <signcorr>

Prints the neoclassical flux from a SFINCS simulation to STDOUT, given a SFINCS input.namelist:
./neoclassical_flux.py <path to a SFINCS input.namelist> <signcorr>

Calculate the ratio of classical/neoclassical and prints to STDOUT:
./nc_c_ratio.py  <path to a SFINCS input.namelist> <signcorr>

Note that the neoclassical_flux.py and nc_c_ratio.py scripts requires that the output specified in the SFINCS input.namelist exists and contains the neoclassical particle flux; the scripts do not actually calculate the neoclassical flux, but merely read it from an already existing SFINCS calculation. The ratio calculation does nothing clever, and is equivalent to running the other two scripts and taking the ratio.

TODOs/improvements:
Signcorr could probably be automatically calculated from the input.namelist.
