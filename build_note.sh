#!/bin/bash

# Install instructions for lammps with flare
# this assumes that you do NOT have mkl libraries and gcc/gxx compilers 
# loaded. If loaded, you may skip the 'ONLY IF NOT LOADED' section below

mamba create --name flare python=3.10
mamba activate flare
mamba install numpy scipy scikit-learn virtualenv psutil pandas tabulate mpi4py Cython
mamba  install conda-forge::libjpeg
mamba  install conda-forge::libpng
mamba  install conda-forge::zlib
# ONLY IF NOT LOADED/AVAILABLE
mamba install conda-forge::gcc
mamba install conda-forge::mkl
mamba install anaconda::mkl-include
mamba  install conda-forge::cxx-compiler

# continue installation instructions here:

# download lammps
git clone -b flare_update https://github.com/jmgoff/lammps_compute_PACE.git
cd lammps_compute_PACE
mkdir build && cd build
cmake ../cmake -DBUILD_MPI=OFF -DBUILD_SHARED_LIBS=ON -DPKG_PYTHON=ON -DPKG_MANYBODY=ON && make -j$(nproc)
make install-python

# navigate to folder above lammps
# download flare
git clone --depth=1 https://github.com/mir-group/flare

cd flare
pip install .
pip install --force-reinstall ase==3.23.0

# after this, try the example in lammps_compute_PACE/flare_lammps_active
