#!/bin/sh

#LATTEQEQPATH=/path/to/latteqeq/
LATTEQEQPATH=/vast/home/zhy/ecp/LATTE-git

 cmake  -DCMAKE_BUILD_TYPE=Release  -DPKG_LATTEQEQ=ON  -DCMAKE_PREFIX_PATH="$LATTEQEQPATH/install" -DPKG_ML-PACE=on -D PKG_MOLECULE=on -DCMAKE_CXX_STANDARD=11   ../cmake

 #copy lammps-user-pace package 
 cp -r ../tmp/lammps-user-pace-v.2022.09.27/ML-PACE/ace-evaluator/* lammps-user-pace-v.2022.09.27/ML-PACE/ace-evaluator/

 make -j


