#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
set -x

rootpath=${PWD}
export GSL_ROOT_DIR=${rootpath}/gsl-2.6/gsl
# 1. build gsl
build_gsl=$rootpath/gsl-2.6/src/build
mkdir $build_gsl
cd $build_gsl
../configure --prefix=${GSL_ROOT_DIR}
make -j 8 #using 8 threads
make install

# 2. build swEOS
build_swEOS=$rootpath/build
mkdir $build_swEOS
cd $build_swEOS
# Detect architecture (sysctl works correctly even under Rosetta)
if [ "$(sysctl -n hw.optional.arm64 2>/dev/null)" = "1" ]; then
  HOST_ARCH=arm64
else
  HOST_ARCH=$(uname -m)
fi
cmake -DCMAKE_OSX_ARCHITECTURES=${HOST_ARCH} ..
make install
