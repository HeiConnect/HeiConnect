#!/bin/bash
if (( $# != 1 )); then
    >&2 echo "Usage: $0 <buildtype:Release/Debug> "
    buildtype=Release
fi

buildtype=$1 # Release or Debug 

git submodule update --init --recursive

NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi

# compile 
rm -rf deploy
rm -rf build
mkdir build
cd build
cmake ../ -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_BUILD_TYPE=${buildtype}
make -j $NCORES
cd ..

mkdir deploy
cp ./build/solver                   deploy/solver
cp ./build/graphchecker             deploy/graphchecker
rm -rf build


