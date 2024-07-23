#!/bin/bash

cd extern/lemon-1.3.1/
rm -rf build
mkdir build
cd build
cmake ../
make 
sudo make install
cd ../../../
