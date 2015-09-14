#!/bin/bash
wget http://gnu.uberglobalmirror.com/gsl/gsl-1.16.tar.gz
mkdir -p gsl
tar xvfz gsl-1.16.tar.gz
cd gsl-1.16
make clean
./configure --disable-shared --disable-dependency-tracking
make
sudo make install
sudo apt-get install pkg-config
