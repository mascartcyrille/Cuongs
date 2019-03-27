#!/bin/sh
NUM_VERSION=9.300.2
DESTDIR=/tmp/my_install/armadillo-$NUM_VERSION
cd /tmp
wget http://sourceforge.net/projects/arma/files/armadillo-$NUM_VERSION.tar.xz
tar xf armadillo-$NUM_VERSION.tar.xz 
cd armadillo-$NUM_VERSION
cmake .
make
make install DESTDIR=$DESTDIR
