#!/usr/bin/env bash

set -e

#TODO: Remove after NGSolve image is rebuild
sudo apt-get install -y python3-pip
pip install scikit-build-core pybind11_stubgen numpy

#HACK: Passing NGSolve variables to NGBem
export NETGENDIR=$HOME/ngsuite/ngsolve-install/bin
export PATH=$NETGENDIR:$PATH

export PYTHONPATH_TMP=`python3 -c "import os.path, sysconfig;print(os.path.relpath(sysconfig.get_path('platlib'), sysconfig.get_path('data')))"`
export PYTHONPATH=$NETGENDIR/../${PYTHONPATH_TMP}:$PATH

export BEMBASEDIR=$HOME/ngbem

mkdir -p $BEMBASEDIR
cd $BEMBASEDIR

git clone https://github.com/Weggler/ngbem.git ngbem-src
cd $BEMBASEDIR/ngbem-src 

git submodule update --init --recursive

mkdir $BEMBASEDIR/ngbem-build
cd $BEMBASEDIR/ngbem-build

cmake ${BEMBASEDIR}/ngbem-src
make -j4
make install

rm -rf $BEMBASEDIR
