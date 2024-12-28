#!/usr/bin/env bash

set -e

export BASEDIR=$HOME/ngsuite

mkdir -p $BASEDIR
cd $BASEDIR

git clone https://github.com/NGSolve/ngsolve.git ngsolve-src
cd $BASEDIR/ngsolve-src

git submodule update --init --recursive

mkdir $BASEDIR/ngsolve-build
mkdir $BASEDIR/ngsolve-install
cd $BASEDIR/ngsolve-build

cmake -DCMAKE_INSTALL_PREFIX=${BASEDIR}/ngsolve-install ${BASEDIR}/ngsolve-src
make -j4
make install

echo "export NETGENDIR=${BASEDIR}/ngsolve-install/bin" >> $HOME/.bashrc
echo "export PATH=\$NETGENDIR:\$PATH" >> $HOME/.bashrc

export PYTHONPATH_TMP=`python3 -c "import os.path, sysconfig;print(os.path.relpath(sysconfig.get_path('platlib'), sysconfig.get_path('data')))"`
echo "export PYTHONPATH=\$NETGENDIR/../${PYTHONPATH_TMP}:\$PATH" >> $HOME/.bashrc

rm -rf $BASEDIR/ngsolve-build
rm -rf $BASEDIR/ngsolve-src
