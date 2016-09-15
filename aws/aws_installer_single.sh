#!/bin/bash

# Exit if something bad happens.
# See http://stackoverflow.com/questions/821396/aborting-a-shell-script-if-any-command-returns-a-non-zero-value
set -eu

sudo yum -y update
sudo yum -y install gcc gcc-c++ openmpi openmpi-devel boost boost-devel scons git

# for ubuntu 13.04 ami for StarCluster
# sudo sed -i -e 's|//.*ubuntu.com/|//old-releases.ubuntu.com/|' /etc/apt/sources.list
# sudo apt-get update
# sudo apt-get install libboost1.53-dev scons libgflags2 libgflags-dev
# # edit SConstruct (for link order)
# # for selecting openmpi (by default, mpiexec is from openmpi and library is mpich)
# sudo update-alternatives --config  mpi

if [ ${PATH:-default} = "default" ]
then
    export PATH=/usr/lib64/openmpi/bin
else
    export PATH=/usr/lib64/openmpi/bin:${PATH}
fi

if [ ${LD_LIBRARY_PATH:-default} = "default" ]
then
    export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:/usr/local/lib
else
    export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:/usr/local/lib:${LD_LIBRARY_PATH}
fi

cat >> ${HOME}/.bashrc <<EOF

if [ ${PATH:-default} = "default" ]
then
    export PATH=/usr/lib64/openmpi/bin
else
    export PATH=/usr/lib64/openmpi/bin:${PATH}
fi

if [ ${LD_LIBRARY_PATH:-default} = "default" ]
then
    export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:/usr/local/lib
else
    export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:/usr/local/lib:${LD_LIBRARY_PATH}
fi
EOF

# gflags
wget https://github.com/gflags/gflags/archive/v2.0.zip
unzip v2.0.zip
cd gflags-2.0/
./configure
make
sudo make install
cd ..

# mp-lamp
# wget https://dl.dropboxusercontent.com/u/1576913/mp-lamp/mp-lamp.tgz
# tar xvzf mp-lamp.tgz
# cd mp-lamp

scons -u

mpiexec -np 2 mp-lamp --item samples/sample_data/sample_item.csv --pos samples/sample_data/sample_expression_over1.csv --a 0.05

echo do \"source \$\{HOME\}/.bashrc\"
