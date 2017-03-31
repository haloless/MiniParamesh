#! /bin/sh


#
echo 'clean lib and app'

make -f make_app clean

#
echo 'build lib and app'
make -f make_app your_app

# 
echo 'test app'
cd your_app
mpirun -np 1 tutor > tutorial_output
cd ..



