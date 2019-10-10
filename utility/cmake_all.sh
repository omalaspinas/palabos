#!/bin/bash
set -e

loc_dir=`pwd`
echo ${loc_dir}

for dir in `find . -name "build"  -type d`
do
    cd ${loc_dir}/${dir}
    cmake ..
done

cd ${loc_dir}

for dir in `find . -name "build"  -type d`
do
    cd ${loc_dir}/${dir}
    make -j4
done

cd ${loc_dir}
