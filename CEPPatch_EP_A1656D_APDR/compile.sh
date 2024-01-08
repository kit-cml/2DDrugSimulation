!#/bin/sh

make clean all
cd bin
./patch_numeric_ORd_$1 $2
cd ..
