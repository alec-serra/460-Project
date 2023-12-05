
rm binary_brute
gfortran -o tree tree.f90 -Wall -Wno-unused-function -finit-local-zero -fno-automatic -g -fbacktrace -fcheck=all 
./tree
python3 tree.py
python3 compare.py
