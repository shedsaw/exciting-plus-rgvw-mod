#!/bin/bash

# Compilers
export MAKE=make
export F90=$TAU/x86_64/bin/tau_f90.sh
export CUDA_PATH=/usr/local/cuda
export NVCC=$CUDA_PATH/bin/nvcc

# TAU environment variables
export TAU_MAKEFILE=$TAU/x86_64/lib/Makefile.tau-pgimpicuda-papi-mpi-pdt-pgi
export TAU_PROFILE=1
export TAU_CALLPATH=1
export TAU_OPTIONS="-optPDTInst -optRevert -optPreProcess -optShared -optLinking=$(LIBS)"
export TAU_METRICS=TIME

# Clean up
make clean
rm *.o *.mod
rm bin/elk-cpu bin/elk-gpu

# Make the binary
cp make.inc.local.tau-pgi.cpu make.inc
make

# CPU-only version
mv src/elk src/elk-cpu

# Compile the necessary codes.
$NVCC -c -g -G cublas_fortran.cu
$F90 -Mpreprocess -c -g cublas_fortran_iso.f90
$F90 -Mpreprocess -g -D_MPI_ -c -I./src/ genmegqblh_cublas.f90
#$F90 -Mpreprocess -g -D_MPI_ -c -I./src/ genvscrn_cublas.f90

# Move the appropriate files over
cp genmegqblh_cublas.o src/addons/expigqr/genmegqblh.o
#cp genvscrn_cublas.o   src/addons/genvscrn.o
cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/addons/expigqr/
cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/addons/
cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/

# Remove main and mod_mpi_grid so they will be recompiled
rm src/main.o
#rm src/mod_mpi_grid.mod # Keep old module for pstop()
rm src/addons/mod_mpi_grid.o

# re-Make the binary.
cp make.inc.local.tau-pgi.gpu make.inc
make

# Copy the hybrid CPU+GPU version
mv src/elk src/elk-gpu

# Keep the two different versions
make install
rm bin/elk
cd bin
ln -s -T ../src/elk-cpu elk-cpu
ln -s -T ../src/elk-gpu elk-gpu
cd ..

