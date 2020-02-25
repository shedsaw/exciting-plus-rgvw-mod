#!/bin/bash

# Compilers
export MAKE=make
export F90=mpifort
export CUDA_PATH=/usr/local/cuda
export NVCC=$CUDA_PATH/bin/nvcc

# Pick your desired compiler here
# Options: ibm (default), pgi, gcc, llvm
COMPILER=ibm

# Copy make.inc
cp make.inc.summit.${COMPILER}.cpu make.inc

# Switch compiler module
case ${COMPILER} in
    ibm)
	module swap xl
	;;
    pgi)
	module swap pgi
	;;
    gcc)
	module swap gcc
	;;
    llvm)
	module swap llvm
	;;
    *)
	echo "Unsupported compiler"
	exit 1
esac

# Load modules
module load essl
module load hdf5

# Clean up
make clean
rm *.o *.mod
rm src/elk-cpu bin/elk-cpu
#rm src/elk-gpu bin/elk-gpu

# Make the binary
make

# CPU-only version
mv src/elk src/elk-cpu

# Compile the necessary codes.
# $NVCC -c -g -G cublas_fortran.cu
# $F90 -cpp -c -g cublas_fortran_iso.f90
# $F90 -cpp -g -D_MPI_ -c -I./src/ genmegqblh_cublas.f90
# #$F90 -cpp -g -D_MPI_ -c -I./src/ genvscrn_cublas.f90

# Move the appropriate files over
# cp genmegqblh_cublas.o src/addons/expigqr/genmegqblh.o
# #cp genvscrn_cublas.o   src/addons/genvscrn.o
# cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/addons/expigqr/
# cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/addons/
# cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/

# Remove main and mod_mpi_grid so they will be recompiled
# rm src/main.o
# #rm src/mod_mpi_grid.mod # Keep old module for pstop()
# rm src/addons/mod_mpi_grid.o

# re-Make the binary.
# cp make.inc.local.gcc.gpu make.inc
# make

# Copy the hybrid CPU+GPU version
# mv src/elk src/elk-gpu

# Keep the two different versions
# make install # Feature isn't cherry-picked yet from cuda-hydra branch
[[ -d bin ]] || mkdir bin
[[ -e bin/elk ]] || rm bin/elk
cd bin
ln -s -T ../src/elk-cpu elk-cpu
#ln -s -T ../src/elk-gpu elk-gpu
cd ..
