#!/bin/bash

# Compilers
export MAKE=make
export F90=mpifort
export CUDA_PATH=/usr/local/cuda
export NVCC=$CUDA_PATH/bin/nvcc

# Pick your desired compiler here
# Options: ibm, pgi, gcc, tau-ibm, tau-pgi, tau-gcc
# TODO: read this as argv[1]
COMPILER=ibm

# Load ESSL and HDF5 modules
module load essl
module load hdf5

# ESSL depends on libxlf90_r
# This function extracts IBM XL compiler paths and saves them to xlvars.sh
getxlvars() {
  module load xl
  cat > summit-xlvars.sh << __EOF__
export OLCF_XL_ROOT=${OLCF_XL_ROOT}
export OLCF_XLF_ROOT=${OLCF_XLF_ROOT}
export OLCF_XLC_ROOT=${OLCF_XLC_ROOT}
export OLCF_XLMASS_ROOT=${OLCF_XLMASS_ROOT}
export OLCF_XLSMP_ROOT=${OLCF_XLSMP_ROOT}
__EOF__
  chmod +x summit-xlvars.sh
}

# PGI's OpenMP implementation relies on GCC's libatomic
# This function extracts the GCC compiler path and saves it to xlvars.sh
getgccvars() {
  module load gcc
  echo "export OLCF_GCC_ROOT=${OLCF_GCC_ROOT}" > summit-gccvars.sh
  chmod +x summit-gccvars.sh
}

# Load compiler module and TAU (if needed)
# Note that TAU 2.29 is botched (doesn't link in MPI libs properly)
# Use tau/2.28.1_patched (but this one doesn't support IBM XL)
# TODO: Resolve ticket #419691
case ${COMPILER} in
  ibm)
    module load xl
    ;;
  pgi)
    getxlvars
    getgccvars
    module load pgi
    source ./summit-xlvars.sh
    source ./summit-gccvars.sh
    ;;
  gcc)
    getxlvars
    module load gcc
    source ./summit-xlvars.sh
    ;;
  llvm)
    echo "Unsupported compiler (TODO: write make.inc.summit.llvm.cpu)"
    exit 1
    #getxlvars
    #module load llvm
    #source ./summit-xlvars.sh
    ;;
  tau-ibm)
    echo "Unsupported compiler (see OLCF ticket #419691)"
    exit 1
    #module load tau/2.29
    #module load xl
    #export TAU_MAKEFILE=${TAU_DIR}/lib/Makefile.tau-xl-papi-mpi-cupti-openmp
    #export USETAU=1
    ;;
  tau-pgi)
    echo "Unsupported compiler (see OLCF ticket #419691)"
    exit 1
    #getxlvars
    #getgccvars
    #module load tau/2.28.1_patched
    #module load pgi
    #source ./summit-xlvars.sh
    #source ./summit-gccvars.sh
    #export TAU_MAKEFILE=${TAU_DIR}/lib/Makefile.tau-pgi-papi-mpi-cupti-pdt-pgi
    #export USETAU=1
    ;;
  tau-gcc)
    echo "Unsupported compiler (TODO: write make.inc.summit.tau-gcc.cpu)"
    exit 1
    #getxlvars
    #module load tau/2.28.1_patched
    #module load gcc
    #source ./summit-xlvars.sh
    #export TAU_MAKEFILE=${TAU_DIR}/lib/Makefile.tau-gnu-papi-gnu-mpi-cupti-pdt-openmp
    #export USETAU=1
    ;;
  tau-llvm)
    echo "Unsupported compiler (TODO: write make.inc.summit.tau-llvm.cpu)"
    exit 1
    ;;
  *)
    echo "Unsupported compiler"
    exit 1
esac

# Copy make.inc
# TODO: Write the unavailable make.inc files
cp make.inc.summit.${COMPILER}.cpu make.inc

# For now, this is always false
# TODO: Resolve ticket #419691
if [ ${USETAU} ]; then
  source ./summit-xlvars.sh
  # Get vars from make.inc
  make lsvars
  source ./libs.sh
  # Apply options
  export TAU_OPTIONS="-optCompInst -optRevert -optTrackIO -optLinking=\"${LIBS}\""
fi

# Clean up
make clean
rm *.o *.mod
rm src/elk-cpu bin/elk-cpu
#rm src/elk-gpu bin/elk-gpu

# Make the binary
make all

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
[[ -d bin ]] || mkdir bin
[[ -e bin/elk ]] || rm bin/elk
cd bin
cp ../src/elk-cpu elk-cpu-${COMPILER}
cd ..
