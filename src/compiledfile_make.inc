# lapack and blas libraries
LIBS=/home/natalie/EL4/publiccode/LAPACK/lapack_intel111.a \
	/home/natalie/EL4/publiccode/LAPACK/blas_intel111.a

# Compilers and flags
F90 = /opt/intel/compiler111/ifort/bin/intel64/ifort
F90FLAGS = -g -debug all -xHOST -fp-model precise -static-intel -check all
LDFLAGS = $(F90FLAGS)
