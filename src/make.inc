# lapack and blas libraries
LIBS=-lblas -llapack

# Compilers and flags
F90 = gfortran
F90FLAGS = -O0 -g
LDFLAGS = $(F90FLAGS)

