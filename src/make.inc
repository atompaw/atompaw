# lapack and blas libraries
LIBS=-L/cea/home/abinit/abinit/APPLIS/g95/linalg -llapack -lblas

# Compilers and flags
F90 = g95
F90FLAGS = -O0 -g
LDFLAGS = $(F90FLAGS)

