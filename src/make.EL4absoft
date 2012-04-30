# lapack and blas libraries
ABSOFT=/opt/absoft-10.1
FC=$(ABSOFT)/bin/f90
LD_LIBRARY_PATH=$ABSOFT/lib64

LIBS =  $(ABSOFT)/lib64/libaU77.a -L$(ABSOFT)/lib64 \
        /opt/absoft-10.1/extras/lapack/precompiled/64bit/liblapack.a \
        /opt/absoft-10.1/extras/lapack/precompiled/64bit/libblas.a
                            
# Compilers and flags
F90 = $(FC)
F90FLAGS = -O2 -march=em64t -mtune=em64t -mcpu=em64t -YEXT_NAMES=LCS -YEXT_SFX="_"
LDFLAGS = $(F90FLAGS)

