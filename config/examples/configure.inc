#!/bin/sh
#
# Configuration script for inc
#

# Allow for VPATH builds
if test -x "./configure"; then
  CFG_SCRIPT="./configure"
elif test -x "../configure"; then
  CFG_SCRIPT="../configure"
else
  echo "Error: could not find configure script"
  exit 1
fi

# Run configure
${CFG_SCRIPT} \
  FC="/opt/intel/compiler91/ifort/bin/ifort" \
  FCFLAGS="-O2 -tpp7 -xW -mp -static" \
  LDFLAGS="${FCFLAGS}" \
  LIBS="/home/natalie/EL4/publiccode/LAPACK/lapack_intel91.a \
    /home/natalie/EL4/publiccode/LAPACK/blas_intel91.a"
