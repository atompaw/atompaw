# -*- Autoconf -*-
#
# M4 macros for AtomPAW
#
# Copyright (C) 2010 Yann Pouillon
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Optimizations
#



# ATP_FC_OPTFLAGS()
# -----------------
#
# Sets optimization flags for the Fortran compiler.
#
AC_DEFUN([ATP_FC_OPTFLAGS],[
  ATP_FCOPTS=""
  ATP_LDOPTS=""
  ATP_LIBS=""

  case "${atp_fc_vendor}" in

    compaq)
      ATP_FCOPTS="-fast"
      ATP_LDOPTS="-fast"
      ATP_LIBS="-lcxml"
      ;;

    g95)
      ATP_FCOPTS="-O2"
      ATP_LDOPTS="-O2"
      case "${with_linalg_flavor}" in
        atlas)
          ATP_LIBS="-llapack -lf77blas -lcblas -latlas"
          ;;
        netlib)
          ATP_LIBS="-llapack -lblas"
          ;;
      esac
      ;;

    gnu)
      ATP_FCOPTS="-O2"
      ATP_LDOPTS="-O2"
      case "${with_linalg_flavor}" in
        atlas)
          ATP_LIBS="-llapack -lf77blas -lcblas -latlas"
          ;;
        netlib)
          ATP_LIBS="-llapack -lblas"
          ;;
      esac
      ;;

    intel)
      ATP_FCOPTS="-O2"
      ATP_LDOPTS="-O2 -Vaxlib"
      case "${with_linalg_flavor}" in
        atlas)
          ATP_LIBS="-llapack -lf77blas -lcblas -latlas"
          ;;
        mkl)
          ATP_LIBS="-lmkl_lapack -lmkl_ia32 -lguide -lpthread"
          ;;
      esac
      ;;

    nag)
      ATP_FCOPTS="-DNAG"
      ATP_LDOPTS=""
      case "${with_linalg_flavor}" in
        atlas)
          ATP_LIBS="-llapack -lf77blas -lcblas -latlas"
          ;;
        netlib)
          ATP_LIBS="-llapack -lblas"
          ;;
      esac
      ;;

    pathscale)
      ATP_FCOPTS="-O3"
      ATP_LDOPTS=""
      ATP_LIBS="-lacml"
      ;;

    pgi)
      ATP_FCOPTS="-fast"
      ATP_LDOPTS=""
      ATP_LIBS="-llapack -lblas"
      ;;

    sun)
      ATP_FCOPTS="-fast -xO3 -xlic_lib=sunperf"
      ATP_LDOPTS="-fast -xO3 -xlic_lib=sunperf"
      ATP_LIBS=""
      ;;

    ibm)
      ATP_FCOPTS="-qfree=f90 -qsuffix=f=f90:cpp=F90 -qarch=auto -O3"
      ATP_LDOPTS="-qfree=f90 -qsuffix=f=f90:cpp=F90 -qarch=auto"
      ATP_LIBS="-llapack -lessl"
      ;;

    *)
      ATP_FCOPTS="-O2"
      ATP_LDOPTS=""
      ATP_LIBS=""
      ;;

  esac

  AC_SUBST(ATP_FCOPTS)
  AC_SUBST(ATP_LDOPTS)
  AC_SUBST(ATP_LIBS)
]) # ATP_FC_OPTFLAGS
