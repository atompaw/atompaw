# -*- Autoconf -*-
#
# M4 macros for AtomPAW (imported from Abinit)
#
# Copyright (C) 2010 Yann Pouillon
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Linear algebra support
#



# ATP_LINALG_CHECK()
# ------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([ATP_LINALG_CHECK],[
  dnl Init
  atp_linalg_ok="unknown"

  dnl Check BLAS and LAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries work])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
      call zhpev
    ]])], [atp_linalg_ok="yes"], [atp_linalg_ok="no"])
  AC_MSG_RESULT([${atp_linalg_ok}])
  AC_LANG_POP([Fortran])
]) # ATP_LINALG_CHECK
