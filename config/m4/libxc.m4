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
# LibXC support
#

# ATP_LIBXC_CHECK()
# ------------------
#
# Checks that the selected libXC libraries properly work.
#
AC_DEFUN([ATP_LIBXC_CHECK],[
  atp_libxc_ok="unknown"

  dnl Check LibXC routine

  AC_LANG_PUSH([Fortran])
  atp_saved_FCFLAGS="${FCFLAGS}"
  FCFLAGS="${atp_saved_FCFLAGS} ${CPPFLAGS}"

# First test with user-defined link flags

  AC_MSG_CHECKING([whether LibXC library works])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use xc_f90_types_m
      use xc_f90_lib_m
      implicit none
      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
      integer :: func_id = 1
      call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
      call xc_f90_func_end(xc_func)
    ]])], [atp_libxc_ok="yes"], [atp_libxc_ok="no"])
  AC_MSG_RESULT([${atp_libxc_ok}])

# Then test adding "-lxcf90" flag (needed from LibXC v2.2)

  if test "${atp_libxc_ok}" = "no"; then
    AC_MSG_CHECKING([whether LibXC library works with -lxcf90 link flag])
    atp_saved_LIBS="${LIBS}"
    LIBS=`echo "${atp_saved_LIBS}" | sed -e 's/-lxc/-lxcf90 -lxc/'`
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use xc_f90_types_m
        use xc_f90_lib_m
        implicit none
        TYPE(xc_f90_pointer_t) :: xc_func
        TYPE(xc_f90_pointer_t) :: xc_info
        integer :: func_id = 1
        call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
        call xc_f90_func_end(xc_func)
      ]])], [atp_libxc_ok="yes"], [atp_libxc_ok="no"])
    AC_MSG_RESULT([${atp_libxc_ok}])
    if test "${atp_libxc_ok}" = "no"; then
      LIBS="${atp_saved_LIBS}"
    fi
  fi

  FCFLAGS="${atp_saved_FCFLAGS}"
  AC_LANG_POP([Fortran])
]) # ATP_LIBXC_CHECK
