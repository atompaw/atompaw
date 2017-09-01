# -*- Autoconf -*-
#
# M4 macros for AtomPAW
#
# Copyright (C) 2017 Marc Torrent
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# FORTRAN 2002 ISO_C_BINDINGS support
#

# ATP_FC_ISO_C_BINDING_CHECK()
# ----------------------------
#
# Checks whether the Fortran compiler provides the intrinsic module ISO_C_BINDING.
#
AC_DEFUN([ATP_FC_ISO_C_BINDING_CHECK],[
  fc_has_iso_c_binding="no"

  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether the Fortran compiler provides the iso_c_binding module])

  dnl Try to compile a simple piece of code using iso_c_binding
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[use iso_c_binding
      implicit none
      integer(c_int) :: ii
      logical :: lbool
      type(c_ptr) :: ptr
      ptr = c_null_ptr
      lbool = c_associated(ptr)
    ]])], [fc_has_iso_c_binding="yes"])
  AC_LANG_POP([Fortran])

  if test "${fc_has_iso_c_binding}" = "yes"; then
    AC_DEFINE([HAVE_FC_ISO_C_BINDING],1,[Define to 1 if your Fortran compiler provides the iso_c_binding module.])
  else
    AC_MSG_ERROR([Fortran compiler does not provide iso_c_binding module. Use a more recent version or a different compiler.])
  fi

  AC_MSG_RESULT(${fc_has_iso_c_binding})
  AC_LANG_POP([Fortran])

]) # ATP_FC_ISO_C_BINDING_CHECK
