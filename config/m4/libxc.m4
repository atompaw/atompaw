# -*- Autoconf -*-
#
# M4 macros for AtomPAW
#
# Copyright (C) 2010-2017 Marc Torrent
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
# Checks that the selected libXC library properly works.
#
AC_DEFUN([ATP_LIBXC_CHECK],[

  atp_libxc_ok="unknown"

  if test "${fc_has_iso_c_binding}" = "no"; then
    atp_libxc_ok="no"
    AC_MSG_ERROR([>>>> LibXC library needs a Fortran 2003 compiler supporting ISO_C_BINDINGS!])
  fi

  atp_saved_FCFLAGS="${FCFLAGS}"
  FCFLAGS="${atp_saved_FCFLAGS} ${CPPFLAGS}"

  AC_MSG_CHECKING([whether LibXC library works])

  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
     [#include "xc.h"
      #include "xc_funcs.h"],
     [[xc_func_type func;
       double exc[5];
       double rho[3] = {0.1, 0.2, 0.3};
       int err,func_id = 1;
       err=xc_func_init(&func, func_id, XC_UNPOLARIZED);
       xc_lda_exc(&func, 3, rho, exc);
     ]])],[atp_libxc_ok="yes"], [atp_libxc_ok="no"])
  AC_LANG_POP([C])

  AC_MSG_RESULT([${atp_libxc_ok}])

  FCFLAGS="${atp_saved_FCFLAGS}"

  if test "${atp_libxc_ok}" = "yes"; then
    AC_DEFINE([HAVE_LIBXC],1,[Define to 1 if you want to use LibXC.])
    AC_MSG_NOTICE([===== LibXC library OK for use])
    AC_MSG_NOTICE([])
  else
    AC_MSG_ERROR([No working LibXC library has been found!])
  fi

]) # ATP_LIBXC_CHECK

# -----------------------------------------------------------------------------

# ATP_LIBXC_SEARCH()
# ------------------
#
# Search were libXC is located.
#
AC_DEFUN([ATP_LIBXC_SEARCH],[

  atp_libxc_found="unknown"

  list_dirs="empty /opt/local /usr/local /opt/etsf /usr / /Library/libxc /usr/local/libxc ${LIBXC_ROOT} ${LIBXC_DIR}"
  list_headers="xc.h xc_funcs.h"
  list_funcs="xc_func_init"

  atp_saved_CPPFLAGS="${CPPFLAGS}"
  atp_saved_LDFLAGS="${LDFLAGS}"
  atp_saved_LIBS="${LIBS}"

  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([===== Looking for LibXC headers and libraries])

  AC_LANG_PUSH([C])

  atp_libxc_headers_ok="no"
  atp_libxc_libs_ok="no"

  if test "x${LIBXC_LIBDIR}" != "x" -a "x${with_libxc_prefix}" = "x"; then
    with_libxc_prefix="${LIBXC_LIBDIR}"
  fi

# First, check with user-defined macros
  if test "x${with_libxc_incs}" != "x" -o "x${with_libxc_libs}" != "x" -o "x${with_libxc_prefix}" != "x"; then
    AC_MSG_NOTICE([...LibXC: trying with command-line options])

    ATP_LIBXC_RESET_HEADERS_CACHE(${list_headers})
    ATP_LIBXC_RESET_HEADERS_CACHE(${list_funcs})

    if test "x${with_libxc_incs}" != "x" -o "x${with_libxc_prefix}" != "x"; then
      CPPFLAGS="";LDFLAGS="";LIBS=""
      if test "x${with_libxc_incs}" != "x"; then
        CPPFLAGS="${with_libxc_incs}"
      else
        CPPFLAGS="-I${with_libxc_prefix}/include"
      fi
      atp_libxc_headers_ok="yes"
      AC_CHECK_HEADERS(${list_headers},[],[atp_libxc_headers_ok="no"])
      if test "${atp_libxc_headers_ok}" = "yes"; then
        with_libxc_incs="${CPPFLAGS}"
      else
        AC_MSG_ERROR([LibXC header files (xc.h, ...) were not found with the specified --with_libxc_incs/--with_libxc_prefix])
      fi
    fi
    if test "x${with_libxc_libs}" != "x" -o "x${with_libxc_prefix}" != "x"; then
      CPPFLAGS="";LDFLAGS="";LIBS=""
      if test "x${with_libxc_libs}" != "x"; then
        LDFLAGS="${with_libxc_libs}"
        LIBS="${with_libxc_libs}"
      else
        LDFLAGS="-L${with_libxc_prefix}/lib -lxc"
        LIBS="-L${with_libxc_prefix}/lib -lxc"
      fi
      atp_libxc_libs_ok="yes"
      AC_SEARCH_LIBS(${list_funcs},[xc],[],[atp_libxc_libs_ok="no"],[-lm])
      if test "${atp_libxc_libs_ok}" = "yes"; then
        with_libxc_libs="${LDFLAGS} -lm"
      else
        AC_MSG_ERROR([LibXC library (libxc.*) was not found with the specified --with_libxc_libs/--with_libxc_prefix option])
      fi
    fi
  fi

# Then, check in a list of directories
  if test "x${with_libxc_incs}" = "x" -o "x${with_libxc_libs}" = "x"; then

    for with_libxc_prefix in ${list_dirs}; do

      ATP_LIBXC_RESET_HEADERS_CACHE(${list_headers})
      ATP_LIBXC_RESET_LIBS_CACHE(${list_funcs})

      CPPFLAGS="";LDFLAGS="";LIBS=""
      if test "x${with_libxc_prefix}" = "xempty"; then
        AC_MSG_NOTICE([...LibXC: trying with current environment])
        CPPFLAGS="";LDFLAGS="-lxc";LIBS="-lxc"
      else
        AC_MSG_NOTICE([...LibXC: trying in ${with_libxc_prefix}])
        if test "x${with_libxc_incs}" != "x"; then
          CPPFLAGS="${with_libxc_incs}"
        else
          CPPFLAGS="-I${with_libxc_prefix}/include"
        fi
        if test "x${with_libxc_libs}" != "x"; then
          LDFLAGS="${with_libxc_libs}"
          LIBS="${with_libxc_libs}"
        else
          LDFLAGS="-L${with_libxc_prefix}/lib -lxc"
          LIBS="-L${with_libxc_prefix}/lib -lxc"
        fi
      fi
      atp_libxc_headers_ok="yes"
      AC_CHECK_HEADERS(${list_headers},[],[atp_libxc_headers_ok="no"])
      if test "${atp_libxc_headers_ok}" = "yes"; then
        atp_libxc_libs_ok="yes"
        AC_SEARCH_LIBS(${list_funcs},[xc],[],[atp_libxc_libs_ok="no"],[-lm])
      fi
      if test "${atp_libxc_headers_ok}" = "yes" -a "${atp_libxc_libs_ok}" = "yes"; then
        with_libxc_incs="${CPPFLAGS}"
        with_libxc_libs="${LDFLAGS}"
        break
      fi

    done

  fi

  AC_LANG_POP([C])

  CPPFLAGS="${atp_saved_CPPFLAGS}"
  LDFLAGS="${atp_saved_LDFLAGS}"
  LIBS="${atp_saved_LIBS}"

  if test "${atp_libxc_headers_ok}" = "yes" -a "${atp_libxc_libs_ok}" = "yes"; then
    atp_libxc_found="yes"
    CPPFLAGS="${CPPFLAGS} ${with_libxc_incs}"
    LDFLAGS="${with_libxc_libs} ${LDFLAGS}"
    LIBS="${with_libxc_libs} ${LIBS}"
    AC_MSG_NOTICE([=> A LibXC library was found on the system])
  else
    atp_libxc_found="no"
    AC_MSG_ERROR([>>>> No LibXC library (libxc.*) was found on the system!])
  fi

]) # ATP_LIBXC_SEARCH

# -----------------------------------------------------------------------------

# ATP_RESET_LIBXC_[HEADERS/LIBS]_CACHE
# ------------------------------------
#
# Helper functions to reset headers/libs cache
#
AC_DEFUN([ATP_LIBXC_RESET_HEADERS_CACHE], [
    AS_FOR([AX_var], [ax_var], [$1], [
        AS_VAR_PUSHDEF([ax_Var], [ac_cv_header_${ax_var}])
        AS_UNSET([ax_Var])
        AS_VAR_POPDEF([ax_Var])
    ])
]) # ATP_LIBXC_RESET_HEADERS_CACHE
AC_DEFUN([ATP_LIBXC_RESET_LIBS_CACHE], [
    AS_FOR([AX_var], [ax_var], [$1], [
        AS_VAR_PUSHDEF([ax_Var], [ac_cv_lib_${ax_var}])
        AS_UNSET([ax_Var])
        AS_VAR_POPDEF([ax_Var])
    ])
]) # ATP_LIBXC_RESET_LIBS_CACHE
