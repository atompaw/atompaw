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
  atp_saved_LDFLAGS="${LDFLAGS}"
  atp_saved_LIBS="${LIBS}"
  FCFLAGS="${atp_saved_FCFLAGS} ${CPPFLAGS}"
  LDFLAGS="${atp_saved_LDFLAGS} -lm"
  LIBS="${atp_saved_LIBS} -lm"

  AC_MSG_CHECKING([whether LibXC library works])

  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
     [#include "xc.h"
      #include "xc_funcs.h"],
     [[
       xc_func_type func;
       double exc[5];
       double rho[3] = {0.1, 0.2, 0.3};
       int err,func_id = 1;
       err=xc_func_init(&func, func_id, XC_UNPOLARIZED);
       xc_lda_exc(&func, 3, rho, exc);
     ]])],[atp_libxc_ok="yes"], [atp_libxc_ok="no"])
  AC_LANG_POP([C])

  if test "${atp_libxc_ok}" = "yes"; then
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
     [[
       use iso_c_binding
       integer(c_int) :: ok
       interface
         integer(c_int) function xc_func_init() bind(C)
         use iso_c_binding, only : c_int
         end function xc_func_init
       end interface
       ok=xc_func_init()
       ]])],[atp_libxc_ok="yes"], [atp_libxc_ok="no"])
    AC_LANG_POP([Fortran])
  fi

  AC_MSG_RESULT([${atp_libxc_ok}])

  FCFLAGS="${atp_saved_FCFLAGS}"
  LDFLAGS="${atp_saved_LDFLAGS}"
  LIBS="${atp_saved_LIBS}"

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

  list_dirs="empty /opt/local /usr/local /opt/etsf /usr / /Library/libxc /usr/local/libxc ${LIBXC_ROOT} ${LIBXC_DIR} ${LIBXC_PREFIX}"
  
  #We shouldnt use this (literals are encouraged)
  list_headers="xc.h xc_funcs.h"
  list_funcs="xc_func_init xc_func_end"

  atp_saved_CPPFLAGS="${CPPFLAGS}"
  atp_saved_LDFLAGS="${LDFLAGS}"
  atp_saved_LIBS="${LIBS}"

  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([===== Looking for LibXC headers and libraries])

  AC_LANG_PUSH([C])

  if test "x${LIBXC_PREFIX}" != "x" -a "x${with_libxc_prefix}" = "x"; then
    with_libxc_prefix="${LIBXC_PREFIX}"
  fi

  atp_libxc_headers_ok="no"
  atp_libxc_libs_ok="no"

# First, check with user-defined macros
  if test "x${with_libxc_incs}" != "x" -o "x${with_libxc_libs}" != "x" -o "x${with_libxc_prefix}" != "x"; then
    AC_MSG_NOTICE([...LibXC: trying with command-line options])

#   Headers
    ATP_LIBXC_RESET_HEADERS_CACHE([${list_headers}])
    ATP_LIBXC_RESET_HEADERS_CACHE([${list_funcs}])
    if test "x${with_libxc_incs}" != "x" -o "x${with_libxc_prefix}" != "x"; then
      CPPFLAGS="";LDFLAGS="";LIBS=""
      if test "x${with_libxc_incs}" != "x"; then
        CPPFLAGS="${with_libxc_incs}"
      else
        CPPFLAGS="-I${with_libxc_prefix}/include"
      fi
      AC_CHECK_HEADERS([xc.h xc_funcs.h],[atp_libxc_headers_ok="yes"],
                                         [atp_libxc_headers_ok="no"])
#     AC_CHECK_HEADERS([${list_headers}],[atp_libxc_headers_ok="yes"],
#                                        [atp_libxc_headers_ok="no"])
      if test "${atp_libxc_headers_ok}" = "yes"; then
        with_libxc_incs="${CPPFLAGS}"
      else
        AC_MSG_ERROR([LibXC header files were not found with the specified --with_libxc_incs/--with_libxc_prefix])
      fi
    fi

#   Libraries
    if test "x${with_libxc_libs}" != "x" -o "x${with_libxc_prefix}" != "x"; then
      CPPFLAGS="";LDFLAGS="";LIBS=""
      if test "x${with_libxc_libs}" != "x"; then
        LDFLAGS="${with_libxc_libs}"
        LIBS="${with_libxc_libs}"
      else
        LDFLAGS="-L${with_libxc_prefix}/lib -lxc"
        LIBS="-L${with_libxc_prefix}/lib -lxc"
      fi
      ATP_LIBXC_SEARCH_LIB_FORTRAN([${list_funcs}],[xc],[atp_libxc_libs_ok="yes"],
                                   [atp_libxc_libs_ok="no"],[-lm])
      if test "${atp_libxc_libs_ok}" = "yes"; then
        with_libxc_libs="${LDFLAGS} -lm"
      else
        if test "x${with_libxc_prefix}" != "x"; then
          LDFLAGS="-L${with_libxc_prefix}/lib64 -lxc"
          LIBS="-L${with_libxc_prefix}/lib64 -lxc"
          ATP_LIBXC_SEARCH_LIB_FORTRAN([${list_funcs}],[xc],[atp_libxc_libs_ok="yes"],
                                       [atp_libxc_libs_ok="no"],[-lm])
          if test "${atp_libxc_libs_ok}" = "yes"; then
            with_libxc_libs="${LDFLAGS} -lm"
          fi
        fi
      fi
      if test "${atp_libxc_libs_ok}" = "no"; then
        AC_MSG_ERROR([LibXC library was not found with the specified --with_libxc_libs/--with_libxc_prefix option])
      fi
    fi
  fi

# Then, check in a list of directories
  if test "${atp_libxc_incs_ok}" = "no" -o "${atp_libxc_libs_ok}" = "no"; then
    if test "x${with_libxc_incs}" = "x" -o "x${with_libxc_libs}" = "x"; then

#     Loop over dirs
	  for with_libxc_prefix in ${list_dirs}; do
		CPPFLAGS="";LDFLAGS="";LIBS=""
		ATP_LIBXC_RESET_HEADERS_CACHE([${list_headers}])
		ATP_LIBXC_RESET_LIBS_CACHE([${list_funcs}])

#       Current environment
		if test "x${with_libxc_prefix}" = "xempty"; then
		  AC_MSG_NOTICE([...LibXC: trying with current environment])

#         Headers
		  CPPFLAGS="";LDFLAGS="-lxc";LIBS="-lxc"
		  AC_CHECK_HEADERS([xc.h xc_funcs.h],[atp_libxc_headers_ok="yes"],
                                             [atp_libxc_headers_ok="no"])
#	      AC_CHECK_HEADERS([${list_headers}],[atp_libxc_headers_ok="yes"],
#                                            [atp_libxc_headers_ok="no"])

#         Libraries
		  if test "${atp_libxc_headers_ok}" = "yes"; then
            ATP_LIBXC_SEARCH_LIB_FORTRAN([${list_funcs}],[xc],[atp_libxc_libs_ok="yes"],
                                         [atp_libxc_libs_ok="no"],[-lm])
		  fi

		else

#         Other directories
		  AC_MSG_NOTICE([...LibXC: trying in ${with_libxc_prefix}])

#         Headers
		  if test "x${with_libxc_incs}" != "x"; then
			CPPFLAGS="${with_libxc_incs}"
		  else
			CPPFLAGS="-I${with_libxc_prefix}/include"
		  fi
		  atp_libxc_headers_ok="yes"
#		  AC_CHECK_HEADERS([${list_headers}],[],[atp_libxc_headers_ok="no"])
		  AC_CHECK_HEADERS([xc.h xc_funcs.h],[],[atp_libxc_headers_ok="no"])

#         Libraries
		  if test "${atp_libxc_headers_ok}" = "yes"; then
			if test "x${with_libxc_libs}" != "x"; then
			  LDFLAGS="${with_libxc_libs}"
			  LIBS="${with_libxc_libs}"
			else
			  LDFLAGS="-L${with_libxc_prefix}/lib -lxc"
			  LIBS="-L${with_libxc_prefix}/lib -lxc"
			fi
            ATP_LIBXC_SEARCH_LIB_FORTRAN([${list_funcs}],[xc],[atp_libxc_libs_ok="yes"],
                                         [atp_libxc_libs_ok="no"],[-lm])
			if test "${atp_libxc_libs_ok}" != "yes"; then
			  if test "x${with_libxc_prefix}" != "x" ; then
			    LDFLAGS="-L${with_libxc_prefix}/lib64 -lxc"
			    LIBS="-L${with_libxc_prefix}/lib64 -lxc"
			    atp_libxc_libs_ok="yes"
                ATP_LIBXC_SEARCH_LIB_FORTRAN([${list_funcs}],[xc],[atp_libxc_libs_ok="yes"],
                                         [atp_libxc_libs_ok="no"],[-lm])
              fi
			fi
		  fi
		fi

#       If found, tore final values of flags
        if test "${atp_libxc_headers_ok}" = "yes" -a "${atp_libxc_libs_ok}" = "yes"; then
	      AC_MSG_NOTICE([=> libXC library was found on the system])
		  with_libxc_incs="${CPPFLAGS}"
		  with_libxc_libs="${LDFLAGS} -lm"
		  break
		fi

	  done

    fi
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

# ATP_LIBXC_SEARCH_LIB_FORTRAN
# ----------------------------
#  (list-funcs, search-lib , [action-if-found],
#   [action-if-not-found], [other-libraries])
#
# Helper function to check a LIB using Fortran calling C
#
AC_DEFUN([ATP_LIBXC_SEARCH_LIB_FORTRAN], [

  AC_MSG_CHECKING([for library containing $1])

  ac_pgm1="" ; ac_pgm2=""
  for fnc in $1; do
    ac_pgm1="${ac_pgm1}
       interface
         integer(c_int) function ${fnc}() bind(C)
         use iso_c_binding, only : c_int
         end function ${fnc}
       end interface";
    ac_pgm2="${ac_pgm2}
       ok=${fnc}()";
  done

  LDFLAGS_SAV="${LDFLAGS}"
  LIBS_SAV="${LIBS}"
  LDFLAGS="${LDFLAGS} $5"
  LIBS="${LIBS} $5"

  AC_LANG_PUSH([Fortran])

  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
     [[
       use iso_c_binding
       integer(c_int) :: ok
       ${ac_pgm1}
       ${ac_pgm2}
       ]])],[lib_fc_ok="yes";$3], [lib_fc_ok="no";$4])

  AC_LANG_POP([Fortran])

  LDFLAGS="${LDFLAGS_SAV}"
  LIBS="${LIBS_SAV}"

  AC_MSG_RESULT([${lib_fc_ok}])

]) # ATP_LIBXC_SEARCH_LIB_FORTRAN


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
