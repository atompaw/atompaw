# -*- Autoconf -*-
#
# M4 macros for AtomPAW (imported from Abinit)
#
# Copyright (C) 2017 Yann Pouillon, Marc Torrent
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

  atp_linalg_ok="unknown"

  AC_MSG_CHECKING([whether Linear Algebra libraries work])

# Check BLAS and LAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
      call zhpev
    ]])], [atp_linalg_ok="yes"], [atp_linalg_ok="no"])
  AC_LANG_POP([Fortran])

  AC_MSG_RESULT([${atp_linalg_ok}])

  if test "${atp_linalg_ok}" = "yes"; then
    AC_MSG_NOTICE([===== Linear Algebra libraries OK for use])
    AC_MSG_NOTICE([])
  else
    AC_MSG_ERROR([No working Linear Algebra libraries have been found!])
  fi

]) # ATP_LINALG_CHECK

# -----------------------------------------------------------------------------

# ATP_LINALG_SEARCH()
# ------------------
#
# Search were Linear Algebra is located.
#
AC_DEFUN([ATP_LINALG_SEARCH],[

  atp_linalg_found="unknown"

  list_dirs="empty /opt/local /usr/local /usr / /Library/lapack /usr/local/lapack ${LINALG_ROOT} ${LINALG_DIR}"
  list_funcs="zhpev zgemm"

  atp_saved_LDFLAGS="${LDFLAGS}"
  atp_saved_LIBS="${LIBS}"

  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([===== Looking for Linear Algebra libraries])

  AC_LANG_PUSH([Fortran])

  atp_linalg_libs_ok="no"

  if test "x${LINALG_LIBDIR}" != "x" -a "x${with_linalg_prefix}" = "x"; then
    with_linalg_prefix="${LINALG_LIBDIR}"
  fi

# First, check with user-defined macros
  if test "x${with_linalg_libs}" != "x" -o "x${with_linalg_prefix}" != "x"; then
    AC_MSG_NOTICE([...Linear Algebra: trying with command-line options])

    ATP_LINALG_RESET_LIBS_CACHE(${list_funcs})

    LDFLAGS="";LIBS=""
    if test "x${with_linalg_libs}" != "x"; then
      LDFLAGS="${with_linalg_libs}"
      LIBS="${with_linalg_libs}"
    else
      LDFLAGS="-L${with_linalg_prefix}/lib -lblas -llapack -lblas"
      LIBS="-L${with_linalg_prefix}/lib -lblas -llapack -lblas"
    fi
    atp_linalg_libs_ok="yes"
    for func in ${list_funcs}; do
      AC_SEARCH_LIBS(${func},[lapack],[],[atp_linalg_libs_ok="no"],[-lblas])
    done
    if test "${atp_linalg_libs_ok}" = "yes"; then
      with_linalg_libs="${LDFLAGS}"
    else
      AC_MSG_ERROR([Linear Algebra libraries (blas, lapack) were not found with the specified --with_linalg_libs/--with_linalg_prefix])
    fi
  fi

# Then, check in a list of directories
  if test "x${with_linalg_libs}" = "x"; then

    for with_linalg_prefix in ${list_dirs}; do

      ATP_LINALG_RESET_LIBS_CACHE(${list_funcs})

      LDFLAGS="";LIBS=""
      if test "x${with_linalg_prefix}" = "xempty"; then
        AC_MSG_NOTICE([...Linear Algebra: trying with current environment])
        LDFLAGS="-lblas -llapack -lblas";LIBS="-lblas -llapack -lblas"
      else
        AC_MSG_NOTICE([...Linear Algebra: trying in ${with_linalg_prefix}])
        if test "x${with_linalg_libs}" != "x"; then
          LDFLAGS="${with_linalg_libs}"
          LIBS="${with_linalg_libs}"
        else
          LDFLAGS="-L${with_linalg_prefix}/lib -lxc"
          LIBS="-L${with_linalg_prefix}/lib -lxc"
        fi
      fi
      atp_linalg_libs_ok="yes"
      for func in ${list_funcs}; do
        AC_SEARCH_LIBS(${func},[lapack],[],[atp_linalg_libs_ok="no"],[-lblas])
      done
      if test "${atp_linalg_libs_ok}" != "yes";then
        if test "x${with_linalg_prefix}" != "xempty"; then
          if test "x${with_linalg_prefix}" != "x"; then
            LDFLAGS="-L${with_linalg_prefix}/lib64 -lxc"
            LIBS="-L${with_linalg_prefix}/lib64 -lxc"
            atp_linalg_libs_ok="yes"
            for func in ${list_funcs}; do
              AC_SEARCH_LIBS(${func},[lapack],[],[atp_linalg_libs_ok="no"],[-lblas])
            done
          fi
        fi
      fi
      if test "${atp_linalg_libs_ok}" = "yes"; then
        with_linalg_libs="${LDFLAGS}"
        break
      fi

    done

  fi

  AC_LANG_POP([Fortran])

  LDFLAGS="${atp_saved_LDFLAGS}"
  LIBS="${atp_saved_LIBS}"

  if test "${atp_linalg_libs_ok}" = "yes"; then
    atp_linalg_found="yes"
    LDFLAGS="${with_linalg_libs} ${LDFLAGS}"
    LIBS="${with_linalg_libs} ${LIBS}"
    AC_MSG_NOTICE([=> Linear Algebra libraries were found on the system])
  else
    atp_linalg_found="no"
    AC_MSG_ERROR([>>>> No Linear Algebra libraries (libblas, liblapack) were found on the system!])
  fi

]) # ATP_LINALG_SEARCH

# -----------------------------------------------------------------------------

# ATP_RESET_LINALG_LIBS_CACHE
# ------- --------------------
#
# Helper function to reset libs cache
#
AC_DEFUN([ATP_LINALG_RESET_LIBS_CACHE], [
    AS_FOR([AX_var], [ax_var], [$1], [
        AS_VAR_PUSHDEF([ax_Var], [ac_cv_lib_${ax_var}])
        AS_UNSET([ax_Var])
        AS_VAR_POPDEF([ax_Var])
    ])
]) # ATP_LINALG_RESET_LIBS_CACHE
