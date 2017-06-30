# -*- Autoconf -*-
#
# Configure script for the AtomPAW package
#
# Copyright (C) 2010-2011 Yann Pouillon, MT
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Autotools startup
#

# Initialize Autoconf
AC_PREREQ(2.62)
AC_INIT([AtomPAW],[4.0.0.14],[natalie@wfu.edu, marc.torrent@cea.fr],[atompaw])
AC_REVISION([Autotools support for AtomPAW])
AC_CONFIG_AUX_DIR(config/gnu)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_SRCDIR(src/atompaw_prog.F90)

# Regenerate the VERSION file
# Note: remember that some editors break long lines
AC_MSG_NOTICE([regenerating the VERSION file])
_AC_SRCDIRS(["."])
cat >"${ac_abs_top_srcdir}/VERSION" <<EOF
!--- Version ${PACKAGE_VERSION} - 12/2016 -- Finite-Nucleus fixed for DFT
EOF
#!--- Version ${PACKAGE_VERSION} - 02/2016 -- XML output improved; --version command line option
#!--- Version ${PACKAGE_VERSION} - 03/2015 -- PBESOL ; SOCCORO output
#!--- Version ${PACKAGE_VERSION} - 08/2014 -- portability corrections; rename f90 in F90
#!--- Version ${PACKAGE_VERSION} - 08/2014 -- added exchange core-valence terms
#!--- Version ${PACKAGE_VERSION} - 06/2013; 1/2014 solver, coretail, xml corr.
#!--- Version ${PACKAGE_VERSION} - 06/2013; 11/2013 solver fixed; 1/2014 coretail

# Initialize Automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.10)
AM_CONFIG_HEADER([config.h])

                    # ------------------------------------ #

#
# Command-line arguments
#

# Declare arguments
AC_ARG_ENABLE([libxc],[
  AS_HELP_STRING([--enable-libxc],
    [Enable the use of LibXC for exchange-correlation])])

AC_ARG_WITH([fc-vendor],[
  AS_HELP_STRING([--with-fc-vendor],
    [Select a particular Fortran compiler vendor])])

AC_ARG_WITH([fc-version],[
  AS_HELP_STRING([--with-fc-version],
    [Select a particular Fortran compiler version])])

AC_ARG_WITH([libxc-incs],[
  AS_HELP_STRING([--with-libxc-incs],
    [LibXC includes for headers and modules])])

AC_ARG_WITH([libxc-libs],[
  AS_HELP_STRING([--with-libxc-libs],
    [LibXC libraries to link with])])

AC_ARG_WITH([linalg-flavor],[
  AS_HELP_STRING([--with-linalg-flavor],
    [Select a particular linear algebra library flavor (atlas, mkl, netlib)])])

AC_ARG_WITH([linalg-libs],[
  AS_HELP_STRING([--with-linalg-libs],
    [Linear algebra libraries to link with])])

# Declare environment variables
AC_ARG_VAR([ATOMPAW_LIBDIR],
  [Directory where to search for linear algebra libraries])

# Set default values of arguments
if test "${with_linalg_flavor}" = ""; then
  with_linalg_flavor="netlib"
  AC_MSG_NOTICE([setting linear algebra flavor to 'netlib'])
fi

                    # ------------------------------------ #

# Check for common programs
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_SED
AC_PROG_AWK
AC_PROG_GREP

# Workaround for the "grep -e" issue on Solaris systems
AC_PROG_EGREP

# Workaround for the wrong path to install-sh on Mac systems
ATP_PROG_MKDIR_P

# Check for other useful programs
AC_CHECK_PROGS(DVIPDF,[dvipdf])
AC_CHECK_PROGS(LATEX,[latex])
AC_CHECK_PROGS(PDFLATEX,[pdflatex])
AC_CHECK_PROGS(PERL,[perl])

                    # ------------------------------------ #

#
# Fortran support
#

# Look for the Fortran compiler
if test "${FC}" != "" -a ! -x "${FC}"; then
  atp_fc_probe=`echo "${FC}" | sed -e 's/ .*//'`
  if test ! -x "${atp_fc_probe}"; then
    AC_PATH_PROG([atp_fc_path],[${atp_fc_probe}])
    if test "${atp_fc_path}" = ""; then
      AC_MSG_ERROR([could not run Fortran compiler "${FC}"])
    fi
  fi
fi
AC_PROG_FC

# Fail if no Fortran compiler is available
if test "${FC}" = ""; then
  AC_MSG_ERROR([no Fortran compiler available])
fi

# Set AtomPAW Fortran parameters
atp_fc_vendor="${with_fc_vendor}"
atp_fc_version="${with_fc_version}"
ATP_PROG_FC

# Set default file extensions
ATP_FC_EXTENSIONS

# Get module file case
ATP_FC_MOD_CASE

# Set default optimizations
ATP_FC_OPTFLAGS

# Check for specific Fortran features
ATP_FC_FEATURES

                    # ------------------------------------ #

#
# Libtool configuration
#

# Initialize Libtool
#   Disable this shared-library feature because it is not convenient
#   to have a script and not an executable for atompaw
LT_INIT([disable-shared])
LT_PREREQ([2.2.2])
LTOBJEXT="lo"
AC_SUBST(LTOBJEXT)

                    # ------------------------------------ #

#
# Library configuration
#

# NOTE: least dependent libraries first

# Linear algebra
if test "${with_linalg_libs}" != ""; then
  LIBS="${with_linalg_libs} ${LIBS}"
else
  if test "${ATP_LIBS}" != ""; then
    LIBS="${ATP_LIBS} ${LIBS}"
  fi
  if test "${ATP_LIBDIR}" != ""; then
    LIBS="-L${ATP_LIBDIR} ${LIBS}"
  fi
fi

# LibXC
LIBXC_DEF_DIR="/opt/etsf"
if test "${enable_libxc}" = "yes"; then
  if test "${with_libxc_incs}" != ""; then
    CPPFLAGS="${CPPFLAGS} ${with_libxc_incs}"
  else
    CPPFLAGS="${CPPFLAGS} -I${LIBXC_DEF_DIR}/include"
  fi
  if test "${with_libxc_libs}" != ""; then
    LIBS="${with_libxc_libs} ${LIBS}"
  else
    LIBS="-L${LIBXC_DEF_DIR}/lib -lxcf90 -lxc ${LIBS}"
  fi
fi

                    # ------------------------------------ #

#
# Final configuration
#

# Note: in order to let the user override the defaults through the
#       command-line, the build flags (FCFLAGS, LDFLAGS) are set only if
#       they are found empty.

if test "${FCFLAGS}" = ""; then
  FCFLAGS="${ATP_FCOPTS}"
  if test "${ac_cv_prog_fc_g}" = "yes"; then
    FCFLAGS="-g ${FCFLAGS}"
  fi
fi
if test "${LDFLAGS}" = ""; then
  LDFLAGS="${ATP_LDOPTS}"
fi
if test "${LDFLAGS}" = ""; then
  LDFLAGS="${ATP_LDOPTS}"
fi

                    # ------------------------------------ #

#
# Final checks
#

# Check that linear algebra libraries properly work
ATP_LINALG_CHECK
if test "${atp_linalg_ok}" != "yes"; then
  AC_MSG_ERROR([no working linear algebra libraries have been found])
fi

# Check that LibXC properly works
if test "${enable_libxc}" = "yes"; then
  ATP_LIBXC_CHECK
  if test "${atp_libxc_ok}" = "yes"; then
    AC_DEFINE([HAVE_LIBXC],1,[Define to 1 if you want to use LibXC.])
    AC_MSG_NOTICE([=== libXC library OK for use])
    AC_MSG_NOTICE([=== generating the src/libxc_names.in file])
    ATP_MK_LIBXC_NAMES
  else
    AC_MSG_ERROR([no working LibXC library has been found])
  fi
fi

                    # ------------------------------------ #

#
# Output configuration
#

AC_MSG_NOTICE([FC       = ${FC}])
AC_MSG_NOTICE([CPPFLAGS = ${CPPFLAGS}])
AC_MSG_NOTICE([FCFLAGS  = ${FCFLAGS}])
AC_MSG_NOTICE([LDFLAGS  = ${LDFLAGS}])
AC_MSG_NOTICE([LIBS     = ${LIBS}])

AC_OUTPUT([
  Makefile
  src/Makefile
  src/pkginfo.f90
])
