# -*- Autoconf -*-
#
# Configure script for the AtomPAW package
#
# Copyright (C) 2010-2017 Yann Pouillon, Marc Torrent
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
AC_INIT([AtomPAW],[4.2.0.0],[natalie@wfu.edu, marc.torrent@cea.fr],[atompaw])
AC_REVISION([Autotools support for AtomPAW])
AC_CONFIG_AUX_DIR(config/gnu)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_SRCDIR(src/atompaw_prog.F90)

# Regenerate the VERSION file
# Note: remember that some editors break long lines
AC_MSG_NOTICE([regenerating the VERSION file])
_AC_SRCDIRS(["."])
cat >"${ac_abs_top_srcdir}/VERSION" <<EOF
!- Version ${PACKAGE_VERSION} - 4/2021 - NH,MT added self-consistent MGGA (beta)
#!- Version ${PACKAGE_VERSION} - 11/2020 - MT,FJ added LDA-1/2 potential
EOF
#!--- Version ${PACKAGE_VERSION} - 11/2020 - MT introduce std I/O variables
#!--- Version ${PACKAGE_VERSION} - 11/2020 - MT change input reading design
#!--- Version ${PACKAGE_VERSION} - 06/2019 - MT,FJ,NH overlap>0, core kinetic
#!--- Version ${PACKAGE_VERSION} - 12/2018 - MT modified abinitinterface pbesol
#!--- Version ${PACKAGE_VERSION} - 09/2018 - NH check for positive xc argument
#!--- Version ${PACKAGE_VERSION} - 07/2018 - MT fixed some libxc issues 
#!--- Version ${PACKAGE_VERSION} - 05/2018 - mod. abinitinterface, other updates
#!--- Version ${PACKAGE_VERSION} - 03/2018 - MT enabled libxc4
#!--- Version ${PACKAGE_VERSION} - 01/2018 - Dirac solver available for graphatom
#!--- Version ${PACKAGE_VERSION} - 09/2017 - LibXC C API now used + libXC detection improved
#!--- Version ${PACKAGE_VERSION} - 12/2016 - Finite-Nucleus fixed for DFT
#!--- Version ${PACKAGE_VERSION} - 02/2016 - XML output improved; --version command line option
#!--- Version ${PACKAGE_VERSION} - 03/2015 - PBESOL ; SOCCORO output
#!--- Version ${PACKAGE_VERSION} - 08/2014 - portability corrections; rename f90 in F90
#!--- Version ${PACKAGE_VERSION} - 08/2014 - added exchange core-valence terms
#!--- Version ${PACKAGE_VERSION} - 01/2014 - solver, coretail, xml corr.
#!--- Version ${PACKAGE_VERSION} - 11/2013 - solver fixed; 1/2014 coretail

# Initialize Automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.10)
AC_CONFIG_HEADERS([config.h])
#AM_CONFIG_HEADER([config.h])

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
    [LibXC includes flags for headers (e.g. -I...)])])

AC_ARG_WITH([libxc-libs],[
  AS_HELP_STRING([--with-libxc-libs],
    [LibXC library flags to link with (e.g. -L... -lxc)])])

AC_ARG_WITH([libxc-prefix],[
  AS_HELP_STRING([--with-libxc-prefix],
    [LibXC installation root directory (e.g. /usr/local, ...)])])

AC_ARG_WITH([linalg-libs],[
  AS_HELP_STRING([--with-linalg-libs],
    [Linear algebra libraries flags to link with (e.g. -L... -lblas -llapack)])])

AC_ARG_WITH([linalg-prefix],[
  AS_HELP_STRING([--with-linalg-prefix],
    [Linear Algebra libraries installation root directory (e.g. /usr/local, ...)])])

# Declare environment variables
AC_ARG_VAR([LINALG_PREFIX],
  [Directory where to search for linear algebra libraries])

# Declare environment variables
AC_ARG_VAR([LIBXC_PREFIX],
  [Directory where to search for libXC library])

                    # ------------------------------------ #

#
# Basic checking 
#

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
ATP_FC_ISO_C_BINDING_CHECK
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
# Library search
#

# NOTE: least dependent libraries first

# Linear Algebra: search for the libraries and check that they properly work
ATP_LINALG_SEARCH
ATP_LINALG_CHECK

# LIBXC: search for the library and check that it properly works
if test "${enable_libxc}" = "yes"; then
  ATP_LIBXC_SEARCH
  ATP_LIBXC_CHECK
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

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/pkginfo.f90
])
AC_OUTPUT
