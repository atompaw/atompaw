# -*- Autoconf -*-
#
# M4 macros for AtomPAW (imported from Abinit)
#
# Copyright (C) 2010 Yann Pouillon
#
# This file is part of the ATOMPAW software package. For license information,
# please see the COPYING file in the top-level directory of the ATOMPAW source
# distribution.
#

#
# Fortran compilers support
#



# _ATP_CHECK_FC_ABSOFT(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the ABSoft Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_ABSOFT],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the ABSoft Fortran compiler])
  fc_info_string=`$1 -V 2> /dev/null`
  atp_result=`echo "${fc_info_string}" | grep '^Pro Fortran'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_ABSOFT],1,[Define to 1 if you are using the ABSOFT Fortran compiler.])
    atp_fc_vendor="absoft"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/Pro Fortran //'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_ABSOFT



# _ATP_CHECK_FC_COMPAQ(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the COMPAQ Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_COMPAQ],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Compaq Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^Compaq Fortran Compiler'`
  atp_result="${fc_info_string}"
  if test "${atp_result}" = ""; then
    fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^HP Fortran Compiler'`
    atp_result="${fc_info_string}"
  fi
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_COMPAQ],1,[Define to 1 if you are using the COMPAQ Fortran compiler.])
    atp_fc_vendor="compaq"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.* V//;s/-.*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_COMPAQ



# _ATP_CHECK_FC_FUJITSU(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Fujitsu Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_FUJITSU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Fujitsu Fortran compiler])
  fc_info_string=`$1 -V 2> /dev/null`
  atp_result=`echo "${fc_info_string}" | grep '^Fujitsu Fortran'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_FUJITSU],1,[Define to 1 if you are using the Fujitsu Fortran compiler.])
    atp_fc_vendor="fujitsu"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.*Driver //;s/ .*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_FUJITSU



# _ATP_CHECK_FC_G95(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the G95 Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_G95],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the G95 Fortran compiler])
  fc_info_string=`$1 --version 2>&1`
  atp_result=`echo "${fc_info_string}" | grep '^G95'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_G95],1,[Define to 1 if you are using the G95 Fortran compiler.])
    atp_fc_vendor="g95"
    atp_fc_version=`echo ${atp_result} | sed -e 's/.*GCC //; s/ .*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_G95



# _ATP_CHECK_FC_GNU(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the GNU Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_GNU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the GNU Fortran compiler])
  fc_info_string=`$1 --version 2>&1`
  atp_result=`echo "${fc_info_string}" | grep '^GNU Fortran'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_GNU],1,[Define to 1 if you are using the GNU Fortran compiler.])
    AC_DEFINE([HAVE_FORTRAN2003],1,[Define to 1 if your Fortran compiler supports Fortran 2003.])
    atp_fc_vendor="gnu"
    atp_fc_version=`echo ${atp_result} | sed -e 's/^[[^(]]*([[^)]]*) //; s/ .*//'`
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_GNU



# _ATP_CHECK_FC_HITACHI(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Hitachi Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_HITACHI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Hitachi Fortran compiler])
  fc_info_string=`$1 -V 2> /dev/null`
  atp_result=`echo "${fc_info_string}" | grep '^Hitachi Fortran'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_HITACHI],1,[Define to 1 if you are using the Hitachi Fortran compiler.])
    atp_fc_vendor="hitachi"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.*Driver //;s/ .*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_HITACHI



# _ATP_CHECK_FC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the IBM XL Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_IBM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the IBM XL Fortran compiler])
  fc_info_string=`$1 -qversion 2>&1 | head -n 1`
  fc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
  atp_result=`echo "${fc_info_string}" | grep 'IBM XL Fortran'`
  if test "${atp_result}" = ""; then
    atp_result=`echo "${fc_info_string}" | grep 'IBM(R) XL Fortran'`
  fi
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
    if test "${fc_garbage}" -gt 50; then
      AC_DEFINE([FC_IBM],1,[Define to 1 if you are using the IBM XL Fortran compiler.])
      atp_fc_vendor="ibm"
      atp_fc_version="unknown"
      atp_result="yes"
    fi
  else
    AC_DEFINE([FC_IBM],1,[Define to 1 if you are using the IBM XL Fortran compiler.])
    atp_fc_vendor="ibm"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.* V//; s/ .*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_IBM



# _ATP_CHECK_FC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified Fortran compiler is the Intel Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_INTEL],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Intel Fortran compiler])
  fc_info_string=`$1 -v -V 2>&1 | sed -e '/^ifc: warning/d'`
  atp_result=`echo "${fc_info_string}" | grep '^Intel(R) Fortran'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_INTEL],1,[Define to 1 if you are using the Intel Fortran compiler.])
    atp_fc_vendor="intel"
    atp_fc_version=`echo "${fc_info_string}" | grep '^Version' | sed -e 's/Version //;s/ .*//;s/ //g' | head -n 1`
    if test "${atp_fc_version}" = ""; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_INTEL



# _ATP_CHECK_FC_MIPSPRO(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the MIPSpro Fortran
# compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_MIPSPRO],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the MIPSpro Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | sed -e '/^$/d'`
  atp_result=`echo "${fc_info_string}" | grep '^MIPSpro'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_MIPSPRO],1,[Define to 1 if you are using the MIPSpro Fortran compiler.])
    atp_fc_vendor="mipspro"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.*Version //'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_MIPSPRO



# _ATP_CHECK_FC_NAG(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the NAGWare Fortran 95
# compiler. If yes, tries to determine its version number and sets the
# atp_fc_vendor and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_NAG],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the NAGWare Fortran 95 compiler])
  fc_info_string=`$1 -v -V 2>&1`
  atp_result=`echo "${fc_info_string}" | grep '^NAG'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_NAG],1,[Define to 1 if you are using the NAGWare Fortran 95 compiler.])
    atp_fc_vendor="nag"
    atp_fc_version=`echo "${fc_info_string}" | sed -e 's/.*Release //;s/[[( ]].*//'`
    if test "${atp_fc_version}" = ""; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_NAG



# _ATP_CHECK_FC_OPEN64(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the Open64
# Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_OPEN64],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
  fc_info_string=`$1 --version 2>&1`
  atp_result=`echo "${fc_info_string}" | grep '^Open64'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_OPEN64],1,[Define to 1 if you are using the Open64 Fortran compiler.])
    atp_fc_vendor="open64"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_OPEN64



# _ATP_CHECK_FC_PATHSCALE(COMPILER)
# ---------------------------------
#
# Checks whether the specified Fortran compiler is the PathScale
# Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_PATHSCALE],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
  fc_info_string=`$1 -version 2>&1`
  atp_result=`echo "${fc_info_string}" | grep '^PathScale'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_PATHSCALE],1,[Define to 1 if you are using the PathScale Fortran compiler.])
    atp_fc_vendor="pathscale"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_PATHSCALE



# _ATP_CHECK_FC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Portland Group
# Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_PGI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Portland Group Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | sed -e '/^$/d'`
  atp_result=`echo "${fc_info_string}" | grep '^pgf9[[05]]' | grep -v 'No files to process'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_PGI],1,[Define to 1 if you are using the Portland Group Fortran compiler.])
    atp_fc_vendor="pgi"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/^pgf9[[05]] //' | sed -e 's/-.*//'`
    if test "${atp_fc_version}" = "${atp_result}"; then
      atp_fc_version="unknown"
    else
      if test "${atp_fc_version}" = "6.0"; then
                AC_DEFINE([FC_PGI6],1,[Define to 1 if you are using the Portland Group Fortran compiler version 6])
      fi
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_PGI



# _ATP_CHECK_FC_SUN(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Sun WorkShop Fortran compiler.
# If yes, tries to determine its version number and sets the atp_fc_vendor
# and atp_fc_version variables accordingly.
#
AC_DEFUN([_ATP_CHECK_FC_SUN],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Sun Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  atp_result=`echo "${fc_info_string}" | grep 'Sun' | grep 'Fortran 95'`
  if test "${atp_result}" = ""; then
    atp_result="no"
    fc_info_string=""
    atp_fc_vendor="unknown"
    atp_fc_version="unknown"
  else
    AC_DEFINE([FC_SUN],1,[Define to 1 if you are using the Sun Fortran compiler.])
    atp_fc_vendor="sun"
    atp_fc_version=`echo "${atp_result}" | sed -e 's/.* Fortran 95 //;s/ .*//'`
    if test "${atp_fc_version}" = "${atp_result}" -o "${atp_fc_version}" = ""; then
      atp_fc_version="unknown"
    fi
    atp_result="yes"
  fi
  dnl AC_MSG_RESULT(${atp_result})
]) # _ATP_CHECK_FC_SUN



 #############################################################################



# ATP_FC_EXTENSIONS()
# -------------------
#
# Sets the default extensions of Fortran source files and modules,
# whenever possible.
#
AC_DEFUN([ATP_FC_EXTENSIONS],[
  dnl Set Fortran module extension
  AX_F90_MODULE_EXTENSION
  if test "${ax_cv_f90_modext}" != ""; then
    MODEXT="${ax_cv_f90_modext}"
  else
    MODEXT="mod"
    AC_MSG_NOTICE([setting Fortran module extension to ".${MODEXT}"])
  fi
  AC_SUBST(MODEXT)

  dnl Change the default Fortran extension for tests
  AC_FC_SRCEXT(F90,[atp_fc_src_ok="yes"],[atp_fc_src_ok="no"])
  if test "${atp_fc_src_ok}" != "yes"; then
    AC_MSG_WARN([Fortran file extension could not be changed])
    AC_MSG_WARN([some advanced Fortran tests may fail])
  fi
]) # ATP_FC_EXTENSIONS



# ATP_FC_MOD_CASE()
# -----------------
#
# Checks whether the Fortran compiler creates upper-case or lower-case
# module files.
#
AC_DEFUN([ATP_FC_MOD_CASE],[
  AC_REQUIRE([ATP_FC_EXTENSIONS])

  dnl Init
  fc_mod_lowercase="yes"
  fc_mod_uppercase="no"
  AC_MSG_NOTICE([determining Fortran module case])

  dnl Compile a dummy module
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([[
    module conftest
    end module conftest
  ]],[],[AC_MSG_FAILURE([unable to compile a simple Fortran module])])
  AC_LANG_POP([Fortran])

  dnl Check module file existence
  if test -f "CONFTEST.${MODEXT}"; then
    fc_mod_lowercase="no"
    fc_mod_uppercase="yes"
  elif test ! -f "conftest.${MODEXT}"; then
    AC_MSG_WARN([conftest.${MODEXT} Fortran module could not be found])
  fi

  dnl Output final outcome
  AC_MSG_CHECKING([whether Fortran modules are upper-case])
  AC_MSG_RESULT([${fc_mod_uppercase}])
]) # ATP_FC_MOD_CASE



# ATP_PROG_FC()
# -------------
#
# Tries to determine which type of Fortran compiler is installed.
#
AC_DEFUN([ATP_PROG_FC],[
  dnl Init
  atp_fc_vendor="${with_fc_vendor}"
  atp_fc_version="${with_fc_version}"

  if test "${atp_fc_vendor}" = ""; then
    atp_fc_vendor="unknown"
  fi
  if test "${atp_fc_version}" = ""; then
    atp_fc_version="unknown"
  fi
  atp_fc_wrap="no"

  dnl Determine Fortran compiler type (the order is important)
  AC_MSG_CHECKING([which type of Fortran compiler we have])

  dnl Get rid of that one as early as possible
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_IBM(${FC})
  fi

  dnl Should be checked before gfortran
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_INTEL(${FC})
  fi

  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_G95(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_GNU(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_PATHSCALE(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_COMPAQ(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_ABSOFT(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_MIPSPRO(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_OPEN64(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_FUJITSU(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_SUN(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_HITACHI(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_NAG(${FC})
  fi
  if test "${atp_fc_vendor}" = "unknown"; then
    _ATP_CHECK_FC_PGI(${FC})
  fi

  dnl Fall back to generic when detection fails
  if test "${atp_fc_vendor}" = "unknown"; then
    atp_fc_vendor="generic"
  fi

  dnl Normalize Fortran compiler version
  if test "${atp_fc_version}" = "unknown"; then
    atp_fc_version="0.0"
  else
    atp_fc_version=`echo ${atp_fc_version} | cut -d. -f1-2`
  fi

  dnl Display final result
  AC_MSG_RESULT([${atp_fc_vendor} ${atp_fc_version}])

  dnl Schedule compiler info for substitution
  AC_SUBST(atp_fc_vendor)
  AC_SUBST(atp_fc_version)
  AC_SUBST(atp_fc_wrap)
]) # ATP_PROG_FC
