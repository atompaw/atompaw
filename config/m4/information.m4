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
# Architecture information
#



# ATP_INIT_CPU_INFO()
# -------------------
#
# Sets architecture-related variables from the information given by the
# specified target. This is a helper for many other AtomPAW macros, that
# should be called quite early in the configure script.
#
# At present, the variables set are:
#
#  * atp_cpu_model  : CPU model, if guessed;
#  * atp_cpu_64bits : whether the CPU is 64 bits or not.
#
AC_DEFUN([ATP_INIT_CPU_INFO],[
  atp_cpu_platform=`echo "${target}" | cut -d- -f2`
  atp_cpu_vendor=""
  atp_cpu_model=""
  atp_cpu_spec=""
  atp_cpu_bits=""
  atp_cpu_64bits=""

  case "${target}" in

    alpha*)
      atp_cpu_vendor="dec"
      atp_cpu_model="${target_cpu}"
      atp_cpu_64bits=`echo "${atp_cpu_model}" | grep '64$'`
      if test "${atp_cpu_64bits}" = ""; then
        atp_cpu_64bits="no"
      else
        atp_cpu_64bits="yes"
      fi
      ;;

    powerpc*)
      atp_cpu_vendor="ibm"
      atp_cpu_model="${target_cpu}"
      atp_cpu_64bits=`echo "${atp_cpu_model}" | grep '64$'`
      if test "${atp_cpu_64bits}" = ""; then
        atp_cpu_64bits="no"
      else
        atp_cpu_64bits="yes"
      fi
      ;;

    i686-*linux*)
      dnl Athlon ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Athlon'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="amd"
          atp_cpu_model="athlon"
          atp_cpu_64bits="no"
        fi
      fi
      dnl Pentium 3 ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Pentium III'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="pentium3"
          atp_cpu_64bits="no"
        fi
      fi
      dnl Pentium 4 ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Pentium(R) 4'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="pentium4"
          atp_cpu_64bits="no"
        fi
      fi
      dnl Pentium 4M ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Pentium(R) M'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="pentium4"
          atp_cpu_64bits="no"
        fi
      fi
      dnl Centrino ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) CPU           T2400'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="centrino"
          atp_cpu_64bits="no"
        fi
      fi
      dnl Pentium CoreDuo ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) CPU           T2050'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="coreduo"
          atp_cpu_64bits="no"
        fi
      fi
      dnl Pentium Core2 ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Core(TM)2 CPU'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="core2"
          atp_cpu_64bits="no"
        fi
      fi
      dnl Unknown
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_vendor="unknown"
        atp_cpu_model="unknown"
        atp_cpu_64bits="unknown"
      fi
      ;;

    ia64-*linux*)
      dnl Itanium 1 ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Itanium' | grep -v 'Itanium 2'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="itanium1"
        fi
      fi
      dnl Itanium 2 ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Itanium 2'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="itanium2"
        fi
      fi
      dnl Unknown
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_vendor="unknown"
        atp_cpu_model="unknown"
      fi
      dnl The processor is anyway 64-bit
      atp_cpu_64bits="yes"
      ;;

    mips*irix*)
      # Get processor type
      atp_cpu_vendor="mips"
      atp_cpu_model=`hinv 2> /dev/null | grep '^CPU: MIPS '`
      if test "${atp_cpu_model}" != ""; then
        atp_cpu_model=`echo "${atp_cpu_model}" | awk '{print tolower($3)}'`
      fi
      atp_cpu_64bits="yes"
      ;;

    x86_64-*linux*)
      dnl Athlon64 ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Athlon'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="amd"
          atp_cpu_model="athlon64"
        fi
      fi
      dnl Opteron ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Opteron'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="amd"
          atp_cpu_model="opteron"
        fi
      fi
      dnl Sempron ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Sempron'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="amd"
          atp_cpu_model="athlon64"
        fi
      fi
      dnl Xeon ?
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) XEON(TM)'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="xeon"
        fi
      fi
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Xeon(R)'`
        if test "${atp_cpu_model}" != ""; then
          atp_cpu_vendor="intel"
          atp_cpu_model="xeon"
        fi
      fi
      dnl Unknown
      if test "${atp_cpu_model}" = ""; then
        atp_cpu_vendor="unknown"
        atp_cpu_model="unknown"
      fi
      dnl The processor is anyway 64-bit
      atp_cpu_64bits="yes"
      ;;

  esac

  dnl Generate CPU identifier
  atp_cpu_spec="${atp_cpu_vendor}_${atp_cpu_model}"

  AC_SUBST(atp_cpu_platform)
  AC_SUBST(atp_cpu_vendor)
  AC_SUBST(atp_cpu_model)
  AC_SUBST(atp_cpu_spec)
  AC_SUBST(atp_cpu_64bits)
  AC_SUBST(atp_cpu_bits)
]) # ATP_INIT_CPU_INFO



# ATP_INIT_OS_INFO()
# ------------------
#
# Sets OS-related variables from the information given by the specified
# target.
#
AC_DEFUN([ATP_INIT_OS_INFO],[
  case "${target_os}" in

    *linux*)
      AC_DEFINE([HAVE_OS_LINUX],1,[Define to 1 if you are using Linux.])
      ;;

    *apple*)
      AC_DEFINE([HAVE_OS_MACOSX],1,[Define to 1 if you are using MacOS X.])
      ;;

  esac
]) # ATP_INIT_OS_INFO



# ATP_INIT_HEADER()
# -----------------
#
# Initializes the contents of the header file produced by Autoheader.
#
AC_DEFUN([ATP_INIT_HEADER],[
  dnl Set top of file ...
  AH_TOP([/*
  * Copyright (C) 2010 Yann Pouillon
  *
  * This file is part of the AtomPAW software package. For license information,
  * please see the COPYING file in the top-level directory of the source
  * distribution.
  *
  */

/* AtomPAW configuration */

#ifndef _ATOMPAW_CONFIG_H
#define _ATOMPAW_CONFIG_H

#ifdef __INTEL_COMPILER
#define FC_INTEL 1
#endif

])

  dnl ... as well as bottom
  AH_BOTTOM([#endif /* _ATOMPAW_CONFIG_H */])
]) # ATP_INIT_HEADER
