# -*- Autoconf -*-
#
# M4 macros for AtomPAW
#
# Copyright (C) 2011 Marc Torrent
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Generate the libxc_names.in file from libXC includes
#


# ATP_MK_LIBXC_NAMES()
# --------------------
#
AC_DEFUN([ATP_MK_LIBXC_NAMES],[

XC_FUNCS_H="xc_funcs.h"

#Search for xc_funcs.h file...
#-----------------------------
#... in libxc_incs:
LIBXC_NAMES_H=""
if test "${with_libxc_incs}" != ""; then
  for ff in $(echo ${with_libxc_incs} | $SED -n "s/\-I/\ /gp") ;do
    if test -f "${ff}/$XC_FUNCS_H";then
      LIBXC_NAMES_H=${ff}"/"${XC_FUNCS_H}
    fi
  done
fi
#... in $INCLUDE:
if test "${LIBXC_NAMES_H}" == ""; then
  if test "${INCLUDE}" != "";then
    for ff in $(echo ${INCLUDE} | $SED -n "s/\:/\ /gp") ;do
      if test -f "${ff}/$XC_FUNCS_H";then
        LIBXC_NAMES_H=${ff}"/"${XC_FUNCS_H}
      fi
    done
  fi
fi
#... in $CPPFLAGS:
if test "${LIBXC_NAMES_H}" == ""; then
  if test "${CPPFLAGS}" != "";then
    for ff in $(echo ${CPPFLAGS} | $SED -n "s/\-I/\ /gp") ;do
      if test -f "${ff}/$XC_FUNCS_H";then
        LIBXC_NAMES_H=${ff}"/"${XC_FUNCS_H}
      fi
    done
  fi
fi
#... in $FCFLAGS:
if test "${LIBXC_NAMES_H}" == ""; then
  if test "${FCFLAGS}" != "";then
    for ff in $(echo ${FCFLAGS} | $SED -n "s/\-I/\ /gp") ;do
      if test -f "${ff}/$XC_FUNCS_H";then
        LIBXC_NAMES_H=${ff}"/"${XC_FUNCS_H}
      fi
    done
  fi
fi
#... in LibXC default:
if test "${LIBXC_NAMES_H}" == ""; then
  if test -f "/opt/etsf/include/$XC_FUNCS_H";then
    LIBXC_NAMES_H="/opt/etsf/include/"${XC_FUNCS_H}
  fi
fi
#...not found:
if test "${LIBXC_NAMES_H}" == ""; then
  AC_MSG_ERROR([could not find "${XC_FUNCS_H}" file])
fi

#Build src/libxc_id.in file
#--------------------------
LIBXC_NAMES_IN=${ac_abs_top_srcdir}"/src/libxc_id.in"

cat <<EOF > ${LIBXC_NAMES_IN}
!!=================================================================
!! NAME    : libxc_id
!! FUNCTION: From a character string, gives the libXC id
!!=================================================================
 integer function libxc_id(xcname)
!------------------------------------------------------------------
 implicit none
 character*(*),intent(in) :: xcname
!------------------------------------------------------------------
 character*50 :: xcstrg
!------------------------------------------------------------------
 libxc_id=0
 xcstrg=trim(xcname)
 if (xcstrg=="") then
   libxc_id=-1
   return
 end if
EOF

if test -f "${LIBXC_NAMES_H}";then
  $SED -n \
  's/.*\#define\ *\([_A-Za-z0-9ÀÁÂÃÄÅÇÑñÇçÈÉÊËÌÍÎÏÒÓÔÕÖØÙÚÛÜÝÑàáâãäåçèéêëìíîïðòóôõöøùúûüýÿñ]*\).*/\ if\ (xcstrg==\"\1\") libxc_id=\1/p' \
  ${LIBXC_NAMES_H} >> ${LIBXC_NAMES_IN}
fi

cat <<EOF >> ${LIBXC_NAMES_IN}
 end function libxc_id

EOF

]) # ATP_MK_LIBXC_NAMES
