#!/bin/sh
#
# Copyright (C) 2010 Yann Pouillon
#
# This file is part of the AtomPAW software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "src/atompaw_prog.F90"; then
  echo "This is not a AtomPAW source tree - aborting now"
  exit 1
fi

# Create possibly missing directories
mkdir -p config/gnu

# Generate M4 includes
echo "Generating aclocal.m4..."
aclocal -I config/m4
echo "done."

# Generate configure auxiliary files
echo "Generating config.h.in..."
autoheader
echo "done."

# Generate configure
echo "Generating configure script..."
autoconf
echo "done."

# Generate libtool scripts
echo "Generating libtool scripts..."
set +e
my_libtoolize="libtoolize"
${my_libtoolize} --version >/dev/null 2>&1
if test "${?}" != "0"; then 
  my_libtoolize="glibtoolize"
fi
${my_libtoolize} --version >/dev/null 2>&1
if test "${?}" != "0"; then 
  echo "Error: could not find a working version of libtoolize" >&2
  exit 1
fi
set -e
${my_libtoolize} --automake --copy --force
echo "done."

# Generate makefile inputs
# Do not use "automake --force-missing", as it overwrites the INSTALL file.
echo "Generating Makefile.in for each directory..."
automake --add-missing --copy
echo "done."
