#!/usr/bin/python

import re
import os
import sys

cfg_opts = dict()

new_opts = dict()
new_opts["F90"] = "FC"
new_opts["F90FLAGS"] = "FCFLAGS"
new_opts["LDFLAGS"] = "LDFLAGS"
new_opts["LIBDIR"] = "ATOMPAW_LIBDIR"
new_opts["LIBS"] = "LIBS"
new_opts["LAPACK_BLAS_LIBS"] = "LINALG_LIBS"
new_opts["BLACS_LIBS"] = "BLACS_LIBS"
new_opts["SCALAPACK_LIBS"] = "SCALAPACK_LIBS"

cfg_template = """#!/bin/sh
#
# Configuration script for @CFG@
#

# Allow for VPATH builds
if test -x "./configure"; then
  CFG_SCRIPT="./configure"
elif test -x "../configure"; then
  CFG_SCRIPT="../configure"
else
  echo "Error: could not find configure script"
  exit 1
fi

@LIB@# Run configure
${CFG_SCRIPT}"""

for old_file in sys.argv[1:]:
  cfg = re.sub("^make\.","",old_file)
  cfg_opts[cfg] = dict()
  opt = ""
  cnt_line = False
  cnt_line_old = False

  for line in file(old_file,"r").readlines():
    line = line.strip()
    line = re.sub("#.*","",line)
    if ( re.search("\\\\$",line) ):
      cnt_line = True
    else:
      cnt_line = False
    if ( (re.match("[^ ]",line)) and (not cnt_line_old) ):
      line = line.split("=")
    else:
      line = [line]

    if ( len(line) > 0 ):
      if ( (len(line) == 1) or cnt_line_old ):
        if ( len(line[0]) > 0 ):
          cfg_opts[cfg][opt] += "\n    " + line[0].strip()
      else:
        if (len(line) > 2 ):
          line[1] = "=".join(line[1:])
        opt = new_opts[line[0].strip()]
        if ( opt in cfg_opts[cfg] ):
          sys.stderr.write("WARNING: multiple definitions of %s\n" % \
            (line[0].strip()))
        line[1] = line[1].strip()
        if ( re.match("^\\$\\(",line[1]) ):
          line_buf = re.sub("\\$\\(([^\\)]*)\\)[^\\$]*","\\1 ",line[1])
          if ( line_buf == line[1] ):
            line_buf = re.sub("\\$\\(([^\\)]*)\\)","\\1",line[1])
          line_tmp = ""
          for var in line_buf.strip().split():
            line_tmp += " ${%s}" % (new_opts[var])
          line[1] = line_tmp.lstrip()
        cfg_opts[cfg][opt] = line[1]

    cnt_line_old = cnt_line

for cfg_key in cfg_opts:
  cfg = cfg_opts[cfg_key]

  la_libs = ""
  if ( "LINALG_LIBS" in cfg ):
    la_libs += " \nLINALG_LIBS=\"%s\"" % (cfg["LINALG_LIBS"])
  if ( "BLACS_LIBS" in cfg ):
    la_libs += " \nBLACS_LIBS=\"%s\"" % (cfg["BLACS_LIBS"])
  if ( "SCALAPACK_LIBS" in cfg ):
    la_libs += " \nSCALAPACK_LIBS=\"%s\"" % (cfg["SCALAPACK_LIBS"])
  if ( la_libs != "" ):
    la_libs = "# Libraries" + la_libs + "\n\n"

  la_opts = ""
  if ( "SCALAPACK_LIBS" in cfg ):
    la_opts += " ${SCALAPACK_LIBS}"
  if ( "BLACS_LIBS" in cfg ):
    la_opts += " ${BLACS_LIBS}"
  if ( "LINALG_LIBS" in cfg ):
    la_opts += " ${LINALG_LIBS}"
  if ( la_opts != "" ):
    la_opts = " \\\n  --with-linalg-libs=\"%s\"" % (la_opts.lstrip())

  fc_opts = ""
  if ( "FC" in cfg ):
    fc_opts += " \\\n  FC=\"%s\"" % (cfg["FC"])
  if ( "FCFLAGS" in cfg ):
    fc_opts += " \\\n  FCFLAGS=\"%s\"" % (cfg["FCFLAGS"])
  if ( "LDFLAGS" in cfg ):
    fc_opts += " \\\n  LDFLAGS=\"%s\"" % (cfg["LDFLAGS"])
  if ( "ATOMPAW_LIBDIR" in cfg ):
    fc_opts += " \\\n  ATOMPAW_LIBDIR=\"%s\"" % (cfg["ATOMPAW_LIBDIR"])
  if ( ("LIBS" in cfg) and (la_opts == "") ):
    fc_opts += " \\\n  LIBS=\"%s\"" % (cfg["LIBS"])

  cfg_data = re.sub("@CFG@",cfg_key,cfg_template)
  cfg_data = re.sub("@LIB@",la_libs,cfg_data)

  cfg_data += la_opts + fc_opts + "\n"

  file("../examples/configure.%s" % (cfg_key),"w").write(cfg_data)

