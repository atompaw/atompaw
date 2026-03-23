#ifndef ATOMPAW_CONFIG_H
#define ATOMPAW_CONFIG_H 1

/* config.h.  Generated from config.h.cmake by cmake.  */

/* Define to 1 according to the detected compiler. */
#cmakedefine FC_ABSOFT 1
#cmakedefine FC_ARM 1
#cmakedefine FC_CRAY 1
#cmakedefine FC_FLANG 1
#cmakedefine FC_G95 1
#cmakedefine FC_GNU 1
#cmakedefine FC_INTEL 1
#cmakedefine FC_INTEL_LLVM 1
#cmakedefine FC_NAG 1
#cmakedefine FC_NVHPC 1
#cmakedefine FC_PGI 1
#cmakedefine FC_PATHSCALE 1
#cmakedefine FC_XL 1
#cmakedefine FC_XL_CLANG 1
#cmakedefine FC_IBM 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H 1

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H 1

/* Define to 1 if you have the <cstdint> header file. */
#cmakedefine HAVE_CSTDINT 1

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H 1

/* Define to 1 if iso_c_binding is available */
#cmakedefine HAVE_FC_ISO_C_BINDING 1

/* Define to 1 if libxc is available */
#cmakedefine HAVE_LIBXC 1

/* Name of package */
#define PACKAGE "@PROJECT_NAME@"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "stephane.schanne@cea.fr;frederic.chateau@cea.fr"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@PROJECT_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@PROJECT_NAME@ @PROJECT_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@PROJECT_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@HOMEPAGE_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@PROJECT_VERSION@"

/* Define to the description of this package. */
#define PACKAGE_DESCRIPTION "@PROJECT_DESCRIPTION@"

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS 1

/* Version number of package */
#define VERSION "@PROJECT_VERSION@"

/* once: ATOMPAW_CONFIG_H */
#endif
