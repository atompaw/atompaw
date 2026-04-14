# ATOMPAW - LibXC library detection
# ---
# Define target  : Libxc::xc
# Define variable: Libxc_LIBRARIES
# Define variable: Libxc_INCLUDE_DIRS
# Define variable: Libxc_VERSION
# ---

set(ENABLE_LIBXC "AUTO" CACHE STRING "Use libxc: AUTO (automatic), ON (force activation), OFF (deactivated)")
set_property(CACHE ENABLE_LIBXC PROPERTY STRINGS "AUTO;ON;OFF")

set(USE_LIBXC FALSE)
if(ENABLE_LIBXC STREQUAL "AUTO" OR ENABLE_LIBXC STREQUAL "ON")
  set(USE_LIBXC TRUE)
endif()

if (USE_LIBXC)
  message(CHECK_START ">>> Detecting libXC")

# 1- Try via CMake find_package
  find_package(Libxc QUIET)
  if (Libxc_FOUND)
    message(STATUS "Libxc found via cmake target")
    set(LIBXC_FOUND TRUE)
    set(Libxc_FOUND_with_cmake TRUE)
  endif()

# 2- Try via PKGCONFIG
  if (NOT LIBXC_FOUND)
    if (PkgConfig_FOUND)
      pkg_check_modules(LIBXC QUIET IMPORTED_TARGET libxc)
    endif()
    if (LIBXC_FOUND)
      message(STATUS "Libxc found via pkg-config: ${LIBXC_VERSION}")
      if (NOT TARGET Libxc::xc)
        add_library(Libxc::xc ALIAS PkgConfig::LIBXC)
      endif()
      set(Libxc_VERSION ${LIBXC_VERSION})
      get_target_property(Libxc_LIBRARIES Libxc::xc INTERFACE_LINK_LIBRARIES)
      get_target_property(Libxc_INCLUDE_DIRS Libxc::xc INTERFACE_INCLUDE_DIRECTORIES)
      if (NOT Libxc_INCLUDE_DIRS OR Libxc_INCLUDE_DIRS STREQUAL "Libxc_INCLUDE_DIRS-NOTFOUND")
        set(Libxc_INCLUDE_DIRS ${LIBXC_INCLUDE_DIRS})
      endif()
      if (NOT Libxc_INCLUDE_DIRS)
        pkg_get_variable(Libxc_INCLUDE_DIRS libxc includedir)
      endif()
      if (NOT Libxc_INCLUDE_DIRS)
        find_path(Libxc_INCLUDE_DIRS NAMES xc.h
                  HINTS ${LIBXC_PREFIX}/include
                  ${LIBXC_LIBRARY_DIRS}/../include
                  PATH_SUFFIXES libxc)
      endif()
      if (NOT Libxc_INCLUDE_DIRS)
        set(LIBXC_FOUND FALSE)
      else()
        set(Libxc_FOUND_with_pkgconfig TRUE)
      endif()
    endif()
  endif()

# 3- Try via environment variables
  if (NOT LIBXC_FOUND)
    if (NOT LIBXC_ROOT)
      set(LIBXC_ROOT "$ENV{LIBXC_ROOT}" CACHE PATH "Libxc root")
    endif()
    if (NOT LIBXC_ROOT)
      set(LIBXC_ROOT "$ENV{LIBXC_DIR}" CACHE PATH "Libxc root")
    endif()
    string(REPLACE ":" ";" _ld_paths "$ENV{LD_LIBRARY_PATH}")
    find_path(LIBXC_INCLUDE_DIR
              NAMES xc.h xc_funcs.h
              HINTS
                "${LIBXC_ROOT}/include"
                "${_ld_paths}/../include"
              PATH_SUFFIXES "" libxc xc
              NO_DEFAULT_PATH)
    find_library(LIBXC_LIBRARY
                 NAMES xc
                 HINTS
                   "${LIBXC_ROOT}/lib"
                   "${LIBXC_ROOT}/lib64"
                   ${_ld_paths}
                 NO_DEFAULT_PATH)
    if (LIBXC_INCLUDE_DIR AND LIBXC_LIBRARY AND NOT TARGET Libxc::xc)
      message(STATUS "Libxc found in environment: ${LIBXC_LIBRARY}")
      add_library(Libxc::xc UNKNOWN IMPORTED)
      set_target_properties(Libxc::xc PROPERTIES
          IMPORTED_LOCATION "${LIBXC_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIR}")
      set(Libxc_INCLUDE_DIRS ${LIBXC_INCLUDE_DIR})
      set(Libxc_LIBRARIES ${LIBXC_LIBRARY})

      # Retrieve Libxc_VERSION
      file(STRINGS "${LIBXC_INCLUDE_DIR}/xc_version.h" _xc_version_line
           REGEX "^#define[ \t]+XC_VERSION[ \t]+\"[0-9.]+\"")
      string(REGEX REPLACE ".*XC_VERSION[ \t]+\"([0-9.]+)\".*" "\\1"
             Libxc_VERSION "${_xc_version_line}")
      if (NOT Libxc_VERSION)
        set(Libxc_VERSION "undefined")
      endif()

      set(Libxc_FOUND_with_environment TRUE)
    endif()
  endif()

# Did we find libXC?
  if (Libxc_FOUND_with_cmake OR Libxc_FOUND_with_pkgconfig OR
      Libxc_FOUND_with_environment)
    set(LIBXC_FOUND_OK TRUE)
  else()
    set(LIBXC_FOUND_OK FALSE)
    set(LIBXC_OK FALSE)
  endif()

# If libxc found, test if it works
  if (LIBXC_FOUND_OK)
    # C program calling LibXC
    file(WRITE ${CMAKE_BINARY_DIR}/tests/test_libxc/test_libxc.c "
#include "xc.h"
int main() {
    xc_func_type func;
    double rho[1] = {1.0};
    double ex[1], vx[1];
    int version_major;
    version_major = (int)XC_VERSION_MAJOR;
    xc_func_init(&func, XC_FUNC_TYPE_LDA, XC_POLARIZATION_NONE, 1, "LDA_X");
    xc_lda(&func, 1, rho, ex, vx);
    xc_func_end(&func);
}
")
    try_run(LIBXC_RUN_RESULT LIBXC_COMPILE_RESULT
            ${CMAKE_BINARY_DIR}/tests/test_libxc
            ${CMAKE_BINARY_DIR}/tests/test_libxc/test_libxc.c
            CMAKE_FLAGS  -DINCLUDE_DIRECTORIES=${Libxc_INCLUDE_DIRS}
            LINK_LIBRARIES ${Libxc_LIBRARIES})
    if (NOT LIBXC_RUN_RESULT)
      set(LIBXC_OK TRUE)
    else()
      set(LIBXC_OK FALSE)
    endif()
  endif()

# If not found and required then fatal error
  if  (ENABLE_LIBXC STREQUAL "ON" AND NOT LIBXC_OK)
    message(FATAL_ERROR "Libxc required (ENABLE_LIBXC=ON) but not available!")
  endif()

  if (LIBXC_OK)
    message(CHECK_PASS "done")
  else()
    message(CHECK_PASS "not found")
  endif()
  
endif()
