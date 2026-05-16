# ATOMPAW - LibXC library detection
# ---
# Define target  : Libxc::xc
# Define variable: Libxc_LIBRARIES
# Define variable: Libxc_INCLUDE_DIRS
# Define variable: Libxc_VERSION
# ---
# Define variable: USE_LIBXC if LibXC is to be used
# Define variable: LIBXC_OK  if LibXC is usable (or not)
# ---
# Define property: ENABLE_LIBXC (default=AUTO)
# ---
# Note: may use PkgConfig_FOUND variable
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
  if (NOT DEFINED Libxc_ROOT AND DEFINED LIBXC_ROOT)
    set(Libxc_ROOT ${LIBXC_ROOT})
  endif()
#  find_package(Libxc QUIET)
  if (Libxc_FOUND OR LIBXC_FOUND)
    if (TARGET Libxc::xc)
      message(STATUS "Libxc found via find_package (config mode)")
      set(LIBXC_FOUND TRUE)
      set(Libxc_FOUND_with_cmake TRUE)
    elseif (LIBXC_INCLUDE_DIRS OR LIBXC_LIBRARIES)
      message(STATUS "Libxc found via find_package (module mode)")
      set(LIBXC_FOUND TRUE)
      set(Libxc_FOUND_with_cmake TRUE)
      if (NOT TARGET Libxc::xc)
        add_library(Libxc::xc UNKNOWN IMPORTED)
        set_target_properties(Libxc::xc PROPERTIES
          IMPORTED_LOCATION             "${LIBXC_LIBRARIES}"
          INTERFACE_LINK_LIBRARIES      "${LIBXC_LIBRARIES}"
          INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIRS}")
      endif()
    endif()
  endif()

# 2- Try via PKGCONFIG
  if (NOT LIBXC_FOUND AND NOT Libxc_FOUND)
    if (PkgConfig_FOUND)
      if (DEFINED LIBXC_ROOT)
        set(ENV{PKG_CONFIG_PATH} "${LIBXC_ROOT}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
      elseif (Libxc_ROOT)
        set(ENV{PKG_CONFIG_PATH} "${Libxc_ROOT}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
      endif()
      pkg_check_modules(LIBXC QUIET IMPORTED_TARGET libxc)
    endif()

    if (LIBXC_FOUND)
      message(STATUS "Libxc found via pkg-config: ${LIBXC_VERSION}")

      # Create target
      if (NOT TARGET Libxc::xc)
        add_library(Libxc::xc INTERFACE IMPORTED GLOBAL)
        if (TARGET PkgConfig::LIBXC)
          target_link_libraries(Libxc::xc INTERFACE PkgConfig::LIBXC)
        endif()
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
                  PATH_SUFFIXES libxc
                  NO_CACHE)
      endif()
      if (NOT Libxc_INCLUDE_DIRS)
        set(LIBXC_FOUND FALSE)
      else()
        set(Libxc_FOUND_with_pkgconfig TRUE)
      endif()
    endif()
  endif()

# 3- Try via environment variables
  if (NOT LIBXC_FOUND AND NOT Libxc_FOUND)
    foreach(_var LIBXC_ROOT LIBXC_PREFIX LIBXC_DIR)
      if (DEFINED ENV{${_var}} AND NOT LIBXC_ROOT)
        set(LIBXC_ROOT "$ENV{${_var}}" CACHE PATH "LIBXC root")
      endif()
    endforeach()
    string(REPLACE ":" ";" _ld_paths "$ENV{LD_LIBRARY_PATH}")
    foreach(_path ${_ld_paths})
      list(APPEND _ld_path_incs "${_path}/../include")
     endforeach()
    find_path(LIBXC_INCLUDE_DIR
              NAMES xc.h xc_funcs.h
              HINTS
                "${LIBXC_ROOT}/include"
                "${_ld_path_incs}"
              PATH_SUFFIXES "" libxc xc
              NO_DEFAULT_PATH
              NO_CACHE)
    find_library(LIBXC_LIBRARY
                 NAMES xc
                 HINTS
                   "${LIBXC_ROOT}/lib"
                   "${LIBXC_ROOT}/lib64"
                   ${_ld_paths}
                 NO_DEFAULT_PATH
                 NO_CACHE)
    if (LIBXC_INCLUDE_DIR AND LIBXC_LIBRARY AND NOT TARGET Libxc::xc)
      message(STATUS "Libxc found in environment: ${LIBXC_LIBRARY}")
      add_library(Libxc::xc UNKNOWN IMPORTED)
      set_target_properties(Libxc::xc PROPERTIES
          IMPORTED_LOCATION             "${LIBXC_LIBRARY}"
          INTERFACE_LINK_LIBRARIES      "${LIBXC_LIBRARY}"
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
#include <xc.h>
#include <xc_funcs.h>
#include <xc_version.h>
    int main() {
      xc_func_type func;
      double rho[1] = {1.0};
      double ex[1], vx[1];
      int version_major;
      version_major = (int)XC_MAJOR_VERSION;
      xc_func_init(&func, XC_LDA_X, XC_UNPOLARIZED);
      xc_lda_exc_vxc(&func, 1, rho, ex, vx);
      xc_func_end(&func);
}
")
    string(REPLACE ";" "\\;" _libxc_incdirs_escaped "${Libxc_INCLUDE_DIRS}")
    try_run(LIBXC_RUN_RESULT LIBXC_COMPILE_RESULT
            ${CMAKE_BINARY_DIR}/tests/test_libxc
            ${CMAKE_BINARY_DIR}/tests/test_libxc/test_libxc.c
            CMAKE_FLAGS  "-DINCLUDE_DIRECTORIES:STRING=${_libxc_incdirs_escaped}"
            LINK_LIBRARIES ${Libxc_LIBRARIES}
            OUTPUT_VARIABLE TRY_OUTPUT)
    if (LIBXC_COMPILE_RESULT AND NOT LIBXC_RUN_RESULT)
      set(LIBXC_OK TRUE)
    else()
      set(LIBXC_OK FALSE)
      message(STATUS "LibXC compilation+execution test failed!")
      if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(FATAL_ERROR "${TRY_OUTPUT}")
      endif()
    endif()
  endif()

# If not found and required then fatal error
  if  (ENABLE_LIBXC STREQUAL "ON" AND NOT LIBXC_OK)
    message(FATAL_ERROR "Libxc required (ENABLE_LIBXC=ON) but not available!")
  endif()

  if (LIBXC_OK)
    message(CHECK_PASS "done")
  else()
    message(CHECK_FAIL "not found")
  endif()
  
endif()
