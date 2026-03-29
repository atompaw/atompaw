# ATOMPAW - Linear Algebra library detection
# ---
# Define target  : BLAS::BLAS
# Define variable: BLAS_LIBRARIES
# Define target  : LAPACK::LAPACK
# Define variable: LAPACK_INCLUDE_DIRS

message(CHECK_START ">>> Detecting BLAS/LAPACK")

# Store LD_LIBRARY_PATH as a list
string(REPLACE ":" ";" _ld_paths "$ENV{LD_LIBRARY_PATH}")

# ----- BLAS -----

# 1- Try BLAS via CMake find_package
find_package(BLAS QUIET)
if (BLAS_FOUND)
  message(STATUS "BLAS found via cmake target")
  set(BLAS_FOUND_with_cmake TRUE)

# 2- Try BLAS via PKGCONFIG
else()
  if (PkgConfig_FOUND)
#    pkg_check_modules(BLAS QUIET IMPORTED_TARGET blas openblas)
  else()
    set(BLAS_FOUND FALSE)
  endif()
  if (BLAS_FOUND)
    message(STATUS "BLAS found via pkg-config")
    if (NOT TARGET BLAS::BLAS)
      add_library(BLAS::BLAS ALIAS PkgConfig::BLAS)
    endif()
    get_target_property(BLAS_LIBRARIES BLAS::BLAS INTERFACE_LINK_LIBRARIES)
    set(BLAS_FOUND_with_pkgconfig TRUE)

# 3- Try BLAS via environment variables
  else()
    if (NOT BLAS_ROOT)
      set(BLAS_ROOT "$ENV{BLAS_ROOT}" CACHE PATH "BLAS root")
    endif()
    if (NOT BLAS_ROOT)
      set(BLAS_ROOT "$ENV{BLAS_DIR}" CACHE PATH "BLAS root")
    endif()
    find_library(BLAS_LIBRARY
                 NAMES blas openblas
                 HINTS
                   "${BLAS_ROOT}/lib"
                   "${BLAS_ROOT}/lib64"
                   "${_ld_paths}"
                 NO_DEFAULT_PATH)
    if (BLAS_LIBRARY AND NOT TARGET BLAS::BLAS)
      message(STATUS "BLAS found in environment: ${BLAS_LIBRARY}")
      add_library(BLAS::BLAS UNKNOWN IMPORTED)
      set_target_properties(BLAS::BLAS PROPERTIES
          IMPORTED_LOCATION "${BLAS_LIBRARY}")
      set(BLAS_LIBRARIES ${BLAS_LIBRARY})
      set(BLAS_FOUND_with_environment TRUE)
    endif()
  endif()
endif()

# ----- LAPACK -----

# 1- Try LAPACK via CMake find_package
find_package(LAPACK QUIET)
if (LAPACK_FOUND)
  message(STATUS "LAPACK found via cmake target")
  set(LAPACK_FOUND_with_cmake TRUE)

# 2- Try LAPACK via PKGCONFIG
else()
  if (PkgConfig_FOUND)
    pkg_check_modules(LAPACK QUIET IMPORTED_TARGET lapack openblas)
  else()
    set(LAPACK_FOUND FALSE)
  endif()
  if (LAPACK_FOUND)
    message(STATUS "LAPACK found via pkg-config")
    if (NOT TARGET LAPACK::LAPACK)
      add_library(LAPACK::LAPACK ALIAS PkgConfig::LAPACK)
    endif()
    get_target_property(LAPACK_LIBRARIES Libxc::xc INTERFACE_LINK_LIBRARIES)
    set(LAPACK_FOUND_with_pkgconfig TRUE)

# 3- Try LAPACK via environment variables
  else()
    if (NOT LAPACK_ROOT)
      set(LAPACK_ROOT "$ENV{LAPACK_ROOT}" CACHE PATH "LAPACK root")
    endif()
    if (NOT LAPACK_ROOT)
      set(LAPACK_ROOT "$ENV{LAPACK_DIR}" CACHE PATH "LAPACK root")
    endif()
    find_library(LAPACK_LIBRARY
                 NAMES lapack openblas
                 HINTS
                   "${LAPACK_ROOT}/lib"
                   "${LAPACK_ROOT}/lib64"
                   "${_ld_paths}"
                 NO_DEFAULT_PATH)
    if (LAPACK_LIBRARY AND NOT TARGET LAPACK::LAPACK)
      message(STATUS "LAPACK found in environment: ${LAPACK_LIBRARY}")
      add_library(LAPACK::LAPACK UNKNOWN IMPORTED)
      set_target_properties(LAPACK::LAPACK PROPERTIES
          IMPORTED_LOCATION "${LAPACK_LIBRARY}")
      set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
      set(LAPACK_FOUND_with_environment TRUE)
    endif()
  endif()
endif()

# ----- TESTS -----

# Did we find BLAS/LAPACK?
if (BLAS_FOUND_with_cmake OR BLAS_FOUND_with_pkgconfig OR
    BLAS_FOUND_with_environment)
  if (LAPACK_FOUND_with_cmake OR LAPACK_FOUND_with_pkgconfig OR
      LAPACK_FOUND_with_environment)
    set(BLAS_LAPACK_FOUND_OK TRUE)
  else()
    set(BLAS_LAPACK_FOUND_OK FALSE)
  endif()
else()
  set(BLAS_LAPACK_FOUND_OK FALSE)
endif()
  
# If BLAS/LAPACK found, test if it works
if (BLAS_LAPACK_FOUND_OK)
  enable_language(Fortran)
  try_run(LINALG_RUN_RESULT LINALG_COMPILE_RESULT
          ${CMAKE_BINARY_DIR}/test_linalg
          ${CMAKE_SOURCE_DIR}/cmake/tests/test_linalg.F90
          LINK_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  if (NOT LINALG_RUN_RESULT)
    set(BLAS_LAPACK_OK TRUE)
  else()
    set(BLAS_LAPACK_OK FALSE)
  endif()
else()
  message(FATAL_ERROR "BLAS+LAPACK detection failed!")
  set(BLAS_LAPACK_OK FALSE)
endif()

if(BLAS_LAPACK_OK)
  message(STATUS "BLAS+LAPACK Fortran OK")
else()
  message(FATAL_ERROR "BLAS+LAPACK test failed: ${TRY_OUTPUT}")
endif()

if (BLAS_LAPACK_OK)
  message(CHECK_PASS "done")
else()
  message(CHECK_PASS "not found")
endif()
