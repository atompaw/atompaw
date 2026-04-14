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
endif()

# 2- Try BLAS via PKGCONFIG
if (NOT BLAS_FOUND)
  if (PkgConfig_FOUND)
    pkg_check_modules(BLAS QUIET IMPORTED_TARGET blas openblas)
  endif()
  if (BLAS_FOUND)
    message(STATUS "BLAS found via pkg-config")
    if (NOT TARGET BLAS::BLAS)
      add_library(BLAS::BLAS ALIAS PkgConfig::BLAS)
    endif()
    get_target_property(BLAS_LIBRARIES BLAS::BLAS INTERFACE_LINK_LIBRARIES)
    if (NOT BLAS_LIBRARIES)
      set(BLAS_FOUND FALSE)
    else()
      set(BLAS_FOUND_with_pkgconfig TRUE)
    endif()
  endif()
endif()

# 3- Try BLAS via environment variables
if (NOT BLAS_FOUND)
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
    set_target_properties(BLAS::BLAS PROPERTIES IMPORTED_LOCATION "${BLAS_LIBRARY}")
    set(BLAS_LIBRARIES ${BLAS_LIBRARY})
    set(BLAS_FOUND_with_environment TRUE)
  endif()
endif()

# ----- LAPACK -----

# 1- Try LAPACK via CMake find_package
find_package(LAPACK QUIET)
if (LAPACK_FOUND)
  message(STATUS "LAPACK found via cmake target")
  set(LAPACK_FOUND_with_cmake TRUE)
endif()

# 2- Try LAPACK via PKGCONFIG
if (NOT LAPACK_FOUND)
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
    if (NOT LAPACK_LIBRARIES)
      set(LAPACK_FOUND FALSE)
    else()
      set(LAPACK_FOUND_with_pkgconfig TRUE)
    endif()
  endif()
endif()

# 3- Try LAPACK via environment variables
if (NOT LAPACK_FOUND)
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
    set_target_properties(LAPACK::LAPACK PROPERTIES IMPORTED_LOCATION "${LAPACK_LIBRARY}")
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
    set(LAPACK_FOUND_with_environment TRUE)
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
  # Fortran program calling BLAS/Lapack
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_linalg/test_linalg.F90 "
program test_linalg
  implicit none
  integer, parameter :: n = 2
  real(8) :: A(2,2), b(2), x(2)
  integer :: ipiv(2), info
  A = reshape([1.0d0, 2.0d0, 3.0d0, 4.0d0], [2,2])
  b = [5.0d0, 6.0d0]
  call dcopy(n, b, 1, x, 1)
    call dgesv(n, 1, A, n, ipiv, x, n, info)
  stop (info == 0)  ! 0=OK, 1=échec
end program test_linalg
")
  try_run(LINALG_RUN_RESULT LINALG_COMPILE_RESULT
          ${CMAKE_BINARY_DIR}/tests/test_linalg
          ${CMAKE_BINARY_DIR}/tests/test_linalg/test_linalg.F90
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
