# ATOMPAW - Linear Algebra library detection
# ---
# Define target  : BLAS::BLAS
# Define variable: BLAS_LIBRARIES
# Define target  : LAPACK::LAPACK
# Define variable: LAPACK_LIBRARIES
# ---
# Define variable: BLAS_LAPACK_OK if BLAS/LAPACK is usable (or not)
# ---
# Note: may use PkgConfig_FOUND variable
# ---
# Adapt this to your needs
set(CHECK_BLAS   TRUE)  # If TRUE, BLAS is checked and BLAS::BLAS is output
set(CHECK_LAPACK TRUE)  # If TRUE, LAPACK is checked and LAPACK::LAPACK is output
# ---

#Customize messages
if (CHECK_BLAS AND CHECK_LAPACK)
  set(_msg "BLAS+LAPACK")
elseif (CHECK_BLAS)
  set(_msg "BLAS")
elseif (CHECK_LAPACK)
  set(_msg "LAPACK")
else()
  set(_msg "NONE")
endif()

if (CHECK_BLAS OR CHECK_LAPACK)
  message(CHECK_START ">>> Detecting ${_msg}")
endif()

# Store LD_LIBRARY_PATH as a list
string(REPLACE ":" ";" _ld_paths "$ENV{LD_LIBRARY_PATH}")

# ----- BLAS -----

if (CHECK_BLAS)

  # 1- Try BLAS via CMake find_package
  find_package(BLAS QUIET)
  if (BLAS_FOUND)
	if (TARGET BLAS::BLAS)
	  message(STATUS "BLAS found via find_package (config mode)")
	  set(BLAS_FOUND_with_cmake TRUE)
	elseif (BLAS_INCLUDE_DIRS OR BLAS_LIBRARIES)
	  message(STATUS "BLAS found via find_package (module mode)")
	  set(BLAS_FOUND_with_cmake TRUE)
	  if (NOT TARGET BLAS::BLAS)
		add_library(BLAS::BLAS UNKNOWN IMPORTED)
		set_target_properties(BLAS::BLAS PROPERTIES
		  IMPORTED_LOCATION             "${BLAS_LIBRARIES}"
		  INTERFACE_LINK_LIBRARIES      "${BLAS_LIBRARIES}"
		  INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}")
	  endif()
	endif()
  endif()

  # 2- Try BLAS via PKGCONFIG
  if (NOT BLAS_FOUND)
	if (PkgConfig_FOUND)
      foreach(_blas_pkg blas lapack openblas atlas)
        if (NOT BLAS_FOUND)
          pkg_check_modules(BLAS QUIET IMPORTED_TARGET ${_blas_pkg})
        endif()
      endforeach()
	else()
	  set(BLAS_FOUND FALSE)
    endif()
	if (BLAS_FOUND)
	  message(STATUS "BLAS found via pkg-config")
	  if (NOT TARGET BLAS::BLAS)
		add_library(BLAS::BLAS INTERFACE IMPORTED GLOBAL)
		if (TARGET PkgConfig::BLAS)
		  target_link_libraries(BLAS::BLAS INTERFACE PkgConfig::BLAS)
		endif()
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
    foreach(_var BLAS_ROOT BLAS_PREFIX BLAS_DIR)
      if (DEFINED ENV{${_var}} AND NOT BLAS_ROOT)
        set(BLAS_ROOT "$ENV{${_var}}" CACHE PATH "BLAS root")
      endif()
    endforeach()
    foreach(_var LAPACK_ROOT LAPACK_PREFIX LAPACK_DIR)
      if (DEFINED ENV{${_var}} AND NOT _lapack_root)
        set(_lapack_root "$ENV{${_var}}" CACHE PATH "BLAS root")
      endif()
    endforeach()
	find_library(BLAS_LIBRARY
				 NAMES blas openblas atlas
				 HINTS
				   "${BLAS_ROOT}/lib"
				   "${BLAS_ROOT}/lib64"
				   "${_lapack_root}/lib"
				   "${_lapack_root}/lib64"
				   "${_ld_paths}"
				 NO_DEFAULT_PATH
				 NO_CACHE)
	if (BLAS_LIBRARY AND NOT TARGET BLAS::BLAS)
	  message(STATUS "BLAS found in environment: ${BLAS_LIBRARY}")
	  add_library(BLAS::BLAS UNKNOWN IMPORTED)
	  set_target_properties(BLAS::BLAS PROPERTIES IMPORTED_LOCATION "${BLAS_LIBRARY}")
	  set(BLAS_LIBRARIES ${BLAS_LIBRARY})
	  set(BLAS_FOUND_with_environment TRUE)
	endif()
  endif()

endif()  # CHECK_BLAS


# ----- LAPACK -----

if (CHECK_LAPACK)

  # 1- Try LAPACK via CMake find_package
  find_package(LAPACK QUIET)
  if (LAPACK_FOUND)
	if (TARGET LAPACK::LAPACK)
	  message(STATUS "LAPACK found via find_package (config mode)")
	  set(LAPACK_FOUND_with_cmake TRUE)
	elseif (LAPACK_INCLUDE_DIRS OR LAPACK_LIBRARIES)
	  message(STATUS "LAPACK found via find_package (module mode)")
	  set(LAPACK_FOUND_with_cmake TRUE)
	  if (NOT TARGET LAPACK::LAPACK)
		add_library(LAPACK::LAPACK UNKNOWN IMPORTED)
		set_target_properties(LAPACK::LAPACK PROPERTIES
		  IMPORTED_LOCATION             "${LAPACK_LIBRARIES}"
		  INTERFACE_LINK_LIBRARIES      "${LAPACK_LIBRARIES}"
		  INTERFACE_INCLUDE_DIRECTORIES "${LAPACK_INCLUDE_DIRS}")
	  endif()
	endif()
  endif()

  # 2- Try LAPACK via PKGCONFIG
  if (NOT LAPACK_FOUND)
	if (PkgConfig_FOUND)
      foreach(_lapack_pkg lapack openblas atlas lapack_atlas)
        if (NOT LAPACK_FOUND)
          pkg_check_modules(LAPACK QUIET IMPORTED_TARGET ${_lapack_pkg})
        endif()
      endforeach()
	else()
	  set(LAPACK_FOUND FALSE)
	endif()
	if (LAPACK_FOUND)
	  message(STATUS "LAPACK found via pkg-config")
	  if (NOT TARGET LAPACK::LAPACK)
		add_library(LAPACK::LAPACK INTERFACE IMPORTED GLOBAL)
		if (TARGET PkgConfig::LAPACK)
		  target_link_libraries(LAPACK::LAPACK INTERFACE PkgConfig::LAPACK)
		endif()
	  endif()
	  get_target_property(LAPACK_LIBRARIES LAPACK::LAPACK INTERFACE_LINK_LIBRARIES)
	  if (NOT LAPACK_LIBRARIES)
		set(LAPACK_FOUND FALSE)
	  else()
		set(LAPACK_FOUND_with_pkgconfig TRUE)
	  endif()
	endif()
  endif()

  # 3- Try LAPACK via environment variables
  if (NOT LAPACK_FOUND)
    foreach(_var LAPACK_ROOT LAPACK_PREFIX LAPACK_DIR)
      if (DEFINED ENV{${_var}} AND NOT LAPACK_ROOT)
        set(LAPACK_ROOT "$ENV{${_var}}" CACHE PATH "LAPACK root")
      endif()
    endforeach()


	find_library(_blas_lib
				 NAMES blas openblas atlas
				 HINTS
				   "${LAPACK_ROOT}/lib"
				   "${LAPACK_ROOT}/lib64"
				   "${_ld_paths}"
				 NO_DEFAULT_PATH
				 NO_CACHE)
	find_library(_lapack_lib
				 NAMES lapack openblas atlas
				 HINTS
				   "${LAPACK_ROOT}/lib"
				   "${LAPACK_ROOT}/lib64"
				   "${_ld_paths}"
				 NO_DEFAULT_PATH
				 NO_CACHE)
    set(LAPACK_LIBRARY ${_blas_lib} ${_lapack_lib})
	if (LAPACK_LIBRARY AND NOT TARGET LAPACK::LAPACK)
	  message(STATUS "LAPACK found in environment: ${LAPACK_LIBRARY}")
	  add_library(LAPACK::LAPACK UNKNOWN IMPORTED)
	  set_target_properties(LAPACK::LAPACK PROPERTIES IMPORTED_LOCATION "${LAPACK_LIBRARY}")
	  set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
	  set(LAPACK_FOUND_with_environment TRUE)
	endif()
  endif()

endif()  # CHECK_LAPACK


# ----- TESTS -----

# Did we find BLAS/LAPACK?
if (BLAS_FOUND_with_cmake OR BLAS_FOUND_with_pkgconfig OR
    BLAS_FOUND_with_environment)
  set(_blas_found TRUE)
elseif (CHECK_BLAS)
  message(STATUS "BLAS not found!")
endif()
if (LAPACK_FOUND_with_cmake OR LAPACK_FOUND_with_pkgconfig OR
    LAPACK_FOUND_with_environment)
  set(_lapack_found TRUE)
elseif (CHECK_LAPACK)
  message(STATUS "LAPACK not found!")
endif()
if (NOT CHECK_BLAS OR _blas_found)
  if (NOT CHECK_LAPACK OR _lapack_found)
    set(BLAS_LAPACK_FOUND_OK TRUE)
  else()
    set(BLAS_LAPACK_FOUND_OK FALSE)
  endif()
else()
  set(BLAS_LAPACK_FOUND_OK FALSE)
endif()

# If BLAS/LAPACK found, test if it works
if (CHECK_BLAS OR CHECK_LAPACK)
  if (BLAS_LAPACK_FOUND_OK)
    enable_language(Fortran)
    if (CHECK_LAPACK)
      set(_comment "")
    else()
      set(_comment "!")
    endif()
    # Fortran program calling BLAS/Lapack
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_linalg/test_linalg.F90 "
program test_linalg
  implicit none
  integer, parameter :: n = 2
  real(8) :: A(2,2), b(2), x(2)
  integer :: ipiv(2), info=0
  A = reshape([1.0d0, 2.0d0, 3.0d0, 4.0d0], [2,2])
  b = [5.0d0, 6.0d0]
  call dcopy(n, b, 1, x, 1)
  ${_comment}call dgesv(n, 1, A, n, ipiv, x, n, info)
  if (info == 0) stop 0
  if (info /= 0) stop 1
end program test_linalg
")
    try_run(LINALG_RUN_RESULT LINALG_COMPILE_RESULT
            ${CMAKE_BINARY_DIR}/tests/test_linalg
            ${CMAKE_BINARY_DIR}/tests/test_linalg/test_linalg.F90
            LINK_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
            OUTPUT_VARIABLE TRY_OUTPUT)
    if (LINALG_COMPILE_RESULT AND NOT LINALG_RUN_RESULT)
      set(BLAS_LAPACK_OK TRUE)
    else()
      set(BLAS_LAPACK_OK FALSE)
    endif()
  else()
    set(BLAS_LAPACK_OK FALSE)
    message(FATAL_ERROR "${_msg} detection failed!")
  endif()

  if(BLAS_LAPACK_OK)
    message(STATUS "${_msg} Fortran OK")
  else()
    message(STATUS "${_msg} compilation+execution test failed!")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      message(FATAL_ERROR "${TRY_OUTPUT}")
    else()
      message(FATAL_ERROR "${_msg} required but not available!")
    endif()
  endif()
endif()

if (CHECK_BLAS OR CHECK_LAPACK)
  if (BLAS_LAPACK_OK)
    message(CHECK_PASS "done")
  else()
    message(CHECK_FAIL "not found")
  endif()
endif()