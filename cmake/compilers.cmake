# ATOMPAW - Additional compiler tests

# --- ISO_C_BINDING check

message(CHECK_START "Checking that Fortran compiler supports ISO C bindings")

# Fortran file
file(WRITE ${CMAKE_BINARY_DIR}/test_iso_c/test_iso_c.f90 "
program test_iso_c
 use, intrinsic :: iso_c_binding
 implicit none
  integer(c_int)   :: i32 = 42_c_int
  type(c_ptr) :: cptr
  print *, 'i32 =', i32
end program test_iso_c
")

# Compilation
try_compile(ISO_C_OK
  ${CMAKE_BINARY_DIR}/test_iso_c
  SOURCES ${CMAKE_BINARY_DIR}/test_iso_c/test_iso_c.f90
  OUTPUT_VARIABLE OUTPUT
)

# Result
if (ISO_C_OK)
  set(HAVE_FC_ISO_C_BINDING 1)
  message(CHECK_PASS "OK")
else()
  message(FATAL_ERROR "Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} doesnt support ISO C Bindings!")
endif()


# --- ABI check - Fortran calling C

message(CHECK_START "Checking that C and Fortran compilers can link")

# C subroutine
file(WRITE ${CMAKE_BINARY_DIR}/test_abi/addone.c "
int addone(int n) {return n + 1;}
")

# Fortran main program
file(WRITE ${CMAKE_BINARY_DIR}/test_abi/test_abi.f90 "
program test_abi
 use iso_c_binding
 implicit none
  interface
   integer function addone(n) bind(C)
    import :: c_int
    integer(c_int), intent(in), value :: n
   end function
  end interface
  if (addone(2) == 3) stop 0
  if (addone(2) /= 3) stop 1
end program test_abi
")

# Compilation
try_compile(ABI_OK
  ${CMAKE_BINARY_DIR}/test_abi
  SOURCES ${CMAKE_BINARY_DIR}/test_abi/addone.c
          ${CMAKE_BINARY_DIR}/test_abi/test_abi.f90
  OUTPUT_VARIABLE OUTPUT
)

# Result
if (ABI_OK)
  message(CHECK_PASS "OK")
else()
  message(FATAL_ERROR "C compiler {CMAKE_C_COMPILER_ID} {CMAKE_C_COMPILER_VERSION} and Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} are ABI-incompatible!")
endif()
