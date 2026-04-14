# ATOMPAW - Additional compiler tests

# --- ISO_C_BINDING check

message(CHECK_START "Checking if Fortran compiler supports ISO C bindings")

# Fortran program using ISO_C_BINDINGS
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_iso_c/test_iso_c.f90 "
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
  ${CMAKE_BINARY_DIR}/tests/test_iso_c
  SOURCES ${CMAKE_BINARY_DIR}/tests/test_iso_c/test_iso_c.f90
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

message(CHECK_START "Checking if C and Fortran compilers can link")

# C subroutine
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_abi/addone.c "
int addone(int n) {return n + 1;}
")

# Fortran program calling C
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_abi/test_abi.f90 "
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
  ${CMAKE_BINARY_DIR}/tests/test_abi
  SOURCES ${CMAKE_BINARY_DIR}/tests/test_abi/addone.c
          ${CMAKE_BINARY_DIR}/tests/test_abi/test_abi.f90
  OUTPUT_VARIABLE OUTPUT
)

# Result
if (ABI_OK)
  message(CHECK_PASS "OK")
else()
  message(FATAL_ERROR "C compiler {CMAKE_C_COMPILER_ID} {CMAKE_C_COMPILER_VERSION} and Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} are ABI-incompatible!")
endif()


# --- FLUSH check

message(CHECK_START "Checking if Fortran compiler supports FLUSH or FLUSH_")

enable_language(Fortran)

# Fortran program using FLUSH
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_flush/test_flush.f90 "
program test_flush
  call flush(1)
end program test_flush
")

# Fortran program using FLUSH_
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_flush/test_flush_.f90 "
program test_flush_
  call flush_(1)
end program test_flush_
")

# Compilation
try_compile(FLUSH_OK
  ${CMAKE_BINARY_DIR}/tests/test_flush
  SOURCES ${CMAKE_BINARY_DIR}/tests/test_flush/test_flush.f90
  OUTPUT_VARIABLE OUTPUT
)
try_compile(FLUSH__OK
  ${CMAKE_BINARY_DIR}/tests/test_flush
  SOURCES ${CMAKE_BINARY_DIR}/tests/test_flush/test_flush_.f90
  OUTPUT_VARIABLE OUTPUT
)

# Result
if (FLUSH_OK)
  set(HAVE_FC_FLUSH 1)
endif()
if (FLUSH__OK)
  set(HAVE_FC_FLUSH_ 1)
endif()

if (FLUSH_OK OR FLUSH__OK)
  message(CHECK_PASS "OK")
else()
  message(CHECK_PASS "not OK")
endif()


# --- ISATTY check

message(CHECK_START "Checking if Fortran compiler supports ISATTY")

enable_language(Fortran)

# Fortran program using ISATTY
file(WRITE ${CMAKE_BINARY_DIR}/tests/test_isatty/test_isatty.f90 "
program test_isatty
  logical :: success
  success=isatty(1)
end program test_isatty
")

# Compilation
try_compile(ISATTY_OK
  ${CMAKE_BINARY_DIR}/tests/test_isatty
  SOURCES ${CMAKE_BINARY_DIR}/tests/test_isatty/test_isatty.f90
  OUTPUT_VARIABLE OUTPUT
)

# Result
if (ISATTY_OK)
  set(HAVE_FC_ISATTY 1)
  message(CHECK_PASS "OK")
else()
  message(CHECK_PASS "not OK")
endif()
