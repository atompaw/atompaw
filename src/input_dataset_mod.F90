!! NAME
!!  input_dataset_mod
!!
!! FUNCTION
!!  This module contains objects and routines used to parse the atompaw
!!  input file. The complete input dataset is stored in the `input_dataset
!!  Fortran data-structure.

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE input_dataset_mod

  IMPLICIT NONE
  
  PRIVATE

!Public functions
 public :: input_dataset_read  ! Initialize input_dataset by reading it from file
 public :: input_dataset_free  ! Destroy input_dataset data-structure

!Public structured datatype
 type,public :: input_dataset_t
   CHARACTER(2) :: atomic_symbol  ! Atomic symbol
   INTEGER      :: atomic_charge  ! Atomic charge     
 end type input_dataset_t

!Public variable containing the complete input dataset
 type(input_dataset_t),public,save :: input_dataset


CONTAINS

!!=================================================================
!! NAME
!!  input_dataset_read
!!
!! FUNCTION
!!  Initialize an input_dataset data-structure by reading it from
!!  a file. If file is omitted, then read from standard input.
!!
!! INPUTS
!!  [inputfile]= name of input file to be read (optional) 
!!
!! OUTPUT
!!  [input_dt]= data-structure containing the complete input file.
!!              If omitted, then the global public `input_dataset`
!!              is used.
!!
!!=================================================================
 SUBROUTINE input_dataset_read(input_dt,inputfile) ! Optional

 character*(*),intent(in),optional :: inputfile
 type(input_dataset_t),intent(out) :: input_dt
 
!---- Local variables

!------------------------------------------------------------------

 END SUBROUTINE input_dataset_read


!!=================================================================
!! NAME
!!  input_dataset_free
!!
!! FUNCTION
!!  Free and destroy a input_dataset data-structure.
!!
!! SIDE EFFECT
!!  input_dt= datstructure to br destroyed
!!
!!=================================================================

 SUBROUTINE input_dataset_free(input_dt)

 type(input_dataset_t),intent(inout) :: input_dt
 
!---- Local variables

!------------------------------------------------------------------

  input_dt%atomic_symbol = ''
  input_dt%atomic_charge = 0

 END SUBROUTINE input_dataset_free

!!=================================================================

END MODULE input_dataset_mod
