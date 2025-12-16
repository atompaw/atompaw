!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the definition of standard I/O streams
!      and some I/O related routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE io_tools

!  Changed 5/23/2022 by NAWH following suggestion by J. Zwanziger

#if defined (__INTEL_COMPILER)
 USE IFPORT
#endif


 IMPLICIT NONE

 PRIVATE

!Public constants - I/O streams (can be modified)
 INTEGER,SAVE,PUBLIC :: STD_IN  = 5
 INTEGER,SAVE,PUBLIC :: STD_OUT = 6
 INTEGER,SAVE,PUBLIC :: STD_ERR = 0

!Public functions
 PUBLIC :: flush_unit     ! Wrapper for the standard flush_unit routine
 PUBLIC :: unit_isatty    ! Wrapper for the Fortran extension isatty function
 PUBLIC :: copy_file      ! Copy a text file into another
 PUBLIC :: compare_files  ! Compare text files

 CONTAINS

!**********************************************************************
!  flush_unit - wrapper for the standard flush_unit routine
!  Available only if the compiler implements this intrinsic procedure.
!**********************************************************************
 SUBROUTINE flush_unit(unit)
   INTEGER,INTENT(in) :: unit
   LOGICAL :: isopen
   IF (unit==-1) RETURN
   INQUIRE(unit=unit,opened=isopen)
#if defined HAVE_FC_FLUSH
   IF (isopen) THEN
     CALL flush(unit)
   ENDIF
#elif defined HAVE_FC_FLUSH_
   IF (isopen) THEN
     CALL flush_(unit)
   ENDIF
#endif
 END SUBROUTINE flush_unit

!**********************************************************************
!  unit_isatty - wrapper for the Fortran extension isatty function
!  Available only if the compiler implements this intrinsic procedure.
!**********************************************************************
 LOGICAL FUNCTION unit_isatty(unit)
   INTEGER,INTENT(in) :: unit
   LOGICAL :: isatty
#if defined HAVE_FC_ISATTY
     unit_isatty=isatty(unit)
#else
     unit_isatty=.true.
#endif    
 END FUNCTION unit_isatty

!**********************************************************************
!  copy_file - copy a text file into another
!  infile:  input file name
!  outfile: output file name
!**********************************************************************
 SUBROUTINE copy_file(infile,outfile)
   CHARACTER(*),INTENT(in) :: infile,outfile
   INTEGER,PARAMETER :: input_unit=77,output_unit=78
   INTEGER :: ios
   CHARACTER(256) :: inputline
   OPEN(input_unit,FILE=trim(infile),FORM='formatted')
   OPEN(output_unit,FILE=trim(outfile),FORM='formatted')
   DO
     READ(input_unit,fmt='(a)',ERR=79,END=79,IOSTAT=ios) inputline
     IF (ios/=0) EXIT
     WRITE(output_unit,fmt='(a)') trim(inputline)
     GOTO 80
79    EXIT
80    CONTINUE
   ENDDO
   CLOSE(output_unit)
   CLOSE(input_unit)
 END SUBROUTINE copy_file

!**********************************************************************
!  compare_files - compare text files
!  file1, file2:  file names
!  line_count: number of lines to compare
!**********************************************************************
 LOGICAL FUNCTION compare_files(file1,file2,line_count)
   CHARACTER(*),INTENT(in) :: file1,file2
   INTEGER,INTENT(in),OPTIONAL :: line_count
   INTEGER,PARAMETER :: unit1=77,unit2=78
   INTEGER :: ios,nline
   CHARACTER(256) :: line1,line2
   compare_files=.true. ; nline=0
   OPEN(unit1,FILE=trim(file1),FORM='formatted',STATUS='old')
   OPEN(unit2,FILE=trim(file2),FORM='formatted',STATUS='replace')
   DO
     READ(unit1,FMT='(a)',ERR=89,END=89,IOSTAT=ios) line1
     READ(unit2,FMT='(a)',ERR=89,END=89,IOSTAT=ios) line2
     IF (ios/=0) EXIT
     IF (trim(line1)/=trim(line2)) THEN
       WRITE(STD_OUT,*) 'Inconsistency in files'
       WRITE(STD_OUT,'(a)') trim(file1)
       WRITE(STD_OUT,'(a)') trim(file2)
       compare_files=.false.
     ENDIF
     nline=nline+1
     IF (present(line_count)) THEN
       IF (nline>line_count) EXIT
     ENDIF
     GOTO 90
89    EXIT
90    CONTINUE
   ENDDO
   CLOSE(unit1)
   CLOSE(unit2)
 END FUNCTION compare_files

END MODULE io_tools
