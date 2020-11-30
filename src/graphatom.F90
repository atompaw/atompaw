  !*************************************************************
  !  program to calculate the self-consistent density functional
  !    atom ground state for atom with atomic number nz
  !************************************************************

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

PROGRAM graphatom

  USE io_tools
  USE GlobalMath
  USE aeatom
  USE atomdata
  USE graphatom_report

  IMPLICIT NONE

  INTEGER :: i,ios,iargc
  CHARACTER(80) :: verbose,token
  LOGICAL :: lotsofoutput=.false.

  if (iargc()>0) then
       call GetArg(1,verbose)
       call UpperCase(verbose)
       if (TRIM(verbose)=='VERBOSE') lotsofoutput=.true.
  endif

  !Read input parameters
  CALL input_dataset_read(echofile='dummy',&
&        read_global_data=.true.,read_elec_data=.true.,&
&        read_coreval_data=.false.,read_basis_data=.false.)

  CALL Init_GlobalConstants()
  CALL SCFatom_Init()
  CALL SCFatom('AE',lotsofoutput)

  CALL Report_Graphatom('AE',Grid,AEOrbit,AEPot,AESCF)

  DO
     WRITE(STD_OUT,*) ' Input 0 or END for plotting results and completing program (once)'
     WRITE(STD_OUT,*) ' Input 1  for changing all-electron configuration (many)'
     WRITE(STD_OUT,*) ' Input 2  for choosing frozencore in current config. (once)'
     WRITE(STD_OUT,*) ' Input 3  for changing frozen-core configuration (many)'

     READ(STD_IN,'(a)',iostat=ios) token
     IF (ios/=0) EXIT

     CALL eliminate_comment(token)

     IF (checkline2(token,"0","END")) THEN
       EXIT
     ELSE
       READ(token,*) i
       IF (i == 1) THEN
         CALL SCFatom('NC',lotsofoutput)
         CALL Report_Graphatom('NC',Grid,AEOrbit,AEPot,AESCF)
       ELSE IF (i == 2) THEN
         CALL SCFatom('SC',lotsofoutput)
         write(std_out,*) ' Finished SC in graphatom '; call flush_unit(std_out)
         CALL Report_Graphatom('FC',Grid,FCOrbit,FCPot,FCSCF)
         write(std_out,*) ' Finished SC report in graphatom '; call flush_unit(std_out)
       ELSE IF (i == 3) THEN
         write(std_out,*) 'before FC in graphtom' ; call flush_unit(std_out)
         CALL SCFatom('FC',lotsofoutput)
         CALL Report_Graphatom('FC',Grid,FCOrbit,FCPot,FCSCF)
       ENDIF
     ENDIF
  ENDDO

  if (scalarrelativistic) CALL deallocate_Scalar_Relativistic
  if (diracrelativistic) CALL deallocate_Dirac_Relativistic
  CALL input_dataset_free()

END PROGRAM graphatom
