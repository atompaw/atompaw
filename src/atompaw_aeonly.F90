  !*************************************************************
  !  atompaw_aeonly
  !  program to calculate the self-consistent density functional
  !    atom ground state for atom with atomic number nz and
  !    given electron configuration with optional frozencore
  !    treatment           9/9/2025     NAWH -- previously "graphatom"
  !    
  !   In this version a required argument is an input file
  !     an optional second argument is verbose which may or may not
  !     actually do something.   NAWH 10/31/2024
  !************************************************************

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

PROGRAM atompaw_aeonly

  USE io_tools
  USE GlobalMath
  USE gridmod
  USE aeatom
  USE atomdata
  USE aeonly_report
  USE radialked
  USE splinesolver
  USE xmlinterface

  IMPLICIT NONE

  INTEGER :: i,ios,iargc
  CHARACTER(80) :: verbose,token
  LOGICAL :: lotsofoutput=.false.
  TYPE(mesh_data_type) :: mesh_data
  CHARACTER(132) ::xcname,file_xml_core

  if (iargc()>0) then
       call GetArg(1,verbose)
       call UpperCase(verbose)
       if (TRIM(verbose)=='VERBOSE') lotsofoutput=.true.
  endif

  CALL Init_GlobalConstants()

  !Read input parameters
  CALL input_dataset_read(echofile='dummy',&
&        read_global_data=.true.,read_elec_data=.true.,&
&        read_coreval_data=.false.,read_basis_data=.false.)

  CALL SCFatom_Init()
  CALL SCFatom('AE',lotsofoutput)

  CALL Report_aeonly('AE',Grid,AEOrbit,AEPot,AESCF)



  DO
     WRITE(STD_OUT,*) ' Input 0 or END for plotting results and completing program (once)'
     WRITE(STD_OUT,*) ' Input 1 or FROZENCORE for choosing frozencore in current config. (once)'
     WRITE(STD_OUT,*) ' Input 2 or PRINTOUTCORE for creating core information in xml format '

     READ(STD_IN,'(a)',iostat=ios) token
     WRITE(STD_OUT,'(a)') 'Token is ', token;call flush(STD_OUT)
     IF (ios/=0) EXIT

     CALL eliminate_comment(token)
     WRITE(STD_OUT,*) 'Token is ', token;call flush(STD_OUT)

     IF (checkline2(token,"0","END")) THEN
       WRITE(STD_OUT,'(a)') 'End of input detected; program ending normally'
       EXIT
     ELSEIF (checkline2(token,"1","FROZENCORE")) THEN
         !Read input parameters
         WRITE(STD_OUT,*) 'In case FROZENCORE '; call flush(STD_OUT)
         CALL input_dataset_read(echofile='dummy',&
&        read_global_data=.false.,read_elec_data=.false.,&
&        read_coreval_data=.true.,read_basis_data=.false.)

         CALL SCFatom('SC',lotsofoutput)
         write(std_out,*) ' Finished SC in aeonly '; call flush_unit(std_out)
         CALL Report_aeonly('FC',Grid,FCOrbit,FCPot,FCSCF)
         write(std_out,*) ' Finished SC report in aeonly '; call flush_unit(std_out)
     ELSEIF (checkline2(token,"2","PRINTOUTCORE")) THEN
         WRITE(STD_OUT,*) 'In case PRINTOUTCORE '; call flush(STD_OUT)
         call build_mesh_data(mesh_data,Grid)
         xcname=exctype;if (have_libxc) call libxc_getshortname(exctype,xcname)
         file_xml_core=TRIM(AEpot%sym)//'.'//TRIM(xcname)//'-paw.corewf.xml'
         call xmlprtcore(file_xml_core,Grid,AEOrbit,AEPot,FC,mesh_data)
     ELSE
         WRITE(STD_OUT,'(a)') 'Unrecognized option -- program ending'
         EXIT
     ENDIF
  ENDDO

  if (scalarrelativistic) CALL deallocate_Scalar_Relativistic
  if (diracrelativistic) CALL deallocate_Dirac_Relativistic
  if (needvtau) CALL deallocate_ked
  if (needvtau) CALL deallocatesplinesolver
  CALL input_dataset_free()

END PROGRAM atompaw_aeonly
