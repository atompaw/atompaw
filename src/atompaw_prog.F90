!*************************************************************
!  Driver program for atompaw program
!
!    Energies always in Rydberg (except when preparing output for
!          datasets to be used in programs that need Hartree units)
!    Lengths always in Bohr
!
!    The following is not currently available:
!      arguments can be SAVEAEATOM (save file name),  OR
!                       LOADAEATOM (load file name)
!      Note:  must have 0 or 2 arguments
!************************************************************

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

PROGRAM atompaw

  USE io_tools
  USE GlobalMath
  USE aeatom
  USE atomdata
  USE atompaw_report
  USE pseudo
  USE abinitinterface
  USE pwscfinterface
  USE xmlinterface
  USE libxc_mod
  USE pkginfo
  USE radialked
  USE splinesolver

  IMPLICIT NONE
  INTEGER :: i,j,iargc
  CHARACTER(120) :: inputfile,outputfile,token
  LOGICAL :: lotsofoutput=.false.    ! can be changed eventually
  LOGICAL :: saveaeatom=.false.,loadaeatom=.false.
  LOGICAL :: OK,scfpaw_done
  INTEGER :: ifen=11

! First write out details on AtomPAW
  WRITE(STD_OUT,*) atp_package,' v',atp_version
  WRITE(STD_OUT,*) 'Compiled for ',atp_target
  WRITE(STD_OUT,*)

  j=iargc()
  if (j>=1) then
    call GetArg(1,token)
    call UpperCase(token)
    if (TRIM(token)=='--VERSION') then
      stop
    else if (j==2.and.TRIM(token)=='SAVEAEATOM') then
      saveaeatom=.true.
      call GetArg(2,outputfile)
      write(std_out,*) 'AEATOM results will be saved to file ',TRIM(outputfile)
    else if (j==2.and.TRIM(token)=='LOADAEATOM') then
      loadaeatom=.true.
      call GetArg(2,inputfile)
      write(std_out,*) 'AEATOM results will be loaded from file ',TRIM(inputfile)
      Inquire(file=TRIM(inputfile),exist=OK)
      If (.not.OK) stop 'Load file does not exist -- program will stop'
    else
      write(std_out,*) 'Argument form not recognized: ', token
      stop
    endif
  else
    write(std_out,*) 'Input/output not saved/dumped in this run'
  endif
  write(std_out,*)

! Read input file
  if (.not.loadaeatom) then
    call input_dataset_read(echofile='dummy')
  else
    call input_dataset_read(echofile='tmp')
    OK=compare_files('tmp','dummy',line_count=3+input_dataset%norbit_mod+input_dataset%norbit)
    if (.not.OK) stop 'Program terminated!'
    CALL copy_file('tmp','dummy')
  endif

  exploremode=.false.

  CALL Init_GlobalConstants()

  If (loadaeatom) then
     CALL Load_AEatom(TRIM(inputfile),&
          Grid,AEOrbit,AEpot,AESCF,FCOrbit,FCpot,FCSCF,FC)
     OPEN (ifen, file=TRIM(AEPot%sym), form='formatted')
     CALL Report_Atomres('AE',Grid,AEOrbit,AEPot,AESCF,ifen)
     CALL Report_Atomres('SC',Grid,FCOrbit,FCPot,FCSCF,ifen)
  else
     CALL SCFatom_Init()
     CALL SCFatom('AE',lotsofoutput)
     OPEN (ifen, file=TRIM(AEPot%sym), form='formatted')
     CALL Report_Atomres('AE',Grid,AEOrbit,AEPot,AESCF,ifen)
     CALL SCFatom('SC',lotsofoutput)
     CALL Report_Atomres('SC',Grid,FCOrbit,FCPot,FCSCF,ifen)
     if (saveaeatom) then
        Call Dump_AEatom(TRIM(outputfile),&
             Grid,AEOrbit,AEpot,AESCF,FCOrbit,FCpot,FCSCF,FC)
     endif
  endif

  if(.not.diracrelativistic) then
    CALL SetPAWOptions1(ifen,Grid)
    CALL InitPAW(PAW,Grid,FCOrbit)
    CALL setbasis(Grid,FCPot,FCOrbit)
    Call setcoretail(Grid,FC%coreden)
    Call setttau(Grid,FC%coretau)
    If (TRIM(FCorbit%exctype)=='HF'.or.TRIM(FCorbit%exctype)=='EXXKLI') PAW%tcore=0
    If (TRIM(FCorbit%exctype)=='EXXKLI') Call fixtcorewfn(Grid,PAW)
    Call SetPAWOptions2(ifen,Grid,FCOrbit,FCPot,OK)
    If (.not.OK) stop 'Stopping due to options failure'
    Call Report_Pseudobasis(Grid,PAW,ifen)

    Call Set_PAW_MatrixElements(Grid,PAW,ifen)
    If (TRIM(FCorbit%exctype)/='HF') then
      CALL logderiv(Grid,FCPot,PAW)
      CALL ftprod(Grid)
    endif

    CALL FindVlocfromVeff(Grid,FCOrbit,PAW)
    CALL Report_Pseudopotential(Grid,PAW)

    CALL SPMatrixElements(Grid,FCPot,FC,PAW)
    CALL Report_pseudo_energies(PAW,6)
    CALL Report_pseudo_energies(PAW,ifen)
  else
    write(STD_OUT,*) 'PAW pseudo routines need more work for dirac case!'
  endif
  CLOSE(ifen)

  scfpaw_done=.false.
  Do
   WRITE(STD_OUT,*);WRITE(STD_OUT,*)
   WRITE(STD_OUT,*) 'Enter 0 or END to end program'
   WRITE(STD_OUT,*) 'Enter 1 or SCFPAW to run SCFPAW'
   WRITE(STD_OUT,*) 'Enter 2 or ABINITOUT to run atompaw2abinit'
   WRITE(STD_OUT,*) 'Enter 3 or UPFOUT or PWSCFOUT to run atompaw2pwscf'
   WRITE(STD_OUT,*) 'Enter 4 or XMLOUT  to run atompaw2xml'
   WRITE(STD_OUT,*) 'Enter 5 or PWPAWOUT  to run atompaw2pwpaw'
   WRITE(STD_OUT,*) 'Enter 6 or SOCORROOUT  to run atompaw2socorro'
   WRITE(STD_OUT,*) 'Enter 10 or EXPLORE to run exploreparms'

   READ(STD_IN,'(a)',iostat=i) token
   if (i/=0) exit

   CALL eliminate_comment(token)

   if (checkline2(token,"0","END")) then
     write(std_out,*) 'End atompaw.' ; exit
   else if (checkline2(token,"1","SCFPAW")) then
     CALL SCFPAW(Grid,PAW)
     scfpaw_done=.true.
   else
     if (scfpaw_done) then
       WRITE(STD_OUT,'(2a)') trim(token),&
 &       ' cannot be executed after SCFPAW -- pgm terminating!'
       EXIT
     else
       if (checkline2(token,"2","ABINITOUT")) then
         CALL Atompaw2Abinit(FCOrbit,FCPot,FCSCF,PAW,FC,Grid)
       else if (checkline2(token,"3","PWSCFOUT").or.checkline2(token,"3","UPFOUT")) then
         CALL Atompaw2Pwscf(Grid,FCPot,FC,PAW,FCOrbit)
       else if (checkline2(token,"4","XMLOUT")) then
         CALL Atompaw2XML(FCOrbit,FCPot,FCSCF,PAW,FC,Grid)
       else if (checkline2(token,"5","PWPAWOUT")) then
         CALL WRITE_ATOMDATA(Grid,FCPot,FCOrbit,FC,PAW)
       else if (checkline2(token,"6","SOCORROOUT")) then
         CALL WRITE_SOCORROATOMDATA(Grid,FCPot,FCOrbit,FC,PAW)
       else if (checkline2(token,"10","EXPLORE")) then
         CALL exploreparms(Grid,FCPot,FC,FCOrbit,PAW)
       else
         WRITE(STD_OUT,'(3a)') 'Option not recognized ',&
 &             trim(token),' -- pgm terminating!'
         CALL flush_unit(std_out)
         EXIT
       endif
     endif
   endif

 Enddo

! Free Memory
  Call DestroyGrid(Grid)
  Call DestroyOrbit(AEOrbit)
  Call DestroyOrbit(FCOrbit)
  Call DestroyPot(AEPot)
  Call DestroyPot(FCPot)
  Call DestroyPAW(PAW)
  Call DestroyFC(FC)
  if (scalarrelativistic) CALL deallocate_Scalar_Relativistic
  if (have_libxc) call libxc_end_func()
  if(needvtau) call Deallocate_ked
  if (usespline) call deallocatesplinesolver
  Call input_dataset_free()

END PROGRAM atompaw
