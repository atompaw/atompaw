PROGRAM atompaw
  !*************************************************************
  !  Driver program for new version of atompaw program
  !      arguments can be SAVEAEATOM (save file name),  OR
  !             LOADAEATOM (load file name)
  !      Note:  must have 0 or 2 arguments
  !************************************************************
  USE GlobalMath
  USE aeatom
  USE atomdata
  USE atompaw_report
  USE pseudo
  USE abinitinterface
  USE pwscfinterface

  IMPLICIT NONE
  INTEGER :: i,j,iargc
  CHARACTER(120) :: inputfile,outputfile,token
  LOGICAL :: lotsofoutput=.false.    ! can be changed eventually
  LOGICAL :: saveaeatom=.false.,loadaeatom=.false.
  LOGICAL :: OK
  INTEGER :: ifinput=10,ifen=11


  j=iargc()
  if (j==2) then
       call GetArg(1,token)
       call UpperCase(token)
       if (TRIM(token)=='SAVEAEATOM') then
           saveaeatom=.true.
           call GetArg(2,outputfile)
           write(6,*) 'AEATOM results will be saved to file ',TRIM(outputfile)
       else if (TRIM(token)=='LOADAEATOM') then
           loadaeatom=.true.
           call GetArg(2,inputfile)
           write(6,*) 'AEATOM results will be loaded from file ',TRIM(inputfile)
           Inquire(file=TRIM(inputfile),exist=OK)
           If (.not.OK) then
              write(6,*) 'Load file does not exist -- program will stop'
              stop
           endif
       else
          write(6,*) 'Argument form not recognized', token
          stop
       endif

  else if (j==0) then
       write(6,*) 'input/output not saved/dumped in this run'
  else
       write(6,*) 'Argument form not recognized', j
       stop
  endif

  exploremode=.false.

  OPEN(ifinput,file='dummy',form='formatted')

  CALL Init_GlobalConstants()

  If (loadaeatom) then
     Call Load_AEatom(TRIM(inputfile),&
        Grid,AEOrbit,AEpot,AESCF,FCOrbit,FCpot,FCSCF,FC,ifinput)
     OPEN (ifen, file=TRIM(AEPot%sym), form='formatted')
     CALL Report_Atomres('AE',Grid,AEOrbit,AEPot,AESCF,ifen)
     CALL Report_Atomres('SC',Grid,FCOrbit,FCPot,FCSCF,ifen)
  else
     CALL Init_SCFatom(ifinput)
     CALL SCFatom('AE',lotsofoutput)

     OPEN (ifen, file=TRIM(AEPot%sym), form='formatted')
     CALL Report_Atomres('AE',Grid,AEOrbit,AEPot,AESCF,ifen)

     CALL SCFatom('SC',lotsofoutput,ifinput)
     CALL Report_Atomres('SC',Grid,FCOrbit,FCPot,FCSCF,ifen)
     if (saveaeatom) then
        Call Dump_AEatom(TRIM(outputfile),&
          Grid,AEOrbit,AEpot,AESCF,FCOrbit,FCpot,FCSCF,FC)
     endif
  endif

  CALL INIT_PAW(ifinput,ifen,Grid)
  Call initpseudopot(Grid) 
  CALL setbasis(Grid,FCPot,FCOrbit,ifinput)
  Call setcoretail(Grid,FC%coreden)
  If (TRIM(FCorbit%exctype)=='HF'.or.TRIM(FCorbit%exctype)=='EXXKLI') then
     PAW%tcore=0
  Endif
  If (TRIM(FCorbit%exctype)=='EXXKLI') Call fixtcorewfn(Grid,PAW)
  Call Set_PAW_Options(ifinput,ifen,Grid,FCOrbit,FCPot,OK)
  If (.not.OK) then
          write(6,*) 'Stopping due to options failure'
          stop
  endif        
  Call Report_Pseudobasis(Grid,PAW,ifen)

  Call Set_PAW_MatrixElements(Grid,PAW)
!  If (TRIM(FCOrbit%exctype)/='HF') CALL logderiv(Grid,FCPot,PAW)
  CALL logderiv(Grid,FCPot,PAW)
  CALL ftprod(Grid)

  CALL FindVlocfromVeff(Grid,FCOrbit,PAW)
  CALL Report_Pseudopotential(Grid,PAW)

  CALL SPMatrixElements(Grid,FCPot,FC,PAW)
  !CALL WRITE_ATOMDATA(Grid,FCPot,FCOrbit,FC,PAW)
  CALL Report_pseudo_energies(PAW,6)
  CALL Report_pseudo_energies(PAW,ifen)

    Do
    WRITE(6,*) 'Enter 0 or END to end program'
       WRITE(6,*) 'Enter 1 or SCFPAW to run SCFPAW'
       WRITE(6,*) 'Enter 2 or ABINITOUT to run atompaw2abinit'
       WRITE(6,*) 'Enter 3 or PWSCFOUT to run atompaw2pwscf'
       WRITE(6,*) 'Enter 4 or PWPAWOUT  to run atompaw2pwpaw'
       WRITE(6,*) 'Enter 5 or XMLOUT  to run atompaw2xml'
       WRITE(6,*) 'Enter 10 or EXPLORE to run exploreparms'
                     
       READ(5,'(a)',iostat=i) token
           if (i/=0) exit
           !if (i==0) then
           if (checkline2(token,"0","END")) then
               exit
           !else if (i==1) then    
           else if (checkline2(token,"1","SCFPAW")) then
               CALL SCFPAW(Grid,PAW)
           !else if (i==2) then    
           else if (checkline2(token,"2","ABINITOUT")) then
               CALL atompaw2abinit(AEOrbit,AEPot,PAW,FC,Grid)
           !else if (i==3) then    
           else if (checkline2(token,"3","PWSCFOUT")) then    
               CALL atompaw2pwscf(Grid,AEPot,FC,PAW,AEOrbit,ifinput)
           else if (checkline2(token,"4","PWPAWOUT")) then    
               CALL WRITE_ATOMDATA(Grid,FCPot,FCOrbit,FC,PAW)
           else if (checkline2(token,"5","XMLOUT")) then    
               CALL WRITE_XML(Grid,FCPot,FCOrbit,FC,PAW) !dummy
           !else if (i==10) then
           else if (checkline2(token,"10","EXPLORE")) then    
                   close(ifinput); close(ifen)
                   CALL exploreparms(Grid,AEPot,FC,AEOrbit,PAW)    
           else 
               write(6,*) 'Option not recognized -- pgm terminating'
               stop
           endif
                                 
    Enddo

!  Free Memory
    call DestroyGrid(Grid)
    Call DestroyOrbit(AEOrbit)
    Call DestroyOrbit(FCOrbit)
    Call DestroyFC(FC)
    Call DestroyPot(AEPot)
    Call DestroyPot(FCPot)
    Call DestroyPAW(PAW)
    !  (have_libxc) call libxc_end()

END PROGRAM atompaw
