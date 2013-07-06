PROGRAM graphatom
  !*************************************************************
  !  program to calculate the self-consistent density functional
  !    atom ground state for atom with atomic number nz
  !************************************************************
  USE GlobalMath
  USE aeatom
  USE atomdata
  USE graphatom_report

  IMPLICIT NONE

  INTEGER :: i,iargc,ifinput
  CHARACTER(80) :: verbose
  LOGICAL :: lotsofoutput=.false.

  if (iargc()>0) then
       call GetArg(1,verbose)
       call UpperCase(verbose)
       if (TRIM(verbose)=='VERBOSE') lotsofoutput=.true.
  endif

  OPEN(ifinput,file='dummy',form='formatted')

  CALL Init_GlobalConstants()
  CALL SCFatom_Init(ifinput)
  CALL SCFatom('AE',lotsofoutput)

  CALL Report_Graphatom('AE',Grid,AEOrbit,AEPot,AESCF)

  CLOSE(ifinput)
  
  DO
     WRITE(6,*) ' Input 0  for plotting results and completing program (once)'
     WRITE(6,*) ' Input 1  for changing all-electron configuration (many)'
     WRITE(6,*) ' Input 2  for choosing frozencore in current config. (once)'
     WRITE(6,*) ' Input 3  for changing frozen-core configuration (many)'


     READ(5,*) i

     IF (i == 0) EXIT
     IF (i == 1) THEN
        CALL SCFatom('NC',lotsofoutput)
        CALL Report_Graphatom('NC',Grid,AEOrbit,AEPot,AESCF)
     ELSE IF (i == 2) THEN
        CALL SCFatom('SC',lotsofoutput)
         write(6,*) ' Finished SC in graphatom '; call flush(6)
        CALL Report_Graphatom('FC',Grid,FCOrbit,FCPot,FCSCF)
         write(6,*) ' Finished SC report in graphatom '; call flush(6)
     ELSE IF (i == 3) THEN
         write(6,*) 'before FC in graphtom' ; call flush(6)
        CALL SCFatom('FC',lotsofoutput)
        CALL Report_Graphatom('FC',Grid,FCOrbit,FCPot,FCSCF)
     ENDIF
  ENDDO

END PROGRAM graphatom
