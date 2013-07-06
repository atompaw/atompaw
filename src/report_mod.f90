MODULE report_mod
  USE atomdata
  USE gridmod
  USE excor

  IMPLICIT NONE

CONTAINS


  SUBROUTINE Total_Energy_Report(SCF,IF)
    TYPE (SCFinfo) ,INTENT(INOUT) :: SCF
    INTEGER, INTENT(IN) :: IF

    WRITE(IF,*)
    WRITE(IF,*) ' Total energies'
    WRITE(IF,*) '    One-electron contribution:  ',SCF%eone
    WRITE(IF,*) '    Kinetic energy contribution:',SCF%ekin
    WRITE(IF,*) '    Coulomb contribution     :  ',SCF%ecoul
    WRITE(IF,*) '    Electrostatic contribution: ',SCF%estatic
    WRITE(IF,*) '    Exch        contribution :  ',SCF%eexc
      SCF%etot=SCF%ekin+SCF%estatic+SCF%eexc
    WRITE(IF,*) '    Ratio Pot/Kin            :  ',&
&                  (SCF%estatic+SCF%eexc)/SCF%ekin
   If (ColleSalvetti) then
    WRITE(IF,*) '    CS correlation           :  ',SCF%oepcs
    SCF%etot=SCF%etot+SCF%oepcs
   endif

    WRITE(IF,*) '    Total                    :  ',SCF%etot

  END SUBROUTINE Total_Energy_Report

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE One_electron_energy_Report(Orbit,IF)
    TYPE (Orbitinfo) ,INTENT(IN) :: Orbit
    INTEGER, INTENT(IN) :: IF

    INTEGER :: io

    WRITE(IF,*)
    WRITE(IF,*) 'Orbital energies'
    WRITE(IF,"(' n  l     occupancy            energy')")
    IF (frozencorecalculation) THEN
       DO io=1,Orbit%norbit
          IF (.NOT.Orbit%iscore(io)) &
&              WRITE(IF,'(i2,1x,i2,4x,1p,2e15.7)') &
&              Orbit%np(io),Orbit%l(io),Orbit%occ(io),Orbit%eig(io)
       ENDDO
    ELSE
       DO io=1,Orbit%norbit
          WRITE(IF,'(i2,1x,i2,4x,1p,2e15.7)') &
&              Orbit%np(io),Orbit%l(io),Orbit%occ(io),Orbit%eig(io)
       ENDDO
    ENDIF
  END SUBROUTINE One_electron_energy_Report

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Core_Electron_Report(Orbit,FC,IF)
    TYPE (Orbitinfo) ,INTENT(IN) :: Orbit
    TYPE (FCinfo) ,INTENT(IN) :: FC
    INTEGER, INTENT(IN) :: IF

    INTEGER :: io

    IF (.NOT.frozencorecalculation) RETURN

    WRITE(IF,*)
    WRITE(IF,*) ' All-Electron core states (zcore)', FC%zcore
    WRITE(IF,"(' n  l     occupancy            energy')")
    DO io=1,Orbit%norbit
       IF (Orbit%iscore(io)) &
&           WRITE(IF,'(i2,1x,i2,4x,1p,2e15.7)') &
&           Orbit%np(io),Orbit%l(io),Orbit%occ(io),Orbit%eig(io)
    ENDDO
  END SUBROUTINE Core_Electron_Report

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Valence_Electron_Report(Orbit,FC,IF)
    TYPE (Orbitinfo) ,INTENT(IN) :: Orbit
    TYPE (FCinfo) ,INTENT(IN) :: FC
    INTEGER, INTENT(IN) :: IF

    INTEGER :: io

    IF (.NOT.frozencorecalculation) RETURN

    WRITE(IF,*)
    WRITE(IF,*) ' All-Electron valence states (zvale)', FC%zvale
    WRITE(IF,*) ' Below are listed the All-Electron valence states'
    WRITE(IF,"(' n  l     occupancy            energy')")
    DO io=1,Orbit%norbit
       IF (.not.Orbit%iscore(io)) &
&           WRITE(IF,'(i2,1x,i2,4x,1p,2e15.7)') &
&           Orbit%np(io),Orbit%l(io),Orbit%occ(io),Orbit%eig(io)
    ENDDO
  END SUBROUTINE Valence_Electron_Report

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Total_FCEnergy_Report(SCF,IF)
    TYPE (SCFinfo) ,INTENT(INOUT) :: SCF
    INTEGER, INTENT(IN) :: IF

    WRITE(IF,*)
    WRITE(IF,*) ' Valence energies'
    WRITE(IF,*) '    Kinetic energy contribution:',SCF%valekin
    WRITE(IF,*) '    Coulomb contribution     :  ',SCF%valecoul
    WRITE(IF,*) '    Exch        contribution :  ',SCF%valeexc

     SCF%evale=SCF%valekin+SCF%valecoul+SCF%valeexc

    WRITE(IF,*) '    Total valence            :  ',SCF%evale

  END SUBROUTINE Total_FCEnergy_Report

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Report SCF AE
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Report_AESCF(Grid,AEPot,AEOrbit)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: AEPot
    TYPE(OrbitInfo), INTENT(IN) :: AEOrbit
    INTEGER i,j,n

    n=Grid%n
    OPEN (unit=1001,file='AE_pot')
    DO i = 1,n
       WRITE(1001,'(1p,7e15.7)') Grid%r(i),AEPot%rv(i) ,&
&           AEPOt%rvx(i),AEPot%rvh(i),AEPot%rvn(i),AEOrbit%den(i)
    ENDDO
    CLOSE(1001)

    OPEN (unit=1001,file='AE_wfn')
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Grid%r(i),&
&           (AEOrbit%wfn(i,j),j=1,AEOrbit%norbit)
    ENDDO
    CLOSE(1001)

  END SUBROUTINE Report_AESCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Report SCF FC
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Report_FCSCF(Grid,FCPOT,FCOrbit)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: FCPot
    TYPE(OrbitInfo), INTENT(IN) :: FCOrbit

    INTEGER i,j,n

    n=Grid%n
    OPEN (unit=2001,file='FC_Pot')
    DO i = 1,n
       WRITE(1001,'(1p,7e15.7)') Grid%r(i),FCPot%rv(i) ,&
&           FCPot%rvx(i),FCPot%rvh(i),FCPot%rvn(i),FCOrbit%den(i)
    ENDDO
    CLOSE(2001)


    OPEN (unit=2001,file='FC_wfn')
    DO i = 1,n
       WRITE(2001,'(1p,50e15.7)') Grid%r(i),&
&           (FCOrbit%wfn(i,j),j=1,FCOrbit%norbit)
    ENDDO
    CLOSE(2001)

  END SUBROUTINE Report_FCSCF

  SUBROUTINE summary_report(ifen,key,Grid,Orbit,Pot,SCF)
    INTEGER, INTENT(IN) :: ifen
    CHARACTER(2), INTENT(IN) :: key
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(SCFInfo), INTENT(IN) :: SCF


    WRITE(ifen,*)
    WRITE(ifen,*) 'Completed calculations for ',TRIM(POT%sym)
    IF (TRIM(Orbit%exctype)=='EXX') THEN
       WRITE(ifen,*) 'Exact exchange calculation'
    ELSEIF (TRIM(Orbit%exctype)=='EXXKLI') THEN
       WRITE(ifen,*) 'Exact exchange KLI calculation'
    ELSEIF (TRIM(Orbit%exctype)=='EXXOCC') THEN
       WRITE(ifen,*) 'Exact exchange OCC calculation'
    ELSEIF (TRIM(Orbit%exctype)=='HF') THEN
       WRITE(ifen,*) 'Hartree Fock calculation'
    ELSE
       CALL Report_EXC(ifen)
    ENDIF

    CALL reportgrid(Grid,ifen)
    IF (scalarrelativistic) THEN
       IF(.NOT.finitenucleus) THEN
          WRITE(ifen,*) 'Scalar relativistic calculation -- point nucleus'
       ELSE
          WRITE(ifen,*) &
&              'Scalar relativistic calculation -- finite (Gaussian) nucleus'
       ENDIF
    ELSE
       WRITE(ifen,*) 'Non-relativistic calculation'
    ENDIF

    IF (key=='AE'.OR.key=='NC') &
&        WRITE(ifen,*) '  AEatom converged in',SCF%iter,' iterations'
    IF (key=='FC'.OR.key=='SC') &
&        WRITE(ifen,*) '  FCatom converged in',SCF%iter,' iterations'
    WRITE(ifen,'(a,f6.2)') '     for nz = ',Pot%nz
    WRITE(ifen,*) '    delta  = ', SCF%delta
    CALL One_electron_energy_Report(Orbit,ifen)
    WRITE(ifen,*)
    WRITE(ifen,*) ' Total energy'
    WRITE(ifen,*) '    Total                    :  ',SCF%etot
    IF (key=='FC'.OR.key=='SC') THEN
       WRITE(ifen,*) '    Valence                  :  ',SCF%evale
    ENDIF

  END SUBROUTINE summary_report

END MODULE report_mod

