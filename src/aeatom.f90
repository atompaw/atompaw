MODULE AEatom
  USE atomdata
  USE excor
  USE exx_mod
  USE general_mod
  USE globalmath
  USE gridmod
  USE hf_mod
  USE ldagga_mod

  IMPLICIT NONE
  REAL(8), PARAMETER, PRIVATE ::small0=1.d-13,rimix=0.5d0,worst=1.d-8
  REAL(8), PARAMETER, PRIVATE ::linrange=50.d0,linh=0.0025d0
  REAL(8), PARAMETER, PRIVATE ::logrange=100.d0,logh=0.020d0,lor00=1.d-5
  INTEGER, PARAMETER, PRIVATE :: niter=1000,mxgrid=20000,mxloop=99
  REAL(8), PARAMETER, PRIVATE :: conv1=4.d13,conv2=3.d13,conv3=2.d13,conv4=1.d13
  INTEGER, PRIVATE :: nps,npp,npd,npf,npg,norbit,nz,n
  REAL(8), PRIVATE :: h, electrons, hadjusted

  TYPE (PotentialInfo) ,TARGET :: AEPot
  TYPE (OrbitInfo) ,TARGET :: AEOrbit
  TYPE (SCFInfo) ,TARGET :: AESCF

  TYPE (PotentialInfo) ,TARGET :: FCPot
  TYPE (OrbitInfo) ,TARGET :: FCOrbit
  TYPE (SCFInfo) ,TARGET :: FCSCF
  TYPE(FCInfo),TARGET :: FC

  TYPE (Gridinfo) ,TARGET :: Grid

  TYPE (PotentialInfo) ,PRIVATE, POINTER :: PotPtr
  TYPE (OrbitInfo) ,PRIVATE, POINTER :: OrbitPtr
  TYPE (SCFInfo) ,PRIVATE, POINTER :: SCFPtr

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Init_SCFatom						!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Init_SCFatom(ifinput)
    INTEGER, INTENT(IN),OPTIONAL :: ifinput

    !  program to calculate the self-consistent density functional
    !    atom ground state for atom with atomic number nz
    !    for self-consistent potential rv


    REAL(8) :: xocc,qf,small,range,zeff,hval,r0
    REAL(8) :: qcal, rescale
    INTEGER :: icount,i,j,k,it,start,np,ierr,gridpoints
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir,nzeff
    INTEGER, ALLOCATABLE :: nl(:,:)
    CHARACTER(128) :: exchangecorrelationandgridline,gridkey,relkey
    CHARACTER(132) :: inputline

    scalarrelativistic=.FALSE.; finitenucleus=.FALSE. ; 
    frozencorecalculation=.FALSE.;setupfrozencore=.false.
    gaussianshapefunction=.FALSE.;besselshapefunction=.FALSE.
    ColleSalvetti=.FALSE.  ; HFpostprocess=.FALSE.
    localizedcoreexchange=.false.

    WRITE(6,*) 'enter atomic symbol and atomic number'
    IF (PRESENT(ifinput)) THEN
       READ(5,'(a)') inputline
       WRITE(ifinput,'(a)') TRIM(inputline)
       READ(inputline,*) AEPot%sym,nz
    ELSE
       READ(5,*) AEPot%sym,nz
    ENDIF

    AEPot%nz=nz
    WRITE(6,*) 'exchange-correlation type -- LDA-PW(default) or GGA-PBE '
    WRITE(6,*) 'optionally (space) "Fixed_Zero" keyword followed by',&
         ' (space) zero_index  if appropriate'
    WRITE(6,*) 'optionally (space) "loggrid" keyword if appropriate'
    WRITE(6,*) 'further optionally (space) n for number of log grid points'
    WRITE(6,*) 'further optionally (space) "scalarrelativistic" keyword'
    WRITE(6,*) 'further optionally (space) "finite-nucleus" keyword'
    READ(5,'(a)') exchangecorrelationandgridline
    IF(PRESENT(ifinput)) WRITE(ifinput,'(a)')&
         TRIM(exchangecorrelationandgridline)
    CALL Uppercase(exchangecorrelationandgridline)
    exchangecorrelationandgridline=TRIM(exchangecorrelationandgridline)

    write(6,*) 'writing exchangecorrelationandgridline',&
      exchangecorrelationandgridline

    READ(unit=exchangecorrelationandgridline,fmt=*) exctype
    write(6,*) exchangecorrelationandgridline
    write(6,*) exctype
    IF (TRIM(exctype)=='EXX'.or.TRIM(exctype)=='EXXKLI'.or.   &
           TRIM(exctype)=='EXXOCC') THEN
       CALL EXX_Input_Settings(exchangecorrelationandgridline)
    ELSE IF (TRIM(exctype)=='EXXCS') THEN
       CALL EXX_Input_Settings(exchangecorrelationandgridline)
       ColleSalvetti=.TRUE.
    ELSE IF (TRIM(exctype)=='HF'.or.TRIM(exctype)=='HFV') THEN
    ELSE
       CALL initexch
    ENDIF
    AEOrbit%exctype=TRIM(exctype)


    i=INDEX(exchangecorrelationandgridline,'LOGGRID')
    IF (i>0) THEN
       gridkey='LOGGRID'
       READ(unit=exchangecorrelationandgridline(i+7:128),fmt=*) gridpoints
    ELSE
       gridkey='LINEAR'
       gridpoints=mxgrid
       hadjusted=linh
    END IF
    i=INDEX(exchangecorrelationandgridline,'LINEARGRID')
    IF (i>0) THEN
       gridkey='LINEAR'
       READ(unit=exchangecorrelationandgridline(i+10:128),fmt=*) gridpoints
       hadjusted=linrange/(gridpoints-1)
    END IF
    write(6,*) 'Grid info ', gridkey,gridpoints

    i=INDEX(exchangecorrelationandgridline,'SCALARRELATIVISTIC')
    IF (i>0) scalarrelativistic=.TRUE.
    i=INDEX(exchangecorrelationandgridline,'FINITE-NUCLEUS')
    IF (i>0) finitenucleus=.TRUE.
    i=INDEX(exchangecorrelationandgridline,'HFPOSTPROCESS')
    IF (i>0) HFpostprocess=.TRUE.
    i=INDEX(exchangecorrelationandgridline,'LOCALIZEDCOREEXCHANGE')
    IF (i>0) localizedcoreexchange=.TRUE.

    write(6,&
     '(" logical variables scalarrelativistic,finitenucleus,HFpostprocess   - ", 3L3)') &
      scalarrelativistic,finitenucleus,HFpostprocess

    WRITE(6,*) 'Calculation for atomic number = ',nz
    WRITE(6,*) 'enter maximum principal quantum numbers for s,p,d,f,g'
    IF(PRESENT(ifinput)) THEN
       READ(5,'(a)') inputline
       WRITE(ifinput,'(a)') TRIM(inputline)
       READ(inputline,*) nps,npp,npd,npf,npg
    ELSE
       READ(5,*) nps,npp,npd,npf,npg
    ENDIF
    IF(nps<0)nps=0
    IF(npp<0)npp=0
    IF(npd<0)npd=0
    IF(npf<0)npf=0
    IF(npg<0)npg=0
    WRITE(6,'(5i4)') nps,npp,npd,npf,npg
    i=MAX(nps,npp,npd,npf,npg)
    !j=nps+npp+npd+npf+npg
    j=nps
    IF(npp>0) j=j+npp-1
    IF(npd>0) j=j+npd-2
    IF(npf>0) j=j+npf-3
    IF(npg>0) j=j+npg-4


    ALLOCATE(nl(i,j),AEOrbit%np(j),AEOrbit%l(j),&
         AEOrbit%eig(j),AEOrbit%occ(j),AEOrbit%iscore(j),stat=k)
    IF (k/=0) THEN
       WRITE(6,*) 'Error in allocation of nl,np...',k
       STOP
    ENDIF
    nl=0;AEOrbit%np=0;AEOrbit%l=0;AEOrbit%eig=0;AEOrbit%occ=0
    AEOrbit%nps=nps;AEOrbit%npp=npp;AEOrbit%npd=npd;AEOrbit%npf=npf
    AEOrbit%npg=npg;AEOrbit%iscore=.false.
    
    icount=0
    IF (nps.GT.0) THEN
       DO is=1,nps
          icount=icount+1
          nl(is,1)=icount
          AEOrbit%occ(icount)=2.d0
          AEOrbit%np(icount)=is
          AEOrbit%l(icount)=0
       ENDDO
    ENDIF
    IF (npp.GT.1) THEN
       DO ip=2,npp
          icount=icount+1
          nl(ip,2)=icount
          AEOrbit%occ(icount)=6.d0
          AEOrbit%np(icount)=ip
          AEOrbit%l(icount)=1
       ENDDO
    ENDIF
    IF (npd.GT.2) THEN
       DO id=3,npd
          icount=icount+1
          nl(id,3)=icount
          AEOrbit%occ(icount)=10.d0
          AEOrbit%np(icount)=id
          AEOrbit%l(icount)=2
       ENDDO
    ENDIF
    IF (npf.GT.3) THEN
       DO jf=4,npf
          icount=icount+1
          nl(jf,4)=icount
          AEOrbit%occ(icount)=14.d0
          AEOrbit%np(icount)=jf
          AEOrbit%l(icount)=3
       ENDDO
    ENDIF
    IF(npg.GT.4) THEN
       DO ig=5,npg
          icount=icount+1
          nl(ig,5)=icount
          AEOrbit%occ(icount)=18.d0
          AEOrbit%np(icount)=ig
          AEOrbit%l(icount)=4
       ENDDO
    ENDIF
    norbit=icount
    AEOrbit%nps=nps
    AEOrbit%npp=npp
    AEOrbit%npd=npd
    AEOrbit%npf=npf
    AEOrbit%npg=npg
    AEOrbit%norbit=norbit
    WRITE(6,*) norbit, ' orbitals will be calculated'
    !
    WRITE(6,*)' Below are listed the default occupations '
    WRITE(6,"(' n  l     occupancy')")
    DO io=1,norbit
       WRITE(6,'(i2,1x,i2,4x,1pe15.7)') &
            AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
    ENDDO
    !
    WRITE(6,*)' enter np l occ for all occupations for all revisions'
    WRITE(6,*)'  enter 0 0 0. to end'

    DO
       IF(PRESENT(ifinput)) THEN
          READ(5,'(a)') inputline
          WRITE(ifinput,'(a)') TRIM(inputline)
          READ(inputline,*) ip,l,xocc
       ELSE
          READ(5,*) ip,l,xocc
       ENDIF
       IF (ip.LE.0) EXIT
       nfix=nl(ip,l+1)
       IF (nfix.LE.0.OR.nfix.GT.norbit) THEN
          WRITE(6,*) 'error in occupations -- ip,l,xocc',                  &
               ip,l,xocc,nfix,norbit
          STOP
       ENDIF
       AEOrbit%occ(nfix)=xocc
    ENDDO

    !
    WRITE(6,*) ' Corrected occupations are: '
    WRITE(6,"(' n  l     occupancy')")
    electrons=0.d0
    DO io=1,norbit
       WRITE(6,'(i2,1x,i2,4x,1pe15.7)')  &
            AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
       electrons=electrons+AEOrbit%occ(io)
    ENDDO
    AEPot%q=electrons
    qf=nz-electrons
    WRITE(6,*)
    WRITE(6,*) 'nuclear charge    = ' , nz
    WRITE(6,*) 'electronic charge = ', electrons
    WRITE(6,*) 'net charge        = ', qf


    IF (TRIM(gridkey)=='LOGGRID') THEN
       hval=logh
       !IF (gridpoints>0) CALL findh(nz,logrange,gridpoints,hval,r0)
       !CALL Init_grid(Grid,hval,logrange,r0)
       IF (gridpoints>0) CALL findh_given_r0(nz,logrange,lor00,gridpoints,hval)
       r0=lor00/nz
       CALL Init_grid(Grid,hval,logrange,r0)
    ELSE
       CALL Init_grid(Grid,hadjusted,linrange)
    ENDIF

    write(6,*) 'Finish Init_SCFatom' ; call flush(6)
    IF (scalarrelativistic) CALL Allocate_Scalar_Relativistic(Grid)

    DEALLOCATE(nl)

  END SUBROUTINE Init_SCFatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  SCFatom                              !!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SCFatom(scftype,lotsofoutput,ifinput)
    IMPLICIT NONE 
    CHARACTER(2),INTENT(IN)::scftype
    LOGICAL, INTENT(IN) :: lotsofoutput
    INTEGER, INTENT(IN), OPTIONAL :: ifinput
    TYPE(SCFInfo) :: SCFPP


    !  program to perform self-consistency loop

    IF(scftype=='AE') THEN
       !------------AE calculation--------------
       frozencorecalculation=.FALSE.
       setupfrozencore=.FALSE.
       PotPtr=>AEPot
       OrbitPtr=>AEOrbit
       SCFPtr=>AESCF

       CALL Init_Orbit()
       CALL Init_Potential()

    ELSEIF(scftype=='NC') THEN
       !------------New Configuration Calculation---------
       frozencorecalculation=.FALSE.
       setupfrozencore=.FALSE.
       PotPtr=>AEPot
       OrbitPtr=>AEOrbit
       SCFPtr=>AESCF

       CALL NC_Init
       Call Init_Potential

    ELSEIF(scftype=='SC') THEN
       !------------Set Core and Valence in current config.---------

       frozencorecalculation=.TRUE.
       setupfrozencore=.TRUE.
       CALL Copy_PotentialInfo(AEPot,FCPot)
       PotPtr=>FCPot

       CALL SC_Init(ifinput)
       OrbitPtr=>FCOrbit
       SCFPtr=>FCSCF

    ELSEIF(scftype=='FC') THEN
       !------------Frozen Core Calculation---------
       frozencorecalculation=.TRUE.
       setupfrozencore=.FALSE.
       PotPtr=>FCPot
       OrbitPtr=>FCOrbit
       SCFPtr=>FCSCF

       CALL FC_Init
       Call Init_Potential

    ENDIF


    IF (TRIM(OrbitPtr%exctype)=='EXX'.OR. &
          TRIM(OrbitPtr%exctype)=='EXXKLI'.OR. &
          TRIM(OrbitPtr%exctype)=='EXXCS') THEN

       CALL EXX_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)
    ELSE IF (TRIM(OrbitPtr%exctype)=='EXXOCC') THEN
       CALL EXXOCC_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)

    ELSE IF (TRIM(OrbitPtr%exctype)=='HF'.or.TRIM(OrbitPtr%exctype)=='HFV') THEN

        write(6,*) 'Just before HF_SCF'; call flush(6)
       CALL HF_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)
        write(6,*) 'Just after HF_SCF'; call flush(6)

    ELSE
       ! LDA or GGA
        write(6,*) 'before ldagga '; call flush(6)
       CALL LDAGGA_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)
    ENDIF

    If (HFpostprocess) then
       write(6,*) "******PostProcessing HF*********"
       call hf_energy_only(Grid,OrbitPtr,PotPtr,SCFPP)
    endif
  END SUBROUTINE SCFatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Init_Orbit				
  !!       From nuclear charge -- generate hydrogenic-like initial wfns
  !!          and densities -- 
  !!          fill AEOrbit%wfn, AEOrbit%eig, and AEOrbit%den and AEOrbit%q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Init_Orbit()
    IMPLICIT NONE
    INTEGER  :: i,io,ir,ip,l,nfix,j,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc
    INTEGER :: initialconfig=0

    if (initialconfig/=0) then
       write(6,*) 'Error in aeatom -- Init_Orbit already called'
       stop
    endif

       n=Grid%n
       j=AEOrbit%norbit

       ALLOCATE(AEOrbit%wfn(n,j))
       AEOrbit%wfn=0; 

       !  calculate initial charge density from hydrogen-like functions
       !    also initial energies
       !
       AEOrbit%wfn=0

       zeff=nz
       DO io=1,norbit
          np=AEOrbit%np(io)
          l=AEOrbit%l(io)
          xocc=AEOrbit%occ(io)
          AEOrbit%eig(io)=-(zeff/(np))**2
          WRITE(6,*) io,np,l,xocc,AEOrbit%eig(io)
          DO ir=1,n
             AEOrbit%wfn(ir,io)=hwfn(zeff,np,l,Grid%r(ir))
             IF (ABS(AEOrbit%wfn(ir,io))<machine_zero)AEOrbit%wfn(ir,io)=0.d0
          ENDDO
          zeff=zeff-0.5d0*xocc
          zeff=MAX(zeff,1.d0)
       ENDDO

       ALLOCATE(AEOrbit%den(n))
       ! check charge and rescale 
       AEOrbit%den(1:n)=0.d0
       DO io=1,norbit
          xocc=AEOrbit%occ(io)
          DO ir=1,n
             AEOrbit%den(ir)=AEOrbit%den(ir)+xocc*(AEOrbit%wfn(ir,io)**2)
          ENDDO
       ENDDO
       qcal=integrator(Grid,AEOrbit%den)
       qf=qcal
       !WRITE(6,*) 'qcal electrons = ',qcal, electrons
       !  rescale density
       rescale=electrons/qcal
       AEOrbit%den(1:n)=AEOrbit%den(1:n)*rescale

       initialconfig=1
        write(6,*) 'completed Init_Orbit '; call flush(6)

  END SUBROUTINE Init_Orbit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE NC_Init
    INTEGER  :: i,io,ir,ip,l,nfix,j,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc
       !----------New Configuration--AE calculation--------------
       ! readin revision and normalize the density

       WRITE(6,*)' enter np l occ for all occupations for all revisions'
       WRITE(6,*)'  enter 0 0 0. to end'

       DO
          READ(5,*) ip,l,xocc
          IF (ip.LE.0) EXIT
          nfix=-100
          DO io=1,norbit
             IF (ip==AEOrbit%np(io).AND.l==AEOrbit%l(io)) THEN
                nfix=io
                EXIT
             ENDIF
          ENDDO
          IF (nfix.LE.0.OR.nfix.GT.norbit) THEN
             WRITE(6,*) 'error in occupations -- ip,l,xocc',                  &
                  ip,l,xocc,nfix,norbit
             STOP
          ENDIF
          AEOrbit%occ(nfix)=xocc
       ENDDO

       WRITE(6,*) ' Corrected occupations are: '
       WRITE(6,"(' n  l     occupancy')")
       electrons=0.d0
       DO io=1,norbit
          WRITE(6,'(i2,1x,i2,4x,1pe15.7)')&
               AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
          electrons=electrons+AEOrbit%occ(io)
       ENDDO
       AEPot%q=electrons
       qf=nz-electrons
       WRITE(6,*)
       WRITE(6,*) 'nuclear charge    = ' , nz
       WRITE(6,*) 'electronic charge = ', electrons
       WRITE(6,*) 'net charge        = ', qf
       !
       !
       !  calculate initial charge density from stored wavefunctions
       !    also initial energies
       !
       AEOrbit%den(1:n)=0.d0
       DO io=1,norbit
          xocc=AEOrbit%occ(io)
          DO ir=1,n
             AEOrbit%den(ir)=AEOrbit%den(ir)+xocc*(AEOrbit%wfn(ir,io)**2)
          ENDDO
       ENDDO
       !
       !  check charge
       !
       qcal=integrator(Grid,AEOrbit%den)
       qf=qcal
       !WRITE(6,*) 'qcal electrons = ',qcal, electrons
       !  rescale density
       rescale=electrons/qcal
       AEOrbit%den(1:n)=AEOrbit%den(1:n)*rescale

  END SUBROUTINE NC_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!         SC_Init -- set core states in current configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SC_Init(ifinput)
    INTEGER :: choosevalence=0
    INTEGER, INTENT(IN), OPTIONAL :: ifinput
    
    INTEGER :: io
    LOGICAL :: noalt=.true.
    If (choosevalence/=0) then
     write(6,*) 'Error in aeatom -- SC_Init already called'
       stop
    endif

          CALL Copy_OrbitInfo(AEOrbit,FCOrbit)
          CALL Set_Valence(ifinput)
         
          Call Core_Electron_Report(FCOrbit,FC,6)
          Call Valence_Electron_Report(FCOrbit,FC,6)

        ! determine valence energy
          CALL Copy_SCFInfo(AESCF,FCSCF)

          IF (TRIM(FCOrbit%exctype)=='EXX') then
             Call Get_FCKinCoul(Grid,FCPot,FCOrbit,FC,FCSCF)
             CALL Get_FCEnergy_EXX(Grid,FCOrbit,FC,FCSCF)
             CALL Set_Vxref(AEPot)   ! probably not a good idea
          ELSEIF (TRIM(FCOrbit%exctype)=='EXXKLI') then
             Call Get_FCKinCoul(Grid,FCPot,FCOrbit,FC,FCSCF,noalt)
             CALL Get_FCEnergy_EXX(Grid,FCOrbit,FC,FCSCF)
          ELSEIF (TRIM(FCOrbit%exctype)=='HF'.or.TRIM(FCOrbit%exctype)=='HFV')&
                     then
             write(6,*) ' Completed Setcore '; call flush(6)
             CALL Get_FCEnergy_HF(Grid,FCPot,FCOrbit,FC,FCSCF)
             CALL Total_FCEnergy_Report(FCSCF,6)
             write(6,*) ' Completed Setcore '; call flush(6)
          ELSE
             Call Get_FCKinCoul(Grid,FCPot,FCOrbit,FC,FCSCF)
             CALL Get_FCEXC(FCSCF)
             CALL Total_FCEnergy_Report(FCSCF,6)
          ENDIF
          !CALL Total_Energy_Report(FCSCF,6)   
              !   For EXX, exchange contribution is not known yet.

       choosevalence=1

  END SUBROUTINE SC_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FC_Init
    INTEGER  :: i,io,ir,ip,l,nfix,j,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc

       WRITE(6,*)' Below are listed the current valence occupations '
       WRITE(6,"(' n  l     occupancy')")
       DO io=1,norbit
          IF (.NOT.FCOrbit%iscore(io)) WRITE(6,'(i2,1x,i2,4x,1pe15.7)') &
               FCOrbit%np(io),FCOrbit%l(io),FCOrbit%occ(io)
       ENDDO
       !
       WRITE(6,*)' enter np l occ for all occupations for all revisions'
       WRITE(6,*)'  enter 0 0 0. to end'

       DO
          READ(5,*) ip,l,xocc
          IF (ip.LE.0) EXIT
          nfix=-100
          DO io=1,norbit
             IF (ip==FCOrbit%np(io).AND.l==FCOrbit%l(io)&
                  .AND.(.NOT.FCOrbit%iscore(io))) THEN
                nfix=io
                EXIT
             ENDIF
          ENDDO
          IF (nfix.LE.0.OR.nfix.GT.norbit) THEN
             WRITE(6,*) 'error in occupations -- ip,l,xocc',                  &
                  ip,l,xocc,nfix,norbit
             STOP
          ENDIF
          FCOrbit%occ(nfix)=xocc
       ENDDO

       !
       WRITE(6,*) ' Corrected occupations are: '
       WRITE(6,"(' n  l     occupancy')")

       FC%zvale=0.d0
       DO io=1,norbit
          If (.not.FCOrbit%iscore(io))then
            WRITE(6,'(i2,1x,i2,4x,1pe15.7)')  &
              FCOrbit%np(io),FCOrbit%l(io),FCOrbit%occ(io)
              FC%zvale=FC%zvale+FCOrbit%occ(io)
          Endif
       ENDDO

       electrons=FC%zvale+FC%zcore

       !
       !  calculate initial charge density from stored wavefunctions
       !    also initial energies
       !
       FC%valeden=0
       DO io=1,norbit
          IF (.NOT.FCOrbit%iscore(io)) THEN
             xocc=FCOrbit%occ(io)
             DO ir=1,n
                FC%valeden(ir)=FC%valeden(ir)+xocc*(FCOrbit%wfn(ir,io)**2)
             ENDDO
          ENDIF
       ENDDO
       !
       !  check charge				!  rescale density

       !
       qcal=integrator(Grid,FC%valeden)
       if (ABS(qcal)<small0.OR.electrons<small0) then
          FC%valeden=0
       else
          rescale=FC%zvale/qcal
          FC%valeden(1:n)=FC%valeden(1:n)*rescale
       endif

    FCOrbit%den=FC%valeden+FC%coreden
       WRITE(6,*)
       WRITE(6,*) 'core charge    = ' , FC%zcore
       WRITE(6,*) 'valence charge = ',  FC%zvale
       WRITE(6,*) 'net charge     = ',  nz-FC%zcore-FC%zvale

  END SUBROUTINE FC_Init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Init_Potential				
  !!    Generate Potptr%rv, rvh, rvx given AEOrbit%wfn, AEOrbit%den, AEOrbit%q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Init_Potential()
    IMPLICIT NONE
    INTEGER  :: i,io,ir,xocc,ip,l,nfix,j,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,ecoul,v0,etxc,eex

    LOGICAL :: success
    INTEGER :: firsttime=0
    CHARACTER(4) :: stuff


    write(6,*) 'In Init_Potential'; call flush(6)
    IF (firsttime==0) THEN

          !------------AE calculation--------------
          ALLOCATE(AEPot%rv(Grid%n), AEPot%rvh(Grid%n),&
               AEPot%rvn(Grid%n),AEPot%rvx(Grid%n),stat=i)
          ! rvh == current hartree potential for this den
          IF ( i /= 0) THEN
             WRITE(6,*) 'InitPotential allocation error ', i,Grid%n
             STOP
          ENDIF
          CALL Get_Nuclearpotential(Grid,AEPot)
   ENDIF
   firsttime=firsttime+1

    CALL poisson(Grid,PotPtr%q,OrbitPtr%den,PotPtr%rvh,ecoul,v0)     
    write(6,*) 'In Init_Potential', PotPtr%q,ecoul; call flush(6)

  END SUBROUTINE Init_Potential



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Set Valence Orbitals
  !      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Set_Valence(ifinput)
    INTEGER, INTENT(IN), OPTIONAL :: ifinput

    CHARACTER(1) :: answer
    INTEGER :: io,n,ok,norbit
    CHARACTER(132) :: inputline

    norbit=FCOrbit%norbit
    n=Grid%n
    ALLOCATE(FC%coreden(n),FC%valeden(n),stat=ok)
    IF (ok /= 0) THEN
       WRITE(6,*) 'Error in FC%coreden allocation',n
       STOP
    ENDIF

    WRITE(6,*) 'for each state enter c for core or v for valence'
    FC%zcore=0; FC%zvale=0
    FC%coreden=0; FC%valeden=0
    DO io=1,norbit
       WRITE(6,'(3i5,1p2e15.7)') io,FCOrbit%np(io),FCOrbit%l(io),&
            FCOrbit%occ(io),FCOrbit%eig(io)
       DO
          IF (PRESENT(ifinput)) THEN
             READ(5,'(a)') inputline
             WRITE(ifinput,'(a)') TRIM(inputline)
             READ(inputline,*) answer
          ELSE
             READ(5,*) answer
          ENDIF
          IF (answer.NE.'c'.AND.answer.NE.'C'.AND.answer.NE.'v'&
               .AND.answer.NE.'V')  THEN
             WRITE(6,*) 'Please input c or v'
          ELSE
             EXIT
          ENDIF
       ENDDO
       IF (answer.EQ.'c'.OR.answer.EQ.'C') FCOrbit%iscore(io)=.TRUE.
       IF (answer.EQ.'v'.OR.answer.EQ.'V') FCOrbit%iscore(io)=.FALSE.
       IF (FCOrbit%iscore(io)) THEN
          FC%zcore=FC%zcore+FCOrbit%occ(io)
          FC%coreden=FC%coreden+FCOrbit%occ(io)*(FCOrbit%wfn(:,io))**2
       ENDIF
       IF (.NOT.FCOrbit%iscore(io)) THEN

          FC%zvale=FC%zvale+FCOrbit%occ(io)
          FC%valeden=FC%valeden+FCOrbit%occ(io)*(FCOrbit%wfn(:,io))**2
       ENDIF
    ENDDO

         write(6,*) 'Returning from SetValence'
         write(6,*) 'Core electrons', FC%zcore,integrator(Grid,FC%coreden) 
         write(6,*) 'Vale electrons', FC%zvale,integrator(Grid,FC%valeden) 

  END SUBROUTINE Set_Valence

  SUBROUTINE dump_aeatom(Fn,Grid,AEOrbit,AEPot,AESCF,FCOrbit,FCPot,FCSCF,FC)
     Character(*), INTENT(IN) :: Fn
     Type(GridInfo), INTENT(IN) :: Grid
     Type(OrbitInfo), INTENT(IN) :: AEOrbit,FCOrbit
     Type(PotentialInfo), INTENT(IN) :: AEPot,FCPot
     Type(SCFInfo), INTENT(IN) :: AESCF,FCSCF
     Type(FCInfo), INTENT(IN) :: FC

     INTEGER, parameter :: ifo=15
     INTEGER :: i,io,n,norbit

     open(ifo,file=Fn,form='formatted',status='replace')
     write(6,*) 'Creating or replacing dump file', Fn

     write(ifo,*) 'ATOMPAW'    ! keyword

    ! grid info
     write(ifo,*) Grid%TYPE,Grid%n,Grid%h
     n=Grid%n
     write(ifo,*)(Grid%r(i),i=1,n)
     if (usingloggrid(Grid)) then
          write(ifo,*)(Grid%drdu(i),Grid%pref(i),Grid%rr02(i),i=1,n)
     endif

    ! atomdata fixed constants
     write(ifo,*)  frozencorecalculation,setupfrozencore,scalarrelativistic ,&
      finitenucleus,gaussianshapefunction,besselshapefunction,ColleSalvetti

    ! Orbit info
     write(ifo,*) AEOrbit%exctype,AEOrbit%nps, AEOrbit%npp, AEOrbit%npd ,&
       AEOrbit%npf, AEOrbit%npg, AEOrbit%norbit
     norbit=AEOrbit%norbit
     write(ifo,*) (AEOrbit%np(io),AEOrbit%l(io),AEOrbit%iscore(io),&
         AEOrbit%eig(io),AEOrbit%occ(io),io=1,norbit)
     write(ifo,*) ((AEOrbit%wfn(i,io),i=1,n),io=1,norbit)
     write(ifo,*) (AEOrbit%den(i),i=1,n)

     write(ifo,*) FCOrbit%exctype,FCOrbit%nps, FCOrbit%npp, FCOrbit%npd ,&
       FCOrbit%npf, FCOrbit%npg, FCOrbit%norbit
     norbit=FCOrbit%norbit
     write(ifo,*) (FCOrbit%np(io),FCOrbit%l(io),FCOrbit%iscore(io),&
         FCOrbit%eig(io),FCOrbit%occ(io),io=1,norbit)
     write(ifo,*) ((FCOrbit%wfn(i,io),i=1,n),io=1,norbit)
     write(ifo,*) (FCOrbit%den(i),i=1,n)

   ! Pot info
      write(ifo,*) AEPot%nz,AEPot%sym,AEPot%q,AEPot%v0,AEPot%v0p
      write(ifo,*) (AEPot%rv(i),AEPot%rvn(i),AEPot%rvh(i),AEPot%rvx(i),i=1,n)

      write(ifo,*) FCPot%nz,FCPot%sym,FCPot%q,FCPot%v0,FCPot%v0p
      write(ifo,*) (FCPot%rv(i),FCPot%rvn(i),FCPot%rvh(i),FCPot%rvx(i),i=1,n)

   ! SCFinfo
      write(ifo,*) AESCF%iter,AESCF%delta,AESCF%eone,AESCF%ekin,&
        AESCF%estatic,AESCF%ecoul,AESCF%eexc,AESCF%oepcs,AESCF%etot,&
        AESCF%valekin,AESCF%valecoul,AESCF%valeexc,AESCF%corekin,AESCF%evale

      write(ifo,*) FCSCF%iter,FCSCF%delta,FCSCF%eone,FCSCF%ekin,&
        FCSCF%estatic,FCSCF%ecoul,FCSCF%eexc,FCSCF%oepcs,FCSCF%etot,&
        FCSCF%valekin,FCSCF%valecoul,FCSCF%valeexc,FCSCF%corekin,FCSCF%evale

   ! FCinfo
      write(ifo,*) FC%zvale,FC%zcore
      write(ifo,*) (FC%coreden(i),FC%valeden(i),i=1,n)

   If (TRIM(AEOrbit%exctype)=='EXX'.or.TRIM(AEOrbit%exctype)=='EXXKLI') &
                   Call EXXdump(Grid,AEOrbit,ifo)
   If (TRIM(AEOrbit%exctype)=='HF'.or.TRIM(AEOrbit%exctype)=='HFV') &
           Call HFdump(Grid,AEOrbit,ifo)

   close(ifo)

   write(6,*) 'Closing dump file'
  END SUBROUTINE dump_aeatom

  SUBROUTINE load_aeatom(Fn,Grid,AEOrbit,AEPot,AESCF,FCOrbit,FCPot,FCSCF,FC,ifd)
     Character(*), INTENT(IN) :: Fn
     Type(GridInfo), INTENT(OUT) :: Grid
     Type(OrbitInfo), INTENT(OUT) :: AEOrbit,FCOrbit
     Type(PotentialInfo), INTENT(OUT) :: AEPot,FCPot
     Type(SCFInfo), INTENT(OUT) :: AESCF,FCSCF
     Type(FCInfo), INTENT(OUT) :: FC
     INTEGER, INTENT(IN) :: ifd

     INTEGER, parameter :: ifo=15
     INTEGER :: i,io,n,norbit,j,k
     CHARACTER(132) :: inputfile,checkfile

     open(ifo,file=Fn,form='formatted',status='old')
     write(6,*) 'Reading dump file ', Fn

     read(ifo,*) inputfile(1:7)
     if (inputfile(1:7)=='ATOMPAW') then
       write(6,*) 'File seems to be fine'
     else
       write(6,*) 'File is not correct -- program will stop'
       stop
     endif

    ! grid info
     read(ifo,*) Grid%TYPE,Grid%n,Grid%h
     n=Grid%n
     Allocate(Grid%r(n))
     read(ifo,*)(Grid%r(i),i=1,n)
     if (usingloggrid(Grid)) then
        Allocate(Grid%drdu(n),Grid%pref(n),Grid%rr02(n))
          read(ifo,*)(Grid%drdu(i),Grid%pref(i),Grid%rr02(i),i=1,n)
     endif

    ! atomdata fixed constants
     read(ifo,*)  frozencorecalculation,setupfrozencore,scalarrelativistic ,&
      finitenucleus,gaussianshapefunction,besselshapefunction,ColleSalvetti

    ! Orbit info
     read(ifo,*) AEOrbit%exctype,AEOrbit%nps, AEOrbit%npp, AEOrbit%npd ,&
       AEOrbit%npf, AEOrbit%npg, AEOrbit%norbit
     norbit=AEOrbit%norbit
     Allocate(AEOrbit%np(norbit),AEOrbit%l(norbit),AEOrbit%iscore(norbit),&
         AEOrbit%eig(norbit),AEOrbit%occ(norbit))
     read(ifo,*) (AEOrbit%np(io),AEOrbit%l(io),AEOrbit%iscore(io),&
         AEOrbit%eig(io),AEOrbit%occ(io),io=1,norbit)
     Allocate(AEOrbit%wfn(n,norbit),AEOrbit%den(n))
     read(ifo,*) ((AEOrbit%wfn(i,io),i=1,n),io=1,norbit)
     read(ifo,*) (AEOrbit%den(i),i=1,n)

     read(ifo,*) FCOrbit%exctype,FCOrbit%nps, FCOrbit%npp, FCOrbit%npd ,&
       FCOrbit%npf, FCOrbit%npg, FCOrbit%norbit
     norbit=FCOrbit%norbit
     Allocate(FCOrbit%np(norbit),FCOrbit%l(norbit),FCOrbit%iscore(norbit),&
         FCOrbit%eig(norbit),FCOrbit%occ(norbit))
     read(ifo,*) (FCOrbit%np(io),FCOrbit%l(io),FCOrbit%iscore(io),&
         FCOrbit%eig(io),FCOrbit%occ(io),io=1,norbit)
     Allocate(FCOrbit%wfn(n,norbit),FCOrbit%den(n))
     read(ifo,*) ((FCOrbit%wfn(i,io),i=1,n),io=1,norbit)
     read(ifo,*) (FCOrbit%den(i),i=1,n)

     exctype=AEOrbit%exctype
     IF(TRIM(exctype)/='EXX'.and.TRIM(exctype)/='EXXCS'.and.&
                  TRIM(exctype)/='EXXKLI') CALL initexch

   ! Pot info
      read(ifo,*) AEPot%nz,AEPot%sym,AEPot%q,AEPot%v0,AEPot%v0p
      Allocate(AEPot%rv(n),AEPot%rvn(n),AEPot%rvh(n),AEPot%rvx(n))
      read(ifo,*) (AEPot%rv(i),AEPot%rvn(i),AEPot%rvh(i),AEPot%rvx(i),i=1,n)

      read(ifo,*) FCPot%nz,FCPot%sym,FCPot%q,FCPot%v0,FCPot%v0p
      Allocate(FCPot%rv(n),FCPot%rvn(n),FCPot%rvh(n),FCPot%rvx(n))
      read(ifo,*) (FCPot%rv(i),FCPot%rvn(i),FCPot%rvh(i),FCPot%rvx(i),i=1,n)

   ! SCFinfo
      read(ifo,*) AESCF%iter,AESCF%delta,AESCF%eone,AESCF%ekin,&
        AESCF%estatic,AESCF%ecoul,AESCF%eexc,AESCF%oepcs,AESCF%etot,&
        AESCF%valekin,AESCF%valecoul,AESCF%valeexc,AESCF%corekin,AESCF%evale

      read(ifo,*) FCSCF%iter,FCSCF%delta,FCSCF%eone,FCSCF%ekin,&
        FCSCF%estatic,FCSCF%ecoul,FCSCF%eexc,FCSCF%oepcs,FCSCF%etot,&
        FCSCF%valekin,FCSCF%valecoul,FCSCF%valeexc,FCSCF%corekin,FCSCF%evale

   ! FCinfo
      read(ifo,*) FC%zvale,FC%zcore
      allocate(FC%coreden(n),FC%valeden(n))
      read(ifo,*) (FC%coreden(i),FC%valeden(i),i=1,n)

   If (TRIM(AEOrbit%exctype)=='EXX'.or.TRIM(AEOrbit%exctype)=='EXXKLI') &
                    Call EXXload(Grid,AEOrbit,ifo)
   If (TRIM(AEOrbit%exctype)=='HF'.or.TRIM(AEOrbit%exctype)=='HFV') &
            Call HFload(Grid,AEOrbit,ifo)

   close(ifo)

   write(6,*) 'Closing load file'

  ! crude checking by comparing input file (unit 5) with dummy file (unit ifd)
   rewind(5); rewind(ifd); 
    Do i=1,3
      read(5,'(a)') inputfile; read(ifd,'(a)') checkfile
      If (TRIM(inputfile)/=TRIM(checkfile)) then
         write(6,*) 'Inconsistency in files'
          write(6,'(a)') inputfile
          write(6,'(a)') checkfile
          stop
      endif
    enddo

    Bigloop: Do
      read(5,'(a)') inputfile; read(ifd,'(a)') checkfile
      read(inputfile,*) i,j,k
      if (i==0) exit Bigloop
      If (TRIM(inputfile)/=TRIM(checkfile)) then
         write(6,*) 'Inconsistency in files'
          write(6,'(a)') inputfile
          write(6,'(a)') checkfile
          stop
      endif
    enddo Bigloop
 
    Do io=1,norbit
      read(5,'(a)') inputfile; read(ifd,'(a)') checkfile
      If (TRIM(inputfile)/=TRIM(checkfile)) then
         write(6,*) 'Inconsistency in files'
          write(6,'(a)') inputfile
          write(6,'(a)') checkfile
          stop
      endif
    enddo

    write(6,*) 'Load completed normally'
  END SUBROUTINE load_aeatom
END  MODULE AEatom
