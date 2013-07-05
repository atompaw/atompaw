MODULE exx_mod
  USE anderson_driver
  USE atomdata
  USE exxdata
  USE fock
  USE general_mod
  USE globalmath
  USE gridmod
  USE hf_mod
  USE report_mod

  IMPLICIT NONE

  TYPE (Gridinfo) ,PRIVATE, POINTER :: Gridwk
  TYPE (PotentialInfo) ,PRIVATE, POINTER :: Potwk
  TYPE (OrbitInfo) ,PRIVATE, POINTER :: Orbitwk
  TYPE (OrbitInfo) ,PRIVATE, POINTER  :: EigOrbitwk
  TYPE (FCInfo) ,PRIVATE, POINTER :: FCwk
  TYPE (SCFInfo) ,PRIVATE, POINTER :: SCFwk

  REAL(8), PRIVATE, ALLOCATABLE :: geqn(:,:)

  TYPE(Anderson_context), PRIVATE :: AC
  LOGICAL, PRIVATE :: verboseoutput

  LOGICAL, PRIVATE :: zeroadjust

CONTAINS

  SUBROUTINE EXX_SCF(scftype,lotsofoutput,Gridin,Orbitin,Potin,FCin,SCFin)
    CHARACTER(2), INTENT(IN) :: scftype
    LOGICAL, INTENT(IN) :: lotsofoutput
    TYPE(GridInfo) ,TARGET:: Gridin
    TYPE(OrbitInfo) ,TARGET:: Orbitin
    TYPE(PotentialInfo) ,TARGET:: Potin
    TYPE(FCInfo) ,TARGET:: FCin
    TYPE(SCFInfo) ,TARGET:: SCFin

    TYPE(OrbitInfo) :: tmpOrbit
    INTEGER :: i,n,io,icount,ok,loop,l,last
    REAL(8) :: en1,x,tol,err
    REAL(8), ALLOCATABLE :: arg(:),starg(:)
    REAL(8) :: errv0=1.d-10,errv,mixHF=0.4d0,mixX=0.2d0
    LOGICAL :: success,done
    INTEGER :: nrestart=8,mxloop1=25,mxloop2=50,mxloop=1000
    INTEGER :: firsttime=0
    LOGICAL :: noalt=.true.

    tol=errv0
    Gridwk=>Gridin
    Orbitwk=>Orbitin
    Potwk=>Potin
    FCwk=>FCin
    SCFwk=>SCFin

    verboseoutput=lotsofoutput

    n=Gridwk%n
    ALLOCATE(arg(n),starg(n))

    IF (firsttime<1) THEN
       io=Orbitwk%norbit
       ALLOCATE(HSZ%psi(n,io),HSZ%shift(n),HSZ%grad(n),HSZ%U(io),&
&           HSZ%Uref(io),HSZ%Ucore(io),HSZ%Uvale(io),HSZ%coreshift(n),&
&           HSZ%rVxKLI(n),HSZ%rDVxKLI(n),&
&           HSZ%LMBD(io,io), HSZ%psiref(n,io),HSZ%rVxref(n),&
&           HSZ%LMBDref(io,io),HSZ%LMBDcore(io,io),HSZ%LMBDvale(io,io),stat=ok)
       IF (ok /=0 ) THEN
          WRITE(6,*) 'Error in allocate HSZdata ', io,n,ok
          STOP
       ENDIF
       HSZ%psi=0
       HSZ%shift=0
       HSZ%grad=0
       HSZ%U=0; HSZ%Ucore=0;HSZ%Uvale=0
       HSZ%LMBD=0;HSZ%LMBDcore=0;HSZ%LMBDvale=0;HSZ%LMBDref=0
       HSZ%grad2=1.d13
       HSZ%rVxref =0
       HSZ%coreshift=0
       HSZ%rVxKLI=0
       HSZ%rDVxKLI=0
       HSZ%matchpoint=Gridwk%n

       !CALL Init_EXX_vx(Gridwk,Orbitwk,Potwk)

       !Potwk%rv=Potwk%rvh+Potwk%rvx
       !CALL zeropot(Gridwk,Potwk%rv,Potwk%v0,Potwk%v0p)
       !Potwk%rv=Potwk%rv+Potwk%rvn

       !arg=Potwk%rv

       !CALL InitAnderson_dr(AC,6,5,n,0.5d0,1.d3,mxloop,1.d-11,1.d-16,.TRUE.)
       !CALL DoAndersonMix(AC,arg,en1,EXX1sub,success)
       !SCFwk%iter=AC%CurIter
       !SCFwk%delta=AC%res
       !CALL FreeAnderson(AC)
       !WRITE(6,*) 'Finished first step ', en1 ,' success = ', success
       !WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta


       ! Initialize using HF
       n=Gridwk%n;io=Orbitwk%norbit
       ALLOCATE(HF%lmbd(io,io),HF%SumY(n,io),HF%lmany(io), &
&           HF%lmap(io,io), HF%emin(io),HF%emax(io),   &
&           HF%CSlmany(io),HF%CSlmap(io,io),stat=ok)
       IF (ok /=0 ) THEN
          WRITE(6,*) 'Error in allocate HFdata ', io,ok
          STOP
       ENDIF
       HF%lmbd=0

       HF%lmax=0
       DO io=1,Orbitwk%norbit
          IF (Orbitwk%l(io)>HF%lmax) &
&              HF%lmax=Orbitwk%l(io)
          HF%emin(io)=-Potwk%nz**2/Orbitwk%np(io)**2
          HF%emax(io)=-1.d-5
       ENDDO

       HF%lmap=0
       DO l=0,HF%lmax
          HF%lmany(l+1)=0
          DO io=1,Orbitwk%norbit
             IF (Orbitwk%l(io)==l) THEN
                HF%lmany(l+1)=HF%lmany(l+1)+1
                HF%lmap(l+1,HF%lmany(l+1))=io
             ENDIF
          ENDDO
       ENDDO

       HF%CSlmap=0
       DO l=0,HF%lmax
          HF%CSlmany(l+1)=0
          x=2*(2*l+1)
          DO io=1,Orbitwk%norbit
             IF (Orbitwk%l(io)==l.AND.ABS(x-Orbitwk%occ(io))<1.d-5) THEN
                HF%CSlmany(l+1)=HF%CSlmany(l+1)+1
                HF%CSlmap(l+1,HF%CSlmany(l+1))=io
             ENDIF
          ENDDO
       ENDDO

       CALL HF_tmp(Gridwk,Orbitwk,Potwk,FCwk,SCFwk)

       errv=1.d10;
       DO loop=1,mxloop2
          WRITE(6,*) 'Init HF ----loop---- ',loop
          CALL  HFIter(mixHF,errv0,errv,success,'ALL')
          IF (success) THEN
             WRITE(6,*) ' wfn iteration converged ', loop
             EXIT
          ENDIF
       ENDDO

       Call KLIVX(Gridwk,Orbitwk,HSZ,Potwk%rvx)
       Potwk%rv=Potwk%rvn+Potwk%rvh+Potwk%rvx
       call Report_EXX_functions('HF0')

       ! estimate Potwk%rvx
       !CALL SetIndex(Orbitwk)
       !last=HSZ%zero_index
       !CALL InitialVx(Gridwk,Orbitwk%wfn(:,last),arg)
       !call NEWVX(Gridwk,Orbitwk,HSZ,arg,Potwk%rvx,done)

       !OPEN(unit=1001,file='VX000',form='formatted')
       !DO i = 1,n
       !    WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), arg(i),Potwk%rvx(i)
       !ENDDO
       !CLOSE(1001)

       firsttime=1
       !CALL Report_EXX_functions('FI')

    ENDIF

    IF (setupfrozencore) THEN
       !HSZ%rVxref=Potwk%rvx    ! should have already been done
       CALL Get_KinCoul(Gridwk,Potwk,Orbitwk,SCFwk,noalt)
       write(6,*) 'setupfrozencore', SCFwk%ekin
       CALL Get_FCKinCoul(Gridwk,Potwk,Orbitwk,FCwk,SCFwk,noalt)
       write(6,*) 'setupfrozencore', SCFwk%ekin
       CALL Get_FCEnergy_EXX(Gridwk,Orbitwk,FCwk,SCFwk)
       CALL Report_EXX_functions('SC')
       CALL Total_FCEnergy_Report(SCFwk,6)
       RETURN
    ENDIF


    arg=Potwk%rvx
    if (TRIM(exctype)=='EXX') then
     write(6,*) 'EXX option is currently broken!!'; call flush(6); stop
     IF (frozencorecalculation.and..not.setupfrozencore) THEN
        WRITE(6,*) 'Beginning FC'
        CALL Init_EXX_vx(Gridwk,Orbitwk,Potwk)
        arg=Potwk%rvx
     ELSE
        WRITE(6,*) 'In  EXX_SCF starting second loop'; CALL flush(6)
        CALL Init_EXX_vx(Gridwk,Orbitwk,Potwk)
        arg=Potwk%rvx
        WRITE(6,*) 'In  EXX_SCF after Init_EXX_vx'
        !CALL Report_EXX_functions('SE')
     ENDIF
    endif

    write(6,*) 'Before KLI loop'
    Do loop=1,mxloop
       !Call EXXKLIiter_innerloop(arg,en1,starg,SCFwk%delta,success,.true.)
       !Call EXXKLIiter(arg,en1,starg,SCFwk%delta,success,.true.)
       Call EXXKLIiter_wfn(loop,mixX,tol,err,success)
       if (success) exit
       success=.false.
       SCFwk%iter=loop
       SCFwk%delta=err
    enddo
       Potwk%rv=Potwk%rvn+Potwk%rvh+Potwk%rvx
       WRITE(6,*) 'Finished KLI iter ', loop,' success = ', success
       WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta
       HSZ%rVxKLI(:)=Potwk%rvx; HSZ%rDVxKLI(:)=0.d0
       if (frozencorecalculation) then
           call Report_EXX_functions('KF')
       else
           call Report_EXX_functions('KA')
       endif
       if(TRIM(exctype)=='EXXKLI') RETURN


    DEALLOCATE(arg,starg)
  END SUBROUTINE EXX_SCF

  SUBROUTINE EXXOCC_SCF(scftype,lotsofoutput,Gridin,Orbitin,Potin,FCin,SCFin)
    CHARACTER(2), INTENT(IN) :: scftype
    LOGICAL, INTENT(IN) :: lotsofoutput
    TYPE(GridInfo) ,TARGET:: Gridin
    TYPE(OrbitInfo) ,TARGET:: Orbitin
    TYPE(PotentialInfo) ,TARGET:: Potin
    TYPE(FCInfo) ,TARGET:: FCin
    TYPE(SCFInfo) ,TARGET:: SCFin

    TYPE(OrbitInfo) :: tmpOrbit
    INTEGER :: i,n,io,icount,ok,loop,l,last
    REAL(8) :: en1,x,tol,err
    REAL(8), ALLOCATABLE :: arg(:),starg(:)
    REAL(8) :: errv0=1.d-10,errv,mixHF=0.4d0,mixX=0.2d0
    LOGICAL :: success,done
    INTEGER :: nrestart=8,mxloop1=25,mxloop2=50,mxloop=1000
    INTEGER :: firsttime=0
    LOGICAL :: noalt=.true.

    tol=errv0
    Gridwk=>Gridin
    Orbitwk=>Orbitin
    Potwk=>Potin
    FCwk=>FCin
    SCFwk=>SCFin

    verboseoutput=lotsofoutput

    n=Gridwk%n
    ALLOCATE(arg(n),starg(n))

    IF (firsttime<1) THEN
       io=Orbitwk%norbit
       ALLOCATE(HSZ%psi(n,io),HSZ%shift(n),HSZ%grad(n),HSZ%U(io),&
&           HSZ%Uref(io),HSZ%Ucore(io),HSZ%Uvale(io),HSZ%coreshift(n),&
&           HSZ%rVxKLI(n),HSZ%rDVxKLI(n),&
&           HSZ%LMBD(io,io), HSZ%psiref(n,io),HSZ%rVxref(n),&
&           HSZ%LMBDref(io,io),HSZ%LMBDcore(io,io),HSZ%LMBDvale(io,io),stat=ok)
       IF (ok /=0 ) THEN
          WRITE(6,*) 'Error in allocate HSZdata ', io,n,ok
          STOP
       ENDIF
       HSZ%psi=0
       HSZ%shift=0
       HSZ%grad=0
       HSZ%U=0; HSZ%Ucore=0;HSZ%Uvale=0
       HSZ%LMBD=0;HSZ%LMBDcore=0;HSZ%LMBDvale=0;HSZ%LMBDref=0
       HSZ%grad2=1.d13
       HSZ%rVxref =0
       HSZ%coreshift=0
       HSZ%rVxKLI=0
       HSZ%rDVxKLI=0
       HSZ%matchpoint=Gridwk%n


       ! Initialize using HF
       n=Gridwk%n;io=Orbitwk%norbit
       ALLOCATE(HF%lmbd(io,io),HF%SumY(n,io),HF%lmany(io), &
&           HF%lmap(io,io), HF%emin(io),HF%emax(io),   &
&           HF%CSlmany(io), HF%CSlmap(io,io),stat=ok)
       IF (ok /=0 ) THEN
          WRITE(6,*) 'Error in allocate HFdata ', io,ok
          STOP
       ENDIF
       HF%lmbd=0

       HF%lmax=0
       DO io=1,Orbitwk%norbit
          IF (Orbitwk%l(io)>HF%lmax) &
&              HF%lmax=Orbitwk%l(io)
          HF%emin(io)=-Potwk%nz**2/Orbitwk%np(io)**2
          HF%emax(io)=-1.d-5
       ENDDO

       HF%lmap=0
       DO l=0,HF%lmax
          HF%lmany(l+1)=0
          DO io=1,Orbitwk%norbit
             IF (Orbitwk%l(io)==l) THEN
                HF%lmany(l+1)=HF%lmany(l+1)+1
                HF%lmap(l+1,HF%lmany(l+1))=io
             ENDIF
          ENDDO
       ENDDO

       HF%CSlmap=0
       DO l=0,HF%lmax
          HF%CSlmany(l+1)=0
          x=2*(2*l+1)
          DO io=1,Orbitwk%norbit
             IF (Orbitwk%l(io)==l.AND.ABS(x-Orbitwk%occ(io))<1.d-5) THEN
                HF%CSlmany(l+1)=HF%CSlmany(l+1)+1
                HF%CSlmap(l+1,HF%CSlmany(l+1))=io
             ENDIF
          ENDDO
       ENDDO

       CALL HF_tmp(Gridwk,Orbitwk,Potwk,FCwk,SCFwk)

       errv=1.d10;
       DO loop=1,mxloop2
          WRITE(6,*) 'Init HF ----loop---- ',loop
          CALL  HFIter(mixHF,errv0,errv,success,'ALL')
          IF (success) THEN
             WRITE(6,*) ' wfn iteration converged ', loop
             EXIT
          ENDIF
       ENDDO

       Call VXOCC(Gridwk,Orbitwk,HSZ,Potwk%rvx)
       Potwk%rv=Potwk%rvn+Potwk%rvh+Potwk%rvx
       call Report_EXX_functions('HF0')


       firsttime=1
       !CALL Report_EXX_functions('FI')

    ENDIF

    IF (setupfrozencore) THEN
       !HSZ%rVxref=Potwk%rvx    ! should have already been done
       CALL Get_KinCoul(Gridwk,Potwk,Orbitwk,SCFwk,noalt)
       write(6,*) 'setupfrozencore', SCFwk%ekin
       CALL Get_FCKinCoul(Gridwk,Potwk,Orbitwk,FCwk,SCFwk,noalt)
       write(6,*) 'setupfrozencore', SCFwk%ekin
       CALL Get_FCEnergy_EXX(Gridwk,Orbitwk,FCwk,SCFwk)
       CALL Report_EXX_functions('SC')
       CALL Total_FCEnergy_Report(SCFwk,6)
       RETURN
    ENDIF


    arg=Potwk%rvx
     IF (frozencorecalculation.and..not.setupfrozencore) THEN
        WRITE(6,*) 'Beginning FC'
        CALL Init_EXX_vx(Gridwk,Orbitwk,Potwk)
        arg=Potwk%rvx
     ELSE
        WRITE(6,*) 'In  EXX_SCF starting second loop'; CALL flush(6)
        CALL Init_EXX_vx(Gridwk,Orbitwk,Potwk)
        arg=Potwk%rvx
        WRITE(6,*) 'In  EXX_SCF after Init_EXX_vx'
        !CALL Report_EXX_functions('SE')
     ENDIF

    write(6,*) 'Before OCC loop'
    Do loop=1,mxloop
       Call EXXOCCiter_wfn(loop,mixX,tol,err,success)
       if (success) exit
       success=.false.
       SCFwk%iter=loop
       SCFwk%delta=err
    enddo
       Potwk%rv=Potwk%rvn+Potwk%rvh+Potwk%rvx
       WRITE(6,*) 'Finished OCC iter ', loop,' success = ', success
       WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta
       HSZ%rVxKLI(:)=Potwk%rvx; HSZ%rDVxKLI(:)=0.d0
       if (frozencorecalculation) then
           call Report_EXX_functions('OF')
       else
           call Report_EXX_functions('OA')
       endif
       if(TRIM(exctype)=='EXXOCC') RETURN


    DEALLOCATE(arg,starg)
  END SUBROUTINE EXXOCC_SCF

  SUBROUTINE EXX_Input_Settings(inputline)
    CHARACTER(128), INTENT(IN) :: inputline

    INTEGER :: i
    HSZ%Fixed_Zero=.FALSE.
    i=INDEX(inputline,'FIXED_ZERO')
    IF (i>0) THEN
       HSZ%Fixed_Zero=.TRUE.
       READ(unit=inputline(i+10:128),fmt=*) HSZ%zero_index
    ENDIF
  END SUBROUTINE EXX_Input_Settings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Init_EXX_vx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Init_EXX_vx(Grid,Orbit,Pot)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot

    REAL(8), POINTER :: rvx(:)
    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: i,j,k,l,n,outer,ok,io
    REAL(8) :: q
    REAL(8), PARAMETER :: threshold=1.d-8

    WRITE(6,*) ' In Init_EXX_vx ' ; CALL flush(6)
    ! determine Fock potential for outer most wavefunction
    CALL SetIndex(Orbit)
    outer=HSZ%zero_index

    WRITE(6,*) 'zero_index ', outer ; CALL flush(6)
    n=Grid%n
    rvx=>Pot%rvx

    ALLOCATE(dum(n), stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error in init_exx_vx ',i,n
       STOP
    ENDIF
    CALL InitialVx(Grid,Orbit%wfn(:,outer),dum)

    WRITE(6,*) ' After InitialVx ' ; CALL flush(6)

    rvx=dum

    DEALLOCATE(dum)
  END SUBROUTINE Init_EXX_vx

  SUBROUTINE InitialVx(Grid,wfn,vl)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: wfn(:)
    REAL(8), INTENT(INOUT) :: vl(:)

    INTEGER :: i,j,k,l,m,n
    REAL(8), ALLOCATABLE :: dum(:)
    REAL(8) :: ecoul,q,v00

    n=Grid%n
    ALLOCATE(dum(n))

    dum=wfn(:)**2
    CALL poisson(Grid,q,dum,vl,ecoul,v00)
    !WRITE(6,*) 'Completed Poisson ', q,ecoul,v00; CALL flush(6)
    vl=-vl

    DEALLOCATE(dum)

  END SUBROUTINE InitialVx

  SUBROUTINE SetIndex(Orbit)
    TYPE(orbitinfo), INTENT(IN) :: Orbit

    INTEGER :: i,j,k,l,n,outer,ok,io
    REAL(8) :: x


    IF (.NOT.HSZ%FIXED_ZERO) THEN
       outer=1;x=Orbit%eig(outer)
       IF (Orbit%norbit>1) THEN
          DO io=2, Orbit%norbit
             IF (Orbit%eig(io)<0.d0.AND.Orbit%occ(io)>1.d-5  &
&                 .AND.Orbit%eig(io)>x) THEN
                outer=io;x=Orbit%eig(outer)
             ENDIF
          ENDDO
       ENDIF

       HSZ%zero_index=outer
    ELSE
       outer=HSZ%zero_index
    ENDIF

    HSZ%lmax=0

    DO io=1,Orbit%norbit
       IF (HSZ%lmax<Orbit%l(io))  HSZ%lmax=Orbit%l(io)
    ENDDO

    WRITE(6,*) 'Returning from SetIndex ', HSZ%zero_index,HSZ%lmax
    IF (frozencorecalculation) THEN
       DO io=1,Orbit%norbit
          IF (Orbit%iscore(io).AND.outer==io) THEN
             WRITE(6,*) 'Error in FC SetIndex -- core state ', io,outer
             STOP
          ENDIF
       ENDDO
    ENDIF


  END SUBROUTINE SetIndex


  SUBROUTINE SetIndexKLI(Orbit)
    TYPE(orbitinfo), INTENT(IN) :: Orbit

    INTEGER :: i,j,k,l,n,outer,ok,io
    REAL(8) :: x,y

    k=FindGridIndex(Gridwk,10.d0)
    IF (.NOT.HSZ%FIXED_ZERO) THEN
       outer=1;x=ABS(Orbit%wfn(k,outer))
       IF (Orbit%norbit>1) THEN
          DO io=2, Orbit%norbit
             y=ABS(Orbit%wfn(k,io))
             IF (y>x.AND.Orbit%occ(io)>1.d-5)  THEN
                outer=io;x=y
             ENDIF
          ENDDO
       ENDIF

       HSZ%zero_index=outer
    ELSE
       outer=HSZ%zero_index
    ENDIF

    WRITE(6,*) 'Returning from SetIndexKLI ', HSZ%zero_index
    IF (frozencorecalculation) THEN
       DO io=1,Orbit%norbit
          IF (Orbit%iscore(io).AND.outer==io) THEN
             WRITE(6,*) 'Error in FC SetIndex -- core state ', io,outer
             STOP
          ENDIF
       ENDDO
    ENDIF


  END SUBROUTINE SetIndexKLI

  SUBROUTINE Setfparms(Grid,Orbit)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit

    REAL(8) :: betaL,q

    CALL SetIndex(Orbit)

    IF (HSZ%zero_index>0) THEN
       betaL=SQRT(ABS(Orbit%eig(HSZ%zero_index)))
    ELSE
       WRITE(6,*) 'Error in Setfparms ', HSZ%zero_index
       STOP
    ENDIF

    HSZ%betaL=ABS(betaL)

    WRITE(6,*) 'Returning from setfparams with ' , HSZ%betaL

  END SUBROUTINE Setfparms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! EXX1sub  -- simple vx version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXX1sub(w,energy,residue,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: residue(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io
    REAL(8),ALLOCATABLE :: dum(:)
    INTEGER:: fcount=0
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot

    ALLOCATE(dum(Gridwk%n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in EXXsub allocation' ,Gridwk%n
       STOP
    ENDIF

    n=Gridwk%n

    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,w,success)
    tmpPot%rv=w

    !if FC core calc , restore core info backinto tmpOrbit

    IF(frozencorecalculation) THEN
       WRITE(6,*) ' Frozencore calculation should not use this routine!'
       STOP
       DO io = 1 , Orbitwk%norbit
          IF(Orbitwk%iscore(io)) THEN
             tmpOrbit%eig(io)=Orbitwk%eig(io)
             tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
          ENDIF
       ENDDO
    ENDIF

    IF (.NOT.success) THEN
       WRITE(6,*) 'Bad luck in Sub'
    ENDIF

    CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)
    CALL Get_Energy_EXX(Gridwk,tmpOrbit,SCFwk%eexc)
    CALL SetIndex(tmpOrbit) ! find index of most extended wfn


    SCFwk%etot=SCFwk%ekin+SCFwk%estatic+SCFwk%eexc

    IF(ColleSalvetti) THEN
       CALL Get_Energy_CS(Gridwk,tmpOrbit,SCFwk%oepcs)
       WRITE(6,*) 'Colle :', SCFwk%oepcs
       SCFwk%etot=SCFwk%etot+SCFwk%oepcs
    ENDIF

    energy=SCFwk%etot
    CALL Total_Energy_Report(SCFwk,6)
    CALL Init_EXX_vx(Gridwk,tmpOrbit,tmpPot)
    dum=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx

    residue=dum-w
    err=DOT_PRODUCT(residue,residue)

    IF (update) THEN
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

       CALL One_electron_energy_Report(Orbitwk,6)
    ENDIF

    CALL DestroyPot(tmpPot)
    CALL DestroyOrbit(tmpOrbit)
    DEALLOCATE (dum)

  END SUBROUTINE  EXX1sub

  SUBROUTINE ApproxVx(Grid,Orbit,last,vl)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    INTEGER, INTENT(IN) :: last
    REAL(8), INTENT(INOUT) :: vl(:)

    INTEGER :: i,j,k,l,m,n
    REAL(8),ALLOCATABLE :: dum(:),vl_cs(:)

    n=Grid%n
    ALLOCATE(dum(n),vl_cs(n))
    CALL InitialVx(Grid,Orbit%wfn(:,last),dum)
    IF (frozencorecalculation) THEN
       CALL Calc_dexdphi_io_v(Gridwk,Orbit,last,vl)
    ELSE
       CALL Calc_dexdphi_io(Gridwk,Orbit,last,vl)
       IF(ColleSalvetti) THEN
          CALL Calc_decdchi_io(Orbitwk,last,vl_cs)
          vl = vl + vl_cs
       ENDIF

    ENDIF
    CALL SetMatchPoint(Grid,Orbit,last)
    m=HSZ%matchpoint
    j=m       ! grid point near last maximum
    !write(6,*) 'ApproxVx: max ', m,Grid%r(m)
    DO i=1,n
       IF (i<j) THEN
          vl(i)=vl(i)/Orbit%wfn(j,last)
       ELSE IF (ABS(Orbit%wfn(i,last))>1.d-36) THEN
          vl(i)=vl(i)/Orbit%wfn(i,last)
       ELSE
          vl(i)=-2.d0
       ENDIF
       !write(6,'(1p,15E15.7)') Grid%r(i),vl(i),Orbit%wfn(i,last)
    ENDDO
    j=m+1
    DO i=m+1,n
       j=i
       IF (ABS(vl(i)-dum(i))<1.d-8) EXIT
    ENDDO
    vl(j:n)=dum(j:n)

    DEALLOCATE(dum,vl_cs)

  END SUBROUTINE ApproxVx


  SUBROUTINE SetMatchPoint(Grid,Orbit,last)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    INTEGER, INTENT(IN) :: last

    INTEGER :: i,j,k,l,m,n

    n=Grid%n
    j=2   ! Find last zero crossing
    DO i=3,n
       IF (Orbit%wfn(i,last)*Orbit%wfn(i-1,last)<0.d0) j=i
    ENDDO

    IF (j>n-1) j=n/10+1
    m=j+1
    DO i=j+1,n   ! Find last maximum
       IF (ABS(Orbit%wfn(i,last))<ABS(Orbit%wfn(i-1,last)) ) THEN
          m=i-1
          EXIT
       ENDIF
    ENDDO
    j=m       ! grid point near last maximum

    HSZ%matchpoint=m+1
    WRITE(6,*) 'SetMatchPoint matchpoint ', m+1,Grid%r(HSZ%matchpoint)
  END SUBROUTINE SetMatchPoint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Calc_dedv
!!!     Accumulate solutions to inhomogeneous equations
!!!        and determine shift
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_dedv(Grid,Pot,Orbit,EigOrbit)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(PotentialInfo), INTENT(IN):: Pot
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    TYPE(OrbitInfo), TARGET, INTENT(IN):: EigOrbit

    REAL(8), POINTER :: phi(:,:),phr(:,:),psi(:,:),r(:),rv(:),U(:),LMBD(:,:)
    REAL(8), ALLOCATABLE :: vl(:),dum(:),rhs(:),shift(:),vs(:),rvx(:),vl_cs(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,io,jo,norbit,last
    INTEGER :: iter
    INTEGER :: counter=0
    LOGICAL :: ok
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,xx,q,ecoul,v00
    REAL(8), PARAMETER :: threshold=1.d-8,tol=1.d-5
    TYPE(Anderson_Context) :: A_DEDV
    REAL(8) :: ModVl,ModVl_CS,TempU

    EigOrbitwk=>EigOrbit
    n=Grid%n

    rv => Pot%rv

    n=Grid%n
    norbit=Orbit%norbit
    phi=>Orbit%wfn
    r=>Grid%r

    IF (counter==0) THEN
       ALLOCATE(geqn(n,norbit))
    ENDIF
    counter=counter+1

    geqn=0
    HSZ%shift=0
    psi=>HSZ%psi
    U=>HSZ%U
    LMBD=>HSZ%LMBD

    psi=0

    ALLOCATE(vl(n),dum(n),rhs(n),shift(n),vs(n),rvx(n),vl_cs(n))

    rvx=Pot%rvx
    last=HSZ%zero_index
    !write(6,*) 'last = ', last; call flush(6)

    ! main contributions
    ! Determine constant coefficients -- U, LMBA, X
    DO io=1,norbit
       occ=Orbit%occ(io)
       IF (occ>threshold) THEN
          IF(frozencorecalculation.AND.Orbit%iscore(io)) THEN
          ELSE
             IF(frozencorecalculation) THEN
                CALL Calc_dexdphi_io_v(Grid,Orbit,io,vl)
             ELSE
                CALL Calc_dexdphi_io(Grid,Orbit,io,vl)
                IF(ColleSalvetti) THEN
                   CALL Calc_decdchi_io(Orbit,io,vl_cs)
                   vl=vl+vl_cs
                ENDIF
             ENDIF
             rhs=0
             rhs(2:n)=(-rvx(2:n)*phi(2:n,io)+vl(2:n))/r(2:n)
             geqn(:,io)=rhs
             vl=rhs
             dum=vl*phi(:,io)
             U(io)=integrator(Gridwk,dum)
             !DO jo=1,norbit
             !   IF (jo/=io) THEN
             !      IF (Orbit%l(jo)==Orbit%l(io)) THEN
             !         dum=vl*EigOrbit%wfn(:,jo)
             !         LMBD(io,jo)=integrator(Gridwk,dum)
             !      ENDIF
             !   ENDIF
             !ENDDO
          ENDIF
       ENDIF
    ENDDO

    shift=0
    ! Determine the g^0 functions
    DO io=1,norbit
       occ=Orbit%occ(io)
       IF (occ>threshold) THEN
          IF(frozencorecalculation.AND.Orbit%iscore(io)) THEN
          ELSE
             rhs=0
             rhs=rhs+geqn(:,io)
             vl=rhs
             dum=vl*phi(:,io)
             term=integrator(Grid,dum)
             rhs=rhs-(term)*phi(:,io)
             !DO jo=1,norbit   !   extra orbit orthogonalization
             !   IF (jo/=io) THEN
             !     IF (Orbitwk%l(jo)==Orbitwk%l(io)) THEN
             !        dum=vl*EigOrbitwk%wfn(:,jo)
             !        term=integrator(Grid,dum)
             !        rhs=rhs-term*EigOrbitwk%wfn(:,jo)
             !     ENDIF
             !  ENDIF
             !ENDDO
             geqn(:,io)=rhs    ! store for possible later use
             CALL Calc_psi(Grid,&
&                 Orbitwk%l(io),Orbitwk%eig(io),rv,phi(:,io),rhs,dum)
             shift=shift+2*occ*phi(:,io)*dum
             psi(:,io)=dum
          ENDIF
       ENDIF
    ENDDO
    geqn(:,last)=geqn(:,last)+U(last)*phi(:,last)    ! ensuring asymptotic form

    !   check orthogonality   (should check for l values)
    !    do io=1,norbit
    !      do jo=1,norbit
    !         write(6,'("Check Ortho ", 2i5,1p,1e15.7)') &
&   !&                 io,jo,overlap(Grid,EigOrbitwk%wfn(:,jo),psi(:,io))
    !      enddo
    !   enddo

    HSZ%shift=shift
    if (setupfrozencore) then
      do io=1,norbit
         if (.not.Orbitwk%iscore(io)) then
            HSZ%Ucore(io)=HSZ%Uref(io)-HSZ%Uvale(io)
            !do jo=1,norbit
            !   if (jo/=io) &
&           !&   HSZ%LMBDcore(io,jo)=HSZ%LMBDref(io,jo)-HSZ%LMBDvale(io,jo)
            !enddo
         endif
      enddo
    else if (frozencorecalculation) then
      do io=1,norbit
         if (.not.Orbitwk%iscore(io)) then
            HSZ%U(io)=HSZ%Ucore(io)+HSZ%Uvale(io)
         endif
      enddo
    endif

    write(6,*) 'Returning from Calc_dedv with Uzero = ', U(HSZ%zero_index)

    DEALLOCATE(vl,dum,rhs,shift,vs,rvx,vl_cs)
    RETURN

    !!!!!!!!!!!no longer used!!!!!!!!!!!
    !IF (MAXVAL(ABS(vs))>0.001d0) THEN
    !   DEALLOCATE(vl,dum,rhs,shift,vs,rvx,vl_cs)
    !   RETURN
    !ENDIF
    !
    !CALL InitAnderson_dr(A_DEDV,6,5,n,0.005d0,1.d3,100,1.d-8,1.d-16,.FALSE.)
    !CALL DoAndersonMix(A_DEDV,vs,xx,vssub,OK)
    !WRITE(6,*) 'Leaving dedv with ', A_DEDV%CurIter, A_DEDV%res,MAXVAL(ABS(vs))
    !CALL FreeAnderson(A_DEDV)
    !
    !DEALLOCATE(vl,dum,rhs,shift,vs,rvx,vl_cs)
  END SUBROUTINE Calc_dedv

  SUBROUTINE Set_Vxref(Pot)
    TYPE(PotentialInfo), INTENT(IN):: Pot

    HSZ%rVxref=Pot%rvx
    HSZ%psiref=HSZ%psi
    HSZ%Uref=HSZ%U   ; HSZ%Ucore=HSZ%Uref
    HSZ%LMBDref=HSZ%LMBD    ;    HSZ%LMBDcore=HSZ%LMBDref

  END SUBROUTINE Set_Vxref

  SUBROUTINE writestuff
    INTEGER :: i,j,k,n,nmap
    INTEGER, SAVE :: counter=0
    CHARACTER(4) :: stuff
    INTEGER, ALLOCATABLE :: mmap(:)

    CALL mkname(counter,stuff)
    OPEN (unit=1001,file='shift'//TRIM(stuff),form='formatted')
    n=Gridwk%n
    IF (frozencorecalculation) THEN
       ALLOCATE(mmap(Orbitwk%norbit))
       j=0
       DO i=1,Orbitwk%norbit
          IF (.NOT.Orbitwk%iscore(i)) THEN
             j=j+1
             mmap(j)=i
          ENDIF
       ENDDO
       nmap=j
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),HSZ%shift(i),&
&              (HSZ%psi(i,mmap(k))*Orbitwk%wfn(i,mmap(k)),k=1,nmap)
       ENDDO
    ELSE
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),HSZ%shift(i) ,&
&              (HSZ%psi(i,j)*Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
       ENDDO
    ENDIF
    CLOSE(1001)
    OPEN (unit=1001,file='pot'//TRIM(stuff),form='formatted')
    DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&              Potwk%rvh(i),Potwk%rvx(i)
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='wfn'//TRIM(stuff),form='formatted')
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&           (Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='psi'//TRIM(stuff),form='formatted')
    IF (frozencorecalculation) THEN
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&              (HSZ%psi(i,mmap(k)),k=1,nmap)
       ENDDO
    ELSE
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&              (HSZ%psi(i,j),j=1,Orbitwk%norbit)
       ENDDO
    ENDIF
    CLOSE(1001)
    counter=counter+1
    IF (ALLOCATED(mmap)) DEALLOCATE(mmap)
  END SUBROUTINE writestuff


  SUBROUTINE Report_EXX_functions(sub)
    CHARACTER(2) :: sub
    INTEGER :: i,j,k,n,nmap,io
    INTEGER, ALLOCATABLE :: mmap(:)

    INTEGER, SAVE :: counter=0
    CHARACTER(4) :: stuff

    WRITE(6,*)
    WRITE(6,*) ' Summary of EXX orbitals '
    WRITE(6,"(' n  l     occupancy            energy              U')")
    DO io=1,Orbitwk%norbit
       WRITE(6,'(i2,1x,i2,4x,1p,3e15.7)') &
&           Orbitwk%np(io),Orbitwk%l(io),Orbitwk%occ(io),Orbitwk%eig(io),HSZ%U(io)
    ENDDO

    CALL mkname(counter,stuff)

    OPEN (unit=1001,file='shift'//sub//TRIM(stuff),form='formatted')
    n=Gridwk%n
    IF (frozencorecalculation) THEN
       ALLOCATE(mmap(Orbitwk%norbit))
       j=0
       DO i=1,Orbitwk%norbit
          IF (.NOT.Orbitwk%iscore(i)) THEN
             j=j+1
             mmap(j)=i
          ENDIF
       ENDDO
       nmap=j
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),HSZ%shift(i) ,&
&              (HSZ%psi(i,mmap(k))*Orbitwk%wfn(i,mmap(k)),k=1,nmap)
       ENDDO
    ELSE
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),HSZ%shift(i) ,&
&              (HSZ%psi(i,j)*Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
       ENDDO
    ENDIF
    CLOSE(1001)
    OPEN (unit=1001,file='pot'//sub//TRIM(stuff),form='formatted')
    DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&              Potwk%rvh(i),Potwk%rvx(i)
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), Orbitwk%den(i),&
&           (Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='psi'//sub//TRIM(stuff),form='formatted')
    IF (frozencorecalculation) THEN
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&              (HSZ%psi(i,mmap(k)),k=1,nmap)
       ENDDO
    ELSE
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&              (HSZ%psi(i,j),j=1,Orbitwk%norbit)
       ENDDO
    ENDIF
    CLOSE(1001)

    IF (ALLOCATED(mmap)) DEALLOCATE(mmap)

    counter=counter+1

  END SUBROUTINE Report_EXX_functions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Get_FCEnergy_EXX                    !!!!
  !!   Valence part of the Fock exchange energy
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_FCEnergy_EXX(Grid,Orbit,FC,SCF)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    TYPE(FCInfo), INTENT(IN) :: FC
    TYPE(SCFInfo), INTENT(INOUT) :: SCF

    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: dum(:)
    REAL(8) :: FCeex, x
    integer :: i,n

    n=Grid%n
    r=>grid%r

    ALLOCATE(dum(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'allocation error in Get_FCEnergy_EXX', i,n
       STOP
    ENDIF

    call Get_Energy_EXX_VC(Grid,Orbit,FCeex)
    write(6,*) 'Frozencore core exchange contribution ', FCeex
    call Get_Energy_EXX_VV(Grid,Orbit,x)
    write(6,*) 'Frozencore valence exchange contribution ', x

    SCF%valeexc=FCeex+x
    SCF%evale=SCF%valekin+SCF%valecoul+SCF%valeexc

   deallocate(dum)

  END SUBROUTINE Get_FCEnergy_EXX

!!!!!!!!!!!!!!!!!!!!!Colle Salvetti routines added by Xiao Xu 10-06-2008


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Colle Salvetti correlation energy
! Still using analytical wavefun and density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Get_Energy_CS(Grid,Orbit,csenergy)
    IMPLICIT NONE
    TYPE(gridinfo),INTENT(IN) ::Grid
    TYPE(orbitinfo),INTENT(INOUT)::Orbit
    REAL(8),INTENT(OUT)::csenergy


    INTEGER :: i,n,norbit,io
    REAL(8),ALLOCATABLE :: RealWfn(:,:),RealDen(:),eta(:),xi(:)
    REAL(8),ALLOCATABLE :: DerivRealWfn(:,:),DerivRealDen(:),&
&        Deriv2RealDen(:),LaplaceRealDen(:)
    REAL(8),ALLOCATABLE :: kernel(:),dum(:)
    REAL(8) :: a ,b ,c ,d ,factor
    REAL(8) :: ECS_1ST , ECS_2ND,ECS_3RD,ECS_4TH,ECS_5TH,ECS
    REAL(8) :: alpha

    a=0.04918d0
    a=a*2
    b=0.132d0
    c=0.2533d0
    d=0.349d0

    n=Grid%n;
    norbit = Orbit%norbit

    !Change norbit for testing
    !norbit = 1
    !alpha =0.1d0
    !DO i=2,n
    !    Orbit%wfn(i,1)=Grid%r(i)*exp(-alpha*Grid%r(i))
    !ENDDO
    !CALL extrapolate(Grid,Orbit%wfn(:,1))

    !Orbit%den(:)=(Orbit%wfn(:,1)**2) !!/(4.0d0*pi*Grid%r(:)**2)

    ALLOCATE(RealWfn(n,norbit),RealDen(n),eta(n),xi(n))
    ALLOCATE(DerivRealWfn(n,norbit),DerivRealDen(n),&
&        Deriv2RealDen(n),LaplaceRealDen(n))
    ALLOCATE(kernel(n))
    ALLOCATE(dum(n))


    !! restore analytical wavefunction and density
    DO io=1,norbit
       DO i=2,n
          RealWfn(i,io) = Orbit%wfn(i,io)/Grid%r(i)
       ENDDO
    ENDDO

    DO i=2,n
       RealDen(i)=Orbit%den(i)/(4.0d0*pi*(Grid%r(i)**2))
    ENDDO


    DO io=1,norbit
       CALL extrapolate(Grid,RealWfn(:,io))
    ENDDO
    CALL extrapolate(Grid,RealDen)

    !! Calc eta and xi
    eta(:)= 1.d0 + d/(RealDen(:)**(1.0d0/3.0d0))
    xi(:) =RealDen(:)**(5.0d0/3.0d0)*EXP(c*(RealDen(:)**(-1.0d0/3.0d0)))*eta(:)
    xi(:) = 1.d0/xi(:)


    ! To deal with eta diverge at long range
    DO i=1,n
       IF(xi(i).NE.xi(i)) THEN
          xi(i)=0
       ENDIF
    ENDDO

    !! derivative of wavefunction , density, and their laplacian
    DO io=1,norbit
       CALL derivative(Grid,RealWfn(:,io),DerivRealWfn(:,io),1,n)
    ENDDO
    CALL derivative(Grid,RealDen,DerivRealDen,1,n)
    CALL derivative(Grid,DerivRealDen,Deriv2RealDen,1,n)

    LaplaceRealDen(:)=2*DerivRealDen(:)/Grid%r(:) +  &
&        Deriv2RealDen(:) !! Blows up at origin
    CALL extrapolate(Grid,LaplaceRealDen)  !! deal with the infinity at origin



    !! Calculate Colle Salvetti Energy
    !! Term 1ECS
    kernel = 0.0d0
    DO io=1,norbit
       kernel(:)=kernel(:) + 2*Orbit%occ(io)*RealDen(:)*(DerivRealWfn(:,io)**2)
    ENDDO
    kernel(:) =- a*b*(0.25d0)*xi(:)*(Grid%r(:)**2)*kernel(:)
    ECS_1ST = integrator(Grid,kernel,1,n)


    !! Term 2
    kernel = 0.0d0
    DO io=1,norbit
       i=Orbit%l(io)
       factor=i*(i+1)

       IF(Orbit%l(io)==0) THEN
          kernel(:)= kernel(:) + 0.0d0
       ELSE
          kernel(:)=kernel(:) + &
&              2*Orbit%occ(io)*RealDen(:)*2*(RealWfn(:,io)**2)*&
&              factor/(Grid%r(:)**2)
       ENDIF
    ENDDO

    CALL extrapolate(Grid, kernel)
    kernel(:) =- a*b*(0.25d0)*xi(:)*(Grid%r(:)**2)*kernel(:)
    ECS_2ND = integrator(Grid,kernel,1,n)


    !! Term 3
    kernel = 0.0d0
    kernel(:) = a*b*xi(:)*(Grid%r(:)**2)*pi*(DerivRealDen(:)**2)
    ECS_3RD = integrator(Grid,kernel,1,n)


    !! Term 4
    kernel = 0.0d0
    kernel(:)=-a*b*xi(:)*(Grid%r(:)**2)*pi*(0.5d0)*RealDen(:)*LaplaceRealDen(:)
    ECS_4TH = integrator(Grid,kernel,1,n)


    !! Term5
    kernel = 0.0d0
    kernel(:) = - a*(Grid%r(:)**2)*(4.0d0)*pi*RealDen(:)/eta(:)
    ECS_5th = integrator(Grid,kernel,1,n)

    !! SUM UP
    ECS = ECS_1ST  + ECS_2ND + ECS_3RD + ECS_4TH + ECS_5TH

    csenergy = ECS

    WRITE(6,*) 'ECS is already Calculated' ,ECS

    DEALLOCATE(RealWfn,RealDen,eta,xi,DerivRealWfn,DerivRealDen,&
&        Deriv2RealDen,LaplaceRealDen,kernel,dum)


  END SUBROUTINE Get_Energy_CS


!!!!!!
!  Calc_decdphi_io
!!!
  SUBROUTINE Calc_decdchi_io(Orbit,io,res)
    TYPE(OrbitInfo),INTENT(IN)::Orbit
    INTEGER,INTENT(IN)::io
    REAL(8),INTENT(OUT) ::res(:)
    INTEGER :: i,n,norbit,IterOrbit
    REAL(8) :: a,b,c,d,ab
    REAL(8) ,ALLOCATABLE :: RealWfn(:,:),DerivRealWfn(:,:),Deriv2RealWfn(:,:)
    REAL(8) ,ALLOCATABLE :: RealDen(:),DerivRealDen(:),Deriv2RealDen(:),&
&        LaplaceRealDen(:)
    REAL(8) ,ALLOCATABLE :: eta(:),xi(:),dxidden(:),Derivxi(:),Deriv2xi(:),&
&        Derivxiden(:),Laplacexi(:)
    REAL(8) ,ALLOCATABLE :: Row1(:),Row2(:),Row3(:),Row4(:),Row5(:),&
&        Row6(:),Row7(:),decdchi(:),RowExtra(:)
    REAL(8) ,ALLOCATABLE :: dum(:),dum1(:),dum2(:)
    REAL(8) :: small=1.d-6
    REAL(8) :: Inte1,Inte2,Inte3,Inte4,Inte5,Inte6,Inte7,InteAll7
    INTEGER :: factor
    REAL(8) :: alpha
    REAL(8) :: XiaoTest


    ! This subroutine was written by Xiao Xu in Aug. 2008;
    !   It is not quite working yet
    WRITE(6,*) 'Colle Salvetti correlation function not yet quite working'
    WRITE(6,*) 'Program will stop.'
    STOP
    a=0.04918d0
    a=a*2
    b=0.132d0
    c=0.2533d0
    d=0.349d0
    ab=a*b

    n = Gridwk%n
    norbit = Orbit%norbit


    ! Fot the testing purpose
    !io=1
    !norbit = 1
    !alpha =0.1d0
    !DO i=2,n
    !    Orbit%wfn(i,1)=Grid%r(i)*exp(-alpha*Grid%r(i))
    !ENDDO
    !CALL extrapolate(Grid,Orbit%wfn(:,1))
    !Orbit%den(:)= (Grid%r(:)**2)*exp(-2*alpha*Grid%r(:))
    !Orbit%l(1) =0


    ALLOCATE(RealWfn(n,norbit),DerivRealWfn(n,norbit),Deriv2RealWfn(n,norbit))
    ALLOCATE(RealDen(n),DerivRealDen(n),Deriv2RealDen(n),LaplaceRealDen(n))
    ALLOCATE(eta(n),xi(n),dxidden(n),Derivxi(n),Deriv2xi(n),Laplacexi(n))
    ALLOCATE(Row1(n),Row2(n),Row3(n),Row4(n),Row5(n),Row6(n),Row7(n),&
&        decdchi(n),RowExtra(n))
    ALLOCATE(Derivxiden(n))
    ALLOCATE(dum(n),dum1(n),dum2(n))

    !-------Wave Function,Density,Eta,Xi Related----------------------------------------
    ! Wavefunction Related
    DO IterOrbit=1,norbit
       DO i=2,n
          RealWfn(i,IterOrbit) = Orbit%wfn(i,IterOrbit)/Gridwk%r(i)
       ENDDO
       CALL extrapolate(Gridwk,RealWfn(:,IterOrbit))
    ENDDO


    DO IterOrbit=1,norbit
       CALL derivative(Gridwk,RealWfn(:,IterOrbit),DerivRealWfn(:,IterOrbit),&
&           1,n)
       CALL derivative(Gridwk,DerivRealWfn(:,IterOrbit),&
&           Deriv2RealWfn(:,IterOrbit),1,n)
    ENDDO

    ! Density Related
    DO i=2,n
       RealDen(i)=Orbit%den(i)/(4.0d0*pi*(Gridwk%r(i)**2))
    ENDDO
    CALL extrapolate(Gridwk,RealDen)

    CALL derivative(Gridwk,RealDen,DerivRealDen,1,n)
    CALL derivative(Gridwk,DerivRealDen,Deriv2RealDen,1,n)
    LaplaceRealDen(:)=2*DerivRealDen(:)/Gridwk%r(:) +  Deriv2RealDen(:)
    !! Blows up at origin
    CALL extrapolate(Gridwk,LaplaceRealDen)   !! deal with the infinity at origin


    ! eta and Xi Related
    eta(:)= 1 + d/(RealDen(:)**(1.0d0/3.0d0))
    xi(:) =RealDen(:)**(5.0d0/3.0d0)*EXP(c*(RealDen(:)**(-1.0d0/3.0d0)))*eta(:)
    xi(:) = 1/xi(:)

    dxidden(:)=(-5.0d0/3.0d0)/RealDen(:)+c/(3.0d0*RealDen(:)**(4.0d0/3.0d0))&
&        + d/(3.0d0*eta(:)*RealDen(:)**(4.0d0/3.0d0))
    dxidden(:)=dxidden(:)*xi(:)

    CALL derivative(Gridwk,xi,Derivxi,1,n)
    CALL derivative(Gridwk,Derivxi,Deriv2xi,1,n)
    Laplacexi(:)=Deriv2xi(:) + 2*Derivxi(:)/Gridwk%r(:)
    CALL extrapolate(Gridwk,Laplacexi)


    dum(:)=xi(:)*RealDen(:)
    CALL derivative(Gridwk,dum,Derivxiden,1,n)

    !------------------------------------------------------------------

    ! Schrodinger equation test :
    !  dum(:)=Gridwk%r(:)*Deriv2RealWfn(:,io) + 2*DerivRealWfn(:,io)
    !  dum1(:)= (Gridwk%r(:)*Orbit%eig(io))*RealWfn(:,io)
    !  dum2(:)=dum(:)-dum1(:)


    ! Row 1
    Row1 = 0.0d0
    Row1(:)= -a/eta(:) - (a*d/3.0d0)*(RealDen(:)**(-1.0d0/3.0d0))/(eta(:)**2)

    Row1(:)= - Orbit%wfn(:,io)*Row1(:)
    Inte1=integrator(Gridwk,Row1,1,n)


    ! Row 2
    Row2 = 0.0d0
    DO IterOrbit=1,norbit
       IF(Orbit%l(IterOrbit) ==0 ) THEN
          Row2(:) = Row2(:) +2*Orbit%occ(IterOrbit)*(1.0d0/(4.0d0*pi))* &
&              RealDen(:)*(DerivRealWfn(:,IterOrbit)**2)
       ELSE
          ! Row2 needs to be accumlated , be careful !
          factor = Orbit%l(IterOrbit)*(Orbit%l(IterOrbit)+1)
          dum(:) = RealWfn(:,IterOrbit)/Gridwk%r(:)
          dum(:) = DerivRealWfn(:,IterOrbit)**2  + (factor*dum(:)**2)
          dum(:) = 2*Orbit%occ(IterOrbit)*(1.0d0/(4.0d0*pi))*RealDen(:)*dum(:)
          Row2(:) = Row2(:) + dum(:)
       ENDIF
    ENDDO

    Row2(:)=-ab*(0.25d0)*dxidden(:)*Orbit%wfn(:,io)*Row2(:)
    CALL extrapolate(Gridwk,Row2)
    Inte2=integrator(Gridwk,Row2,1,n)

    ! Row3
    Row3 = 0.0d0
    Row3(:) = - DerivRealDen(:)**2 + 0.5d0*RealDen(:)*LaplaceRealDen(:)
    Row3(:) = -ab*(0.25d0)*dxidden(:)*orbit%wfn(:,io)*Row3(:)
    Inte3 = integrator(Gridwk,Row3,1,n)

    ! Row4
    Row4 = 0.0d0
    Row4(:) = Gridwk%r(:)*DerivRealWfn(:,io)*Derivxiden(:)
    Row4(:) = ab*(0.5d0)*Row4(:)
    Inte4 = integrator(Gridwk,Row4,1,n)

    !  New Row5  combines part of the Row4
    Row5 = 0.0d0
    ! IF(Orbit%l(io) == 0 ) THEN
    !    Row5(:) = 0.0d0
    ! ELSE
    !    factor = orbit%l(io)*(orbit%l(io)+1)
    !    Row5(:)  = -ab*(0.5d0)*xi(:)*RealDen(:)*factor/(Gridwk%r(:)**2)
    ! ENDIF

    Row5(:)= ( Potwk%rv(:)  - Gridwk%r(:)*Orbit%eig(io) )*RealWfn(:,io)
    Row5(:)=ab*(0.50d0)*xi(:)*RealDen(:)*Row5(:)


    CALL extrapolate(Gridwk,Row5)

    dum(:)=Row4(:)+Row5(:)
    XiaoTest=integrator(Gridwk,dum)

    WRITE(*,*) "XiaoTest" , XiaoTest

    ! Row6
    Row6 = 0.0d0
    Row6(:) = derivxi(:)*derivRealDen(:)+xi(:)*LaplaceRealDen(:)
    Row6(:) = -ab*(0.25d0)*Row6(:)*orbit%wfn(:,io)
    Inte6=integrator(Gridwk,Row6,1,n)

    ! Row7
    Row7 = 0.0d0
    Row7(:) = -ab*(0.125d0)*(Laplacexi(:)*RealDen(:)+&
&        Derivxi(:)*DerivRealDen(:)+xi(:)*LaplaceRealDen(:))
    Row7(:) = -ab*(0.125d0)*xi(:)*LaplaceRealDen(:) + Row7(:)
    Row7(:) = Row7(:)*orbit%wfn(:,io)

    Inte7= integrator(Gridwk, Row7,1,n)


    ! Row Extra
    RowExtra = 0.0d0
    DO IterOrbit=1,norbit
       IF(Orbit%l(IterOrbit) ==0 ) THEN
          RowExtra(:) = RowExtra(:) +&
&              2*Orbit%occ(IterOrbit)*(1.0d0/(4.0d0*pi))* &
&              RealDen(:)*(DerivRealWfn(:,IterOrbit)**2)
       ELSE
          ! RowExtra needs to be accumlated , be careful !
          factor = Orbit%l(IterOrbit)*(Orbit%l(IterOrbit)+1)
          dum(:) = RealWfn(:,IterOrbit)/Gridwk%r(:)
          dum(:) = DerivRealWfn(:,IterOrbit)**2  + (factor*dum(:)**2)
          dum(:) = 2*Orbit%occ(IterOrbit)*(1.0d0/(4.0d0*pi))*dum(:)
          RowExtra(:) = RowExtra(:) + dum(:)
       ENDIF
    ENDDO

    RowExtra(:)=-RowExtra(:)*ab*(0.25d0)*xi(:)*orbit%wfn(:,io)
    CALL extrapolate(Gridwk, RowExtra)

    ! Sum up
    decdchi(:) = Row1(:) +  Row2(:) +  Row3(:) +  Row4(:) + &
&        Row5(:) +  Row6(:) +  Row7(:) + RowExtra(:)
    InteAll7 = integrator(Gridwk,decdchi,1,n)


    ! For Test Simple Case purpose Xiao
    ! dum(:)= -exp(-alpha*Gridwk%r(:))*(Gridwk%r(:)**2)
    ! dum(:)= -Gridwk%r(:)*Orbit%wfn(:,1);
    ! dum(:)=decdchi(:)*dum(:)
    ! XiaoTest=integrator(Gridwk,dum)



    res = decdchi

    DEALLOCATE(RealWfn,DerivRealWfn,Deriv2RealWfn, &
&        RealDen,DerivRealDen,Deriv2RealDen,LaplaceRealDen, &
&        eta,xi,dxidden,Derivxi,Deriv2xi,Laplacexi, &
&        Row1,Row2,Row3,Row4,Row5,Row6,Row7,decdchi,RowExtra,&
&        Derivxiden,dum,dum1,dum2)

  END SUBROUTINE Calc_decdchi_io

  !*************************************************************
  ! subroutine invert_numerov(Grid,l,mn,energy,zeroval,st,rv,res)
  !    subroutine similar to inhomo_bound_numerov but return inverted
  !        matrix of Numerov Coefficients -- not including r=0 point
  !        Also return bp== b(mn) coefficient for special boundary use
  !*************************************************************
  SUBROUTINE invert_numerov(Grid,l,mn,energy,rv,res,bp)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,mn
    REAL(8), INTENT(IN) :: energy,rv(:)
    REAL(8), INTENT(OUT) :: res(:,:),bp

    REAL(8), ALLOCATABLE :: aa(:),a(:),b(:),c(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale,zeroval,st
    INTEGER :: i,j,k,n

    ALLOCATE(a(mn-1),b(mn-1),c(mn-1),p(mn-1),aa(mn-1),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error inhomo_multi_numerov ', mn,i
       STOP
    ENDIF

    ! assume wfn(r) ~ r**(l+1) for r->0

    angm=l*(l+1)
    aa=0; res=0
    h=Grid%h;    h2=h*h
    zeroval=0
    IF (l==0)  zeroval=rv(1)
    IF (l==1)  zeroval=2
    st=rv(1)/2; st=st/(l+1)
    DO i=1,mn-1
       aa(i)=rv(i+1)/Grid%r(i+1)-energy+angm/(Grid%r(i+1)**2)
    ENDDO
    bp=rv(mn+1)/Grid%r(mn+1)-energy+angm/(Grid%r(mn+1)**2)

    xx=zeroval
    IF (usingloggrid(Grid)) THEN
       aa(1:mn-1)=0.25d0+Grid%rr02(2:mn)*aa(1:mn-1)
       bp=0.25d0+Grid%rr02(mn+1)*bp
       xx=Grid%rr02(1)*xx/Grid%pref(1)
    ENDIF
    bp=1.2d0-0.1d0*h2*bp

    res=0

    DO k=1,mn-1
       a(1:mn-1)=aa(1:mn-1)
       c=0; c(k)=1

       b=-2.4d0-h2*a   ;
       b(1)=b(1)-0.1d0*h2*xx/((Grid%r(2)**(l+1))*(1.d0+st*Grid%r(2)))
       a=1.2d0-0.1d0*h2*a
       p=0;p(2:mn-1)=a(1:mn-2)
       a(1:mn-2)=a(2:mn-1);a(mn-1)=0

       CALL thomas(mn-1,p,b,a,c)

       res(1:mn-1,k)=c(1:mn-1)
    ENDDO

    IF (usingloggrid(Grid)) THEN
       DO i=1,mn-1
          res(1:mn-1,i)=res(1:mn-1,i)*(Grid%pref(2:mn))
       ENDDO
       bp=bp*(Grid%pref(mn+1))
    ENDIF

    DEALLOCATE(a,b,c,p,aa)

  END SUBROUTINE invert_numerov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Version of VXOCC which just returns rvx given input wfn's
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE VXOCC(Grid,Orbit,HSZ,rvx)
     TYPE(GridInfo), INTENT(IN) :: Grid
     TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
     TYPE(HSZInfo), INTENT(INOUT) :: HSZ
     REAL(8), INTENT(INOUT) :: rvx(:)   ! new Vx*r

     INTEGER :: q,i,j,k,l,n,io,jo,ip,jp,last,mdim,adim,lastindex,norbit
     REAL(8) :: x,y,z
     REAL(8), ALLOCATABLE :: dum1(:),dum2(:),arg(:),der1(:),der2(:)
     REAL(8), ALLOCATABLE :: M(:,:),MM(:,:),ux(:),vx(:),Rp(:)
     INTEGER, ALLOCATABLE :: ind1(:),ind2(:),indm1(:),indm2(:)
     REAL(8), parameter :: scalefac=0.001d0,Rfix=20.d0
     INTEGER :: jcount=0
     CHARACTER(4) :: stuff

     n=Grid%n
     ALLOCATE(dum1(n),dum2(n),arg(n),der1(n),der2(n))

     mdim=norbit*norbit
     ALLOCATE(M(mdim,mdim),MM(mdim,mdim),ux(mdim),vx(mdim),&
&             Rp(mdim),ind1(mdim),ind2(mdim),indm1(mdim),indm2(mdim))

     last=0
     z=-1.d13
     Orbit%den=0
     Do io=1,norbit
        ! check normalization and sign of wfn
        arg=Orbit%wfn(:,io)**2
        x=sqrt(Integrator(Grid,arg))
        Orbit%wfn(:,io)=Orbit%wfn(:,io)/x
        CALL ADJUSTSIGN(Orbit%wfn(:,io),3)
        x=ABS(Orbit%wfn(n-3,io))
        if (Orbit%occ(io)>0.0001d0.and.x>z) then
             last=io
             z=x
        endif
        Orbit%den=Orbit%den+Orbit%occ(io)*(Orbit%wfn(:,io)**2)
     Enddo
     write(6,*) 'VXOCC ', jcount,last,z
     CALL SetIndex(Orbit)
     last=HSZ%zero_index
     write(6,*) 'Reset last = ', last

  ! calculate Xi(r) and <psi_q|X_p>
     adim=0;dum1=0;arg=0;ux=0;ind1=0;ind2=0;indm1=0;indm2=0

        Do io=1,norbit
           if (Orbit%occ(io)>0.00001d0) then
              call Calc_dexdphi_io(Grid,Orbit,io,der1);Orbit%X(:,io)=der1
              dum1=dum1+Orbit%occ(io)*der1*Orbit%wfn(:,io)
              do jo=1,norbit
                 if (Orbit%occ(jo)>0.00001d0.and.Orbit%l(io)==Orbit%l(jo)) then
                    adim=adim+1
                    ind1(adim)=io; ind2(adim)=jo
                    arg=Orbit%wfn(:,jo)*der1
                    arg(1)=0.d0; arg(2:n)=arg(2:n)/Grid%r(2:n)
                    ux(adim)=integrator(Grid,arg)
                    write(6,'("ux ", 3i5,1p,1e15.7)') io, jo, adim, ux(adim)
                     call flush(6)
                    if (jo==last.and.io==last) then
                       lastindex=adim
                       write(6,*) 'Last index = ', lastindex,io,jo
                    endif
                  endif
               enddo
            endif
         enddo

        ! Linear algebra setup
        dum1(1)=0
        do q=2,n
           if (Orbit%den(q)>machine_zero) then
              dum1(q)=dum1(q)/(Grid%r(q)*Orbit%den(q))
           else
              dum1(q)=0.d0
           endif
        enddo

        Rp=0; M=0; MM=0;k=0
        Do i=1,adim
           if (i/=lastindex) then
              k=k+1
              io=ind1(i);jo=ind2(i)
              indm1(k)=io; indm1(k)=jo
              arg=dum1*Orbit%wfn(:,io)*Orbit%wfn(:,jo)
              Rp(k)=integrator(Grid,arg)
              l=0; MM(k,k)=1.d0
              do j=1,adim
                 ip=ind1(j);jp=ind2(j)
                 if (j/=lastindex) then
                    l=l+1
                    arg=Orbit%wfn(:,io)*Orbit%wfn(:,jo)* &
&                      Orbit%wfn(:,ip)*Orbit%wfn(:,jp)*Orbit%occ(ip)
                   do q=2,n
                    if (Orbit%den(q)>machine_zero) then
                       arg(q)=arg(q)/(Orbit%den(q))
                    else
                       arg(q)=0.d0
                    endif
                   enddo
                   x=integrator(Grid,arg)
                   Rp(k)=Rp(k)-x*ux(j)
                   MM(k,l)=MM(k,l)-x
                 endif
              enddo
            endif
         enddo

    !Write(6,*) 'MM matrix'
    !Do k=1,adim-1
    !   write(6,'(1p,6e12.5)') (MM(k,l),l=1,adim-1)
    !enddo

    !write(6,*) 'Rp '
    !   write(6,'(1p,6e12.5)') (Rp(l),l=1,adim-1)

        call linsol(MM,Rp,adim-1,mdim,mdim,mdim)
        i=0;arg=0;HSZ%U=0
        Do k=1,adim
           if (k/=lastindex) then
              i=i+1
              vx(k)=Rp(i)
           else
              vx(k)=ux(k)
           endif
           io=ind1(k);jo=ind2(k)
           write(6,'(3i5," ux ", 1p,1e15.7,"  vx  ",1p,1e15.7)')k,io,jo,ux(k),vx(k)
           arg=arg+Orbit%occ(io)*(vx(k)-ux(k))*&
&                  (Orbit%wfn(:,io)*Orbit%wfn(:,jo))
           If (io==jo) HSZ%U(io)=vx(io)-ux(io)
        enddo
        arg(1)=0
        do i=2,n
           if (Orbit%den(i)>machine_zero) then
              arg(i)=arg(i)/Orbit%den(i)
           endif
        enddo
        rvx=(dum1+arg)*Grid%r

    !   call mkname(jcount,stuff)
    !   OPEN(unit=1001,file='VX'//TRIM(stuff),form='formatted')
    !   DO i = 1,n
    !      WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
    !&       rvx(i),dum1(i),arg(i),Orbit%den(i),&
    !&                 (Orbit%wfn(i,last))**2
    !   ENDDO
    !   CLOSE(1001)

       jcount=jcount+1

     DEALLOCATE(dum1,dum2,arg,der1,der2,M,MM,ux,vx,Rp,ind1,ind2,indm1,indm2)

  END SUBROUTINE VXOCC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  EXXOCCiter_wfn -- version using wfn mixing
!!!!   w is rvx
!!!!     On input it is assumed that Orbitwk contains current iteration wfn
!!!!     On output it is assumed that Orbitwk contains next iteration wfn
!!!!     Potwk%rvh and Potwk%rvx are calculated for input Orbitwk
!!!!     It is assumed that input wfn's are orthonormal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXXOCCiter_wfn(counter,mix,tol,err,success)
    INTEGER, INTENT(IN) :: counter
    REAL(8), INTENT(IN) :: mix,tol
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success

    INTEGER :: i,j,k,n,io,l,iter,last,lastchance
    REAL(8) :: x,y,z,fac,energy,val
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot
    INTEGER:: fcount=0
    REAL(8),save :: previouserror =1.d13
    CHARACTER(4) :: stuff
    LOGICAL :: OK,noalt=.true.,fix=.true.

    success=.FALSE. ; err=0.d0
    n=Gridwk%n

    lastchance=0
    fac=1.d0
    if (Potwk%nz<=29.d0+tol.and.Potwk%nz>=21.d0-tol) fac=0.5d0
    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in EXXOCCiter allocation' ,n,i
       STOP
    ENDIF

    write(6,*) 'Starting EXXOCCiter_wfn'; call flush(6)
    write(6,*) 'Initial eigenenergies'
    call One_electron_energy_Report(Orbitwk,6)

    Orbitwk%den=0     ! update density
    Do io=1,Orbitwk%norbit
       Orbitwk%den=Orbitwk%den+Orbitwk%occ(io)*(Orbitwk%wfn(:,io)**2)
    enddo

    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    ! replace   tmpPot%rvh   tmpPot%rvx
    call poisson(Gridwk,x,Orbitwk%den,tmpPot%rvh,y,z)
    Call VXOCC(Gridwk,Orbitwk,HSZ,tmpPot%rvx)

    rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
    tmpPot%rv=rv
    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,rv,OK)

    !if FC core calc , restore core info backinto tmpOrbit

    IF(frozencorecalculation) THEN
       write(6,*) 'Frozencore case'
       DO io = 1 , Orbitwk%norbit
          IF(Orbitwk%iscore(io)) THEN
             tmpOrbit%eig(io)=Orbitwk%eig(io)
             tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
          ENDIF
       ENDDO
    ENDIF

    write(6,*) 'EXXOCCIter ', fcount
    DO io = 1 , Orbitwk%norbit
       write(6,'(3i10,2x,1p,2e15.7)') io,tmpOrbit%np(io),tmpOrbit%l(io),&
&           tmpOrbit%occ(io),tmpOrbit%eig(io)
    ENDDO

    IF (.NOT.OK) THEN
       If (lastchance==0) then
           lastchance=1
           Call VXOCC(Gridwk,Orbitwk,HSZ,tmpPot%rvx)
           rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
           tmpPot%rv=rv
           CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,rv,OK)

            Write(6,*) 'Redo OCC '
             IF(frozencorecalculation) THEN
                write(6,*) 'Frozencore case'
                DO io = 1 , Orbitwk%norbit
                   IF(Orbitwk%iscore(io)) THEN
                      tmpOrbit%eig(io)=Orbitwk%eig(io)
                      tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
                   ENDIF
                ENDDO
             ENDIF

             write(6,*) 'EXXOCCIter ', fcount
             DO io = 1 , Orbitwk%norbit
                write(6,'(3i10,2x,1p,2e15.7)') io,tmpOrbit%np(io),&
&                    tmpOrbit%l(io),tmpOrbit%occ(io),tmpOrbit%eig(io)
             ENDDO

       ELSE
             WRITE(6,*) 'Bad luck in EXXIter'
             Potwk%rv=tmpPot%rv
             Potwk%rvh=tmpPot%rvh
             Potwk%rvx=tmpPot%rvx

             Orbitwk%wfn=tmpOrbit%wfn
             Orbitwk%eig=tmpOrbit%eig
             Orbitwk%den=tmpOrbit%den

             CALL Report_EXX_functions('ER')
             STOP
       ENDIF
    ENDIF

    !CALL mkname(fcount,stuff)
    !OPEN (unit=1001,file='occdiff.'//TRIM(stuff),form='formatted')
    !DO i=1,n
    !   WRITE(1001,'(1p,60E15.7)') Gridwk%r(i),Potwk%rvh(i),tmpPot%rvh(i),&
    !&       Potwk%rvx(i),tmpPot%rvx(i)
    !ENDDO
    !CLOSE(1001)

    !OPEN (unit=1001,file='occwfn.'//TRIM(stuff),form='formatted')
    !DO i=1,n
    !   WRITE(1001,'(1p,60E15.7)') Gridwk%r(i), &
    !&        (Orbitwk%wfn(i,io),io=1,Orbitwk%norbit),(tmpOrbit%wfn(i,io),io=1,Orbitwk%norbit)
    !ENDDO
    !CLOSE(1001)

    err=0
    do io=1,Orbitwk%norbit
       If(.not.frozencorecalculation.or.&
&           frozencorecalculation.and..not.Orbitwk%iscore(io)) then
          !CALL ADJUSTSIGN(Orbitwk%wfn(:,io),3)
          !CALL ADJUSTSIGN(tmpOrbit%wfn(:,io),3)
          dum=(Orbitwk%wfn(:,io)-tmpOrbit%wfn(:,io))**2
          err=err+Orbitwk%occ(io)*Integrator(Gridwk,dum)
       endif
    enddo

    ! update wfn if tolerance not satisfied
    IF (err>tol.OR.ABS(err-previouserror)>tol) THEN
       if (counter==40) then     ! reset
          write(6,*) '***Reset wfn***'
          DO io=1,Orbitwk%norbit
             If(.not.frozencorecalculation.or.&
&              frozencorecalculation.and..not.Orbitwk%iscore(io)) then
               Orbitwk%wfn(:,io)=tmpOrbit%wfn(:,io)
               Orbitwk%eig(io)=tmpOrbit%eig(io)
             endif
          ENDDO
       else
          val=1.d0-mix*fac
          WRITE(6,*) 'mixing wfns ', mix*fac,err
          DO io=1,Orbitwk%norbit
             If(.not.frozencorecalculation.or.&
&              frozencorecalculation.and..not.Orbitwk%iscore(io)) then
               Orbitwk%wfn(:,io)=val*Orbitwk%wfn(:,io)+&
&                        mix*fac*tmpOrbit%wfn(:,io)
               dum=(Orbitwk%wfn(:,io))**2
               x=sqrt(Integrator(Gridwk,dum))
               Orbitwk%wfn(:,io)=Orbitwk%wfn(:,io)/x
               Orbitwk%eig(io)=tmpOrbit%eig(io)       ! only true at convergence
             endif
          ENDDO
          !CALL ORTHONORMALIZE(Gridwk,Orbitwk)    !  normalize only
          write(6,*) 'OCCIter tol',fcount,tol,err,previouserror
       endif
    ELSE
       success=.TRUE.
       write(6,*) 'Calculation has converged with ', err,tol
    ENDIF

    Orbitwk%den=0
    Do io=1,Orbitwk%norbit
       Orbitwk%den=Orbitwk%den+Orbitwk%occ(io)*(Orbitwk%wfn(:,io)**2)
    Enddo

    Potwk%rv=tmpPot%rv
    Potwk%rvh=tmpPot%rvh
    Potwk%rvx=tmpPot%rvx

  !  calculate total energy for this new set; note that eone may be off

    Call Get_KinCoul(Gridwk,tmpPot,Orbitwk,SCFwk,noalt)

    write(6,*) 'Get_KinCoul', SCFwk%ekin


    If (frozencorecalculation) &
&           Call Get_FCKinCoul(Gridwk,tmpPot,Orbitwk,FCwk,SCFwk,noalt)

    write(6,*) 'Get_FCKinCoul', SCFwk%ekin

    CALL One_electron_energy_Report(Orbitwk,6)
    CALL Get_Energy_EXX(Gridwk,Orbitwk,SCFwk%eexc)
    SCFwk%etot=SCFwk%ekin+SCFwk%estatic+SCFwk%eexc
    IF (frozencorecalculation) THEN
       CALL Get_FCEnergy_EXX(Gridwk,Orbitwk,FCwk,SCFwk)
       energy=SCFwk%evale
       CALL Total_FCEnergy_Report(SCFwk,6)
       CALL Total_Energy_Report(SCFwk,6)
    ELSE
       energy=SCFwk%etot
       CALL Total_Energy_Report(SCFwk,6)
    ENDIF

    WRITE(6,'("ASiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err

    fcount=fcount+1
    previouserror=err
    DEALLOCATE (dum,res,rv)
    CALL DestroyPot(tmpPot)
    CALL DestroyOrbit(tmpOrbit)

  END SUBROUTINE   EXXOCCiter_wfn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  EXXKLIiter
!!!!   w is rvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXXKLIiter(w,energy,grad,err,success,verbose)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: verbose

    INTEGER :: i,j,k,n,io,l,iter,last
    REAL(8) :: x,fac
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot
    INTEGER:: fcount=0
    REAL(8), PARAMETER :: tol=1.d-10,mix=0.2d0
    CHARACTER(4) :: stuff
    LOGICAL :: OK

    success=.TRUE. ; err=0.d0
    n=Gridwk%n

    fac=1.d0
    if (Potwk%nz>=21.d0-tol.and.Potwk%nz<=29.d0+tol) fac=0.5d0
    HSZ%matchpoint=n
    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in EXXKLIiter allocation' ,n,i
       STOP
    ENDIF

    Potwk%rvx=w
    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    rv=Potwk%rvn+Potwk%rvh+Potwk%rvx
    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,rv,OK)

    !if FC core calc , restore core info backinto tmpOrbit

    IF(frozencorecalculation) THEN
       write(6,*) 'Frozencore case'
       DO io = 1 , Orbitwk%norbit
          IF(Orbitwk%iscore(io)) THEN
             tmpOrbit%eig(io)=Orbitwk%eig(io)
             tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
          ENDIF
       ENDDO
    ENDIF

    write(6,*) 'EXXKLIIter ', fcount
    DO io = 1 , Orbitwk%norbit
       write(6,'(3i10,2x,1p,2e15.7)') io,tmpOrbit%np(io),tmpOrbit%l(io),&
&           tmpOrbit%occ(io),Orbitwk%eig(io)
    ENDDO

    IF (.NOT.OK) THEN
       WRITE(6,*) 'Bad luck in EXXIter'
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

       CALL Report_EXX_functions('ER')
       STOP
    ENDIF

    tmpPot%rv=rv
    if (frozencorecalculation) then
      CALL Get_FCKinCoul(Gridwk,tmpPot,tmpOrbit,FCwk,SCFwk)
    else
      CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)
    endif
    dum=tmpPot%rvh-Potwk%rvh
    err=err+Dot_Product(dum,dum)
    write(6,*) 'HartreeIter ', fcount,err
    call simplemix(fac*mix,Potwk%rvh,tmpPot%rvh,dum)
    tmpPot%rvh=dum

    CALL One_electron_energy_Report(tmpOrbit,6)

    CALL Get_Energy_EXX(Gridwk,tmpOrbit,SCFwk%eexc)
    CALL SetIndex(tmpOrbit) ! find index of most extended wfn
    SCFwk%etot=SCFwk%ekin+SCFwk%estatic+SCFwk%eexc
    IF (frozencorecalculation) THEN
       CALL Get_FCEnergy_EXX(Gridwk,tmpOrbit,FCwk,SCFwk)
       energy=SCFwk%evale
       CALL Total_FCEnergy_Report(SCFwk,6)
    ELSE
       energy=SCFwk%etot
       CALL Total_Energy_Report(SCFwk,6)
    ENDIF

    call NEWVX(Gridwk,tmpOrbit,HSZ,Potwk%rvx,tmpPot%rvx,OK)
    res=tmpPot%rvx-Potwk%rvx
    x=Dot_Product(res,res)
    write(6,*) 'VxIter ', fcount,x; call flush(6)
    err=err+x;
    if (.not.OK) then
           err=err+100      ! penalty for bad tail
           !   reset potential to tmpPot%rvx
           w=tmpPot%rvx
           grad=tmpPot%rvx
    else
          call simplemix(fac*mix,Potwk%rvx,tmpPot%rvx,res)
          tmpPot%rvx=res; w=res
          grad=res
    endif
    success=OK

   ! CALL mkname(fcount,stuff)
   ! OPEN (unit=1001,file='klidiff.'//TRIM(stuff),form='formatted')
   ! DO i=1,n
   !    WRITE(1001,'(1p,60E15.7)') Gridwk%r(i),Potwk%rvh(i),tmpPot%rvh(i),&
   !&        Potwk%rvx(i),tmpPot%rvx(i)
   ! ENDDO
   ! Close(1001)

       tmpPot%rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

    WRITE(6,'("ASiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err
    WRITE(6,*) 'zero_index', fcount,HSZ%zero_index,HSZ%U(HSZ%zero_index)

    fcount=fcount+1
    DEALLOCATE (dum,res,rv)
    CALL DestroyPot(tmpPot)
    CALL DestroyOrbit(tmpOrbit)

  END SUBROUTINE   EXXKLIiter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  EXXKLIiter_innerloop
!!!!   w is rvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXXKLIiter_innerloop(w,energy,grad,err,success,verbose)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: verbose

    INTEGER :: i,j,k,n,io,l,iter,last
    REAL(8) :: x,fac
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot
    INTEGER:: fcount=0
    REAL(8), PARAMETER :: tol=1.d-10,mix=0.2d0
    CHARACTER(4) :: stuff
    LOGICAL :: OK

    success=.TRUE. ; err=0.d0   ; fac=1.d0
    if (HSZ%FIXED_ZERO) then
       fac=0.05d0
       write(6,*) '**** FIXED_ZERO calculation -- mix reduced to ', mix*fac
    endif
    n=Gridwk%n

    HSZ%matchpoint=n
    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in EXXKLIiter allocation' ,n,i
       STOP
    ENDIF

    Call Find_Hartree(w)
    Potwk%rvx=w
    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    rv=Potwk%rvn+Potwk%rvh+Potwk%rvx

    write(6,*) 'EXXKLIIter ', fcount
    DO io = 1 , Orbitwk%norbit
       write(6,'(3i10,2x,1p,2e15.7)') io,tmpOrbit%np(io),tmpOrbit%l(io),&
&           tmpOrbit%occ(io),Orbitwk%eig(io)
    ENDDO

    CALL One_electron_energy_Report(tmpOrbit,6)

    CALL Get_Energy_EXX(Gridwk,tmpOrbit,SCFwk%eexc)
    CALL SetIndex(tmpOrbit) ! find index of most extended wfn
    SCFwk%etot=SCFwk%ekin+SCFwk%estatic+SCFwk%eexc
    IF (frozencorecalculation) THEN
       CALL Get_FCKinCoul(Gridwk,Potwk,tmpOrbit,FCwk,SCFwk)
       CALL Get_FCEnergy_EXX(Gridwk,tmpOrbit,FCwk,SCFwk)
       energy=SCFwk%evale
       CALL Total_FCEnergy_Report(SCFwk,6)
    ELSE
       energy=SCFwk%etot
       CALL Total_Energy_Report(SCFwk,6)
    ENDIF

    call NEWVX(Gridwk,tmpOrbit,HSZ,Potwk%rvx,tmpPot%rvx,OK)
    res=tmpPot%rvx-Potwk%rvx
    x=Dot_Product(res,res)
    write(6,*) 'VxIter ', fcount,x; call flush(6)
    err=err+x;
    if (.not.OK) err=err+100  ! penalty for bad tail

    call simplemix(fac*mix,Potwk%rvx,tmpPot%rvx,res)
    tmpPot%rvx=res; w=res
    grad=res
    success=OK

   ! CALL mkname(fcount,stuff)
   ! OPEN (unit=1001,file='klidiff.'//TRIM(stuff),form='formatted')
   ! DO i=1,n
   !    WRITE(1001,'(1p,60E15.7)') Gridwk%r(i),Potwk%rvh(i),tmpPot%rvh(i),&
   !&        Potwk%rvx(i),tmpPot%rvx(i)
   ! ENDDO
   ! CLOSE(1001)

       tmpPot%rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

    WRITE(6,'("ASiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err
    WRITE(6,*) 'zero_index', fcount,HSZ%zero_index,HSZ%U(HSZ%zero_index)

    fcount=fcount+1
    DEALLOCATE (dum,res,rv)
    CALL DestroyOrbit(tmpOrbit)
    CALL DestroyPot(tmpPot)

  END SUBROUTINE   EXXKLIiter_innerloop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  EXXKLIiter_wfn -- version using wfn mixing
!!!!   w is rvx
!!!!     On input it is assumed that Orbitwk contains current iteration wfn
!!!!     On output it is assumed that Orbitwk contains next iteration wfn
!!!!     Potwk%rvh and Potwk%rvx are calculated for input Orbitwk
!!!!     It is assumed that input wfn's are orthonormal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXXKLIiter_wfn(counter,mix,tol,err,success)
    INTEGER, INTENT(IN) :: counter
    REAL(8), INTENT(IN) :: mix,tol
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success

    INTEGER :: i,j,k,n,io,l,iter,last,lastchance
    REAL(8) :: x,y,z,fac,energy,val
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot
    INTEGER:: fcount=0
    REAL(8),save :: previouserror =1.d13
    CHARACTER(4) :: stuff
    LOGICAL :: OK,noalt=.true.,fix=.true.

    success=.FALSE. ; err=0.d0
    n=Gridwk%n

    lastchance=0
    fac=1.d0
    if (Potwk%nz<=29+tol.and.Potwk%nz>=21-tol) fac=0.5d0
    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in EXXKLIiter allocation' ,n,i
       STOP
    ENDIF

    write(6,*) 'Starting EXXKLIiter_wfn'; call flush(6)
    write(6,*) 'Initial eigenenergies'
    call One_electron_energy_Report(Orbitwk,6)

    Orbitwk%den=0     ! update density
    Do io=1,Orbitwk%norbit
       Orbitwk%den=Orbitwk%den+Orbitwk%occ(io)*(Orbitwk%wfn(:,io)**2)
    enddo

    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    ! replace   tmpPot%rvh   tmpPot%rvx
    call poisson(Gridwk,x,Orbitwk%den,tmpPot%rvh,y,z)
    Call KLIVX(Gridwk,Orbitwk,HSZ,tmpPot%rvx)

    rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
    tmpPot%rv=rv
    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,rv,OK)

    !if FC core calc , restore core info backinto tmpOrbit

    IF(frozencorecalculation) THEN
       write(6,*) 'Frozencore case'
       DO io = 1 , Orbitwk%norbit
          IF(Orbitwk%iscore(io)) THEN
             tmpOrbit%eig(io)=Orbitwk%eig(io)
             tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
          ENDIF
       ENDDO
    ENDIF

    write(6,*) 'EXXKLIIter ', fcount
    DO io = 1 , Orbitwk%norbit
       write(6,'(3i10,2x,1p,2e15.7)') io,tmpOrbit%np(io),tmpOrbit%l(io),&
&           tmpOrbit%occ(io),tmpOrbit%eig(io)
    ENDDO

    IF (.NOT.OK) THEN
       If (lastchance==0) then
           lastchance=1
           Call KLIVX(Gridwk,Orbitwk,HSZ,tmpPot%rvx,fix)
           rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
           tmpPot%rv=rv
           CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,rv,OK)

            Write(6,*) 'Redo KLI '
             IF(frozencorecalculation) THEN
                write(6,*) 'Frozencore case'
                DO io = 1 , Orbitwk%norbit
                   IF(Orbitwk%iscore(io)) THEN
                      tmpOrbit%eig(io)=Orbitwk%eig(io)
                      tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
                   ENDIF
                ENDDO
             ENDIF

             write(6,*) 'EXXKLIIter ', fcount
             DO io = 1 , Orbitwk%norbit
                write(6,'(3i10,2x,1p,2e15.7)') io,tmpOrbit%np(io),&
&                    tmpOrbit%l(io),tmpOrbit%occ(io),tmpOrbit%eig(io)
             ENDDO

       ELSE
             WRITE(6,*) 'Bad luck in EXXIter'
             Potwk%rv=tmpPot%rv
             Potwk%rvh=tmpPot%rvh
             Potwk%rvx=tmpPot%rvx

             Orbitwk%wfn=tmpOrbit%wfn
             Orbitwk%eig=tmpOrbit%eig
             Orbitwk%den=tmpOrbit%den

             CALL Report_EXX_functions('ER')
             STOP
       ENDIF
    ENDIF

   ! CALL mkname(fcount,stuff)
   ! OPEN (unit=1001,file='klidiff.'//TRIM(stuff),form='formatted')
   ! DO i=1,n
   !    WRITE(1001,'(1p,60E15.7)') Gridwk%r(i),Potwk%rvh(i),tmpPot%rvh(i),&
   !&         Potwk%rvx(i),tmpPot%rvx(i)
   ! ENDDO
   ! CLOSE(1001)

   ! OPEN (unit=1001,file='kliwfn.'//TRIM(stuff),form='formatted')
   ! DO i=1,n
   !    WRITE(1001,'(1p,60E15.7)') Gridwk%r(i), &
   !&         (Orbitwk%wfn(i,io),io=1,tmpOrbit%norbit),(tmpOrbit%wfn(i,io),io=1,tmpOrbitnorbit)
   ! ENDDO
   ! CLOSE(1001)

    err=0
    do io=1,Orbitwk%norbit
       If(.not.frozencorecalculation.or.&
&           frozencorecalculation.and..not.Orbitwk%iscore(io)) then
          !CALL ADJUSTSIGN(Orbitwk%wfn(:,io),3)
          !CALL ADJUSTSIGN(tmpOrbit%wfn(:,io),3)
          dum=(Orbitwk%wfn(:,io)-tmpOrbit%wfn(:,io))**2
          err=err+Orbitwk%occ(io)*Integrator(Gridwk,dum)
       endif
    enddo

    ! update wfn if tolerance not satisfied
    IF (err>tol.OR.ABS(err-previouserror)>tol) THEN
       if (counter==40) then     ! reset
          write(6,*) '***Reset wfn***'
          DO io=1,Orbitwk%norbit
             If(.not.frozencorecalculation.or.&
&              frozencorecalculation.and..not.Orbitwk%iscore(io)) then
               Orbitwk%wfn(:,io)=tmpOrbit%wfn(:,io)
               Orbitwk%eig(io)=tmpOrbit%eig(io)
             endif
          ENDDO
       else
          val=1.d0-mix*fac
          WRITE(6,*) 'mixing wfns ', mix*fac,err
          DO io=1,Orbitwk%norbit
             If(.not.frozencorecalculation.or.&
&              frozencorecalculation.and..not.Orbitwk%iscore(io)) then
               Orbitwk%wfn(:,io)=val*Orbitwk%wfn(:,io)+&
&                        mix*fac*tmpOrbit%wfn(:,io)
               dum=(Orbitwk%wfn(:,io))**2
               x=sqrt(Integrator(Gridwk,dum))
               Orbitwk%wfn(:,io)=Orbitwk%wfn(:,io)/x
               Orbitwk%eig(io)=tmpOrbit%eig(io)       ! only true at convergence
             endif
          ENDDO
          !CALL ORTHONORMALIZE(Gridwk,Orbitwk)    !  normalize only
          write(6,*) 'KLIIter tol',fcount,tol,err,previouserror
       endif
    ELSE
       success=.TRUE.
       write(6,*) 'Calculation has converged with ', err,tol
    ENDIF

    Orbitwk%den=0
    Do io=1,Orbitwk%norbit
       Orbitwk%den=Orbitwk%den+Orbitwk%occ(io)*(Orbitwk%wfn(:,io)**2)
    Enddo

    Potwk%rv=tmpPot%rv
    Potwk%rvh=tmpPot%rvh
    Potwk%rvx=tmpPot%rvx


  !  calculate total energy for this new set; note that eone may be off

    Call Get_KinCoul(Gridwk,tmpPot,Orbitwk,SCFwk,noalt)

    write(6,*) 'Get_KinCoul', SCFwk%ekin


    If (frozencorecalculation) &
&           Call Get_FCKinCoul(Gridwk,tmpPot,Orbitwk,FCwk,SCFwk,noalt)

    write(6,*) 'Get_FCKinCoul', SCFwk%ekin

    CALL One_electron_energy_Report(Orbitwk,6)
    CALL Get_Energy_EXX(Gridwk,Orbitwk,SCFwk%eexc)
    SCFwk%etot=SCFwk%ekin+SCFwk%estatic+SCFwk%eexc
    IF (frozencorecalculation) THEN
       CALL Get_FCEnergy_EXX(Gridwk,Orbitwk,FCwk,SCFwk)
       energy=SCFwk%evale
       CALL Total_FCEnergy_Report(SCFwk,6)
       CALL Total_Energy_Report(SCFwk,6)
    ELSE
       energy=SCFwk%etot
       CALL Total_Energy_Report(SCFwk,6)
    ENDIF

    WRITE(6,'("ASiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err

    fcount=fcount+1
    previouserror=err
    DEALLOCATE (dum,res,rv)
    CALL DestroyOrbit(tmpOrbit)
    CALL DestroyPot(tmpPot)

  END SUBROUTINE   EXXKLIiter_wfn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Version of KLIVX which just returns rvx given input wfn's
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE KLIVX(Grid,Orbit,HSZ,rvx,fix)
     TYPE(GridInfo), INTENT(IN) :: Grid
     TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
     TYPE(HSZInfo), INTENT(INOUT) :: HSZ
     REAL(8), INTENT(INOUT) :: rvx(:)   ! new Vx*r
     LOGICAL, OPTIONAL :: fix

     INTEGER :: i,j,k,n,io,jo,norbit,last
     REAL(8) :: x,y,z
     REAL(8), ALLOCATABLE :: dum1(:),dum2(:),arg(:),der1(:),der2(:)
     REAL(8), ALLOCATABLE :: M(:,:),MM(:,:),ux(:),vx(:),Rp(:)
     REAL(8), parameter :: scalefac=0.001d0,Rfix=20.d0
     INTEGER :: jcount=0
     CHARACTER(4) :: stuff

     n=Grid%n
     ALLOCATE(dum1(n),dum2(n),arg(n),der1(n),der2(n))

     norbit=Orbit%norbit
     ALLOCATE(M(norbit,norbit),MM(norbit-1,norbit-1),ux(norbit),vx(norbit),&
&             Rp(norbit-1))

     last=0
     z=-1.d13
     Orbit%den=0
     Do io=1,norbit
        ! check normalization and sign of wfn
        arg=Orbit%wfn(:,io)**2
        x=sqrt(Integrator(Grid,arg))
        Orbit%wfn(:,io)=Orbit%wfn(:,io)/x
        CALL ADJUSTSIGN(Orbit%wfn(:,io),3)
        x=ABS(Orbit%wfn(n-3,io))
        if (Orbit%occ(io)>0.0001d0.and.x>z) then
             last=io
             z=x
        endif
        Orbit%den=Orbit%den+Orbit%occ(io)*(Orbit%wfn(:,io)**2)
     Enddo

     write(6,*) 'KLIVX ', jcount,last,z
     CALL SetIndex(Orbit)
     last=HSZ%zero_index
     write(6,*) 'Reset last = ', last

     If (present(fix)) then
        CALL InitialVx(Grid,Orbit%wfn(:,last),rvx)
     else
        dum1=0;  arg=0; M=0;
        Do io=1,norbit
           call Calc_dexdphi_io(Grid,Orbit,io,der1);Orbit%X(:,io)=der1
           dum1=dum1+Orbit%occ(io)*der1*Orbit%wfn(:,io)
           der1=der1*Orbit%wfn(:,io)
           der1(1)=0.d0; der1(2:n)=der1(2:n)/Grid%r(2:n)
           ux(io)=integrator(Grid,der1)
           write(6,*) 'io const' , io,ux(io) ; call flush(6)
           do jo=1,norbit
              arg=0
              arg(2:n)=(Orbit%wfn(2:n,io)*Orbit%wfn(2:n,jo))**2/Orbit%den(2:n)
              M(io,jo)=integrator(Grid,arg)
           enddo
        Enddo

        dum1(1)=0; dum2=0.d0

        do i=2,n
           if (Orbit%den(i)>machine_zero) then
              dum1(i)=dum1(i)/Orbit%den(i)/Grid%r(i)
           else
             dum1(i)=0
           endif
        enddo
        call extrapolate(Grid,dum1)

        MM=0
        i=0
        do io=1,norbit
           if(io/=last) then
             i=i+1;j=0
             MM(i,i)=1.d0
             arg=(dum1+dum2)*(Orbit%wfn(:,io)**2)
             Rp(i)=integrator(Grid,arg)
             do jo=1,norbit
                if(jo/=last)  then
                  j=j+1
                  Rp(i)=Rp(i)-M(io,jo)*Orbit%occ(jo)*ux(jo)
                  MM(i,j)=MM(i,j)-M(io,jo)*Orbit%occ(jo)
                endif
             enddo
           endif
        enddo

        call linsol(MM,Rp,norbit-1,norbit,norbit,norbit-1)
        i=0;arg=0
        Do io=1,norbit
           if (io/=last) then
              i=i+1
              vx(io)=Rp(i)
           else
              vx(io)=ux(io)
           endif
           write(6,'(i5," ux ", 1p,1e15.7,"  vx  ",1p,1e15.7)') io,ux(io),vx(io)
           arg=arg+Orbit%occ(io)*(vx(io)-ux(io))*(Orbit%wfn(:,io)**2)
           HSZ%U(io)=vx(io)-ux(io)
        enddo
        arg(1)=0
        do i=2,n
           if (Orbit%den(i)>machine_zero) then
              arg(i)=arg(i)/Orbit%den(i)
           endif
        enddo
        call extrapolate(Grid,arg)
        rvx=(dum1+dum2+arg)*Grid%r

      ! call mkname(jcount,stuff)
      ! OPEN(unit=1001,file='VX'//TRIM(stuff),form='formatted')
      ! DO i = 1,n
      !    WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
      !&     rvx(i),dum1(i),dum2(i),arg(i),Orbit%den(i),&
      !&               (Orbit%wfn(i,last))**2
      ! ENDDO
      ! CLOSE(1001)

   Endif
       jcount=jcount+1

     DEALLOCATE(dum1,dum2,arg,der1,der2,M,MM,ux,vx,Rp)

  END SUBROUTINE KLIVX

  SUBROUTINE Find_Hartree(rvx)
    REAL(8), INTENT(INOUT) :: rvx(:)

    TYPE(Anderson_Context) :: A_H
    REAL(8), ALLOCATABLE :: arg(:)
    REAL(8) :: energy
    LOGICAL :: OK,update

    ALLOCATE(arg(Gridwk%n))
    CALL InitAnderson_dr(A_H,6,5,Gridwk%n,0.5d0,1.d3,200,1.d-8,1.d-12,.FALSE.)
    arg=Potwk%rvh
    Potwk%rvx=rvx
    CALL DoAndersonMix(A_H,arg,energy,Hsub,OK)

    WRITE(6,*) ' Find_Hartree completed --  iter, delta = ',&
&        A_H%CurIter, A_H%res

    CALL FreeAnderson(A_H)
    DEALLOCATE(arg)
  END SUBROUTINE Find_Hartree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   In Hsub w is the Hartree potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE  Hsub(w,energy,residue,err,OK,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: residue(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: OK
    LOGICAL, INTENT(IN)  :: update

    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot
    REAL(8), ALLOCATABLE :: rv(:)
    INTEGER :: i,j,k,n,io

    ALLOCATE(rv(Gridwk%n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in Hsub allocation' ,Gridwk%n,i
       STOP
    ENDIF

    WRITE(6,*) 'Starting Hartree '
    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    rv=tmpPot%rvn+w+tmpPot%rvx

    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,rv,OK)
    IF (.NOT.OK) THEN
       WRITE(6,*) 'Bad luck in Hsub'
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

       CALL Report_EXX_functions('ER')
       STOP
    ENDIF


    tmpPot%rv=rv

    !if FC core calc , restore core info backinto tmpOrbit

    IF(frozencorecalculation) THEN
       DO io = 1 , Orbitwk%norbit
          IF(Orbitwk%iscore(io)) THEN
             tmpOrbit%eig(io)=Orbitwk%eig(io)
             tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
          ENDIF
       ENDDO
    ENDIF

    CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)

    residue=tmpPot%rvh-w
    err=DOT_PRODUCT(residue,residue)
    ! tmpPot%rvh has been updated; find new tmpPot%rv
    energy=SCFwk%ecoul

    IF (update) THEN
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

       !Call One_electron_energy_Report(Orbitwk,6)
    ENDIF

    !call writestuff
    CALL DestroyOrbit(tmpOrbit)
    DEALLOCATE (rv,tmpPot%rvn,tmpPot%rv,tmpPot%rvh,tmpPot%rvx)

  END SUBROUTINE Hsub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  EXXIter
!!!!   w is rvx
!!!!         Use shift to update w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXXiter(w,energy,grad,err,success,verbose)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: verbose

    INTEGER :: i,j,k,n,io,l,iter,last
    REAL(8) :: x
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot
    INTEGER:: fcount=0
    TYPE (OrbitInfo), TARGET :: eigOrbit
    REAL(8), PARAMETER :: tol=1.d-10,mix=0.5d0,scaletol=0.01d0
    CHARACTER(4) :: stuff
    LOGICAL :: OK

    success=.TRUE. ; err=0.d0
    n=Gridwk%n

    HSZ%matchpoint=n
    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in EXXIter allocation' ,n,i
       STOP
    ENDIF

    Potwk%rvx=w
    Potwk%rv=Potwk%rvn+Potwk%rvh+Potwk%rvx

    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    rv=Potwk%rvn+Potwk%rvh+Potwk%rvx
    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,rv,OK)
    IF (.NOT.OK) THEN
       WRITE(6,*) 'Bad luck in EXXIter'
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

       CALL Report_EXX_functions('ER')
       STOP
    ENDIF

    tmpPot%rv=rv
    CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)
    dum=tmpPot%rvh-Potwk%rvh
    err=err+Dot_Product(dum,dum)
    write(6,*) 'HartreeIter ', fcount,err
    call simplemix(mix,Potwk%rvh,tmpPot%rvh,dum)
    tmpPot%rvh=dum

    !write(6,*) 'Before one electron report'; call flush(6)
    CALL One_electron_energy_Report(tmpOrbit,6)
    !write(6,*) 'After one electron report'; call flush(6)

    CALL   CopyOrbit(tmpOrbit,EigOrbit)
    CALL SetIndex(tmpOrbit)
    !write(6,*) 'Before Calc_dedv report'; call flush(6)
    CALL Calc_dedv(Gridwk,Potwk,tmpOrbit,EigOrbit)

    !call mkname(fcount,stuff)
    !OPEN(unit=1001,file='psiii'//TRIM(stuff),form='formatted')
    !DO i = 1,n
    !   WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
    !&       (HSZ%psi(i,io),io=1,Orbitwk%norbit)
    !ENDDO
    !CLOSE(1001)

    !OPEN(unit=1001,file='shifti'//TRIM(stuff),form='formatted')
    !DO i = 1,n
    !   WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), HSZ%shift(i),&
    !&       (HSZ%psi(i,io)*tmpOrbit%wfn(i,io),io=1,Orbitwk%norbit)
    !ENDDO
    !CLOSE(1001)

    CALL Get_Energy_EXX(Gridwk,tmpOrbit,SCFwk%eexc)
    CALL SetIndex(tmpOrbit) ! find index of most extended wfn
    SCFwk%etot=SCFwk%ekin+SCFwk%estatic+SCFwk%eexc
    IF (frozencorecalculation) THEN
       CALL Get_FCKinCoul(Gridwk,Potwk,tmpOrbit,FCwk,SCFwk)
       CALL Get_FCEnergy_EXX(Gridwk,tmpOrbit,FCwk,SCFwk)
       energy=SCFwk%evale
       CALL Total_FCEnergy_Report(SCFwk,6)
    ELSE
       energy=SCFwk%etot
       CALL Total_Energy_Report(SCFwk,6)
    ENDIF

    grad=-HSZ%shift
    x=Dot_Product(grad,grad)
    write(6,*) 'VxIter ', fcount,x; call flush(6)
    err=err+x
    x=sqrt(x)
    if (x>scaletol) then
           x=scaletol/x
       else
           x=1.d0
    endif
    w=w+mix*x*grad
    tmpPot%rvx=w


    !CALL mkname(fcount,stuff)
    !OPEN (unit=1001,file='hfwfni.'//TRIM(stuff),form='formatted')
    !DO i=1,n
    !   WRITE(1001,'(1p,60E15.7)') Gridwk%r(i),Orbitwk%den(i),tmpOrbit%den(i),&
    !&      (Orbitwk%wfn(i,io),&
    !&       tmpOrbit%wfn(i,io),io=1,Orbitwk%norbit)
    !ENDDO
    !CLOSE(1001)

    !OPEN (unit=1001,file='poti.'//TRIM(stuff),form='formatted')
    !DO i=1,n
    !   WRITE(1001,'(1p,60E15.7)') Gridwk%r(i),Potwk%rvh(i),tmpPot%rvh(i),&
    !&      Potwk%rvx(i),tmpPot%rvx(i),grad(i)
    !ENDDO
    !CLOSE(1001)

       tmpPot%rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

    WRITE(6,'("ASiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err
    WRITE(6,*) 'zero_index', fcount,HSZ%zero_index,HSZ%U(HSZ%zero_index)

    fcount=fcount+1
    CALL DestroyOrbit(eigOrbit)
    CALL DestroyOrbit(tmpOrbit)
    DEALLOCATE(dum,res,rv)
    DEALLOCATE (tmpPot%rvn,tmpPot%rv,tmpPot%rvh,tmpPot%rvx)

  END SUBROUTINE   EXXiter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  AdjustV   -- simplified version
!!!!   Return rvx = (w+cons*r)*theta+rvx0*(1-theta)
!!!!     theta=phiL(r)/phiL(rm) for r>rm and 1 for r<rm
!!!!   So that outer orbital constraint is satisfied
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE AdjustV(Grid,Orbit,w)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    REAL(8), INTENT(INOUT) :: w(:)

    INTEGER :: i,j,k,n,io,last,l,point,begin
    REAL(8),ALLOCATABLE :: rvx(:),rvx0(:),f(:),arg1(:),arg2(:),arg3(:),dum(:)
    REAL(8), POINTER :: r(:)
    REAL(8) :: energy,val,x,Numerator,beta,rc
    REAL(8), PARAMETER :: smll=1.d-9,lsmll=1.d-7
    INTEGER, PARAMETER :: many=5
    INTEGER, SAVE :: matchpt
    LOGICAL :: success
    INTEGER :: counter=0
    CHARACTER(4) :: stuff

    !write(6,*) 'AdjustV not used in this verison'
    !stop
    CALL SetIndex(Orbit)
    last=HSZ%zero_index
    if (ABS(HSZ%U(last))<=smll) return
    n=Grid%n
    r=>Grid%r
    ALLOCATE(rvx0(n),rvx(n),f(n),arg1(n),arg2(n),arg3(n),dum(n))

    rvx=0

    CALL ApproxVx(Grid,Orbit,last,rvx0)
    f=0
    f(:)=Orbit%wfn(:,last)
    if (f(n)<0.d0) f(:)=-f(:)
    do i=n,2,-1
       j=i
       if (f(i-1)<f(i)) exit
    enddo
    point=j
    write(6,*) 'Simple AdjustV ', point,Grid%r(point) ; call flush(6)
    f(1:j-1)=1.d0
    f(j:n)=f(j:n)/f(j)
    HSZ%matchpoint=point
    matchpt=point
    rc=Grid%r(point)
    dum=w(:)*f(:)+rvx0(:)*(1.d0-f(:))
    arg1=dum*(Orbitwk%wfn(:,last)**2)
    arg2=f*(Orbitwk%wfn(:,last)**2)

    CALL Calc_dexdphi_io(Gridwk,Orbitwk,last,arg3)

    arg1=arg3*Orbitwk%wfn(:,last)-arg1
    arg1(1)=0; arg1(2:n)=arg1(2:n)/Grid%r(2:n)
    val=integrator(Grid,arg1)/integrator(Gridwk,arg2)

    WRITE(6,*) 'AdjustV -- ', val

    rvx=0
    rvx=dum+val*Grid%r*f

    write(6,*) 'AdjustV diff:',SUM(ABS(rvx-w))

   ! CALL mkname(counter,stuff)
   ! OPEN (unit=1001,file='adjust'//TRIM(stuff),form='formatted')
   ! DO i=1,n
   !    WRITE(1001,'(1p,15E15.7)') r(i),w(i),rvx(i),rvx0(i),f(i)
   ! ENDDO
   ! CLOSE(1001)
    w=rvx

    counter=counter+1
    DEALLOCATE(rvx0,rvx,f,arg1,arg2,arg3,dum)

  END SUBROUTINE AdjustV

   SUBROUTINE simplemix(mix,v0,v1,vmix)
      REAL(8), INTENT(IN) :: mix
      REAL(8), INTENT(IN) :: v0(:),v1(:)
      REAL(8), INTENT(OUT) :: vmix(:)

      vmix=v0*(1.d0-mix)+v1*mix
   END SUBROUTINE simplemix

  SUBROUTINE NEWVX(Grid,Orbit,HSZ,rvx0,rvx,success)
     TYPE(GridInfo), INTENT(IN) :: Grid
     TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
     TYPE(HSZInfo), INTENT(INOUT) :: HSZ
     REAL(8), INTENT(IN) :: rvx0(:)     ! original Vx*r
     REAL(8), INTENT(INOUT) :: rvx(:)   ! new Vx*r
     LOGICAL, INTENT(OUT) :: success

     INTEGER :: i,j,k,n,io,jo,norbit,last
     REAL(8) :: x,y,z
     REAL(8), ALLOCATABLE :: dum1(:),dum2(:),arg(:),der1(:),der2(:)
     REAL(8), ALLOCATABLE :: M(:,:),MM(:,:),ux(:),vx(:),Rp(:)
     REAL(8), parameter :: scalefac=0.001d0,Rfix=20.d0
     INTEGER :: jcount=0
     CHARACTER(4) :: stuff

     n=Grid%n
     ALLOCATE(dum1(n),dum2(n),arg(n),der1(n),der2(n))

     norbit=Orbit%norbit
     ALLOCATE(M(norbit,norbit),MM(norbit-1,norbit-1),ux(norbit),vx(norbit),&
&             Rp(norbit-1))

     CALL SetIndex(Orbit)
     last=HSZ%zero_index

     Orbit%den=0
     Do io=1,norbit
        Orbit%den=Orbit%den+Orbit%occ(io)*(Orbit%wfn(:,io)**2)
     Enddo

     dum1=0;  arg=0; M=0;
     Do io=1,norbit
        call Calc_dexdphi_io(Grid,Orbit,io,der1);Orbit%X(:,io)=der1
        dum1=dum1+Orbit%occ(io)*der1*Orbit%wfn(:,io)
        der1=der1*Orbit%wfn(:,io)
        der1(1)=0.d0; der1(2:n)=der1(2:n)/Grid%r(2:n)
        ux(io)=integrator(Grid,der1)
        write(6,*) 'io const' , io,ux(io) ; call flush(6)
        do jo=1,norbit
           arg=0
           arg(2:n)=(Orbit%wfn(2:n,io)*Orbit%wfn(2:n,jo))**2/Orbit%den(2:n)
           M(io,jo)=integrator(Grid,arg)
        enddo
     Enddo

     dum1(1)=0; dum2=0.d0

     do i=2,n
        if (Orbit%den(i)>machine_zero) then
           dum1(i)=dum1(i)/Orbit%den(i)/Grid%r(i)
        else
          dum1(i)=0
        endif
     enddo
     call extrapolate(Grid,dum1)

     MM=0
     i=0
     do io=1,norbit
        if(io/=last) then
          i=i+1;j=0
          MM(i,i)=1.d0
          arg=(dum1+dum2)*(Orbit%wfn(:,io)**2)
          Rp(i)=integrator(Grid,arg)
          do jo=1,norbit
             if(jo/=last)  then
               j=j+1
               Rp(i)=Rp(i)-M(io,jo)*Orbit%occ(jo)*ux(jo)
               MM(i,j)=MM(i,j)-M(io,jo)*Orbit%occ(jo)
             endif
          enddo
        endif
     enddo

     call linsol(MM,Rp,norbit-1,norbit,norbit,norbit-1)
     i=0;arg=0
     Do io=1,norbit
        if (io/=last) then
           i=i+1
           vx(io)=Rp(i)
        else
           vx(io)=ux(io)
        endif
        write(6,'(i5," ux ", 1p,1e15.7,"  vx  ",1p,1e15.7)') io,ux(io),vx(io)
        arg=arg+Orbit%occ(io)*(vx(io)-ux(io))*(Orbit%wfn(:,io)**2)
        HSZ%U(io)=vx(io)-ux(io)
     enddo
     arg(1)=0
     do i=2,n
        if (Orbit%den(i)>machine_zero) then
           arg(i)=arg(i)/Orbit%den(i)
        endif
     enddo
     call extrapolate(Grid,arg)
     rvx=(dum1+dum2+arg)*Grid%r

     der2=rvx
     CALL Fixrvx(Grid,Orbit,der2,rvx,success)

   ! call mkname(jcount,stuff)
   ! OPEN(unit=1001,file='VX'//TRIM(stuff),form='formatted')
   ! DO i = 1,n
   !    WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
   !&        rvx0(i),rvx(i),der2(i),der1(i),dum1(i),dum2(i),arg(i),Orbit%den(i),&
   !&                  (Orbit%wfn(i,last))**2
   ! ENDDO
   ! CLOSE(1001)
    jcount=jcount+1

    !if (FULL) stop
     DEALLOCATE(dum1,dum2,arg,der1,der2,M,MM,ux,vx,Rp)

  END SUBROUTINE NEWVX

    SUBROUTINE Fixrvx(Grid,Orbit,rvxin,rvxout,success)
      Type(GridInfo), INTENT(IN) :: Grid
      Type(OrbitInfo), INTENT(IN) :: Orbit
      Real(8), INTENT(IN) :: rvxin(:)
      Real(8), INTENT(OUT) :: rvxout(:)
      Logical, INTENT(OUT) :: success

      INTEGER :: i,j,k,n,last
      INTEGER :: jcount=0
      REAL(8) :: x,y,z
      REAL(8), parameter :: Rfix=20.0d0, err=0.01d0

      success=.true.
      n=Grid%n
      CALL SetIndex(Orbit)
      last=HSZ%zero_index

      rvxout=0
      CALL InitialVx(Grid,Orbit%wfn(:,last),rvxout)
      x=ABS(rvxin(n)-rvxout(n))
      write(6,*) 'KLI conv ', jcount, x
      if (x>err) then
         !open(1001,file='badidea',form='formatted')
         success=.false.; write(6,*) 'Warning KLI not converging'
         j=FindGridIndex(Grid,Rfix)
         Do i=1,j
            !write(1001,'(1p,60E15.7)') Grid%r(i),rvxin(i),rvxout(i)
            rvxout(i)=rvxin(i)
         Enddo
         Do i=j+1,n
            z=Pi*(Grid%r(i)-Grid%r(j))/(Grid%r(n)-Grid%r(j))
            y=1.d0-(sin(z)/z)**8
            !write(1001,'(1p,60E15.7)') Grid%r(i),rvxin(i),rvxout(i),&
            !&        rvxout(i)*y+rvxin(i)*(1-y),y
            rvxout(i)=rvxout(i)*y+rvxin(i)*(1-y)
         Enddo
         !close(1001)
         !stop
      else
         rvxout=rvxin
      endif
      jcount=jcount+1
    END SUBROUTINE Fixrvx

  SUBROUTINE NEWAdjustVX(Grid,Orbit,HSZ,rvx0,rvx)
     TYPE(GridInfo), INTENT(IN) :: Grid
     TYPE(OrbitInfo), INTENT(IN) :: Orbit
     TYPE(HSZInfo), INTENT(IN) :: HSZ
     REAL(8), INTENT(IN) :: rvx0(:)     ! original Vx*r
     REAL(8), INTENT(INOUT) :: rvx(:)   ! new Vx*r
     INTEGER :: i,j,k,n,io,jo,norbit,last
     REAL(8) :: x,y,z
     REAL(8), ALLOCATABLE :: dum1(:),dum2(:),arg(:),der1(:),der2(:)
     REAL(8), ALLOCATABLE :: M(:,:),MM(:,:),ux(:),vx(:),Rp(:)
     REAL(8), parameter :: scalefac=0.001d0
     INTEGER :: jcount=0
     CHARACTER(4) :: stuff

     n=Grid%n
     ALLOCATE(dum1(n),dum2(n),arg(n),der1(n),der2(n))

     norbit=Orbit%norbit
     ALLOCATE(M(norbit,norbit),MM(norbit-1,norbit-1),ux(norbit),vx(norbit),&
&             Rp(norbit-1))

     CALL SetIndex(Orbit)
     last=HSZ%zero_index


     dum1=0;  arg=0; M=0
     Do io=1,norbit
        call Calc_dexdphi_io(Grid,Orbit,io,der1)
        dum1=dum1+Orbit%occ(io)*der1*Orbit%wfn(:,io)
        der1=der1*Orbit%wfn(:,io)
        der1(1)=0.d0; der1(2:n)=der1(2:n)/Grid%r(2:n)
        ux(io)=integrator(Grid,der1)
        write(6,*) 'io const' , io,ux(io) ; call flush(6)
        do jo=1,norbit
           arg=0
           arg(2:n)=(Orbit%wfn(2:n,io)*Orbit%wfn(2:n,jo))**2/Orbit%den(2:n)
           M(io,jo)=integrator(Grid,arg)
        enddo
     Enddo

     dum1(1)=0; dum2=0.d0
     do io=1,norbit
        dum2=dum2-Orbit%occ(io)*geqn(:,io)*Orbit%wfn(:,io)
     enddo
  ! correct for incorrect asymptote:
  !  dum2=dum2-Orbit%occ(last)*HSZ%U(last)*Orbit%wfn(:,last)**2  ! already done

     do i=2,n
        if (Orbit%den(i)>machine_zero) then
           dum1(i)=dum1(i)/Orbit%den(i)/Grid%r(i)
           dum2(i)=dum2(i)/Orbit%den(i)
        else
          dum1(i)=0
        endif
     enddo
     call extrapolate(Grid,dum1)
     call extrapolate(Grid,dum2)

     MM=0
     i=0
     do io=1,norbit
        if(io/=last) then
          i=i+1;j=0
          MM(i,i)=1.d0
          arg=(dum1+dum2)*(Orbit%wfn(:,io)**2)
          Rp(i)=integrator(Grid,arg)
          do jo=1,norbit
             if(jo/=last)  then
               j=j+1
               Rp(i)=Rp(i)-M(io,jo)*Orbit%occ(jo)*ux(jo)
               MM(i,j)=MM(i,j)-M(io,jo)*Orbit%occ(jo)
             endif
          enddo
        endif
     enddo

     call linsol(MM,Rp,norbit-1,norbit,norbit,norbit-1)
     i=0;arg=0
     Do io=1,norbit
        if (io/=last) then
           i=i+1
           vx(io)=Rp(i)
        else
           vx(io)=ux(io)
        endif
        write(6,'(i5," ux ", 1p,1e15.7,"  vx  ",1p,1e15.7)') io,ux(io),vx(io)
        arg=arg+Orbit%occ(io)*(vx(io)-ux(io))*(Orbit%wfn(:,io)**2)
     enddo
     arg(1)=0
     do i=2,n
        if (Orbit%den(i)>machine_zero) then
           arg(i)=arg(i)/Orbit%den(i)
        endif
     enddo
     call extrapolate(Grid,arg)
     rvx=(dum1+dum2+arg)*Grid%r

    !call mkname(jcount,stuff)
    !OPEN(unit=1001,file='AVX'//TRIM(stuff),form='formatted')
    !DO i = 1,n
    !   WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
    !&       rvx0(i),rvx(i),dum1(i),dum2(i),arg(i),Orbit%den(i)
    !ENDDO
    !CLOSE(1001)
    jcount=jcount+1

    !if (FULL) stop
     DEALLOCATE(dum1,dum2,arg,der1,der2,M,MM,ux,vx,Rp)

  END SUBROUTINE NEWAdjustVX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Program to fix origin of radial function affected by
!     uncontrolled derivatives at origin
!     Fit to polynomial of order power
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE fixorigin(Grid,fin,gin,fout,gout)
       Type(GridInfo), INTENT(IN) :: Grid
       Real(8), INTENT(IN) :: fin(:),gin(:)
       REAL(8), INTENT(INOUT) :: fout(:),gout(:)

       INTEGER, parameter ::  nmin=10,orders=3,dimm=5
       REAL(8), parameter ::  rmin=0.001d0
       REAL(8), ALLOCATABLE :: M(:,:),X(:)
       INTEGER :: i,j,k,n,n0,n1,s1,s2,s3
       REAL(8) :: rm,xx

       fout=fin
       gout=gin
    !   check to see if correction is needed
       if (SUM(ABS(gout(1:5)))/SUM(ABS(gout(5:10)))>1.d0) then

           n=Grid%n
           n0=5
           n1=n0+nmin;j=0
           xx=10**orders
           do i=n0+1,n
              j=j+1
              if (Grid%r(i)/Grid%r(n0)>xx) exit
           enddo
           n1=MAX(n1,j)

           ALLOCATE(M(dimm,dimm),X(dimm))
           M=0;X=0

           rm=Grid%r(n1)
           write(6,*) 'fixorigin == n0,n1 = ', n0,n1,rm

           do i=1,dimm
              do k=n0,n1
                 X(i)=X(i)+fin(k)*((Grid%r(k)/rm)**(i-1))
              enddo
              do j=1,dimm
                 do k=n0,n1
                    M(i,j)=M(i,j)+((Grid%r(k)/rm)**(i+j-2))
                 enddo
              enddo
           enddo

           s1=size(M,1);s2=size(M,2);s3=size(X)
           call linsol(M,X,dimm,s1,s2,s3)
           write(6,*) 'Returning from linsol in fixorigin'
             do i=1,dimm
                write(6,*) i,X(i)
             enddo

           do k=1,n1
              fout(k)=0;gout(k)=0
              do i=1,dimm
                 fout(k)=fout(k)+X(i)*((Grid%r(k)/rm)**(i-1))
                 if (i>1) then
                    gout(k)=gout(k)+(i-1)*X(i)*((Grid%r(k)/rm)**(i-2))/rm
                 endif
              enddo
           enddo

           !do k=1,n1+20
           !   write(6,'(1p,20e15.7)') Grid%r(k),fin(k),fout(k),gin(k),gout(k)
           !enddo

           DEALLOCATE(M,X)
       ENDIF
   END SUBROUTINE fixorigin



  SUBROUTINE EXXdump(Grid,Orbit,ifo)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    INTEGER, INTENT(IN) :: ifo

    INTEGER :: i,io,jo,n,norbit

    n=Grid%n;norbit=Orbit%norbit

    WRITE(ifo,*) HSZ%zero_index,HSZ%lmax,HSZ%betaL,HSZ%grad2,HSZ%matchpoint
    WRITE(ifo,*)((HSZ%psi(i,io),i=1,n),io=1,norbit)
    WRITE(ifo,*)((HSZ%psiref(i,io),i=1,n),io=1,norbit)
    WRITE(ifo,*) (HSZ%U(io),(HSZ%LMBD(jo,io),jo=1,norbit),io=1,norbit)
    WRITE(ifo,*) (HSZ%Ucore(io),(HSZ%LMBDcore(jo,io),jo=1,norbit),io=1,norbit)
    WRITE(ifo,*) (HSZ%Uvale(io),(HSZ%LMBDvale(jo,io),jo=1,norbit),io=1,norbit)
    WRITE(ifo,*) (HSZ%Uref(io),io=1,norbit)
    WRITE(ifo,*) (HSZ%rVxref(i), &
&             HSZ%shift(i),HSZ%coreshift(i),i=1,n)

  END SUBROUTINE EXXdump

  SUBROUTINE EXXload(Grid,Orbit,ifo)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    INTEGER, INTENT(IN) :: ifo

    INTEGER :: i,io,jo,n,norbit

    n=Grid%n;norbit=Orbit%norbit

    READ(ifo,*) HSZ%zero_index,HSZ%lmax,HSZ%betaL,HSZ%grad2,HSZ%matchpoint
    ALLOCATE(HSZ%psi(n,norbit),HSZ%U(norbit),HSZ%LMBD(norbit,norbit),&
&        HSZ%shift(n),HSZ%rVxref(n),&
&        HSZ%psiref(n,norbit),HSZ%Ucore(norbit),HSZ%LMBDcore(norbit,norbit),&
&        HSZ%Uvale(norbit),HSZ%LMBDvale(norbit,norbit),HSZ%Uref(norbit),&
&        HSZ%coreshift(n))
    READ(ifo,*)((HSZ%psi(i,io),i=1,n),io=1,norbit)
    READ(ifo,*)((HSZ%psiref(i,io),i=1,n),io=1,norbit)
    READ(ifo,*) (HSZ%U(io),(HSZ%LMBD(jo,io),jo=1,norbit),io=1,norbit)
    READ(ifo,*) (HSZ%Ucore(io),(HSZ%LMBDcore(jo,io),jo=1,norbit),io=1,norbit)
    READ(ifo,*) (HSZ%Uvale(io),(HSZ%LMBDvale(jo,io),jo=1,norbit),io=1,norbit)
    READ(ifo,*) (HSZ%Uref(io),io=1,norbit)
    READ(ifo,*) (HSZ%rVxref(i),&
&             HSZ%shift(i),HSZ%coreshift(i),i=1,n)

  END SUBROUTINE EXXload

END MODULE exx_mod
