MODULE exx_mod
  USE anderson_driver
  USE atomdata
  USE exxdata
  USE fock
  USE general_mod
  USE globalmath
  USE gridmod
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
    INTEGER :: n,io,icount,ok,loop
    REAL(8) :: en1,x
    REAL(8), ALLOCATABLE :: arg(:),starg(:)
    REAL(8) :: errv0=1.d-7,errv
    LOGICAL :: success,done
    INTEGER :: nrestart=5,mxloop1=12,mxloop2=50,mxloop=200
    INTEGER :: firsttime=0

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
&           HSZ%LMBD(io,io),HSZ%rVxcore(n),HSZ%psicore(n,io),&
&           HSZ%psivale(n,io),HSZ%psiref(n,io),HSZ%rVxref(n),HSZ%rVxvale(n),&
&           HSZ%LMBDref(io,io),HSZ%LMBDcore(io,io),HSZ%LMBDvale(io,io),stat=ok)
       IF (ok /=0 ) THEN
          WRITE(6,*) 'Error in allocate HSZdata ', io,n,ok
          STOP
       ENDIF
       HSZ%psi=0;HSZ%psicore=0;HSZ%psivale=0
       HSZ%shift=0
       HSZ%grad=0
       HSZ%U=0; HSZ%Ucore=0;HSZ%Uvale=0
       HSZ%LMBD=0;HSZ%LMBDcore=0;HSZ%LMBDvale=0;HSZ%LMBDref=0
       HSZ%grad2=1.d13
       HSZ%rVxref =0
       HSZ%rVxcore=0
       HSZ%rVxvale=0
       HSZ%coreshift=0
       HSZ%matchpoint=Gridwk%n

       CALL Init_EXX_vx(Gridwk,Orbitwk,Potwk)

       Potwk%rv=Potwk%rvh+Potwk%rvx
       CALL zeropot(Gridwk,Potwk%rv,Potwk%v0,Potwk%v0p)
       Potwk%rv=Potwk%rv+Potwk%rvn

       arg=Potwk%rv

       CALL InitAnderson_dr(AC,6,5,n,0.5d0,1.d3,mxloop1,1.d-11,1.d-16,.TRUE.)
       CALL DoAndersonMix(AC,arg,en1,EXX1sub,success)
       SCFwk%iter=AC%CurIter
       SCFwk%delta=AC%res
       CALL FreeAnderson(AC)
       WRITE(6,*) 'Finished first step ', en1 ,' success = ', success
       WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta

       firsttime=1
       !CALL Report_EXX_functions('FI')

    ENDIF

    IF (setupfrozencore) THEN
       !HSZ%rVxref=Potwk%rvx    ! should have already been done
       CALL Init_EXX_core(Gridwk,Orbitwk,arg,done)
       If (.not.done) then
           CALL InitAnderson_dr(AC,6,5,n,0.005d0,1.d3,mxloop1,&
&                 1.d-11,1.d-16,.TRUE.)
           CALL DoAndersonMix(AC,arg,en1,VXcoreSub,success)
           WRITE(6,*) 'Finished setupfrozencore ',AC%CurIter,AC%res
           CALL FreeAnderson(AC)
           do icount=1,nrestart/2
              call AdjustVcore(Gridwk,arg)
              CALL InitAnderson_dr(AC,6,5,n,0.005d0,1.d3,mxloop1,&
&                 1.d-11,1.d-16,.TRUE.)
              CALL DoAndersonMix(AC,arg,en1,VXcoreSub,success)
              WRITE(6,*) 'Finished setupfrozencore # ',icount ,AC%CurIter,AC%res
              CALL FreeAnderson(AC)
           enddo
           starg=arg; HSZ%rVxcore=arg
   !  check result by running frozencore calculation -- arg now Vxvale
           CALL Init_EXX_vx(Gridwk,Orbitwk,Potwk)
           arg=Potwk%rvx
           loop=mxloop1; errv=errv0;  zeroadjust=.TRUE.
           CALL SetMatchPoint(Gridwk,Orbitwk,HSZ%zero_index)
           CALL InitAnderson_dr(AC,6,5,n,0.005d0,1.d3,loop,errv,1.d-16,.TRUE.)
           CALL DoAndersonMix(AC,arg,en1,EXX3Sub,success)
           SCFwk%iter=AC%CurIter
           SCFwk%delta=AC%res
           CALL FreeAnderson(AC)
           WRITE(6,*) 'Finished Anderson Mix FC', en1 ,' success = ', success
           WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta
           x=MAXVAL(ABS(starg(1:HSZ%corerange)-(HSZ%rVxref(1:HSZ%corerange)-&
&                    arg(1:HSZ%corerange))))
           write(6,*) 'frozencore and set core max diff ' , x
           if (ABS(x)<0.001d0) then
              write(6,*) 'Frozencore consistent with setcore '
              arg=HSZ%rVxref-arg
           else
              starg=HSZ%rVxref-arg; HSZ%rVxcore=starg
              DO icount=1,nrestart
                  loop=mxloop1; errv=errv0;   zeroadjust=.TRUE.
                  CALL InitAnderson_dr(AC,6,5,n,0.005d0,1.d3,loop,&
&                    errv,1.d-16,.TRUE.)
                  CALL AdjustV(Gridwk,arg,icount)
                  CALL DoAndersonMix(AC,arg,en1,EXX3Sub,success)
                  SCFwk%iter=AC%CurIter
                  SCFwk%delta=AC%res
                  CALL FreeAnderson(AC)
                  WRITE(6,*) 'Finished Anderson Mix ',&
&                    icount,en1 ,' success = ', success
                  WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta
                  x=MAXVAL(ABS(starg(1:HSZ%corerange)-&
&                   (HSZ%rVxref(1:HSZ%corerange)-arg(1:HSZ%corerange))))
                  write(6,*) 'frozencore and set core max diff ' , x
                  if (ABS(x)<0.001d0) exit
                  starg=HSZ%rVxref-arg
              ENDDO
              arg=HSZ%rVxref-arg
           endif
        endif
       HSZ%rVxcore=arg
       HSZ%rVxvale=HSZ%rVxref-HSZ%rVxcore
       CALL Calc_dedv_v(Gridwk,Potwk,Orbitwk,HSZ%rVxvale)
       CALL Report_EXX_functions(scftype)
       CALL Reportandset_coreEXX_functions
       CALL Get_FCKinCoul(Gridwk,Potwk,Orbitwk,FCwk,SCFwk)
       CALL Get_FCEnergy_EXX(Gridwk,Orbitwk,FCwk,SCFwk)
       DEALLOCATE(arg,starg)
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
    loop=mxloop1; errv=errv0;  zeroadjust=.FALSE.
    CALL InitAnderson_dr(AC,6,5,n,0.1d0,1.d3,loop,errv,1.d-16,.TRUE.)
    !CALL InitAnderson_dr(AC,6,0,n,0.1d0,1.d3,loop,errv,1.d-16,.TRUE.)
    CALL DoAndersonMix(AC,arg,en1,EXX3Sub,success)
    SCFwk%iter=AC%CurIter
    SCFwk%delta=AC%res
    CALL FreeAnderson(AC)
    WRITE(6,*) 'Finished Anderson Mix second', en1 ,' success = ', success
    WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta

    DO icount=1,nrestart
       IF (frozencorecalculation) then
          x=ABS(HSZ%Uvale(HSZ%zero_index))
       ELSE
          x=ABS(HSZ%U(HSZ%zero_index))
       ENDIF
       IF (x<1.d-6) EXIT
       loop=mxloop1; errv=errv0;   zeroadjust=.TRUE.
       CALL InitAnderson_dr(AC,6,5,n,0.005d0,1.d3,loop,errv,1.d-16,.TRUE.)
       CALL AdjustV(Gridwk,arg,icount)
       CALL DoAndersonMix(AC,arg,en1,EXX3Sub,success)
       SCFwk%iter=AC%CurIter
       SCFwk%delta=AC%res
       CALL FreeAnderson(AC)
       WRITE(6,*) 'Finished Anderson Mix ', icount,en1 ,' success = ', success
       WRITE(6,*) 'Convergence', SCFwk%iter,SCFwk%delta
    ENDDO
       If (frozencorecalculation) then
           CALL Get_FCKinCoul(Gridwk,Potwk,Orbitwk,FCwk,SCFwk)
           CALL Get_FCEnergy_EXX(Gridwk,Orbitwk,FCwk,SCFwk)
       Endif

    CALL Report_EXX_functions(scftype)

    DEALLOCATE(arg,starg)
  END SUBROUTINE EXX_SCF


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

    !WRITE(6,*) 'Returning from SetIndex ', HSZ%zero_index,HSZ%lmax
    IF (frozencorecalculation) THEN
       DO io=1,Orbit%norbit
          IF (Orbit%iscore(io).AND.outer==io) THEN
             WRITE(6,*) 'Error in FC SetIndex -- core state ', io,outer
             STOP
          ENDIF
       ENDDO
    ENDIF


  END SUBROUTINE SetIndex

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
!!!!  AdjustV   -- simplified version
!!!!   Return rvx = (w+cons*r)*theta+rvx0*(1-theta)
!!!!     theta=tanh((r/rc)^2)
!!!!   So that outer orbital constraint is satisfied
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE AdjustV(Grid,w,icase)
    TYPE(GridInfo) :: Grid
    REAL(8), INTENT(INOUT) :: w(:)
    INTEGER, INTENT(IN) :: icase

    INTEGER :: i,j,k,n,io,last,l,point,begin
    REAL(8),ALLOCATABLE :: rvx(:), rvx0(:),f(:,:),arg(:),Coeff(:),m(:,:)
    REAL(8),ALLOCATABLE ::arg_cs(:)
    REAL(8), POINTER :: r(:),psi(:,:),U(:)
    REAL(8) :: energy,val,x,Numerator,beta,rc
    REAL(8), PARAMETER :: smll=1.d-9
    INTEGER, PARAMETER :: many=5
    INTEGER, SAVE :: matchpt
    LOGICAL :: success
    INTEGER :: counter=0
    CHARACTER(4) :: stuff

    n=Grid%n
    r=>Grid%r
    ALLOCATE(rvx0(n),rvx(n),f(n,many),arg(n),Coeff(many),m(many,many),&
&        arg_cs(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in AdjustV allocation' ,n,i
       STOP
    ENDIF

    if (frozencorecalculation) then
       psi=>HSZ%psivale
       U=>HSZ%Uvale
    else
       psi=>HSZ%psi
       U=>HSZ%U
    endif
    WRITE(6,*) 'In AdjustV '; CALL flush(6)
    IF (.NOT.setupfrozencore) CALL Find_Hartree(w)

    rvx=0
    CALL SetIndex(Orbitwk)
    last=HSZ%zero_index

    CALL ApproxVx(Gridwk,Orbitwk,last,rvx0)
    !arg=Orbitwk%wfn(:,last)*psi(:,last)
    !j=HSZ%matchpoint
    !DO i=HSZ%matchpoint,FindGridindex(Gridwk,10.d0)
    !   IF (arg(i)*arg(i-1)<0.d0) THEN
    !      j=i+2
    !      EXIT
    !   ENDIF
    !ENDDO
    !WRITE(6,*) 'AdjustV: ', j,Gridwk%r(j),HSZ%matchpoint,Gridwk%r(HSZ%matchpoint)
    !HSZ%matchpoint=MAX(j,HSZ%matchpoint)
    !IF (icase>1.AND.matchpt>HSZ%matchpoint) THEN
    !   point=matchpt+5
    !   WRITE(6,'( "AdjustV point ",i5,1p,1e15.7,2x,i5,1p,1e15.7)')&
    !&       point,Gridwk%r(point)
    !ELSE
    !   beta=SQRT(ABS(Orbitwk%eig(last)))
    !   arg=w+U(last)*Gridwk%r-rvx0
    !   i=HSZ%matchpoint
    !   j=FindGridindex(Gridwk,3.d0/beta)
    !   point=MINLOC(ABS(arg(i:j)),1)+i-1
    !   WRITE(6,'( "AdjustV point ",i5,1p,1e15.7,2x,i5,1p,1e15.7)')&
    !&       point,Gridwk%r(point),i,Gridwk%r(i),j,Gridwk%r(j)
    !ENDIF
    write(6,*) 'In adjustv before beta ',Orbitwk%eig(last); call flush(6)
    f=0
    f(:,1)=Orbitwk%wfn(:,last)
    if (f(n,1)<0.d0) f(:,1)=-f(:,1)
    do i=n,2,-1
       j=i
       if (f(i-1,1)<f(i,1)) exit
    enddo
    point=j
    write(6,*) 'Simple AdjustV ', point,Gridwk%r(point) ; call flush(6)
    f(1:j-1,1)=1.d0
    f(j:n,1)=f(j:n,1)/f(j,1)
    f(:,2)=1.d0-f(:,1)
    HSZ%matchpoint=point
    matchpt=point
    rc=Gridwk%r(point)
    f(:,3)=f(:,1)*(Orbitwk%wfn(:,last)**2)
    f(:,4)=(w(:)*f(:,1)+rvx0(:)*f(:,2))*(Orbitwk%wfn(:,last)**2)

    IF (frozencorecalculation) THEN
       CALL Calc_dexdphi_io_v(Gridwk,Orbitwk,last,arg)
    ELSE
       CALL Calc_dexdphi_io(Gridwk,Orbitwk,last,arg)
       ! if Colle Salvetti
       IF(ColleSalvetti==.TRUE.) THEN
          CALL Calc_decdchi_io(Orbitwk,last,arg_cs)
          arg = arg + arg_cs
       ENDIF

    ENDIF

    arg=arg*Orbitwk%wfn(:,last)-f(:,4); f(:,5)=arg
    arg(1)=0; arg(2:n)=arg(2:n)/Gridwk%r(2:n)
    val=integrator(Gridwk,arg)/integrator(Gridwk,f(:,3))

    WRITE(6,*) 'AdjustV -- ', val

    rvx=0
    rvx=(w(:)+val*Gridwk%r(:))*f(:,1)+rvx0(:)*(f(:,2))

    write(6,*) 'AdjustV diff:',SUM(ABS(rvx-w))

    CALL mkname(counter,stuff)
    OPEN (unit=1001,file='adjust'//TRIM(stuff),form='formatted')
    DO i=1,n
       WRITE(1001,'(1p,15E15.7)') r(i),w(i),rvx(i),rvx0(i),&
&           Orbitwk%den(i),Orbitwk%wfn(i,last),f(i,1:5)
    ENDDO
    CLOSE(1001)
    w=rvx

    counter=counter+1
    DEALLOCATE(rvx0,rvx,f,arg,Coeff,m)

  END SUBROUTINE AdjustV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  AdjustV
!!!!   Return rvx = rvx0 + SUM_i Coef_i*f_i
!!!!   So that outer orbital constraint is satisfied
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE AdjustV_old(Grid,w,icase)
    TYPE(GridInfo) :: Grid
    REAL(8), INTENT(INOUT) :: w(:)
    INTEGER, INTENT(IN) :: icase

    INTEGER :: i,j,k,n,io,last,l,point,begin
    REAL(8),ALLOCATABLE :: rvx(:), rvx0(:),f(:,:),arg(:),Coeff(:),m(:,:)
    REAL(8),ALLOCATABLE ::arg_cs(:)
    REAL(8), POINTER :: r(:),psi(:,:),U(:)
    REAL(8) :: energy,val,x,Numerator,beta
    REAL(8), PARAMETER :: smll=1.d-9
    INTEGER, PARAMETER :: many=5
    INTEGER, SAVE :: matchpt
    LOGICAL :: success
    INTEGER :: counter=0
    CHARACTER(4) :: stuff

    n=Grid%n
    r=>Grid%r
    ALLOCATE(rvx0(n),rvx(n),f(n,many),arg(n),Coeff(many),m(many,many),&
&        arg_cs(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in AdjustV allocation' ,n,i
       STOP
    ENDIF

    if (frozencorecalculation) then
       psi=>HSZ%psivale
       U=>HSZ%Uvale
    else
       psi=>HSZ%psi
       U=>HSZ%U
    endif
    !WRITE(6,*) 'In AdjustV '; CALL flush(6)
    IF (.NOT.setupfrozencore) CALL Find_Hartree(w)

    rvx=0
    CALL SetIndex(Orbitwk)
    last=HSZ%zero_index

    CALL ApproxVx(Gridwk,Orbitwk,last,rvx0)
    arg=Orbitwk%wfn(:,last)*psi(:,last)
    j=HSZ%matchpoint
    DO i=HSZ%matchpoint,FindGridindex(Gridwk,10.d0)
       IF (arg(i)*arg(i-1)<0.d0) THEN
          j=i+2
          EXIT
       ENDIF
    ENDDO
    WRITE(6,*) 'AdjustV: ', j,Gridwk%r(j),HSZ%matchpoint,Gridwk%r(HSZ%matchpoint)
    HSZ%matchpoint=MAX(j,HSZ%matchpoint)
    IF (icase>1.AND.matchpt>HSZ%matchpoint) THEN
       point=matchpt+5
       WRITE(6,'( "AdjustV point ",i5,1p,1e15.7,2x,i5,1p,1e15.7)')&
&           point,Gridwk%r(point)
    ELSE
       beta=SQRT(ABS(Orbitwk%eig(last)))
       arg=w+U(last)*Gridwk%r-rvx0
       i=HSZ%matchpoint
       j=FindGridindex(Gridwk,3.d0/beta)
       point=MINLOC(ABS(arg(i:j)),1)+i-1
       WRITE(6,'( "AdjustV point ",i5,1p,1e15.7,2x,i5,1p,1e15.7)')&
&           point,Gridwk%r(point),i,Gridwk%r(i),j,Gridwk%r(j)
    ENDIF
    HSZ%matchpoint=point
    matchpt=point
    begin=point-many/2-1
    f=0.d0
    f(begin:n,2)=(rvx0(begin:n)+2.d0)/r(begin:n)
    f(begin:n,2)=f(begin:n,2)/f(begin,2)
    IF(many>2) THEN
       DO i=3,many
          f(begin:n,i)=f(begin:n,i-1)/r(begin:n)
          f(begin:n,i)=f(begin:n,i)/f(begin,i)
       ENDDO
    ENDIF

    m=0; Coeff=0

    m(1:many-1,1)=-r(begin:begin+many-2)    ! constant term * r
    DO i=0,many-2
       DO j=2,many
          m(i+1,j)=f(begin+i,j)
       ENDDO
       Coeff(i+1)=w(begin+i)-rvx0(begin+i)
    ENDDO

    IF (frozencorecalculation) THEN
       CALL Calc_dexdphi_io_v(Gridwk,Orbitwk,last,arg)
    ELSE
       CALL Calc_dexdphi_io(Gridwk,Orbitwk,last,arg)
       ! if Colle Salvetti
       IF(ColleSalvetti==.TRUE.) THEN
          CALL Calc_decdchi_io(Orbitwk,last,arg_cs)
          arg = arg + arg_cs
       ENDIF

    ENDIF
    arg=arg*Orbitwk%wfn(:,last); arg(1)=0; arg(2:n)=arg(2:n)/r(2:n)
    Numerator=integrator(Grid,arg)
    arg(1:point)=w(1:point)*(Orbitwk%wfn(1:point,last)**2)
    arg(1)=0; arg(2:point)=arg(2:point)/r(2:point)
    Numerator=Numerator-integrator(Grid,arg,1,point)
    arg(point:n)=rvx0(point:n)*(Orbitwk%wfn(point:n,last)**2)
    arg(point:n)=arg(point:n)/r(point:n)
    Numerator=Numerator-integrator(Grid,arg,point,n)
    Coeff(many)=Numerator

    arg(1:point)=Orbitwk%wfn(1:point,last)**2
    m(many,1)=integrator(Grid,arg,1,point)
    DO j=2,many
       arg=0
       arg(point:n)=(f(point:n,j)/r(point:n))*(Orbitwk%wfn(point:n,last)**2)
       m(many,j)=integrator(Grid,arg,point,n)
    ENDDO

    CALL SolveAXeqB(many,m,Coeff,1.d3)
    WRITE(6,*) 'AdjustV -- ', Coeff

    rvx=0
    rvx(1:point)=w(1:point)+Coeff(1)*r(1:point)
    rvx(point+1:n)=rvx0(point+1:n)
    DO j=2,many
       rvx(point+1:n)=rvx(point+1:n)+Coeff(j)*f(point+1:n,j)
    ENDDO

    !write(6,*) 'AdjustV diff:',SUM(ABS(rvx-w))

    CALL mkname(counter,stuff)
    OPEN (unit=1001,file='adjust'//TRIM(stuff),form='formatted')
    DO i=1,n
       WRITE(1001,'(1p,15E15.7)') r(i),w(i),rvx(i),rvx0(i),&
&           Orbitwk%den(i),Orbitwk%wfn(i,last),f(i,2:many)
    ENDDO
    CLOSE(1001)
    w=rvx

    counter=counter+1
    DEALLOCATE(rvx0,rvx,f,arg,Coeff,m)

  END SUBROUTINE AdjustV_old


!!!!  AdjustVcore
!!!!   For input w == rVxcore, return adjusted function
!!!!   So that outer orbital constraint is satisfied
!!!!      <ux|phi(outer)>=<phi|w|phi>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE AdjustVcore(Grid,w)
    TYPE(GridInfo) :: Grid
    REAL(8), INTENT(INOUT) :: w(:)

    INTEGER :: i,j,k,n,io,last
    REAL(8),ALLOCATABLE :: rvx(:), rvx0(:),dum(:),ux(:)
    REAL(8), POINTER :: r(:),phi(:)
    REAL(8) :: energy,val,x,fix,beta
    REAL(8), PARAMETER :: smll=1.d-9

    n=Grid%n
    r=>Grid%r
    ALLOCATE(rvx0(n),rvx(n),dum(n),ux(n))

    rvx=0
    CALL SetIndex(Orbitwk)
    last=HSZ%zero_index
    CALL Calc_dexdphi_io_c(Gridwk,Orbitwk,last,ux)
    phi=>Orbitwk%wfn(:,last)
    beta=SQRT(ABS(Orbitwk%eig(last)))
    dum=0; dum(2:n)=ux(2:n)*phi(2:n)/r(2:n)
    val=integrator(Grid,dum)
    dum=0; dum(2:n)=w(2:n)*(phi(2:n)**2)/r(2:n)
    x=integrator(Grid,dum)
    write(6,*) 'AdjustVcore : ', val,x, val-x
    if (ABS(val-x) >smll) then
        rvx0=exp(-2*beta*r); dum=rvx0*(phi**2)
        fix=(val-x)/integrator(Grid,dum)
        w=w+fix*r*rvx0
        dum=0; dum(2:n)=w(2:n)*(phi(2:n)**2)/r(2:n)
        x=integrator(Grid,dum)
        write(6,*) 'AdjustVcore -- fix: ', fix,val,x
    endif
    DEALLOCATE(rvx0,rvx,dum,ux)

  END SUBROUTINE AdjustVcore

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
          IF(Orbitwk%iscore(io)==.TRUE.) THEN
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

    IF(ColleSalvetti==.TRUE.) THEN
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
       IF(ColleSalvetti==.TRUE.) THEN
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
    IF (frozencorecalculation) rv=rv+HSZ%rVxcore

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
          IF(Orbitwk%iscore(io)==.TRUE.) THEN
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
    CALL DestroyPot(tmpPot)
    CALL DestroyOrbit(tmpOrbit)
    DEALLOCATE (rv)

  END SUBROUTINE Hsub


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  EXX3Sub									!!!!
!!!!   w is rvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXX3Sub(w,energy,grad,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io,l,iter,last
    REAL(8) :: x
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    INTEGER:: fcount=0
    TYPE (OrbitInfo), TARGET :: eigOrbit
    REAL(8), PARAMETER :: tol=1.d-10
    CHARACTER(4) :: stuff

    success=.TRUE.
    n=Gridwk%n

    HSZ%matchpoint=n
    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in EXX3sub allocation' ,n,i
       STOP
    ENDIF

    Potwk%rvx=w
    dum=w

    IF (zeroadjust.AND.AC%CurIter>10) THEN
       !  write(6,*) 'CALLING ADJUST V'; call flush(6)
       !CALL AdjustV(Gridwk,w)
       CALL Find_Hartree(w)
    ELSE
       CALL Find_Hartree(w)
    ENDIF


    rv=Potwk%rv

    CALL One_electron_energy_Report(Orbitwk,6)

    CALL   CopyOrbit(Orbitwk,EigOrbit)
    IF (frozencorecalculation) CALL Updatewfn(Gridwk,Potwk,EigOrbit,rv,success)

    CALL SetIndex(Orbitwk)
    CALL Calc_dedv(Gridwk,Potwk,Orbitwk,EigOrbit)

    CALL Get_Energy_EXX(Gridwk,Orbitwk,SCFwk%eexc)
    CALL SetIndex(Orbitwk) ! find index of most extended wfn
    SCFwk%etot=SCFwk%ekin+SCFwk%estatic+SCFwk%eexc
    IF(ColleSalvetti==.TRUE.) THEN
       CALL Get_Energy_CS(Gridwk,Orbitwk,SCFwk%oepcs)
       WRITE(6,*) 'Colle :', SCFwk%oepcs
       SCFwk%etot=SCFwk%etot+SCFwk%oepcs
    ENDIF

    IF (frozencorecalculation) THEN
       CALL Get_FCKinCoul(Gridwk,Potwk,Orbitwk,FCwk,SCFwk)
       CALL Get_FCEnergy_EXX(Gridwk,Orbitwk,FCwk,SCFwk)
       energy=SCFwk%evale
       CALL Total_FCEnergy_Report(SCFwk,6)
    ELSE
       energy=SCFwk%etot
       CALL Total_Energy_Report(SCFwk,6)
    ENDIF
    last=HSZ%zero_index
    energy=energy
    grad(2:n)=-HSZ%shift(2:n)/Gridwk%r(2:n)
    CALL extrapolate(Gridwk,grad)
    HSZ%grad2=DOT_PRODUCT(grad,grad)
    HSZ%grad=grad
    k=HSZ%matchpoint
    err=DOT_PRODUCT(grad(1:k),grad(1:k))

    IF (verboseoutput) CALL writestuff

    WRITE(6,'("ASiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err,HSZ%grad2
    WRITE(6,*) 'zero_index', fcount,HSZ%zero_index,HSZ%U(HSZ%zero_index)

    !CALL mkname(fcount,stuff)
    !OPEN (unit=1001,file='junk'//TRIM(stuff),form='formatted')
    !do i=1,n
    !   write(1001,'(1p,15E15.7)') Gridwk%r(i),w(i),dum(i)
    !enddo
    !close(1001)


    fcount=fcount+1
    CALL DestroyOrbit(eigOrbit)
    DEALLOCATE (dum,res,rv)

  END SUBROUTINE  EXX3Sub

  SUBROUTINE Calc_g3(Grid,l,eig,rv,rhs,tol,nn,MN)
    TYPE(GridInfo), INTENT(IN):: Grid
    INTEGER, INTENT(IN) :: l,nn
    REAL(8), INTENT(IN) :: eig
    REAL(8), INTENT(IN) :: rv(:)
    REAL(8), INTENT(IN) :: rhs(:)
    REAL(8), INTENT(IN) :: tol
    REAL(8), INTENT(OUT) :: MN(:,:)

    INTEGER :: i,j,k,n,ok,LWORK
    REAL(8) :: scale
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),d(:)
    REAL(8), ALLOCATABLE :: u(:,:),vt(:,:),s(:),work(:)
    INTEGER, ALLOCATABLE :: iwork(:)

    n=Grid%n;   LWORK=8*nn*nn
    ALLOCATE (a(n),b(n),c(n),d(n),u(nn,nn),vt(nn,nn),s(nn),&
&        work(LWORK),IWORK(8*nn) )

    CALL inhomo_numerov_coeff(Grid,l,nn+1,eig,rv,rhs,a,b,c,d)

    MN=0
    DO i=1,nn
       MN(i,i)=a(i+1)
    ENDDO
    DO i=1,nn-1
       MN(i,i+1)=b(i+2)
    ENDDO
    DO i=2,nn
       MN(i,i-1)=b(i)
    ENDDO

    CALL  dgesdd('A',nn,nn,MN,nn,s,u,nn,vt,nn,work,LWORK,IWORK,ok)
    WRITE(6,*) 'dgesdd completed for MN with info = ', ok ; call flush(6)

    MN=0; scale=s(1)*tol
    DO i=1,nn
          WRITE(6,*) 'g3 res:  i s = ', i,s(i) ; call flush(6)
       IF (s(i)>scale) THEN
          DO j=1,nn
             DO k=1,nn
                MN(j,k)=MN(j,k) + vt(i,j)*u(k,i)/s(i)
             ENDDO
          ENDDO
       ELSE
          WRITE(6,*) 'i s = ', i,s(i); call flush(6)
       ENDIF
    ENDDO

    u=MN
    vt=0
    DO i=1,nn
       vt(i,i)=c(i+1)
    ENDDO
    DO i=1,nn-1
       vt(i,i+1)=d(i+2)
    ENDDO
    DO i=2,nn
       vt(i,i-1)=d(i)
    ENDDO
    MN=0
    DO i=1,nn
       DO j=1,nn
          DO k=1,nn
             MN(i,j)=MN(i,j)+u(i,k)*vt(k,j)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE (a,b,c,d,u,vt,s,work,IWORK )
  END SUBROUTINE Calc_g3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    For exponentially decreasing function f(i), find index
  !       such that f(i)<val for i>index
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GetRange(f,val,index)
    REAL(8), INTENT(IN) :: f(:),val
    INTEGER, INTENT(OUT) :: index

    INTEGER :: n,i

    n=SIZE(f)
    index=n
    DO i=1,n
       IF (ABS(f(n-i+1))>ABS(val))  EXIT
       index=n-i+1
    ENDDO

    IF (index==1.OR.index==n) THEN
       WRITE(6,*) 'Strange result in GetRange', n,index,val,f(1),f(n)
       STOP
    ENDIF

  END SUBROUTINE GetRange

  SUBROUTINE Calc_Vxcore(Grid,Pot,Orbit,Vxcore)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    REAL(8), INTENT(INOUT) :: Vxcore(:)

    INTEGER :: i,j,k,n,ic,io,rangeindex,norbit,ncore,nn
    REAL(8) :: occ
    REAL(8), PARAMETER :: smallwfn=1.d-12,threshold=1.d-8,tol=1.d-8
    REAL(8), ALLOCATABLE :: coreshift(:),accumshift(:),g1(:,:),rhs(:),vl(:)
    REAL(8), ALLOCATABLE :: MN(:,:),MNaccum(:,:), st(:,:)
    REAL(8), POINTER :: rv(:)

    n=Grid%n; norbit=Orbit%norbit

    Vxcore=0

    rangeindex=0; k=0;
    DO io=1,norbit
       IF(Orbit%iscore(io)) THEN
          CALL GetRange(Orbit%wfn(:,io),smallwfn,i)
          IF (rangeindex<i) rangeindex=i
          k=k+1
       ENDIF
    ENDDO
    ncore=k;nn=rangeindex-1

    IF (ncore<1) THEN
       WRITE(6,*) ' No core states found  -- Vxcore == 0'
       RETURN
    ENDIF

    WRITE(6,*) 'In Calc_Vxcore rangeindex = ', rangeindex,Grid%r(rangeindex)
    rangeindex=n
    WRITE(6,*) 'Reset rangeindex ', rangeindex

    ALLOCATE(coreshift(n), accumshift(n),g1(n,norbit),st(n,norbit),rhs(n),vl(n))
    ALLOCATE(MN(nn,nn),MNaccum(nn,nn))

    coreshift(:)=0
    DO io=1,norbit
       IF (Orbit%iscore(io)) THEN
          coreshift=coreshift+2*Orbit%occ(io)*Orbit%wfn(:,io)*HSZ%psi(:,io)
       ENDIF
    ENDDO

    accumshift=coreshift;   g1=0 ;    MNaccum=0
    rv=>Pot%rv
    DO io=1,norbit
       IF (.NOT.Orbit%iscore(io)) THEN
          occ=Orbit%occ(io)
          IF (occ>threshold) THEN
             rhs=0; vl=0
             CALL Calc_dexdphi_io_c(Grid,Orbit,io,vl)
             rhs(2:n)=vl(2:n)/Grid%r(2:n)
             st(:,io)=rhs(:)
             CALL inhomo_numerov_SVD(Grid,Orbit%l(io),rangeindex,&
&                 Orbit%eig(io),tol,rv,rhs,vl)
             g1(:,io)=vl
             accumshift=accumshift+occ*vl*Orbit%wfn(:,io)
             rhs(2:n)=Orbit%wfn(2:n,io)/Grid%r(2:n)
             CALL extrapolate(Grid,rhs)
             CALL Calc_g3(Grid,&
&                 Orbit%l(io),Orbit%eig(io),rv,rhs,tol,nn,MN)
             DO j=1,nn
                DO i=1,nn
                   MNaccum(i,j)=MNaccum(i,j)+occ*Orbit%wfn(i+1,io)*MN(i,j)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDDO

    Vxcore=0; Vxcore(2:rangeindex)=accumshift(2:rangeindex)
    CALL SolveAXeqB(nn,MNaccum,Vxcore(2:rangeindex),1.d0/tol)

    OPEN (unit=1001,file='extra',form='formatted')
    DO i=1,n
       WRITE(1001,'(1p,15E15.7)') Gridwk%r(i),coreshift(i), &
&           accumshift(i),Vxcore(i),&
&           (g1(i,j),j=1,norbit),(st(i,j),j=1,norbit)
    ENDDO
    CLOSE(1001)


    DEALLOCATE(coreshift,accumshift,g1,rhs,vl,MN,MNaccum)

  END SUBROUTINE Calc_Vxcore


  SUBROUTINE Init_EXX_core(Grid,Orbit,rvcore,done)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    REAL(8), INTENT(INOUT) :: rvcore(:)
    Logical, INTENT(OUT) :: done

    INTEGER :: i,j,k,n,ic,io,rangeindex,norbit,ncore,nn
    REAL(8) :: occ,accumocc,rm
    REAL(8), PARAMETER :: smallwfn=1.d-8,threshold=1.d-8,tol=1.d-8
    REAL(8), ALLOCATABLE :: rhs(:),vl(:),tmp(:,:)
    REAL(8), POINTER :: coreshift(:)

    done=.false.
    n=Grid%n; norbit=Orbit%norbit
    HSZ%Ucore=HSZ%U; HSZ%Uvale=0
    HSZ%LMBDcore=HSZ%LMBD; HSZ%LMBDvale=0
    rvcore=0

    rangeindex=0; k=0;
    DO io=1,norbit
       IF(Orbit%iscore(io)) THEN
          CALL GetRange(Orbit%wfn(:,io),smallwfn,i)
          IF (rangeindex<i) rangeindex=i
          k=k+1
       ENDIF
    ENDDO
    ncore=k;HSZ%corerange=rangeindex

    IF (ncore<1) THEN
       WRITE(6,*) ' No core states found  -- Vxcore == 0'
       RETURN
    ENDIF

    WRITE(6,*) 'In Calc_Vxcore rangeindex = ', rangeindex,Grid%r(rangeindex)

    ALLOCATE(rhs(n),vl(n),tmp(n,norbit))

    coreshift=>HSZ%coreshift
    coreshift(:)=0 ; HSZ%psicore=0; HSZ%psivale=0
    DO io=1,norbit
       IF (Orbit%iscore(io)) THEN
          coreshift=coreshift+2*Orbit%occ(io)*Orbit%wfn(:,io)*HSZ%psi(:,io)
          HSZ%psicore(:,io)=HSZ%psi(:,io)
       ENDIF
    ENDDO

    !   Check to see if this is a case of single valence shell
    ic=0
    Do io=1,norbit
       IF (.not.Orbit%iscore(io)) THEN
          IF (Orbit%occ(io)>threshold) THEN
             ic=ic+1
             if (ic==1) CALL InitialVx(Grid,Orbit%wfn(:,io),rvcore)
             HSZ%psicore(:,io)=HSZ%psi(:,io)
          ENDIF
       ENDIF
    Enddo

    If (ic==1) then
       rvcore=HSZ%rVxref-rvcore
       done=.true.
    else
       rvcore=0; accumocc=0; k=0
       DO io=1,norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             occ=Orbit%occ(io)
             IF (occ>threshold) THEN
                rhs=0; vl=0; accumocc=accumocc+occ
                CALL Calc_dexdphi_io_c(Grid,Orbit,io,vl)
                rvcore=rvcore+occ*vl
                k=k+1; tmp(:,k)=vl
             ENDIF
          ENDIF
       ENDDO

       rvcore=rvcore/accumocc


       accumocc=0
       do i=rangeindex+1,n
          accumocc=accumocc+ABS(rvcore(i))
          rvcore(i)=0
       enddo

    !  newreset
       rm=Grid%r(rangeindex)
       rvcore=0
       rvcore(2:rangeindex)=HSZ%rVxref(2:rangeindex)*&
&          (sin(Pi*Grid%r(2:rangeindex)/rm)/(Pi*Grid%r(2:rangeindex)/rm))**2

       call AdjustVcore(Grid,rvcore)

       write(6,*) 'Init_rvcore range, error ', Grid%r(rangeindex),accumocc
   endif
    OPEN (unit=1001,file='exSetCore',form='formatted')
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Grid%r(i),(tmp(i,j),j=1,k)
       ENDDO
    CLOSE(1001)


    DEALLOCATE(rhs,vl,tmp)

  END SUBROUTINE Init_EXX_core


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Calc_dedv_v
!!!     Solve  inhomogeneous equations for valence electrons only
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_dedv_v(Grid,Pot,Orbit,rvxvale)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(PotentialInfo), INTENT(IN):: Pot
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    REAL(8), INTENT(IN):: rvxvale(:)

    REAL(8), POINTER :: phi(:,:),phr(:,:),psi(:,:),r(:),rv(:),U(:),LMBD(:,:)
    REAL(8), ALLOCATABLE :: vl(:),dum(:),rhs(:),shift(:),vs(:),rvx(:),st(:),ss(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,io,jo,norbit,last
    INTEGER :: iter
    INTEGER :: counter=0
    LOGICAL :: ok
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,xx,q,ecoul,v00,energy
    REAL(8), PARAMETER :: threshold=1.d-8,tol=1.d-12

    n=Grid%n

    rv => Pot%rv

    n=Grid%n
    norbit=Orbit%norbit
    phi=>Orbit%wfn
    r=>Grid%r

    psi=>HSZ%psivale
    U=>HSZ%Uvale

    psi=0

    ALLOCATE(vl(n),dum(n),rhs(n),shift(n),vs(n),rvx(n),st(n),ss(n))

    rvx=rvxvale

    DO io=1,norbit
       occ=Orbit%occ(io)
       l=Orbit%l(io)
       energy=Orbit%eig(io)
       IF (occ>threshold) THEN
          IF(.not.Orbit%iscore(io)) THEN
             CALL Calc_dexdphi_io_v(Grid,Orbit,io,vl)
             rhs=0
             rhs(2:n)=(-rvx(2:n)*phi(2:n,io)+vl(2:n))/r(2:n)
             st(2:n)=(vl(2:n))/r(2:n)
             ss(2:n)=(rvx(2:n)*phi(2:n,io))/r(2:n)
             vl=rhs
             dum=vl*phi(:,io)
             U(io)=integrator(Gridwk,dum)
             rhs=rhs-U(io)*phi(:,io)
             CALL Calc_psi(Grid,&
&                 Orbitwk%l(io),Orbitwk%eig(io),rv,phi(:,io),rhs,dum)
             psi(:,io)=dum
             HSZ%psicore(:,io)=HSZ%psiref(:,io)-psi(:,io)
             HSZ%Ucore(io)=HSZ%Uref(io)-U(io)
             counter=counter+1
             vl=psi(:,io); dum=psi(:,io)
             call inhomo_numerov_SVD_bv(Gridwk,l,300,energy,tol,rv,rhs,vl)
             call inhomo_numerov_SVD_bv(Gridwk,l,400,energy,tol,rv,rhs,dum)
             do i=1,n
                write(200+counter,'(1p,10e15.7)') r(i),phi(i,io),rhs(i),&
&                         st(i),ss(i),psi(i,io),vl(i),dum(i)
             enddo
          ENDIF
      ENDIF
    ENDDO

    write(6,*) 'Returning from Calc_dedv with Uzero = ', U(HSZ%zero_index)

    DEALLOCATE(vl,dum,rhs,shift,vs,rvx,st,ss)

  END SUBROUTINE Calc_dedv_v

  SUBROUTINE Calc_dedv_c(Grid,Pot,Orbit,rvxcore)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    REAL(8), INTENT(IN) :: rvxcore(:)

    INTEGER :: i,j,k,n,ic,io,jo,rangeindex,norbit,ncore,nn
    REAL(8) :: occ
    REAL(8), PARAMETER :: smallwfn=1.d-8,threshold=1.d-8,tol=1.d-8
    REAL(8), ALLOCATABLE :: rhs(:),vl(:),dum(:)
    REAL(8), POINTER :: rv(:),phi(:,:)

    n=Grid%n; norbit=Orbit%norbit

    HSZ%shift=HSZ%coreshift


    ALLOCATE(rhs(n),vl(n),dum(n))

    rv=>Pot%rv;  phi=>Orbit%wfn
    DO io=1,norbit
       IF (.NOT.Orbit%iscore(io)) THEN
          occ=Orbit%occ(io)
          IF (occ>threshold) THEN
             rhs=0; vl=0;
             CALL Calc_dexdphi_io_c(Grid,Orbit,io,vl)
             rhs(2:n)=(-rvxcore(2:n)*phi(2:n,io)+vl(2:n))/Grid%r(2:n)
             vl=rhs*phi(:,io)
             HSZ%Ucore(io)=overlap(Grid,phi(:,io),rhs)
             HSZ%Uvale(io)=HSZ%Uref(io)-HSZ%Ucore(io)
             rhs=rhs-HSZ%Ucore(io)*phi(:,io)
             write(6,*) 'Ucore ',io,HSZ%Ucore(io)
             !DO jo=1,norbit
             !   IF (jo/=io) THEN
             !      IF (Orbit%l(jo)==Orbit%l(io)) THEN
             !         HSZ%LMBDcore(io,jo)=overlap(Grid,Orbit%wfn(:,jo),rhs)
             !         HSZ%LMBDvale(io,jo)=HSZ%LMBDref(io,jo)-HSZ%LMBDcore(io,jo)
             !         rhs=rhs-HSZ%LMBDcore(io,jo)*Orbit%wfn(:,jo)
             !         write(6,*) 'LMBD ', io,jo,HSZ%LMBDcore(io,jo)
             !      ENDIF
             !   ENDIF
             !ENDDO
             CALL Calc_psi(Grid,&
&                 Orbit%l(io),Orbit%eig(io),rv,phi(:,io),rhs,vl)
             HSZ%shift=HSZ%shift+2*occ*phi(:,io)*vl
             HSZ%psicore(:,io)=vl
             HSZ%psivale(:,io)=HSZ%psiref(:,io)-vl
          ENDIF
       ENDIF
    ENDDO

    write(6,*) 'Returning from Calc_dedv_c with Uzero = ',&
&            HSZ%Ucore(HSZ%zero_index)
   deallocate(vl,dum,rhs)
  END SUBROUTINE Calc_dedv_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  VXcoreSub
!!!!    Subroutine to find Vxcore after SCF calculation
!!!!   w is r*Vxcore
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE VXcoreSub(w,energy,grad,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io,l,iter,last
    REAL(8) :: x
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    INTEGER:: fcount=0
    REAL(8), PARAMETER :: tol=1.d-10
    CHARACTER(4) :: stuff

    success=.TRUE.
    n=Gridwk%n

    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in VXcoresub allocation' ,n,i
       STOP
    ENDIF

    CALL Calc_dedv_c(Gridwk,Potwk,Orbitwk,w)

    grad(2:n)=-HSZ%shift(2:n)/Gridwk%r(2:n)
    CALL extrapolate(Gridwk,grad)
    HSZ%grad2=DOT_PRODUCT(grad,grad)
    energy=HSZ%grad2
    HSZ%grad=grad
    k=HSZ%corerange
    err=DOT_PRODUCT(grad(1:k),grad(1:k))

    WRITE(6,'("Coreiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err,HSZ%grad2
    grad(k+1:n)=0    ! force core potential to be zero outside range

    !CALL mkname(fcount,stuff)
    !OPEN (unit=1001,file='junk'//TRIM(stuff),form='formatted')
    !do i=1,n
    !   write(1001,'(1p,25E15.7)') Gridwk%r(i),w(i),grad(i),&
    !&      (HSZ%psicore(i,io)*Orbitwk%wfn(i,io),io=1,Orbitwk%norbit)
    !enddo
    !close(1001)


    fcount=fcount+1
    DEALLOCATE (dum,res,rv)

  END SUBROUTINE VXcoreSub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  VXvaleSub
!!!!    Subroutine to find Vxvale after SCF calculation
!!!!   w is r*Vxvale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE VXvaleSub(w,energy,grad,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io,l,iter,last
    REAL(8) :: x
    REAL(8),ALLOCATABLE :: dum(:),res(:),rv(:)
    INTEGER:: fcount=0
    REAL(8), PARAMETER :: tol=1.d-10
    CHARACTER(4) :: stuff

    success=.TRUE.
    n=Gridwk%n

    HSZ%matchpoint=n
    ALLOCATE(dum(n),res(n),rv(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in VXvalesub allocation' ,n,i
       STOP
    ENDIF

    Potwk%rvx=w
    dum=w

    IF (zeroadjust.AND.AC%CurIter>10) THEN
       !  write(6,*) 'CALLING ADJUST V'; call flush(6)
       ! CALL AdjustV(Gridwk,w)
    ENDIF

    CALL SetIndex(Orbitwk)
    CALL Calc_dedv(Gridwk,Potwk,Orbitwk,Orbitwk)

    last=HSZ%zero_index
    grad(2:n)=-HSZ%shift(2:n)/Gridwk%r(2:n)
    CALL extrapolate(Gridwk,grad)
    HSZ%grad2=DOT_PRODUCT(grad,grad)
    energy=HSZ%grad2
    HSZ%grad=grad
    k=HSZ%matchpoint
    err=DOT_PRODUCT(grad(1:k),grad(1:k))

    WRITE(6,'("Valeiter",i5,1p,1e20.12,1p,2e15.7)')fcount,energy,err,HSZ%grad2
    WRITE(6,*) 'zero_index', fcount,HSZ%zero_index,HSZ%U(HSZ%zero_index)

    !CALL mkname(fcount,stuff)
    !OPEN (unit=1001,file='junk'//TRIM(stuff),form='formatted')
    !do i=1,n
    !   write(1001,'(1p,15E15.7)') Gridwk%r(i),w(i),grad(i)
    !enddo
    !close(1001)


    fcount=fcount+1
    DEALLOCATE (dum,res,rv)

  END SUBROUTINE VXvaleSub

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
    if (frozencorecalculation) then
       psi=>HSZ%psivale
       U=>HSZ%Uvale
       LMBD=>HSZ%LMBDvale
    else
       psi=>HSZ%psi
       U=>HSZ%U
       LMBD=>HSZ%LMBD
    endif

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
                IF(ColleSalvetti==.TRUE.) THEN
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
             rhs=rhs-term*phi(:,io)
             !DO jo=1,norbit   !   extra orbit orthogonalization
             !   IF (jo/=io) THEN
             !     IF (Orbitwk%l(jo)==Orbitwk%l(io)) THEN
             !        dum=vl*EigOrbitwk%wfn(:,jo)
             !        term=integrator(Grid,dum)
             !        rhs=rhs-term*EigOrbitwk%wfn(:,jo)
             !     ENDIF
             !  ENDIF
             !ENDDO
             geqn(:,io)=rhs    ! store for possible late use
             CALL Calc_psi(Grid,&
&                 Orbitwk%l(io),Orbitwk%eig(io),rv,phi(:,io),rhs,dum)
             shift=shift+2*occ*phi(:,io)*dum
             psi(:,io)=dum
          ENDIF
       ENDIF
    ENDDO

    !   check orthogonality   (should check for l values)
    !    do io=1,norbit
    !      do jo=1,norbit
    !         write(6,'("Check Ortho ", 2i5,1p,1e15.7)') &
    !&                io,jo,overlap(Grid,EigOrbitwk%wfn(:,jo),psi(:,io))
    !      enddo
    !   enddo

    HSZ%shift=shift
    if (setupfrozencore) then
      do io=1,norbit
         if (.not.Orbitwk%iscore(io)) then
            HSZ%psicore(:,io)=HSZ%psiref(:,io)-psi(:,io)
            HSZ%Ucore(io)=HSZ%Uref(io)-HSZ%Uvale(io)
            !do jo=1,norbit
            !   if (jo/=io) &
            !&   HSZ%LMBDcore(io,jo)=HSZ%LMBDref(io,jo)-HSZ%LMBDvale(io,jo)
            !enddo
         endif
      enddo
    else if (frozencorecalculation) then
      do io=1,norbit
         if (.not.Orbitwk%iscore(io)) then
            HSZ%psi(:,io)=HSZ%psicore(:,io)+HSZ%psi(:,io)
            HSZ%U(io)=HSZ%Ucore(io)+HSZ%Uvale(io)
            !do jo=1,norbit
            !   if (jo/=io) &
            !&   HSZ%LMBD(io,jo)=HSZ%LMBDcore(io,jo)+HSZ%LMBDvale(io,jo)
            !enddo
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
       IF (frozencorecalculation) THEN
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&              Potwk%rvh(i),Potwk%rvx(i),HSZ%rVxcore(i)
       ELSE
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&              Potwk%rvh(i),Potwk%rvx(i)
       ENDIF
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
       IF (frozencorecalculation) THEN
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&              Potwk%rvh(i),Potwk%rvx(i),HSZ%rVxcore(i)
       ELSE
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&              Potwk%rvh(i),Potwk%rvx(i)
       ENDIF
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
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

  SUBROUTINE Reportandset_coreEXX_functions
    INTEGER :: i,j,k,n,nmap,io,jo
    INTEGER, ALLOCATABLE :: mmap(:)

    ALLOCATE(mmap(Orbitwk%norbit))
       j=0
       DO i=1,Orbitwk%norbit
          IF (.NOT.Orbitwk%iscore(i)) THEN
             j=j+1
             mmap(j)=i
          ENDIF
       ENDDO
       nmap=j

    WRITE(6,*)
    WRITE(6,*) ' Summary of EXX valence orbitals for frozencore'
    WRITE(6,"(' n  l     occupancy            energy              U           UC     UV')")
    DO jo=1,nmap
       io=mmap(jo)
       WRITE(6,'(i2,1x,i2,4x,1p,5e15.7)') &
&           Orbitwk%np(io),Orbitwk%l(io),Orbitwk%occ(io),Orbitwk%eig(io),&
&                HSZ%Uref(io), HSZ%Ucore(io),HSZ%Uvale(io)
    ENDDO
    !WRITE(6,*) '   Summary of LMBD values '
    !DO jo=1,nmap
    !   io=mmap(jo)
    !   WRITE(6,*) ' io, l: ',io,Orbitwk%l(io)
    !        do j=1,Orbitwk%norbit
    !           if (j/=io.and.Orbitwk%l(io)==Orbitwk%l(j)) &
    !&              WRITE(6,'(2i4,3x,1p,3e15.7)') j,Orbitwk%l(j),&
    !&               HSZ%LMBDref(io,j),HSZ%LMBDcore(io,j),HSZ%LMBDvale(io,j)
    !        enddo
    !ENDDO

    OPEN (unit=1001,file='shiftSetCore',form='formatted')
    n=Gridwk%n
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),HSZ%shift(i),HSZ%coreshift(i),&
&              (HSZ%psicore(i,mmap(k))*Orbitwk%wfn(i,mmap(k)),k=1,nmap), &
&              (HSZ%psivale(i,mmap(k))*Orbitwk%wfn(i,mmap(k)),k=1,nmap)
       ENDDO
    CLOSE(1001)

    OPEN (unit=1001,file='potSetCore',form='formatted')
    DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&              Potwk%rvh(i),HSZ%rVxvale(i),HSZ%rVxcore(i)
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='psiSetCore',form='formatted')
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&              (HSZ%psicore(i,mmap(k)),k=1,nmap),&
&              (HSZ%psivale(i,mmap(k)),k=1,nmap)
       ENDDO
    CLOSE(1001)

    IF (ALLOCATED(mmap)) DEALLOCATE(mmap)

  END SUBROUTINE Reportandset_coreEXX_functions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Get_FCEnergy_EXX                                  !!!!
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
       WRITE(6,*) 'allocation error in Get_Energy_EXX', i,n
       STOP
    ENDIF

    dum=0; dum(2:n)=HSZ%rvxcore(2:n)*FC%valeden(2:n)/r(2:n)
    FCeex=integrator(Grid,dum)
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

  SUBROUTINE Calc_VXvale

    INTEGER :: i,j,k,n,l,nvale,io,ip,norbit,corematch
    INTEGER, ALLOCATABLE :: map(:)
    REAL(8) :: x,energy,zeroval,term,bp,term1,term2,term3,hfac
    REAL(8), PARAMETER :: threshold=1.d-8
    REAL(8), ALLOCATABLE :: RR(:,:),LM(:,:),RV(:),M(:,:),P(:,:)

    HSZ%rVxvale=0
    HSZ%rVxcore=HSZ%rVxref

    n=Gridwk%n; hfac=0.1d0*((Gridwk%h)**2)
    norbit=Orbitwk%norbit
    nvale=0
    DO io=1,norbit
       IF (.NOT.Orbitwk%iscore(io) .AND. Orbitwk%occ(io)>threshold) THEN
          nvale=nvale+1
       ENDIF
    ENDDO

    IF (nvale<=0) THEN
       WRITE(6,*) ' No valence orbitals found -- returning from Calc_VXvale'
       RETURN
    ELSE
       ALLOCATE(map(nvale),RR(n,nvale))
    ENDIF

    ip=0
    DO io=1,norbit
       IF (.NOT.Orbitwk%iscore(io) .AND. Orbitwk%occ(io)>threshold) THEN
          ip=ip+1
          map(ip)=io
       ENDIF
    ENDDO

    IF(nvale==1) THEN

       CALL  Calc_expot_io(Gridwk,Orbitwk,map(1),HSZ%rVxvale)
       HSZ%rVxcore=HSZ%rVxref-HSZ%rVxvale
       DEALLOCATE(map,RR)
       RETURN
    ENDIF

    ! Evaluate RR functions  !   not quite clear about 1/r
    RR=0
    DO ip=1,nvale
       CALL Calc_dexdphi_io_v(Gridwk,Orbitwk,map(ip),RR(:,ip))
       RR(1,ip)=0;   RR(2:n,ip)=RR(2:n,ip)/Gridwk%r(2:n)
       RR(:,ip)=-RR(:,ip)+HSZ%U(map(ip))*Orbitwk%wfn(:,map(ip))
    ENDDO

    OPEN(1001,file='RR', form='formatted')
    DO i=1,n
       WRITE(1001,'(1p,20E15.7)') Gridwk%r(i),(RR(i,ip),ip=1,nvale)
    ENDDO
    CLOSE(1001)


    ! Determine core match point
    x=1000000.d0
    DO io=1,norbit
       IF(Orbitwk%iscore(io)) THEN
          IF (ABS(Orbitwk%eig(io))<x) x=ABS(Orbitwk%eig(io))
       ENDIF
    ENDDO

    WRITE(6,*) 'Most extended core state is ', x
    x=18.42d0/(2*SQRT(x))
    WRITE(6,*) 'Estimated core range is ', x

    corematch=FindGridIndex(Gridwk,x)
    WRITE(6,*) 'corematch = ', corematch

    k=corematch      !  note that r(2) maps to index 1

    ALLOCATE(LM(k-1,k-1),RV(k-1),M(k-1,k-1),P(k-1,k-1))

    LM=0; RV=0
    DO ip=1,nvale
       io=map(ip)
       l=Orbitwk%l(io)
       energy=Orbitwk%eig(io)

       M=0; P=0
       CALL invert_numerov(Gridwk,l,k,energy,Potwk%rv,M,bp)
       WRITE(6,*) 'bp = ',bp
       DO i=1,k-1
          DO j=2,k-2
             P(i,j)=10*M(i,j)+M(i,j-1)+M(i,j+1)
          ENDDO
          P(i,1)=10*M(i,1)+M(i,2)
          P(i,k-1)=10*M(i,k-1)+M(i,k-2)
       ENDDO

       x=Orbitwk%occ(io); term1=-bp*HSZ%psi(k+1,io)/hfac
       term2=Orbitwk%wfn(k+1,io)*&
&           HSZ%rVXref(k+1)/Gridwk%r(k+1)
       term3=RR(k+1,ip)
       WRITE(6,*) 'term1,term2,term3 = ', term1,term2,term3
       DO i=1,k-1
          DO j=1,k-1
             LM(i,j)=LM(i,j)+x*Orbitwk%wfn(i+1,io)*Orbitwk%wfn(j+1,io)*P(i,j)
             RV(i)=RV(i)-x*Orbitwk%wfn(i+1,io)*RR(j+1,ip)*P(i,j)
          ENDDO
          RV(i)=RV(i)-x*Orbitwk%wfn(i+1,io)*M(i,k-1)*(term1+term2+term3)
       ENDDO
    ENDDO

    OPEN(1001,file='stuff',form='formatted')
    DO i=1,k-1
       WRITE(1001,'(1p,10E15.7)') Gridwk%r(i+1),RV(i),LM(i,i)
    ENDDO
    CLOSE(1001)

    CALL linsol(LM,RV,k-1)
    HSZ%rVXvale=HSZ%rVXref
    HSZ%rVXvale(2:k)=RV(1:k-1)*Gridwk%r(2:k)
    HSZ%rVXcore=HSZ%rVXref-HSZ%rVXvale

    DEALLOCATE(map,RR,LM,RV,M,P)

  END SUBROUTINE Calc_VXvale
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

  SUBROUTINE EXXdump(Grid,Orbit,ifo)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    INTEGER, INTENT(IN) :: ifo

    INTEGER :: i,io,jo,n,norbit

    n=Grid%n;norbit=Orbit%norbit

    WRITE(ifo,*) HSZ%zero_index,HSZ%lmax,HSZ%betaL,HSZ%grad2,HSZ%matchpoint,&
&       HSZ%corerange
    WRITE(ifo,*)((HSZ%psi(i,io),i=1,n),io=1,norbit)
    WRITE(ifo,*)((HSZ%psicore(i,io),i=1,n),io=1,norbit)
    WRITE(ifo,*)((HSZ%psivale(i,io),i=1,n),io=1,norbit)
    WRITE(ifo,*)((HSZ%psiref(i,io),i=1,n),io=1,norbit)
    WRITE(ifo,*) (HSZ%U(io),(HSZ%LMBD(jo,io),jo=1,norbit),io=1,norbit)
    WRITE(ifo,*) (HSZ%Ucore(io),(HSZ%LMBDcore(jo,io),jo=1,norbit),io=1,norbit)
    WRITE(ifo,*) (HSZ%Uvale(io),(HSZ%LMBDvale(jo,io),jo=1,norbit),io=1,norbit)
    WRITE(ifo,*) (HSZ%Uref(io),io=1,norbit)
    WRITE(ifo,*) (HSZ%rVxref(i),HSZ%rVxcore(i),HSZ%rVxvale(i),&
&             HSZ%shift(i),HSZ%coreshift(i),i=1,n)

  END SUBROUTINE EXXdump

  SUBROUTINE EXXload(Grid,Orbit,ifo)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    INTEGER, INTENT(IN) :: ifo

    INTEGER :: i,io,jo,n,norbit

    n=Grid%n;norbit=Orbit%norbit

    READ(ifo,*) HSZ%zero_index,HSZ%lmax,HSZ%betaL,HSZ%grad2,HSZ%matchpoint,&
&       HSZ%corerange
    ALLOCATE(HSZ%psi(n,norbit),HSZ%U(norbit),HSZ%LMBD(norbit,norbit),&
&        HSZ%shift(n),HSZ%rVxref(n),HSZ%rVxcore(n),HSZ%rVxvale(n),&
&        HSZ%psicore(n,norbit),HSZ%psivale(n,norbit),HSZ%psicore(n,norbit),&
&        HSZ%psiref(n,norbit),HSZ%Ucore(norbit),HSZ%LMBDcore(norbit,norbit),&
&        HSZ%Uvale(norbit),HSZ%LMBDvale(norbit,norbit),HSZ%Uref(norbit),&
&        HSZ%coreshift(n))
    READ(ifo,*)((HSZ%psi(i,io),i=1,n),io=1,norbit)
    READ(ifo,*)((HSZ%psicore(i,io),i=1,n),io=1,norbit)
    READ(ifo,*)((HSZ%psivale(i,io),i=1,n),io=1,norbit)
    READ(ifo,*)((HSZ%psiref(i,io),i=1,n),io=1,norbit)
    READ(ifo,*) (HSZ%U(io),(HSZ%LMBD(jo,io),jo=1,norbit),io=1,norbit)
    READ(ifo,*) (HSZ%Ucore(io),(HSZ%LMBDcore(jo,io),jo=1,norbit),io=1,norbit)
    READ(ifo,*) (HSZ%Uvale(io),(HSZ%LMBDvale(jo,io),jo=1,norbit),io=1,norbit)
    READ(ifo,*) (HSZ%Uref(io),io=1,norbit)
    READ(ifo,*) (HSZ%rVxref(i),HSZ%rVxcore(i),HSZ%rVxvale(i),&
&             HSZ%shift(i),HSZ%coreshift(i),i=1,n)

  END SUBROUTINE EXXload

END MODULE exx_mod
