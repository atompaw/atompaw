Module fock
  Use atomdata
  Use gridmod
  Use globalmath

  IMPLICIT NONE

 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Get_Energy_EXX					!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_Energy_EXX(Grid,Orbit,eex)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    REAL(8), INTENT(OUT) :: eex

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),dum(:),ddum(:)
    REAL(8), ALLOCATABLE :: Fnu(:,:),Lnu(:,:),Snu(:),Mnunup(:,:),Cnu(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,io,jo,ok,norbit,last,nodes
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,kappa
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=orbit%norbit
    phi=>orbit%wfn
    r=>grid%r

    ALLOCATE(wfp(n),dum(n),ddum(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'allocation error in Get_Energy_EXX', i,n
       STOP
    ENDIF

    eex=0
    DO io=1,norbit
       occ=orbit%occ(io)
       IF (occ>threshold) THEN
          li=orbit%l(io);ni=orbit%np(io)
          DO jo=1,norbit
             occj=orbit%occ(jo)
             IF (occj>threshold) THEN
                lj=orbit%l(jo);nj=orbit%np(jo)
                lmax=li+lj
                lmin=ABS(li-lj)
                wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
                DO ll=lmin,lmax,2
                   CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                   IF (wgt>threshold) THEN
                      wgt=wgt*(2*ll+1)   ! because of apoisson convention
                      CALL apoisson(Grid,ll,n,wfp,dum)
                      ddum(2:n)=dum(2:n)*phi(2:n,jo)*phi(2:n,io)/r(2:n)
                      ddum(1)=0
                      eex=eex-wgt*integrator(Grid,ddum)/2
                   ENDIF
                ENDDO
             ENDIF
          ENDDO   !jo

       ENDIF
    ENDDO

    DEALLOCATE(wfp,dum,ddum)

  END SUBROUTINE Get_Energy_EXX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Get_Energy_EXX_VV       !!!!
  !!   Valence - valence part of the Fock exchange energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_Energy_EXX_VV(Grid,Orbit,eex)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    REAL(8), INTENT(OUT) :: eex

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),dum(:),ddum(:)
    REAL(8), ALLOCATABLE :: Fnu(:,:),Lnu(:,:),Snu(:),Mnunup(:,:),Cnu(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,io,jo,ok,norbit,last,nodes
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,kappa,fac,FCeex
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=orbit%norbit
    phi=>orbit%wfn
    r=>grid%r

    ALLOCATE(wfp(n),dum(n),ddum(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'allocation error in Get_Energy_EXX_FC', i,n
       STOP
    ENDIF


    eex=0
    DO io=1,norbit
       occ=Orbit%occ(io)
       IF (occ>threshold.AND.(.NOT.Orbit%iscore(io))) THEN
          li=Orbit%l(io);ni=Orbit%np(io)
          DO jo=1,norbit
             occj=orbit%occ(jo);fac=1
             IF (occj>threshold.AND.(.NOT.Orbit%iscore(jo))) THEN
                lj=orbit%l(jo);nj=orbit%np(jo)
                IF (Orbit%iscore(jo)) fac=2
                lmax=li+lj
                lmin=ABS(li-lj)
                wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
                DO ll=lmin,lmax,2
                   CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                   IF (wgt>threshold) THEN
                      wgt=wgt*(2*ll+1)*fac   ! because of apoisson convention
                      CALL apoisson(Grid,ll,n,wfp,dum)
                      ddum(2:n)=dum(2:n)*phi(2:n,jo)*phi(2:n,io)/r(2:n)
                      ddum(1)=0
                      eex=eex-wgt*integrator(Grid,ddum)/2
                   ENDIF
                ENDDO
             ENDIF
          ENDDO   !jo

       ENDIF
    ENDDO

    DEALLOCATE(wfp,dum,ddum)

  END SUBROUTINE Get_Energy_EXX_VV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Get_Energy_EXX_VC       !!!!
  !!   Valence - core part of the Fock exchange energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_Energy_EXX_VC(Grid,Orbit,eex)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(orbitinfo), INTENT(IN) :: Orbit
    REAL(8), INTENT(OUT) :: eex

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),dum(:),ddum(:)
    REAL(8), ALLOCATABLE :: Fnu(:,:),Lnu(:,:),Snu(:),Mnunup(:,:),Cnu(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,io,jo,ok,norbit,last,nodes
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,kappa,fac,FCeex
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=orbit%norbit
    phi=>orbit%wfn
    r=>grid%r

    ALLOCATE(wfp(n),dum(n),ddum(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'allocation error in Get_Energy_EXX_FC', i,n
       STOP
    ENDIF


    eex=0
    DO io=1,norbit
       occ=Orbit%occ(io)
       IF (occ>threshold.AND.(.NOT.Orbit%iscore(io))) THEN
          li=Orbit%l(io);ni=Orbit%np(io)
          DO jo=1,norbit
             occj=orbit%occ(jo);fac=1
             IF (occj>threshold.AND.(Orbit%iscore(jo))) THEN
                lj=orbit%l(jo);nj=orbit%np(jo)
                IF (Orbit%iscore(jo)) fac=2
                lmax=li+lj
                lmin=ABS(li-lj)
                wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
                DO ll=lmin,lmax,2
                   CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                   IF (wgt>threshold) THEN
                      wgt=wgt*(2*ll+1)*fac   ! because of apoisson convention
                      CALL apoisson(Grid,ll,n,wfp,dum)
                      ddum(2:n)=dum(2:n)*phi(2:n,jo)*phi(2:n,io)/r(2:n)
                      ddum(1)=0
                      eex=eex-wgt*integrator(Grid,ddum)/2
                   ENDIF
                ENDDO
             ENDIF
          ENDDO   !jo

       ENDIF
    ENDDO

    DEALLOCATE(wfp,dum,ddum)

  END SUBROUTINE Get_Energy_EXX_VC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  EXXwgt -- written by Xiao Xu   05-02-07  following Talman & Shadwick
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE EXXwgt(occi,occj,ni,li,nj,lj,ll,wgt)
    REAL(8), INTENT(IN) :: occi,occj
    INTEGER, INTENT(IN) :: li,lj,ll,ni,nj
    REAL(8), INTENT(OUT) :: wgt

    INTEGER ::BigJ,J1,J2,J3,VV
    REAL(8) ::ThreeJ

    wgt=0

    J1=li;
    J2=lj;
    J3=ll;
    BigJ=J1+J2+J3;
    VV=2*(2*li+1);

    ! Following the 3J symbol notes
    IF(MOD(BigJ,2) == 1)  THEN
       ThreeJ=0;
    ELSE
       ThreeJ=(-1)**(BigJ/2);
       ThreeJ=  ThreeJ* ( factorial(BigJ - &
&           2*J1)*factorial(BigJ - 2*J2)*factorial(BigJ - &
&           2*J3)/factorial(BigJ +1) )**(1.0d0/2.0d0);
       ThreeJ = ThreeJ*factorial(BigJ/2)/( &
&           factorial(BigJ/2-J1) *factorial(BigJ/2-J2)*factorial(BigJ/2 - J3));
    ENDIF

    ThreeJ=ThreeJ*ThreeJ

    ! Occ factors
    IF (nj /= ni)   wgt = 0.5d0*occi*occj*ThreeJ;
    IF( (nj == ni).AND.(ll /= 0))   &
&        wgt = 0.5d0*ThreeJ*occi*(occj-1)*VV/(VV-1);
    IF((nj == ni).AND.(ll == 0))    wgt = occi;

  END SUBROUTINE EXXwgt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine to calculate Condon-Shortley integral
!
!   RL = \int_0^\infty dr dr' r<^L     w1(r)*w2(r)*w3(r')*w4(r')
!                             --
!                             r>^(L+1)
!   Integration is performed over entire radial range
!     No error checking is done
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CondonShortley(Grid,L,w1,w2,w3,w4,RL)
    Type(GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: L
    REAL(8), INTENT(IN) :: w1(:),w2(:),w3(:),w4(:)
    REAL(8), INTENT(OUT) :: RL

    INTEGER :: n,i
    REAL(8), allocatable :: dum(:),arg(:)

    n=Grid%n
    allocate(arg(n),dum(n))
         arg(:)=w1(:)*w2(:)
         CALL apoisson(Grid,L,n,arg,dum)
         dum(2:n)=dum(2:n)*w3(2:n)*w4(2:n)/Grid%r(2:n)
         dum(1)=0
         RL=integrator(Grid,dum)*(2*L+1)   ! because of apoisson convention

    deallocate(arg,dum)
  END  SUBROUTINE CondonShortley

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc_dexdphi_io(Grid,Orbit,io,res)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_dexdphi_io(Grid,Orbit,io,res)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: res(:)

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,jo,ok,norbit,ko
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=Orbit%norbit
    ALLOCATE(wfp(n),vl(n),dum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_dexdphi_io allocation error:', n,norbit,ok
       STOP
    ENDIF

    phi=>Orbit%wfn
    r=>Grid%r

    res=0.d0
    occ=orbit%occ(io)
    li=Orbit%l(io);ni=Orbit%np(io)
    DO jo=1,norbit
       occj=Orbit%occ(jo); vl=0.d0
       IF (occj>threshold) THEN
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL apoisson(Grid,ll,n,wfp,dum)
                vl(1:n)=vl(1:n)-wgt*dum(1:n)*phi(1:n,jo)/occ
                !write(6,*) 'In vx', io,jo,wgt/((2*ll+1)*occ)
             ENDIF
          ENDDO
       ENDIF
       res=res+vl
       ! testing
       !vl(1)=0;vl(2:n)=vl(2:n)/Grid%r(2:n)
       !do ko=1,norbit
       !   if(Orbit%l(ko)==li) then
       !     write(6,'("vx sum", 4i4,2x,1p,1e15.7)') &
       !&           io,jo,ko,jo,overlap(Grid,Orbit%wfn(:,ko),vl)
       !   endif
       !enddo
    ENDDO   !jo


    DEALLOCATE(dum,wfp,vl)
  END SUBROUTINE Calc_dexdphi_io

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc_Xp(Grid,Orbit,ni,li,wfn,res)
  !      calculate X_p for possibly
  !        unoccupied state with angular moment li, stored in wfn(:)
  !         Note: ni should not be equal to real shell index
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_Xp(Grid,Orbit,ni,li,wfn,res,lng)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    INTEGER, INTENT(IN) :: ni,li
    REAL(8), INTENT(IN) :: wfn(:)
    REAL(8), INTENT(OUT) :: res(:)
    INTEGER, OPTIONAL, INTENT(IN) :: lng

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:)
    INTEGER :: i,j,k,n,m,ll,lj,lmin,lmax,jo,ok,norbit,ko
    INTEGER :: nu,nup,nj,last
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=Orbit%norbit
    ALLOCATE(wfp(n),vl(n),dum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_Xp allocation error:', n,norbit,ok
       STOP
    ENDIF

    last=n
    if (present(lng)) last=min(last,lng)
    phi=>Orbit%wfn
    r=>Grid%r

    res=0.d0
    occ=1.d0

    DO jo=1,norbit
       occj=Orbit%occ(jo); vl=0.d0
       IF (occj>threshold) THEN
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:last)=wfn(1:last)*phi(1:last,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,ni,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL apoisson(Grid,ll,last,wfp,dum)
                vl(1:last)=vl(1:last)-wgt*dum(1:last)*phi(1:last,jo)/occ
             ENDIF
          ENDDO
       res(1:last)=res(1:last)+vl(1:last)
       ENDIF
    ENDDO   !jo


    DEALLOCATE(dum,wfp,vl)
  END SUBROUTINE Calc_Xp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc_dexdphi_io_v(Grid,Orbit,io,res)
!!  valence-valence contribution for use in frozencore calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_dexdphi_io_v(Grid,Orbit,io,res)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: res(:)

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,jo,ok,norbit
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    if(.not.frozencorecalculation) then
       write(6,*) 'In Calc_dexdphi_io_v  but not frozencore calculation'
       stop
    endif

    if(Orbit%iscore(io)) then
       write(6,*) 'In Calc_dexdphi_io_v but not valence orbit ', io
       stop
    endif

    n=Grid%n
    norbit=Orbit%norbit
    ALLOCATE(wfp(n),vl(n),dum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_dexdphi_io_v allocation error:', n,norbit,ok
       STOP
    ENDIF

    phi=>Orbit%wfn
    r=>Grid%r

    occ=orbit%occ(io)
    vl=0; li=Orbit%l(io);ni=Orbit%np(io)
    DO jo=1,norbit
      If(.not.Orbit%iscore(jo)) then
       occj=Orbit%occ(jo)
       IF (occj>threshold) THEN
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL apoisson(Grid,ll,n,wfp,dum)
                vl(1:n)=vl(1:n)-wgt*dum(1:n)*phi(1:n,jo)/occ
             ENDIF
          ENDDO
       ENDIF
      Endif
    ENDDO   !jo

    res=vl

    DEALLOCATE(dum,wfp,vl)
  END SUBROUTINE Calc_dexdphi_io_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calc_dexdphi_io_c(Grid,Orbit,io,res)
!!  valence-core contribution for use in frozencore calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_dexdphi_io_c(Grid,Orbit,io,res)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: res(:)

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,jo,ok,norbit
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    if(.not.frozencorecalculation) then
       write(6,*) 'In Calc_dexdphi_io_c  but not frozencore calculation'
       stop
    endif

    if(Orbit%iscore(io)) then
       write(6,*) 'In Calc_dexdphi_io_c but not valence orbit ', io
       stop
    endif

    n=Grid%n
    norbit=Orbit%norbit
    ALLOCATE(wfp(n),vl(n),dum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_dexdphi_io_c allocation error:', n,norbit,ok
       STOP
    ENDIF

    phi=>Orbit%wfn
    r=>Grid%r

    occ=orbit%occ(io)
    vl=0; li=Orbit%l(io);ni=Orbit%np(io)
    DO jo=1,norbit
      If(Orbit%iscore(jo)) then
       occj=Orbit%occ(jo)
       IF (occj>threshold) THEN
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL apoisson(Grid,ll,n,wfp,dum)
                vl(1:n)=vl(1:n)-wgt*dum(1:n)*phi(1:n,jo)/occ
             ENDIF
          ENDDO
       ENDIF
      Endif
    ENDDO   !jo

    res=vl

    DEALLOCATE(dum,wfp,vl)
  END SUBROUTINE Calc_dexdphi_io_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc_expot_io(Grid,Orbit,io,res)
!        Determine exchange potential contribution from single orbital
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_expot_io(Grid,Orbit,io,res)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: res(:)

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,jo,ok,norbit
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=Orbit%norbit
    ALLOCATE(wfp(n),vl(n),dum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_expot_io allocation error:', n,norbit,ok
       STOP
    ENDIF

    phi=>Orbit%wfn
    r=>Grid%r

    res=0
    occ=orbit%occ(io)
    if(occ<threshold) return
    vl=0; li=Orbit%l(io);ni=Orbit%np(io)
    jo=io
       occj=Orbit%occ(jo)
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL apoisson(Grid,ll,n,wfp,dum)
                vl(1:n)=vl(1:n)-wgt*dum(1:n)/occ
             ENDIF
          ENDDO

    res=vl

    DEALLOCATE(dum,wfp,vl)
  END SUBROUTINE Calc_expot_io

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  calc_psi
  !     Solve inhomogeneous differential equation to determine psi_io
  !         of the form (H-e_io) psi_io = rhs
  !         H=T+rv/r
  !            psi_io is orthogonalized to phi_io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_psi(Grid,l,eig,rv,phi,rhs,psi)
    TYPE(GridInfo), INTENT(IN):: Grid
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: eig
    REAL(8), INTENT(IN) :: rv(:),phi(:)
    REAL(8), INTENT(INOUT) :: rhs(:)
    REAL(8), INTENT(OUT) :: psi(:)

    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:),ddum(:)
    INTEGER :: i,j,k,n,ok
    REAL(8) :: occ,term,term1,zeroval,energy,tol
    REAL(8), PARAMETER :: tol0=1.d-12,threshold=1.d-8
    INTEGER :: counter=1

    tol=tol0
    n=Grid%n
    ALLOCATE(wfp(n),vl(n),dum(n),ddum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_psi allocation error:', n,ok
       STOP
    ENDIF

    r=>Grid%r

    psi=0

    !!  check orthogonality
    !dum=rhs*phi
    !term=integrator(Grid,dum)
    !rhs=rhs-term*phi

    !CALL inhomo_bound_numerov(Grid,l,n,eig,rv,rhs,psi)

    !call filter(n,psi,machine_zero)
    !! orthogonalize to  phi(:)
    !dum=phi(:)*psi(:); ddum=phi(:)**2
    !term=integrator(Grid,dum)
    !term1=integrator(Grid,ddum)
    !psi(:)=psi(:)-term*phi(:)/term1
    !!dum=phi(:,io)*psi(:)
    !!write(6,*) 'ortho ', term,term1,integrator(Grid,dum)

    !!do jo=1,norbit
    !!   dum=EigOrbit%wfn(:,jo)*psi(:)
    !!   write(6,*) 'ortho' , io,jo,integrator(Grid,dum)
    !!enddo

    !CALL inhomo_numerov_SVD(Grid,l,n,eig,tol,rv,rhs,psi,phi)

    psi=0.d0
    CALL inhomo_bound_numerov(Grid,l,n,eig,rv,rhs,psi,phi)
    call filter(n,psi,machine_zero)

    !do i=1,n
    !   write(500,'(1P20e15.7)') Grid%r(i),phi(i),rhs(i),psi(i),wfp(i)
    !enddo
    !stop
    DEALLOCATE(wfp,vl,dum,ddum)

  END SUBROUTINE Calc_psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc_dexdphi_io_diag_byphi(Grid,Orbit,io,res)
!     Version for HF equations -- diagonal contribution
!         without phi factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_dexdphi_io_diag_byphi(Grid,Orbit,io,res)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: res(:)

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,jo,ok,norbit
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=Orbit%norbit
    ALLOCATE(wfp(n),vl(n),dum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_dexdphi_io allocation error:', n,norbit,ok
       STOP
    ENDIF

    phi=>Orbit%wfn
    r=>Grid%r

    occ=orbit%occ(io)
    vl=0; li=Orbit%l(io);ni=Orbit%np(io)
    DO jo=io,io
       occj=Orbit%occ(jo)
       IF (occj>threshold) THEN
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL apoisson(Grid,ll,n,wfp,dum)
                vl(1:n)=vl(1:n)-wgt*dum(1:n)/occ
             ENDIF
          ENDDO
       ENDIF
    ENDDO   !jo

    res=vl

    DEALLOCATE(dum,wfp,vl)
  END SUBROUTINE Calc_dexdphi_io_diag_byphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc_dexdphi_io_offdiag(Grid,Orbit,io,res)
!     Version for HF equations -- off-diagonal contribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_dexdphi_io_offdiag(Grid,Orbit,io,res)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(OrbitInfo), INTENT(IN):: Orbit
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: res(:)

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,jo,ok,norbit
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=Orbit%norbit
    ALLOCATE(wfp(n),vl(n),dum(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_dexdphi_io allocation error:', n,norbit,ok
       STOP
    ENDIF

    phi=>Orbit%wfn
    r=>Grid%r

    occ=orbit%occ(io)
    vl=0; li=Orbit%l(io);ni=Orbit%np(io)
    DO jo=1,norbit
       occj=Orbit%occ(jo)
       IF (occj>threshold.and.jo/=io) THEN
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL apoisson(Grid,ll,n,wfp,dum)
                vl(1:n)=vl(1:n)-wgt*dum(1:n)/occ
             ENDIF
          ENDDO
       ENDIF
    ENDDO   !jo

    res=vl

    DEALLOCATE(dum,wfp,vl)
  END SUBROUTINE Calc_dexdphi_io_offdiag

End module fock
