Module exx_pseudo
   USE anderson_driver
   USE atomdata
   USE exx_mod
   USE fock
   USE globalmath
   USE gridmod
   USE pseudodata
   USE pseudo_sub

   IMPLICIT NONE

   INTEGER, PRIVATE :: S_n,S_nocc,S_nbase,S_irc,S_constr
   INTEGER, ALLOCATABLE, TARGET, PRIVATE  :: S_l(:),S_bmap(:),S_omap(:)
   REAL(8), ALLOCATABLE, TARGET, PRIVATE  :: S_gvv(:,:),S_tgvv(:,:)
   REAL(8), ALLOCATABLE, TARGET, PRIVATE  :: S_W(:,:),S_X(:,:),S_U(:),S_E(:)
   REAL(8), ALLOCATABLE, TARGET, PRIVATE  :: S_Vx(:,:),S_shift(:)
   REAL(8), ALLOCATABLE, TARGET, PRIVATE  :: S_Vxvale(:),S_tVxvale(:)
   REAL(8), ALLOCATABLE, TARGET, PRIVATE  :: S_dExdtphi(:,:),S_UOphi(:,:)

   TYPE(Anderson_context), PRIVATE :: AC
   TYPE(GridInfo), POINTER, PRIVATE :: Gridwk
   TYPE(PseudoInfo), POINTER, PRIVATE :: PAW
 !!!! Note: S_omap(:) maps occupied valence states to atomdata indices
 !!!!       S_bmap(:) maps occupied valence states to PAW basis indices


   CONTAINS

    SUBROUTINE EXXKLI_pseudoVx(Grid,PAW,rtvx)
        TYPE(GridInfo), INTENT(IN) :: Grid
        TYPE(PseudoInfo), INTENT(INOUT) :: PAW
        REAL(8), INTENT(INOUT) :: rtvx(:)

        INTEGER :: i,j,k,l,n,io,jo
        REAL(8), ALLOCATABLE :: dum1(:),dum2(:),dum3(:),den(:),tden(:)
        REAL(8) :: occ

        n=Grid%n
        Allocate(dum1(n),dum2(n),dum3(n),den(n),tden(n))

        rtvx=0;dum1=0;dum2=0;den=0;tden=0
        do io=1,PAW%TOCCWFN%norbit
           occ=PAW%TOCCWFN%occ(io)
           den=den+occ*PAW%OCCwfn%wfn(:,io)**2
           tden=tden+occ*PAW%tOCCwfn%wfn(:,io)**2
           call Calc_dexdphi_io(Grid,PAW%OCCwfn,io,dum3);PAW%OCCwfn%X(:,io)=dum3
               ! not really needed -- just for checking
           call Calc_tdexdphi_io(Grid,PAW,io,dum3)
               PAW%tOCCwfn%X(:,io)=dum3
           dum1=dum1+occ*(PAW%OCCwfn%X(:,io)*PAW%OCCwfn%wfn(:,io) &
&                         + (PAW%OCCwfn%wfn(:,io)**2)*HSZ%U(io)*Grid%r)
           dum2=dum2+occ*(PAW%tOCCwfn%X(:,io)*PAW%tOCCwfn%wfn(:,io) &
&                         + (PAW%tOCCwfn%wfn(:,io)**2)*HSZ%U(io)*Grid%r)
        enddo

        rtvx(1)=0; dum3=0
        do i=2,n
           if (den(i)>machine_zero) then
              dum3(i)=dum1(i)/den(i)
           else
              dum3(i)=0
           endif
           if (tden(i)>machine_zero) then
              rtvx(i)=dum2(i)/tden(i)
           else
              rtvx(i)=0
           endif
        enddo

       OPEN(1001,file='checkpseudovx',form='formatted')
          do i=1,n
             write(1001,'(1p,16e15.7)') Grid%r(i),dum1(i),dum2(i),dum3(i),&
&                rtvx(i),den(i),tden(i)
          enddo
       CLOSE(1001)

       OPEN(1001,file='checkagain',form='formatted')
          do i=1,n
             write(1001,'(1p,50e15.7)') Grid%r(i),(PAW%OCCWFN%wfn(i,io),&
&                PAW%tOCCWFN%wfn(i,io),io=1,PAW%TOCCWFN%norbit)
          enddo
       CLOSE(1001)


        Deallocate(dum1,dum2,dum3,den,tden)
    END SUBROUTINE EXXKLI_pseudoVx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calc_texpot_io(Grid,Orbit,io,res)
!        Determine exchange potential contribution from single pseudo orbital
!        Appropriate for reference configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_tdexdphi_io(Grid,PAW,io,res)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: res(:)

    REAL(8), POINTER :: r(:)
    TYPE(OrbitInfo), POINTER :: Orbit, tOrbit
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),arg(:),dum(:),dum1(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,jo,ok
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    write(6,*) 'Calc_tecpot_io ',io; call flush(6)
    n=Grid%n
    Orbit=>PAW%OCCwfn
    tOrbit=>PAW%tOCCwfn
    ALLOCATE(wfp(n),vl(n),arg(n),dum(n),dum1(n),stat=ok)
    IF (ok/=0) THEN
       WRITE(6,*) 'calc_expot_io allocation error:', n,Orbit%norbit,ok
       STOP
    ENDIF

    r=>Grid%r

    res=0
    occ=orbit%occ(io)
    if(occ<threshold) return
    li=Orbit%l(io);ni=Orbit%np(io)
    DO jo=1,Orbit%norbit
       occj=Orbit%occ(jo);vl=0
          lj=Orbit%l(jo);nj=Orbit%np(jo)
          IF(occj>threshold) THEN
             lmax=li+lj
             lmin=ABS(li-lj)
             wfp(1:n)=tOrbit%wfn(1:n,io)*tOrbit%wfn(1:n,jo)
             arg(1:n)=Orbit%wfn(1:n,io)*Orbit%wfn(1:n,jo)-wfp(1:n)
             DO ll=lmin,lmax,2
                CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                IF (wgt>threshold) THEN
                   wgt=wgt*(2*ll+1)   ! because of apoisson convention
                   dum=(r**ll)*arg;term1=integrator(Grid,dum)
                   dum1=wfp+term1*PAW%g(:,ll+1)
                   CALL apoisson(Grid,ll,n,dum1,dum)
                   vl(1:n)=vl(1:n)-wgt*dum(1:n)*tOrbit%wfn(1:n,jo)/occ
                ENDIF
             ENDDO
           ENDIF
        res=res+vl
     ENDDO


    DEALLOCATE(dum,wfp,arg,vl)
  END SUBROUTINE Calc_tdexdphi_io

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Assumes tildeg(r) = \sum_n=1^N C_n r^{lead+n+l}

  SUBROUTINE EXX_pseudoVx(Grid,Orbit,PAWin,trvx)
     Type(Gridinfo), TARGET,  INTENT(IN) :: Grid
     Type(Orbitinfo), INTENT(IN) :: Orbit
     Type(PseudoInfo), TARGET, INTENT(IN) :: PAWin
     REAL(8), INTENT(OUT) :: trvx(:)

     LOGICAL :: success
     INTEGER :: i,j,k,l,m,n,nbase,io,ip,iq,ib,ic,nocc,kmax,kmin,lp,lq,irc,ind
     INTEGER :: ii,jj,kk,loop,outer,p
     INTEGER :: nequations,nterms
     INTEGER, POINTER :: map(:)
     INTEGER, PARAMETER :: nshift=4,nmatch=3,lead=4
     REAL(8) :: occp,occq,wgt,xxx,yyy,r1
     REAL(8), ALLOCATABLE :: trvx0(:),d(:),d1(:),d2(:),d3(:),v(:),arg(:),F2(:,:)
     REAL(8), ALLOCATABLE :: F(:),rho(:)
     INTEGER, PARAMETER :: mxloop=100
     REAL(8), PARAMETER :: threshold=1.d-8,CondNo=1.d12
     REAL(8) :: svdcut=1.d-8
     LOGICAL :: parameterfileexists
     CHARACTER(120) :: inputline

     Gridwk=>Grid
     PAW=>PAWin
     n=Grid%n; S_n=n;  irc=PAW%irc; S_irc=irc

     Inquire(file='in_parameters',exist=parameterfileexists)

     If(parameterfileexists) then
         Open(1002,file='in_parameters')
            do
              READ(1002,'(a)') inputline
              Call Uppercase (inputline)
              If (TRIM(inputline)=='<SVDCUTOFF>') exit
            enddo
            READ(1002,*)  svdcut
            WRITE(6,*) 'Modified SVD Cutoff ', svdcut
         close(1002)
     endif

     trvx=0

!     nocc=0; nbase=PAW%nbase
!       do io=1,nbase
!          if (PAW%occ(io)>threshold) nocc=nocc+1
!       enddo
!
!     ALLOCATE(S_bmap(nocc),S_omap(nocc),S_E(nocc),S_l(nocc),S_U(nocc),&
!&          S_W(n,nocc),S_X(nbase,nocc),S_Vx(nbase,nocc), S_shift(n), &
!&          S_dExdtphi(n,nocc),S_UOphi(n,nocc),S_Vxvale(n),S_tVxvale(n))
!     ALLOCATE(trvx0(n),d(n),d1(n),d2(n),d3(n),v(n),arg(irc),F2(n,nocc),&
!&           F(n),rho(n))
!
!     map=>S_bmap
!     S_nocc=nocc; S_nbase=nbase;
!     S_bmap=0; S_omap=0
!
!     io=0
!     do ib=1,nbase
!        if (PAW%occ(ib)>threshold) then
!           io=io+1
!           map(io)=ib
!        endif
!     enddo
!
!     Do ip=1,nocc
!        S_E(ip)=PAW%eig(map(ip))
!        S_l(ip)=PAW%l(map(ip))
!        i=0
!        do io=1,Orbit%norbit
!           if (ABS(Orbit%eig(io)-S_E(ip))<1.d-8) then
!              i=io; S_omap(ip)=io; write(6,*) 'ip io map', ip,io ; exit
!           endif
!        enddo
!        if (i==0) then
!           write(6,*) 'Error:  No match found for ',ip,S_E(ip)
!           stop
!        endif
!    Enddo
!
!    allocate(S_gvv(n,nocc),S_tgvv(n,nocc))
!    do ip=1,nocc
!       S_gvv(:,ip)=HSZ%psivale(:,S_omap(ip))
!       S_tgvv(:,ip)=HSZ%psivale(:,S_omap(ip))
!    enddo
!
!
!     S_W=0;S_X=0
!     F2=0;F=0;rho=0;
!     Do ip=1,nocc
!        occp=PAW%occ(map(ip))
!        lp=PAW%l(map(ip))
!        Do iq=1,nocc
!           occq=PAW%occ(map(iq))
!           lq=PAW%l(map(iq))
!           kmax=lp+lq
!           kmin=ABS(lp-lq)
!           d1=PAW%otphi(:,map(ip))*PAW%otphi(:,map(iq))
!           d3=PAW%ophi(:,map(ip))*PAW%ophi(:,map(iq))
!           Do k=kmin,kmax,2
!              CALL EXXwgt(occp,occq,map(ip),lp,map(iq),lq,k,wgt)
!              IF (wgt>threshold) THEN
!                 wgt=(wgt*(2*k+1))/occp   ! because of apoisson convention
!                 CALL CompensationRHO(Grid,PAW,k,map(ip),map(iq),d2)
!                 d=d1+d2
!                 CALL apoisson(Grid,k,n,d,trvx0)
!                 S_W(:,ip)=S_W(:,ip)-wgt*trvx0(:)*PAW%otphi(:,map(iq))
!                 CALL apoisson(Grid,k,n,d3,v)
!                 F2(:,ip)=F2(:,ip)-wgt*v(:)*PAW%ophi(:,map(iq))
!                 do ib=1,nbase
!                    if (PAW%l(ib)==lp) then
!                       d=PAW%ophi(:,ib)*PAW%ophi(:,map(iq))*v - &
!&                            PAW%otphi(:,ib)*PAW%otphi(:,map(iq))*trvx0
!                       d(1)=0; d(2:n)=d(2:n)/Grid%r(2:n)
!                       S_X(ib,ip)=S_X(ib,ip)-wgt*integrator(Grid,d,1,irc)
!                       write(6,*) 'Check ',ib,ip,integrator(Grid,d,1,irc),&
!&                        integrator(Grid,d)
!                    endif
!                 enddo
!              ENDIF
!           ENDDO   !k
!        ENDDO      !iq
!
!
!
!      ! (re)calculate U factor
!        d=(F2(:,ip)-HSZ%rVxvale(:)*PAW%ophi(:,map(ip)))*PAW%ophi(:,map(ip))
!        d(1)=0;d(2:n)=d(2:n)/Grid%r(2:n)
!        S_U(ip)=integrator(Grid,d)
!        write(6,*) 'Checking U ', ip,map(ip),integrator(Grid,d),&
!&          HSZ%Uvale(S_omap(ip))
!        d=d-HSZ%Uvale(S_omap(ip))*PAW%ophi(:,map(ip))**2
!        write(6,*) 'FC rhs',ip,integrator(Grid,d),integrator(Grid,d,1,irc)
!        d=S_W(:,ip)*PAW%otphi(:,map(ip))
!        d(1)=0;d(2:n)=d(2:n)/Grid%r(2:n)
!        write(6,*) 'Checking dtEx ', ip,map(ip),&
!&            integrator(Grid,d)+S_X(map(ip),ip),S_X(map(ip),ip)
!        d=F2(:,ip)*PAW%ophi(:,map(ip))
!        d(1)=0;d(2:n)=d(2:n)/Grid%r(2:n)
!        write(6,*) 'Checking dEx ', ip,map(ip),integrator(Grid,d)
!        d=S_W(:,ip)*PAW%otphi(:,map(ip))
!        d(1)=0;d(2:n)=d(2:n)/Grid%r(2:n)
!        F(:)=F(:)+occp*d(:)
!        rho(:)=rho(:)+occp*PAW%otphi(:,map(ip))**2
!        do iq=1,nocc
!           if (iq==ip) then
!              F(:)=F(:)-occp*HSZ%Uvale(S_omap(ip))*PAW%otphi(:,map(ip))**2
!           else
!              !if (S_l(iq)==S_l(ip)) then
!              !F(:)=F(:)-occp*HSZ%LMBDvale(S_omap(ip),S_omap(iq))&
!&             !   *PAW%otphi(:,map(ip))*PAW%otphi(:,map(iq))
!              !endif
!           endif
!        enddo
!     ENDDO
!  ! accumulate dExdtphi
!  S_dExdtphi=0 ; S_UOphi=0
!    do ip=1,nocc
!       S_dExdtphi(2:n,ip)=S_W(2:n,ip)/Grid%r(2:n)
!       S_UOphi(:,ip)=PAW%otphi(:,map(ip))
!       do ib=1,nbase
!          if (PAW%l(ib)==S_l(ip)) then
!              S_dExdtphi(:,ip)=S_dExdtphi(:,ip)+S_X(ib,ip)*PAW%otp(:,ib)
!              S_UOphi(:,ip)=S_UOphi(:,ip)+PAW%Oij(ib,map(ip))*PAW%otp(:,ib)
!          endif
!       enddo
!       S_UOphi(:,ip)=HSZ%Uvale(S_omap(ip))*S_UOphi(:,ip)
!    enddo
!
!  !  Calculated matrix element <\phi(ib)|VxVale|\phi(ip)>
!    S_Vx=0
!    do ip=1,nocc
!       lp=PAW%l(S_bmap(ip))  ;
!       do ib=1,nbase
!          if (lp==PAW%l(ib)) THEN
!             d=HSZ%rVxVale(:)*PAW%ophi(:,ib)*PAW%ophi(:,S_bmap(ip))
!             d(1)=0; d(2:n)=d(2:n)/Gridwk%r(2:n)
!             S_Vx(ib,ip)=integrator(Gridwk,d,1,irc)
!          endif
!          write(6,'(" vx ",2i5,1p,1e15.7)') ib,ip,S_Vx(ib,ip)
!       enddo
!    enddo
!
!    ! Note that S_W is r*(desired function)
!
!     Open(1001,file='testF1',form='formatted')
!     do i=1,n
!        write(1001,'(1p,50e15.7)') Grid%r(i),HSZ%rVxvale(i),&
!&            (S_W(i,ip),F2(i,ip),ip=1,nocc)
!     enddo
!     close(1001)
!
!     Do ip=1,nocc
!        do ib=1,nbase
!           write(6,*) 'X  ', ip,ib,S_X(ib,ip)
!        enddo
!     enddo
!
!   call Init_trvx(Grid,irc,HSZ%rVxvale,trvx)
!   S_Vxvale=HSZ%rVxvale; S_tVxvale=trvx
!
!     Open(1001,file='trvx0',form='formatted')
!     do i=1,n
!        write(1001,'(1p,50e15.7)') Grid%r(i),trvx(i),HSZ%rVxvale(i),&
!&                        HSZ%rVxcore(i)
!     enddo
!     close(1001)
!
!   arg(1:irc)=trvx(1:irc)
!   CALL InitAnderson_dr(AC,6,5,irc,0.1d0,1.d3,100,&
!&                 1.d-8,1.d-16,.TRUE.)
!   CALL DoAndersonMix(AC,arg,xxx,tVXsub1,success)
!   trvx(1:irc)=arg(1:irc)
!      WRITE(6,*) 'Finished trvx iter # ',AC%CurIter,AC%res
!   CALL FreeAnderson(AC)
!
!     Open(1001,file='trvx',form='formatted')
!     do i=1,n
!        write(1001,'(1p,50e15.7)') Grid%r(i),trvx(i),HSZ%rVxvale(i),&
!&                        HSZ%rVxcore(i)
!     enddo
!     close(1001)
!
!     Deallocate(trvx0,d,d1,d2,d3,v,F2,arg,rho,F)
!
  END SUBROUTINE EXX_pseudoVx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ib is the PAW basis index of the outer
!     most orbital for initializing tVx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE InitialVx_pseudo(Grid,ib,trvx0)
    TYPE(GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: ib
    REAL(8), INTENT(OUT) :: trvx0(:)

    INTEGER :: lp,kmax,kmin,k,n
    REAL(8) :: occp,wgt
    REAL(8), ALLOCATABLE :: d(:),d1(:),d2(:)
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    allocate(d(n),d1(n),d2(n))

        trvx0=0
        occp=PAW%occ(ib)
        lp=PAW%l(ib)
           kmax=2*lp
           kmin=0
           d1=(PAW%otphi(:,ib)**2)
           Do k=kmin,kmax,2
              CALL EXXwgt(occp,occp,ib,lp,ib,lp,k,wgt)
              IF (wgt>threshold) THEN
                 wgt=(wgt*(2*k+1))/occp   ! because of apoisson convention
                 CALL CompensationRHO(Grid,PAW,k,ib,ib,d)
                 d=d1+d
                 CALL apoisson(Grid,k,n,d,d2)
                 trvx0=trvx0-wgt*d2
              ENDIF
           ENDDO   !k

     deallocate(d,d1,d2)
  END SUBROUTINE InitialVx_pseudo

  SUBROUTINE Init_trvx(Grid,p,trvxin,trvxout)
    Type(GridInfo),INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: p       ! smooth potential for r<r(irc+1)
    REAL(8), INTENT(IN) ::trvxin(:)
    REAL(8), INTENT(OUT) ::trvxout(:)

    INTEGER :: n,i
    REAL(8) :: r1u1,r2u2,r3u3,d0rc, d1rc,d2rc,rc,x
    REAL(8) :: u1,u2,u3
    REAL(8), ALLOCATABLE :: d1(:),d2(:)

    n=Grid%n
    ALLOCATE(d1(n),d2(n))

    CALL derivative(Grid,trvxin,d1)
    CALL derivative(Grid,d1,d2)

    d0rc=trvxin(p)
    d1rc=d1(p)
    d2rc=d2(p)
    rc=Grid%r(p)

    r2u2 = (-d0rc + rc*d1rc -(1.0d0/3.0d0)*rc*rc*d2rc)*3
    r3u3 = (rc*d1rc -d0rc - r2u2)*(1.0d0/2.0d0)
    r1u1 = d0rc - r2u2 - r3u3

    u1 = r1u1/rc;
    u2 = r2u2/(rc*rc)
    u3 = r3u3/(rc**3)

    trvxout=trvxin

    do i=1,p
       x = Grid%r(i)
       trvxout(i) = x*u1 + x*x*u2 + x*x*x*u3
    enddo

     open(unit=2001,file='rvxsmooth.txt')
     n=Grid%n
     do i=1,n
           write(2001,'(1p,6e15.7)') Grid%r(i), trvxin(i),trvxout(i)
     enddo
     close(2001)

  END SUBROUTINE Init_trvx

  SUBROUTINE CompensationRHO(Grid,PAW,k,alpha,beta,rho)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER,INTENT(IN) :: k,alpha,beta
    REAL(8),INTENT(INOUT) :: rho(:)

    REAL(8) , ALLOCATABLE :: Temp1(:),Temp2(:)
    INTEGER :: n,irc
    REAL(8) :: n1
    REAL(8) ,POINTER :: r(:)

    n=Grid%n
    ALLOCATE(Temp1(n),Temp2(n))
    r=>Grid%r
    irc = PAW%irc

    Temp1 = PAW%ophi(:,alpha)*PAW%ophi(:,beta) - &
&      PAW%otphi(:,alpha)*PAW%otphi(:,beta)
    Temp1 = Temp1 * (r(:)**k)
    n1 = integrator(Grid,Temp1,1,irc)
    rho(:)=(r(:)**(k+2))*PAW%hatshape(:)
    Temp1=(r(:)**(k))*rho(:)
    rho(:) =n1*rho(:)/integrator(Grid,Temp1)

  END SUBROUTINE CompensationRHO


!!!!!!!
!  tVXsub0  -- version without projector constraints
!!!!!!!
  SUBROUTINE tVXsub0(w,en,grad,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: en
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io,l,iter,last,irc,ip,nocc,nbase,many,ib,ic,jo,p
    REAL(8) :: x,y,energy,occ,boundaryv
    REAL(8), POINTER :: rv(:)
    REAL(8),ALLOCATABLE :: shift(:),dvx(:,:),tg(:,:),tu(:),rhs(:),tw(:)
    INTEGER:: fcount=0
    REAL(8), PARAMETER :: tol=1.d-10, wgt0=0.01d0
    CHARACTER(4) :: stuff

    success=.TRUE.
    irc=S_irc; nocc=S_nocc; n=S_n; nbase=PAW%nbase; p=irc+1

    rv=>PAW%rVeff

    allocate(shift(n),tg(n,nocc),tu(n),rhs(n),tw(n),&
&    dvx(nbase,nocc))
    shift=0;grad=0

   !  Matrix element of \tilde{V}_x = w (input)
    do ip=1,nocc
       dvx(:,ip)=0; l=PAW%l(S_bmap(ip))  ;
       do ib=1,nbase
          if (l==PAW%l(ib)) THEN
          tu(1:irc)=w(1:irc)*PAW%otphi(1:irc,ib)*PAW%otphi(1:irc,S_bmap(ip))
          tu(1)=0; tu(2:irc)=tu(2:irc)/Gridwk%r(2:irc)
          dvx(ib,ip)=S_Vx(ib,ip)-integrator(Gridwk,tu,1,irc)
          endif
          write(6,'(" dvx ",2i5,1p,1e15.7)') ib,ip,dvx(ib,ip)
       enddo
    enddo
    call flush(6)
    tw=S_tVxvale
    tw(1:irc)=w(1:irc)
    shift=0 ;
    do ip=1,nocc
       occ=PAW%occ(S_bmap(ip))
       energy=S_E(ip)
       l=S_l(ip)
       rhs=0;
       rhs(2:n)=S_dExdtphi(2:n,ip)  &
&            -tw(2:n)*PAW%otphi(2:n,S_bmap(ip))/Gridwk%r(2:n)
       rhs(2:n)=rhs(2:n)-S_UOphi(2:n,ip)

       do ib=1,nbase
          if (l==PAW%l(ib)) then
             rhs(:)=rhs(:)-dvx(ib,ip)*PAW%otp(:,ib)
          endif
       enddo

      write(6,*)' Check smooth rhs ip ', ip,  &
&                overlap(Gridwk,rhs,PAW%otphi(:,S_bmap(ip)),1,irc),&
&                overlap(Gridwk,rhs,PAW%otphi(:,S_bmap(ip)))


      call inhomo_numerov_SVD_bv(Gridwk,l,irc+1,energy,tol,rv,rhs,S_tgvv(:,ip))

      ! check orthogonalities

      do ib=1,nbase
         if (l==PAW%l(ib)) then
            write(6,*) ' ip ib <g|p> = ', ip,ib ,&
&             overlap(Gridwk,S_tgvv(:,ip),PAW%otp(:,ib))
            write(6,*) ' ip ib <g|tphi> = ', &
&             overlap(Gridwk,S_tgvv(:,ip),PAW%otphi(:,ib))
         endif
      enddo
   enddo   !ip

   shift=0;tu=0
   do ip=1,nocc
      occ=PAW%occ(S_bmap(ip))

      shift(:)=shift(:)+2*occ*S_tgvv(:,ip)*PAW%otphi(:,S_bmap(ip))
   enddo

    grad=0
    grad(2:irc)=-shift(2:irc)/Gridwk%r(2:irc)
    err=DOT_PRODUCT(grad(1:irc),grad(1:irc))
    en=err
    S_shift=shift

    WRITE(6,'("PAWiter",i5,1p,1e20.12,1p,2e15.7)')fcount,en

    call mkname(fcount,stuff)
    open(1001,file='pseudo.'//TRIM(stuff),form='formatted')
    do i=1,n
       write(1001,'(1p,15e15.7)') Gridwk%r(i),tw(i),shift(i),&
&                 (S_tgvv(i,ip),ip=1,nocc),&
&                 (S_tgvv(i,ip)*PAW%tphi(i,S_bmap(ip)),ip=1,nocc)
    enddo
    close(1001)

    fcount=fcount+1
    deallocate(shift,tg,tu,rhs,tw,dvx)

  END SUBROUTINE tVXsub0

!!!!!!!
!  tVXsub1
!!!!!!!
  SUBROUTINE tVXsub1(w,en,grad,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: en
    REAL(8), INTENT(OUT) :: grad(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io,l,iter,last,irc,ip,nocc,nbase,many,ib,ic,jo,p
    REAL(8) :: x,y,energy,occ,boundaryv
    REAL(8), POINTER :: rv(:)
    REAL(8),ALLOCATABLE :: shift(:),dvx(:,:),tg(:,:),tu(:),rhs(:,:),tw(:)
    INTEGER:: fcount=0
    REAL(8), PARAMETER :: tol=1.d-10, wgt0=0.01d0
    CHARACTER(4) :: stuff

    success=.TRUE.
    irc=S_irc; nocc=S_nocc; n=S_n; nbase=PAW%nbase; p=irc+1

    rv=>PAW%rVeff

    allocate(shift(n),tg(n,nocc),tu(n),rhs(n,nbase+1),tw(n),&
&    dvx(nbase,nocc))
    shift=0;grad=0

   !  Matrix element of \tilde{V}_x = w (input)
    do ip=1,nocc
       dvx(:,ip)=0; l=PAW%l(S_bmap(ip))  ;
       do ib=1,nbase
          if (l==PAW%l(ib)) THEN
          tu(1:irc)=w(1:irc)*PAW%otphi(1:irc,ib)*PAW%otphi(1:irc,S_bmap(ip))
          tu(1)=0; tu(2:irc)=tu(2:irc)/Gridwk%r(2:irc)
          dvx(ib,ip)=S_Vx(ib,ip)-integrator(Gridwk,tu,1,irc)
          endif
          write(6,'(" dvx ",2i5,1p,1e15.7)') ib,ip,dvx(ib,ip)
       enddo
    enddo
    call flush(6)
    tw=S_tVxvale
    tw(1:irc)=w(1:irc)
    shift=0 ;
    do ip=1,nocc
       occ=PAW%occ(S_bmap(ip))
       energy=S_E(ip)
       l=S_l(ip)
       rhs=0;
       rhs(2:n,1)=S_dExdtphi(2:n,ip)  &
&            -tw(2:n)*PAW%otphi(2:n,S_bmap(ip))/Gridwk%r(2:n)
       rhs(2:n,1)=rhs(2:n,1)-S_UOphi(2:n,ip)

       many=1
       do ib=1,nbase
          if (l==PAW%l(ib)) then
             rhs(:,1)=rhs(:,1)-dvx(ib,ip)*PAW%otp(:,ib)
             many=many+1
             rhs(:,many)=PAW%otp(:,ib)
          endif
       enddo

      S_tgvv(:,ip)=0
      If (many<=2) then
         call inhomo_numerov_SVD_bv(Gridwk,l,n-1,energy,tol,rv,rhs(:,1),S_tgvv(:,ip))
      else
         call  multisolv(S_bmap(ip),l,many,energy,rv,rhs(:,1:many),S_tgvv(:,ip))
      endif

      ! Remove homogeneous solution
      x=0; y=0
      write(6,*) 'Fix tail ', ip,S_bmap(ip)
      do i=irc+1,n
         x=x+PAW%otphi(i,S_bmap(ip))**2
         y=y+(S_tgvv(i,ip)-S_gvv(i,ip))*PAW%otphi(i,S_bmap(ip))
      enddo
         write(6,*) 'Adjust tail value ', ip ,y,x
         S_tgvv(:,ip)=S_tgvv(:,ip)-(y/x)*PAW%otphi(:,S_bmap(ip))

      do ib=1,nbase
         if (l==PAW%l(ib)) then
            write(6,*) ' ip ib <g|p> = ', ip,ib ,&
&             overlap(Gridwk,S_tgvv(:,ip),PAW%otp(:,ib))
            write(6,*) ' ip ib <g|tphi> = ', &
&             overlap(Gridwk,S_tgvv(:,ip),PAW%otphi(:,ib))
         endif
      enddo
   enddo   !ip

   shift=0;tu=0
   do ip=1,nocc
      occ=PAW%occ(S_bmap(ip))
      shift(:)=shift(:)+2*occ*S_tgvv(:,ip)*PAW%otphi(:,S_bmap(ip))
   enddo

    grad=0
    grad(2:irc)=-shift(2:irc)/Gridwk%r(2:irc)
    err=DOT_PRODUCT(grad(1:irc),grad(1:irc))
    en=err
    S_shift=shift

    WRITE(6,'("PAWiter",i5,1p,1e20.12,1p,2e15.7)')fcount,en

    call mkname(fcount,stuff)
    open(1001,file='pseudo.'//TRIM(stuff),form='formatted')
    do i=1,n
       write(1001,'(1p,15e15.7)') Gridwk%r(i),tw(i),shift(i),&
&         (S_gvv(i,ip),S_tgvv(i,ip),PAW%otphi(i,S_bmap(ip)),ip=1,nocc)
    enddo
    close(1001)

    fcount=fcount+1
    deallocate(shift,tg,tu,rhs,tw,dvx)

  END SUBROUTINE tVXsub1

   SUBROUTINE FindVX(mxloop,n,threshold,arg,trvx)
      Integer, INTENT(IN) :: mxloop,n
      REAL(8), INTENT(IN) :: threshold
      REAL(8), INTENT(INOUT) :: arg(:)
      REAL(8), INTENT(INOUT) :: trvx(:)

      LOGICAL :: success
      REAL(8) :: xxx
      INTEGER :: loop,i
      INTEGER, parameter :: mxl=50
      REAL(8), parameter :: mix=0.2d0
      CHARACTER(4) :: nm
      trvx=arg
      Do loop=1,mxloop*10
         CALL InitAnderson_dr(AC,6,2,n,0.001d0,1.d1,mxl,1.d-6,1.d-6,.true.)
         CALL DoAndersonMix(AC,arg,xxx,tVXsub0,success)
         write(6,*) 'Finished EXX_PseudoVx',AC%CurIter,AC%res
         CALL FreeAnderson(AC)

         trvx=arg
         arg=0
         arg(2:n)=-S_shift(2:n)/Gridwk%r(2:n)
         xxx=DOT_PRODUCT(arg(1:n),arg(1:n))
         write(6,*) 'tVX loop ', loop,xxx
           call mkname(loop,nm)
              Open(1001,file='anal'//TRIM(nm),form='formatted')
                 do i=1,n
                    write(1001,'(1p,20e15.7)') Gridwk%r(i),trvx(i),S_shift(i),&
&                            trvx(i)+mix*arg(i)
                 enddo
              close(1001)
         if (xxx<threshold) then
            write(6,*) ' tVX loop converged with ', loop,xxx
            exit
         endif
         trvx=trvx+mix*arg
         arg=trvx
      enddo

   END SUBROUTINE FindVX


  SUBROUTINE multisolv(in,l,mult,energy,rv,rhs,res)
    INTEGER, INTENT(IN):: in,l,mult
    REAL(8), INTENT(IN):: energy,rv(:),rhs(:,:)
    REAL(8), INTENT(INOUT):: res(:)

    integer :: i,io,jo,ib,many,n,nbase
    integer, ALLOCATABLE :: map(:)
    REAL(8), ALLOCATABLE :: tu(:),tw(:,:),M(:,:),MM(:,:),MMM(:,:),c(:),cc(:)
    REAL(8), PARAMETER :: tol=1.d-10

    n=S_n; res=0
    nbase=PAW%nbase
    allocate(map(nbase),M(nbase,nbase),MM(nbase,nbase),MMM(nbase,nbase),&
&      tw(n,nbase+1),tu(n),c(nbase),cc(nbase))

    many=0; map=0
    do ib=1,nbase
       if (l==PAW%l(ib)) then
          many=many+1; map(many)=ib;
       endif
    enddo

    if (mult-1/=many) then
         write(6,*) 'Error in multisolv ', mult,many
         stop
    endif

    call inhomo_numerov_SVD_bvm(Gridwk,l,n-1,mult,energy,tol,rv,rhs(:,1:mult),&
&             tw(:,1:mult))

    MMM=0; M=0
    do io=1,many
       do jo=1,many
          tu(:)=PAW%otp(:,map(io))*tw(:,jo+1)
          MMM(io,jo)=integrator(Gridwk,tu)
          M(io,jo)=PAW%dij(map(io),map(jo))-&
&                energy*PAW%oij(map(io),map(jo))
          enddo
       enddo
       MM=0
       do io=1,many
          do jo=1,many
             do ib=1,many
                MM(io,jo)=MM(io,jo)+M(io,ib)*MMM(ib,jo)
             enddo
          enddo
          MM(io,io)=MM(io,io)+1
      enddo

      c=0;cc=0
      do io=1,many
         tu(:)=PAW%otp(:,map(io))*tw(:,1)
         c(io)=integrator(Gridwk,tu)
      enddo
      do io=1,many
         do jo=1,many
            cc(io)=cc(io)-M(io,jo)*c(jo)
         enddo
      enddo

      MMM=MM
      call linsol(MMM,cc,many,nbase,nbase,nbase)

      write(6,*) ' linsol results '
      do io=1,many
         write(6,*) io,cc(io)
      enddo

      res=tw(:,1)
      do io=1,many
         res(:)=res(:)+cc(io)*tw(:,io+1)
      enddo

    deallocate(map,M,MM,MMM,tw,tu,c,cc)

  END SUBROUTINE multisolv

  SUBROUTINE Calc_w(Grid,PAW,io,w)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(OUT) :: w(:,:)     ! w(1:n,1:nbase)

    INTEGER :: i,j,l,n,nbase
    REAL(8) :: eig,zeroval,term
    REAL(8),POINTER :: rv(:),r(:)
    REAL(8) ,ALLOCATABLE :: rhs(:),temp(:)

    eig=PAW%eig(io)
    l=PAW%l(io)        !<-- lp
    rv=>PAW%rveff
    n=Grid%n
    r=>Grid%r
    ALLOCATE(rhs(n),temp(n))
    rhs =0
    w=0
    temp=0

    nbase=PAW%nbase
    DO i=1,nbase
       IF (l==PAW%l(i))THEN ! <-- lp = lk
          rhs = PAW%otp(:,i)
          CALL inhomo_bound_numerov(Grid,l,n,eig,rv,rhs,temp)
          w(:,i)=temp(:)
       ENDIF
    ENDDO

    OPEN(unit=2001,file='w.txt')
    DO i=1,n
       WRITE(2001,'(1p,50e15.7)') r(i), (w(i,j),j=1,nbase)
    ENDDO
    CLOSE(2001)

    DEALLOCATE(rhs,temp)
  END SUBROUTINE Calc_w

  SUBROUTINE Calc_u(Grid,PAW,io,rhs,u)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(IN) :: rhs(:)
    REAL(8), INTENT(OUT) :: u(:)

    INTEGER :: i,j,l,n
    REAL(8) :: eig,zeroval,term
    REAL(8),POINTER :: rv(:),r(:)
    REAL(8) ,ALLOCATABLE :: temp(:)

    eig=PAW%eig(io)
    l=PAW%l(io)
    rv=>PAW%rveff
    n=Grid%n
    r=>Grid%r

    ALLOCATE(temp(n))
    temp=0;u=0

    CALL inhomo_bound_numerov(Grid,l,n,eig,rv,rhs,u)

    OPEN(unit=2001,file='uio.txt')
    DO i=1,n
       WRITE(2001,'(1p,6e15.7)') r(i), u(i)
    ENDDO
    CLOSE(2001)

    Deallocate(temp)
  END SUBROUTINE Calc_u

END Module exx_pseudo
