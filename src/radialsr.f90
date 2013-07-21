!*****************************************************************
! Module for solving scalar relativistic radial equations
!   Uses program adapted by Marc Torrent and Francois Jollet from
!      USPS pgm of David Vanderbilt based on two coupled first order
!      differential equations
!      Previous version, based on second order differential equation
!        from formalism of Shadwick, Talman, and Norman, Comp. Phys. Comm.
!      54, 95-102 (1989)  found to be unstable
!   09-16-06  NAWH
!*****************************************************************
MODULE radialsr
  USE GlobalMath
  USE gridmod
  USE atomdata

  IMPLICIT NONE

  REAL(8), parameter :: inverse_fine_structure=137.03599911d0
  Real(8), private :: gamma,c1,c2,MA,MB
  Real(8), private, allocatable :: ww(:),jj(:)
     ! jj stores (r+(alpha/2)**2*(E*r-rv) == r*M(r)
     ! ww stores kappa*(kappa+1)/(r**2*M(r)) - (E - V(r))

CONTAINS

!******************************************************************
! Allocate_scalar_relativistic
!******************************************************************
 Subroutine Allocate_scalar_relativistic(Grid)
   Type(GridInfo), INTENT(IN) :: Grid

   INTEGER :: n,i

   n=Grid%n
   allocate(ww(n),jj(n), stat=i)
   if (i/=0) then
      write(6,*)  'Allocate_scalar_relativistic: error in allocation ',i,n
      stop
   endif

 End subroutine Allocate_scalar_relativistic

!******************************************************************
! Deallocate_scalar_relativistic
!******************************************************************
 Subroutine deallocate_scalar_relativistic

    deallocate(ww,jj)

 end subroutine deallocate_scalar_relativistic

!*******************************************************************
!  Subroutine Azeroexpand(Grid,Pot,l,energy)
!      If finitenucleus==.true. assumes potential for r--> has form
!         -2*Z*erf(r/RR)/r, where RR is a nuclear size parameter
!      Otherwise, assumes nuclear potential is -2*Z/r
!*******************************************************************
 Subroutine Azeroexpand(Grid,Pot,l,energy,nr)
   Type(GridInfo), INTENT(IN) :: Grid
   Type(PotentialInfo), INTENT(IN) :: Pot
   Integer, INTENT(IN) :: l
   Real(8), INTENT(IN) :: energy
   Integer, optional, INTENT(IN) :: nr

   Integer :: i,j,k,n
   Real(8) :: nz,xx,yy,angm,alpha2,balpha2

   n=Grid%n
   if (present(nr)) n=min(n,nr)

   nz=Pot%nz
    ww=0; jj=0;
   balpha2=inverse_fine_structure**2
   alpha2=1.d0/balpha2
   jj(1:n)=(Grid%r(1:n) + &
&       0.25d0*alpha2*(energy*Grid%r(1:n)-Pot%rv(1:n)))
   angm=l*(l+1)
   ww(2:n)=(Pot%rv(2:n)/Grid%r(2:n)-energy) &
&     + angm/(Grid%r(2:n)*jj(2:n))
   ww(1)=0

   if (.not.finitenucleus) then
      gamma=sqrt(angm+1.d0-alpha2*nz**2)
      c1=-((2.d0*balpha2)/(nz*(2*gamma+1)))*((1-gamma)*&
&          (1+0.25d0*alpha2*(energy-Pot%v0))+alpha2*nz**2+&
&          0.5d0*((nz*alpha2)**2)*(energy-Pot%v0))
      xx=2*nz*(1+alpha2*(energy-Pot%v0)+3*(0.25*alpha2*(energy-Pot%v0))**2)&
&      +0.25d0*alpha2*Pot%v0p*(alpha2*nz*nz+2*(gamma-1))
      yy=(1+gamma)*(1+0.25d0*alpha2*(energy-Pot%v0))+alpha2*nz*nz+&
&          0.5d0*((alpha2*nz)**2)*(energy-Pot%v0)
      c2=-(xx+yy*c1)/(2*(gamma+1)*alpha2*nz)
      !write(6,*) 'Azeroexpand: ', gamma,c1,c2
      MA=0; MB=0

   else  ! version for finite nuclear size
       gamma=0
       MA=1.d0+0.25d0*alpha2*(energy-Pot%v0)
       MB=0.25d0*alpha2*Pot%v0p
       c1=-MB*l/(2*MA*(l+1))
       c2=(-MA*MA*(energy-Pot%v0)*MB*c1*(l+1))/(MA*(4*l+6))
   endif

  end subroutine Azeroexpand

!*******************************************************************
! SUBROUTINE wfnsrinit(Grid,l,wfn,lwfn,istart)
!*******************************************************************
  SUBROUTINE wfnsrinit(Grid,l,wfn,lwfn,istart)
   ! returns the solution of the scalar relativistic equations near r=0
   !  using power series expansion
   Type(GridInfo), INTENT(IN) :: Grid
   INTEGER, INTENT(IN) :: l
    REAL(8),INTENT(INOUT) :: wfn(:),lwfn(:)
    INTEGER, INTENT(OUT) :: istart

    REAL(8) :: rr,M
    INTEGER :: i,j,n

    wfn=0; lwfn=0

    istart=6
    do i=1,istart
       rr=Grid%r(i+1)
       if (.not.finitenucleus) then
          wfn(i+1)=1+rr*(c1+rr*c2)
          lwfn(i+1)=(gamma-1)+rr*(c1*gamma+rr*c2*(gamma+1))
          wfn(i+1)=wfn(i+1)*(rr**gamma)
          lwfn(i+1)=lwfn(i+1)*(rr**gamma)/jj(i+1)

       else   ! finite nucleus case
          M=MA-MB*rr
          wfn(i+1)=(1+rr*(c1+rr*c2))*(rr**(l+1))
          lwfn(i+1)=(l+rr*((l+1)*c1+rr*(l+2)*c2))*(rr**(l+1))/M
       endif

    enddo

  End SUBROUTINE wfnsrinit

 subroutine wfnsrasym(Grid,wfn,lwfn,energy,iend)
  ! returns the solution of the scalar relativistic equations near r=inf
  !  using exp(-x*r) for upper component
   Type(GridInfo), INTENT(IN) :: Grid
    REAL(8),INTENT(INOUT) :: wfn(:),lwfn(:)
    REAL(8), INTENT(IN) :: energy
    INTEGER, INTENT(OUT) :: iend

    REAL(8) :: rr,x,m
    INTEGER :: i,j,n

    if (energy>0.d0) then
       write(6,*) 'Error in wfnsrasym -- energy > 0', energy
       stop
    endif

    wfn=0; lwfn=0
    n=Grid%n

    m=1+0.25d0*energy/(inverse_fine_structure**2)
    x=sqrt(-m*energy)
    !write(6,*) ' in wfnsrasym with x = ',x
    iend=5
    do i=n-iend,n
       wfn(i)=exp(-x*(Grid%r(i)-Grid%r(n-iend)))
       lwfn(i)=-wfn(i)*(x+1.d0/Grid%r(i))/m
    enddo
  end subroutine wfnsrasym

  !**********************************************************************
  !      subroutine unboundsr(Grid,Pot,nr,l,energy,wfn,nodes)
  !  pgm to solve radial scalar relativistic equation for unbound states
  !    at energy 'energy' and at angular momentum l
  !
  !    with potential rv/r, given in uniform linear or log mesh of n points
  !   assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
  !
  !  nz=nuclear charge
  !
  !  Does not use Noumerov algorithm -- but uses coupled first-order
  !       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
  !
  ! also returns node == number of nodes for calculated state
  !************************************************************************
  SUBROUTINE unboundsr(Grid,Pot,nr,l,energy,wfn,nodes)
    TYPE(GridInfo), INTENT(IN)  :: Grid
    TYPE(PotentialInfo), INTENT(IN)  :: Pot
    INTEGER, INTENT(IN) :: nr,l
    REAL(8), INTENT(IN) :: energy
    REAL(8), INTENT(INOUT) :: wfn(:)
    INTEGER, INTENT(INOUT) :: nodes

    INTEGER :: n,i,j,k,ierr,istart
    REAL(8) :: scale
    REAL(8), allocatable :: lwfn(:),zz(:,:,:),yy(:,:)

    n=Grid%n
    IF (nr > n) THEN
       WRITE(6,*) 'Error in unboundsr -- nr > n', nr,n
       STOP
    ENDIF

    call Azeroexpand(Grid,Pot,l,energy,nr)

    allocate(lwfn(nr),zz(2,2,nr),yy(2,nr),stat=ierr)
       if (ierr/=0) then
          write(6,*) ' allocation error in unboundsr ', nr,ierr
          stop
       endif

    lwfn=0;zz=0;yy=0;

    call wfnsrinit(Grid,l,wfn,lwfn,istart)
    call prepareforcfdsol(Grid,1,istart,nr,wfn,lwfn,yy,zz)
    call cfdsol(Grid,zz,yy,istart,nr)
    call getwfnfromcfdsol(1,nr,yy,wfn)
    nodes=countnodes(2,nr,wfn)
    !
    ! normalize to unity within integration range
    !
    !call filter(nr,wfn,machine_zero)     !NH votes to remove
    scale=1.d0/overlap(Grid,wfn(1:nr),wfn(1:nr),1,nr)
    scale=SIGN(SQRT(scale),wfn(nr-2))
    wfn(1:nr)=wfn(1:nr)*scale

    deallocate(lwfn,yy,zz)

  END SUBROUTINE unboundsr

!******************************************************************
!  SUBROUTINE boundsr(Grid,Pot,eig,wfn,l,nroot,emin,ierr,success)
!******************************************************************
  SUBROUTINE boundsr(Grid,Pot,eig,wfn,l,nroot,emin,ierr,success)
    !  pgm to solve radial scalar relativistic equation for nroot bound state
    !    energies and wavefunctions for angular momentum l
    !    with potential rv/r, given in uniform linear or log mesh of n points
    !  nz=nuclear charge
    !  emin=is estimate of lowest eigenvalue; used if nz=0
    !     otherwise, set to the value of -(nz/(l+1))**2
    !
    !  It is assumed that the wavefunction has np-l-1 nodes, where
    !    np is the principle quantum number-- np=1,2,..nroot
    !
    !  Does not use Noumerov algorithm -- but uses coupled first-order
    !       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
    !
    !  Corrections are also needed for r>n*h, depending on:
    !         e0 (current guess of energy eigenvalue
    !         the extrapolated value of rv == r * v
    !
    ! ierr=an nroot digit number indicating status of each root
    !   a digit of 1 indicates success in converging root
    !              2 indicates near success in converging root
    !              9 indicates that root not found
    !
    ! first check how many roots expected =  ntroot (returned as argument)
    !
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    REAL(8), INTENT(INOUT) :: eig(:),wfn(:,:)
    INTEGER, INTENT(IN) :: l,nroot
    INTEGER, INTENT(INOUT) :: ierr
    REAL(8), INTENT(INOUT) :: emin
    LOGICAL, INTENT(INOUT) :: success

    REAL(8), PARAMETER :: convre=1.d-10,vlrg=1.d30
    INTEGER, PARAMETER :: niter=1000

    REAL(8), POINTER :: rv(:)
    REAL(8), ALLOCATABLE :: p1(:),p2(:),dd(:)
    INTEGER :: n
    REAL(8) :: nz,h,v0,v0p
    REAL(8) :: err,convrez,energy,zeroval
    REAL(8) :: scale,emax,best,rout
    REAL(8) :: arg,r,r2,veff,pppp1,rin,dele,x,rvp1,pnp1,bnp1
    INTEGER :: iter,i,j,k,node,match,mxroot,ntroot,ir,iroot
    INTEGER :: least,many,ifac,istart,iend
    LOGICAL :: ok
    REAL(8), allocatable :: lwfn(:),zz(:,:,:),yy(:,:)
    ! integer :: icount=0

    n=Grid%n
    h=Grid%h
    ALLOCATE(p1(n),p2(n),dd(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) ' Error in boundsr allocation ',i,n
       STOP
    ENDIF

    success=.true.
    allocate(lwfn(n),zz(2,2,n),yy(2,n),stat=i)
       if (i/=0) then
          write(6,*) ' allocation error in boundsr ', n,i
          stop
       endif


    nz=Pot%nz
    v0=Pot%v0
    v0p=Pot%v0p
    rv=>Pot%rv
    err=n*nz*(h**4)
    convrez=convre
    IF (nz>0.001d0) convrez=convre*nz
    ierr=0

    WRITE(6,*) 'z , l = ',nz,l
    ! check how many roots expected by integration outward at
    !   energy = 0
    energy = 0
    call Azeroexpand(Grid,Pot,l,energy)
    lwfn=0;zz=0;yy=0;
    call wfnsrinit(Grid,l,p1,lwfn,istart)
    !
    !start outward integration
    call prepareforcfdsol(Grid,1,istart,n,p1,lwfn,yy,zz)
    call cfdsol(Grid,zz,yy,istart,n)
    call getwfnfromcfdsol(1,n,yy,p1)
    node=countnodes(2,n,p1)

    WRITE(6,*) ' nodes at e=0  ', node

    mxroot=node+1
    ntroot=node
    IF (mxroot.LT.nroot) THEN
       WRITE(6,*)'error in boundsr - for l = ',l
       WRITE(6,*) nroot,' states requested but only',mxroot,' possible'
       DO ir=mxroot+1,nroot
          ierr=ierr+9*(10**(ir-1))
       ENDDO
       success=.false.
    ENDIF
    mxroot=min0(mxroot,nroot)
    !
    IF (nz.EQ.0) energy=-ABS(emin)
    IF (nz.NE.0) energy=-1.1d0*(nz/(l+1.d0))**2
    emin=energy-err
    emax=0.d0

    DO iroot=1,mxroot
       best=1.d10; dele=1.d10
       energy=emin+err
       IF (energy.LT.emin) energy=emin
       IF (energy.GT.emax) energy=emax
       ok=.FALSE.
       !write(6,*) 'iter,iroot,energy',iter,iroot,energy
       !write(6,*) 'emin,max',emin,emax
       BigIter: DO iter=1,niter
          !write(6,*) 'In iter with energy', iter,energy,niter,l,iroot
          !  start inward integration
          !  start integration at n
          call Azeroexpand(Grid,Pot,l,energy)
          ! find classical turning point
          call ClassicalTurningPoint(Grid,Pot%rv,l,energy,match)
          match=max(match,10); match=min(match,n-20)
          call wfnsrasym(Grid,p2,lwfn,energy,iend)
          call prepareforcfdsol(Grid,n-iend,n,n,p2,lwfn,yy,zz)
          call cfdsol(Grid,zz,yy,n-iend,match)
          call getwfnfromcfdsol(match,n,yy,p2)
          match=match+6
          rin=Gfirstderiv(Grid,match,p2)/p2(match)

          call wfnsrinit(Grid,l,p1,lwfn,istart)
          call prepareforcfdsol(Grid,1,istart,n,p1,lwfn,yy,zz)
          call cfdsol(Grid,zz,yy,istart,match+6)
          call getwfnfromcfdsol(1,match+6,yy,p1)
          node= countnodes(2,match+6,p1)

          !icount=icount+1
          !  do i=1,match+6
          !    write(100+icount,'(1p,2e15.7)') Grid%r(i),p1(i)
          !  enddo
          rout=Gfirstderiv(Grid,match,p1)/p1(match)
          ! check whether node = (iroot-1)
          !   not enough nodes -- raise energy
          IF (node.LT.iroot-1) THEN
             emin=MAX(emin,energy)-err
             energy=emax-(emax-energy)*ranx()
             ifac=9
             !   too many nodes -- lower energy
          ELSEIF (node.GT.iroot-1) THEN
             IF (energy.LE.emin) THEN
                ierr=ierr+9*(10**(iroot-1))
                WRITE(6,*) 'boundsr error -- emin too high',l,nz,emin,energy
                do i=2,n
                   write(999,'(1p,4e15.7)') Grid%r(i),jj(i)/Grid%r(i),ww(i),Pot%rv(i)
                enddo
                STOP
             ENDIF
             emax=MIN(emax,energy+err)
             energy=emin+(energy-emin)*ranx()
             !   correct number of nodes -- estimate correction
          ELSEIF (node.EQ.iroot-1) THEN
             DO j=1,match
                p1(j)=p1(j)/p1(match)
                         !write(6,*) 'j,p1',j,p1(j)
             ENDDO
             DO j=match,n
                p1(j)=p2(j)/p2(match)
                         !write(6,*) 'j,p2',j,p1(j)
             ENDDO
             scale=1.d0/overlap(Grid,p1,p1)
             dele=(rout-rin)*scale
                  !write(6,*) 'energy,dele,scale',energy,dele,scale
             x=ABS(dele)
             IF (x.LT.best) THEN
                scale=SQRT(scale)
                p1(1:n)=p1(1:n)*scale
                call filter(n,p1,machine_zero)
                wfn(1:n,iroot)=p1(1:n)
                eig(iroot)=energy
                !write(6,*) 'root',l,iroot,eig(iroot),emin,emax
                best=x
             ENDIF
             IF (ABS(dele).LE.convrez) THEN
                !write(6,*) 'iter with dele' , iter,dele
                ok=.TRUE.
                !  eigenvalue found
                ierr=ierr+10**(iroot-1)
                IF (iroot+1.LE.mxroot) THEN
                   emin=energy+err
                   emax=0
                   energy=(emin+emax)/2
                   IF (energy.LT.emin) energy=emin
                   IF (energy.GT.emax) energy=emax
                   best=1.d10
                ENDIF
                EXIT BigIter
             ENDIF
             IF (ABS(dele).GT.convrez) THEN
                energy=energy+dele
                ! if energy is out of range, pick random energy in correct range
                IF (emin-energy.GT.convrez.OR.energy-emax.GT.convrez)         &
&                    energy=emin+(emax-emin)*ranx()
                ifac=2
             ENDIF
          ENDIF
       ENDDO BigIter !iter
       IF (.NOT.ok) THEN
          success=.false.     
          ierr=ierr+ifac*(10**(iroot-1))
          WRITE(6,*) 'no convergence in boundsr',iroot,l,dele,energy
          WRITE(6,*) ' best guess of eig, dele = ',eig(iroot),best
          IF (iroot.LT.mxroot) THEN
             DO ir=iroot+1,mxroot
                ierr=ierr+9*(10**(ir-1))
             ENDDO
          ENDIF
        ! reset wfn with hydrogenic form
        j=iroot+l+1
        wfn(:,iroot)=0
        x=(j)*sqrt(abs(eig(iroot)))
        do i=2,n
           wfn(i,iroot)=hwfn(x,j,l,Grid%r(i))
        enddo
       ENDIF
    ENDDO !iroot

    DEALLOCATE(p1,p2,dd,lwfn,yy,zz)
             write(6,*) 'returning from boundsr -- ierr=',ierr
  END SUBROUTINE Boundsr

  subroutine scalarrelativisticturningpt(Grid,least,turningpoint)
     Type(GridInfo), INTENT(IN) :: Grid
     Integer, INTENT(IN) :: least
     Integer, INTENT(OUT) :: turningpoint

    integer :: i,n

    n=Grid%n

    turningpoint=n
    do i=n,least,-1
      if (ww(i)<0.d0) exit
    enddo
    turningpoint=i

    !write(6,*) 'Found turning point at ', turningpoint, Grid%r(turningpoint)

  End subroutine scalarrelativisticturningpt

  subroutine prepareforcfdsol(Grid,i1,i2,n,wfn,lwfn,yy,zz)
     Type(gridinfo), INTENT(IN) :: Grid
     INTEGER, INTENT(IN) :: i1,i2,n
     REAL(8), INTENT(IN) :: wfn(:),lwfn(:)
     REAL(8), INTENT(OUT) :: yy(:,:),zz(:,:,:)

     INTEGER :: i

      yy=0;zz=0
      yy(1,i1:i2)=wfn(i1:i2)
      yy(2,i1:i2)=lwfn(i1:i2)

      do  i=2,n
       zz(1,1,i)=1.d0/Grid%r(i)
       zz(1,2,i)=jj(i)/Grid%r(i)
       zz(2,2,i)=-1.d0/Grid%r(i)
       zz(2,1,i)=ww(i)
      enddo

   end subroutine prepareforcfdsol

   subroutine getwfnfromcfdsol(start,finish,yy,wfn)
      INTEGER, INTENT(IN) :: start,finish
      REAL(8), INTENT(IN) :: yy(:,:)
      REAL(8), INTENT(INOUT) :: wfn(:)

      INTEGER :: i

      wfn=0
      do i=start,finish
         wfn(i)=yy(1,i)
      enddo
   end subroutine getwfnfromcfdsol


END MODULE radialsr
