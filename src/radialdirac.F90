!*****************************************************************
! Module for solving Dirac  relativistic radial equations
!   Uses program adapted by Marc Torrent and Francois Jollet from
!      USPS pgm of David Vanderbilt based on two coupled first order
!      differential equations
!      internally calulates G(r) -->wfn  and cF(r); at return,
!          lwfn = cF/c
!    Note:   at the moment, does not work for finite nuclear models
!      Uses structures from radialsr.F90 code  -- may need to 
!        modify Dzeroexpand
!   01-07-18  NAWH
!*****************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains the following active subroutines:
!     Allocate_Dirac_relativistic, deallocate_Dirac_relativistic
!        Dzeroexpand, wfnDinit, wfnDasym, unboundD, boundD,
!        prepareforcfdsolD, getwfnfromcfdsolD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE radialdirac

  USE io_tools
  USE GlobalMath
  USE gridmod
  USE atomdata

  IMPLICIT NONE

  !REAL(8), parameter :: inverse_fine_structure=137.035999139d0
  !  moved to globalmath
  Real(8), private :: s,A0,B0,A1,B1
  Real(8), private, allocatable :: ww(:),jj(:)
     ! ww stores  (E - V(r))
     ! jj stores  1 + 0.25d0*alpha**2*(E - V(r))

CONTAINS

!******************************************************************
! Allocate_Dirac_relativistic
!******************************************************************
 Subroutine Allocate_Dirac_relativistic(Grid)
   Type(GridInfo), INTENT(IN) :: Grid

   INTEGER :: n,i

   n=Grid%n
   allocate(ww(n),jj(n), stat=i)
   if (i/=0) then
      write(std_out,*)  'Allocate_Dirac_relativistic: error in allocation ',i,n
      stop
   endif

 End subroutine Allocate_Dirac_relativistic

!******************************************************************
! Deallocate_Dirac_relativistic
!******************************************************************
 Subroutine deallocate_Dirac_relativistic

    deallocate(ww,jj)

 end subroutine deallocate_Dirac_relativistic

!*******************************************************************
!  Subroutine Dzeroexpand(Grid,Pot,kappa,energy)
!      If finitenucleus==.true. assumes potential is non-singular
!          at origin and Pot%v0 and Pot%v0p are properly set
!          -- actually not programmed yet
!      Otherwise, assumes nuclear potential is -2*Z/r
!*******************************************************************
 Subroutine Dzeroexpand(Grid,Pot,kappa,energy,nr)
   Type(GridInfo), INTENT(IN) :: Grid
   Type(PotentialInfo), INTENT(IN) :: Pot
   Integer, INTENT(IN) :: kappa
   Real(8), INTENT(IN) :: energy
   Integer, optional, INTENT(IN) :: nr

   Integer :: i,j,k,n
   Real(8) :: nz,xx,yy,angm,alpha2,balpha2
   Real(8) :: x,y,z

   n=Grid%n
   if (present(nr)) n=min(n,nr)

   nz=Pot%nz
    ww=0; jj=0;
   balpha2=inverse_fine_structure**2
   alpha2=1.d0/balpha2
   ww(2:n)=energy-Pot%rv(2:n)/Grid%r(2:n)
   jj(2:n)=(1.d0 + 0.25d0*alpha2*ww(2:n))

   if (.not.finitenucleus) then
      s=sqrt(kappa*kappa-alpha2*nz**2)
      A0=1.d0
      B0=2.d0*(s+kappa)*balpha2/nz
      z=2*s+1
      x=alpha2*(nz**2)
      y=4*alpha2+energy-Pot%v0
      A1=(4*alpha2*x+y*(s+kappa-2*x))/(2*nz*z)
      y=2*alpha2+energy-Pot%v0
      B1=-(y*(2*(s+kappa)+energy-Pot%v0))/z

   else  ! version for finite nuclear size
           write(std_out,*) 'Dirac case not yet programmed for finite nucleus'
           stop
   endif        

  end subroutine Dzeroexpand

!*******************************************************************
! SUBROUTINE wfnDinit(Grid,kappa,wfn,lwfn,istart)
!*******************************************************************
  SUBROUTINE wfnDinit(Grid,kappa,wfn,lwfn,istart)
   ! returns the solution of the Dirac relativistic equations near r=0
   !  using power series expansion
   Type(GridInfo), INTENT(IN) :: Grid
   INTEGER, INTENT(IN) :: kappa
    REAL(8),INTENT(INOUT) :: wfn(:),lwfn(:)
    INTEGER, INTENT(OUT) :: istart

    REAL(8) :: rr,M
    INTEGER :: i,j,n

    wfn=0; lwfn=0

    !write(std_out,*) 'Entering wfnDinit with s =', kappa,s
    istart=6
    do i=1,istart
       rr=Grid%r(i+1)
       if (.not.finitenucleus) then
          wfn(i+1)=(rr**s)*(A0+A1*rr)
          lwfn(i+1)=(rr**s)*(B0+B1*rr)

       else   ! finite nucleus case
               write(std_out,*) 'Dirac case not programmed for finite nucleus'
               STOP
       endif

    enddo

  End SUBROUTINE wfnDinit

 subroutine wfnDasym(Grid,wfn,lwfn,energy,iend)
  ! returns the solution of the Dirac relativistic equations near r=inf
  !  using exp(-x*r) for upper component
   Type(GridInfo), INTENT(IN) :: Grid
    REAL(8),INTENT(INOUT) :: wfn(:),lwfn(:)
    REAL(8), INTENT(IN) :: energy
    INTEGER, INTENT(OUT) :: iend

    REAL(8) :: rr,x,m,qx
    INTEGER :: i,j,n

    if (energy>0.d0) then
       write(std_out,*) 'Error in wfnDasym -- energy > 0', energy
       stop
    endif

    wfn=0; lwfn=0; 
    n=Grid%n


    m=1.d0+0.25d0*energy/(inverse_fine_structure**2)
    x=sqrt(-m*energy)
    rr=energy/x
    iend=5
    do i=n-iend,n
       wfn(i)=exp(-x*(Grid%r(i)-Grid%r(n-iend)))
       lwfn(i)=rr*wfn(i)
    enddo
  end subroutine wfnDasym

  !**********************************************************************
  !      subroutine unboundD(Grid,Pot,nr,kappa,energy,wfn,lwfn,nodes)
  !  pgm to solve radial Dirac relativistic equation for unbound states
  !    at energy 'energy' and at spin-orbit parameter kappa
  !
  !    with potential rv/r, given in uniform linear or log mesh of n points
  !
  !  nz=nuclear charge
  !
  !  Does not use Noumerov algorithm -- but uses coupled first-order
  !       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
  !
  ! also returns node == number of nodes for calculated state
  !************************************************************************
  SUBROUTINE unboundD(Grid,Pot,nr,kappa,energy,wfn,lwfn,nodes)
    TYPE(GridInfo), INTENT(IN)  :: Grid
    TYPE(PotentialInfo), INTENT(IN)  :: Pot
    INTEGER, INTENT(IN) :: nr,kappa
    REAL(8), INTENT(IN) :: energy
    REAL(8), INTENT(INOUT) :: wfn(:),lwfn(:)
    INTEGER, INTENT(INOUT) :: nodes

    INTEGER :: n,i,j,k,ierr,istart
    REAL(8) :: scale
    REAL(8), allocatable :: zz(:,:,:),yy(:,:)

    n=Grid%n
    IF (nr > n) THEN
       WRITE(STD_OUT,*) 'Error in unboundD -- nr > n', nr,n
       STOP
    ENDIF

    call Dzeroexpand(Grid,Pot,kappa,energy,nr)

    allocate(zz(2,2,nr),yy(2,nr),stat=ierr)
       if (ierr/=0) then
          write(std_out,*) ' allocation error in unboundD ', nr,ierr
          stop
       endif

    wfn=0;lwfn=0;zz=0;yy=0;

    call wfnDinit(Grid,kappa,wfn,lwfn,istart)
    call prepareforcfdsolD(Grid,1,istart,nr,kappa,wfn,lwfn,yy,zz)
    call cfdsol(Grid,zz,yy,istart,nr)
    call getwfnfromcfdsolD(1,nr,yy,wfn,lwfn)
    nodes=countnodes(2,nr,wfn)
    !
    ! normalize to unity within integration range
    !
    ! change back to original lower component
    scale=0.5d0/inverse_fine_structure
    lwfn(1:nr)=scale*lwfn(1:nr)
    scale=overlap(Grid,wfn(1:nr),wfn(1:nr),1,nr)
    scale=scale+overlap(Grid,lwfn(1:nr),lwfn(1:nr),1,nr)
    scale=SIGN(SQRT(1.d0/scale),wfn(nr-2))
    wfn(1:nr)=wfn(1:nr)*scale
    lwfn(1:nr)=lwfn(1:nr)*scale

    deallocate(yy,zz)

  END SUBROUTINE unboundD

!******************************************************************
!  SUBROUTINE boundD(Grid,Pot,eig,wfn,lwfn,kappa,nroot,emin,ierr,success)
!******************************************************************
  SUBROUTINE boundD(Grid,Pot,eig,wfn,lwfn,kappa,nroot,emin,ierr,success)
    !  pgm to solve radial Dirac relativistic equation for nroot bound state
    !    energies and wavefunctions for spin orbit parameter kappa
    !    with potential rv/r, given in uniform linear or log mesh of n points
    !  nz=nuclear charge
    !  emin=is estimate of lowest eigenvalue; used if nz=0
    !     otherwise, set to the value of -(nz/(l+1))**2
    !
    !  It is assumed that the wavefunction has np-l-1 nodes, where
    !    np is the principal quantum number-- np=1,2,..nroot
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
    REAL(8), INTENT(INOUT) :: eig(:),wfn(:,:),lwfn(:,:)
    INTEGER, INTENT(IN) :: kappa,nroot
    INTEGER, INTENT(INOUT) :: ierr
    REAL(8), INTENT(INOUT) :: emin
    LOGICAL, INTENT(INOUT) :: success

    REAL(8), PARAMETER :: convre=1.d-10,vlrg=1.d30
    REAL(8), PARAMETER :: ftr=0.5d0/inverse_fine_structure
    INTEGER, PARAMETER :: niter=1000

    REAL(8), POINTER :: rv(:)
    REAL(8), ALLOCATABLE :: p1(:),lp1(:),p2(:),lp2(:),dd(:)
    INTEGER :: n
    REAL(8) :: nz,h,v0,v0p
    REAL(8) :: err,convrez,energy,zeroval
    REAL(8) :: scale,emax,best,rout
    REAL(8) :: arg,r,r2,veff,pppp1,rin,dele,x,rvp1,pnp1,bnp1
    INTEGER :: iter,i,j,k,node,match,mxroot,ntroot,ir,iroot,l
    INTEGER :: least,many,ifac,istart,iend
    LOGICAL :: ok
    REAL(8), allocatable :: zz(:,:,:),yy(:,:)
    ! integer :: icount=0

    n=Grid%n
    h=Grid%h
    if (kappa<0)  l=-kappa-1
    if (kappa>0)  l=kappa        

    ALLOCATE(p1(n),lp1(n),p2(n),lp2(n),dd(n),stat=i)
    IF (i/=0) THEN
       WRITE(STD_OUT,*) ' Error in boundD allocation ',i,n
       STOP
    ENDIF

    success=.true.
    allocate(zz(2,2,n),yy(2,n),stat=i)
       if (i/=0) then
          write(std_out,*) ' allocation error in boundD ', n,i
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

    WRITE(STD_OUT,*) 'z ,kappa, l = ',nz,kappa,l
    ! check how many roots expected by integration outward at
    !   energy = 0
    energy = 0
    call Dzeroexpand(Grid,Pot,kappa,energy)
    zz=0;yy=0;
    call wfnDinit(Grid,l,p1,lp1,istart)
    !
    !start outward integration
    call prepareforcfdsolD(Grid,1,istart,n,kappa,p1,lp1,yy,zz)
    call cfdsol(Grid,zz,yy,istart,n)
    call getwfnfromcfdsolD(1,n,yy,p1,lp1)
    node=countnodes(2,n,p1)

    WRITE(STD_OUT,*) ' nodes at e=0  ', node

    mxroot=node+1
    ntroot=node
    IF (mxroot.LT.nroot) THEN
       WRITE(STD_OUT,*)'error in boundD - for l = ',l
       WRITE(STD_OUT,*) nroot,' states requested but only',mxroot,' possible'
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
       !write(std_out,*) 'iter,iroot,energy',iter,iroot,energy
       !write(std_out,*) 'emin,max',emin,emax
       BigIter: DO iter=1,niter
          !write(std_out,*) 'In iter with energy', iter,energy,niter,l,iroot
          !  start inward integration
          !  start integration at n
          call Dzeroexpand(Grid,Pot,kappa,energy)
          ! find classical turning point
          call ClassicalTurningPoint(Grid,Pot%rv,l,energy,match)
          match=max(match,10); match=min(match,n-20)
          call wfnDasym(Grid,p2,lp2,energy,iend)
          call prepareforcfdsolD(Grid,n-iend,n,n,kappa,p2,lp2,yy,zz)
          call cfdsol(Grid,zz,yy,n-iend,match)
          call getwfnfromcfdsolD(match,n,yy,p2,lp2)
          match=match+6
          rin=lp2(match)/p2(match)

          call wfnDinit(Grid,l,p1,lp1,istart)
          call prepareforcfdsolD(Grid,1,istart,n,kappa,p1,lp1,yy,zz)
          call cfdsol(Grid,zz,yy,istart,match+6)
          call getwfnfromcfdsolD(1,match+6,yy,p1,lp1)
          node= countnodes(2,match+6,p1)

          !icount=icount+1
          !  do i=1,match+6
          !    write(100+icount,'(1p,2e15.7)') Grid%r(i),p1(i)
          !  enddo
          rout=lp1(match)/p1(match)
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
                WRITE(STD_OUT,*) 'boundD error -- emin too high',l,nz,emin,energy
                do i=2,n
                   write(999,'(1p,4e15.7)') Grid%r(i),jj(i),ww(i),Pot%rv(i)
                enddo
                STOP
             ENDIF
             emax=MIN(emax,energy+err)
             energy=emin+(energy-emin)*ranx()
             !   correct number of nodes -- estimate correction
          ELSEIF (node.EQ.iroot-1) THEN
             DO j=1,match
                p1(j)=p1(j)/p1(match)
                lp1(j)=ftr*lp1(j)/p1(match)
                   !if (j>match-3)   write(std_out,*) 'j,p1',j,p1(j),lp1(j)
             ENDDO
             DO j=match,n
                p1(j)=p2(j)/p2(match)
                lp1(j)=ftr*lp2(j)/p2(match)
                  !if (j<match+3)  write(std_out,*) 'j,p2',j,p1(j),lp1(j)
             ENDDO
             scale=overlap(Grid,p1,p1)+overlap(Grid,lp1,lp1)
             dele=(rout-rin)/scale
                  !write(std_out,*) 'energy,dele,scale',energy,dele,scale
             x=ABS(dele)
             IF (x.LT.best) THEN
                scale=1.d0/SQRT(scale)
                p1(1:n)=p1(1:n)*scale
                lp1(1:n)=lp1(1:n)*scale
                call filter(n,p1,machine_zero)
                call filter(n,lp1,machine_zero)
                wfn(1:n,iroot)=p1(1:n)
                lwfn(1:n,iroot)=lp1(1:n)
                eig(iroot)=energy
                !write(std_out,*) 'root',l,iroot,eig(iroot),emin,emax
                best=x
             ENDIF
             IF (ABS(dele).LE.convrez) THEN
                !write(std_out,*) 'iter with dele' , iter,dele
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
          WRITE(STD_OUT,*) 'no convergence in boundD',iroot,l,dele,energy
          WRITE(STD_OUT,*) ' best guess of eig, dele = ',eig(iroot),best
          IF (iroot.LT.mxroot) THEN
             DO ir=iroot+1,mxroot
                ierr=ierr+9*(10**(ir-1))
             ENDDO
          ENDIF
        ! reset wfn with hydrogenic form
        j=iroot+l+1
        wfn(:,iroot)=0
        lwfn(:,iroot)=0
        x=(j)*sqrt(abs(eig(iroot)))
        do i=2,n
           call dirachwfn(j,kappa,x,Grid%r(i),arg,wfn(i,iroot),lwfn(i,iroot))
        enddo
       ENDIF
    ENDDO !iroot

    DEALLOCATE(p1,lp1,p2,lp2,dd,yy,zz)
             write(std_out,*) 'returning from boundD -- ierr=',ierr
  END SUBROUTINE BoundD

  subroutine prepareforcfdsolD(Grid,i1,i2,n,kappa,wfn,lwfn,yy,zz)
     Type(gridinfo), INTENT(IN) :: Grid
     INTEGER, INTENT(IN) :: i1,i2,n,kappa
     REAL(8), INTENT(IN) :: wfn(:),lwfn(:)
     REAL(8), INTENT(OUT) :: yy(:,:),zz(:,:,:)

     INTEGER :: i

      yy=0;zz=0
      yy(1,i1:i2)=wfn(i1:i2)
      yy(2,i1:i2)=lwfn(i1:i2)

      do  i=2,n
       zz(1,1,i)=-kappa/Grid%r(i)
       zz(1,2,i)=jj(i)
       zz(2,2,i)=kappa/Grid%r(i)
       zz(2,1,i)=-ww(i)
      enddo

   end subroutine prepareforcfdsolD

   subroutine getwfnfromcfdsolD(start,finish,yy,wfn,lwfn)
      INTEGER, INTENT(IN) :: start,finish
      REAL(8), INTENT(IN) :: yy(:,:)
      REAL(8), INTENT(INOUT) :: wfn(:),lwfn(:)

      INTEGER :: i

      wfn=0
      do i=start,finish
         wfn(i)=yy(1,i)
         lwfn(i)=yy(2,i)
      enddo
   end subroutine getwfnfromcfdsolD


END MODULE radialDirac
