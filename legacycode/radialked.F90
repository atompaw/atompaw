!*****************************************************************
!  Previous version of radialked -- no longer used
! Module for solving modified Kohn-Sham radial equations
!   For the case that the exchange-correlation functional depends
!     on the kinetic energy density
!   Uses program adapted by Marc Torrent and Francois Jollet from
!      USPS pgm of David Vanderbilt based on two coupled first order
!      differential equations
!   08-15-19  NAWH
!*****************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains the following active subroutines:
!        ALLOCATE_KED, DEALLOCATE_KED, Set_Pot,
!        wfnkedinit, wfnkedasym, unboundked, boundked,
!        ClassicalTurningked, setupforcfdsol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE radialked

  USE io_tools
  USE GlobalMath
  USE gridmod
  USE atomdata

  IMPLICIT NONE
    REAL(8), private, allocatable :: oneplusvtau(:),dvtaudr(:)
    REAL(8), private :: qq        ! ionic charge (qq=0 for neutral)
    REAL(8), private :: zxc       !rvx(1)

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate_ked
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine Allocate_ked(Grid)
     Type(GridInfo), INTENT(IN) :: Grid
     ALLOCATE( oneplusvtau(Grid%n),dvtaudr(Grid%n))
  End subroutine Allocate_ked

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEAllocate_ked
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine Deallocate_ked
     DEALLOCATE( oneplusvtau,dvtaudr)
 End subroutine DEAllocate_ked

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set_Pot
!    Version with polynomial fitting of vtau for 0<=r<=0.001
!       and corresponding reseting of oneplusvtau and dvtaudr in that range
!       NAWH   4/6/2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Subroutine Set_Pot(Grid,Pot)
     Type(GridInfo), INTENT(IN) :: Grid
     TYPE(PotentialInfo), INTENT(IN) :: Pot

     INTEGER :: i,j,k,n,nfit
     REAL(8),parameter ::  smallr=0.001d0
     INTEGER, parameter :: order=4
     REAL(8) :: A(4,4),X(4),B(4),xx,xxx


     n=Grid%n

!  check for possible ionic charge
      qq=-Pot%rv(n)/2
      if(qq<0.001d0) qq=0


     oneplusvtau=0.d0;   dvtaudr=0.d0

!!!! Hopefully the following code is no longer needed     
!     nfit=0
!     do i=1,n
!         if (Grid%r(i).le.smallr) then
!             nfit=i
!         else
!             exit
!         endif     
!     enddo    
!     !x write(std_out,*) 'In Set_Pot   nfit ', nfit
!     if (nfit.le.5) then
!        write(std_out,*) 'Program stopping because nfit too small '
!        stop
!     endif     
!
!     A=0.d0;B=0.d0;X=0.d0;
!     do i=1,nfit   
!        do j=1,order
!          if (j==1) then
!            xx=1
!          else
!            xx=(Grid%r(i)/smallr)**(j-1)
!          endif  
!          B(j)=B(j)+Pot%vtau(i)*xx
!          do k=1,j
!             if (k==1) then
!               xxx=xx
!             else
!               xxx=xx*(Grid%r(i)/smallr)**(k-1)
!             endif  
!             A(j,k)=A(j,k)+xxx
!             if (k<j) A(k,j)=A(k,j)+xxx
!          enddo   
!        enddo
!     enddo  
!
!     call SolveAXeqB(order,A,B)
!
!     !x write(std_out,*) 'completed SolveAXeqB ',B(1:order)
!     k=max(10,nfit/2)
!     !x write(std_out,*) 'Resetting vtau for first ',k,'  points'
!
!     xx=1.d0/smallr
!     do j=2,order
!        B(j)=B(j)*xx
!        xx=xx/smallr
!     enddo   
!     !x write(std_out,*) 'mod B: ',B(1:order)
!     do i=1,nfit
!        xx=B(1);xxx=0.d0
!        do j=2,order
!           xx=xx+B(j)*(Grid%r(i))**(j-1)
!           xxx=xxx+(j-1)*B(j)*(Grid%r(i))**(j-2)
!        enddo
!        !x write(std_out,'(i10,1p,5e15.7)')i,Grid%r(i),pot%vtau(i),xx,(pot%vtau(i)-xx),xxx
!        if(i.le.k) then 
!           pot%vtau(i)=xx
!           dvtaudr(i)=xxx
!        endif   
!     enddo   


     oneplusvtau=1.d0+Pot%vtau
     call derivative(Grid,Pot%vtau,dvtaudr)
     zxc=Pot%rvx(1)

!     write(std_out,*) 'In set_pot ', oneplusvtau(1),dvtaudr(1),zxc
!     !x do i=1,k+8
!        !x write(std_out,'(i10,1p,3e15.7)') i,Grid%r(i),oneplusvtau(i),dvtaudr(i)
!     !x enddo

   END Subroutine Set_Pot

!*******************************************************************
! SUBROUTINE wfnkedinit(Grid,l,nz,v0,energy,wfn,lwfn,istart)
!   returns the solution of the modified KS equations near r=0
!   using power series expansion assuming including vtau contributions
!   wfn=P(r)   lwfn=(1+vtau)*dP/dr
!   P(r)~~(r**(l+1))*(1+c1*r)
!   Assumes v(r) ~~ -2*nz/r+zxc/r   for r-->0
!   Assumes vtau(r) -- t0 +t1*r  for r-->0
!*******************************************************************
  SUBROUTINE wfnkedinit(Grid,l,nz,v0,energy,wfn,lwfn,istart)
   Type(GridInfo), INTENT(IN) :: Grid
   INTEGER, INTENT(IN) :: l,nz
   REAL(8), INTENT(IN) :: v0,energy
    REAL(8),INTENT(INOUT) :: wfn(:),lwfn(:)
    INTEGER, INTENT(OUT) :: istart

    REAL(8) :: rr,c1,c2,t1,t0
    INTEGER :: i,j,n

    t0=oneplusvtau(1);t1=dvtaudr(1)
    wfn=0; lwfn=0
    c1=-(2*nz-zxc+l*t1)/(2*(l+1)*(oneplusvtau(1)))

    !x write(std_out,*)' wfnkedinit for l ',l,energy;call flush_unit(std_out)
    istart=6
    do i=1,istart
       rr=Grid%r(i)
       wfn(i)=1+rr*c1
       lwfn(i)=(rr**l)*((l+1)*wfn(i)+rr*(c1))
       lwfn(i)=(t0+t1*rr)*lwfn(i)
       wfn(i)=wfn(i)*(rr**(l+1))
       !x write(std_out,'(i10,1p,10e15.7)') i,rr,wfn(i)/rr,lwfn(i),t0+t1*rr-oneplusvtau(i)
    enddo

  End SUBROUTINE wfnkedinit


!*******************************************************************
! SUBROUTINE wfnkedinitinhomo(Grid,l,proj,wfn,lwfn,istart)
!   returns the solution of the modified KS equations near r=0
!   using power series expansion including vtau contributions
!   wfn=P(r)   lwfn=(1+vtau)*dP/dr
!   P(r)~~(r**(l+1))*(c2*r**2)
!   Assumes v(r) ~~ v0   for r-->0
!   Assumes vtau(r) -- t0 +t1*r  for r-->0
!*******************************************************************
  SUBROUTINE wfnkedinitinhomo(Grid,l,proj,wfn,lwfn,istart)
   Type(GridInfo), INTENT(IN) :: Grid
   INTEGER, INTENT(IN) :: l
   REAL(8), INTENT(IN) :: proj(:)
    REAL(8),INTENT(INOUT) :: wfn(:),lwfn(:)
    INTEGER, INTENT(OUT) :: istart

    REAL(8) :: rr,c1,c2,p0,t0,t1
    INTEGER :: i,j,n

    !!!write(std_out,*) 'in wfnkedinitinhomo ', l; call flush_unit(std_out)
    !!!write(std_out,*) 'May not be correct '
    wfn=0; lwfn=0
    t0=oneplusvtau(1);t1=dvtaudr(1)
   !!  assume proj(r) -- r**(l=1)*p0  for r-->0;  determine p0
    lwfn(2:8)=proj(2:8)/Grid%r(2:8)
    call extrapolate(Grid,lwfn)
    p0=lwfn(1)
    lwfn=0
    !write(std_out,*) 'p0   ',p0; call flush_unit(std_out)
    c2=-p0/(t0*(4*l+6))
    istart=6
    do i=1,istart
       rr=Grid%r(i)
       wfn(i)=c2*(rr**(l+3))
       lwfn(i)=(l+3)*c2*(rr**(l+2))
       lwfn(i)=(t0+t1*rr)*lwfn(i)
       !x write(std_out,*) i,rr,wfn(i),lwfn(i); call flush_unit(std_out)
    enddo

    !x write(std_out,*) 'finished init'; call flush_unit(std_out)
  End SUBROUTINE wfnkedinitinhomo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form taken from wfnend in numerov_mod.F90
!  Find asymptotic form of wavefunction
!    assuming equations has form
!   (- d^2  +l(l+1) -2q   +b^2 )
!   (  ---   ------ ---        )  P(r) = 0
!   (  dr^2    r^2    r        )
!       where b^2=-energy
!
!        P(r) =
!        exp(-b*(r-rN))*r^(q/b)(1+(l*(l+1)-q/b)*(q/b-1)/(2*b)/r
!        + ...)
!      only first term is kept
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine wfnkedasym(Grid,wfn,lwfn,energy,iend)
  ! returns the solution of the modified KS equations near r=inf
  !  using only leading term
   Type(GridInfo), INTENT(IN) :: Grid
    REAL(8),INTENT(INOUT) :: wfn(:),lwfn(:)
    REAL(8), INTENT(IN) :: energy
    INTEGER, INTENT(OUT) :: iend

    REAL(8) :: rr,x,qx,qbb,b,cn,term,fac,rN
    INTEGER :: i,j,n


    if (energy>0.d0) then
       write(std_out,*) 'Error in wfnkedasym -- energy > 0', energy
       stop
    endif

    wfn=0; lwfn=0; 
    n=Grid%n

    b=oneplusvtau(n)
    if (abs(1.d0-oneplusvtau(n))> 1.d-3) then
        write(std_out,*) '*****WARNING***** Long range vtau --',&
&        oneplusvtau(n-3),oneplusvtau(n-2),oneplusvtau(n-1),oneplusvtau(n) 
    endif
    x=sqrt(-energy/b)
    qx=qq/b     !  Possible net ionic charge
    qx=(qx/x)
    iend=5
    do i=n-iend,n
       wfn(i)=ddexp(-x*(Grid%r(i)-Grid%r(n-iend)))
       if (qx>0.d0) then
               rr=(Grid%r(i)/Grid%r(n-iend))**qx
               wfn(i)=wfn(i)*rr
       endif        
       lwfn(i)=oneplusvtau(i)*wfn(i)*(qx/Grid%r(i)-x)

    enddo
  end subroutine wfnkedasym

  !**********************************************************************
  !      subroutine unboundked(Grid,Pot,nr,l,energy,wfn,nodes)
  !  pgm to solve radial kinetic energy derivative equations for unbound states
  !    at energy 'energy' and at angular momentum l
  !
  !    with potential rv/r, given in uniform linear or log mesh of n points
  !   assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
  !
  !  nz=nuclear chargedd
  !
  !  Does not use Noumerov algorithm -- but uses coupled first-order
  !       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
  !
  ! also returns node == number of nodes for calculated state
  !************************************************************************
  SUBROUTINE unboundked(Grid,Pot,nr,l,energy,wfn,nodes)
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
       write(std_out,*) 'Error in unboundked -- nr > n', nr,n
       STOP
    ENDIF

    call Set_Pot(Grid,Pot)

    allocate(lwfn(nr),zz(2,2,nr),yy(2,nr),stat=ierr)
       if (ierr/=0) then
          write(std_out,*) ' allocation error in unboundked ', nr,ierr
          stop
       endif

    lwfn=0;zz=0;yy=0;

    call wfnkedinit(Grid,l,Pot%nz,Pot%v0,energy,wfn,lwfn,istart)
    call setupforcfdsol(Grid,Pot%rv,1,istart,nr,l,energy,wfn,lwfn,yy,zz)
    call cfdsoliter(Grid,zz,yy,istart,nr)
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

  END SUBROUTINE unboundked

  
  !**********************************************************************
  !      subroutine unboundkedinhomo(Grid,Pot,proj,nr,l,energy,wfn)
  !  pgm to solve radial kinetic energy derivative equations for unbound states
  !    at energy 'energy' and at angular momentum l
  !
  !    with potential rv/r, given in uniform linear or log mesh of n points
  !   assuming proj(r)=C*r**(l+1)*polynomial(r) for r==0;
  !
  !  nz=nuclear chargedd
  !
  !  Does not use Noumerov algorithm -- but uses coupled first-order
  !       equations from David Vanderbilt, Marc Torrent, and Francois Jollet
  !
  ! also returns node == number of nodes for calculated state
  !************************************************************************
  SUBROUTINE unboundkedinhomo(Grid,Pot,proj,nr,l,energy,wfn)
    TYPE(GridInfo), INTENT(IN)  :: Grid
    TYPE(PotentialInfo), INTENT(IN)  :: Pot
    REAL(8), INTENT(IN) :: proj(:)
    INTEGER, INTENT(IN) :: nr,l
    REAL(8), INTENT(IN) :: energy
    REAL(8), INTENT(INOUT) :: wfn(:)

    INTEGER :: n,i,j,k,ierr,istart
    REAL(8) :: scale
    REAL(8), allocatable :: lwfn(:),zz(:,:,:),yy(:,:),ff(:,:)

    !write(std_out,*) 'in unboundkedinhomo ' ; call flush_unit(std_out)

    n=Grid%n
    IF (nr > n) THEN
       write(std_out,*) 'Error in unboundked -- nr > n', nr,n
       STOP
    ENDIF

    call Set_Pot(Grid,Pot)

    allocate(lwfn(nr),zz(2,2,nr),yy(2,nr),ff(2,nr),stat=ierr)
       if (ierr/=0) then
          write(std_out,*) ' allocation error in unboundked ', nr,ierr
          stop
       endif

       !x write(std_out,*) ' after Set_Pot '; call flush_unit(std_out)
    lwfn=0;zz=0;yy=0;ff=0

 
    call wfnkedinitinhomo(Grid,l,proj,wfn,lwfn,istart)
    call setupforcfdsol(Grid,Pot%rv,1,istart,nr,l,energy,wfn,lwfn,yy,zz)
    ff(2,:)=-proj(:)
    !!!write(std_out,*) 'before inhomcfdsol '; call flush_unit(std_out)
    call inhomocfdsol(Grid,zz,yy,ff,istart,nr)
    !x write(std_out,*) 'after inhomcfdsol '; call flush_unit(std_out)
    call getwfnfromcfdsol(1,nr,yy,wfn)
    !
    !

    !!!write(std_out,*) 'completed unboundkedhomo'; call flush_unit(std_out)
    deallocate(lwfn,yy,zz,ff)
    !!!write(std_out,*) 'completed unboundkedhomo'; call flush_unit(std_out)

  END SUBROUTINE unboundkedinhomo

!******************************************************************
!  SUBROUTINE Boundked(Grid,Pot,eig,wfn,l,nroot,emin,ierr,success)
!******************************************************************
  SUBROUTINE Boundked(Grid,Pot,eig,wfn,l,nroot,emin,ierr,success)
    !  pgm to solve  meta gga radial equation for nroot bound state
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
       write(std_out,*) ' Error in boundked allocation ',i,n
       STOP
    ENDIF

    success=.true.
    allocate(lwfn(n),zz(2,2,n),yy(2,n),stat=i)
       if (i/=0) then
          write(std_out,*) ' allocation error in boundked ', n,i
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

    write(std_out,*) 'z , l = ',nz,l; call flush_unit(std_out)
    ! check how many roots expected by integration outward at
    !   energy = 0
    energy = 0
    call Set_Pot(Grid,Pot)
    p1=0;lwfn=0;zz=0;yy=0;
    call wfnkedinit(Grid,l,Pot%nz,Pot%v0,energy,p1,lwfn,istart)
    !x write(std_out,*) ' after wfnkedinit ', istart
    !
    !start outward integration
    call setupforcfdsol(Grid,Pot%rv,1,istart,n,l,energy,p1,lwfn,yy,zz)
    write(std_out,*) ' after setupfocfdsol ', istart,energy
    call cfdsoliter(Grid,zz,yy,istart,n)
    !x do i=1,istart+5
    !x    write(std_out,'(i10,1p20e15.7)') i,Grid%r(i),yy(1,i),zz(2,1,i)
    !x enddo

    call getwfnfromcfdsol(1,n,yy,p1)
    !x write(std_out,*) 'afterwfnfromcfdsol';call flush_unit(std_out)
    !open(1001,file='initwfn',form='formatted')
    ! do i=1,n
    !     write(1001,*) i,yy(1,i),p1(i)
    !     enddo
    ! close(1001)
    node=countnodes(2,n,p1)

    write(std_out,*) ' nodes at e=0  ', node

    mxroot=node+1
    ntroot=node
    IF (mxroot.LT.nroot) THEN
       write(std_out,*)'error in boundked - for l = ',l
       write(std_out,*) nroot,' states requested but only',mxroot,' possible'
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
       !!!energy=emin+err
       energy=eig(iroot)     !  assume good initial guess
       IF (energy.LT.emin) energy=emin
       IF (energy.GT.emax) energy=emax
       ok=.FALSE.
       !write(std_out,*) 'iter,iroot,energy',iter,iroot,energy
       !write(std_out,*) 'emin,max',emin,emax
       BigIter: DO iter=1,niter
          !write(std_out,*) 'In iter with energy', iter,energy,niter,l,iroot
          !  start inward integration
          !  start integration at n
          ! find classical turning point
          !call ClassicalTurningked(Grid,Pot%rv,l,energy,match)
          call ClassicalTurningpoint(Grid,Pot%rv,l,energy,match)
          match=max(match,10); match=min(match,n-20)
          call wfnkedasym(Grid,p2,lwfn,energy,iend)
          call setupforcfdsol(Grid,Pot%rv,n-iend,n,n,l,energy,p2,lwfn,yy,zz)
          call cfdsoliter(Grid,zz,yy,n-iend,match)
          call getwfnfromcfdsol(match,n,yy,p2)
          match=match+6
          rin=Gfirstderiv(Grid,match,p2)/p2(match)

          call wfnkedinit(Grid,l,Pot%nz,Pot%v0,energy,p1,lwfn,istart)
          call setupforcfdsol(Grid,Pot%rv,1,istart,n,l,energy,p1,lwfn,yy,zz)
          call cfdsoliter(Grid,zz,yy,istart,match+6)
          call getwfnfromcfdsol(1,match+6,yy,p1)
          node= countnodes(2,match+6,p1)

          rout=Gfirstderiv(Grid,match,p1)/p1(match)
          ! check whether node = (iroot-1)
          !   not enough nodes -- raise energy
          IF (node.LT.iroot-1) THEN
             emin=MAX(emin,energy)-err
             energy=emax-(emax-energy)*ranx()
             ifac=9
             !   too many nodes -- lower energy
          ELSEIF (node.GT.iroot-1) THEN
             IF (energy.LT.emin) THEN
                ierr=ierr+9*(10**(iroot-1))
                write(std_out,*) 'boundked error -- emin too high',node,iroot-1,l,nz,emin,energy
                do i=2,n
                   write(999,'(1p,4e15.7)') Grid%r(i),Pot%rv(i),p1(i)
                enddo
                STOP
             ENDIF
             emax=MIN(emax,energy+err)
             energy=emin+(energy-emin)*ranx()
             !   correct number of nodes -- estimate correction
          ELSEIF (node.EQ.iroot-1) THEN
             DO j=1,match
                p1(j)=p1(j)/p1(match)
                         !write(std_out,*) 'j,p1',j,p1(j)
             ENDDO
             DO j=match,n
                p1(j)=p2(j)/p2(match)
                         !write(std_out,*) 'j,p2',j,p1(j)
             ENDDO
             scale=1.d0/overlap(Grid,p1,p1)
             dele=(rout-rin)*scale
                  write(std_out,*) 'energy,dele,scale',energy,dele,scale
             x=ABS(dele)
             IF (x.LT.best) THEN
                scale=1.d0/overlap(Grid,p1,p1)
                scale=SQRT(scale)
                p1(1:n)=p1(1:n)*scale
                call filter(n,p1,machine_zero)
                wfn(1:n,iroot)=p1(1:n)
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
          write(std_out,*) 'no convergence in boundked',iroot,l,dele,energy
          write(std_out,*) ' best guess of eig, dele = ',eig(iroot),best
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
             write(std_out,*) 'returning from boundked -- ierr=',ierr
  END SUBROUTINE Boundked

  SUBROUTINE ClassicalTurningked(Grid,rv,l,energy,turningpoint)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: rv(:)
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: energy
    INTEGER, INTENT(OUT) :: turningpoint

    INTEGER :: i,n
    REAL(8), ALLOCATABLE :: v(:)

    n=Grid%n
    ALLOCATE(v(n), stat=i)
    IF (i /= 0) THEN
       write(std_out,*) 'Allocation error in ClassicalTurningked ', i,n
       STOP
    ENDIF

    v=0
    v(2:n)=rv(2:n)/Grid%r(2:n)+&
&                l*(l+1)*(oneplusvtau(2:n))/(Grid%r(2:n)**2)

    turningpoint=n
    DO i=n,2,-1
       IF (v(i)<energy) EXIT
    ENDDO
    turningpoint=i
    turningpoint=MIN(turningpoint,FindGridIndex(Grid,10.0d0))

    write(std_out,*) 'Found turning point at ', turningpoint, Grid%r(turningpoint)

    DEALLOCATE(v)

  END SUBROUTINE ClassicalTurningked
  subroutine setupforcfdsol(Grid,rv,i1,i2,n,l,energy,wfn,lwfn,yy,zz)
     Type(gridinfo), INTENT(IN) :: Grid
     INTEGER, INTENT(IN) :: i1,i2,n,l
     REAL(8), INTENT(IN) :: energy
     REAL(8), INTENT(IN) :: wfn(:),lwfn(:),rv(:)
     REAL(8), INTENT(INOUT) :: yy(:,:),zz(:,:,:)

     INTEGER :: i
     REAL(8) :: x

      x=l*(l+1)
      yy=0;zz=0
      yy(1,i1:i2)=wfn(i1:i2)
      yy(2,i1:i2)=lwfn(i1:i2)

      do  i=1,n
       zz(1,2,i)=1.d0/oneplusvtau(i)
       if(i==1) then
         zz(2,1,i)=0.d0
       else  
         zz(2,1,i)=oneplusvtau(i)*x/(Grid%r(i)*Grid%r(i))+&
&            dvtaudr(i)/Grid%r(i)+(rv(i)/Grid%r(i)-energy)
       endif  
      enddo

   end subroutine setupforcfdsol


END MODULE radialked
