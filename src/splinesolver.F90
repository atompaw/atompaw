!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    This module contains the following active subroutines:
!      initsplinesolver,deallocatesplinesolver,initpotforsplinesolver
!        Boundsplinesolver
! 07/11/2021 -- NAWH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE splinesolver

  USE io_tools
  USE GlobalMath
  USE gridmod
  USE interpolation_mod
  USE search_sort
  USE atomdata
  USE excor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  For solving eigenvalue problem, need local logarithmic
  !    grid -- grids  which has tuned local r0 and h
  !    The eigensolver does not use r=0 point
  !     Has dimensions ns=n-1
  !    Local arrays srv,svtau,sden,stau,etc refer to grids
  !    Output interpolated to universal grid of calculation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
     INTEGER, private :: ns      ! # spline nodes, not including origin
     REAL(8), private :: r0,h
     TYPE(gridinfo), private :: Grids     ! private grid but includes origin
     REAL(8), private, allocatable :: u(:),pref(:),rr1(:),rr2(:)
     REAL(8), private, allocatable :: srv(:),svtau(:),soneplusvt(:),sdvt(:)

     !special extra private variable for MGGA case (needvtau==.true.)
     TYPE(gridinfo), private ::Gridf   ! private fine linear grid
     REAL(8), private, allocatable :: fvtau(:),fdvtaudr(:),frvx(:)
     REAL(8), private, allocatable :: fden(:),ftau(:)

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Program uses two grids -- "universal grid" Grid
!   and modified grid Grids, with fixed Grids%n=401 and ns=400
!     and Grids%r0=0.1 which are found to work well for splinesolver
!   Internally need to omit origin and so ns=Grids%n-1
!   For MGGA case (needvtau=.true.) also need fine linear grid
!      Gridf
!  Setup local private grid  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE initsplinesolver(Grid)
      Type(GridInfo), INTENT(IN) :: Grid       !Input universal grid

      INTEGER :: i,nf
      REAL(8) :: hf,x

      r0=0.1d0    ! 0.1 seems to be a good idea
      ns=400      ! seems to be good enough
      h=log(1.d0+Grid%r(Grid%n)/r0)/ns
      call initgridwithn(Grids,2,ns+1,r0,h)  !local loggrid

      write(std_out,*) 'initsplinesolve ', r0,ns,h;call flush_unit(std_out)

      allocate(u(ns+1),pref(ns+1),rr1(ns+1),rr2(ns+1))
      allocate(srv(ns+1),svtau(ns+1),sdvt(ns+1))
      allocate(soneplusvt(ns+1))
      do i=1,ns+1
         u(i)=(i-1)*h
         !r(i)=r0*(exp(u(i))-1.d0)      !grids%r
         pref(i)=exp(0.5d0*u(i))
         rr1(i)=((Grids%r(i)+r0))
         rr2(i)=((Grids%r(i)+r0)**2)
      enddo    

      if(needvtau) then    ! set up fine linear grid
          nf=20001
          hf=Grid%r(Grid%n)/(nf-1)    
          call initgridwithn(Gridf,1,nf,0.d0,hf)  !local fine linear grid
          allocate(fvtau(nf),fdvtaudr(nf),frvx(nf),fden(nf),ftau(nf))
          fvtau=0.d0;fdvtaudr=0.d0;frvx=0.d0;fden=0.d0;ftau=0.d0
      endif    

      soneplusvt=1.d0;
      svtau=0.d0;sdvt=0.d0
   END SUBROUTINE initsplinesolver

   SUBROUTINE deallocatesplinesolver
      deallocate(u,pref,rr1,rr2)
      deallocate(srv,svtau,soneplusvt)
      call destroygrid(Grids)

      if(needvtau) then
          deallocate(fvtau,fdvtaudr,frvx,fden,ftau)
      endif
   END SUBROUTINE deallocatesplinesolver   

   SUBROUTINE initpotforsplinesolver(Grid,Pot,den,tau)   
      Type(GridInfo), INTENT(IN) :: Grid
      Type(Potentialinfo), INTENT(INOUT) :: Pot   !Universal grid
      REAL(8), INTENT(IN) :: den(:),tau(:)  !Universal grid

      INTEGER :: i,j,k,l,m,n
      INTEGER :: icount=0
      REAL(8), allocatable :: dum(:),dum1(:)
      REAL(8) :: etxc,eexc

      n=Grid%n
      
      write(std_out,*) 'initpot ', n,ns;call flush_unit(std_out)

      If (needvtau) then
          allocate(dum(ns+1),dum1(ns+1))
          call interpfunc(n,Grid%r,den,Gridf%n,Gridf%r,fden)     
          call interpfunc(n,Grid%r,tau,Gridf%n,Gridf%r,ftau)     
          call exch(Gridf,fden,frvx,etxc,eexc,tau=ftau,vtau=fvtau)
         write(std_out,*) 'called exch from splinesolver ', eexc
          call nderiv(Gridf%h,fvtau,fdvtaudr,Gridf%n,i)
        ! write(7000+icount,*)'#r fden  ftau   fvxc   fvtau     fdvtaudr'
        ! do i=1,Gridf%n
        !     write(7000+icount,'(1p,20e15.7)')Gridf%r(i),&
!&              fden(i),ftau(i),frvx(i),fvtau(i),fdvtaudr(i)   
        ! enddo    
        ! close(7000+icount)
          !rv
          dum=Pot%rvn+Pot%rvh    !  presumably these are smooth
          call interpfunc(n,Grid%r,dum,ns+1,Grids%r,srv)
          call interpfunc(Gridf%n,Gridf%r,frvx,Grids%n,Grids%r,dum1)
          srv=srv+dum1
          call interpfunc(Gridf%n,Gridf%r,fvtau,Grids%n,Grids%r,svtau)
          call interpfunc(Gridf%n,Gridf%r,fdvtaudr,Grids%n,Grids%r,sdvt)
          soneplusvt=1.d0+svtau
         ! write(9000+icount,*)'#r srv  svxc   svtau     sdvt    sdvtdu'
         !  do i=1,Grids%n
         !     write(9000+icount,'(1p,20e15.7)')Grids%r(i),&
! &               srv(i),dum1(i),svtau(i),sdvt(i),sdvt(i)*Grids%drdu(i)
         ! enddo    
         ! close(9000+icount)
         ! icount=icount+1
          sdvt=sdvt*Grids%drdu      ! needed in algorithm
      
          deallocate(dum,dum1)
      else    !  finegrid not needed        
          call interpfunc(n,Grid%r,Pot%rv,ns+1,Grids%r,srv)  
          soneplusvt=1.d0;
          svtau=0.d0;sdvt=0.d0
      endif

  END   SUBROUTINE initpotforsplinesolver

  SUBROUTINE Boundsplinesolver(Grid,l,neig,eig,wfn,otau,OK)
      Type(GridInfo), INTENT(IN) :: Grid
      INTEGER, INTENT(IN) :: l,neig
      REAL(8), INTENT(INOUT) :: eig(:),wfn(:,:),otau(:,:)
      LOGICAL, INTENT(OUT) :: OK
          
      real(8), allocatable :: A(:,:),B(:,:),G(:,:),D(:),DL(:),DU(:),work(:)
      real(8), allocatable :: vl(:,:),vr(:,:),wr(:),wi(:)
      integer, allocatable :: lut(:)
      REAL(8), allocatable :: P(:),MP(:)
      REAL(8), allocatable :: dum(:),dP(:)
      REAL(8), allocatable :: F1(:,:),F2(:,:),S(:),E(:),V(:)
      integer :: i,j,m,n, info,lwork,nu
      real(8) :: ol,x
      integer :: icount=0

      write(std_out,*) 'entering boundspline with l,neig ', l, neig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Note that although the universal Grid and the local Grids
!     have the same range, they differ by the number of points and
!     the r0 parameter.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     write(std_out,*) 'in splinesolver ';call flush_unit(std_out)

      OK=.false.
      n=ns+1       !Grids%n   within this routine only
      nu=Grid%n     ! Universal grid
      ol=l*(l+1)
      lwork=ns**2
      allocate(A(ns,ns),B(ns,ns),G(ns,ns),D(ns),DL(ns),DU(ns),work(lwork))
      allocate(vl(ns,ns),vr(ns,ns),wr(ns),wi(ns),lut(ns))
      allocate(P(n),MP(n),dum(nu),dP(nu))
      allocate(F1(ns,ns),F2(ns,ns),S(n),E(n),V(n))

      D(1:ns)=2.d0
      DL(1:ns-1)=0.5d0
      DU(1:ns-1)=0.5d0
   !correct first row assuming wfn=r^(l+1)*(W0+r*W1)
      D(1)=0.5d0*(5.d0+2*l)/(1.d0+l)
      B=0.d0
      do i=1,ns
         B(i,i)=1.d0
      enddo   

      call dgtsv(ns,ns,DL,D,DU,B,ns,info)
      write(std_out,*) 'Completed dgtsv with info = ', info;call flush_unit(std_out)

      A =0.d0
      do i=1,ns
      A(i,i)=-2.d0
      enddo
      do i=1,ns-1
      A(i,i+1)=1.d0
      A(i+1,i)=1.d0
      enddo   

      ! correct first row values
      A(1,1)=-0.5d0*(l+4.d0);  
      A=3*A/(h**2)

      G=0.d0
      G=MATMUL(B,A)      !  G stores transformation to find M=G*y

      ! Calculate full matrix
      A=0;S=0;E=0;V=0
      ! These arrays  have the full range  1..n
      S=-soneplusvt/rr2
      E=-sdvt/rr2
      V=0.5d0*E-0.25d0*S     !Including only non diverging term
      do i=2,n
        V(i)=V(i)+soneplusvt(i)*ol/(Grids%r(i)**2)+sdvt(i)/(Grids%r(i)*rr1(i)) &
&               +srv(i)/Grids%r(i)
      enddo
      F1=0.d0;F2=0.d0
      do i=1,ns
         F1(i,i)=S(i+1)-E(i+1)*h/3
         F1(i,i+1)=-E(i+1)*h/6       
         F2(i,i)=V(i+1)-E(i+1)/h
         F2(i,i+1)=E(i+1)/h
      ENDDO   

      A=MATMUL(F1,G)+F2

      call dgeev('N','V',ns,A,ns,wr,wi,vl,ns,vr,ns,work,lwork,info)
      write(std_out,*) 'dgeev completed with info = ',info;call flush_unit(std_out)

      if (info/=0) then
          OK=.false.
          return
       endif    
      call insertion_sort(wr,lut,.true.)      

      write(std_out,*) 'Results  for l = ', l, (wr(lut(i)),i=1,10);call  flush_unit(std_out)

      write(std_out,*)  'enumeration for neig solutions ', neig;call flush_unit(std_out)
      do m=1,neig
          write(std_out,'(1p,2e20.8)') wr(lut(m)),wi(lut(m));call flush_unit(std_out)


         !!!! Perhaps not needed?? 
         ! If (wr(lut(m))>0.d0.or.(abs(wi(lut(m)))>1.d-8)) then
         !     write(std_out,*) 'Problem in splinesolver program stopping '
         !     write(std_out,*) 'eig - ', wr(lut(m)),wi(lut(m))
         !     stop
         ! endif    
          eig(m)=wr(lut(m))

          D=vr(:,lut(m))     !Q
          DL=MATMUL(G,D)     !MQ

          ! now extend grid to r=0   Still Q and MQ
          P=0;MP=0
          P(2:ns+1)=D(1:ns)
          MP(2:ns+1)=DL(1:ns)
          if(l==0) then
            dum=0      
            do i=2,10      
               dum(i)=P(i)/Grids%r(i)
            enddo   
            call extrapolate(Grids,dum) 
            x=1.d0/(S(1)-E(1)*h/3)
            MP(1)=(E(1)*(MP(2)*h/6-P(2))-(sdvt(1)/rr1(1)+srv(1))*dum(1))*x
             write(std_out,*) 'MP(1) for l=0',MP(1),dum(1)
          endif
          if(l==1) then
            dum=0      
            do i=2,10      
               dum(i)=P(i)/(Grids%r(i)**2)
            enddo   
            call extrapolate(Grids,dum) 
            x=1.d0/(S(1)-E(1)*h/3)
            MP(1)=(E(1)*(MP(2)*h/6-P(2))-(2*soneplusvt(1))*dum(1))*x
             write(std_out,*) 'MP(1) for l=1',MP(1),dum(1)
          endif


          MP=pref*(MP-0.25d0*P)/rr2
          P=pref*P


         ! write(std_out,*) ' Writing out some initial P and MP '
         ! do i=1,15
         !     write(std_out,'(1p,3e15.7)') Grids%r(i),P(i),MP(i)
         ! enddo    
          call specialinterp(l,n,Grids%r,P,MP,Grid%n,Grid%r,wfn(:,m),dP(:))
          dum=0.d0; dum(2:nu)=wfn(2:nu,m)/Grid%r(2:nu)
          call extrapolate(Grid,dum)
          do i=1,nu
             otau(i,m)= (dP(i)-dum(i))**2+l*(l+1)*(dum(i))**2
          enddo   

          x=overlap(Grid,wfn(:,m),wfn(:,m))
          write(std_out,*) 'overlap integral ', x
          x=1.d0/x
          otau(:,m)=x*otau(:,m)
          x=sqrt(x)
          wfn(:,m)=x*wfn(:,m)    ! should be normalized now
         ! call taufromwfn(Grid,wfn(:,m),l,dum)    ! testing

         ! do i=1,Grid%n
         !    write(500+icount,'(1p,4e16.8)') Grid%r(i),wfn(i,m),otau(i,m),dum(i)
         ! enddo   

         ! close(500+icount)
         ! icount=icount+1
      enddo    
      OK=.true.
      deallocate(A,B,G,D,DL,DU,work)
      deallocate(vl,vr,wr,wi,lut,P,MP,dP,dum)
      deallocate(F1,F2,S,E,V)

  end SUBROUTINE Boundsplinesolver
END MODULE splinesolver
