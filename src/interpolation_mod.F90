!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Module for interpolation of functions between different grids
!    Based on DeBoor's cubic splines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!     interpfunc, cubspl (from de Boor) ,specialinterp,checkrange
!     binterpfunc -- a special version of interpfunc that assumes
!     cubspl is called first -- used for setting boundary values
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

Module interpolation_mod

   Use io_tools

   Implicit none

   CONTAINS

   SUBROUTINE interpfunc(nin,rin,fin,nout,rout,fout)
     INTEGER, INTENT(IN) :: nin,nout
     REAL(8), INTENT(IN) :: rin(:),fin(:),rout(:)
     REAL(8), INTENT(INOUT) :: fout(:)

     REAL(8), ALLOCATABLE :: c(:,:)
     INTEGER :: i,j,nacount=0
     REAL(8) :: x
     LOGICAL :: leftin,lefton,leftout,rightin,righton,rightout

     !  check if grid is within interpolation range
     call checkrange(rin(1),rout(1),leftin,lefton,leftout)
     call checkrange(rout(nin),rin(nin),rightin,righton,rightout)
     if (leftout.or.rightout) then
      write(std_out,*) 'Grid error in interpfunc ',rin(1),rout(1),rin(nin),rout(nout)
      stop
     endif

     allocate(c(4,nin))

     c=0;
     c(1,:)=fin(:)
     call cubspl(rin,c,nin,0,0)

     fout=0

     do i=1,nout
        do j=1,nin-1
          call checkrange(rin(j),rout(i),leftin,lefton,leftout)
          call checkrange(rout(i),rin(j+1),rightin,righton,rightout)
           if ((leftin.or.lefton).and.(rightin.or.righton)) then
              x=rout(i)-rin(j)
              fout(i)=c(1,j)+x*(c(2,j)+x*(c(3,j)+x*c(4,j)/3)/2)
              exit
           endif
        enddo
     enddo

     !nacount=nacount+1

     !do i=1,nout
     !   write(300+nacount,'(1p,5e15.7)') rout(i),fout(i)
     !enddo

     !do i=1,nin
     !   write(500+nacount,'(1p,5e15.7)') rin(i),fin(i)
     !enddo

     !close(300+nacount)
     !close(500+nacount)

     deallocate(c)
   END SUBROUTINE interpfunc

   SUBROUTINE binterpfunc(c,nin,rin,fin,nout,rout,fout)
     REAL(8), INTENT(IN) :: c(:,:)   !input from cubicspl
     INTEGER, INTENT(IN) :: nin,nout
     REAL(8), INTENT(IN) :: rin(:),fin(:),rout(:)
     REAL(8), INTENT(INOUT) :: fout(:)

     INTEGER :: i,j
     REAL(8) :: x
     LOGICAL :: leftin,lefton,leftout,rightin,righton,rightout

     !  check if grid is within interpolation range
     call checkrange(rin(1),rout(1),leftin,lefton,leftout)
     call checkrange(rout(nin),rin(nin),rightin,righton,rightout)
     if (leftout.or.rightout) then
      write(std_out,*) 'Grid error in binterpfunc ',rin(1),rout(1),rin(nin),rout(nout)
      stop
     endif

     fout=0

     do i=1,nout
        do j=1,nin-1
          call checkrange(rin(j),rout(i),leftin,lefton,leftout)
          call checkrange(rout(i),rin(j+1),rightin,righton,rightout)
           if ((leftin.or.lefton).and.(rightin.or.righton)) then
              x=rout(i)-rin(j)
              fout(i)=c(1,j)+x*(c(2,j)+x*(c(3,j)+x*c(4,j)/3)/2)
              exit
           endif
        enddo
     enddo

   END SUBROUTINE binterpfunc

!! from webpage
!!   http://pages.cs.wisc.edu/~deboor/pgs/cubspl.f
!!   Transformed into fortran 90 and REAL(8)
      subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
!  from  * a practical guide to splines *  by c. de boor
!     ************************  input  ***************************
!     n = number of data points. assumed to be .ge. 2.
!     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
!        data points. tau is assumed to be strictly increasing.
!     ibcbeg, ibcend = boundary condition indicators, and
!     c(2,1), c(2,n) = boundary condition information. specifically,
!        ibcbeg = 0  means no boundary condition at tau(1) is given.
!           in this case, the not-a-knot condition is used, i.e. the
!           jump in the third derivative across tau(2) is forced to
!           zero, thus the first and the second cubic polynomial pieces
!           are made to coincide.)
!        ibcbeg = 1  means that the slope at tau(1) is made to equal
!           c(2,1), supplied by input.
!        ibcbeg = 2  means that the second derivative at tau(1) is
!           made to equal c(2,1), supplied by input.
!        ibcend = 0, 1, or 2 has analogous meaning concerning the
!           boundary condition at tau(n), with the additional infor-
!           mation taken from c(2,n).
!     ***********************  output  **************************
!     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
!        of the cubic interpolating spline with interior knots (or
!        joints) tau(2), ..., tau(n-1). precisely, in the interval
!        (tau(i), tau(i+1)), the spline f is given by
!           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
!        where h = x - tau(i). the function program *ppvalu* may be
!        used to evaluate f or its derivatives from tau,c, l = n-1,
!        and k=4.
      integer ibcbeg,ibcend,n,   i,j,l,m
      real(8) c(4,n),tau(n),   divdf1,divdf3,dtau,g
!****** a tridiagonal linear system for the unknown slopes s(i) of
!  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
!  ination, with s(i) ending up in c(2,i), all i.
!     c(3,.) and c(4,.) are used initially for temporary storage.
      l = n-1
!compute first differences of tau sequence and store in c(3,.). also,
!compute first divided difference of data and store in c(4,.).
      do  m=2,n
         c(3,m) = tau(m) - tau(m-1)
         c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
      enddo
!construct first equation from the boundary condition, of the form
!             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
      if (ibcbeg-1 <0)                  go to 11
      if (ibcbeg-1==0)                  go to 15
      if (ibcbeg-1 >0)                  go to 16

   11 if (n .gt. 2)                     go to 12
!     no condition at left end and n = 2.
      c(4,1) = 1.d0
      c(3,1) = 1.d0
      c(2,1) = 2.d0*c(4,2)
                                        go to 25
!     not-a-knot condition at left end and n .gt. 2.
   12 c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) =((c(3,2)+2.d0*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))/c(3,1)
                                        go to 19
!     slope prescribed at left end.
   15 c(4,1) = 1.d0
      c(3,1) = 0.d0
                                        go to 18
!     second derivative prescribed at left end.
   16 c(4,1) = 2.d0
      c(3,1) = 1.d0
      c(2,1) = 3.d0*c(4,2) - c(3,2)/2.d0*c(2,1)
   18 if(n .eq. 2)                      go to 25
!  if there are interior knots, generate the corresp. equations and car-
!  ry out the forward pass of gauss elimination, after which the m-th
!  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
   19 do m=2,l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.d0*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
         c(4,m) = g*c(3,m-1) + 2.d0*(c(3,m) + c(3,m+1))
      enddo
!construct last equation from the second boundary condition, of the form
!           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
!     if slope is prescribed at right end, one can go directly to back-
!     substitution, since c array happens to be set up just right for it
!     at this point.
      if (ibcend-1 <0)                  go to 21
      if (ibcend-1==0)                  go to 30
      if (ibcend-1 >0)                  go to 24

   21 if (n .eq. 3 .and. ibcbeg .eq. 0) go to 22
!     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
!     left end point.
      g = c(3,n-1) + c(3,n)
      c(2,n) = ((c(3,n)+2.d0*g)*c(4,n)*c(3,n-1) &
&                 + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g = -g/c(4,n-1)
      c(4,n) = c(3,n-1)
                                        go to 29
!     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
!     knot at left end point).
   22 c(2,n) = 2.d0*c(4,n)
      c(4,n) = 1.d0
                                        go to 28
!     second derivative prescribed at right endpoint.
   24 c(2,n) = 3.d0*c(4,n) + c(3,n)/2.d0*c(2,n)
      c(4,n) = 2.d0
                                        go to 28
   25 continue
      if (ibcend-1 <0)                  go to 26
      if (ibcend-1==0)                  go to 30
      if (ibcend-1 >0)                  go to 24

   26 if (ibcbeg .gt. 0)                go to 22
!     not-a-knot at right endpoint and at left endpoint and n = 2.
      c(2,n) = c(4,n)
                                        go to 30
   28 g = -1.d0/c(4,n-1)
!complete forward pass of gauss elimination.
   29 c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
!carry out back substitution
   30 j = l
   40    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         if (j .gt. 0)                  go to 40
!****** generate cubic coefficients in each interval, i.e., the deriv.s
!  at its left endpoint, from value and slope at its endpoints.
      do  i=2,n
         dtau = c(3,i)
         divdf1 = (c(1,i) - c(1,i-1))/dtau
         divdf3 = c(2,i-1) + c(2,i) - 2.d0*divdf1
         c(3,i-1) = 2.d0*(divdf1 - c(2,i-1) - divdf3)/dtau
         c(4,i-1) = (divdf3/dtau)*(6.d0/dtau)
      enddo
 END subroutine cubspl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!  Subroutine to use input from Ahlberg spline interpolation
!     for node points rin(1..nin) and function values yin and
!     second derivative values Min to interpolate to radial grid
!     rout with values yout.    For the range 0 \le r \le rin(1)
!     it is assumed that yout=r**(l+1)(W0+W1*r) where l denotes
!     the angular momentum of the function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine specialinterp(l,nin,rin,yin,MMin,nout,rout,yout,ypout)
    integer, intent(IN) :: l,nin,nout
    real(8), intent(IN) :: rin(:),yin(:),MMin(:),rout(:)
    real(8), intent(INOUT) :: yout(:),ypout(:)

    real(8), allocatable :: h(:),C(:,:)
    real(8) :: W0,W1,x
    integer :: i,j,k
    logical :: ok
     LOGICAL :: leftin,lefton,leftout,rightin,righton,rightout

    allocate(h(nin),C(4,nin))

    h=0
    do i=2,nin
       h(i)=rin(i)-rin(i-1)
    enddo  
    !h(1) is hopefully not used
    C=0
    C(1,1:nin)=yin(1:nin)
    C(3,1:nin)=MMin(1:nin)
    do i=1,nin-1
       C(4,i)=(MMin(i+1)-MMin(i))/h(i+1)
       C(2,i)=((yin(i+1)-yin(i))/h(i+1)-(MMin(i+1)+2*MMin(i))*h(i+1)/6)
    enddo   

    yout=0;ypout=0
     !  check if grid is within interpolation range
     call checkrange(rin(1),rout(1),leftin,lefton,leftout)
     call checkrange(rout(nin),rin(nin),rightin,righton,rightout)
     if (leftout.or.rightout) then
      write(std_out,*) 'Grid error in specialint ',rin(1),rout(1),rin(nin),rout(nout)
      stop
     endif 

    do i=1,nout
       do j=1,nin-1
          call checkrange(rin(j),rout(i),leftin,lefton,leftout)
          call checkrange(rout(i),rin(j+1),rightin,righton,rightout)
           if ((leftin.or.lefton).and.(rightin.or.righton)) then
             x=rout(i)-rin(j)
             yout(i)=c(1,j)+x*(c(2,j)+x*(c(3,j)+x*c(4,j)/3)/2)
             ypout(i)=c(2,j)+x*(c(3,j)+0.5d0*x*c(4,j))
             exit
          endif
       enddo
    enddo

    deallocate(h,C)
  end subroutine specialinterp    

  subroutine  checkrange(rin,rout,inrange,onrange,outofbounds)
      LOGICAL, INTENT(OUT) :: inrange,onrange,outofbounds
      REAL(8), INTENT(IN) :: rin,rout    
      REAL(8), parameter :: tol=1.d-7

      inrange=.false.;onrange=.false.;outofbounds=.false.
      if (rout>=rin) then
        inrange=.true.
        return
      endif  
      if (abs(rout-rin).le.tol) then
        onrange=.true.
        return
      endif
      outofbounds=.true.      
  end subroutine  checkrange
END Module interpolation_mod
