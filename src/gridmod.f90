MODULE gridmod
  USE globalmath

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: lineargrid=1  ! r(i)=h*(i-1)
  INTEGER, PARAMETER, PRIVATE :: loggrid=2     ! r(i)=r0*(exp(h*(i-1))-1)

  TYPE GridInfo
     INTEGER :: TYPE
     INTEGER :: n
     INTEGER :: ishift
     REAL(8) :: h
     REAL(8), POINTER :: r(:)
     REAL(8), POINTER :: drdu(:)    ! for loggrid -- dr/du
     REAL(8), POINTER :: pref(:)    ! for loggrid -- r0*exp(u/2)
     REAL(8), POINTER :: rr02(:)    ! for loggrid -- (r+r0)**2
  END TYPE GridInfo

CONTAINS
  !**********************************************************************
  ! function usingloggrid(Grid)
  !**********************************************************************
  FUNCTION usingloggrid(Grid)
    LOGICAL :: usingloggrid
    TYPE(GridInfo), INTENT(IN) :: Grid

    usingloggrid=.FALSE.
    IF (Grid%type==loggrid) usingloggrid=.TRUE.
  END FUNCTION usingloggrid

  !**********************************************************************
  FUNCTION overint(n,h,f1,icorr)
    !   function to calculate the integral of one vectors f1
    !      using simpsons rule assuming a regular grid with
    !      spacing of h and n total points
    !
    !      icorr: optional parameter: used only when n is even
    !             if icorr<0,  a trapezoidal correction is applied
    !                          at the start of interval
    !             if icorr>=0, a trapezoidal correction is applied
    !                          at the end of interval
    !             default (if missing) is icorr=0
    REAL(8) :: overint
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: h,f1(:)
    INTEGER, OPTIONAL :: icorr


    REAL(8),PARAMETER :: tol=1.D-14
    INTEGER :: i,j,istart,m

    overint=0

    !Eliminate zeros at end of interval
    i=n;DO WHILE(ABS(f1(i))<machine_zero.AND.i>2);i=i-1;ENDDO
    m=MIN(i+1,n)

    IF (m<=1) THEN
       RETURN
    ELSEIF (m==2) THEN
       overint=(f1(1)+f1(2))*(h/2)   ! Trapezoidal rule
       RETURN
    ENDIF

    istart=1
    IF (PRESENT(icorr)) THEN
       IF (icorr<0.AND.MOD(m,2)==0) istart=2
    ENDIF
    overint=f1(istart)+4*f1(istart+1)+f1(istart+2)
    j=((m-istart)/2)*2+istart
    IF (j>=istart+4) THEN
       DO i=istart+4,j,2
          overint=overint+f1(i-2)+4*f1(i-1)+f1(i)
       ENDDO
    ENDIF
    overint=overint*(h/3)
    IF (m>j) overint=overint+(f1(j)+f1(m))*(h/2)
    IF (istart==2) overint=overint+(f1(1)+f1(2))*(h/2)
    RETURN
  END FUNCTION overint

  !*******************************************************************
  ! function integrator(Grid,arg)
  !*******************************************************************
  FUNCTION integrator(Grid,arg,str,fin)
    REAL(8) :: integrator
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: arg(:)
    INTEGER, INTENT(IN), OPTIONAL :: str,fin

    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: i,n,i1,i2

    n=Grid%n
    i1=1;i2=n

    IF (PRESENT(str).AND.PRESENT(fin)) THEN
       i1=str; i2=fin; n=i2-i1+1
    ENDIF

    SELECT CASE(Grid%type)
    CASE default
       WRITE(6,*) 'Error in integrator -- grid ', Grid%type
       STOP
    CASE(lineargrid)
       integrator=overint(n,Grid%h,arg(i1:i2))
    CASE(loggrid)
       ALLOCATE(dum(i1:i2),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in integrator -- allocation ', Grid%n,i
          STOP
       ENDIF
       dum(i1:i2)=arg(i1:i2)*Grid%drdu(i1:i2)
       integrator=overint(n,Grid%h,dum(i1:i2),-1)
       DEALLOCATE(dum)
    END SELECT

  END FUNCTION integrator

  !**********************************************************************
  SUBROUTINE INTWGT(Grid,wgt,str,fin,icorr)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(INOUT) :: wgt(:)
    INTEGER, INTENT(IN), OPTIONAL :: str,fin,icorr


    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: i,j,n,i1,i2,istart

    n=Grid%n
    i1=1;i2=n

    IF (PRESENT(str).AND.PRESENT(fin)) THEN
       i1=str; i2=fin; n=i2-i1+1
    ENDIF

    wgt=0.d0

    IF (n<=1) THEN
       RETURN
    ELSEIF (n==2) THEN
       wgt(i1)=Grid%h/2;   wgt(i2)=Grid%h/2  ! Trapezoidal rule
       if (Grid%type==loggrid) wgt(i1:i2)=wgt(i1:i2)*Grid%drdu(i1:i2)
       RETURN
    ENDIF

    istart=i1
    IF (PRESENT(icorr)) THEN
       IF (icorr<0.AND.MOD(n,2)==0) istart=i1+1
    ENDIF
    wgt(istart)=1; wgt(istart+1)=4; wgt(istart+2)=1
    j=((n-istart)/2)*2+istart
    IF (j>=istart+4) THEN
       DO i=istart+4,j,2
          wgt(i-2)=wgt(i-2)+1
          wgt(i-1)=4
          wgt(i)=1
       ENDDO
    ENDIF
    wgt=wgt*Grid%h/3
    IF (n>j) then
         wgt(j)=wgt(j)+Grid%h/2; wgt(n)=Grid%h/2
    ENDIF

    IF (istart>i1) then
       wgt(i1)=Grid%h/2; wgt(istart)=wgt(istart)+Grid%h/2
    ENDIF
    if (Grid%type==loggrid) wgt(i1:i2)=wgt(i1:i2)*Grid%drdu(i1:i2)

  END SUBROUTINE INTWGT

  !*****************************************************************
  FUNCTION overlap(Grid,f1,f2,str,fin)
    !   function to calculate the overlap between two vectors f1 and f2
    REAL(8) :: overlap
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: f1(:),f2(:)
    INTEGER, INTENT(IN), OPTIONAL :: str,fin

    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: i,n,i1,i2

    n=Grid%n
    i1=1;i2=n
    IF (PRESENT(str).AND.PRESENT(fin)) THEN
       i1=str; i2=fin; n=i2-i1+1
    ENDIF


    ALLOCATE(dum(i1:i2),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in overlap allocation ', n,i
       STOP
    ENDIF
    dum(1:n)=f1(i1:i2)*f2(i1:i2)
    overlap=integrator(Grid,dum(1:n),1,n)
    DEALLOCATE(dum)
  END FUNCTION overlap
  !*****************************************************

  SUBROUTINE nderiv(h,y,z,ndim,ierr)
    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(INOUT) :: ierr
    REAL(8) , INTENT(IN) :: h,y(:)
    REAL(8) , INTENT(INOUT) :: z(:)
    !     subroutine ddet5(h,y,z,ndim)
    !      ssp routine modified by nawh 6/8/76

    REAL(8) :: hh,yy,a,b,c
    INTEGER :: i

    ierr=-1
    IF (ndim.LT.5) RETURN
    !        prepare differentiation loop
    hh=.08333333333333333d0/h
    yy=y(ndim-4)
    b=hh*(-25.d0*y(1)+48.d0*y(2)-36.d0*y(3)+16.d0*y(4)-3.d0*y(5))
    c=hh*(-3.d0*y(1)-10.d0*y(2)+18.d0*y(3)-6.d0*y(4)+y(5))

    !
    !        start differentiation loop
    DO  i=5,ndim
       a=b
       b=c
       c=hh*(y(i-4)-y(i)+8.d0*(y(i-1)-y(i-3)))
       z(i-4)=a
    ENDDO
    !        end of differentiation loop
    !
    !        normal exit
    a=hh*(-yy+6.d0*y(ndim-3)-18.d0*y(ndim-2)+10.d0*y(ndim-1)          &
&        +3.d0*y(ndim))
    z(ndim)=hh*(3.d0*yy-16.d0*y(ndim-3)+36.d0*y(ndim-2)               &
&        -48.d0*y(ndim-1)+25.d0*y(ndim))
    z(ndim-1)=a
    z(ndim-2)=c
    z(ndim-3)=b
    !
    ierr=0
    RETURN
  END SUBROUTINE nderiv

  !**********************************************************************
  ! subroutine derivative(Grid,f,dfdr)
  !*********************************************************************
  SUBROUTINE derivative(Grid,f,dfdr,begin,bend)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: f(:)
    REAL(8), INTENT(OUT) :: dfdr(:)
    INTEGER, OPTIONAL, INTENT(IN) :: begin,bend

    INTEGER :: i,n,i1,i2

    i1=1;i2=Grid%n;n=i2-i1+1
    IF (PRESENT(begin).OR.PRESENT(bend)) THEN
       IF (begin>=1.AND.bend<= Grid%n) THEN
          i1=begin;i2=bend;n=i2-i1+1
       ELSE
          WRITE(6,*) 'Error in derivative', begin,bend,Grid%n
          STOP
       ENDIF
    ENDIF


    SELECT CASE(Grid%type)
    CASE default
       WRITE(6,*) 'Error in derivative -- grid ', Grid%type
       STOP
    CASE(lineargrid)
       CALL nderiv(Grid%h,f(i1:i2),dfdr(i1:i2),n,i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in derivative -nderiv problem', i
          STOP
       ENDIF
    CASE(loggrid)
       CALL nderiv(Grid%h,f(i1:i2),dfdr(i1:i2),n,i)
       IF (i/=0) THEN
          WRITE(6,*) 'Error in derivative -nderiv problem', i
          STOP
       ENDIF
       dfdr(i1:i2)=dfdr(i1:i2)/Grid%drdu(i1:i2)
    END SELECT

  END SUBROUTINE derivative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! simplederiv(Grid,f,dfdr,begin,bend)
  !   low order formula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE simplederiv(Grid,f,dfdr,begin,bend)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: f(:)
    REAL(8), INTENT(OUT) :: dfdr(:)
    INTEGER, OPTIONAL, INTENT(IN) :: begin,bend

    INTEGER :: i,n,i1,i2
    REAL(8) :: HH

    i1=1;i2=Grid%n;n=i2-i1+1
    IF (PRESENT(begin).OR.PRESENT(bend)) THEN
       IF (begin>=1.AND.bend<= Grid%n) THEN
          i1=begin;i2=bend;n=i2-i1+1
          IF (n<3) THEN
             WRITE(6,*) 'Error in simplederiv -- n too small',n,i1,i2
             STOP
          ENDIF
       ELSE
          WRITE(6,*) 'Error in simplederive', begin,bend,Grid%n
          STOP
       ENDIF
    ENDIF

    HH=0.5d0/Grid%h
    dfdr(i1)=HH*(-3*f(i1)+4*f(i1+1)-f(i1+2))
    DO i=i1,i2-2
       dfdr(i+1)=HH*(f(i+2)-f(i))
    ENDDO
    dfdr(i2)=HH*(3*f(i2)-4*f(i2-1)+f(i2-2))

    IF (Grid%type==loggrid) THEN
       dfdr(i1:i2)=dfdr(i1:i2)/Grid%drdu(i1:i2)
    ENDIF

  END SUBROUTINE simplederiv

  SUBROUTINE laplacian(Grid,l,wfn,del,fin)
    TYPE(GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: wfn(:)
    REAL(8), INTENT(INOUT) :: del(:)
    INTEGER, INTENT(IN),OPTIONAL :: fin

    INTEGER :: OK,i,lfac,n
    REAL(8), ALLOCATABLE :: dum1(:),dum2(:)

    n=Grid%n
    IF (PRESENT(fin)) n=fin
    ALLOCATE(dum1(n),dum2(n),stat=ok)
    IF (ok /= 0) THEN
       WRITE(6,*) 'Error in laplace allocation',n
       STOP
    ENDIF
    lfac=l*(l+1)
    del=0
    DO i=2,n
       dum2(i)=wfn(i)/Grid%r(i)
       del(i)=-lfac*dum2(i)
    ENDDO
    CALL extrapolate(Grid,dum2)
    CALL derivative(Grid,dum2,dum1,1,n)
    del(:)=del(:)+2*Grid%r(:)*dum1(:)
    CALL derivative(Grid,dum1,dum2,1,n)
    del(:)=del(:)+(Grid%r(:)**2)*dum2(:)

    DEALLOCATE(dum1,dum2)
  END SUBROUTINE laplacian

  !******************************************************
  ! SUBROUTINE poisson(Grid,q,den,rv,ecoul,v00)
  !*****************************************************
  SUBROUTINE poisson(Grid,q,den,rv,ecoul,v00)
    !  use Numerov algorithm to solve poisson equation
    !  den(n) is electron density * (4*pi*r**2)
    !  rv(n) is returned as electrostatic potential * r
    !  ecoul is the coulomb interaction energy

    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN):: den(:)
    REAL(8), INTENT(INOUT) :: rv(:),ecoul,q
    REAL(8), OPTIONAL, INTENT(OUT) :: v00

    REAL(8), ALLOCATABLE :: a(:),b(:)
    REAL(8) :: sd,sl,h,h2
    INTEGER :: i,n

    n=Grid%n
    h=Grid%h

    rv=0.d0
    q=integrator(Grid,den)
    write(6,*) 'In poisson q = ', q; call flush(6)

    ALLOCATE(a(n),b(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in allocating arrays in poisson',i,n
       STOP
    ENDIF

    IF (Grid%type==lineargrid) THEN

       sd=2
       sl=-1
       a(1)=0
       DO i=2,n
          a(i)=h*den(i)/(6*(i-1))
       ENDDO
       rv(1)=0
       rv(2)=10*a(2)+a(3)
       DO i=3,n-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       rv(n)=10*a(n)+a(n-1)+2*q

    ELSEIF (Grid%type==loggrid) THEN

       sd=2+10*h*h/48
       sl=-1+h*h/48
       a(1)=0
       h2=h*h
       DO i=2,n
          a(i)=h2*Grid%rr02(i)*den(i)/(6*Grid%r(i))/Grid%pref(i)
       ENDDO
       rv(1)=0
       rv(2)=10*a(2)+a(3)
       DO i=3,n-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       !   last term is boundary value at point n+1
       rv(n)=10*a(n)+a(n-1)-2*q*sl/(Grid%pref(n)*EXP(h/2))

    ENDIF

    CALL conthomas(n,sl,sd,rv)
    IF (Grid%type==loggrid) rv=rv*Grid%pref

    !
    !  calculate ecoul
    !
    DO i=2,n
       a(i)=den(i)*rv(i)/Grid%r(i)
    ENDDO
    a(1)=0
    ecoul=0.5d0*integrator(Grid,a)
    !WRITE(6,*) ' from poisson: ecoul = ',ecoul

    IF (PRESENT(v00)) THEN
       a=0
       a(2:n)=den(2:n)/Grid%r(2:n)
       v00=2*integrator(Grid,a)
    ENDIF
    DEALLOCATE(a,b)
    !
  END SUBROUTINE poisson

  !*****************************************************************
  !Alternative form of poisson solver written by Marc Torrent 6/9/06
  !  works well for loggrid
  !******************************************************************
  SUBROUTINE poisson_marc(Grid,q,den,rv,ecoul)
    !  use Numerov algorithm to solve poisson equation
    !  den(n) is electron density * (4*pi*r**2)
    !  rv(n) is returned as electrostatic potential * r
    !  ecoul is the coulomb interation energy

    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN):: den(:)
    REAL(8), INTENT(INOUT) :: rv(:),ecoul,q

    REAL(8), ALLOCATABLE :: aa(:),bb(:),cc(:),dd(:)
    REAL(8) :: sd,sl,h,h2
    INTEGER :: i,n,ir,jr

    n=Grid%n
    h=Grid%h

    rv=0.d0
    q=integrator(Grid,den)

    ALLOCATE(aa(n),bb(n),cc(n),dd(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in allocating arrays in poisson',i,n
       STOP
    ENDIF

    DO jr=n,2,-1
       ir=n-jr+1
       aa(ir)=den(jr)*Grid%drdu(jr)
       bb(ir)=den(jr)*Grid%drdu(jr)/Grid%r(jr)
    END DO

    cc=0.d0
    cc(5)=aa(n-4);cc(4)=aa(n-3);cc(3)=aa(n-2);cc(2)=aa(n-1)
    call extrapolate(Grid,cc);aa(n)=cc(1)
    cc(5)=bb(n-4);cc(4)=bb(n-3);cc(3)=bb(n-2);cc(2)=bb(n-1)
    call extrapolate(Grid,cc);bb(n)=cc(1)

    cc(1)=0.d0;dd(1)=0.d0
    DO ir=3,n,2
       cc(ir)  =cc(ir-2)+h/3.d0*(aa(ir-2)+4.d0*aa(ir-1)+aa(ir))
       cc(ir-1)=cc(ir-2)+h/3.d0*(1.25d0*aa(ir-2)+2.0d0*aa(ir-1)-0.25d0*aa(ir))
       dd(ir)  =dd(ir-2)+h/3.d0*(bb(ir-2)+4.d0*bb(ir-1)+bb(ir))
       dd(ir-1)=dd(ir-2)+h/3.d0*(1.25d0*bb(ir-2)+2.d0*bb(ir-1)-0.25d0*bb(ir))
    END DO
    IF (MOD(n,2)==0) THEN
       cc(n)=cc(n-2)+h/3.d0*(aa(n-2)+4.d0*aa(n-1)+aa(n))
       dd(n)=dd(n-2)+h/3.d0*(bb(n-2)+4.d0*bb(n-1)+bb(n))
    END IF

    rv(1)=0.d0
    DO ir=2,n
       jr=n-ir+1
       rv(ir)=2.d0*(dd(jr)*Grid%r(ir)+(cc(n)-cc(jr)))
    END DO

    !  calculate ecoul
    aa(1)=0.d0
    do i=2,n
     aa(i)=den(i)*rv(i)/Grid%r(i)
    end do
    ecoul=0.5d0*integrator(Grid,aa)

    DEALLOCATE(aa,bb,cc,dd)

  END SUBROUTINE poisson_marc

  !********************************************************************
  !  use Numerov algorithm to solve poisson equation
  !  for angularly dependent charge distribution of angular momentum l
  !  den(n) is electron density * (4*pi*r**2) appropriate for l
  !  rv(n) is returned as electrostatic potential * r
  !  a(n), b(n), and c(n) are work arrays
  !********************************************************************

  SUBROUTINE apoisson(Grid,l,irc,den,rv)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,irc
    REAL(8), INTENT(IN) :: den(:)
    REAL(8), INTENT(OUT) :: rv(:)

    INTEGER :: i,j
    REAL(8) :: angm,r,q,h,h2
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:)


    ALLOCATE(a(irc),b(irc),c(irc),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in allocating arrays in apoisson ',i,irc
       STOP
    ENDIF

    b=den(1:irc)*(Grid%r(1:irc))**L
    q=integrator(Grid,b,1,irc)/(2*l+1)
    h=Grid%h
    WRITE(6,*) 'check l multipole',l,q
    IF (Grid%type==lineargrid) THEN
       a=0;b=0; c=0; rv=0
       DO i=2,irc
          a(i)=0.2d0*h*den(i)/((i-1))
       ENDDO
       DO i=2,irc-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       !  set up Tridiagonal equations
       a=0; angm=l*(l+1)
       DO i=2,irc-1
          IF (i>2) a(i)=-1.2d0+0.1d0*angm/((i-2)**2)
          c(i)=-1.2d0+0.1d0*angm/((i)**2)
          b(i)=2.4d0+angm/((i-1)**2)
       ENDDO
       a(2)=0
       IF (l==1) b(2)=b(2)+0.2d0
       rv(irc)=2*q/(Grid%r(irc)**l)
       rv(irc-1)=rv(irc-1)-c(irc-1)*rv(irc)
       c(irc-1)=0;rv(1)=0
       CALL thomas(irc-2,a(2:irc-1),b(2:irc-1),c(2:irc-1),rv(2:irc-1))
    ELSEIF (Grid%type==loggrid) THEN
       a=0;b=0; c=0; rv=0; h2=h*h
       DO i=2,irc
          a(i)=0.2d0*h2*Grid%rr02(i)*den(i)/(Grid%pref(i)*Grid%r(i))
       ENDDO
       DO i=2,irc-1
          rv(i)=10*a(i)+a(i+1)+a(i-1)
       ENDDO
       !  set up Tridiagonal equations
       a=0; angm=l*(l+1)
       DO i=2,irc-1
          IF (i>2) a(i)=-1.2d0+&
&              0.1d0*h2*(0.25d0+angm*Grid%rr02(i-1)/Grid%r(i-1)**2)
          c(i)=-1.2d0+0.1d0*h2*(0.25d0+angm*Grid%rr02(i+1)/Grid%r(i+1)**2)
          b(i)=2.4d0+h2*(0.25d0+angm*Grid%rr02(i)/Grid%r(i)**2)
       ENDDO
       a(2)=0
       !IF (l==1) b(2)=b(2)+0.2d0*Grid%rr02(1)/(Grid%r(2)**2)  Slight error
       IF (l==1) b(2)=b(2)+0.2d0*h2*Grid%rr02(1)/(Grid%r(2)**2)
       rv(irc)=2*q/(Grid%r(irc)**l)/Grid%pref(irc)
       rv(irc-1)=rv(irc-1)-c(irc-1)*rv(irc)
       c(irc-1)=0;rv(1)=0
       CALL thomas(irc-2,a(2:irc-1),b(2:irc-1),c(2:irc-1),rv(2:irc-1))
       rv(1:irc)=Grid%pref(1:irc)*rv(1:irc)
    ENDIF
    !
    DEALLOCATE(a,b,c)
  END SUBROUTINE apoisson

  !******************************************************************
  ! pgm to determine r=0 form of potential assuming that
  !   nuclear contribution (-2*nz) is not yet included
  !******************************************************************

  SUBROUTINE zeropot(Grid,rv,v0,v0p)
    ! extrapolate potential to value at r=0
    TYPE (GridInfo), INTENT(IN):: Grid
    REAL(8), INTENT(IN) :: rv(:)    ! Note: rv(1) corresponds to r=0
    REAL(8), INTENT(OUT) :: v0,v0p

    REAL(8) :: tmp(15),tmp1(15)
    INTEGER :: i,n

    tmp(2:15)=rv(2:15)/Grid%r(2:15)
    CALL extrapolate(Grid,tmp(1:15))
    v0=tmp(1)

    CALL derivative(Grid,tmp(1:15),tmp1(1:15),2,15)

    CALL extrapolate(Grid,tmp1(1:15))
    v0p=tmp1(1)
  END SUBROUTINE zeropot

  SUBROUTINE extrapolate(Grid,v)
    ! extrapolate array v to r=0 at v(1)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(INOUT) :: v(:)  ! assume v(2),v(3)...  given

    !v(1)=v(4)+3.d0*(v(2)-v(3))                         !second order formula
    !v(1))=4.d0*v(2)-6.d0*v(3)+4.d0*v(4)-v(5)           ! third order formula
    v(1)=5.d0*v(2)-10.d0*v(3)+10.d0*v(4)-5.d0*v(5)+v(6) ! fourth order formula

  END SUBROUTINE extrapolate

  !*************************************************************
  ! subroutine forward_numerov(Grid,l,many,energy,rv,zeroval,wfn,nodes)
  !*************************************************************
  SUBROUTINE forward_numerov(Grid,l,many,energy,rv,zeroval,wfn,nodes)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,many
    REAL(8), INTENT(IN) :: energy,zeroval,rv(:)
    ! zeroval == lim r-->0 A(r)*P(r)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(1) and wfn(2) given
    ! on ouput wfn(i) given for
    !      i<=many
    INTEGER, INTENT(OUT) :: nodes

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    ALLOCATE(a(many),b(many),p(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error forward_numerov ', many,i
       STOP
    ENDIF

    p(1)=wfn(1)
    p(2)=wfn(2)
    xx=zeroval
    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=2,many
       a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(1:2)=wfn(1:2)/Grid%pref(1:2)
       xx=Grid%rr02(1)*xx/Grid%pref(1)
       a=0.25d0+Grid%rr02(1:many)*a
    ENDIF

    b=2.4d0+h2*a
    a=1.2d0-0.1d0*h2*a
    p(3)=(b(2)*p(2)+0.1d0*h2*xx)/a(3)

    nodes=0
    DO i=4,many
       p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2))/a(i)
       IF (p(i)*p(i-1) < 0.d0) nodes=nodes+1
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(1:i)=scale*p(1:i)
       ENDIF
    ENDDO

    wfn(1:many)=p(1:many)

    IF (Grid%type==loggrid) THEN
       wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE forward_numerov

  !*************************************************************
  ! subroutine shifted_forward_numerov(Grid,many,istart,ww,wfn,nodes)
  !     designed for use with scalar relativistic case in which
  !       l and energy information is stored in ww
  !*************************************************************
  SUBROUTINE shifted_forward_numerov(Grid,many,istart,ww,wfn,nodes)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: many,istart
    REAL(8), INTENT(IN) :: ww(:)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(1:5) given
    ! on ouput wfn(i) given for
    !      i<=many
    INTEGER, INTENT(OUT) :: nodes

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    IF (istart>many) THEN
       WRITE(6,*) 'shifted_forward_numerov:  error istart many',istart,many
       STOP
    ENDIF
    ALLOCATE(a(many),b(many),p(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error forward_numerov ', many,i
       STOP
    ENDIF

    p(1:istart)=wfn(1:istart)
    a=0
    h=Grid%h;    h2=h*h
    DO i=2,many
       a(i)=ww(i)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(1:istart)=wfn(1:istart)/Grid%pref(1:istart)
       a=0.25d0+Grid%rr02(1:many)*a
    ENDIF

    b=2.4d0+h2*a
    a=1.2d0-0.1d0*h2*a

    nodes=0
    DO i=istart+1,many
       p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2))/a(i)
       IF (p(i)*p(i-1) < 0.d0) nodes=nodes+1
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(1:i)=scale*p(1:i)
       ENDIF
    ENDDO

    wfn(1:many)=p(1:many)

    IF (Grid%type==loggrid) THEN
       wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE shifted_forward_numerov
  !**********************************************************************
  !  pgm to integrate outward the radial schroedinger inhomogeneous equation
  !    at energy 'energy' and at angular momentum l
  !    with potential smooth rv/r
  !        proj/r == inhomogeneous function
  !     It is assumed that for r~0, proj~~(r**(l+1))*(c0+r**2*c2+...)
  !       and wfn~C*r**(l+3)*polynomial(r) for r~0;
  !     Works well for small ranges, but diverges at long range
  !  uses Noumerov algorithm
  !*************************************************************************
  !*************************************************************
  ! subroutine inhomogeneous_numerov(Grid,l,many,energy,rv,proj,wfn)
  !*************************************************************
  SUBROUTINE inhomogeneous_numerov(Grid,l,many,energy,rv,proj,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,many
    REAL(8), INTENT(IN) :: energy,rv(:),proj(:)
    REAL(8), INTENT(OUT) :: wfn(:)
    ! initial values of wfn determined from proj(r=0)
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    INTEGER :: i,j,k,n
    INTEGER :: count=0

    ALLOCATE(a(many),b(many),c(many),p(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error forward_numerov ', many,i
       STOP
    ENDIF

     !Note for NAWH:  Need to check these initial values....

    a=0
    a(2:8)=proj(2:8)/(Grid%r(2:8)**(l+1))
    CALL extrapolate(Grid,a)
    !write(6,'("extrapolate ",1p,9e15.7)') a(1:8)
    wfn=0
    wfn(2)=-a(1)*(Grid%r(2)**(l+3))/(4*l+6.d0)

    p(1)=wfn(1)
    p(2)=wfn(2)
    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=2,many
       a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    b(1:many)=0.1d0*h2*proj(1:many)
    IF (Grid%type==loggrid) THEN
       p(1:2)=wfn(1:2)/Grid%pref(1:2)
       a=0.25d0+Grid%rr02(1:many)*a
       b=Grid%rr02(1:many)*b/Grid%pref(1:many)
    ENDIF

    c=0
    DO i=2,many-1
       c(i)=10*b(i)+b(i-1)+b(i+1)
    ENDDO

    b=2.4d0+h2*a
    a=1.2d0-0.1d0*h2*a
    p(3)=(b(2)*p(2)-c(2))/a(3)

    DO i=4,many
       p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2)-c(i-1))/a(i)
    ENDDO

    wfn(1:many)=p(1:many)

    IF (Grid%type==loggrid) THEN
       wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
    ENDIF

    DEALLOCATE(a,b,c,p)

  END SUBROUTINE inhomogeneous_numerov

  !*************************************************************
  ! subroutine inhomo_bound_numerov(Grid,l,mn,energy,rv,rhs,wfn)
  !      Assumes that wfn(mn+1)==0
  !      Method can be unstable due to singularity of matrix inversion
  !        introduced by homogeneous solution
  !*************************************************************
  SUBROUTINE inhomo_bound_numerov(Grid,l,mn,energy,rv,rhs,wfn,phi)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,mn
    REAL(8), INTENT(IN) :: energy,rv(:),rhs(:)
    REAL(8), INTENT(OUT) :: wfn(:)
    REAL(8), INTENT(IN), OPTIONAL :: phi(:)

    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale,zeroval,st
    INTEGER :: i,j,k,n

    ALLOCATE(a(mn-1),b(mn-1),c(mn-1),p(mn-1),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error inhomo_bound_numerov ', mn,i
       STOP
    ENDIF

    ! assume wfn(r) ~ r**(l+1) for r->0 for homogeneous solution
    !

    zeroval=0

    IF (l==0)  zeroval=rv(1)
    IF (l==1)  zeroval=2
    st=rv(1)/2; st=st/(l+1)

    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=1,mn-1
       a(i)=rv(i+1)/Grid%r(i+1)-energy+angm/(Grid%r(i+1)**2)
    ENDDO
    b(1:mn-1)=0.1d0*h2*rhs(2:mn)
    xx=zeroval
    IF (Grid%type==loggrid) THEN
       a(1:mn-1)=0.25d0+Grid%rr02(2:mn)*a(1:mn-1)
       b(1:mn-1)=Grid%rr02(2:mn)*b(1:mn-1)/Grid%pref(2:mn)
       xx=Grid%rr02(1)*xx/Grid%pref(1)
    ENDIF

    c=0; c(1)=10*b(1)+b(2); c(mn-1)=10*b(mn-1)+b(mn-2)
    DO i=2,mn-2
       c(i)=10*b(i)+b(i-1)+b(i+1)
    ENDDO

    b=-2.4d0-h2*a   ;
    ! zero value extrapolation
    b(1)=b(1)-0.1d0*h2*xx/((Grid%r(2)**(l+1))*(1.d0+st*Grid%r(2)))
    a=1.2d0-0.1d0*h2*a
    p=0;p(2:mn-1)=a(1:mn-2)
    a(1:mn-2)=a(2:mn-1);a(mn-1)=0

    CALL thomas(mn-1,p,b,a,c)
    wfn=0
    wfn(2:mn)=c(1:mn-1)

    IF (Grid%type==loggrid) THEN
       wfn(1:mn)=wfn(1:mn)*Grid%pref(1:mn)
    ENDIF

    if (present(phi)) then   !  orthogonalize wfn to phi
       wfn=wfn-phi*overlap(Grid,phi,wfn)/overlap(Grid,phi,phi)
    endif

    DEALLOCATE(a,b,c,p)

  END SUBROUTINE inhomo_bound_numerov

  !*************************************************************
  ! subroutine inhomo_numerov_coeff(Grid,l,mn,energy,rv,rhs,a,b)
  !      Returns coefficients to solve
  ![   d^2     l(l+1)   rv(r)         ]
  ![ - ---  +  ------ + ----  -energy ] g(r) =  rhs(r)*(unknown func(r))
  ![   dr^2      r^2      r           ]
  !
  !  coefficients a(i), b(i), etc. indexed to Grid%r(i)
  !  It is assumed that g(mn)=0 and rhs(mn)=0
  !*************************************************************
  SUBROUTINE inhomo_numerov_coeff(Grid,l,mn,energy,rv,rhs,a,b,c,d)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,mn
    REAL(8), INTENT(IN) :: energy,rv(:),rhs(:)
    REAL(8), INTENT(OUT) :: a(:),b(:),c(:),d(:)
    REAL(8) :: xx,angm,h,h2,scale,zeroval,st
    INTEGER :: i,j,k,n

    ! assume g(r) ~ r**(l+3) for r->0

    !zeroval=0

    !IF (l==0)  zeroval=rv(1)
    !IF (l==1)  zeroval=2
    !st=rv(1)/2; st=st/(l+1)

    angm=l*(l+1)
    a=0; b=0; c=0;   d=0
    h=Grid%h;    h2=h*h
    DO i=2,mn
       b(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    d(2:mn)=0.1d0*h2*rhs(2:mn)
    !xx=zeroval
    IF (Grid%type==loggrid) THEN
       b(2:mn)=0.25d0+Grid%rr02(2:mn)*b(2:mn)
       d(2:mn)=Grid%rr02(2:mn)*d(2:mn)/Grid%pref(2:mn)
       !xx=Grid%rr02(1)*xx/Grid%pref(1)
    ENDIF

    a(2:mn)=-2.4d0-h2*b(2:mn)   ;
    ! zero value extrapolation       ! not used
    !a(2)=a(2)-0.1d0*h2*xx/((Grid%r(2)**(l+1))*(1.d0+st*Grid%r(2)))
    b(2:mn)=1.2d0-0.1d0*h2*b(2:mn)
    c(2:mn)=10*d(2:mn)
    If (Grid%type==loggrid) then
         a(2:mn)=a(2:mn)/Grid%pref(2:mn)
         b(2:mn)=b(2:mn)/Grid%pref(2:mn)
    endif

  END SUBROUTINE inhomo_numerov_coeff

  !*************************************************************
  ! subroutine inhomo_numerov_SVD(Grid,l,nr,energy,rv,rhs,sol)
  !      Returns sol(r) for equation
  ![   d^2     l(l+1)   rv(r)         ]
  ![ - ---  +  ------ + ----  -energy ] sol(r) =  -rhs(r)
  ![   dr^2      r^2      r           ]
  !
  !  It is assumed that sol(nr)=0 and rhs(nr)=0; in general nr=n
  !*************************************************************
  SUBROUTINE inhomo_numerov_SVD(Grid,l,nr,energy,tol,rv,rhs,sol,phi)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,nr
    REAL(8), INTENT(IN) :: energy,tol,rv(:),rhs(:)
    REAL(8), INTENT(OUT) :: sol(:)
    REAL(8), INTENT(IN), OPTIONAL :: phi(:)

    REAL(8) :: xx,angm,h,h2,scale,zeroval,st
    INTEGER :: i,j,k,n,ok,nn,LWORK
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),d(:)
    REAL(8), ALLOCATABLE :: MN(:,:),u(:,:),vt(:,:),s(:),work(:)
    INTEGER, ALLOCATABLE :: iwork(:)

    n=Grid%n;  nn=nr-1; LWORK=8*nn*nn
    ALLOCATE (a(n),b(n),c(n),d(n),u(nn,nn),vt(nn,nn),s(nn),&
&        MN(nn,nn),work(LWORK),IWORK(8*nn) )

    ! assume rhs(r) ~ r**(l+3) for r->0

    zeroval=0

    IF (l==0)  zeroval=rv(1)
    IF (l==1)  zeroval=2
    st=rv(1)/2; st=st/(l+1)

    angm=l*(l+1)
    a=0; b=0; c=0;   d=0
    h=Grid%h;    h2=h*h
    DO i=2,nr
       b(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    d(2:nr)=0.1d0*h2*rhs(2:nr)
    xx=zeroval
    IF (Grid%type==loggrid) THEN
       b(2:nr)=0.25d0+Grid%rr02(2:nr)*b(2:nr)
       d(2:nr)=Grid%rr02(2:nr)*d(2:nr)/Grid%pref(2:nr)
       xx=Grid%rr02(1)*xx/Grid%pref(1)
    ENDIF

    a(2:nr)=-2.4d0-h2*b(2:nr)   ;
    ! zero value extrapolation      ! not used
    a(2)=a(2)-0.1d0*h2*xx/((Grid%r(2)**(l+1))*(1.d0+st*Grid%r(2)))
    b(2:nr)=1.2d0-0.1d0*h2*b(2:nr)
    c=0; c(2)=10*d(2)+d(3); c(nr)=10*d(nr)+d(nr-1);
    do i=3,nr-1
       c(i)=10*d(i)+d(i-1)+d(i+1)
    enddo
    If (Grid%type==loggrid) then
         a(2:nr)=a(2:nr)/Grid%pref(2:nr)
         b(2:nr)=b(2:nr)/Grid%pref(2:nr)
    endif

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
    !WRITE(6,*) 'dgesdd completed for MN with info = ', ok

    sol=0; scale=tol
    DO i=1,nn
          !WRITE(6,*) 'i s = ', i,s(i)
       IF (s(i)>scale) THEN
          xx=Dot_Product(c(2:nr),u(:,i))/s(i)
          sol(2:nr)=sol(2:nr)+xx*vt(i,:)
       ELSE
          WRITE(6,*) 'i s = ', i,s(i)
       ENDIF
    ENDDO

    if (present(phi)) then   !  orthogonalize sol to phi
       sol=sol-phi*overlap(Grid,phi,sol)/overlap(Grid,phi,phi)
    endif

    DEALLOCATE (a,b,c,d,u,vt,s,work,IWORK,MN )
  END SUBROUTINE inhomo_numerov_SVD

  !*************************************************************
  ! subroutine inhomo_numerov_SVD_bv(Grid,l,nr,energy,tol,rv,rhs,sol,phi)
  !      Returns sol(r) for equation
  ![   d^2     l(l+1)   rv(r)         ]
  ![ - ---  +  ------ + ----  -energy ] sol(r) =  rhs(r)
  ![   dr^2      r^2      r           ]
  !
  !  It is assumed that sol(nr+1)=is given and sol(r) is determined
  !      for r<Grid%r(nr)
  !    rhs(nr+1) is accessed also
  !*************************************************************
  SUBROUTINE inhomo_numerov_SVD_bv(Grid,l,nr,energy,tol,rv,rhs,sol,phi)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,nr
    REAL(8), INTENT(IN) :: energy,tol,rv(:),rhs(:)
    REAL(8), INTENT(INOUT) :: sol(:)
    REAL(8), INTENT(IN), OPTIONAL :: phi(:)

    REAL(8) :: xx,angm,h,h2,scale,zeroval,st
    INTEGER :: i,j,k,n,ok,nn,LWORK
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),d(:)
    REAL(8), ALLOCATABLE :: MN(:,:),u(:,:),vt(:,:),s(:),work(:)
    INTEGER, ALLOCATABLE :: iwork(:)

    n=Grid%n;  nn=nr-1; LWORK=8*nn*nn
    ALLOCATE (a(n),b(n),c(n),d(n),u(nn,nn),vt(nn,nn),s(nn),&
&        MN(nn,nn),work(LWORK),IWORK(8*nn) )

    If (nr>n-1) then
       write(6,*) 'Error in inhomo_numerov_SVD_bv ', n,nr
       stop
    endif

    IF (l==0)  zeroval=rv(1)
    IF (l==1)  zeroval=2
    st=rv(1)/2; st=st/(l+1)

    angm=l*(l+1)
    a=0; b=0; c=0;   d=0
    h=Grid%h;    h2=h*h
    DO i=2,n
       b(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    d(2:n)=0.1d0*h2*rhs(2:n)
    xx=zeroval
    IF (Grid%type==loggrid) THEN
       b(2:n)=0.25d0+Grid%rr02(2:n)*b(2:n)
       d(2:n)=Grid%rr02(2:n)*d(2:n)/Grid%pref(2:n)
       xx=Grid%rr02(1)*xx/Grid%pref(1)
    ENDIF

    a(2:n)=-2.4d0-h2*b(2:n)   ;
    ! zero value extrapolation
    a(2)=a(2)-0.1d0*h2*xx/((Grid%r(2)**(l+1))*(1.d0+st*Grid%r(2)))
    b(2:n)=1.2d0-0.1d0*h2*b(2:n)
    c=0; c(2)=10*d(2)+d(3);
    do i=3,nr
       c(i)=10*d(i)+d(i-1)+d(i+1)
    enddo
    If (Grid%type==loggrid) then
         a(2:n)=a(2:n)/Grid%pref(2:n)
         b(2:n)=b(2:n)/Grid%pref(2:n)
    endif
    c(nr)=c(nr)-b(nr+1)*sol(nr+1)       ! boundary value

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
    !WRITE(6,*) 'dgesdd completed for MN with info = ', ok

    sol(1:nr)=0; scale=tol
    DO i=1,nn
       IF (s(i)>scale) THEN
          xx=Dot_Product(c(2:nr),u(:,i))/s(i)
          sol(2:nr)=sol(2:nr)+xx*vt(i,:)
          !WRITE(6,*) 'i s = ', i,s(i)
       ELSE
          WRITE(6,*) 'i s = ', i,s(i)
       ENDIF
    ENDDO

    if (present(phi)) then   !  orthogonalize sol to phi
       sol=sol-phi*overlap(Grid,phi,sol)/overlap(Grid,phi,phi)
    endif

    DEALLOCATE (a,b,c,d,u,vt,s,work,IWORK,MN )
  END SUBROUTINE inhomo_numerov_SVD_bv


  !*************************************************************
  ! subroutine inhomo_numerov_SVD_bvm(Grid,l,nr,mult,energy,tol,rv,rhs,sol)
  !      Returns sol(r) for equation
  ![   d^2     l(l+1)   rv(r)         ]
  ![ - ---  +  ------ + ----  -energy ] sol(r) =  rhs(r,i)
  ![   dr^2      r^2      r           ]
  !
  !    for k=1,2, .. mult
  !  It is assumed that sol(nr+1)=is given and sol(r) is determined
  !      for r<Grid%r(nr)
  !    rhs(nr+1) is accessed also
  !*************************************************************
  SUBROUTINE inhomo_numerov_SVD_bvm(Grid,l,nr,mult,energy,tol,rv,rhs,sol)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,nr,mult
    REAL(8), INTENT(IN) :: energy,tol,rv(:),rhs(:,:)
    REAL(8), INTENT(INOUT) :: sol(:,:)

    REAL(8) :: xx,angm,h,h2,scale,zeroval,st
    INTEGER :: i,j,k,n,ok,nn,LWORK
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:,:),d(:,:)
    REAL(8), ALLOCATABLE :: MN(:,:),u(:,:),vt(:,:),s(:),work(:)
    INTEGER, ALLOCATABLE :: iwork(:)

    n=Grid%n;  nn=nr-1; LWORK=8*nn*nn
    ALLOCATE (a(n),b(n),c(n,mult),d(n,mult),u(nn,nn),vt(nn,nn),s(nn),&
&        MN(nn,nn),work(LWORK),IWORK(8*nn) )

    If (nr>n-1) then
       write(6,*) 'Error in inhomo_numerov_SVD_bv ', n,nr
       stop
    endif

    IF (l==0)  zeroval=rv(1)
    IF (l==1)  zeroval=2
    st=rv(1)/2; st=st/(l+1)

    angm=l*(l+1)
    a=0; b=0; c=0;   d=0
    h=Grid%h;    h2=h*h
    DO i=2,n
       b(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    d(2:n,:)=0.1d0*h2*rhs(2:n,:)
    xx=zeroval
    IF (Grid%type==loggrid) THEN
       b(2:n)=0.25d0+Grid%rr02(2:n)*b(2:n)
       do i=1,mult
          d(2:n,i)=Grid%rr02(2:n)*d(2:n,i)/Grid%pref(2:n)
       enddo
       xx=Grid%rr02(1)*xx/Grid%pref(1)
    ENDIF

    a(2:n)=-2.4d0-h2*b(2:n)   ;
    ! zero value extrapolation
    a(2)=a(2)-0.1d0*h2*xx/((Grid%r(2)**(l+1))*(1.d0+st*Grid%r(2)))
    b(2:n)=1.2d0-0.1d0*h2*b(2:n)
    c=0; c(2,:)=10*d(2,:)+d(3,:);
    do i=3,nr
       c(i,:)=10*d(i,:)+d(i-1,:)+d(i+1,:)
    enddo
    If (Grid%type==loggrid) then
         a(2:n)=a(2:n)/Grid%pref(2:n)
         b(2:n)=b(2:n)/Grid%pref(2:n)
    endif
    c(nr,:)=c(nr,:)-b(nr+1)*sol(nr+1,:)       ! boundary value

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
    !WRITE(6,*) 'dgesdd completed for MN with info = ', ok

    sol(1:nr,:)=0; scale=tol
    DO i=1,nn
       IF (s(i)>scale) THEN
         do j=1,mult
          xx=Dot_Product(c(2:nr,j),u(:,i))/s(i)
          sol(2:nr,j)=sol(2:nr,j)+xx*vt(i,:)
          !WRITE(6,*) 'i s = ', i,s(i)
         enddo
       ELSE
          WRITE(6,*) 'i s = ', i,s(i)
       ENDIF
    ENDDO

    DEALLOCATE (a,b,c,d,u,vt,s,work,IWORK,MN )
  END SUBROUTINE inhomo_numerov_SVD_bvm


  !*************************************************************
  ! subroutine inhomo_bv_numerov(Grid,l,mn,energy,boundaryv,rv,rhs,wfn)
  !      Assumes that wfn(mn+1)==boundaryv
  !!!!!!!!!!!!!  Not yet checked!!!!!
  !*************************************************************
  SUBROUTINE inhomo_bv_numerov(Grid,l,mn,energy,boundaryv,rv,rhs,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,mn
    REAL(8), INTENT(IN) :: energy,boundaryv,rv(:),rhs(:)
    REAL(8), INTENT(OUT) :: wfn(:)
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale,zeroval,st,bp,bv
    INTEGER :: i,j,k,n

    ALLOCATE(a(mn-1),b(mn-1),c(mn-1),p(mn-1),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error inhomo_bound_numerov ', mn,i
       STOP
    ENDIF

    ! assume wfn(r) ~ r**(l+1) for r->0

    zeroval=0

    IF (l==0)  zeroval=rv(1)
    IF (l==1)  zeroval=2
    st=rv(1)/2; st=st/(l+1)

    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=1,mn-1
       a(i)=rv(i+1)/Grid%r(i+1)-energy+angm/(Grid%r(i+1)**2)
    ENDDO
    bp=rv(mn+1)/Grid%r(mn+1)-energy+angm/(Grid%r(mn+1)**2)
    bv=boundaryv
    b(1:mn-1)=0.1d0*h2*rhs(2:mn)
    xx=zeroval
    IF (Grid%type==loggrid) THEN
       a(1:mn-1)=0.25d0+Grid%rr02(2:mn)*a(1:mn-1)
       bp=0.25d0+Grid%rr02(mn+1)*bp
       b(1:mn-1)=Grid%rr02(2:mn)*b(1:mn-1)/Grid%pref(2:mn)
       xx=Grid%rr02(1)*xx/Grid%pref(1)
       bv=bv/Grid%pref(mn+1)
    ENDIF

    c=0; c(1)=10*b(1)+b(2); c(mn-1)=10*b(mn-1)+b(mn-2)
    c(mn-1)=c(mn-1)-(1.2d0-0.1d0*h2*bp)*bv
    DO i=2,mn-2
       c(i)=10*b(i)+b(i-1)+b(i+1)
    ENDDO

    b=-2.4d0-h2*a   ;
    ! zero value extrapolation
    b(1)=b(1)-0.1d0*h2*xx/((Grid%r(2)**(l+1))*(1.d0+st*Grid%r(2)))
    a=1.2d0-0.1d0*h2*a
    p=0;p(2:mn-1)=a(1:mn-2)
    a(1:mn-2)=a(2:mn-1);a(mn-1)=0

    CALL thomas(mn-1,p,b,a,c)
    wfn=0
    wfn(2:mn)=c(1:mn-1)

    IF (Grid%type==loggrid) THEN
       wfn(1:mn)=wfn(1:mn)*Grid%pref(1:mn)
    ENDIF

    DEALLOCATE(a,b,c,p)

  END SUBROUTINE inhomo_bv_numerov

  !*********************************************************
  ! subroutine backward_numerov(Grid,l,match,energy,rv,wfn)
  !*********************************************************
  SUBROUTINE backward_numerov(Grid,l,match,energy,rv,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,match
    REAL(8), INTENT(IN) :: energy,rv(:)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(n-1) and wfn(n) given
    ! on ouput wfn(i) given for
    !      start<=i<=n

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    n=Grid%n
    ALLOCATE(a(n),b(n),p(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error backward_numerov ', n,i
       STOP
    ENDIF

    p(n)=wfn(n)
    p(n-1)=wfn(n-1)
    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=match,n
       a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(n-1:n)=p(n-1:n)/Grid%pref(n-1:n)
       a(match:n)=0.25d0+Grid%rr02(match:n)*a(match:n)
    ENDIF

    b(match:n)=2.4d0+h2*a(match:n)
    a(match:n)=1.2d0-0.1d0*h2*a(match:n)

    DO i=n-2,match,-1
       p(i)=(b(i+1)*p(i+1)-a(i+2)*p(i+2))/a(i)
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(i:n)=scale*p(i:n)
       ENDIF
    ENDDO

    wfn(match:n)=p(match:n)

    IF (Grid%type==loggrid) THEN
       wfn(match:n)=wfn(match:n)*Grid%pref(match:n)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE backward_numerov

  ! Subroutine from David Vanderbilt's USPS code, modified by Marc
  !     Torrent and Francois Jollet, further modified by NAWH
  !===========================================================================
  !      subroutine cfdsol(zz,yy,jj1,jj2,mesh)
  !===========================================================================

  !     routine for solving coupled first order differential equations
  !
  !      d yy(x,1)
  !      ---------   =  zz(x,1,1) * yy(x,1) + zz(x,1,2) * yy(2,1)
  !         dx
  !
  !      d yy(x,2)
  !      ---------   =  zz(x,2,1) * yy(x,1) + zz(x,2,2) * yy(2,1)
  !         dx
  !
  !
  !     using fifth order predictor corrector algorithm
  !
  !     routine integrates from jj1 to jj2 and can cope with both cases
  !     jj1 < jj2 and jj1 > jj2.  first five starting values of yy must
  !     be provided by the calling program.

  SUBROUTINE cfdsol(Grid,zz,yy,jj1,jj2)
    TYPE(gridinfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN):: zz(:,:,:)
    REAL(8), INTENT(INOUT):: yy(:,:)
    INTEGER, INTENT(IN)  :: jj1,jj2


    REAL(8):: fa(0:5),fb(0:5),abp(1:5),amc(0:4)
    INTEGER :: isgn,i,j,ip,mesh
    REAL(8):: arp,brp
    REAL(8), ALLOCATABLE :: tmpz(:,:,:)
    REAL(8), PARAMETER :: verylarge=1.d30
    REAL(8) :: scale

    mesh=SIZE(yy(2,:))
    !write (6,*) ' in cdfdol with mesh jj1,j22', mesh, jj1,jj2
    IF (SIZE(zz(2,2,:))/=mesh) THEN
       WRITE(6,*) 'cfdsol error - incompatible arrays', mesh,SIZE(zz(2,2,:))
       STOP
    ENDIF
    isgn = ( jj2 - jj1 ) / iabs( jj2 - jj1 )
    IF ( isgn .EQ. + 1 ) THEN
       IF ( jj1 .LE. 5 .OR. jj2 .GT. mesh ) THEN
          WRITE(6,10) isgn,jj1,jj2,mesh
          CALL EXIT(1)
       ENDIF
    ELSEIF ( isgn .EQ. - 1 ) THEN
       IF ( jj1 .GE. ( mesh - 4 ) .OR. jj2 .LT. 1 ) THEN
          WRITE(6,10) isgn,jj1,jj2,mesh
          CALL EXIT(1)
       ENDIF
    ELSE
       WRITE(6,10) isgn,jj1,jj2,mesh
    ENDIF

10  FORMAT(' ***error in subroutine difsol',/,&
         &' isgn =',i2,' jj1 =',i5,' jj2 =',i5,' mesh =',i5,&
         &' are not allowed')

    ALLOCATE(tmpz(2,2,mesh))
    tmpz=zz

    DO i=1,2
       DO j=1,2
          tmpz(i,j,:)=tmpz(i,j,:)*Grid%h
          if (Grid%TYPE==loggrid) tmpz(i,j,1:mesh)=tmpz(i,j,1:mesh)*Grid%drdu(1:mesh)
       ENDDO
    ENDDO

    abp(1) = 1901.d0 / 720.d0
    abp(2) = -1387.d0 / 360.d0
    abp(3) = 109.d0 / 30.d0
    abp(4) = -637.d0 / 360.d0
    abp(5) = 251.d0 / 720.d0
    amc(0) = 251.d0 / 720.d0
    amc(1) = 323.d0 / 360.d0
    amc(2) = -11.d0 / 30.d0
    amc(3) = 53.d0 / 360.d0
    amc(4) = -19.d0 / 720.d0

    DO j = 1,5
       ip = jj1 - isgn * j
       fa(j) = tmpz(1,1,ip) * yy(1,ip) + tmpz(1,2,ip) * yy(2,ip)
       fb(j) = tmpz(2,1,ip) * yy(1,ip) + tmpz(2,2,ip) * yy(2,ip)
    ENDDO

    DO j = jj1,jj2,isgn
       arp = yy(1,j-isgn)
       brp = yy(2,j-isgn)
       IF (ABS(arp)>verylarge.OR.brp>verylarge) THEN
          scale=1.d0/(ABS(arp)+ABS(brp))
          arp=arp*scale
          brp=brp*scale
          fa(:)=fa(:)*scale; fb(:)=fb(:)*scale
          yy=yy*scale
       ENDIF
       DO  i = 1,5
          arp = arp + DBLE(isgn) * abp(i) * fa(i)
          brp = brp + DBLE(isgn) * abp(i) * fb(i)
       ENDDO

       fa(0) = tmpz(1,1,j) * arp + tmpz(1,2,j) * brp
       fb(0) = tmpz(2,1,j) * arp + tmpz(2,2,j) * brp
       yy(1,j) = yy(1,j-isgn)
       yy(2,j) = yy(2,j-isgn)
       DO  i = 0,4,1
          yy(1,j) = yy(1,j) + DBLE(isgn) * amc(i) * fa(i)
          yy(2,j) = yy(2,j) + DBLE(isgn) * amc(i) * fb(i)
       ENDDO

       DO i = 5,2,-1
          fa(i) = fa(i-1)
          fb(i) = fb(i-1)
       ENDDO
       fa(1) = tmpz(1,1,j) * yy(1,j) + tmpz(1,2,j) * yy(2,j)
       fb(1) = tmpz(2,1,j) * yy(1,j) + tmpz(2,2,j) * yy(2,j)
    ENDDO

    DEALLOCATE(tmpz)
  END SUBROUTINE cfdsol


  !*********************************************************
  ! subroutine mod_backward_numerov(Grid,match,ww,wfn)
  !    version modified for scalar relativistic case when
  !     l and energy terms represented in array ww
  !*********************************************************
  SUBROUTINE mod_backward_numerov(Grid,match,ww,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: match
    REAL(8), INTENT(IN) :: ww(:)
    REAL(8), INTENT(INOUT) :: wfn(:)    ! on input wfn(n-1) and wfn(n) given
    ! on ouput wfn(i) given for
    !      start<=i<=n

    REAL(8), ALLOCATABLE :: a(:),b(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    REAL(8), PARAMETER :: vlarg=1.d30
    INTEGER :: i,j,k,n

    n=Grid%n
    ALLOCATE(a(n),b(n),p(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error backward_numerov ', n,i
       STOP
    ENDIF

    p(n)=wfn(n)
    p(n-1)=wfn(n-1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=match,n
       a(i)=ww(i)
    ENDDO
    IF (Grid%type==loggrid) THEN
       p(n-1:n)=p(n-1:n)/Grid%pref(n-1:n)
       a(match:n)=0.25d0+Grid%rr02(match:n)*a(match:n)
    ENDIF

    b(match:n)=2.4d0+h2*a(match:n)
    a(match:n)=1.2d0-0.1d0*h2*a(match:n)

    DO i=n-2,match,-1
       p(i)=(b(i+1)*p(i+1)-a(i+2)*p(i+2))/a(i)
       !renormalize if necessary
       scale=ABS(p(i))
       IF (scale > vlarg) THEN
          scale=1.d0/scale
          p(i:n)=scale*p(i:n)
       ENDIF
    ENDDO

    wfn(match:n)=p(match:n)

    IF (Grid%type==loggrid) THEN
       wfn(match:n)=wfn(match:n)*Grid%pref(match:n)
    ENDIF

    DEALLOCATE(a,b,p)

  END SUBROUTINE mod_backward_numerov
  !******************************************************************
  !  subroutine kinetic(Grid,wfn,l,ekin)
  !       calculates expectation value of kinetic energy for wfn
  !        with orbital angular momentum l
  !        wfn == r*radialwfn in Schroedinger Equation
  !        assumes wfn=(constant)*r^(l+1) at small r
  !*****************************************************************
  SUBROUTINE kinetic(Grid,wfn,l,ekin)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: wfn(:)
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(OUT) :: ekin

    REAL(8), ALLOCATABLE :: dfdr(:),arg(:)
    INTEGER :: i,n

    n=Grid%n

    ALLOCATE(dfdr(n),arg(n),stat=i)

    CALL derivative(Grid,wfn,dfdr)

    arg=0
    DO i=2,n
       arg(i)=wfn(i)/Grid%r(i)
    ENDDO

    DO i=1,n
       arg(i)=(dfdr(i))**2+(l*(l+1))*(arg(i))**2
    ENDDO

    ekin=integrator(Grid,arg)

    DEALLOCATE(dfdr,arg)
  END SUBROUTINE kinetic

  !******************************************************************
  !  subroutine kinetic_ij(Grid,wfn1,wfn2,l,ekin)
  !       calculates matrix element  of kinetic energy for wfn1 and wfn2
  !        with orbital angular momentum l
  !        wfn == r*radialwfn in Schroedinger Equation
  !        assumes wfn=(constant)*r^(l+1) at small r
  !*****************************************************************
  SUBROUTINE kinetic_ij(Grid,wfn1,wfn2,l,ekin,last)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: wfn1(:),wfn2(:)
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(OUT) :: ekin
    INTEGER, INTENT(IN), OPTIONAL :: last

    REAL(8), ALLOCATABLE :: dfdr1(:),dfdr2(:),arg1(:),arg2(:)
    INTEGER :: i,n

    n=Grid%n
    if (present(last)) n=last
    if (ABS(wfn1(n)*wfn2(n))> 1.d-8) then
       write(6,*) 'WARNING in kinetic_ij ',ABS(wfn1(n)*wfn2(n))
    endif
    ALLOCATE(dfdr1(n),arg1(n),dfdr2(n),arg2(n),stat=i)

    CALL derivative(Grid,wfn2,dfdr2,1,n)
    CALL derivative(Grid,wfn1,dfdr1,1,n)

    arg1=0; arg2=0
    DO i=2,n
       arg1(i)=wfn1(i)/Grid%r(i)
       arg2(i)=wfn2(i)/Grid%r(i)
    ENDDO

    DO i=1,n
       arg1(i)=(dfdr1(i)*dfdr2(i))+(l*(l+1))*(arg1(i)*arg2(i))
    ENDDO

    ekin=integrator(Grid,arg1,1,n)

    DEALLOCATE(dfdr1,arg1,dfdr2,arg2)
  END SUBROUTINE kinetic_ij

  !******************************************************************
  !  subroutine deltakinetic_ij(Grid,wfn1,wfn2,twfn1,twfn2,l,ekin,last)
  !       calculates difference matrix element  of kinetic energy for
  !         all electron functions wfn1 and wfn2
  !         and pseudo functions twfn1 and twfn2
  !        with orbital angular momentum l
  !        wfn == r*radialwfn in Schroedinger Equation
  !        assumes wfn=(constant)*r^(l+1) at small r
  !*****************************************************************
  SUBROUTINE deltakinetic_ij(Grid,wfn1,wfn2,twfn1,twfn2,l,ekin,last)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: wfn1(:),wfn2(:),twfn1(:),twfn2(:)
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(OUT) :: ekin
    INTEGER, INTENT(IN), OPTIONAL :: last

    REAL(8), ALLOCATABLE :: dfdr1(:),dfdr2(:),arg1(:),arg2(:)
    REAL(8), ALLOCATABLE :: tdfdr1(:),tdfdr2(:),targ1(:),targ2(:)
    INTEGER :: i,n

    n=Grid%n
    if (present(last)) n=last
    ALLOCATE(dfdr1(n),arg1(n),dfdr2(n),arg2(n),stat=i)
    ALLOCATE(tdfdr1(n),targ1(n),tdfdr2(n),targ2(n),stat=i)

    CALL derivative(Grid,wfn1,dfdr1,1,n)
    CALL derivative(Grid,wfn2,dfdr2,1,n)
    CALL derivative(Grid,twfn1,tdfdr1,1,n)
    CALL derivative(Grid,twfn2,tdfdr2,1,n)

    arg1=0; arg2=0; targ1=0; targ2=0
    DO i=2,n
       arg1(i)=wfn1(i)/Grid%r(i)
       arg2(i)=wfn2(i)/Grid%r(i)
       targ1(i)=twfn1(i)/Grid%r(i)
       targ2(i)=twfn2(i)/Grid%r(i)
    ENDDO

    DO i=1,n
       arg1(i)=(dfdr1(i)*dfdr2(i)-tdfdr1(i)*tdfdr2(i))&
&           +(l*(l+1))*(arg1(i)*arg2(i)-targ1(i)*targ2(i))
    ENDDO

    ekin=integrator(Grid,arg1,1,n)

    DEALLOCATE(dfdr1,arg1,dfdr2,arg2,tdfdr1,targ1,tdfdr2,targ2)
  END SUBROUTINE deltakinetic_ij
  !******************************************************************
  !******************************************************************
  !  subroutine altkinetic(Grid,wfn,energy,rv,ekin)
  !       calculates expectation value of kinetic energy for wfn
  !        with orbital wfn by integrating
  !          int(wfn**2 * (energy-rv/r), r=0..rmax)
  !*****************************************************************
  SUBROUTINE altkinetic(Grid,wfn,energy,rv,ekin)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: wfn(:),rv(:),energy
    REAL(8), INTENT(OUT) :: ekin

    REAL(8), ALLOCATABLE :: arg(:)
    INTEGER :: i,n

    n=Grid%n

    ALLOCATE(arg(n),stat=i)

    arg=0
    DO i=2,n
       arg(i)=(wfn(i)**2)*(energy-rv(i)/Grid%r(i))
    ENDDO

    ekin=integrator(Grid,arg)

    DEALLOCATE(arg)
  END SUBROUTINE altkinetic

  !****************************************************************
  ! function FindGridIndex(Grid,rpoint)
  !****************************************************************
  FUNCTION FindGridIndex(Grid,rpoint)
    INTEGER :: FindGridIndex
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: rpoint

    REAL(8) :: r0

    FindGridIndex=0
    IF (Grid%type==lineargrid) THEN
       FindGridIndex=rpoint/Grid%h+1
       IF (Grid%h*(FindGridIndex-1)<rpoint-1.d-10) FindGridIndex=FindGridIndex+1
    ELSEIF (Grid%type==loggrid) THEN
       r0=Grid%drdu(1)
       FindGridIndex=LOG(rpoint/r0+1)/Grid%h+1
       IF (r0*EXP(Grid%h*(FindGridIndex-1))<rpoint-1.d-10) &
&           FindGridIndex=FindGridIndex+1
    ENDIF
  END FUNCTION FindGridIndex

  !****************
  ! finite difference second derivative
  !   based on 5 point formula
  !   Ref. Engeln-Mullges & Uhlig (1996)
  !*****************
  FUNCTION secondderiv(index,f,h)
    REAL(8) :: secondderiv
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: f(:),h

    INTEGER :: i,n

    n=SIZE(f)
    secondderiv=0

    if (index==1.and.n>=5) THEN
       secondderiv=(70*f(1)-208*f(2)+228*f(3)-112*f(4)+22*f(5))/(24*h*h)
    else if (index==2.and.n>=5) THEN
       secondderiv=(22*f(1)-40*f(2)+12*f(3)+8*f(4)-2*f(5))/(24*h*h)
    else if (index>2.and.index<=n-2) THEN
       secondderiv=-(f(index-2)+f(index+2))/12 + &
&        4*(f(index-1)+f(index+1))/3 - 5*f(index)/2
       secondderiv=secondderiv/(h*h)
    else if (index>=5.and.index==n-1)   THEN
       secondderiv=(-2*f(n-4)+8*f(n-3)+12*f(n-2)-40*f(n-1)+22*f(n))/(24*h*h)
    else if (index>=5.and.index==n)   THEN
       secondderiv=(22*f(n-4)-112*f(n-3)+228*f(n-2)-208*f(n-1)+70*f(n))/(24*h*h)
    else
       WRITE(6,*) 'Error in secondderiv', index, SIZE(f)
       STOP
    ENDIF

  END FUNCTION secondderiv


  !****************
  ! finite difference first derivative
  !*****************
  FUNCTION firstderiv(index,f,h)
    REAL(8) :: firstderiv
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: f(:),h

    INTEGER :: n

    n=SIZE(f)
    firstderiv=0

    if (index==1.and.n>=5) THEN
       firstderiv=(-25*f(1)+48*f(2)-36*f(3)+16*f(4)-3*f(5))/(12*h)
    else if (index==2.and.n>=5) THEN
       firstderiv=(-3*f(1)-10*f(2)+18*f(3)-6*f(4)+f(5))/(12*h)
    else if (index>2.and.index<=n-2) THEN
       firstderiv=(f(index-2)-8*f(index-1)+8*f(index+1)-f(index+2))/(12*h)
    else if (index>=5.and.index==n-1)   THEN
       firstderiv=(-f(n-4)+6*f(n-3)-18*f(n-2)+10*f(n-1)+3*f(n))/(12*h)
    else if (index>=5.and.index==n)   THEN
       firstderiv=(3*f(n-4)-16*f(n-3)+36*f(n-2)-48*f(n-1)+25*f(n))/(12*h)
    else
       WRITE(6,*) 'Error in firstderiv', index, SIZE(f)
       STOP
    ENDIF
  END FUNCTION firstderiv


  !*****************
  ! Second derivative for general grid
  !*****************

  FUNCTION Gsecondderiv(Grid,index,g)
    REAL(8) :: Gsecondderiv
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: g(:)

    REAL(8), ALLOCATABLE :: f(:)
    INTEGER :: i,n

    Gsecondderiv=0

    IF (Grid%type==lineargrid) THEN
       Gsecondderiv=secondderiv(index,g,Grid%h)
    ELSEIF  (Grid%type==loggrid) THEN
       Gsecondderiv=(secondderiv(index,g,Grid%h)&
&            -firstderiv(index,g,Grid%h))/Grid%rr02(index)
    ENDIF

  END FUNCTION Gsecondderiv


  !*****************
  ! First  derivative for general grid
  !*****************

  FUNCTION Gfirstderiv(Grid,index,g)
    REAL(8) :: Gfirstderiv
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: index
    REAL(8), INTENT(IN) :: g(:)

    REAL(8), ALLOCATABLE :: f(:)
    INTEGER :: i,n

    Gfirstderiv=0

    IF (Grid%type==lineargrid) THEN
       Gfirstderiv=firstderiv(index,g,Grid%h)
    ELSEIF  (Grid%type==loggrid) THEN
       Gfirstderiv=firstderiv(index,g,Grid%h)/Grid%drdu(index)
    ENDIF

  END FUNCTION Gfirstderiv


  !*****************************************************************
  ! subroutine reportgrid(Grid,unit)
  !*****************************************************************
  SUBROUTINE reportgrid(Grid,unit)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: unit

    IF (Grid%type==lineargrid) THEN
       WRITE(unit,*) ' Radial integration grid is linear '
       WRITE(unit,'(" h = ", 1p,1e15.7,"   n = ",i9," rmax = ", 1p,1e15.7)')&
&       Grid%h,Grid%n,Grid%r(Grid%n)
    ELSEIF (Grid%type==loggrid) THEN
       WRITE(unit,*) ' Radial integration grid is logarithmic '
       WRITE(unit,&
&       '("r0 = ",1p,1E15.7," h = ", 1p,1e15.7,"   n = ",i9," rmax = ", 1p,1e15.7)')&
&       Grid%drdu(1), Grid%h, Grid%n,Grid%r(Grid%n)
    ENDIF
  END SUBROUTINE reportgrid


  !******************************************************************
  ! function gridindex(Grid,r)
  !******************************************************************
  FUNCTION gridindex(Grid,r)
    INTEGER :: gridindex
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: r

    gridindex=0
    IF (Grid%type==lineargrid) THEN
       gridindex=r/Grid%h +0.1d0 +1
    ELSEIF (Grid%type==loggrid) THEN
       gridindex=LOG(1.d0+r/Grid%drdu(1))/Grid%h +0.1d0 +1
    ENDIF

  END FUNCTION gridindex

  !*********************************************************************
  !  subroutine findh(Z,range,n,hval,r0)
  !    find hval for fixed number of input grid points n in loggrid case
  !    assumes form r(i)=(h/Z)*(exp(h*(i-1))-1);   r0=h/Z
  !*********************************************************************
   SUBROUTINE findh(Z,range,n,hval,r0)
     REAL(8), INTENT(IN) :: Z
     INTEGER, INTENT(IN) :: n
     REAL(8), INTENT(IN) :: range
     REAL(8), INTENT(INOUT) :: hval,r0

     REAL(8) :: h0,dh,f,df
     INTEGER :: i,j,k
     INTEGER, parameter :: iter=1000
     REAL(8), parameter :: eps=1.e-15
     LOGICAL :: success

     h0=hval
     success=.false.
     do i=1,iter
        f=LOG(Z*range/h0+1.d0)/h0
        df=-f/h0-(Z*range/h0**3)/(Z*range/h0+1.d0)
        dh=(n-1-f)/df
        if (ABS(dh)< eps) then
           success=.true.
           exit
        endif
        if (h0+dh<0.d0) then
           h0=h0/2
        else
           h0=h0+dh
        endif
      enddo

      if (.not.success) then
        write(6,*) 'Warning in findh -- dh > eps ', dh,h0
      endif
      hval=h0
      r0=hval/Z

  end subroutine findh


  !*********************************************************************
  !  subroutine findh_given_r0(Z,range,r0,n,hval)
  !    find hval for fixed number of input grid points n in loggrid case
  !    assumes form r(i)=(r0/Z)*(exp(h*(i-1))-1);
  !*********************************************************************
   SUBROUTINE findh_given_r0(Z,range,r0,n,hval)
     INTEGER, INTENT(IN) :: n
     REAL(8), INTENT(IN) :: Z,range,r0
     REAL(8), INTENT(INOUT) :: hval

     REAL(8) :: h0,dh,f,df
     INTEGER :: i,j,k

     hval=log((Z*range/r0 + 1.d0))/(n-1)

  end subroutine findh_given_r0

  !*********************************************************************
  !  subroutine findh_worse(Z,range,n,hval,r0)
  !    find hval for fixed number of input grid points n in loggrid case
  !    assumes form r(i)=(h/(Z**1/3))*(exp(h*(i-1))-1); r0=(h/(Z**1/3)
  !       Note: this choice of r0 does poor job for some integrals
  !                    than r0=h/Z
  !*********************************************************************
  SUBROUTINE findh_worse(Z,range,n,hval,r0)
    INTEGER, INTENT(IN) :: Z,n
    REAL(8), INTENT(IN) :: range
    REAL(8), INTENT(INOUT) :: hval,r0

    REAL(8) :: h0,dh,f,df,zz
    INTEGER :: i,j,k
    INTEGER, PARAMETER :: iter=1000
    REAL(8), PARAMETER :: eps=1.e-15
    LOGICAL :: success

    h0=hval;zz=z; zz=zz**(1.d0/3.d0)
    success=.FALSE.
    DO i=1,iter
       f=LOG(zz*range/h0+1.d0)/h0
       df=-f/h0-(zz*range/h0**3)/(zz*range/h0+1.d0)
       dh=(n-1-f)/df
       IF (ABS(dh)< eps) THEN
          success=.TRUE.
          EXIT
       ENDIF
       IF (h0+dh<0.d0) THEN
          h0=h0/2
       ELSE
          h0=h0+dh
       ENDIF
    ENDDO

    IF (.NOT.success) THEN
       WRITE(6,*) 'Warning in findh -- dh > eps ', dh,h0
    ENDIF
    hval=h0; r0=h0/zz

  END SUBROUTINE findh_worse



  !******************************************************************
  ! subroutine initgrid(Grid,h,range,r0)
  !******************************************************************
  SUBROUTINE InitGrid(Grid,h,range,r0)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    REAL(8), INTENT(IN) :: range
    REAL(8), INTENT(IN) :: h
    REAL(8), OPTIONAL, INTENT(IN) :: r0

    INTEGER :: i,n

    IF (PRESENT(r0)) THEN

       Grid%type=loggrid
       n=LOG(range/r0+1)/h+1
       Grid%ishift=5
       IF (r0*(EXP(h*(n-1))-1)<range-1.d-5) n=n+1
       Grid%h=h
       Grid%n=n
       WRITE(6,*) 'InitGrid: -- logarithmic ',n, h,range,r0
       ALLOCATE(Grid%r(n),Grid%drdu(n),Grid%pref(n),Grid%rr02(n),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Allocation error in initgrid ', n,i
          STOP
       ENDIF
       DO i=1,n
          Grid%r(i)=r0*(EXP(Grid%h*(i-1))-1)
          Grid%drdu(i)=r0*EXP(Grid%h*(i-1))
          Grid%pref(i)=r0*EXP(Grid%h*(i-1)/2.d0)
          Grid%rr02(i)=(Grid%r(i)+r0)**2
       ENDDO
    ELSE
       Grid%type=lineargrid
       n=range/h+1
       Grid%ishift=25
       IF (h*(n-1)<range-1.d-5) n=n+1
       Grid%n=n
       Grid%h=h
       WRITE(6,*) 'InitGrid: -- linear  ', n,h,range
       ALLOCATE(Grid%r(n),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Allocation error in initgrid ', n,i
          STOP
       ENDIF
       DO i=1,n
          Grid%r(i)=(Grid%h*(i-1))
          Grid%drdu(i)=1.d0
       ENDDO
       NULLIFY(Grid%pref)
       NULLIFY(Grid%rr02)
    ENDIF

  END SUBROUTINE InitGrid

  !******************************************************************
  ! subroutine destroygrid(Grid)
  !******************************************************************
  SUBROUTINE DestroyGrid(Grid)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    IF (ASSOCIATED(Grid%r)) DEALLOCATE(Grid%r)
    IF (ASSOCIATED(Grid%drdu)) DEALLOCATE(Grid%drdu)
    IF (ASSOCIATED(Grid%pref)) DEALLOCATE(Grid%pref)
    IF (ASSOCIATED(Grid%rr02)) DEALLOCATE(Grid%rr02)
  END SUBROUTINE DestroyGrid

  !******************************************************************
  ! subroutine nullifygrid(Grid)
  !******************************************************************
  SUBROUTINE NullifyGrid(Grid)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    NULLIFY(Grid%r)
    NULLIFY(Grid%drdu)
    NULLIFY(Grid%pref)
    NULLIFY(Grid%rr02)
  END SUBROUTINE NullifyGrid

  SUBROUTINE ClassicalTurningPoint(Grid,rv,l,energy,turningpoint)
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
       WRITE(6,*) 'Allocation error in ClassicalTurningPoint ', i,n
       STOP
    ENDIF

    v=0
    v(2:n)=rv(2:n)/Grid%r(2:n)+l*(l+1)/(Grid%r(2:n)**2)

    turningpoint=n
    DO i=n,2,-1
       IF (v(i)<energy) EXIT
    ENDDO
    turningpoint=i
    turningpoint=MIN(turningpoint,FindGridIndex(Grid,10.0d0))

    !write(6,*) 'Found turning point at ', turningpoint, Grid%r(turningpoint)

    DEALLOCATE(v)

  END SUBROUTINE ClassicalTurningPoint

   INTEGER function countnodes(start,finish,wfn,filter)
      INTEGER, INTENT(IN) :: start,finish
      REAL(8), INTENT(IN) :: wfn(:)
      REAL(8), INTENT(IN), OPTIONAL :: filter

      INTEGER :: i,nodes
      nodes=0

      do i=start+1,finish
         if (wfn(i)*wfn(i-1)<0.d0) nodes=nodes+1
         if (PRESENT(filter).and.(abs(wfn(i))+abs(wfn(i+1)))<filter) exit
      enddo

      countnodes=nodes
   end function countnodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  On input:  many functions stored in wfn are orthonormal
  !     On output, if f is linearly independent, it is orthnormalized and
  !       stored as wfn(:,many+1) many -> many+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE gramschmidt(Grid,many,wfn,f)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: many
    REAL(8), INTENT(INOUT) :: wfn(:,:),f(:)

    REAL(8) :: h,x,norm
    INTEGER :: i,j,k,l,n,irc

    k=Size(f)

    norm=overlap(Grid,f,f,1,k)
    IF (many>0) THEN
       DO i=1,many
          x=overlap(Grid,f,wfn(:,i),1,k)
          f=f-x*wfn(:,i)
       ENDDO
       norm=overlap(Grid,f,f,1,k)
    ENDIF
    IF(norm>machine_zero) THEN
       f=f/SQRT(norm)
       wfn(:,many+1)=f
    ELSE
       write(6,*) '!!!! warning in Gram Schmidt !!!!!!', norm,many
       stop
    ENDIF

  END SUBROUTINE gramschmidt


    SUBROUTINE Milne(Grid,many,energy,l,ZZ,ve,x,wfn)
      TYPE (GridInfo), INTENT(IN) :: Grid
      INTEGER, INTENT(IN) :: l,many,ZZ
      REAL(8), INTENT(IN) :: energy,ve(:),x(:)   ! ve==Hartree, x=Fock
      REAL(8), INTENT(INOUT) :: wfn(:)

      REAL(8), allocatable :: z(:),y(:),zp(:),yp(:),A(:),B(:),C(:)
      REAL(8) :: c0,c1,c2,predy,corry,predz,corrz,err,diff
      INTEGER :: i,j,k,m,n,it,last
      INTEGER, parameter :: maxit=10

      diff=machine_precision*1000
      n=Grid%n
      if (many.gt.n.or.many.lt.6) then
            write(6,*) 'Error in Milne == n>many',n,many
            stop
      endif

      wfn=0
      allocate(z(many),y(many),zp(many),yp(many),A(many),B(many),C(many))

   ! Note that must not evaluate effective potential terms at r=0
      z=0;zp=0;y=0;yp=0;A=0;B=0;C=0
      c0=1.d0
      c1=-float(ZZ)/(l+1)
      C(2:many)=-x(2:many)/(Grid%r(2:many)**(l+1))
      CALL extrapolate(Grid,C)
      c2=(C(1)+ve(1)-energy + 2*float(ZZ**2)/(l+1))/(4*l+6)
      write(6,*) 'coeff', c0,c1,c2
      do i=1,6
        y(i)=c0+Grid%r(i)*(c1+Grid%r(i)*c2)
        yp(i)=c1+2*c2*Grid%r(i)
        z(i)=yp(i)
        zp(i)=2*c2
        write(6,'("init",1p,5e16.7)') Grid%r(i),y(i),yp(i),z(i),zp(i)
      enddo

      write(6,*) 'Starting iterations with h = ', Grid%h
      If(Grid%type==loggrid) then
        zp(1:6)= Grid%drdu(1:6)*(yp(1:6)+Grid%drdu(1:6)*2*c2)
        yp(1:6)= Grid%drdu(1:6)*yp(1:6)
        z(1:6)=yp(1:6)
        A(2:many)=2*(l+1)*Grid%drdu(2:many)/Grid%r(2:many)-1.d0
        B(2:many)=(energy+2*ZZ/Grid%r(2:many)-ve(2:many))*Grid%rr02(2:many)
        C(1:many)=C(1:many)*Grid%rr02(1:many)
      Else
        A(2:many)=2*(l+1)/Grid%r(2:many)
        B(2:many)=(energy+2*ZZ/Grid%r(2:many)-ve(2:many))
      ENDIF

      Do i=6,many-1
            predy=y(i-5)+0.3d0*Grid%h*(26.d0*yp(i-2) &
&               + 11.d0*(yp(i)+yp(i-4))-14.d0*(yp(i-1)+yp(i-3)))
            y(i+1)=predy
            predz=z(i-5)+0.3d0*Grid%h*(26.d0*zp(i-2) &
&               + 11.d0*(zp(i)+zp(i-4))-14.d0*(zp(i-1)+zp(i-3)))
            z(i+1)=predz
            yp(i+1)=z(i+1)
            zp(i+1)=C(i+1)-A(i+1)*z(i+1)-B(i+1)*y(i+1)
         Do it=1,maxit
            last=it
            corry=y(i-3)+0.04444444444444444d0*Grid%h*(12*yp(i-1) &
&               +7.d0*(yp(i+1)+yp(i-3))+32.d0*(yp(i)+yp(i-2)))
            corrz=z(i-3)+0.04444444444444444d0*Grid%h*(12*zp(i-1) &
&               +7.d0*(zp(i+1)+zp(i-3))+32.d0*(zp(i)+zp(i-2)))
            !write(800,'(i5,1p,8e16.7)') it,Grid%r(i), y(i+1),corry,z(i+1),corrz
            err= abs(corry-y(i+1))+abs(corrz-z(i+1))
            if (err<diff) exit
            z(i+1)=corrz
            y(i+1)=corry
            yp(i+1)=z(i+1)
            zp(i+1)=C(i+1)-A(i+1)*z(i+1)-B(i+1)*y(i+1)
        Enddo
        !Write(6,*) 'Completed PC ', i,last
        If (last==maxit) write(6,*) 'Warning from Pred-Corr',i,err
    Enddo
    Do i=1,many
       wfn(i)=(Grid%r(i)**(l+1))*y(i)
    Enddo
       wfn(1:many)=wfn(1:many)/wfn(many)

      deallocate(z,y,zp,yp,A,B,C)
   END SUBROUTINE

  !**********************************************************************
  !  pgm to integrate outward the radial schroedinger inhomogeneous equation
  !    at energy 'energy' and at angular momentum l
  !    with potential smooth rv/r   and inhomogeous contribution rhs
  !     Assumes wfn(i) is known for 1 <= i <= istart
  !  uses Noumerov algorithm
  !*************************************************************************
  !*************************************************************
  ! subroutine midrange_numerov(Grid,l,istart,many,energy,rv,rhs,wfn)
  !*************************************************************
  SUBROUTINE midrange_numerov(Grid,l,istart,many,energy,rv,rhs,wfn)
    TYPE (GridInfo), INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,istart,many
    REAL(8), INTENT(IN) :: energy,rv(:),rhs(:)
    REAL(8), INTENT(INOUT) :: wfn(:)
    ! initial values of wfn determined from proj(r=0)
    REAL(8), ALLOCATABLE :: a(:),b(:),c(:),p(:)
    REAL(8) :: xx,angm,h,h2,scale
    INTEGER :: i,j,k,n
    INTEGER :: count=0

    ALLOCATE(a(many),b(many),c(many),p(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error midrange_numerov ', many,i
       STOP
    ENDIF

     !Note   values of a,b,c for i<istart are not used

     if (istart<3) then
         write(6,*) 'Error in midrange_numerov -- istart < 3 ', istart
         stop
     elseif (istart>many-1) then
         write(6,*) 'Error in midrange_numerov -- istart > many ', istart
         stop
     endif


    p(1:istart)=wfn(1:istart)
    angm=l*(l+1)
    a=0
    h=Grid%h;    h2=h*h
    DO i=2,many
       a(i)=rv(i)/Grid%r(i)-energy+angm/(Grid%r(i)**2)
    ENDDO
    b(1:many)=-0.1d0*h2*rhs(1:many)
    IF (Grid%type==loggrid) THEN
       p(1:istart)=wfn(1:istart)/Grid%pref(1:istart)
       a=0.25d0+Grid%rr02(1:many)*a
       b=Grid%rr02(1:many)*b/Grid%pref(1:many)
    ENDIF

    c=0
    DO i=2,many-1
       c(i)=10*b(i)+b(i-1)+b(i+1)
    ENDDO

    b=2.4d0+h2*a
    a=1.2d0-0.1d0*h2*a

    DO i=istart+1,many
       p(i)=(b(i-1)*p(i-1)-a(i-2)*p(i-2)-c(i-1))/a(i)
    ENDDO

    wfn(1:many)=p(1:many)

    IF (Grid%type==loggrid) THEN
       wfn(1:many)=wfn(1:many)*Grid%pref(1:many)
    ENDIF

    DEALLOCATE(a,b,c,p)

  END SUBROUTINE midrange_numerov

END MODULE gridmod
