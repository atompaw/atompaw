MODULE GlobalMath

  IMPLICIT NONE

  REAL(8) :: pi , machine_precision , machine_zero , machine_infinity

  REAL(8), PRIVATE :: minlog,maxlog,minexp,maxexp
  REAL(8), PRIVATE :: minlogarg,maxlogarg,minexparg,maxexparg

CONTAINS

  !****************************************************
  SUBROUTINE Init_GlobalConstants()

    INTEGER :: i
    REAL(8)    :: tmp , a1,a2,a3

    ! Calculate machine accuracy
    machine_precision = 0
    a1 = 4.d0/3.d0
    DO WHILE (machine_precision == 0.d0)
       a2 = a1 - 1.d0
       a3 = a2 + a2 + a2
       machine_precision = ABS(a3 - 1.d0)
    ENDDO

    machine_zero= machine_precision**4
    machine_infinity = 1.d0/machine_zero

    pi = ACOS(-1.d0)

    minlogarg=machine_precision; minlog=LOG(minlogarg)
    maxlogarg=1.d0/machine_precision; maxlog=LOG(maxlogarg)
    minexparg=LOG(machine_precision);  minexp=0.d0
    maxexparg=-LOG(machine_precision);  maxexp=EXP(maxexparg)

    RETURN
  END SUBROUTINE Init_GlobalConstants

  !****************************************************
  FUNCTION ddlog(arg)
    REAL(8) :: arg,ddlog

    IF (arg>maxlogarg) THEN
       ddlog=maxlog
    ELSE IF (arg<minlogarg) THEN
       ddlog=minlog
    ELSE
       ddlog=LOG(arg)
    ENDIF

    RETURN
  END FUNCTION ddlog

  !****************************************************
  FUNCTION ddexp(arg)
    REAL(8) :: arg,ddexp

    IF (arg>maxexparg) THEN
       ddexp=maxexp
    ELSE IF (arg<minexparg) THEN
       ddexp=minexp
    ELSE
       ddexp=EXP(arg)
    ENDIF

    RETURN
  END FUNCTION ddexp

  !**********************************************
  SUBROUTINE sphbes(lval,ndim,x)
    INTEGER, INTENT(IN)    :: lval
    INTEGER, INTENT(IN)    :: ndim
    REAL(8),    INTENT(INOUT) :: x(ndim)

    INTEGER :: it,lll,l,n,il
    REAL(8) :: z,bs,bs0,bs1,xx,z2,arg,sum
    REAL(8), PARAMETER :: tol=1.e-11

    DO it=1,ndim
       z=x(it)

       IF(z.LE.0.5) THEN ! series expansion for small arguments
          IF (z.LE.tol) THEN ! - special treatment if z = 0
             bs=0.d0
             IF(lval.EQ.0) bs=1.d0

          ELSE
             lll=2*lval+1
             xx=1.d0
             DO l=1,lll,2
                xx=xx*l
             END DO

             z2=0.5d0*z*z
             lll=lll+2
             arg=-z2/lll
             sum=1.d0+arg

             DO  n=2,500
                lll=lll+2
                arg=-arg*z2/(n*lll)
                sum=sum+arg
                IF (ABS(arg).LE.tol) EXIT
             END DO

             bs=(z**lval)*sum/xx
          END IF
       ELSE        ! -- trigometric form
          bs=SIN(z)/z

          IF (lval.GT.0) THEN
             bs0=bs
             bs=(bs0-COS(z))/z
             IF (lval.GT.1) THEN
                bs1=bs

                DO  il=2,lval
                   bs=(2*il-1)*bs1/z-bs0
                   bs0=bs1
                   bs1=bs
                END DO
             END IF
          END IF
       END IF

       IF (ABS(bs)<1D-50) bs = 0
       x(it)=bs

    END DO

    RETURN
  END  SUBROUTINE sphbes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE bessandneumfunc(l,val,jl,jlp,nl,nlp)
!      returns spherical bessel and neumann functions and
!        their derivatives;  actually jl== j_l(val)*val
!                                     nl== n_l(val)*val
!             jlp and nlp are their derivatives;
!          programmed so far only for l<=2
    SUBROUTINE bessandneumfunc(l,val,jl,jlp,nl,nlp)
       INTEGER, INTENT(IN) :: l
       REAL(8), INTENT(IN) :: val
       REAL(8), INTENT(OUT) :: jl,jlp,nl,nlp

       REAL(8) :: x

       if (val<1.d-8) then
          write(6,*) 'error in bessandneumfunc   -- argument should  be > 0',val
          stop
       endif

       if (l==0) then
          jl=SIN(val) ; nl= -COS(val)
          jlp=COS(val) ; nlp= SIN(val)
       elseif (l==1) then
          jl=SIN(val)/val-COS(val);   nl=-COS(val)/val-SIN(val)
          jlp=SIN(val)*(1-1.d0/val**2) + COS(val)/val
          nlp=-COS(val)*(1-1.d0/val**2) + SIN(val)/val
       elseif (l==2) then
          jl=-SIN(val)*(1-3.d0/val**2)-3*COS(val)/val
          nl=COS(val)*(1-3.d0/val**2)-3*SIN(val)/val
          jlp=3*SIN(val)*(1-2.d0/val**2)/val-COS(val)*(1-6.d0/val**2)
          nlp=-3*COS(val)*(1-2.d0/val**2)/val-SIN(val)*(1-6.d0/val**2)
       else
          write(6,*) 'Error in bessandneumfunc ;  '&
&             ,'Only l<=2 has been programmed so far '
          stop
       endif

       write(6,*) 'Check Bessel Wronskian : ',l,jl*nlp-nl*jlp

    END SUBROUTINE bessandneumfunc
  !**********************************************************************

  FUNCTION ranx()
    REAL(8) :: ranx
    INTEGER, PARAMETER :: konst=125
    INTEGER  :: m=100001
    m=m*konst
    m=m-2796203*(m/2796203)
    ranx=m/2796203.d0
    RETURN
  END FUNCTION ranx

  SUBROUTINE shift4(v1,v2,v3,v4,NEW)
    REAL(8), INTENT(IN) :: NEW
    REAL(8), INTENT(INOUT) :: v1,v2,v3,v4

    v1=v2
    v2=v3
    v3=v4
    v4=NEW

    RETURN
  END SUBROUTINE shift4

  SUBROUTINE mkname(i,stuff)
    !  subroutine to take an integer (.le. 4 digits) and return it
    !   in the form of a character string
    CHARACTER(4) stuff
    INTEGER i,i1,i10,i100,i1000
    stuff='?'
    IF (i.GT.9999) i=MOD(i,10000)
    i1000=i/1000
    i100=(i-1000*i1000)/100
    i10=(i-1000*i1000-100*i100)/10
    i1=(i-1000*i1000-100*i100-10*i10)
    IF (i.GE.1000) THEN
       stuff=CHAR(i1000+48)//CHAR(i100+48)//CHAR(i10+48)//CHAR(i1+48)
       RETURN
    ENDIF
    IF (i.GE.100) THEN
       stuff=CHAR(i100+48)//CHAR(i10+48)//CHAR(i1+48)
       RETURN
    ENDIF
    IF (i.GE.10) THEN
       stuff=CHAR(i10+48)//CHAR(i1+48)
       RETURN
    ENDIF
    IF (i.GE.0) stuff=CHAR(i1+48)
    RETURN
  END SUBROUTINE mkname

  SUBROUTINE linsol(a,b,kk,la,ra,lb,det)
    INTEGER, INTENT(IN) :: kk,la,ra,lb
    REAL(8), INTENT(INOUT) :: a(la,ra),b(lb)
    REAL(8), OPTIONAL, INTENT(OUT) :: det

    REAL(8) :: d,s,r
    INTEGER :: kkm,i,j,k,l,ipo,n,kmo

    d = 1.d0
    if (kk>min(la,ra,lb)) then
       write(6,*) 'Dimension error in linsol ', la,ra,lb,kk
       stop
    endif
    kkm=kk-1
    IF (kkm == 0) THEN
       b(1)=b(1)/a(1,1)
    ELSE IF (kkm > 0) THEN
       DO i=1, kkm
          s = 0.0d0
          l=i
          DO j=i,kk
             r=ABS(a(j,i))
             IF(r >  s) THEN
                s=r
                l=j
             ENDIF
          ENDDO
          IF(l /= i) THEN
             DO j=i,kk
                s=a(i,j)
                a(i,j)=a(l,j)
                a(l,j)=s
             ENDDO
             s=b(i)
             b(i)=b(l)
             b(l)=s
             d = -d
          ENDIF
          IF (a(i,i) /= 0.0d0) THEN
             ipo=i+1
             DO j=ipo,kk
                IF (a(j,i) /= 0.0d0) THEN
                   s=a(j,i)/a(i,i)
                   a(j,i) = 0.0d0
                   DO k=ipo,kk
                      a(j,k)=a(j,k)-a(i,k)*s
                   ENDDO
                   b(j)=b(j)-b(i)*s
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       DO i=1,kk
          d=d*a(i,i)
       ENDDO
       kmo=kk-1
       b(kk)=b(kk)/a(kk,kk)
       DO i=1,kmo
          n=kk-i
          DO j=n,kmo
             b(n)=b(n)-a(n,j+1)*b(j+1)
          ENDDO
          b(n)=b(n)/a(n,n)
       ENDDO
    ENDIF
    !write(6,*) 'determinant from linsol ' , d
    IF(ABS(d).LT.1.d-10) WRITE(6,*) '**warning from linsol --',&
&        'determinant too small --',d
    If (present(det)) det=d
  END SUBROUTINE linsol

  SUBROUTINE minverse(a,kk,la,ra)
    INTEGER, INTENT(IN) :: kk,la,ra
    REAL(8), INTENT(INOUT) :: a(la,ra)

    REAL(8) :: d,s,r
    REAL(8), allocatable :: ai(:,:)
    INTEGER :: kkm,i,j,k,l,ipo,n,kmo

    if (kk>min(la,ra)) then
       write(6,*) 'Dimension error in minverse ', la,ra,kk
       stop
    endif
    ALLOCATE(ai(kk,kk),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Allocation error in minverse',i,kk
       STOP
    ENDIF

    ai=0
    DO i=1,kk
       ai(i,i)=1
    ENDDO

    d = 1.0d0
    kkm=kk-1
    IF (kkm == 0) THEN
       ai(1,1)=ai(1,1)/a(1,1)
    ELSE IF (kkm > 0) THEN
       DO i=1, kkm
          s = 0.0d0
          DO j=i,kk
             r=ABS(a(j,i))
             IF(r >  s) THEN
                s=r
                l=j
             ENDIF
          ENDDO
          IF(l /= i) THEN
             DO j=i,kk
                s=a(i,j)
                a(i,j)=a(l,j)
                a(l,j)=s
             ENDDO
             DO k=1,kk
                s=ai(i,k)
                ai(i,k)=ai(l,k)
                ai(l,k)=s
             ENDDO
             d = -d
          ENDIF
          IF (a(i,i) /= 0.0d0) THEN
             ipo=i+1
             DO j=ipo,kk
                IF (a(j,i) /= 0.0d0) THEN
                   s=a(j,i)/a(i,i)
                   a(j,i) = 0.0d0
                   DO k=ipo,kk
                      a(j,k)=a(j,k)-a(i,k)*s
                   ENDDO
                   DO k=1,kk
                      ai(j,k)=ai(j,k)-ai(i,k)*s
                   ENDDO
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       DO i=1,kk
          d=d*a(i,i)
       ENDDO
       kmo=kk-1
       DO k=1,kk
          ai(kk,k)=ai(kk,k)/a(kk,kk)
          DO i=1,kmo
             n=kk-i
             DO j=n,kmo
                ai(n,k)=ai(n,k)-a(n,j+1)*ai(j+1,k)
             ENDDO
             ai(n,k)=ai(n,k)/a(n,n)
          ENDDO
       ENDDO
    ENDIF
    a=0.d0
    do i=1,kk
       do j=1,kk
          a(i,j)=ai(i,j)
       enddo
    enddo
    DEALLOCATE(ai)

    IF(ABS(d).LT.1.d-10) WRITE(6,*) '**warning from linsol --',&
&        'determinant too small --',d
  END SUBROUTINE minverse

  FUNCTION factorial(n)
    REAL(8) :: factorial
    INTEGER, INTENT(IN) :: n
    INTEGER :: i
    factorial=1
    IF (n.LT.2) RETURN
    DO i=2,n
       factorial=factorial*i
    ENDDO
  END FUNCTION factorial

  !*************************************************************************
  ! FUNCTION hwfn(z,np,l,r)
  !
  ! function to calculate the radial H wfn for nuclear charge z
  !          (note in this version z is real and need not be integral)
  !                                            principal qn   np
  !                                            orbital qn     l
  !   r*(radial H wfn) is returned
  !*************************************************************************
  FUNCTION hwfn(z,np,l,r)
    REAL(8) :: hwfn
    REAL(8), INTENT(IN) :: z,r
    INTEGER, INTENT(IN) :: np,l
    INTEGER :: node,k
    REAL(8) :: scale,rho,pref,term,sum
    node=np-l-1
    scale=2.d0*z/np
    rho=scale*r
    pref=scale*SQRT(scale*factorial(np+l)/(2*np*factorial(node)))
    term=(rho**l)/factorial(2*l+1)
    sum=term
    IF (node.GT.0) THEN
       DO k=1,node
          term=-term*(node-k+1)*rho/(k*(2*l+1+k))
          sum=sum+term
       ENDDO
    ENDIF
    hwfn=r*pref*EXP(-0.5d0*rho)*sum
    !     write(6,*) 'r,hwfn=',r,hwfn
  END FUNCTION hwfn

  !***************************************************************************
  ! subroutine filter(n,func,small)
  !***************************************************************************
  SUBROUTINE filter(n,func,small)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: func(:)
    REAL(8), INTENT(IN) :: small

    INTEGER :: i
    DO i=1,n
       IF (ABS(func(i)).LT.small) func(i)=0.d0
    ENDDO
  END SUBROUTINE filter

  !******************************************************************************
  !  From:
  !  File : misc.f90
  !    by : Alan Tackett
  !    on : 10/17/95
  !   for : Misc general purpose functions
  !
  !  This module contains misc. general purpose routines that don't fit
  !  in another library. Below is a list of routines contained :
  !
  !  PrintDate(Unit, Text)
  !      Prints the date and time to the specified unit along with the TEXT
  !
  !*****************************************************************************

  !******************************************************************************
  !
  !  PrintDate - Prints the date to the specified unit with TEXT prepended.
  !
  !    Unit - Output unit
  !    Text - Text for prepending. All trailing blanks are removed.
  !
  !******************************************************************************

  SUBROUTINE PrintDate(Unit, Text)
    INTEGER,       INTENT(IN) :: Unit
    CHARACTER*(*), INTENT(IN) :: Text

    CHARACTER*10 :: DateStr, TimeStr
    CHARACTER*50 :: FmtStr

    CALL DATE_AND_TIME(DateStr, TimeStr)

    FmtStr=' ' // DateStr(5:6) // '/' // DateStr(7:8) // '/' // DateStr(1:4) // &
&        ', ' // TimeStr(1:2) // ':' // TimeStr(3:4) // ':' // &
&        TimeStr(5:10)

    WRITE(Unit, '(a,a)')  TRIM(Text), TRIM(FmtStr)

    RETURN
  END SUBROUTINE PrintDate

!******************************************************************************
!
!  PrintDateStr - Prints the date to the specified string
!
!******************************************************************************

Subroutine PrintDateStr(text)
  Character*(*), Intent(out) :: Text

  Character*10 :: DateStr, TimeStr
  Character*50 :: FmtStr

  Call Date_And_Time(DateStr, TimeStr)

  FmtStr=' ' // DateStr(5:6) // '/' // DateStr(7:8) // '/' // DateStr(1:4) // &
       &       ', ' // TimeStr(1:2) // ':' // TimeStr(3:4) // ':' // &
       &         TimeStr(5:10)

  Write(text, '(a)')  Trim(FmtStr)

  Return
End Subroutine


  FUNCTION stripchar(inputchar)
    CHARACTER(132) :: stripchar
    CHARACTER*(*), INTENT(IN) :: inputchar

    INTEGER :: i,j,n

    n=LEN(inputchar)

    DO i=1,132
       stripchar(i:i)=''
    ENDDO

    j=0
    DO i=1,n
       IF (inputchar(i:i) /= '') THEN
          j=j+1
          stripchar(j:j)=inputchar(i:i)
       ENDIF
    ENDDO

  END FUNCTION stripchar

  SUBROUTINE ConvertChar(inchar,outn)
    CHARACTER*(*),INTENT(IN) :: inchar
    INTEGER, INTENT(OUT) :: outn

    INTEGER :: i, j, k, n, fac

    n=LEN(inchar)

    fac=1; outn=0
    DO i=n,1,-1
       IF (inchar(i:i)==''.OR.inchar(i:i)=="-") EXIT
       j=ICHAR(inchar(i:i))-48
       outn=outn+fac*j
       fac=fac*10
    ENDDO

    outn=outn
    WRITE(6,*) 'exiting ConvertChar', inchar,outn
  END SUBROUTINE ConvertChar

  SUBROUTINE extractword(wordindex,stringin,stringout)
    INTEGER, INTENT(IN) :: wordindex
    CHARACTER(*), INTENT(IN) :: stringin
    CHARACTER(*), INTENT(OUT) :: stringout

    INTEGER :: i,j,n,str,fin,icount

    stringout=''
    n=LEN(stringin)
    i=INDEX(stringin,'!');IF (i==0) i=n
    j=INDEX(stringin,'#');IF (j==0) j=n
    n=MIN(i,j,n)
    str=1;fin=n

    DO icount=1,MAX(1,wordindex-1)
       DO i=str,n
          IF (stringin(i:i)/=' ') EXIT
       ENDDO
       str=i
       IF (n>str) THEN
          DO i=str+1,n
             IF(stringin(i:i)==' ') EXIT
          ENDDO
          fin=i
       ENDIF
       IF (wordindex>2) THEN
          IF (fin<n) THEN
             str=fin+1
          ELSE
             EXIT
          ENDIF
       ENDIF
    ENDDO

    IF (wordindex>1) THEN
       IF (fin>=n) RETURN
       DO i=fin+1,n
          IF (stringin(i:i)/=' ') EXIT
       ENDDO
       str=i
       IF (n>str) THEN
          DO i=str+1,n
             IF(stringin(i:i)==' ') EXIT
          ENDDO
          fin=i
       ENDIF
    ENDIF

    stringout=stringin(str:fin)

  END SUBROUTINE extractword


  !*****************************************************************************
  !  !!! Routine written by Alan Tackett
  ! UpperCase - Converts a string to Upper Case
  !
  ! str     - String to convert
  !
  !*****************************************************************************

  SUBROUTINE UpperCase(str)
    CHARACTER(*), INTENT(INOUT) :: str

    INTEGER  :: i, j, k

    j = LEN(Str)

    DO i=1, j
       k = IACHAR(str(i:i))
       IF ((k>96) .AND. (k<123)) str(i:i) = ACHAR(k-32)
    END DO

    RETURN
  END SUBROUTINE UpperCase


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Logical function  checkline2
!     Returns .true. if inputline matchs either in1 or in2
!     Note: in1 and in2 should either be integers or uppercase characters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LOGICAL FUNCTION checkline2(inputline,in1,in2)
     CHARACTER(*), INTENT(IN) :: inputline
     CHARACTER(*), INTENT(IN) :: in1,in2

     INTEGER :: i,len1,len2
     CHARACTER(120) :: inputline_u

     inputline_u=trim(inputline)
     call UpperCase(inputline_u)
     len1=len(trim(in1));len2=len(trim(in2))

     checkline2=(inputline_u(1:len1)==trim(in1).or.inputline_u(1:len2)==trim(in2))

     RETURN
  END FUNCTION checkline2

  !*****************************************************************
  ! subroutine conthomas(n,o,d,sol)
  !*****************************************************************
  ! use Thomas's algorithm for inverting matrix
  !   Dale U. von Rosenberg, "Methods for the Numerical Solution of
  !     Partial Differential Equations,
  !        Am. Elsevier Pub., 1969, pg. 113
  !   On input, sol contains the RHS of the equation
  !   On ouput, sol contains the solution of the equation
  !    Equation:  o*sol(i-1)+d*sol(i)+o*sol(i+1) = RHS(i)
  !      sol(1)==sol(n+1)==0
  !    simplified version for constant tridiagonal terms --

  SUBROUTINE conthomas(n,o,d,sol)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: o,d
    REAL(8), INTENT(INOUT) :: sol(:)

    REAL(8), ALLOCATABLE :: a(:),b(:)
    REAL(8) :: ss2
    INTEGER :: i

    ALLOCATE(a(n),b(n),stat=i)
    IF (i /= 0) THEN
       WRITE(6,*) 'Thomas: allocation error ', i,n
       STOP
    ENDIF

    a(2)=d
    ss2=o*o
    DO i=3,n
       a(i)=d-ss2/a(i-1)
    ENDDO
    b(2)=sol(2)/d
    DO i=3,n
       b(i)=(sol(i)-o*b(i-1))/a(i)
    ENDDO
    sol(n)=b(n)
    DO i=n-1,2,-1
       sol(i)=b(i)-o*sol(i+1)/a(i)
    ENDDO
    sol(1)=0
    DEALLOCATE(a,b)
  END SUBROUTINE conthomas

  !  subroutine to use Thomas's algorithm to solve tridiagonal matrix
  !   Dale U. von Rosenberg, "Methods for the Numerical Solution of
  !     Partial Differential Equations, Am. Elsevier Pub., 1969, pg. 113
  !      a(i)*u(i-1)+b(i)*u(i)+c(i)*u(i+1)=d(i)
  !      a(1)=c(n)==0
  !      upon return, d(i)=u(i), a,b, & c modified
  SUBROUTINE thomas(n,a,b,c,d)
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(INOUT) :: a(:),b(:),c(:),d(:)

    INTEGER :: i
    IF (n.LT.1) THEN
       WRITE(6,*) '***error in thomas -- n = ',n
       STOP
    ELSE IF (n.EQ.1) THEN
       d(1)=d(1)/b(1)
       RETURN
    ELSE
       DO i=2,n
          b(i)=b(i)-a(i)*c(i-1)/b(i-1)
       ENDDO
       d(1)=d(1)/b(1)
       DO i=2,n
          d(i)=(d(i)-a(i)*d(i-1))/b(i)
       ENDDO
       DO i=n-1,1,-1
          d(i)=d(i)-c(i)*d(i+1)/b(i)
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE thomas

  FUNCTION intjl(rc,q,ff,l)
    ! Compute the integral from rc to inf of
    !     r**2*f(r)*jl(q*r)
    ! assuming large q and small f(inf)
    ! Expanded up to q**(-5)

    INTEGER :: l
    REAL(8) :: intjl,q,rc,ff(4)
    REAL(8) :: cqr,sqr,qinv,qinv2
    cqr=COS(q*rc);sqr=SIN(q*rc)
    qinv=1.0/q;qinv2=qinv*qinv
    IF(l.EQ.0) THEN
       intjl=cqr*qinv2*(ff(1)-ff(3)*qinv2)&
            &        -sqr*qinv2*qinv*(ff(2)-ff(4)*qinv2)
    ELSE IF(l.EQ.1) THEN
       intjl=sqr*qinv2*(ff(1)-qinv2*(ff(3)+ff(2)/rc-ff(1)/(rc*rc)))&
            &        +cqr*qinv2*qinv*(ff(2)+ff(1)/rc-qinv2*(ff(4)&
            &                         +ff(3)/rc-2.0*(ff(2)-ff(1)/rc)/(rc*rc)))
    ELSE IF(l.EQ.2) THEN
       intjl=cqr*qinv2*(-ff(1)+qinv2*(ff(3)+3.0*ff(2)/rc))&
            &        +sqr*qinv2*qinv*(ff(2)+3.0*ff(1)/rc&
            &                  -qinv2*(ff(4)+3.0*(ff(3)-ff(2)/rc)/rc))
    ELSE IF(l.EQ.3) THEN
       intjl=ff(1)*(cqr*(-6.d0*qinv**3/rc-3.d0*qinv**5/rc**3)+&
            &         sqr*(9.d0*qinv**4/rc**2-qinv2)) +&
            &         ff(2)*(cqr*(3.d0*qinv**5/rc**2-qinv**3)+&
            &         sqr*6.d0*qinv**4/rc) +&
            &         ff(3)*(6.d0*cqr*qinv**5/rc+sqr*qinv**4)+&
            &         ff(4)*cqr*qinv**5
    ELSE
       STOP 'Error in intjl: l too large !'
    ENDIF

  END FUNCTION intjl

  !******************************************************************
  ! subroutine jbessel(bes,besp,bespp,ll,order,xx)
  !    Spherical bessel function and derivatives
  !******************************************************************
  SUBROUTINE jbessel(bes,besp,bespp,ll,order,xx)

    INTEGER,INTENT(IN) :: ll,order
    REAL(8),INTENT(IN) :: xx
    REAL(8),INTENT(OUT) :: bes,besp,bespp

    INTEGER,PARAMETER :: imax=40
    REAL(8),PARAMETER :: prec=1.d-15
    INTEGER :: ii,il
    REAL(8) :: besp1,fact,factp,factpp,jn,jnp,jnpp,jr,xx2,xxinv


    IF (order>2) STOP "Wrong order in jbessel !"

    IF (ABS(xx)<prec) THEN
       bes=0.d0;IF (ll==0) bes=1.d0
       IF (order>=1) THEN
          besp=0.d0;IF (ll==1) besp=1.d0/3.d0
       ENDIF
       IF (order==2) THEN
          bespp=0.d0
          IF (ll==0) bespp=-1.d0/3.d0
          IF (ll==2) bespp=2.d0/15.d0
       ENDIF
       RETURN
    ENDIF

    xxinv=1.d0/xx

    IF (xx<1.d0) THEN
       xx2=0.5d0*xx*xx
       fact=1.D0;DO il=1,ll;fact=fact*xx/DBLE(2*il+1);ENDDO
       jn=1.D0;jr=1.D0;ii=0
       DO WHILE(ABS(jr)>=prec.AND.ii<imax)
          ii=ii+1;jr=-jr*xx2/DBLE(ii*(2*(ll+ii)+1))
          jn=jn+jr
       ENDDO
       bes=jn*fact
       IF (ABS(jr)>prec) STOP 'Error: Bessel function did not converge !'
       IF (order>=1) THEN
          factp=fact*xx/DBLE(2*ll+3)
          jnp=1.D0;jr=1.D0;ii=0
          DO WHILE(ABS(jr)>=prec.AND.ii<imax)
             ii=ii+1;jr=-jr*xx2/DBLE(ii*(2*(ll+ii)+3))
             jnp=jnp+jr
          ENDDO
          besp=-jnp*factp+jn*fact*xxinv*DBLE(ll)
          IF (ABS(jr)>prec) STOP 'Error: 1st der. of Bessel function did not converge !'
       ENDIF
       IF (order==2) THEN
          factpp=factp*xx/DBLE(2*ll+5)
          jnpp=1.D0;jr=1.D0;ii=0
          DO WHILE(ABS(jr)>=prec.AND.ii<imax)
             ii=ii+1;jr=-jr*xx2/DBLE(ii*(2*(ll+ii)+5))
             jnpp=jnpp+jr
          ENDDO
          besp1=-jnpp*factpp+jnp*factp*xxinv*DBLE(ll+1)
          IF (ABS(jr)>prec) STOP 'Error: 2nd der. of Bessel function did not converge !'
       ENDIF
    ELSE
       jn =SIN(xx)*xxinv
       jnp=(-COS(xx)+jn)*xxinv
       DO il=2,ll+1
          jr=-jn+DBLE(2*il-1)*jnp*xxinv
          jn=jnp;jnp=jr
       ENDDO
       bes=jn
       IF (order>=1) besp =-jnp+jn *xxinv*DBLE(ll)
       IF (order==2) besp1= jn -jnp*xxinv*DBLE(ll+2)
    ENDIF

    IF (order==2) bespp=-besp1+besp*ll*xxinv-bes*ll*xxinv*xxinv

  END SUBROUTINE jbessel


  !******************************************************************
  ! subroutine solvbes(root,alpha,l,nq)
  !    Find nq first roots of instrinsic equation:
  !                            alpha.jl(Q) + beta.Q.djl/dr(Q) = 0
  !******************************************************************
  SUBROUTINE solvbes(root,alpha,beta,ll,nq)

    INTEGER,INTENT(IN) :: ll,nq
    REAL(8),INTENT(IN) :: alpha,beta
    REAL(8),INTENT(OUT) :: root(nq)

    REAL(8),PARAMETER :: dh=1.D-1, tol=1.D-14

    INTEGER :: nroot
    REAL(8) :: dum,y1,y2,jbes,jbesp,qq,qx,hh

    qq=dh;nroot=0

    DO WHILE (nroot<nq)
       CALL jbessel(jbes,jbesp,dum,ll,1,qq)
       y1=alpha*jbes+beta*qq*jbesp
       qq=qq+dh
       CALL jbessel(jbes,jbesp,dum,ll,1,qq)
       y2=alpha*jbes+beta*qq*jbesp

       DO WHILE (y1*y2>=0.D0)
          qq=qq+dh
          CALL jbessel(jbes,jbesp,dum,ll,1,qq)
          y2=alpha*jbes+beta*qq*jbesp
       ENDDO

       hh=dh;qx=qq
       DO WHILE (hh>tol)
          hh=0.5D0*hh
          IF (y1*y2<0) THEN
             qx=qx-hh
          ELSE
             qx=qx+hh
          ENDIF
          CALL jbessel(jbes,jbesp,dum,ll,1,qx)
          y2=alpha*jbes+beta*qx*jbesp
       ENDDO
       nroot=nroot+1
       root(nroot)=qx

    ENDDO

  END SUBROUTINE solvbes


  !******************************************************************
  ! subroutine shapebes(al,ql,ll,rc)
  !    Find al and ql parameters for a "Bessel" shape function:
  !    Shape(r)=al1.jl(ql1.r)+al2.jl(ql2.r)
  !      such as Shape(r) and 2 derivatives are zero at r=rc
  !              Intg_0_rc[Shape(r).r^(l+2).dr]=1
  !******************************************************************
  SUBROUTINE shapebes(al,ql,ll,rc)

    INTEGER,INTENT(IN) :: ll
    REAL(8),INTENT(IN) :: rc
    REAL(8),INTENT(OUT) :: al(2),ql(2)

    INTEGER :: i
    REAL(8) :: alpha,beta,det,qr,jbes,jbesp,jbespp,amat(2,2),bb(2)

    alpha=1.D0;beta=0.D0
    CALL solvbes(ql,alpha,beta,ll,2)
    ql(1:2)=ql(1:2)/rc

    DO i=1,2
       qr=ql(i)*rc
       CALL jbessel(jbes,jbesp,jbespp,ll,1,qr)
       amat(1,i)=jbesp*ql(i)
       CALL jbessel(jbes,jbesp,jbespp,ll+1,0,qr)
       amat(2,i)=jbes*rc**(ll+2)/ql(i)  !  Intg_0_rc[jl(qr).r^(l+2).dr]
    ENDDO
    bb(1)=0.d0;bb(2)=1.d0

    det=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
    al(1)=(amat(2,2)*bb(1)-amat(1,2)*bb(2))/det
    al(2)=(amat(1,1)*bb(2)-amat(2,1)*bb(1))/det

  END SUBROUTINE shapebes


  !******************************************************************
  SUBROUTINE CALERF(ARG,RESULT,JINT)
    ! ------------------------------------------------------------------
    !
    ! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
    !   for a real argument  x.  It contains three FUNCTION type
    !   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
    !   and one SUBROUTINE type subprogram, CALERF.  The calling
    !   statements for the primary entries are:
    !
    !                   Y=ERF(X)     (or   Y=DERF(X)),
    !
    !                   Y=ERFC(X)    (or   Y=DERFC(X)),
    !   and
    !                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
    !
    !   The routine  CALERF  is intended for internal packet use only,
    !   all computations within the packet being concentrated in this
    !   routine.  The function subprograms invoke  CALERF  with the
    !   statement
    !
    !          CALL CALERF(ARG,RESULT,JINT)
    !
    !   where the parameter usage is as follows
    !
    !      Function                     Parameters for CALERF
    !       call              ARG                  Result          JINT
    !
    !     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
    !     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
    !     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
    !
    !   The main computation evaluates near-minimax approximations
    !   from "Rational Chebyshev approximations for the error function"
    !   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
    !   transportable program uses rational functions that theoretically
    !   approximate  erf(x)  and  erfc(x)  to at least 18 significant
    !   decimal digits.  The accuracy achieved depends on the arithmetic
    !   system, the compiler, the intrinsic functions, and proper
    !   selection of the machine-dependent constants.
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! Explanation of machine-dependent constants
    !
    !   XMIN   = the smallest positive floating-point number.
    !   XINF   = the largest positive finite floating-point number.
    !   XNEG   = the largest negative argument acceptable to ERFCX;
    !            the negative of the solution to the equation
    !            2*exp(x*x) = XINF.
    !   XSMALL = argument below which erf(x) may be represented by
    !            2*x/sqrt(pi)  and above which  x*x  will not underflow.
    !            A conservative value is the largest machine number X
    !            such that   1.0 + X = 1.0   to machine precision.
    !   XBIG   = largest argument acceptable to ERFC;  solution to
    !            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
    !            W(x) = exp(-x*x)/[x*sqrt(pi)].
    !   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
    !            machine precision.  A conservative value is
    !            1/[2*sqrt(XSMALL)]
    !   XMAX   = largest acceptable argument to ERFCX; the minimum
    !            of XINF and 1/[sqrt(pi)*XMIN].
    !
    !   Approximate values for some important machines are:
    !
    !                          XMIN       XINF        XNEG     XSMALL
    !
    !  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
    !  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
    !  IEEE (IBM/XT,
    !    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
    !  IEEE (IBM/XT,
    !    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
    !  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
    !  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
    !  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
    !  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
    !
    !
    !                          XBIG       XHUGE       XMAX
    !
    !  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
    !  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
    !  IEEE (IBM/XT,
    !    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
    !  IEEE (IBM/XT,
    !    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
    !  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
    !  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
    !  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
    !  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
    !
    !*******************************************************************
    !*******************************************************************
    !
    ! Error returns
    !
    !  The program returns  ERFC = 0      for  ARG .GE. XBIG;
    !
    !                       ERFCX = XINF  for  ARG .LT. XNEG;
    !      and
    !                       ERFCX = 0     for  ARG .GE. XMAX.
    !
    !
    ! Intrinsic functions required are:
    !
    !     ABS, AINT, EXP
    !
    !
    !  Author: W. J. Cody
    !          Mathematics and Computer Science Division
    !          Argonne National Laboratory
    !          Argonne, IL 60439
    !
    !  Latest modification: March 19, 1990
    !
    !------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER I,JINT
    !S    REAL
    DOUBLE PRECISION ::  &
&        A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI,  &
&        TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL,&
&        Y,YSQ,ZERO
    DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
    LOGICAL :: goto300
    !------------------------------------------------------------------
    !  Mathematical constants
    !------------------------------------------------------------------
    !S    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/,
    !S   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/,
    !S   2     SIXTEN/16.0E0/
    DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/, &
&        SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/,     &
&        SIXTEN/16.0D0/
    !------------------------------------------------------------------
    !  Machine-dependent constants
    !------------------------------------------------------------------
    !S    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/,
    !S   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/
    DATA XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/,  &
&        XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/
    !------------------------------------------------------------------
    !  Coefficients for approximation to  erf  in first interval
    !------------------------------------------------------------------
    !S    DATA A/3.16112374387056560E00,1.13864154151050156E02,
    !S   1       3.77485237685302021E02,3.20937758913846947E03,
    !S   2       1.85777706184603153E-1/
    !S    DATA B/2.36012909523441209E01,2.44024637934444173E02,
    !S   1       1.28261652607737228E03,2.84423683343917062E03/
    DATA A/3.16112374387056560D00,1.13864154151050156D02, &
&        3.77485237685302021D02,3.20937758913846947D03,   &
         1.85777706184603153D-1/
    DATA B/2.36012909523441209D01,2.44024637934444173D02,  &
&        1.28261652607737228D03,2.84423683343917062D03/
    !------------------------------------------------------------------
    !  Coefficients for approximation to  erfc  in second interval
    !------------------------------------------------------------------
    !S    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
    !S   1       6.61191906371416295E01,2.98635138197400131E02,
    !S   2       8.81952221241769090E02,1.71204761263407058E03,
    !S   3       2.05107837782607147E03,1.23033935479799725E03,
    !S   4       2.15311535474403846E-8/
    !S    DATA D/1.57449261107098347E01,1.17693950891312499E02,
    !S   1       5.37181101862009858E02,1.62138957456669019E03,
    !S   2       3.29079923573345963E03,4.36261909014324716E03,
    !S   3       3.43936767414372164E03,1.23033935480374942E03/
    DATA C/5.64188496988670089D-1,8.88314979438837594D0, &
&       6.61191906371416295D01,2.98635138197400131D02, &
&        8.81952221241769090D02,1.71204761263407058D03, &
&        2.05107837782607147D03,1.23033935479799725D03, &
&        2.15311535474403846D-8/
    DATA D/1.57449261107098347D01,1.17693950891312499D02, &
&        5.37181101862009858D02,1.62138957456669019D03,  &
&        3.29079923573345963D03,4.36261909014324716D03,  &
&        3.43936767414372164D03,1.23033935480374942D03/
    !------------------------------------------------------------------
    !  Coefficients for approximation to  erfc  in third interval
    !------------------------------------------------------------------
    !S    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
    !S   1       1.25781726111229246E-1,1.60837851487422766E-2,
    !S   2       6.58749161529837803E-4,1.63153871373020978E-2/
    !S    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
    !S   1       5.27905102951428412E-1,6.05183413124413191E-2,
    !S   2       2.33520497626869185E-3/
    DATA P/3.05326634961232344D-1,3.60344899949804439D-1, &
&        1.25781726111229246D-1,1.60837851487422766D-2, &
&        6.58749161529837803D-4,1.63153871373020978D-2/
    DATA Q/2.56852019228982242D00,1.87295284992346047D00, &
&        5.27905102951428412D-1,6.05183413124413191D-2, &
&        2.33520497626869185D-3/
    !------------------------------------------------------------------
    Goto300 = .FALSE.
    X = ARG
    Y = ABS(X)
    IF (Y .LE. THRESH) THEN
       !------------------------------------------------------------------
       !  Evaluate  erf  for  |X| <= 0.46875
       !------------------------------------------------------------------
       YSQ = ZERO
       IF (Y .GT. XSMALL) YSQ = Y * Y
       XNUM = A(5)*YSQ
       XDEN = YSQ
       DO I = 1, 3
          XNUM = (XNUM + A(I)) * YSQ
          XDEN = (XDEN + B(I)) * YSQ
       END DO
       RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
       IF (JINT .NE. 0) RESULT = ONE - RESULT
       IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
       RETURN        !*** goto 800
       !------------------------------------------------------------------
       !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
       !------------------------------------------------------------------
    ELSE IF (Y .LE. FOUR) THEN
       XNUM = C(9)*Y
       XDEN = Y
       DO I = 1, 7
          XNUM = (XNUM + C(I)) * Y
          XDEN = (XDEN + D(I)) * Y
       END DO
       RESULT = (XNUM + C(8)) / (XDEN + D(8))
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
       !------------------------------------------------------------------
       !  Evaluate  erfc  for |X| > 4.0
       !------------------------------------------------------------------
    ELSE
       Goto300 = .FALSE.
       RESULT = ZERO
       IF (Y .GE. XBIG) THEN
          IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GOTO300 = .TRUE.
          IF ((Y .GE. XHUGE) .AND. (.NOT. Goto300)) THEN
             RESULT = SQRPI / Y
             GOTO300 = .TRUE.
          END IF
       END IF

       IF (.NOT. Goto300) THEN
          YSQ = ONE / (Y * Y)
          XNUM = P(6)*YSQ
          XDEN = YSQ
          DO I = 1, 4
             XNUM = (XNUM + P(I)) * YSQ
             XDEN = (XDEN + Q(I)) * YSQ
          END DO
          RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
          RESULT = (SQRPI -  RESULT) / Y
          IF (JINT .NE. 2) THEN
             YSQ = AINT(Y*SIXTEN)/SIXTEN
             DEL = (Y-YSQ)*(Y+YSQ)
             RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
          END IF
       END IF
    END IF
    !------------------------------------------------------------------
    !  Fix up for negative argument, erf, etc.
    !------------------------------------------------------------------
    IF (JINT .EQ. 0) THEN                    !*** Goto300 ***
       RESULT = (HALF - RESULT) + HALF
       IF (X .LT. ZERO) RESULT = -RESULT
    ELSE IF (JINT .EQ. 1) THEN
       IF (X .LT. ZERO) RESULT = TWO - RESULT
    ELSE
       IF (X .LT. ZERO) THEN
          IF (X .LT. XNEG) THEN
             RESULT = XINF
          ELSE
             YSQ = AINT(X*SIXTEN)/SIXTEN
             DEL = (X-YSQ)*(X+YSQ)
             Y = EXP(YSQ*YSQ) * EXP(DEL)
             RESULT = (Y+Y) - RESULT
          END IF
       END IF
    END IF
    RETURN           !*** 800
    !---------- Last card of CALERF ----------
  END SUBROUTINE CALERF
  !S    REAL FUNCTION ERF(X)
  DOUBLE PRECISION FUNCTION DERF(X)
    !--------------------------------------------------------------------
    !
    ! This subprogram computes approximate values for erf(x).
    !   (see comments heading CALERF).
    !
    !   Author/date: W. J. Cody, January 8, 1985
    !
    !--------------------------------------------------------------------
    INTEGER JINT
    !S    REAL             X, RESULT
    DOUBLE PRECISION :: X, RESULT
    !------------------------------------------------------------------
    JINT = 0
    CALL CALERF(X,RESULT,JINT)
    !S    ERF = RESULT
    DERF = RESULT
    RETURN
    !---------- Last card of DERF ----------
  END FUNCTION DERF
  !S    REAL FUNCTION ERFC(X)
  DOUBLE PRECISION FUNCTION DERFC(X)
    !--------------------------------------------------------------------
    !
    ! This subprogram computes approximate values for erfc(x).
    !   (see comments heading CALERF).
    !
    !   Author/date: W. J. Cody, January 8, 1985
    !
    !--------------------------------------------------------------------
    INTEGER JINT
    !S    REAL             X, RESULT
    DOUBLE PRECISION :: X, RESULT
    !------------------------------------------------------------------
    JINT = 1
    CALL CALERF(X,RESULT,JINT)
    !S    ERFC = RESULT
    DERFC = RESULT
    RETURN
    !---------- Last card of DERFC ----------
  END FUNCTION DERFC
  !S    REAL FUNCTION ERFCX(X)
  DOUBLE PRECISION FUNCTION DERFCX(X)
    !------------------------------------------------------------------
    !
    ! This subprogram computes approximate values for exp(x*x) * erfc(x).
    !   (see comments heading CALERF).
    !
    !   Author/date: W. J. Cody, March 30, 1987
    !
    !------------------------------------------------------------------
    INTEGER JINT
    !S    REAL             X, RESULT
    DOUBLE PRECISION :: X, RESULT
    !------------------------------------------------------------------
    JINT = 2
    CALL CALERF(X,RESULT,JINT)
    !S    ERFCX = RESULT
    DERFCX = RESULT
    RETURN
    !---------- Last card of DERFCX ----------
  END FUNCTION DERFCX

  SUBROUTINE SolveAXeqB(n,A,B,conditionNo)
    ! General purpose AX=B solver for A, B real.
    ! On return, B stores X
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: A(:,:)
    REAL(8), INTENT(INOUT) :: B(:)
    REAL(8), INTENT(IN), OPTIONAL :: conditionNo

    REAL(8), ALLOCATABLE :: C(:,:),U(:,:),VT(:,:),X(:)
    REAL(8), ALLOCATABLE :: WORK(:)
    REAL(8), ALLOCATABLE :: S(:)
    REAL(8), PARAMETER :: rtol=1.d-9
    INTEGER :: i,j,LWORK
    REAL(8) :: xx,tol
    REAL(8), PARAMETER :: one=1,zero=0

    IF (n == 1) THEN
       B(1)=B(1)/A(1,1)
       RETURN
    ENDIF
    LWORK=MAX(200,n*n)

    ALLOCATE(C(n,n),X(n),U(n,n),VT(n,n),WORK(LWORK),S(n),STAT=i)
    IF (i /= 0) THEN
       WRITE(6,*) 'SolveAXeqB: allocation ', n,Lwork,i
       STOP
    ENDIF

    tol=rtol
    IF (PRESENT(conditionNo)) tol=1.d0/conditionNo
    !WRITE(6,*) ' In SolveAXeqB ', n
    !CALL flush(6)

    C(1:n,1:n)=A(1:n,1:n)
    CALL DGESVD('A','A',n,n,C,n,S,U,n,VT,n,WORK,LWORK,i)

    !DO i=1,n
    !   WRITE(6,*) i,S(i)
    !   CALL flush(6)
    !ENDDO

    tol=tol*S(1)
    X=0
    DO i=1,n
       write(6,*) 'Solver, tol ',i,S(i),tol
       IF (S(i)>tol) THEN
          xx=DOT_PRODUCT(U(1:n,i),B(1:n))/S(i)
          X(1:n)=X(1:n)+xx*(VT(i,1:n))
       ENDIF
    ENDDO

    B=X

    DEALLOCATE(C,X,U,VT,WORK,S)
  END SUBROUTINE SolveAXeqB

  SUBROUTINE SolveAXeqBM(n,A,B,many)
    ! General purpose AX=B solver for A, B real.
    !  Modified version using "many" largest singular values
    ! On return, B stores X
    INTEGER, INTENT(IN) :: n,many
    REAL(8), INTENT(IN) :: A(:,:)
    REAL(8), INTENT(INOUT) :: B(:)

    REAL(8), ALLOCATABLE :: C(:,:),U(:,:),VT(:,:),X(:)
    REAL(8), ALLOCATABLE :: WORK(:)
    REAL(8), ALLOCATABLE :: S(:)
    REAL(8), PARAMETER :: rtol=1.d-9
    INTEGER :: i,j,LWORK
    REAL(8) :: xx,tol
    REAL(8), PARAMETER :: one=1,zero=0

    If (many<1.or.many>n) then
       Write(6,*) 'Error in SOlveAXeqBM -- n, many ', n, many
       stop
    endif
    IF (n == 1) THEN
       B(1)=B(1)/A(1,1)
       RETURN
    ENDIF
    LWORK=MAX(200,n*n)

    ALLOCATE(C(n,n),X(n),U(n,n),VT(n,n),WORK(LWORK),S(n),STAT=i)
    IF (i /= 0) THEN
       WRITE(6,*) 'SolveAXeqB: allocation ', n,Lwork,i
       STOP
    ENDIF


    C(1:n,1:n)=A(1:n,1:n)
    CALL DGESVD('A','A',n,n,C,n,S,U,n,VT,n,WORK,LWORK,i)


    X=0
    DO i=1,n
       if (i<=many) then
         write(6,*) 'Used SV  ',i,S(i)
          xx=DOT_PRODUCT(U(1:n,i),B(1:n))/S(i)
          X(1:n)=X(1:n)+xx*(VT(i,1:n))
       else
         write(6,*) 'UnUsed SV  ',i,S(i)
       ENDIF
    ENDDO

    B=X

    DEALLOCATE(C,X,U,VT,WORK,S)
  END SUBROUTINE SolveAXeqBM

  SUBROUTINE SolveAXeqB_eig(n,A,B,conditionNo)
    !  AX=B solver for A, B real and A symmetric
    ! On return, B stores X
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: A(:,:)
    REAL(8), INTENT(INOUT) :: B(:)
    REAL(8), INTENT(IN), OPTIONAL :: conditionNo

    REAL(8), ALLOCATABLE :: C(:,:),U(:,:),X(:)
    REAL(8), ALLOCATABLE :: WORK(:)
    REAL(8), ALLOCATABLE :: S(:)
    REAL(8), PARAMETER :: rtol=1.d-9
    INTEGER :: i,j,LWORK
    REAL(8) :: xx,tol
    REAL(8), PARAMETER :: one=1,zero=0

    IF (n == 1) THEN
       B(1)=B(1)/A(1,1)
       RETURN
    ENDIF
    LWORK=MAX(200,n*n)

    ALLOCATE(C(n,n),X(n),WORK(LWORK),S(n),STAT=i)
    IF (i /= 0) THEN
       WRITE(6,*) 'SolveAXeqB: allocation ', n,Lwork,i
       STOP
    ENDIF

    tol=rtol
    IF (PRESENT(conditionNo)) tol=1.d0/conditionNo
    !WRITE(6,*) ' In SolveAXeqB ', n
    !CALL flush(6)

    C(1:n,1:n)=A(1:n,1:n)
    CALL DSYEV('V','U',n,C,n,S,WORK,LWORK,i)

    !DO i=1,n
    !   WRITE(6,*) i,S(i)
    !   CALL flush(6)
    !ENDDO

    tol=tol*MAXVAL(S(1:n))
    X=0
    DO i=1,n
       IF (S(i)>tol) THEN
          xx=DOT_PRODUCT(C(1:n,i),B(1:n))/S(i)
          X(1:n)=X(1:n)+xx*(C(1:n,i))
       ENDIF
    ENDDO

    B=X

    DEALLOCATE(C,X,WORK,S)
  END SUBROUTINE SolveAXeqB_eig

  SUBROUTINE SVDINVERSE(m,n,tol,MN,RINVMN)      !MN(m,n)   RINVMN(n,m)
    INTEGER, INTENT(IN) :: n,m
    REAL(8), INTENT(IN) :: tol
    REAL(8), INTENT(INOUT) :: MN(:,:),RINVMN(:,:)

    INTEGER :: i,j,k,kk,LWORK,LIWORK,ok
    REAL(8), ALLOCATABLE :: u(:,:),vt(:,:),s(:),work(:)
    INTEGER, ALLOCATABLE :: iwork(:)

    LWORK=8*n*m; LIWORK=8*MAX(n,m); kk=MIN(n,m)
    allocate (u(m,m),vt(n,n),s(kk),&
&         work(LWORK),IWORK(LIWORK) )

    s=0.d0;  u=0.d0;  vt=0.d0
    call  dgesdd('A',m,n,MN,m,s,u,m,vt,n,work,LWORK,IWORK,ok)
    write(6,*) 'dgesdd completed with info = ', ok
    !CALL DGESVD('A','A',m,n,MN,m,S,U,m,VT,n,WORK,LWORK,i)
    RINVMN=0
    do i=1,kk
       write(6,*) 'i s = ', i,s(i)
       if (s(i)>tol) then
       do j=1,n
          do k=1,m
             RINVMN(j,k)=RINVMN(j,k) + vt(i,j)*u(k,i)/s(i)
          enddo
       enddo
       endif
    enddo

    deallocate (u,vt,s,work,IWORK )
  END SUBROUTINE SVDINVERSE

  FUNCTION ASSOCIATEDLAGUERRE(n,x)
     REAL(8) :: ASSOCIATEDLAGUERRE
     REAL(8), INTENT(IN) :: x
     INTEGER, INTENT(IN) :: n

     REAL(8) :: fac,term,acc
     INTEGER :: i

     if (n<0.or.n>1000) then
        write(6,*) 'ERROR in ASSOCIATEDLAGUERRE -- n = ', n,x
        stop
     endif
     fac=sqrt(float((n+1)*(n+2)))/2; acc=1
     if (n>0) then
        term=1
        do i=1,n
           term=term*x*(-n+i-1)/(3.d0+i-1)/i
           acc=acc+term
        enddo
     endif
     ASSOCIATEDLAGUERRE=acc*fac
  END FUNCTION ASSOCIATEDLAGUERRE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Given integer l1, l2, l3  return 3-j value
!
!           ( l1  l2  l3 )^2
!           (  0   0   0 )
!    Based on subroutine written by Xiao Xu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION THRJ2(l1,l2,l3)
   REAL(8) :: THRJ2
   INTEGER, INTENT(IN) :: l1, l2, l3

   INTEGER ::BigJ,J1,J2,J3

    J1=l1;
    J2=l2;
    J3=l3;
    BigJ=J1+J2+J3;

    ! Following the 3J symbol notes
    IF(MOD(BigJ,2) == 1)  THEN
       THRJ2=0;
    ELSE
       THRJ2 = factorial(BigJ/2)/( &
&           factorial(BigJ/2-J1) *factorial(BigJ/2-J2)*factorial(BigJ/2 - J3));
       THRJ2 = THRJ2*THRJ2
       THRJ2=  THRJ2* ( factorial(BigJ - 2*J1)*factorial(BigJ - 2*J2)* &
&             factorial(BigJ - 2*J3)/factorial(BigJ +1) )
    ENDIF

  END FUNCTION THRJ2

END MODULE GlobalMath

