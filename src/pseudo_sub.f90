MODULE pseudo_sub
  USE globalmath
  USE gridmod
  USE pseudodata
  USE excor
  USE fock

  IMPLICIT NONE

CONTAINS

  !***************************************************************
  ! SUBROUTINE besselps(lmax,Grid,Pot)
  !  Creates screened pseudopotential by simply pseudizing the
  !    AE potential with a l=0 spherical Bessel function:
  !                                     Vps(r) = a.sin(qr)/r
  !***************************************************************
  SUBROUTINE besselps(Grid,Pot,PAW)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW

    INTEGER :: i,irc,l,n
    REAL(8) :: e,rc,alpha,beta,vv,vvp,AA,QQ,xx(1)
    REAL(8),ALLOCATABLE ::  VNC(:),wfn(:)
    REAL(8),POINTER :: r(:),rv(:)
    CHARACTER(132) :: line

    n=Grid%n
    r=>Grid%r
    rv=>Pot%rv
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc

    vv=rv(irc);vvp=Gfirstderiv(Grid,irc,rv)

    alpha=1.D0-rc*vvp/vv;beta=1.D0
    call solvbes(xx,alpha,beta,0,1);QQ=xx(1)
    AA=vv/sin(QQ);QQ=QQ/rc

    PAW%rveff(1)=0.d0
    PAW%rveff(irc+1:n)=rv(irc+1:n)
    do i=2,irc
     PAW%rveff(i)=AA*sin(QQ*r(i))
    enddo

  END SUBROUTINE besselps


    !***************************************************************
    ! SUBROUTINE EvaluateP
    !   Inverts 4x4 matrix used  by kerker subroutine
    !***************************************************************
    SUBROUTINE EvaluateP(m,A,B,C,D,coef)
      INTEGER, INTENT(IN) :: m(4)
      REAL(8), INTENT(IN) :: A,B,C,D
      REAL(8), INTENT(OUT) ::  coef(4)

      REAL(8) :: t(4,4)
      INTEGER :: i,n

      t=0
      Coef(1)=A; Coef(2)=B;  Coef(3)=C;    Coef(4)=D
      t(1,1:4)=1
      t(2,1:4)=m(1:4)
      DO i=1,4
         t(3,i)=m(i)*(m(i)-1)
      ENDDO
      DO i=1,4
         t(4,i)=m(i)*(m(i)-1)*(m(i)-2)
      ENDDO
      n=4
      CALL linsol(t,Coef,n,4,4,4)
    END SUBROUTINE EvaluateP

    !***************************************************************
    ! SUBROUTINE EvaluateTp
    !   Inverts 5x5 matrix used  by troullier subroutine
    !***************************************************************
    SUBROUTINE EvaluateTp(l,A,B,C,D,F,coef)
      INTEGER, INTENT(IN) :: l
      REAL(8), INTENT(IN) :: A,B,C,D,F
      REAL(8), INTENT(OUT) ::  coef(6)

      REAL(8) :: t(6,6),coef10,old
      REAL(8), PARAMETER :: small=1.e-10
      INTEGER :: i,n,iter
      INTEGER, PARAMETER :: niter=1000

      old=-1.e30; Coef10=-1; iter=-1
      DO WHILE (iter < niter .AND. ABS(old-coef10)> small)
         iter=iter+1
         t=0
         Coef(1)=A-Coef10; Coef(2)=B-2*Coef10;  Coef(3)=C-2*Coef10;
         Coef(4)=D;    Coef(5)=F
         Coef(6)=-Coef10**2
         DO i=1,6
            t(1,i)=1
            t(2,i)=2*i
            t(3,i)=2*i*(2*i-1)
            t(4,i)=2*i*(2*i-1)*(2*i-2)
            t(5,i)=2*i*(2*i-1)*(2*i-2)*(2*i-3)
         ENDDO
         t(6,1)=2*Coef10;  t(6,2)=2*l+5

         n=6
         CALL linsol(t,Coef,n,6,6,6)

         old=Coef10; Coef10=Coef10+Coef(1)
         !WRITE(6,'("EvaluateTp: iter",i5,1p,2e15.7)') iter,Coef(1),Coef10
         !WRITE(6,'("Coef: ",1p,6e15.7)')Coef10,(Coef(i),i=2,6)
         Coef(1)=Coef10
      ENDDO

      IF (iter >= niter) THEN
         WRITE(6,*) 'Error in EvaluateTP -- no convergence'
         STOP
      ENDIF
    END SUBROUTINE EvaluateTp

  !*******************************************************************
  !  function to calculated <wfn|O|wfn> for smooth paw wavefunction
  !*******************************************************************
  FUNCTION sepnorm(Grid,PAW,nr,l,wfn)
    REAL(8) :: sepnorm
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: nr,l
    REAL(8), INTENT(IN) :: wfn(:)

    INTEGER :: n,ib,ic,nbase,irc
    REAL(8) :: h
    REAL(8), ALLOCATABLE :: b(:)

    n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc
    ALLOCATE(b(nbase),stat=ib)
    IF (ib/=0) THEN
       WRITE(6,*) 'Error in sepnorm allocation', nbase,ib
       STOP
    ENDIF

    IF (nr<irc) THEN
       WRITE(6,*) 'Error in sepnorm -- nr < irc'
       STOP
    ENDIF
    sepnorm=overlap(Grid,wfn,wfn,1,nr)
    b=0
    DO ib=1,nbase
       IF (l==PAW%l(ib)) b(ib)=overlap(Grid,wfn,PAW%otp(:,ib),1,nr)
    ENDDO
    DO ib=1,nbase
       DO ic=1,nbase
          sepnorm=sepnorm+b(ib)*PAW%oij(ib,ic)*b(ic)
       ENDDO
    ENDDO

    DEALLOCATE(b)
  END FUNCTION sepnorm

  !*******************************************************************
  !  function to calculated <wfn1|O|wfn2> for smooth paw functions
  !*******************************************************************
  FUNCTION genoverlap(Grid,PAW,nr,l,wfn1,wfn2)
    REAL(8) :: genoverlap
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: nr,l
    REAL(8), INTENT(IN) :: wfn1(:),wfn2(:)

    INTEGER :: n,ib,ic,nbase,irc,p
    REAL(8) :: h
    REAL(8), ALLOCATABLE :: a(:),b(:)

    n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc;p=irc+1
    ALLOCATE(a(nbase),b(nbase),stat=ib)
    IF (ib/=0) THEN
       WRITE(6,*) 'Error in genoverlap allocation', nbase,ib
       STOP
    ENDIF

    IF (nr<irc) THEN
       WRITE(6,*) 'Error in genoverlap -- nr < irc'
       STOP
    ENDIF
    genoverlap=overlap(Grid,wfn1,wfn2,1,nr)
    a=0; b=0
    DO ib=1,nbase
       IF (l==PAW%l(ib)) THEN
          a(ib)=overlap(Grid,wfn1,PAW%otp(:,ib))
          b(ib)=overlap(Grid,wfn2,PAW%otp(:,ib))
       ENDIF
    ENDDO
    DO ib=1,nbase
       DO ic=1,nbase
          genoverlap=genoverlap+a(ib)*PAW%oij(ib,ic)*b(ic)
       ENDDO
    ENDDO

    DEALLOCATE(a,b)
  END FUNCTION genoverlap


  ! generalized Graham-Schmidt orthogonalization of wfn1 to wfn2
  ! on output <wfn1|O|wfn2>=0
  SUBROUTINE genOrthog(Grid,PAW,nr,l,wfn1,wfn2)
    !REAL(8) :: genoverlap
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: nr,l
    REAL(8), INTENT(INOUT) :: wfn1(:)
    REAL(8), INTENT(IN) :: wfn2(:)

    REAL(8) :: x,y

    !   Orthogonalize
    x=genoverlap(Grid,PAW,nr,l,wfn1,wfn2)
    y=genoverlap(Grid,PAW,nr,l,wfn2,wfn2)

    WRITE(6,*) 'overlap ', x,y
    wfn1(:)=wfn1(:)-(x/y)*wfn2(:)

  END SUBROUTINE genOrthog

  !**************************************************************
  ! subroutine hatpotL
  !   Calculates potential associated with L component
  !    of unit hat density
  !**************************************************************
  SUBROUTINE hatpotL(Grid,PAW,l,vhat)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(OUT) :: vhat(:)

    INTEGER :: n,irc,i
    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: den(:),a(:)
    REAL(8) :: h,con
    REAL(8) :: qr,jbes1,jbes2,dum1,dum2,al(2),ql(2)

    n=Grid%n
    h=Grid%h
    r=>Grid%r

    irc=PAW%irc

    ALLOCATE(den(n),a(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in hatpotL allocation',n,i
       STOP
    ENDIF


    IF (besselshapefunction) THEN
       CALL shapebes(al,ql,l,PAW%rc_shap)
       DO i=1,PAW%irc_shap
          qr=ql(1)*r(i);CALL jbessel(jbes1,dum1,dum2,l,0,qr)
          qr=ql(2)*r(i);CALL jbessel(jbes2,dum1,dum2,l,0,qr)
          den(i)=(al(1)*jbes1+al(2)*jbes2)*r(i)**2
       ENDDO
       IF (n>PAW%irc_shap) den(PAW%irc_shap+1:n)=0.d0
    ELSE
       DO i=1,n
          den(i)=(r(i)**l)*PAW%hatden(i)
       ENDDO
       a(1:n)=den(1:n)*(r(1:n)**l)
       con=integrator(Grid,a,1,PAW%irc_shap)
       den=den/con
    ENDIF

    a(1:n)=den(1:n)*(r(1:n)**l)
    write(6,*) 'hatpotL check l ', l,integrator(Grid,a)

    vhat=0
    CALL apoisson(Grid,l,n,den,vhat(1:n))

    ! apoisson returns vhat*r
    !DO i=1,n
    !   WRITE (78+l,'(i5,1p,5e15.7)') i,Grid%r(i),den(i),vhat(i)
    !ENDDO
    vhat(2:n)=vhat(2:n)/r(2:n)
    CALL extrapolate(Grid,vhat)

    DEALLOCATE(den,a)
  END SUBROUTINE hatpotL

  !**************************************************************
  ! subroutine hatL
  !   Calculates density associated with L component
  !    normalized to unity
  !**************************************************************
  SUBROUTINE hatL(Grid,PAW,l,dhat)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(OUT) :: dhat(:)

    INTEGER :: n,irc,i
    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: den(:),a(:)
    REAL(8) :: h,con
    REAL(8) :: qr,jbes1,jbes2,dum1,dum2,al(2),ql(2)

    n=Grid%n
    h=Grid%h
    r=>Grid%r

    irc=PAW%irc

    ALLOCATE(den(n),a(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in hatL allocation',irc,i
       STOP
    ENDIF

    IF (besselshapefunction) THEN
       CALL shapebes(al,ql,l,PAW%rc_shap)
       DO i=1,PAW%irc_shap
          qr=ql(1)*r(i);CALL jbessel(jbes1,dum1,dum2,l,0,qr)
          qr=ql(2)*r(i);CALL jbessel(jbes2,dum1,dum2,l,0,qr)
          den(i)=(al(1)*jbes1+al(2)*jbes2)*r(i)**2
       ENDDO
       IF (n>PAW%irc_shap) den(PAW%irc_shap+1:n)=0.d0
    ELSE
       DO i=1,n
          den(i)=(r(i)**l)*PAW%hatden(i)
       ENDDO
       a(1:n)=den(1:n)*(r(1:n)**l)
       con=integrator(Grid,a,1,PAW%irc_shap)
       den=den/con
    ENDIF

    dhat=0

    dhat(1:n)=den(1:n)

    DEALLOCATE(den,a)
  END SUBROUTINE hatL

  !***********************************************************************88
  ! on input: f1(i) and f2(i) are radial wfn * r for angular momentum l
  ! on input: t1(i) and t2(i) are smooth radial wfn * r for angular momentum l
  !   for r > rc, f1=t1, f2=t2
  ! on output: qqqq is difference overlap matrix element
  !   qqqq=<f1|f2>-<t1|t2>
  !***********************************************************************88
  SUBROUTINE dqij(Grid,PAW,ib,ic,qqqq)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: ib,ic
    REAL(8), INTENT(OUT) :: qqqq

    INTEGER :: n,i,ok,irc
    REAL(8) :: h
    REAL(8), ALLOCATABLE :: dum(:)

    qqqq=0
    IF (PAW%l(ib)/=PAW%l(ic)) RETURN
    n=Grid%n; h=Grid%h;  irc=PAW%irc
    ALLOCATE(dum(n),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in dqij allocation', n,ok
       STOP
    ENDIF
    DO i=1,n
       dum(i)=PAW%ophi(i,ib)*PAW%ophi(i,ic)-PAW%otphi(i,ib)*PAW%otphi(i,ic)
    ENDDO
    qqqq=integrator(Grid,dum,1,irc)

    DEALLOCATE(dum)
  END SUBROUTINE dqij

  !***********************************************************************
  ! SUBROUTINE dtij
  ! on input: f1(i) and f2(i) are radial wfn * r for angular momentum l
  ! on input: t1(i) and t2(i) are smooth radial wfn * r for angular momentum l
  !   for r > rc, f1=t1, f2=t2
  ! on output: tij is difference kinetic energy matrix element in Rydberg units
  !   tij =<f1|T|f2>-<t1|T|t2>
  !************************************************************************
  SUBROUTINE dtij(Grid,PAW,ib,ic,tij)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: ib,ic
    REAL(8), INTENT(OUT) :: tij

    INTEGER :: n,i,ok,l,irc
    REAL(8) :: angm
    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: dum(:),del1(:),del2(:),tdel1(:),tdel2(:)

    tij=0
    IF (PAW%l(ib)/=PAW%l(ic)) RETURN
    n=Grid%n;  r=>Grid%r;  l=PAW%l(ib);  irc=PAW%irc
    ALLOCATE(dum(n),del1(n),tdel1(n),del2(n),tdel2(n),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in dtij allocation', n,ok
       STOP
    ENDIF
    CALL derivative(Grid,PAW%ophi(:,ib),del1)
    CALL derivative(Grid,PAW%ophi(:,ic),del2)
    CALL derivative(Grid,PAW%otphi(:,ib),tdel1)
    CALL derivative(Grid,PAW%otphi(:,ic),tdel2)
    dum=0 ;   angm=l*(l+1)
    DO i=1,irc
       dum(i)=del1(i)*del2(i)-tdel1(i)*tdel2(i)
    ENDDO
    del1=0;del2=0;tdel1=0;tdel2=0
    del1(2:irc)=PAW%ophi(2:irc,ib)/Grid%r(2:irc)
    del2(2:irc)=PAW%ophi(2:irc,ic)/Grid%r(2:irc)
    tdel1(2:irc)=PAW%otphi(2:irc,ib)/Grid%r(2:irc)
    tdel2(2:irc)=PAW%otphi(2:irc,ic)/Grid%r(2:irc)
    DO i=1,irc
       dum(i)=dum(i)+angm*(del1(i)*del2(i)-tdel1(i)*tdel2(i))
    ENDDO
    tij=integrator(Grid,dum,1,irc)

    DEALLOCATE(dum,del1,del2,tdel1,tdel2)
  END SUBROUTINE dtij

  !***********************************************************************
  ! SUBROUTINE altdtij
  ! on input: f1(i) and f2(i) are radial wfn * r for angular momentum l
  ! on input: t1(i) and t2(i) are smooth radial wfn * r for angular momentum l
  !   for r > rc, f1=t1, f2=t2
  ! on output: tij is difference kinetic energy matrix element in Rydberg units
  !   tij =<f1|T|f2>-<t1|T|t2>
  !************************************************************************
  SUBROUTINE altdtij(Grid,PAW,ib,ic,tij)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: ib,ic
    REAL(8), INTENT(OUT) :: tij

    INTEGER :: n,i,ok,l,irc
    REAL(8) :: angm
    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: dum(:),tdel1(:),tdel2(:)

    tij=0
    IF (PAW%l(ib)/=PAW%l(ic)) RETURN
    n=Grid%n;  r=>Grid%r;  l=PAW%l(ib);  irc=PAW%irc
    ALLOCATE(dum(n),tdel1(n),tdel2(n),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in dtij allocation', n,ok
       STOP
    ENDIF
    dum=0
    DO i=2,irc
       !dum(i)=(PAW%eig(ic)-AEPot%rv(i)/Grid%r(i))*PAW%ophi(i,ib)*PAW%ophi(i,ic)
       dum(i)=PAW%ophi(i,ib)*PAW%Kop(i,ic)
    ENDDO
    CALL derivative(Grid,PAW%otphi(:,ic),tdel1)
    CALL derivative(Grid,tdel1,tdel2)
    angm=l*(l+1)
    DO i=2,irc
       dum(i)=dum(i)+PAW%otphi(i,ib)*(tdel2(i)-&
&           angm*PAW%otphi(i,ic)/(Grid%r(i)**2))
    ENDDO
    tij=integrator(Grid,dum,1,irc)

    DEALLOCATE(dum,tdel1,tdel2)
  END SUBROUTINE altdtij

  SUBROUTINE dvij(Grid,PAW,FC,nz,ib,ic,vij)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    TYPE(FCInfo), INTENT(IN) :: FC
    INTEGER, INTENT(IN) :: ib,ic
    REAL(8), INTENT(IN) :: nz
    REAL(8), INTENT(OUT) :: vij

    INTEGER :: n,i,ok,irc
    REAL(8) :: h,en,q,qt
    REAL(8), ALLOCATABLE :: dum(:),d1(:)
    REAL(8), POINTER :: r(:)

    vij=0
    IF (PAW%l(ib)/=PAW%l(ic)) RETURN
    n=Grid%n; h=Grid%h;  r=>Grid%r;  irc=PAW%irc
    ALLOCATE(dum(n),d1(n),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in dvij allocation', n,ok
       STOP
    ENDIF

    q=integrator(Grid,FC%coreden)
    WRITE(6,*) 'core electrons ',q,FC%zcore
    CALL poisson(Grid,q,FC%coreden,dum,en)
    dum=dum-2*nz
    qt=integrator(Grid,PAW%tcore)
    WRITE(6,*) 'coretail electrons ',qt
    CALL poisson(Grid,qt,PAW%tcore,d1,en)
    dum(1)=0
    DO i=2,irc
       dum(i)=PAW%ophi(i,ib)*PAW%ophi(i,ic)*dum(i)/r(i)-&
&           PAW%otphi(i,ib)*PAW%otphi(i,ic)*(PAW%vloc(i)+d1(i)/r(i))
    ENDDO
    vij=integrator(Grid,dum,1,irc)

    DEALLOCATE(dum)
  END SUBROUTINE dvij

  SUBROUTINE avij(Grid,PAW,ib,ic,vij)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: ib,ic
    REAL(8), INTENT(OUT) :: vij

    INTEGER :: n,i,ok,irc
    REAL(8) :: h,en,q,qt
    REAL(8), ALLOCATABLE :: dum(:),d1(:)
    REAL(8), POINTER :: r(:)

    vij=0
    IF (PAW%l(ib)/=PAW%l(ic)) RETURN
    n=Grid%n; h=Grid%h;  r=>Grid%r;  irc=PAW%irc
    ALLOCATE(dum(n),d1(n),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in dvij allocation', n,ok
       STOP
    ENDIF

    dum(1)=0
    DO i=2,irc
       dum(i)=(PAW%ophi(i,ib)*PAW%ophi(i,ic)*PAW%AErefrv(i)-&
&           PAW%otphi(i,ib)*PAW%otphi(i,ic)*PAW%rveff(i))/r(i)
    ENDDO
    vij=integrator(Grid,dum,1,irc)

    DEALLOCATE(dum)
  END SUBROUTINE avij

  !********************************************************
  ! SUBROUTINE calcwij
  !  subroutine to accumulate the wij coefficience for an input
  !     smooth wavefunction twfn and occupancy and l
  !********************************************************
  SUBROUTINE calcwij(Grid,PAW,l,occ,twfn,wij)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: twfn(:),occ
    REAL(8), INTENT(INOUT) :: wij(:,:)

    INTEGER :: n,i,ok,irc,ib,ic,nbase
    REAL(8) :: h
    REAL(8), ALLOCATABLE :: bm(:)

    n=Grid%n; h=Grid%h ;  irc=PAW%irc

    nbase=PAW%nbase

    ALLOCATE(bm(nbase))

    bm=0
    DO ib=1,nbase
       IF (l==PAW%l(ib)) bm(ib)=overlap(Grid,PAW%otp(:,ib),twfn,1,irc)
       !IF (l==PAW%l(ib)) write(6,*) 'accum wij',l,ib,bm(ib)
    ENDDO
    DO ib=1,nbase
       DO ic=1,nbase
          wij(ib,ic)=wij(ib,ic) + occ*bm(ib)*bm(ic)
       ENDDO
    ENDDO
    DEALLOCATE(bm)
  END SUBROUTINE calcwij

  SUBROUTINE FillHat(Grid,PAW)
    TYPE(GridInfo) , INTENT(IN):: Grid
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW

    INTEGER :: ll,n,l

    ll=MAXVAL(PAW%TOCCWFN%l(:)); ll=MAX(ll,PAW%lmax); ll=2*ll
    n=Grid%n

    ALLOCATE(PAW%g(n,ll+1))

    DO l=0,ll
       CALL hatL(Grid,PAW,l,PAW%g(:,l+1))
    ENDDO

  END SUBROUTINE FillHat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Get_Energy_EXX_smoothpseudo             !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_Energy_EXX_pseudo(Grid,PAW,eex)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(pseudoinfo), INTENT(IN) :: PAW
    REAL(8), INTENT(OUT) :: eex

    REAL(8), POINTER :: phi(:,:),r(:)
    REAL(8), ALLOCATABLE :: wfp(:),dum(:),ddum(:),ch(:),ar(:)
    REAL(8), ALLOCATABLE :: Fnu(:,:),Lnu(:,:),Snu(:),Mnunup(:,:),Cnu(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,io,jo,ok,norbit,last,nodes
    INTEGER :: nu,nup,ni,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,kappa,test
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=PAW%TOCCwfn%norbit
    phi=>PAW%TOCCwfn%wfn
    r=>grid%r

    ALLOCATE(wfp(n),dum(n),ddum(n),ch(n),ar(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'allocation error in Get_Energy_EXX', i,n
       STOP
    ENDIF

    eex=0;test=0
    DO io=1,norbit
       occ=PAW%TOCCwfn%occ(io)
       IF (occ>threshold.AND.(.NOT.PAW%TOCCwfn%iscore(io))) THEN
          li=PAW%TOCCwfn%l(io);ni=PAW%TOCCwfn%np(io)
          DO jo=1,norbit
             occj=PAW%TOCCwfn%occ(jo)
             IF (occj>threshold.AND.(.NOT.PAW%TOCCwfn%iscore(jo))) THEN
                lj=PAW%TOCCwfn%l(jo);nj=PAW%TOCCwfn%np(jo)
                lmax=li+lj
                lmin=ABS(li-lj)
                wfp(1:n)=phi(1:n,io)*phi(1:n,jo)
                ar(1:n)=PAW%OCCwfn%wfn(1:n,io)*PAW%OCCwfn%wfn(1:n,jo)
                DO ll=lmin,lmax,2
                   CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                   IF (wgt>threshold) THEN
                      wgt=wgt*(2*ll+1)   ! because of apoisson convention
                      call Calc_Moment(Grid,PAW,phi(1:n,io),phi(1:n,jo),&
&                             li,lj,ll,dum)
                      ddum=wfp+dum
                      CALL apoisson(Grid,ll,n,ddum,dum)
                      CALL apoisson(Grid,ll,n,ar,ch)
                      ddum(2:n)=(dum(2:n)*ddum(2:n))/r(2:n)
                      ddum(1)=0
                      eex=eex-wgt*integrator(Grid,ddum)/2
                      ch(2:n)=(ch(2:n)*ar(2:n))/r(2:n); ch(1)=0
                      test=test-wgt*integrator(Grid,ch)/2
                   ENDIF
                ENDDO
             ENDIF
          ENDDO   !jo

       ENDIF
    ENDDO

    write(6,*) 'Get_Energy_EXX_pseudo', eex,test
    DEALLOCATE(wfp,dum,ddum,ch,ar)

  END SUBROUTINE Get_Energy_EXX_pseudo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Calc_Moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Calc_Moment(Grid,PAW,wfn1,wfn2,l1,l2,ll,m)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(pseudoinfo), INTENT(IN) :: PAW
    REAL(8), INTENT(IN) :: wfn1(:),wfn2(:)
    INTEGER, INTENT(IN) :: l1,l2,ll
    REAL(8), INTENT(INOUT) :: m(:)

    REAL(8) :: dot1,dot2
    INTEGER :: i,j,k,n,nbasis,ib,jb
    REAL(8), PARAMETER :: threshold=1.d-8
    REAL(8), ALLOCATABLE :: dum(:)

    nbasis=PAW%nbase
    n=Grid%n
    ALLOCATE(dum(n))

    m=0
    DO ib=1,nbasis
       if (PAW%l(ib)==l1) then
           dot1=overlap(Grid,wfn1,PAW%otp(:,ib)); write(6,*)'dot1',dot1
           do jb=1,nbasis
              if (PAW%l(jb)==l2) then
                 dot2=overlap(Grid,wfn2,PAW%otp(:,jb)); write(6,*)'dot2',dot2
                 m=m+dot1*dot2*PAW%mLij(ib,jb,ll+1)*PAW%g(:,ll+1)
              endif
           enddo
        endif
     Enddo

     dum=m*(Grid%r**ll)
     write(6,*)'Check moment',l1,l2,ll,integrator(Grid,dum)
     deallocate(dum)
  END  SUBROUTINE Calc_Moment


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Get_Energy_EXX_onecenter_pseudo             !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_Energy_EXX_pseudo_one(Grid,PAW,eex)
    TYPE(gridinfo), INTENT(IN) :: Grid
    TYPE(pseudoinfo), INTENT(IN) :: PAW
    REAL(8), INTENT(OUT) :: eex

    REAL(8), POINTER :: phi(:,:),r(:)
    INTEGER :: i,j,k,n,m,l,ll,li,lj,lmin,lmax,io,jo,ok,norbit,last,nodes
    INTEGER :: nu,nup,ni,nj,ib,jb,kb,lb,nbasis
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,kappa,doti,dotj,dotk,dotl
    REAL(8), PARAMETER :: threshold=1.d-8

    n=Grid%n
    norbit=PAW%TOCCwfn%norbit
    phi=>PAW%TOCCwfn%wfn
    r=>grid%r
    nbasis=PAW%nbase


    eex=0
    DO io=1,norbit
       occ=PAW%TOCCwfn%occ(io)
       IF (occ>threshold.AND.(.NOT.PAW%TOCCwfn%iscore(io))) THEN
          li=PAW%TOCCwfn%l(io);ni=PAW%TOCCwfn%np(io)
          DO jo=1,norbit
             occj=PAW%TOCCwfn%occ(jo)
             IF (occj>threshold.AND.(.NOT.PAW%TOCCwfn%iscore(jo))) THEN
                lj=PAW%TOCCwfn%l(jo);nj=PAW%TOCCwfn%np(jo)
                lmax=li+lj
                lmin=ABS(li-lj)
                DO ib=1,nbasis
                   if (PAW%l(ib)==li) then
                    doti=overlap(Grid,PAW%otp(:,ib),phi(:,io))
                    Do jb=1,nbasis
                       If (PAW%l(jb)==lj) then
                        dotj=overlap(Grid,PAW%otp(:,jb),phi(:,jo))
                       Do kb=1,nbasis
                        If (PAW%l(kb)==lj) then
                          dotk=overlap(Grid,PAW%otp(:,kb),phi(:,jo))
                         Do lb=1,nbasis
                          If (PAW%l(lb)==li) then
                              dotl=overlap(Grid,PAW%otp(:,lb),phi(:,io))
                             term=doti*dotj*dotk*dotl
                              DO ll=lmin,lmax,2
                                 CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                                 IF (wgt>threshold) THEN
                                   wgt=wgt*(2*ll+1) !Not in DR
                                   eex=eex-wgt*term*PAW%DR(ib,jb,kb,lb,ll+1)/2
                                 endif
                              enddo  !ll
                           endif
                          enddo  !ib
                         endif
                        enddo  !kb
                       endif
                      enddo  !jb
                     endif
                    enddo  !ib
                  endif
                 enddo    !jo
                endif
               enddo     !io

       write(6,*) 'Get_Energy_EXX_pseudo_one val-val: ', eex

       If   (PAW%ncoreshell>0) then
         do k=1,PAW%ncoreshell
            do ib=1,nbasis
               do jb=1,nbasis
                  if (PAW%l(ib)==PAW%l(jb)) then
                      eex=eex-PAW%wij(ib,jb)*PAW%DRVC(k,ib,jb)
                  endif
               enddo
            enddo
         enddo
         write(6,*) 'Get_Energy_EXX_pseudo_one corrected: ', eex
        endif
  END SUBROUTINE Get_Energy_EXX_pseudo_one



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Fill stored matrix elements and calculate atomic energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SPMatrixElements(Grid,Pot,FC,PAW)
    TYPE(GridInfo) , INTENT(IN):: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(FCInfo), INTENT(IN) :: FC
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW

    INTEGER :: io,jo,ko,lo,ll,lmax,i,j,k,l,n,norbit,nbase,ib,jb,kb,lb,irc
    INTEGER :: llmin, llmax,lp,lr
    REAL(8), ALLOCATABLE :: arg(:),dum(:),rh(:),trh(:),d(:),td(:),wij(:,:)
    REAL(8) :: x,y,z,rr,occ
    REAL(8) :: Qcore,tQcore
    LOGICAL :: even

    n=Grid%n ; irc=PAW%irc; lr=min(irc+20,n) ; nbase=PAW%nbase
    ll=MAXVAL(PAW%TOCCWFN%l(:)); ll=MAX(ll,PAW%lmax); ll=2*ll
    ALLOCATE(arg(n),dum(n),rh(n),trh(n),d(n),td(n))
    arg=0;   dum=0
    ALLOCATE(PAW%mLij(nbase,nbase,ll+1),PAW%DR(nbase,nbase,nbase,nbase,ll+1))
    ALLOCATE(wij(nbase,nbase))
    PAW%mLij=0.d0;PAW%DR=0.d0

    CALL poisson(Grid,Qcore,PAW%core,dum,y)
    PAW%rVf=-2*Pot%nz+dum

    CALL poisson(Grid,tQcore,PAW%tcore,dum,y)
        PAW%rtVf=Grid%r*PAW%vloc + dum + &
&       (-Pot%nz + Qcore - tQcore)*PAW%hatpot
        !DO l=0,ll
        !      CALL hatL(Grid,PAW,l,PAW%g(:,l+1))
        !ENDDO

    OPEN(1001,file='rvf',form='formatted')
    DO i=1,n
      WRITE(1001,'(1p,30e15.7)') Grid%r(i),PAW%rVf(i),PAW%rtVf(i),&
&       PAW%hatpot(i), PAW%hatden(i),(PAW%g(i,l+1),l=0,ll)
    ENDDO
    CLOSE(1001)

    arg=PAW%tcore+(-Pot%nz + Qcore - tQcore)*PAW%hatden
    CALL poisson(Grid,y,arg,dum,PAW%Eaion)
    write(6,*) 'Eaion  ', y, PAW%Eaion
    arg=dum*PAW%hatden; arg(1)=0; arg(2:n)=arg(2:n)/Grid%r(2:n)
    PAW%Eaionhat=integrator(Grid,arg)
    write(6,*) 'Eaionhat  ', PAW%Eaionhat

    DO ib=1,nbase
       l=PAW%l(ib)
       DO jb=1,nbase
         IF (PAW%l(jb)==l) THEN
           !CALL kinetic_ij(Grid,PAW%ophi(:,ib),PAW%ophi(:,jb),l,x,lr)
           !CALL kinetic_ij(Grid,PAW%otphi(:,ib),PAW%otphi(:,jb),l,y,lr)
            If (scalarrelativistic) then
               call altdtij(Grid,PAW,ib,jb,x)
            Else
              CALL deltakinetic_ij(Grid,PAW%ophi(:,ib),PAW%ophi(:,jb), &
&                  PAW%otphi(:,ib),PAW%otphi(:,jb),l,x,PAW%irc)
            Endif
           WRITE(6,'(" Kinetic ", 3i5, 1p,3e15.7)') ib,jb,l,x
           PAW%Kij(ib,jb)=x
           dum=PAW%ophi(:,ib)*PAW%ophi(:,jb)*PAW%rVf - &
&          PAW%otphi(:,ib)*PAW%otphi(:,jb)*PAW%rtVf
           dum(2:n)=dum(2:n)/Grid%r(2:n)
           dum(1)=0
           PAW%VFij(ib,jb)=integrator(Grid,dum,1,lr)
           WRITE(6,'(" Fix pot ", 3i5, 1p,3e15.7)') ib,jb,l,PAW%VFij(ib,jb)
         ENDIF
       ENDDO
    ENDDO

    ! Average equivalent terms
    DO ib=1,nbase
      DO jb=ib,nbase
        IF(jb>ib) THEN
          x=PAW%Kij(ib,jb); y=PAW%Kij(jb,ib)
          x=0.5d0*(x+y)
          PAW%Kij(ib,jb)=x; PAW%Kij(jb,ib)=x
          x=PAW%VFij(ib,jb); y=PAW%VFij(jb,ib)
          x=0.5d0*(x+y)
          PAW%VFij(ib,jb)=x; PAW%VFij(jb,ib)=x
        ENDIF
      ENDDO
    ENDDO

    DO ib=1,nbase
      DO jb=1,nbase
        llmin=ABS(PAW%l(ib)-PAW%l(jb))
        llmax=PAW%l(ib)+PAW%l(jb)
        DO l=llmin,llmax,2
          DO i=1,irc
            rr=Grid%r(i)
            dum(i)=(rr**l)*(PAW%ophi(i,ib)*PAW%ophi(i,jb) &
&               -PAW%otphi(i,ib)*PAW%otphi(i,jb))
          ENDDO
          PAW%mLij(ib,jb,l+1)=integrator(Grid,dum,1,irc)
          WRITE(6,'("mLij ",3i5,1p,1e15.7)') ib,jb,l,PAW%mLij(ib,jb,l+1)
        ENDDO
      ENDDO
    ENDDO

    ! Average equivalent terms
    DO ib=1,nbase
      DO jb=ib,nbase
        IF (jb>ib) THEN
          llmin=ABS(PAW%l(ib)-PAW%l(jb))
          llmax=PAW%l(ib)+PAW%l(jb)
          DO l=llmin,llmax,2
            x=PAW%mLij(ib,jb,l+1); y=PAW%mLij(jb,ib,l+1)
            x=0.5d0*(x+y)
            PAW%mLij(ib,jb,l+1)=x; PAW%mLij(jb,ib,l+1)=x
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    WRITE(6,*) 'DR matrix elements '
    DO ib=1,nbase
      DO jb=1,nbase
        d=PAW%ophi(:,ib)*PAW%ophi(:,jb)
        td=PAW%otphi(:,ib)*PAW%otphi(:,jb)
        llmin=ABS(PAW%l(ib)-PAW%l(jb))
        llmax=PAW%l(ib)+PAW%l(jb)
        DO l=llmin,llmax,2
          CALL apoisson(Grid,l,lr,d,rh)
          arg=td+PAW%mLij(ib,jb,l+1)*PAW%g(:,l+1)
          CALL apoisson(Grid,l,lr,arg,trh)
          DO kb=1,nbase
            DO lb=1,nbase
              lp=llmax+PAW%l(kb)+PAW%l(lb)
              even=.FALSE.
              IF (2*(lp/2)==lp) even=.TRUE.
              IF (even.AND.l.GE.ABS(PAW%l(kb)-PAW%l(lb)).AND. &
&                 l.LE.(PAW%l(kb)+PAW%l(lb)))  THEN
                arg=PAW%otphi(:,kb)*PAW%otphi(:,lb)+&
&               PAW%mLij(kb,lb,l+1)*PAW%g(:,l+1)
                dum(1:lr)=PAW%ophi(1:lr,kb)*PAW%ophi(1:lr,lb)*rh(1:lr)&
&                 - arg(1:lr)*trh(1:lr)
                dum(1)=0; dum(2:lr)=dum(2:lr)/Grid%r(2:lr)
                PAW%DR(ib,jb,kb,lb,l+1)=integrator(Grid,dum,1,lr)
                !if (abs(PAW%DR(ib,jb,kb,lb,l+1))>100.d0) then
                !   write(6,*) 'problem', ib,jb,kb,lb,l+1,PAW%DR(ib,jb,kb,lb,l+1),&
                !&       PAW%mLij(kb,lb,l+1)
                !   open(unit=1001,file='stuff',form='formatted')
                !   do i=1,lr
                !      write(1001,'(1p,20e15.7)') &
                !&      Grid%r(i),PAW%ophi(i,kb),PAW%ophi(i,lb), &
                !&        PAW%otphi(i,kb),PAW%otphi(i,lb),PAW%g(i,l+1),rh(i),trh(i)
                !   enddo
                !   close(1001)
                !   stop
                ! endif
                WRITE(6,'(5i5,1p,1e20.10)') ib,jb,kb,lb,l, &
&                     PAW%DR(ib,jb,kb,lb,l+1)
              ENDIF
            ENDDO  !lb
          ENDDO  !kb
        ENDDO  !l
      ENDDO  !jb
    ENDDO  !ib

    ! Average equivalent terms
    DO ib=1,nbase
      DO jb=ib,nbase
        llmin=ABS(PAW%l(ib)-PAW%l(jb))
        llmax=PAW%l(ib)+PAW%l(jb)
        DO l=llmin,llmax,2
          DO kb=1,nbase
            DO lb=kb,nbase
              lp=llmax+PAW%l(kb)+PAW%l(lb)
              even=.FALSE.
              IF (2*(lp/2)==lp) even=.TRUE.
              IF (even.AND.l.GE.ABS(PAW%l(kb)-PAW%l(lb)).AND. &
&                 l.LE.(PAW%l(kb)+PAW%l(lb)))  THEN
                arg(1)=PAW%DR(ib,jb,kb,lb,l+1)
                arg(2)=PAW%DR(ib,jb,lb,kb,l+1)
                arg(3)=PAW%DR(jb,ib,kb,lb,l+1)
                arg(4)=PAW%DR(jb,ib,lb,kb,l+1)
                x=0.25d0*(arg(1)+arg(2)+arg(3)+arg(4))
                PAW%DR(ib,jb,kb,lb,l+1)=x
                PAW%DR(ib,jb,lb,kb,l+1)=x
                PAW%DR(jb,ib,kb,lb,l+1)=x
                PAW%DR(jb,ib,lb,kb,l+1)=x
              ENDIF
            ENDDO  !lb
          ENDDO  !kb
        ENDDO  !l
      ENDDO   !jb
    ENDDO    !ib

    If (TRIM(PAW%exctype)=="EXXKLI".OR.TRIM(PAW%exctype)=="HF") &
&     Call COREVAL_EXX(Grid,PAW)

    IF (TRIM(PAW%exctype)=="HF") THEN
      norbit=PAW%TOCCWFN%norbit
      ALLOCATE(PAW%mLic(nbase,norbit,ll+1),PAW%DRC(nbase,nbase,norbit,ll+1),&
&     PAW%mLcc(norbit,norbit,ll+1),PAW%DRCC(nbase,norbit,norbit,ll+1),&
&     PAW%DRCjkl(norbit,nbase,nbase,nbase,ll+1),&
&     PAW%Dcj(norbit,nbase),PAW%TOCCWFN%lqp(norbit,norbit))
      PAW%mLic=0.d0; PAW%DRC=0.d0; PAW%mLcc=0.d0; PAW%DRCC=0.d0; PAW%Dcj=0.d0
      PAW%DRCjkl=0.d0; PAW%TOCCWFN%lqp=0.d0
      lr=MAX(lr,PAW%coretailpoints)
      write(6,*) 'lr for core treatment ', lr
      DO io=1,norbit
        IF (PAW%TOCCWFN%iscore(io)) THEN
          DO jb=1,nbase
            llmin=ABS(PAW%TOCCWFN%l(io)-PAW%l(jb))
            llmax=PAW%TOCCWFN%l(io)+PAW%l(jb)
            DO l=llmin,llmax,2
              DO i=1,n
                rr=Grid%r(i)
                dum(i)=(rr**l)*(PAW%OCCWFN%wfn(i,io)*PAW%ophi(i,jb) &
&                     -PAW%TOCCWFN%wfn(i,io)*PAW%otphi(i,jb))
              ENDDO
              PAW%mLic(jb,io,l+1)=integrator(Grid,dum,1,lr)
              WRITE(6,'("mLic ",3i5,1p,1e15.7)') jb,io,l,PAW%mLic(jb,io,l+1)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      DO io=1,norbit
        IF (PAW%TOCCWFN%iscore(io)) THEN
          DO jo=1,norbit
            IF (PAW%TOCCWFN%iscore(jo)) THEN
              llmin=ABS(PAW%TOCCWFN%l(io)-PAW%TOCCWFN%l(jo))
              llmax=PAW%TOCCWFN%l(io)+PAW%TOCCWFN%l(jo)
              DO l=llmin,llmax,2
                DO i=1,n
                  rr=Grid%r(i)
                  dum(i)=(rr**l)*&
&                    (PAW%OCCWFN%wfn(i,io)*PAW%OCCWFN%wfn(i,jo) &
&                    -PAW%TOCCWFN%wfn(i,io)*PAW%TOCCWFN%wfn(i,jo))
                ENDDO
                PAW%mLcc(jo,io,l+1)=integrator(Grid,dum,1,lr)
                WRITE(6,'("mLcc ",3i5,1p,1e15.7)') jo,io,l,&
&               PAW%mLcc(jo,io,l+1)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      WRITE(6,*) 'DRC matrix elements '
      DO io=1,norbit
        IF (PAW%TOCCWFN%iscore(io)) THEN
          DO ib=1,nbase
            d=PAW%ophi(:,ib)*PAW%OCCWFN%wfn(:,io)
            td=PAW%otphi(:,ib)*PAW%TOCCWFN%wfn(:,io)
            llmin=ABS(PAW%l(ib)-PAW%TOCCWFN%l(io))
            llmax=PAW%l(ib)+PAW%TOCCWFN%l(io)
            DO l=llmin,llmax,2
              CALL apoisson(Grid,l,lr,d,rh)
              arg=td+PAW%mLic(ib,io,l+1)*PAW%g(:,l+1)
              CALL apoisson(Grid,l,n,arg,trh)
              DO jb=1,nbase
                lp=llmax+PAW%l(jb)+PAW%TOCCWFN%l(io)
                even=.FALSE.
                IF (2*(lp/2)==lp) even=.TRUE.
                IF (even.AND.l.GE.ABS(PAW%l(jb)-PAW%TOCCWFN%l(io)).AND. &
&                   l.LE.(PAW%l(jb)+PAW%TOCCWFN%l(io)))  THEN
                  arg=PAW%otphi(:,jb)*PAW%TOCCWFN%wfn(:,io) +&
&                     PAW%mLic(jb,io,l+1)*PAW%g(:,l+1)
                  dum=PAW%ophi(:,jb)*PAW%OCCWFN%wfn(:,io)*rh - arg*trh
                  dum(1)=0; dum(2:n)=dum(2:n)/Grid%r(2:n)
                      PAW%DRC(ib,jb,io,l+1)=integrator(Grid,dum,1,lr)
                  WRITE(6,'(4i5,1p,1e20.10)') ib,jb,io,l, &
&                     PAW%DRC(ib,jb,io,l+1)
                ENDIF
              ENDDO  !jb
            ENDDO  !l
          ENDDO    !ib
        ENDIF
      ENDDO    !io

      ! Average equivalent terms
      DO io=1,norbit
        IF (PAW%TOCCWFN%iscore(io)) THEN
          DO ib=1,nbase
            llmin=ABS(PAW%l(ib)-PAW%TOCCWFN%l(io))
            llmax=PAW%l(ib)+PAW%TOCCWFN%l(io)
            DO l=llmin,llmax,2
              DO jb=ib,nbase
                lp=llmax+PAW%l(jb)+PAW%TOCCWFN%l(io)
                even=.FALSE.
                IF (2*(lp/2)==lp) even=.TRUE.
                IF (even.AND.l.GE.ABS(PAW%l(jb)-PAW%TOCCWFN%l(io)).AND. &
&                   l.LE.(PAW%l(jb)+PAW%TOCCWFN%l(io)))  THEN
                  arg(1)=PAW%DRC(ib,jb,io,l+1)
                  arg(2)=PAW%DRC(jb,ib,io,l+1)
                  x=0.5d0*(arg(1)+arg(2))
                  PAW%DRC(ib,jb,io,l+1)=x
                  PAW%DRC(jb,ib,io,l+1)=x
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      WRITE(6,*) 'DRCC matrix elements '
      DO io=1,norbit
        IF (PAW%TOCCWFN%iscore(io)) THEN
          DO jo=1,norbit
            IF (PAW%TOCCWFN%iscore(jo)) THEN
              d=PAW%OCCWFN%wfn(:,jo)*PAW%OCCWFN%wfn(:,io)
              td=PAW%TOCCWFN%wfn(:,jo)*PAW%TOCCWFN%wfn(:,io)
              llmin=ABS(PAW%TOCCWFN%l(jo)-PAW%TOCCWFN%l(io))
              llmax=PAW%TOCCWFN%l(jo)+PAW%TOCCWFN%l(io)
              DO l=llmin,llmax,2
                CALL apoisson(Grid,l,lr,d,rh)
                arg=td+PAW%mLcc(jo,io,l+1)*PAW%g(:,l+1)
                CALL apoisson(Grid,l,lr,arg,trh)
                DO jb=1,nbase
                  lp=llmax+PAW%l(jb)+PAW%TOCCWFN%l(io)
                  even=.FALSE.
                  IF (2*(lp/2)==lp) even=.TRUE.
                  IF (even.AND.l.GE.ABS(PAW%l(jb)-PAW%TOCCWFN%l(io)) &
&                    .AND.l.LE.(PAW%l(jb)+PAW%TOCCWFN%l(io)))  THEN
                    arg=PAW%otphi(:,jb)*PAW%TOCCWFN%wfn(:,io) +&
&                       PAW%mLic(jb,io,l+1)*PAW%g(:,l+1)
                    dum=PAW%ophi(:,jb)*PAW%OCCWFN%wfn(:,io)*rh - arg*trh
                    dum(1)=0; dum(2:n)=dum(2:n)/Grid%r(2:n)
                    PAW%DRCC(jb,jo,io,l+1)=integrator(Grid,dum,1,lr)
                    WRITE(6,'(4i5,1p,1e20.10)') jb,jo,io,l,PAW%DRCC(jb,jo,io,l+1)
                  ENDIF
                ENDDO  !jb
              ENDDO  !l
            ENDIF
          ENDDO    !jo
        ENDIF
      ENDDO    !io

      WRITE(6,*) 'Dcj matrix elements '
      DO io=1,norbit
        IF (PAW%TOCCWFN%iscore(io)) THEN
          l=PAW%TOCCWFN%l(io)
          DO jb=1,nbase
            IF (PAW%l(jb)==l) THEN
              !CALL kinetic_ij(Grid,PAW%OCCWFN%wfn(:,io),&
              !&           PAW%ophi(:,jb),l,x,lr)
              !CALL kinetic_ij(Grid,PAW%TOCCWFN%wfn(:,io),&
              !&       PAW%otphi(:,jb),l,y,lr)
              CALL deltakinetic_ij(Grid,PAW%OCCWFN%wfn(:,io),PAW%ophi(:,jb), &
&                  PAW%TOCCWFN%wfn(:,io),PAW%otphi(:,jb),l,x,PAW%irc)
              !WRITE(6,'(" Kinetic ", 3i5, 1p,3e15.7)') io,jb,l,x,y,x-y
              PAW%Dcj(io,jb)=x
              dum=PAW%OCCWFN%wfn(:,io)*PAW%ophi(:,jb)*PAW%rVf - &
&             PAW%TOCCWFN%wfn(:,io)*PAW%otphi(:,jb)*PAW%rtVf
              dum(2:n)=dum(2:n)/Grid%r(2:n)
              dum(1)=0
              x=integrator(Grid,dum,1,lr)
              !WRITE(6,'(" Fix pot ", 3i5, 1p,3e15.7)') io,jb,l,x
              PAW%Dcj(io,jb)=PAW%Dcj(io,jb)+x
              WRITE(6,'(3i5,1p,e15.7)') io,jb,l,PAW%Dcj(io,jb)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      WRITE(6,*) 'DRCjkl matrix elements '
      DO io=1,norbit
        IF (PAW%TOCCWFN%iscore(io)) THEN
          DO jb=1,nbase
            d=PAW%OCCWFN%wfn(:,io)*PAW%ophi(:,jb)
            td=PAW%TOCCWFN%wfn(:,io)*PAW%otphi(:,jb)
            llmin=ABS(PAW%TOCCWFN%l(io)-PAW%l(jb))
            llmax=PAW%TOCCWFN%l(io)+PAW%l(jb)
            DO l=llmin,llmax,2
              CALL apoisson(Grid,l,lr,d,rh)
              arg=td+PAW%mLic(jb,io,l+1)*PAW%g(:,l+1)
              CALL apoisson(Grid,l,lr,arg,trh)
              DO kb=1,nbase
                DO lb=1,nbase
                  lp=llmax+PAW%l(kb)+PAW%l(lb)
                  even=.FALSE.
                  IF (2*(lp/2)==lp) even=.TRUE.
                  IF (even.AND.l.GE.ABS(PAW%l(kb)-PAW%l(lb)).AND. &
&                     l.LE.(PAW%l(kb)+PAW%l(lb)))  THEN
                    arg=PAW%otphi(:,kb)*PAW%otphi(:,lb)+&
&                   PAW%mLij(kb,lb,l+1)*PAW%g(:,l+1)
                    dum=PAW%ophi(:,kb)*PAW%ophi(:,lb)*rh - arg*trh
                    dum(1)=0; dum(2:n)=dum(2:n)/Grid%r(2:n)
                    PAW%DRCjkl(io,jb,kb,lb,l+1)=integrator(Grid,dum,1,lr)
                    WRITE(6,'(5i5,1p,1e20.10)') io,jb,kb,lb,l, &
&                   PAW%DRCjkl(io,jb,kb,lb,l+1)
                  ENDIF
                ENDDO  !lb
              ENDDO  !kb
            ENDDO  !l
          ENDDO   !jb
        ENDIF
      ENDDO    !io
    ENDIF

    WRITE(6,*) ' Completed SPMatrixElements'; CALL flush(6)

!   Calculate atomic energy from PAW matrix elements
    PAW%tkin=0; PAW%tion=0; PAW%tvale=0;PAW%txc=0;PAW%Ea=0
    PAW%Etotal=0;PAW%Eaxc=0;PAW%den=0; PAW%tden=0
    norbit=PAW%TOCCWFN%norbit
    DO io=1,norbit
      if (.NOT.PAW%TOCCWFN%iscore(io)) then
         l=PAW%TOCCWFN%l(io) ; occ=PAW%TOCCWFN%occ(io)
         CALL kinetic(Grid,PAW%TOCCWFN%wfn(:,io),l,x)
         PAW%tkin=PAW%tkin+occ*x
         PAW%den=PAW%den+occ*(PAW%OCCWFN%wfn(:,io))**2
         PAW%tden=PAW%tden+occ*(PAW%TOCCWFN%wfn(:,io))**2
      endif   
    ENDDO
    write(6,*) 'smooth kinetic ', PAW%tkin

    arg=PAW%den+FC%coreden-PAW%tden-PAW%tcore
    x=-Pot%nz+integrator(Grid,arg,1,irc)
    write(6,*) 'q00 for atom ', x
    arg=PAW%tden+PAW%tcore+x*PAW%hatden
    write(6,*) ' Total charge check ', integrator(Grid,arg)
    call poisson(Grid,x,arg,dum,y)
    write(6,*) ' Smooth coulomb ', y
    PAW%tvale=PAW%tkin+y
    arg=PAW%vloc*PAW%tden
    PAW%tion=integrator(Grid,arg,1,irc)
    write(6,*) ' Vloc energy ', PAW%tion
    PAW%tvale=PAW%tvale+PAW%tion

    IF (TRIM(PAW%exctype)=="HF".or.TRIM(PAW%exctype)=="EXX".or.&
&       TRIM(PAW%exctype)=="EXXKLI") THEN
      !CALL Get_Energy_EXX(Grid,PAW%TOCCWFN,x) ! not correct need moment
      CALL Get_Energy_EXX_pseudo(Grid,PAW,x)
      write(6,*) 'Warning: does not include core contributions'
    ELSE
      arg=PAW%tden+PAW%tcore
      CALL exch(Grid,arg,dum,y,x)
    ENDIF
    write(6,*) 'Smooth exchange-correlation contribution ', x
    PAW%txc=x   ; PAW%tvale=PAW%tvale+PAW%txc

!   one center terms
    arg=PAW%den-PAW%tden
    x=integrator(Grid,arg,1,irc)
    write(6,*) 'valence Q', x
    PAW%Ea=-PAW%Eaion-x*PAW%Eaionhat

    wij=0
    do io=1,norbit
      l=PAW%TOCCWFN%l(io); occ=PAW%TOCCWFN%occ(io)
      if (occ>1.d-8.and..NOT.PAW%TOCCWFN%iscore(io)) &
&      call calcwij(Grid,PAW,l,occ,PAW%TOCCWFN%wfn(:,io),wij)      
    enddo

    
    do ib=1,nbase
      do jb=1,nbase
        PAW%wij(lb,jb)=wij(ib,jb)
      enddo
    enddo
    DO ib=1,nbase
      do jb=1,nbase
        if (PAW%l(ib)==PAW%l(jb)) then
          PAW%Ea=PAW%Ea+PAW%wij(ib,jb)*(PAW%Kij(ib,jb)+PAW%VFij(ib,jb))
          x=0
          do kb=1,nbase
            do lb=1,nbase
              if (PAW%l(kb)==PAW%l(lb)) then
                x=x+PAW%wij(kb,lb)*PAW%DR(ib,jb,kb,lb,1)
              endif
            enddo
          enddo
          PAW%Ea=PAW%Ea+0.5d0*x*PAW%wij(ib,jb)
        endif
      enddo
    enddo

    write(6,*) 'Ea up to exchange-correlation terms ', PAW%Ea

    IF (TRIM(PAW%exctype)=="HF".or.TRIM(PAW%exctype)=="EXX".or.&
&       TRIM(PAW%exctype)=="EXXKLI") THEN
      CALL Get_Energy_EXX_pseudo_one(Grid,PAW,x)
      write(6,*) 'one-center HF exchange', x
      PAW%Ea=PAW%Ea+x; PAW%Eaxc=x
    ELSE
      arg=PAW%tden+PAW%tcore
      CALL exch(Grid,arg,dum,y,x,irc)
      arg=PAW%den+FC%coreden
      CALL exch(Grid,arg,dum,y,z,irc)
      write(6,*) ' one center xc ', z,x,z-x
      PAW%Ea=PAW%Ea+z-x; PAW%Eaxc=z-x
    ENDIF

    PAW%Etotal=PAW%tvale+PAW%Ea
    write(6,*) 'Energy terms ', PAW%tvale, PAW%Ea, PAW%Etotal

    PAW%mesh_size=PAW%irc+Grid%ishift
    PAW%coretailpoints=MAX(PAW%coretailpoints,PAW%mesh_size)


    DEALLOCATE(arg,dum,rh,trh,d,td,wij)

  END SUBROUTINE SPMatrixElements


  SUBROUTINE COREVAL_EXX(Grid,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER :: i,k, io,ib,ic,nbase,li,lj,l
      REAL(8) :: accum,term
      REAL(8), ALLOCATABLE :: f(:)

      nbase=PAW%nbase
      k=0            ! core-valence terms
      do io=1,PAW%OCCWFN%norbit
        if (PAW%OCCWFN%iscore(io)) then
          li=PAW%OCCWFN%l(io)
          do ib=1,PAW%nbase
            do ic=1,PAW%nbase
              if (PAW%l(ib)==PAW%l(ic)) k=k+1
            enddo
          enddo
        endif
      enddo
      PAW%ncoreshell=k
      if(k>0) then
        ALLOCATE (PAW%DRVC(k,nbase,nbase),f(k)); PAW%DRVC=0
        !!!!WRITE(ifatompaw,'("  COREVAL_LIST   ",i10)') k
        !!!!WRITE(ifatompaw,'("  COREVAL_R      ")')
        i=0;f=0
        do io=1,PAW%OCCWFN%norbit
          if (PAW%OCCWFN%iscore(io)) then
            i=i+1
            li=PAW%OCCWFN%l(io)
            do ib=1,PAW%nbase
              lj=PAW%l(ib)
              do ic=1,PAW%nbase
                if (PAW%l(ib)==PAW%l(ic)) then
                  f(i)=0
                  do l=abs(li-lj),(li+lj),2
                    call EXXwgt(1.d0,1.d0,1,li,2,lj,l,accum)
                    call CondonShortley(Grid,l,PAW%OCCWFN%wfn(:,io), &
&                          PAW%ophi(:,ib),PAW%OCCWFN%wfn(:,io), &
&                          PAW%ophi(:,ic),term)
                    f(i)=f(i)+2*accum*term !EXXwgt returns 1/2*3J
                    write(6,*) 'core-val CondonShortley',&
&                              i,li,lj,l,2*accum,term
                  enddo
                  !!!!!WRITE(ifatompaw,'(3i10,1p,1e25.17)') i, ib,ic,f(i)
                  PAW%DRVC(i,ib,ic)=f(i)
                endif
              enddo   !ic
            enddo   !ib
          endif
        enddo   !norbit

        DEALLOCATE(f)
      endif

   END SUBROUTINE COREVAL_EXX

END MODULE pseudo_sub
