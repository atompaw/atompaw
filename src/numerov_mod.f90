MODULE Numerov_mod
  USE gridmod
  USE blockdavidson_mod

  IMPLICIT NONE

  REAL(8), ALLOCATABLE, PRIVATE :: hn(:),hd(:),on(:),od(:)
  INTEGER, PRIVATE :: firsttime=0
  REAL(8), PRIVATE :: qq
  REAL(8), PRIVATE, PARAMETER :: tailparm=12.d0,small=1.d-5
  TYPE(GridInfo), PRIVATE, POINTER :: GridPtr

CONTAINS

  SUBROUTINE BoundNumerov(Grid,rv,v0,v0p,nz,l,nroot,Eig,Psi,BDsolve,success)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: rv(:),v0,v0p
    INTEGER, INTENT(IN) :: nz,l,nroot
    REAL(8), INTENT(INOUT) :: Eig(:), Psi(:,:)
    LOGICAL, INTENT(IN) :: BDsolve
    LOGICAL, INTENT(INOUT) :: success

    INTEGER, PARAMETER :: repeat=4
    INTEGER :: i,j,k,n,many,count
    REAL(8), ALLOCATABLE :: vec(:,:),f(:),dum(:),e(:)
    REAL(8) :: x

  If (BDsolve) then  
    success=.false.
    CALL initBoundNumerov(Grid)
    CALL startBoundNumerov(Grid,rv,l)

    n=Grid%n; k=n-2 ; many=4*nroot
    ALLOCATE(vec(k,many),f(k),dum(k),e(many),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in BoundNumerov ', i,n,k
       STOP
    ENDIF

    vec=0
    count=0
    DO i=1,nroot
       f=0
       f(1:n-2)=Psi(2:n-1,i)
       CALL gramschmidt(Grid,count,vec,f)
       count=count+1
    ENDDO
    !WRITE(6,*) 'Starting with ',count,'  orthonormal functions '
    IF (count<1) THEN
       WRITE(6,*) 'Error in initial wavefunctions -- stopping '
       STOP
    ENDIF

    i=1
    DO
       IF (count>=many) EXIT
       IF (i>count) THEN
          WRITE(6,*) 'Error in initial wavefunctions -- i>count ', i,count
          STOP
       ENDIF
       dum=vec(:,i)
       IF (usingloggrid(Grid)) dum(1:n-2)=dum(1:n-2)/Grid%pref(2:n-1)
       CALL tridiagmult(hn,hd,dum,f)
       IF (usingloggrid(Grid)) f(1:n-2)=f(1:n-2)*Grid%pref(2:n-1)
       CALL gramschmidt(Grid,count,vec,f)
       i=i+1; count=count+1
    ENDDO

    IF (usingloggrid(Grid)) THEN
       DO j=1,many
          vec(1:n-2,j)=vec(1:n-2,j)/Grid%pref(2:n-1)
       ENDDO
    ENDIF

    CALL BlockDavidson(nroot,many,vec,e,success,repeat,hvov,multres)

    DO j=1,nroot
       If (e(j)<-1.d-5) then
          Eig(j)=e(j)
       elseif (e(j)>=-1.d-5.and.j==1) then
          Eig(j)=-REAL(nz**2)/(j+l)**2   
       elseif (e(j)>=-1.d-5.and.j>1) then   
          Eig(j)=(Eig(j-1)*(nroot-j+1))/nroot
       endif   
       Psi(1,j)=0;  Psi(n,j)=0
       !WRITE(6,*) 'BD Eigenvalue -- ', j, Eig(j)
       !call flush(6)
       IF (usingloggrid(Grid)) THEN
          DO i=2,n-1
             Psi(i,j)=vec(i-1,j)*Grid%pref(i)
          ENDDO
       ELSE
          DO i=2,n-1
             Psi(i,j)=vec(i-1,j)
          ENDDO
       ENDIF
       x=Integrator(Grid,Psi(:,j)*Psi(:,j))
       !WRITE(6,*) 'normalization = ', x
       Psi(:,j)=Psi(:,j)/SQRT(x)
    ENDDO

                                                                                    DEALLOCATE(vec,f,dum,e)
    CALL endBoundNumerov

    If (.not.success) then
       write(6,*) 'program faltering due to failure of BlockDavidson'
       write(6,*) 'l = ',l
       write(6,*) 'Eig = ', Eig(1:nroot)

    endif
  Endif     !BDsolve

    write(6,*) 'Before newboundsch',l,nroot, Eig(1:nroot); call flush(6)
    CALL newboundsch(Grid,rv,v0,v0p,l,nroot,Eig,Psi,success)
    write(6,*) 'After newboundsch',l,nroot, Eig(1:nroot); call flush(6)

    ! adjust sign
    Do j=1,nroot
       if (Psi(3,j)<0.d0) Psi(:,j)=-Psi(:,j)
    Enddo

  END SUBROUTINE BoundNumerov

  SUBROUTINE initBoundNumerov(Grid)
    TYPE(GridInfo), TARGET, INTENT(IN) :: Grid

    INTEGER :: n,i

    IF(firsttime==0) THEN
       n=Grid%n
       ALLOCATE(hn(n),hd(n),on(n),od(n),stat=i)
       IF (i/=0) THEN
          WRITE(6,*) 'Allocation error in initBoundNumerov', i,n
          STOP
       ENDIF
       firsttime=1
       GridPtr=>Grid
    ENDIF
  END SUBROUTINE initBoundNumerov

  SUBROUTINE endBoundNumerov
    DEALLOCATE(hn,hd,on,od)
    firsttime=0
  END SUBROUTINE endBoundNumerov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   fin(1:n-2), fout(1:n-2) map to grid points 2:n-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE tridiagmult(nn,dd,fin,fout)
    REAL(8), INTENT(IN) :: nn(:),dd(:),fin(:)
    REAL(8), INTENT(OUT) :: fout(:)

    INTEGER :: i,j,k,n

    n=GridPtr%n
    fout(1:n-2)=dd(2:n-1)*fin(1:n-2)

    DO i=1,n-3
       fout(i)=fout(i)+nn(i+2)*fin(i+1)
       fout(i+1)=fout(i+1)+nn(i+1)*fin(i)
    ENDDO

  END SUBROUTINE tridiagmult

  SUBROUTINE startBoundNumerov(Grid,rv,l)
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: rv(:)
    INTEGER, INTENT(IN) :: l

    INTEGER:: i,j,k,n
    REAL(8):: h,x,h2,angm

    angm=l*(l+1);h2=(Grid%h**2);n=Grid%n
    x=0
    IF (l==0) x=rv(1)/Grid%r(2)
    IF (l==1) x=2.d0/(Grid%r(2)**2)

    hd=0.d0 ; od=1.d0
    DO i=2,n
       hd(i)=rv(i)/Grid%r(i)+angm/(Grid%r(i)**2)
    ENDDO


    IF (usingloggrid(Grid)) THEN
       hd=Grid%rr02*hd +0.25d0
       x=Grid%rr02(1)*x
       od=Grid%rr02*od
    ENDIF
    hn=-1.2d0+0.1d0*h2*hd
    hd=2.4d0+h2*hd
    hd(2)=hd(2)+0.1d0*h2*x
    on=0.1d0*od
    od=od
    hd=hd/h2
    hn=hn/h2

    qq=-rv(n)/2
    IF (qq<small) qq=0.d0

  END SUBROUTINE startBoundNumerov

  SUBROUTINE hvov(vin,hv,ov,ee)
    REAL(8), INTENT(IN) :: vin(:)
    REAL(8), INTENT(OUT) :: hv(:),ov(:),ee

    INTEGER :: n,i
    REAL(8) :: x

    hv=0.d0; ov=0.d0
    CALL tridiagmult(hn,hd,vin,hv)
    CALL tridiagmult(on,od,vin,ov)

    x=DOT_PRODUCT(vin,ov)
    IF (ABS(x)<1.d-9) x=1
    ee=DOT_PRODUCT(vin,hv)/x

  END SUBROUTINE hvov

  FUNCTION multres(v1,v2)
    REAL(8) :: multres
    REAL(8), INTENT(IN) :: v1(:),v2(:)

    INTEGER :: n,i

    multres=DOT_PRODUCT(v1,v2)

  END FUNCTION multres

  !******************************************************************
  !  SUBROUTINE newboundsch(Grid,rv,v0,v0p,l,nroot,Eig,Psi,ok)
  !******************************************************************
  SUBROUTINE newboundsch(Grid,rv,v0,v0p,l,nroot,Eig,Psi,ok)
    !  pgm to solve radial schroedinger equation for nroot bound state
    !    energies and wavefunctions for angular momentum l
    !    with potential rv/r
    !
    !   Assymptotic from of wfn:
    !     Psi(r) = const*(r**q/kappa)*EXP(-kappa*r), where eig=-kappa**2
    !  uses Noumerov algorithm
    !
    !  For l=0,1 corrections are needed to approximate wfn(r=0)
    !     These depend upon:
    !         e0 (current guess of energy eigenvalue)
    !         l,nz
    !         v(0) == v0 electronic potential at r=0
    !         v'(0) == v0p derivative of electronic potential at r=0
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
    !
    TYPE(GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: rv(:),v0,v0p
    INTEGER, INTENT(IN) :: l,nroot
    REAL(8), INTENT(INOUT) :: Eig(:), Psi(:,:)
    LOGICAL, INTENT(OUT) :: ok

    REAL(8), PARAMETER :: convre=1.d-10,vlrg=1.d30
    INTEGER, PARAMETER :: niter=1000

    REAL(8), ALLOCATABLE :: p1(:),p2(:),dd(:)
    INTEGER :: nz,n,ierr
    REAL(8) :: h,q
    REAL(8) :: err,convrez,energy,zeroval,zz
    REAL(8) :: scale,emin,emax,best,rout,ppp
    REAL(8) :: arg,r,r2,veff,pppp1,rin,dele,x,rvp1,pnp1,bnp1
    INTEGER :: iter,i,j,k,node,match,mxroot,ntroot,ir,iroot
    INTEGER :: least,many,ifac

    ifac=0
    n=Grid%n
    h=Grid%h
    ALLOCATE(p1(n),p2(n),dd(n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) ' Error in newboundsch allocation ',i,n
       STOP
    ENDIF

    nz=-(rv(1)-0.1d0)/2;zz=nz
    q=-rv(n)/2
    IF (q<0.001d0) q=0.d0
    err=n*nz*(h**4)
    convrez=convre
    IF (nz.GT.0) convrez=convre*nz
    !     write(6,*) 'expected error = ',err
    ierr=0

    !
    emin=MAX(-REAL((nz)**2)/(l+1)**2-0.5d0,Eig(1)-1.d0)
    !write(6,*) 'in boundsch --' , emin, nz, l
    emax=0.d0

    !WRITE(6,'("starting boundsch with eigenvalues -- ",1p,20e15.7)') &
    !&    Eig(1:nroot)

    DO iroot=1,nroot
       best=1.d10; dele=1.d10
       energy=Eig(iroot)
       !write(6,*) 'starting iroot energy' ,iroot, energy
       IF (energy.LT.emin) energy=emin
       IF (energy.GT.emax) energy=emax
       ok=.FALSE.
       !write(6,*) 'iter,iroot,energy',iter,iroot,energy
       !write(6,*) 'emin,max',emin,emax
       BigIter: DO iter=1,niter
          !write(6,*) 'In iter with energy', iter,energy,niter,l,iroot
          !  start inward integration
          !  start integration at n
          ! find classical turning point
          CALL ClassicalTurningPoint(Grid,rv,l,energy,match)
          match=MAX(5,match)
          match=MIN(n-15,match)
          ppp=SQRT(ABS(-energy))
          p2=0
          !p2(n)=(Grid%r(n)**(qq/ppp))
          !p2(n-1)=(Grid%r(n-1)**(qq/ppp))*EXP(-ppp*(Grid%r(n-1)-Grid%r(n)))
          p2(n)=wfnend(l,energy,Grid%r(n),Grid%r(n))
          p2(n-1)=wfnend(l,energy,Grid%r(n-1),Grid%r(n))
          CALL backward_numerov(Grid,l,match,energy,rv,p2)
          match=match+6
          CALL derivative(Grid,p2,dd,match-5,match+5)
          rin=dd(match)/p2(match)
          ! write(6,*) ' match point = ',match,rin,p2(match)
          !  start outward integration
          !    correct behavior near r=0
          ! initialize p1
          p1=0
          p1(2)=wfninit(zz,l,v0,v0p,energy,Grid%r(2))
          zeroval=0
          IF (l==0) zeroval=-2*nz
          IF (l==1) zeroval=2

          CALL forward_numerov(Grid,l,match+6,energy,rv,zeroval,p1,node)

          CALL derivative(Grid,p1,dd,match-5,match+5)
          rout=dd(match)/p1(match)
          !write(6,*) 'node,match,rin,rout',node,(iroot-1),match,rin,rout
          ! check whether node = (iroot-1)
          !   not enough nodes -- raise energy
          IF (node.LT.iroot-1) THEN
             emin=MAX(emin,energy)-1.d-5
             energy=emax-(emax-energy)*ranx()
             ifac=9
             !   too many nodes -- lower energy
          ELSEIF (node.GT.iroot-1) THEN
             IF (energy.LT.emin) THEN
                ierr=ierr+9*(10**(iroot-1))
                WRITE(6,*) 'newboundsch error -- emin too high',l,nz,emin,energy
                RETURN
             ENDIF
             emax=MIN(emax,energy+1.d-5)
             energy=emin+(energy-emin)*ranx()
             !write(6,*) 'energy reset to ', energy
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
                Psi(1:n,iroot)=p1(1:n)
                Eig(iroot)=energy
                !write(6,*) 'root',l,iroot,Eig(iroot),emin,emax
                best=x
             ENDIF
             IF (ABS(dele).LE.convrez) THEN
                !WRITE(6,*) 'iter with dele' , iter,dele
                ok=.TRUE.
                !  eigenvalue found
                ierr=ierr+10**(iroot-1)
                IF (iroot+1.LE.nroot) THEN
                   emin=energy+1.d-5
                   emax=0
                   energy=(emin+emax)/2
                   IF (energy.LT.emin) energy=emin
                   IF (energy.GT.emax) energy=emax
                   best=1.d10
                ENDIF
                EXIT BigIter
             ENDIF
             IF (ABS(dele).GT.convrez) THEN
                !write(6,*) 'iter with dele' , iter,dele
                energy=energy+dele
                ! if energy is out of range, pick random energy in correct range
                IF (emin-energy.GT.convrez.OR.energy-emax.GT.convrez)         &
&                    energy=emin+(emax-emin)*ranx()
                ifac=2
                !write(6,*) 'continuing with iter dele', iter,dele
             ENDIF
          ENDIF
       ENDDO BigIter !iter
       IF (.NOT.ok) THEN
          ierr=ierr+ifac*(10**(iroot-1))
          WRITE(6,*) 'no convergence in newboundsch',iroot,l,dele,energy
          WRITE(6,*) ' best guess of eig, dele = ',Eig(iroot),best
          IF (iroot.LT.nroot) THEN
             DO ir=iroot+1,nroot
                ierr=ierr+9*(10**(ir-1))
             ENDDO
          ENDIF
          ! reset wfn with hydrogenic form
          j=iroot+l+1
          Psi(:,iroot)=0
          ppp=(j)*SQRT(ABS(Eig(iroot)))
          DO i=2,n
             Psi(i,iroot)=hwfn(ppp,j,l,Grid%r(i))
          ENDDO
       ENDIF
    ENDDO !iroot

    !WRITE(6,'("finish boundsch with eigenvalues -- ",1p,20e15.7)') &
    !&    Eig(1:nroot)
    DEALLOCATE(p1,p2,dd)
    !WRITE(6,*) 'returning from newboundsch -- ierr=',ierr

  END SUBROUTINE newboundsch

  !*******************************************************************
  !  FUNCTION wfninit(nz,l,v0,v0p,energy,r)
  !*******************************************************************
  FUNCTION wfninit(nz,l,v0,v0p,energy,r)
    ! returns the solution of the Schroedinger equation near r=0
    !  using power series expansion
    REAL(8) :: wfninit
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: nz,v0,v0p,energy,r

    REAL(8) :: c1,c2,c3

    c1=-REAL(nz)/(l+1.d0)
    c2=((v0-energy)-2*nz*c1)/(4*l+6.d0)
    c3=(v0p+(v0-energy)*c1-2*nz*c2)/(6*l+12.d0)

    wfninit=(r**(l+1))*(1+r*(c1+r*(c2+r*c3)))

    !write(6,*) 'In wfninit', r,v0,v0p, energy,wfninit

  END FUNCTION wfninit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FUNCTION wfnend(l,energy,r,rN)
  !
  !  Find asymptotic form of wavefunction
  !    assuming equations has form
  !   (- d^2  +l(l+1) -2q   +b^2 )
  !   (  ---   ------ ---        )  P(r) = 0
  !   (  dr^2    r^2    r        )
  !       where b^2=-energy
  !
  !        P(r) = exp(-b*(r-rN))*r^(q/b)(1+(l*(l+1)-q/b)*(q/b-1)/(2*b)/r + ...)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION wfnend(l,energy,r,rN)
    REAL(8) :: wfnend
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: energy,r,rN

    REAL(8) :: qbb,b,wfn,cn,term,fac
    INTEGER :: i
    INTEGER, PARAMETER :: last=5

    IF (energy>=0.d0) THEN
       wfnend=0
       RETURN
    ENDIF

    b=SQRT(-energy)
    qbb=qq/b
    !write(6,*) ' qbb = ', qbb
    cn=l*(l+1)
    fac=EXP(-b*(r-rN))*(r**qbb)
    term=1.d0;   wfn=0
    DO i=1,last
       wfn=wfn+term
       IF (i<last) THEN
          term=-term*((qbb-i+1)*(qbb-i)-cn)/(2*b*i)/r
       ENDIF

    ENDDO
    wfnend=fac*wfn
  END FUNCTION wfnend

  !**********************************************************************
  !      subroutine unboundsch(Grid,rv,v0,v0p,nr,l,energy,wfn,nodes)
  !  pgm to solve radial schroedinger equation for unbound states
  !    at energy 'energy' and at angular momentum l
  !
  !    with potential rv/r, given in uniform mesh of n points
  !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
  !                               p((n+1)*h)=0
  !  nz=nuclear charge
  !
  !  uses Noumerov algorithm
  !
  !  For l=0,1 corrections are needed to approximate wfn(r=0)
  !     These depend upon:
  !         e0 (current guess of energy eigenvalue)
  !         l,nz
  !         v(0) == v0 electronic potential at r=0
  !         v'(0) == v0p derivative of electronic potential at r=0
  !
  ! also returns node == number of nodes for calculated state
  !************************************************************************
  SUBROUTINE unboundsch(Grid,rv,v0,v0p,nr,l,energy,wfn,nodes)
    TYPE(GridInfo), INTENT(IN)  :: Grid
    REAL(8), INTENT(IN) :: rv(:),v0,v0p
    INTEGER, INTENT(IN) :: nr,l
    REAL(8), INTENT(IN) :: energy
    REAL(8), INTENT(INOUT) :: wfn(:)
    INTEGER, INTENT(INOUT) :: nodes

    INTEGER :: n,nz,i,j,k,ierr
    REAL(8) :: zeroval,scale,zz

    n=Grid%n
    IF (nr > n) THEN
       WRITE(6,*) 'Error in unboundsch -- nr > n', nr,n
       STOP
    ENDIF
    nz=-(rv(1)-0.1d0)/2;zz=nz

    ! initialize wfn
    wfn=0
    wfn(2)=wfninit(zz,l,v0,v0p,energy,Grid%r(2))
    zeroval=0
    if (l==0) zeroval=-2*nz
    if (l==1) zeroval=2

    call forward_numerov(Grid,l,nr,energy,rv,zeroval,wfn,nodes)
    !
    ! normalize to unity within integration range
    !
    scale=1.d0/overlap(Grid,wfn(1:nr),wfn(1:nr),1,nr)
    scale=SIGN(SQRT(scale),wfn(nr-2))
    wfn(1:nr)=wfn(1:nr)*scale

  END SUBROUTINE unboundsch

END MODULE Numerov_mod
