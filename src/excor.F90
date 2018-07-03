!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!     initexch, Report_EXC, pwldafunc, pbefunc, radialexcpbe, exch
!  This module contains the following active functions:
!     Logofterm, Aofec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE excor

  USE Tools
  USE globalmath
  USE gridmod
  USE libxc_mod
  USE atomdata, only : scalarrelativistic,diracrelativistic

  IMPLICIT NONE

  CHARACTER(132) :: exctype
  INTEGER, PRIVATE, SAVE :: itype
  INTEGER, PRIVATE, PARAMETER :: LDA_PW=14
  INTEGER, PRIVATE, PARAMETER :: GGA_PBE=16
  INTEGER, PRIVATE, PARAMETER :: GGA_PBESOL=18
  INTEGER, PRIVATE, PARAMETER :: LIBXC=-1

  ! Parameters for the Perdew-Wang (PRB 45,13244 (1992)) LDA correlation
  REAL(8), PRIVATE, PARAMETER :: AA=0.0310907d0
  REAL(8), PRIVATE, PARAMETER :: a1=0.21370d0
  REAL(8), PRIVATE, PARAMETER :: b1=7.59570d0
  REAL(8), PRIVATE, PARAMETER :: b2=3.58760d0
  REAL(8), PRIVATE, PARAMETER :: b3=1.63820d0
  REAL(8), PRIVATE, PARAMETER :: b4=0.49294d0

  ! Parameters from PBE paper (and in program from web)
  !  Perdew, Burke, and Ernzerhof (PRL 77, 3865 (1996))
  !  http://chem.ps.uci.edu/~kieron/dft/pubs/PBE.asc
  REAL(8), PRIVATE, PARAMETER :: kappa= 0.804d0
  REAL(8), PRIVATE, PARAMETER :: muorig = 0.2195149727645171d0
  REAL(8), PRIVATE, PARAMETER :: betorig = 0.06672455060314922d0
  REAL(8), PRIVATE, PARAMETER :: gamm = 0.03109069086965489503494086371273d0

  ! Parameters from new PBEsol
  !      (arXiv:0711.0156v2 [cond-mat.mtrl-sci] 22 Feb 2008)
  !   http://chem.ps.uci.edu/~kieron/dft/pubs/PBEsol.html
  REAL(8), PRIVATE, PARAMETER :: musol = 0.123456790123456d0
  REAL(8), PRIVATE, PARAMETER :: betsol = 0.046d0

  REAL(8), PRIVATE, SAVE :: mu,beta,betabygamm

CONTAINS

  !***********************************************************
  !* Initialization
  !***********************************************************

  SUBROUTINE initexch
    integer :: id(2)=(/0,0/)
    !  choose form of exchange-correlation potential
    CALL Uppercase(exctype)
    WRITE(6,*) exctype; call flush_unit(6)

    SELECT CASE(TRIM(exctype))
    CASE default
      IF (have_libxc) THEN
        call libxc_getid_fromInput(exctype,id)
        call libxc_init_func(id,1)
        itype = LIBXC
        WRITE(6,*) 'Using Libxc -- ', TRIM(exctype) ; call flush_unit(6)
      ELSE
        WRITE(6,*) 'ERROR: cannot understand Exchange-Correlation functional from input file!'
        WRITE(6,*) '       possible issue: atompaw not compiled with libXC.'
        stop
        !itype = LDA_PW
        !exctype='LDA-PW'
        !WRITE(6,*) 'Perdew-Wang correlation'
      END IF
    CASE('LDA-PW')
       itype = LDA_PW
       WRITE(6,*) 'Perdew-Wang correlation'
    CASE('GGA-PBE')
       itype = GGA_PBE
       WRITE(6,*) 'Perdew-Burke-Ernzerhof GGA'
       mu=muorig
       beta=betorig
       betabygamm=beta/gamm
    CASE('GGA-PBESOL')
       itype = GGA_PBESOL
       WRITE(6,*) 'Perdew-Burke-Ernzerhof modified (PBEsol) GGA'
       mu=musol
       beta=betsol
       betabygamm=beta/gamm
    END SELECT
    WRITE(6,*)

  END SUBROUTINE initexch

  SUBROUTINE Report_EXC(if)
     INTEGER, INTENT(IN) :: if

     SELECT CASE(itype)
     CASE default
        write(if,*) 'Incorrect exchange-correlation set up'
     CASE(LDA_PW)
        write(if,*) 'Perdew-Wang LDA -- PRB 45, 13244 (1992)'
     CASE(GGA_PBE)
        write(if,*) 'Perdew-Burke-Ernzerhof GGA -- PRL 77, 3865 (1996)'
     CASE(GGA_PBESOL)
        write(if,*) 'Perdew-Burke-Ernzerhof modified (PBEsol) GGA'
        write(if,*) 'http://chem.ps.uci.edu/~kieron/dft/pubs/PBEsol.html'
     CASE(LIBXC)   
        write(if,*)  'Using Libxc -- ',TRIM(exctype)   
     END SELECT

  END SUBROUTINE Report_EXC

  !**********************************************************
  ! Subroutine to calculate the LDA exchange correlation functionals
  !   using the form of Perdew and Wang (PRB 45, 13244 (1992)
  !   assuming no spin polarization
  !  Inside this routine, energies are in Hartree units

  SUBROUTINE pwldafunc(den,exc,vxc)
    REAL(8), INTENT(IN) :: den    !density
    REAL(8), INTENT(OUT) :: exc,vxc

    ! Variables depending on den
    REAL(8) :: n,g,kf,rs,ks

    REAL(8) :: ex,ec,pprs,decdrs
    REAL(8) :: term

    n=den
    IF (n < machine_zero)  THEN
       exc=0.d0; vxc=0.d0
       RETURN
    ENDIF

    kf=(3.d0*(pi**2)*n)**0.3333333333333333333333333333d0
    rs=(3.d0/(4.d0*pi*n))**0.3333333333333333333333333333d0
    ks=SQRT(4.d0*kf/pi)

    ex=-3.d0*kf/(4.d0*pi)
    pprs=SQRT(rs)*(b1+b3*rs)+rs*(b2+b4*rs)
    !ec=-2.d0*AA*(1.d0+a1*rs)*ddlog(1.d0+1.d0/(2.d0*AA*pprs))
    term=Logofterm(1.d0/(2.d0*AA*pprs))
    ec=-2.d0*AA*(1.d0+a1*rs)*term

    exc=ex+ec

    !decdrs=-(2.d0*AA*a1)*ddlog(1.d0+1.d0/(2*AA*pprs)) &
    decdrs=-(2.d0*AA*a1)*term &
&        +((1.d0+a1*rs)*((b1+3*b3*rs)/(2.d0*SQRT(rs))+&
&        b2+2*b4*rs))/(pprs*(pprs+1.d0/(2.d0*AA)))

    vxc = (4.d0/3.d0)*ex+ec-(decdrs*rs)/3.d0

    IF ((ABS(exc).GT.1.d65).OR.(ABS(vxc).GT.1.d65)) THEN
       WRITE(6,*) 'Problem in PW',n,rs,ec
    ENDIF

    exc=2*exc; vxc=2*vxc      ! change to Rydberg units

    RETURN
  END SUBROUTINE pwldafunc
  !**********************************************************
  ! Subroutine to calculate the exchange correlation functionals
  !   using the form of Perdew, Burke, and Ernzerhof (PRL 77, 3865 (1996))
  !   assuming no spin polarization
  !  Inside this routine, energies are in Hartree units

  SUBROUTINE pbefunc(den,grad,fxc,dfxcdn,dfxcdgbg)
    REAL(8), INTENT(IN) :: den,grad    !density, magnitude of grad(density)
    REAL(8), INTENT(OUT) :: fxc,dfxcdn,dfxcdgbg

    ! Variables depending on den & grad
    REAL(8) :: n,g,kf,rs,ks,s,t

    REAL(8) :: ex,ec,Fx,H,A,pprs,ppt,At2,dFds,dHdt,decdrs,dHdrs,dHdA,dAdrs
    REAL(8) :: term,dHdtbg,dFdsbg

    n=den
    IF (n < machine_zero)  THEN
       fxc=0.d0; dfxcdn=0.d0; dfxcdgbg=0.d0
       RETURN
    ENDIF
    g=grad
    IF (g < machine_zero) g=machine_zero

    kf=(3.d0*(pi**2)*n)**0.3333333333333333333333333333d0
    rs=(3.d0/(4.d0*pi*n))**0.3333333333333333333333333333d0
    ks=SQRT(4.d0*kf/pi)
    s=g/(2.d0*kf*n)
    t=g/(2.d0*ks*n)
    IF (s*s > machine_infinity .or. t*t > machine_infinity)  THEN
       fxc=0.d0; dfxcdn=0.d0; dfxcdgbg=0.d0
       RETURN
    ENDIF

    ex=-3.d0*kf/(4.d0*pi)
    pprs=SQRT(rs)*(b1+b3*rs)+rs*(b2+b4*rs)
    !ec=-2.d0*AA*(1.d0+a1*rs)*ddlog(1.d0+1.d0/(2.d0*AA*pprs))
    term=Logofterm(1.d0/(2.d0*AA*pprs))
    ec=-2.d0*AA*(1.d0+a1*rs)*term
    Fx=1.d0+kappa -kappa/(1.d0+(mu/kappa)*s*s)
    A=Aofec(ec)
    At2=A*t*t
    ppt=(1.d0+At2*(1.d0+At2))
    !H=gamm*ddlog(1.d0+(betabygamm)*(t*t)*((1.d0+At2)/ppt))
    H=gamm*Logofterm((betabygamm)*(t*t)*((1.d0+At2)/ppt))

    fxc=n*(ex*Fx+ec+H)

    dFds = (2.d0*mu*s)/(1.d0+(mu/kappa)*(s**2))**2
    dFdsbg = ((2.d0*mu)/(1.d0+(mu/kappa)*(s**2))**2)/(2.d0*kf*n)
    dHdt = (2.d0*t*beta*gamm*(1.d0+2.d0*At2))/&
&        ((gamm*ppt+beta*t*t*(1.d0+At2))*ppt)
    dHdtbg = ((2.d0*beta*gamm*(1.d0+ &
&        2.d0*At2))/((gamm*ppt+beta*t*t*(1.d0+At2))*ppt))/(2.d0*ks*n)
    !decdrs=-(2.d0*AA*a1)*ddlog(1.d0+1.d0/(2*AA*pprs)) &
    decdrs=-(2.d0*AA*a1)*term &
&        +((1.d0+a1*rs)*((b1+3*b3*rs)/(2.d0*SQRT(rs))+ &
&          b2+2*b4*rs))/(pprs*(pprs+1.d0/(2.d0*AA)))
    dHdA=((2.d0+At2)*(At2*t*t*t*t*beta*gamm))/&
&        ((gamm*ppt+beta*t*t*(1.d0+At2))*ppt)
    dAdrs=-ddexp(-ec/gamm)*A*A*decdrs/beta
    dHdrs=dHdA*dAdrs

    dfxcdn = (4.d0/3.d0)*ex*(Fx-dFds*s)+ec-(decdrs*rs)/3.d0+H-(dHdrs*rs)/3.d0 &
&        - (7.d0/6.d0)*dHdt*t

    dfxcdgbg = ex*dFdsbg/(2.d0*kf) + dHdtbg/(2.d0*ks)

    IF ((ABS(fxc).GT.1.d65).OR.(ABS(dfxcdn).GT.1.d65).OR.&
&             (ABS(dfxcdgbg).GT.1.d65)) THEN
       WRITE(6,*) 'Problem in PBE',n,g,rs,s,t,ec,A,H
    ENDIF


    RETURN
  END SUBROUTINE pbefunc

  !*******************************************************************
  !
  ! Function Logofterm -- needed to take care of behavior for small term
  !  Evaluates log(1+term)

  FUNCTION Logofterm(term)

    REAL(8) :: term, Logofterm

    IF (ABS(term)>machine_precision) THEN
       Logofterm=ddlog(1.d0+term)
    ELSE
       Logofterm=term
    ENDIF

    RETURN
  END FUNCTION Logofterm

  !*******************************************************************
  !
  ! Function Aofec -- needed to take care of behavior for small ec

  FUNCTION Aofec(ec)

    REAL(8) :: ec, Aofec

    IF (ABS(ec)>machine_precision) THEN
       Aofec=betabygamm/(ddexp(-ec/gamm)-1.d0)
    ELSEIF (ABS(ec)>machine_zero) THEN
       Aofec=beta/(-ec)
    ELSE
       Aofec=-beta*DSIGN(machine_infinity,ec)
    ENDIF

    RETURN
  END FUNCTION Aofec

  !********************************************************************
  !
  ! Subroutine radialexcpbe
  !   Density(:) input on a uniform radial mesh of Npts
  !   Grid%r(:) input mesh points
  !   Exc - output integrated exchange correlation energy   -- in Rydberg units
  !   vxc(:) -- output exchange correlation potential       -- in Rydberg units
  !
  !********************************************************************

  SUBROUTINE radialexcpbe(Grid,density,Exc,vxc,fin)

    IMPLICIT NONE

    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: density(:)
    REAL(8), INTENT(OUT) :: Exc, vxc(:)
    INTEGER, INTENT(IN), OPTIONAL :: fin

    INTEGER :: i,j , ierr,Npts
    REAL(8),ALLOCATABLE :: gradient(:),gradmag(:),gxc(:),dgxcdr(:),fxc(:)
    REAL(8) :: dfxcdn,dfxcdgbg,sgn
    REAL(8) :: h

    !rewind(201); rewind(202)
    Npts=Grid%n
    IF (PRESENT(fin)) Npts=fin
    ALLOCATE(gradient(Npts),gradmag(Npts),gxc(Npts),dgxcdr(Npts), &
&        fxc(Npts),stat=i)
    IF (i /=0) THEN
       WRITE(6,*) 'error in radialexcpbe allocation ', Npts,i
       STOP
    ENDIF

    ! if (.not.scalarrelativistic) then
    CALL derivative(Grid,density(1:Npts),gradient(1:Npts),1,Npts)
    ! else
    !    CALL simplederiv(Grid,density(1:Npts),gradient(1:Npts),1,Npts)
    ! endif

    gradmag=ABS(gradient)

    DO i=1,Npts
       CALL  pbefunc(density(i),gradmag(i),fxc(i),dfxcdn,dfxcdgbg)
       vxc(i)=dfxcdn
       gxc(i)=dfxcdgbg*gradient(i)

       !write(201,'(i5,1p,6e15.7)') i,Grid%r(i),density(i),gradient(i),vxc(i),gxc(i)
    ENDDO

    ! if (.not.scalarrelativistic) then
    CALL derivative(Grid,gxc(1:Npts),dgxcdr(1:Npts),1,Npts)
    ! else
    !    CALL simplederiv(Grid,gxc(1:Npts),dgxcdr(1:Npts),1,Npts)
    ! endif

    DO i=2,Npts
       fxc(i)=2*fxc(i)*4*pi*(Grid%r(i)**2)  !2* changes from Har to Ryd
       vxc(i)=2*vxc(i)-2*dgxcdr(i)-4*gxc(i)/Grid%r(i)  ! Correction thanks
       ! to Marc Torrent and Francois Jollet
    ENDDO
    fxc(1)=0
    CALL extrapolate(Grid,vxc)

    !Do i=1,Npts
    !  write(202,'(i5,1p,6e15.7)') i,Grid%r(i),density(i),fxc(i),vxc(i),dgxcdr(i)
    !enddo

    Exc = integrator(Grid,fxc,1,Npts)

    RETURN
  END SUBROUTINE radialexcpbe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mggascanfunc contains only x part and is not working
!x  subroutine mggascanfunc(Grid,tmpd,tmpt,eexc,vxc,vtau)
!x     IMPLICIT NONE
!x     TYPE(GridInfo), INTENT(IN) :: Grid
!x     REAL(8), INTENT(IN) :: tmpd(:),tmpt(:)
!x     REAL(8), INTENT(OUT) :: eexc,vxc(:),vtau(:)
!x !   For eexc=int(Gxc)
!x !           vxc=partial Gxc/partial n
!x !           gxc=partial Gxc/partial abs(grad(n))/abs(grad(n))
!x !           vtau=partial Gxc/partial tau
!x 
!x   ! parameters set in PRL 115, 036402 (2015)
!x     REAL(8), PARAMETER :: sch0x=1.174d0  
!x     REAL(8), PARAMETER :: sca1=4.9479d0
!x     REAL(8), PARAMETER :: scc1x=0.667d0
!x     REAL(8), PARAMETER :: scc2x=0.8d0
!x     REAL(8), PARAMETER :: scdx=1.24d0
!x     REAL(8), PARAMETER :: sck1=0.065d0
!x     REAL(8), PARAMETER :: scmuAK=0.1234567901d0
!x     REAL(8), PARAMETER :: scb1=0.1566320775d0
!x     REAL(8), PARAMETER :: scb2=0.1208304597d0
!x     REAL(8), PARAMETER :: scb3=0.5d0
!x     REAL(8), PARAMETER :: scb4=0.1218315102d0
!x 
!x     REAL(8), PARAMETER :: tol=1.d-64
!x     INTEGER :: i,j,k,l,m,n
!x     REAL(8), allocatable :: grad(:),exc(:),deriv(:),arg(:),gxc(:),stuff(:)
!x     REAL(8) :: den,term1,g,a,s,s2,x,h1x,gx,fx,ffx,abrev,dfxda,dh1xda
!x     REAL(8) :: dffxda,dffxds,abrev1,dh1xds,term2
!x     REAL(8) :: dsdn,dsdg,dadn,dadg,dadt,dterm1dn
!x 
!x     n=Grid%n
!x     vxc=0.d0;vtau=0.d0;eexc=0.d0
!x     allocate( grad(n),exc(n),deriv(n),arg(n),gxc(n),stuff(n))
!x     grad=0.d0;exc=0.d0;deriv=0.d0;arg=0.d0;gxc=0.d0
!x     call derivative(Grid,tmpd,grad,1,n)
!x !!!Testing  -- explicitly H atom  
!x     grad(1:n)=-2.d0*tmpd(1:n)
!x !!!Testing    
!x     do i=1,n
!x        den=max(tmpd(i),tol)
!x        g=abs(grad(i))   
!x        term1=-0.7385587660d0*(den**(4.d0/3.d0))
!x        dterm1dn=-0.7385587660d0*(4.d0/3.d0)*(den**(1.d0/3.d0))
!x        term2=0.1616204597d0*(den**(-4.d0/3.d0))
!x        s=term2*g
!x        dsdn=-(4.d0/3.d0)*0.1616204597d0*g*(den**(-7.d0/3.d0))
!x        dsdg=term2
!x        s2=s*s
!x        a=0.3482823063d0*(tmpt(i)-(g**2)/(8*den))*(den**(-5.d0/3.d0))
!x        dadt=0.3482823063d0*(den**(-5.d0/3.d0))
!x        dadn=-(5.d0/3.d0)*0.3482823063d0*(tmpt(i))*(den**(-8.d0/3.d0)) &
!x &      + (8.d0/3.d0)*0.3482823063d0*(g**2)*(den**(-11.d0/3.d0))/8 
!x        dadg=0.3482823063d0*(g)*(den**(-8.d0/3.d0))/4
!x        x=scmuAK*s2*(1+(scb4*s2/scmuAK)*exp(-scb4*s2/scmuAK)) &
!x &             * (scb1*s+scb2*(1.d0-a)*exp(-scb3*(1.d0-a)**2))**2
!x        abrev1=exp(-scb4*s2/scmuAK)
!x        abrev=exp(-scb3*(1.d0-a)**2)
!x        h1x=1.d0+sck1-sck1/(1.d0*x/sck1)
!x        gx=1.d0
!x        if(s>0.d0)gx=gx-exp(-sca1/sqrt(s))
!x        fx=0.d0
!x        if(a<1.d0)  fx=exp(-scc1x*a/(1.d0-a))
!x        if(a>1.d0)  fx=-scdx*exp(-scc2x/(a-1.d0))
!x        ffx=(h1x+fx*(sch0x-h1x))*gx
!x        exc(i)=term1*ffx
!x        dfxda=0.d0
!x        if(a<1.d0)  dfxda=-scc1x*fx/(1.d0-a)**2
!x        if(a>1.d0)  dfxda=scc2x*fx/(1.d0-a)**2
!x        dh1xda= 2.d0*abrev*scb2*(scb1*s2+scb2*(1.d0-a)*abrev)*&
!x &               (-1.d0+2*scb3*(1.d0-a)**2)/(1.d0+x/sck1)**2
!x        dffxda=gx*(dfxda*(sch0x-h1x)+dh1xda*(1.d0-fx))
!x        vtau(i)=term1*dffxda*dadt
!x        dh1xds=(2*scmuAK*s*abrev1*(((1.d0+scb4*s2/scmuAK)*(1.d0-2*scb4*s2)) &
!x &       +2*s2*scb4)+4*scb1*s*(scb1*s2+scb2*(1.d0-a)*abrev))/(1.d0+x/sck1)**2               
!x        dffxds=0.d0
!x        if(s>0.d0)gxc(i)=dffxds&
!x &          -0.5d0**sca1*(exp(-sca1/sqrt(s))/(s**1.5d0))*(h1x+fx*(sch0x-h1x))
!x        dffxds=dffxds+dh1xds*gx*(1.d0-fx)
!x        arg(i)=term1*dffxds*dsdg   !  g term alone
!x        gxc(i)=arg(i)+term1*dffxda*dadg
!x        vxc(i)=dterm1dn*ffx &
!x &       +term1*dffxds*dsdn+term1*dffxda*dadn
!x        write(101,'(1p,50e15.7)') &
!x &       Grid%r(i),den,s,a,x,exc(i),vtau(i),arg(i),gxc(i),vxc(i)
!x 
!x     enddo    
!x     open(1001,file='pexc',form='formatted')
!x       arg=gxc*grad
!x       call derivative(Grid,arg,deriv,1,n)
!x       stuff=vxc
!x       vxc(2:n)=vxc(2:n)-deriv(2:n)-2*arg(2:n)/Grid%r(2:n)
!x       call extrapolate(Grid,vxc)
!x       do i=1,n
!x         write(1001,'(1p,50e15.7)')Grid%r(i),tmpd(i),grad(i),tmpt(i),exc(i),vtau(i),stuff(i),gxc(i),arg(i),deriv(i),vxc(i)
!x       enddo  
!x    close(1001)
!x     open(1001,file='veff',form='formatted')
!x       arg=log(1.d0+vtau)
!x       call derivative(Grid,arg,deriv,1,n)
!x       do i=1,n
!x        write(1001,'(1p,50e15.7)')Grid%r(i),tmpd(i),exc(i),deriv(i),vxc(i)/(1.d0+vtau(i))
!x        enddo
!x     close(1001)   
!x    stop   
!x ! note at the moment vxc is not rvxc   and still in Hartree units 
!x !   only x included -- c not included yet
!x     deallocate(grad,exc,deriv,arg,gxc)
!x 
!x end   subroutine mggascanfunc
!x 

  !*******************************************************************
  SUBROUTINE exch(Grid,den,rvxc,etxc,eexc,fin,v0,v0p,needvtau,tau,vtau)
    !  calculate exchange correlation potentials and energies
    !    for density functional theory from electron density
    !  den(n) is electron density * (4*pi*r**2)
    !  rvxc(n) is returned as vxc * r
    !  eexc is the total exchange energy (int(den*exc))
    !  etxc is eexc - int(den*vxc)
    !  fin (optional) is integer range of densities and potentials
    !  v0  (optional) is extrapolated value of vxc for r=0
    !  v0p  (optional) is extrapolated value of dvxc/dr for r=0
    !  needvtau (optional) is logical .true. if mgga
    !  tau(n) (optional) is kinetic energy density * (4*pi*r**2)
    !  vtau(n) (optional) is kinetic energy contribution for mgga
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: den(:)
    REAL(8), INTENT(INOUT) :: rvxc(:),etxc,eexc
    INTEGER, INTENT(IN), OPTIONAL :: fin
    REAL(8), INTENT(OUT), OPTIONAL :: v0,v0p
    LOGICAL, INTENT(OUT), OPTIONAL :: needvtau
    REAL(8), INTENT(IN),OPTIONAL :: tau(:)
    REAL(8), INTENT(INOUT),OPTIONAL :: vtau(:)
    REAL(8), ALLOCATABLE ::  tmpd(:),tmpv(:),dum(:),lrho(:)
            !  lrho(n)  is Laplacian of the density for mgga
    REAL(8), ALLOCATABLE :: exci(:),dfxcdgbg(:),gxc(:),dgxcdr(:)
    REAL(8), ALLOCATABLE :: grad(:),gradmag(:)
    REAL(8), ALLOCATABLE :: tmpt(:),tmpl(:),dgxcdl(:)
    REAL(8) :: fpi
    INTEGER :: i,n,i1,i2
    REAL(8) :: r,r2,rho,exc,vxc

    !write(6,*) 'In exch ', itype   ; call flush_unit(6)
    n=Grid%n
    IF (PRESENT(fin)) n=fin
    fpi=4*pi
    !write(6,*) 'In exch ', n   ; call flush_unit(6)
    rvxc=0;etxc=0;eexc=0
    if (PRESENT(v0)) v0=0
    if (PRESENT(v0p)) v0p=0
    if (PRESENT(needvtau)) needvtau=.false.
    if (PRESENT(vtau)) vtau=0.d0
    If (itype==GGA_PBE.or.itype==GGA_PBESOL) then !!!!!!!PBE form!!!!!!
       ALLOCATE(tmpd(n),tmpv(n))

       tmpd=0
       DO i=2,n
          tmpd(i)=den(i)/(fpi*(Grid%r(i)**2))
       ENDDO
       CALL extrapolate(Grid,tmpd)

       IF (PRESENT(fin)) THEN
          CALL radialexcpbe(Grid,tmpd,eexc,tmpv,fin)
       ELSE
          CALL radialexcpbe(Grid,tmpd,eexc,tmpv)
       ENDIF

       !WRITE(6,*) 'eexc',eexc

       IF (PRESENT(v0).AND.PRESENT(v0p)) THEN
          CALL derivative(Grid,tmpv,tmpd,1,15)
          v0=tmpv(1)
          v0p=tmpd(1)
       ENDIF


       DO i=1,n
          rvxc(i)=tmpv(i)*Grid%r(i)
          tmpv(i)=tmpv(i)*den(i)
       ENDDO

       etxc=eexc-integrator(Grid,tmpv(1:n),1,n)
       DEALLOCATE(tmpd,tmpv)

    ELSE IF (itype==LDA_PW) then !!! ! Perdew-Wang LDA !!!!
       !     write(6,*) 'LDA -- '; call flush_unit(6)
       ALLOCATE(tmpd(n),tmpv(n),dum(n))
       tmpd=0;tmpv=0;rvxc=0;dum=0
       DO i=2,n
          r=Grid%r(i)
          r2=r*r
          tmpd(i)=den(i)/(fpi*r2)
       ENDDO
       CALL extrapolate(Grid,tmpd)
       DO i=1,n
          CALL pwldafunc(tmpd(i),exc,vxc)
          tmpd(i)=den(i)*(exc-vxc)
          tmpv(i)=den(i)*exc
          rvxc(i)=Grid%r(i)*vxc
          IF (PRESENT(v0).AND.PRESENT(v0p)) THEN
             IF (i==1) v0=vxc
             dum(i)=vxc
          ENDIF
       ENDDO
       !  calculate exchange-correlation contribution to the potential
       !       etxc=dsum(n,a,1)*h
       !       eexc=dsum(n,b,1)*h
       etxc=integrator(Grid,tmpd(1:n),1,n)
       eexc=integrator(Grid,tmpv(1:n),1,n)
       IF (PRESENT(v0).AND.PRESENT(v0p)) THEN
          CALL derivative(Grid,dum,tmpd,1,15)
          v0p=tmpd(1)
       ENDIF
       DEALLOCATE(tmpd,tmpv,dum)
    ELSE IF (itype==LIBXC.and.have_libxc) then !!!!! External LibXC library !!!!!
       allocate(tmpd(n),tmpv(n),exci(n),tmpt(n))
       tmpd=0.d0; tmpv=0.d0; exci=0.d0;tmpt=0.d0
       tmpd(2:n)=den(2:n)/(fpi*(Grid%r(2:n)**2))
       call extrapolate(Grid,tmpd)
       if (libxc_islda()) then
        call libxc_getvxc(n,exci,tmpv,1,tmpd)
       elseif (libxc_isgga()) then
        allocate(grad(n),gradmag(n),gxc(n),dgxcdr(n),dfxcdgbg(n))
        grad=0.d0;gradmag=0.d0;dgxcdr=0.d0;dfxcdgbg=0.d0
        call derivative(Grid,tmpd,grad,1,n)
        gradmag=ABS(grad)
        call libxc_getvxc(n,exci,tmpv,1,tmpd,grho=gradmag,vxcgr=dfxcdgbg)
        gxc(1:n)=dfxcdgbg(1:n)*grad(1:n)
        call derivative(Grid,gxc,dgxcdr,1,n)
        tmpv(2:n)=tmpv(2:n)-dgxcdr(2:n)-2.d0*gxc(2:n)/Grid%r(2:n)
        call extrapolate(Grid,tmpv)
        deallocate(grad,gradmag,gxc,dgxcdr,dfxcdgbg)
        !  not sure why/if this is needed
        do i=1,n
           if (isnan(tmpv(i))) tmpv(i)=0.d0
           if (isnan(exci(i))) exci(i)=0.d0
        enddo   
       elseif (libxc_ismgga()) then
        write(6,*) '  atompaw not yet available for mgga -- stop '
        call flush_unit(6)
        stop       
        if(PRESENT(needvtau)) needvtau=.true.       
        allocate(grad(n),gradmag(n),gxc(n),dgxcdr(n),dfxcdgbg(n))
        grad=0.d0;gradmag=0.d0;gxc=0.d0;dgxcdr=0.d0;dfxcdgbg=0.d0
        call derivative(Grid,tmpd,grad,1,n)
        gradmag=ABS(grad)
        tmpt(2:n)=tau(2:n)/(fpi*(Grid%r(2:n)**2))
        call extrapolate(Grid,tmpt)
        if(libxc_needs_laplacian()) then
          allocate(dum(n),lrho(n))
          dum=0.d0;lrho=0
          call derivative(Grid,grad,dum,1,n)
          lrho(2:n)=dum(2:n)+2.d0*grad(2:n)/Grid%r(2:n)
          call extrapolate(Grid,lrho); dum=0.d0       
          call libxc_getvxc(n,exci,tmpv,1,tmpd,grho=gradmag,lrho=lrho, &
&            tau=tmpt,vxcgr=dfxcdgbg,vxclrho=dum,vxctau=vtau)
          gxc(1:n)=dfxcdgbg(1:n)*grad(1:n)
          call derivative(Grid,gxc,dgxcdr,1,n)
          tmpv(2:n)=tmpv(2:n)-dgxcdr(2:n)-2.d0*gxc(2:n)/Grid%r(2:n)
          call derivative(Grid,dum,gxc,1,n)  ! add Laplacian contributions
          call derivative(Grid,gxc,dgxcdr,1,n)
          tmpv(2:n)=tmpv(2:n)+dgxcdr(2:n)+2.d0*gxc(2:n)/Grid%r(2:n)
          call extrapolate(Grid,tmpv)
          deallocate(grad,gradmag,gxc,dgxcdr,dfxcdgbg,dum,lrho)
        else   ! no Laplacian contribution
          call libxc_getvxc(n,exci,tmpv,1,tmpd,grho=gradmag, &
&            tau=tmpt,vxcgr=dfxcdgbg,vxctau=vtau)
!          open(1001,file='expot',form='formatted')
!          do i=1,n
!             write(1001,'(1p,50e15.7)')Grid%r(i),tmpd(i),gradmag(i),tmpt(i),exci(i),tmpv(i),dfxcdgbg(i),vtau(i)
!             enddo
!          close(1001)
!          stop   
          gxc(1:n)=dfxcdgbg(1:n)*grad(1:n)
          call derivative(Grid,gxc,dgxcdr,1,n)
          tmpv(2:n)=tmpv(2:n)-dgxcdr(2:n)-2.d0*gxc(2:n)/Grid%r(2:n)
          call extrapolate(Grid,tmpv)
          deallocate(grad,gradmag,gxc,dgxcdr,dfxcdgbg)
        endif
       else
           write(6,*) 'unknown libxc family -- need to work harder '
           stop        
       end if
       rvxc=0.d0
       rvxc(1:n)=tmpv(1:n)*Grid%r(1:n)
       exci(1:n)=exci(1:n)*tmpd(1:n)*fpi*Grid%r(1:n)**2
       eexc=integrator(Grid,exci,1,n)
       if (present(v0).and.present(v0p)) then
        call derivative(Grid,tmpv,tmpd,1,15)
        v0=tmpv(1);v0p=tmpd(1)
       endif
       tmpv(1:n)=tmpv(1:n)*den(1:n)
       etxc=eexc-integrator(Grid,tmpv(1:n),1,n)
       WRITE(6,*) 'etxc,eexc = ',etxc,eexc;call flush_unit(6)
       open(1001,file='expot',form='formatted')
       do i=1,n
          write(1001,'(1p,20e15.7)')Grid%r(i),rvxc(i),den(i),exci(i)
          enddo   
       deallocate(tmpd,tmpv,exci,tmpt)
    else

       WRITE(6,*) 'Warning (EXCOR): ', itype,' no results returned !'
       STOP
    END if

    !write(6,*) 'Exiting exch '; call flush_unit(6)
  END SUBROUTINE exch

END MODULE excor
