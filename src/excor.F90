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

  USE io_tools
  USE Tools
  USE globalmath
  USE gridmod
  USE libxc_mod
  USE r2scanmod
  USE atomdata, only : scalarrelativistic,diracrelativistic,needvtau

  IMPLICIT NONE

  CHARACTER(132),PUBLIC,SAVE :: exctype
  INTEGER, PRIVATE, SAVE :: itype
  INTEGER, PRIVATE, PARAMETER :: LDA_PW=14
  INTEGER, PRIVATE, PARAMETER :: GGA_PBE=16
  INTEGER, PRIVATE, PARAMETER :: GGA_PBESOL=18
  INTEGER, PRIVATE, PARAMETER :: MGGA_R2SCAN_001=13001
  INTEGER, PRIVATE, PARAMETER :: MGGA_R2SCAN_01=1301
  INTEGER, PRIVATE, PARAMETER :: LIBXC=-1
  INTEGER, PRIVATE, PARAMETER :: NO_XC=-100

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
    WRITE(STD_OUT,*) exctype; call flush_unit(std_out)

    SELECT CASE(TRIM(exctype))
    CASE default
      IF (have_libxc) THEN
        call libxc_getid_fromInput(exctype,id)
        call libxc_init_func(id,1)
        itype = LIBXC
        WRITE(STD_OUT,*) 'Using Libxc -- ', TRIM(exctype) ; call flush_unit(std_out)
        if(needvtau.and.(.not.libxc_needs_tau())) then
          WRITE(STD_OUT,*) 'Problem with XC functional choice -- need KE based mgga form for vtau '   
          WRITE(STD_OUT,*) '    Program stopping '
          stop
        endif    
      ELSE
        WRITE(STD_OUT,*) 'ERROR: cannot understand Exchange-Correlation functional from input file!'
        WRITE(STD_OUT,*) '       possible issue: atompaw not compiled with libXC.'
        stop
        !itype = LDA_PW
        !exctype='LDA-PW'
        !WRITE(STD_OUT,*) 'Perdew-Wang correlation'
      END IF
    CASE('LDA-PW')
       itype = LDA_PW
       WRITE(STD_OUT,*) 'Perdew-Wang correlation'
    CASE('GGA-PBE')
       itype = GGA_PBE
       WRITE(STD_OUT,*) 'Perdew-Burke-Ernzerhof GGA'
       mu=muorig
       beta=betorig
       betabygamm=beta/gamm
    CASE('GGA-PBESOL')
       itype = GGA_PBESOL
       WRITE(STD_OUT,*) 'Perdew-Burke-Ernzerhof modified (PBEsol) GGA'
       mu=musol
       beta=betsol
       betabygamm=beta/gamm
    CASE ('MGGA-R2SCAN-001')
       itype = MGGA_R2SCAN_001    
       WRITE(STD_OUT,*) 'R2SCAN MGGA with eta=0.001'
       call r2scaninit(0.001d0)
       needvtau=.true.
    CASE ('MGGA-R2SCAN-01')
       itype = MGGA_R2SCAN_01    
       WRITE(STD_OUT,*) 'R2SCAN MGGA with eta=0.01'
       call r2scaninit(0.01d0)
       needvtau=.true.
    CASE ('HF')
       itype = NO_XC
       WRITE(STD_OUT,*) 'No XC'
    END SELECT
    WRITE(STD_OUT,*)

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
     CASE(MGGA_R2SCAN_001)
        write(if,*) 'Furness et al. J. Phys. Chem. Lett. 2020, 11, 8208−8215 with eta=0.001'
        write(if,*)  'https://dx.doi.org/10.1021/acs.jpclett.0c02405'
        write(if,*)  'Using modified KS equations for tau dependence'
     CASE(MGGA_R2SCAN_01)
        write(if,*) 'Furness et al. J. Phys. Chem. Lett. 2020, 11, 8208−8215 with eta=0.01'
        write(if,*)  'https://dx.doi.org/10.1021/acs.jpclett.0c02405'
        write(if,*)  'Using modified KS equations for tau dependence'
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
       WRITE(STD_OUT,*) 'Problem in PW',n,rs,ec
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
       WRITE(STD_OUT,*) 'Problem in PBE',n,g,rs,s,t,ec,A,H
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
       WRITE(STD_OUT,*) 'error in radialexcpbe allocation ', Npts,i
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

  !*******************************************************************
  SUBROUTINE exch(Grid,den,rvxc,etxc,eexc,fin,v0,v0p,tau,vtau)
    !  calculate exchange correlation potentials and energies
    !    for density functional theory from electron density
    !  den(n) is electron density * (4*pi*r**2)
    !  rvxc(n) is returned as vxc * r
    !  eexc is the total exchange energy (int(den*exc))
    !  etxc is eexc - int(den*vxc)
    !  fin (optional) is integer range of densities and potentials
    !  v0  (optional) is extrapolated value of vxc for r=0
    !  v0p  (optional) is extrapolated value of dvxc/dr for r=0
    !  tau(n) (optional) is kinetic energy density * (4*pi*r**2)
    !    Note that tau(n) is always included but only used from some
    !      mgga formulations
    !  vtau(n) (optional) is kinetic energy contribution for mgga
    !  vtau is stored only if needvtau==.true.
    !  needvtau is a global variable which is true for mgga and
    !      explicit treatment of kinetic energy density tau

    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(IN) :: den(:)
    REAL(8), INTENT(INOUT) :: rvxc(:),etxc,eexc
    INTEGER, INTENT(IN), OPTIONAL :: fin
    REAL(8), INTENT(OUT), OPTIONAL :: v0,v0p
    REAL(8), INTENT(IN),OPTIONAL :: tau(:)
    REAL(8), INTENT(INOUT),OPTIONAL :: vtau(:)
    REAL(8), ALLOCATABLE ::  tmpd(:),tmpv(:),dum(:),lrho(:)
            !  lrho(n)  is Laplacian of the density for mgga
    REAL(8), ALLOCATABLE :: exci(:),dfxcdgbg(:),gxc(:),dgxcdr(:)
    REAL(8), ALLOCATABLE :: grad(:),gradmag(:),sigma(:)
    REAL(8), ALLOCATABLE :: tmpt(:),tmpl(:),dgxcdl(:)
    REAL(8), ALLOCATABLE :: tmpvt(:),dum1(:)
    REAL(8), ALLOCATABLE :: dexcdn(:),dexcds(:)
    REAL(8) :: fpi,A
    INTEGER :: i,n,i1,i2
    REAL(8) :: r,r2,rho,exc,vxc

    !write(std_out,*) 'In exch ', itype   ; call flush_unit(std_out)
    n=Grid%n
    IF (PRESENT(fin)) n=fin
    fpi=4*pi
    !write(std_out,*) 'In exch ', n   ; call flush_unit(std_out)
    rvxc=0;etxc=0;eexc=0
    if (PRESENT(v0)) v0=0
    if (PRESENT(v0p)) v0p=0
    if (needvtau.and.(PRESENT(vtau).eqv..false.)) then
       write(6,*) 'exch:  Inconsistency in mgga ', needvtau,PRESENT(vtau)
       write(6,*) 'exch:  stopping execution '
       call flush_unit(6)
       stop
    endif

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

       !WRITE(STD_OUT,*) 'eexc',eexc

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
       !     write(std_out,*) 'LDA -- '; call flush_unit(std_out)
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
       
    ELSE IF (itype==MGGA_R2SCAN_001.or.itype==MGGA_R2SCAN_01) then !!r2scan 
       allocate(tmpd(n),tmpv(n),exci(n),tmpt(n),dum(n),dum1(n))
       tmpd=0.d0; tmpv=0.d0; exci=0.d0;tmpt=0.d0
       dum1=0.d0
       tmpd(2:n)=den(2:n)/(fpi*(Grid%r(2:n)**2))
       call extrapolate(Grid,tmpd)
       allocate(grad(n),sigma(n),dexcdn(n),dexcds(n))
       grad=0.d0;sigma=0.d0;dexcdn=0.d0;dexcds=0.d0
       !call derivative(Grid,tmpd,grad,1,n)
       do i=1,n
          dum1(i)=ddlog(tmpd(i))
       enddo   
       call derivative(Grid,dum1,grad,1,n)
       grad(1:n)=grad(1:n)*tmpd(1:n)     ! perhaps more accurate???
       sigma=grad**2
       !   Prepare kinetic energy input tau -- used for most mgga
       tmpt(2:n)=tau(2:n)/(fpi*(Grid%r(2:n)**2))
       call extrapolate(Grid,tmpt)
       ! convert to Hartree units
       tmpt=0.5d0*tmpt

       do i=1,n
         call r2scanfun(tmpd(i),grad(i),tmpt(i),&
&                 exci(i),vtau(i),dexcdn(i),dexcds(i))
       enddo
       ! convert to Rydberg units
       exci=2*exci
       !!!wrong -- dimensionless vtau=2*vtau
       dexcdn=2*dexcdn
       dexcds=4*dexcds       !extra factor of two due to sigma=grad**2
       rvxc=0.d0
       rvxc=dexcdn
       dum(1:n)=dexcds(1:n)
       call derivative(Grid,dum,dum1,1,n)
       dum(2:n)=2*dum(2:n)/(Grid%r(2:n))
       call extrapolate(Grid,dum)
       rvxc=rvxc-dum1-dum
       !IF (PRESENT(v0).AND.PRESENT(v0p)) THEN
       !  v0=rvxc(1)
       !  call derivative(Grid,rvxc,dum,1,15)
       !  v0p=dum(1)
       !ENDIF  
       dum(1:n)=exci(1:n)*fpi*(Grid%r(1:n)**2)
       eexc=integrator(Grid,dum(1:n),1,n)
       dum(1:n)=rvxc(1:n)*den(1:n)
       etxc=eexc-integrator(Grid,dum(1:n),1,n)
      ! rvxc(1:n)=rvxc(1:n)*Grid%r(1:n)
        rvxc(1:n)=(dexcdn(1:n)-dum1(1:n))*Grid%r(1:n)-2*dexcds(1:n)
       deallocate(grad,sigma,dexcdn,dexcds)
       deallocate(tmpd,tmpv,exci,tmpt,dum,dum1)
    ELSE IF (itype==LIBXC.and.have_libxc) then !!!!! External LibXC library !!!!!
       !write(std_out,*)    'libxc  case in exch -- islda isgga ismgga', libxc_islda(),libxc_isgga(), libxc_ismgga()
!! 6/9/2023 NAWH changed the order of the if statements to accomodate the possibility of multiple
!!       types for exchange and correlation which can be accomodated properly by the order: ismgga, isgga, islda
       allocate(tmpd(n),tmpv(n),exci(n),tmpt(n))
       tmpd=0.d0; tmpv=0.d0; exci=0.d0;tmpt=0.d0
       tmpd(2:n)=den(2:n)/(fpi*(Grid%r(2:n)**2))
       call extrapolate(Grid,tmpd)

       if (libxc_ismgga()) then
        !write(std_out,*) 'found libxc_ismgga '; call flush_unit(std_out)
        allocate(grad(n),gradmag(n),gxc(n),dgxcdr(n),dfxcdgbg(n))
        grad=0.d0;gradmag=0.d0;gxc=0.d0;dgxcdr=0.d0;dfxcdgbg=0.d0
        tmpv=0.d0
        do i=1,n
           tmpv(i)=ddlog(tmpd(i))
        enddo   
        call derivative(Grid,tmpv,grad,1,n)
        grad(1:n)=grad(1:n)*tmpd(1:n)     !  perhaps more accurate???
        tmpv=0.d0
        gradmag=ABS(grad)
        allocate(tmpl(n),dgxcdl(n),tmpvt(n),dum(n),dum1(n))
        tmpl=0.d0; dgxcdl=0.d0;tmpvt=0.d0;dum=0.d0;dum1=0.d0
        !   Prepare kinetic energy input tau -- used for most mgga
        if(libxc_needs_tau()) then
           write(STD_OUT,*) 'Needs tau '; call flush_unit(STD_OUT)
           tmpt(2:n)=tau(2:n)/(fpi*(Grid%r(2:n)**2))
           call extrapolate(Grid,tmpt)
        end if
        if(libxc_needs_laplacian()) then
           write(STD_OUT,*) 'Needs laplacian '; call flush_unit(STD_OUT)
          call derivative(Grid,grad,tmpl,1,n)
          tmpl(2:n)=tmpl(2:n)+2.d0*grad(2:n)/Grid%r(2:n)
          call extrapolate(Grid,tmpl)
        endif  
!!!     Note: libxc_getvxc in libxc_mod.F90 takes care of Ry<->Ha conversion 
        call libxc_getvxc(n,exci,tmpv,1,tmpd,grho=gradmag,lrho=tmpl, &
&            tau=tmpt,vxcgr=dfxcdgbg,vxclrho=dgxcdl,vxctau=tmpvt)
        gxc(1:n)=dfxcdgbg(1:n)*grad(1:n)
        call derivative(Grid,gxc,dgxcdr,1,n)
        tmpv(2:n)=tmpv(2:n)-dgxcdr(2:n)-2.d0*gxc(2:n)/Grid%r(2:n)
        if(libxc_needs_laplacian()) then
           write(STD_OUT,*) 'Needs laplacian '; call flush_unit(STD_OUT)
          dum=0.d0      
          call derivative(Grid,dgxcdl,dum,1,n)
          dum1=0.d0      
          call derivative(Grid,dum,dum1,1,n)
          tmpv(2:n)=tmpv(2:n)+dum1(2:n)+2.d0*dum(2:n)/Grid%r(2:n)
        endif  
        call extrapolate(Grid,tmpv)
        if (needvtau) vtau=tmpvt
        deallocate(grad,gradmag,gxc,dgxcdr,dfxcdgbg)
        deallocate(tmpl,dgxcdl,tmpvt,dum,dum1)
        do i=1,n
           if (isnan(tmpv(i)).or.abs(tmpv(i)).lt.machine_zero) tmpv(i)=0.d0
           if (isnan(exci(i)).or.abs(exci(i)).lt.machine_zero) exci(i)=0.d0
           if (isnan(vtau(i)).or.abs(vtau(i)).lt.machine_zero) vtau(i)=0.d0
           if (isnan(tmpd(i)).or.abs(tmpd(i)).lt.machine_zero) tmpd(i)=0.d0
        enddo

       elseif (libxc_isgga()) then
        allocate(grad(n),gradmag(n),gxc(n),dgxcdr(n),dfxcdgbg(n))
        grad=0.d0;gradmag=0.d0;dgxcdr=0.d0;dfxcdgbg=0.d0
        do i=1,n
           tmpv(i)=ddlog(tmpd(i))
        enddo   
        call derivative(Grid,tmpv,grad,1,n)
        grad(1:n)=grad(1:n)*tmpd(1:n)     !  perhaps more accurate???
        gradmag=ABS(grad)
        call libxc_getvxc(n,exci,tmpv,1,tmpd,grho=gradmag,vxcgr=dfxcdgbg)
        gxc(1:n)=dfxcdgbg(1:n)*grad(1:n)
        call derivative(Grid,gxc,dgxcdr,1,n)
        tmpv(2:n)=tmpv(2:n)-dgxcdr(2:n)-2.d0*gxc(2:n)/Grid%r(2:n)
        call extrapolate(Grid,tmpv)
        deallocate(grad,gradmag,gxc,dgxcdr,dfxcdgbg)
        do i=1,n
           if (isnan(tmpv(i)).or.abs(tmpv(i)).lt.machine_zero) tmpv(i)=0.d0
           if (isnan(exci(i)).or.abs(exci(i)).lt.machine_zero) exci(i)=0.d0
           if (isnan(tmpd(i)).or.abs(tmpd(i)).lt.machine_zero) tmpd(i)=0.d0
        enddo

       elseif (libxc_islda()) then
        call libxc_getvxc(n,exci,tmpv,1,tmpd)
        do i=1,n
           if (isnan(tmpv(i)).or.abs(tmpv(i)).lt.machine_zero) tmpv(i)=0.d0
           if (isnan(exci(i)).or.abs(exci(i)).lt.machine_zero) exci(i)=0.d0
           if (isnan(tmpd(i)).or.abs(tmpd(i)).lt.machine_zero) tmpd(i)=0.d0
        enddo

       else
           write(std_out,*) 'unknown libxc family -- need to work harder '
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
       WRITE(STD_OUT,*) 'etxc,eexc = ',etxc,eexc;call flush_unit(std_out)
       deallocate(tmpd,tmpv,exci,tmpt)
    else

       WRITE(STD_OUT,*) 'Warning (EXCOR): ', itype,' no results returned !'
       STOP
    END if

    !write(std_out,*) 'Exiting exch '; call flush_unit(std_out)
  END SUBROUTINE exch

END MODULE excor
