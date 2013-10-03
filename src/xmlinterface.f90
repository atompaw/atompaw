!! NAME
!!  XMLInterface
!!
!! FUNCTION
!!  This module contains routines used to convert PAW atomic dataset into a XML file
!!
!! COPYRIGHT
!! Copyright (C) 2013-2013 ATOMPAW group (NHolzwarth, MTorrent, FJollet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt.

Module XMLInterface

 use atomdata
 use GlobalMath
 use gridmod
 use excor
 use pseudo
 use pkginfo

 implicit none

 private

 public :: Atompaw2XML

!!=================================================================
!! CONSTANTS
!!=================================================================

!Nb of bytes related to subtypes
 integer, parameter :: i1b=selected_int_kind(2)
 integer, parameter :: i2b=selected_int_kind(4)
 integer, parameter :: i4b=selected_int_kind(9)
 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0d0,1.0d0))
 integer, parameter :: lgt=kind(.true.)

!Default lengths
 integer, parameter :: fnlen=132    ! maximum length of file name variables
 integer, parameter :: strlen=32000 ! maximum length of input string

!UNIX unit numbers : standard input, standard output
 integer, parameter :: std_in=5,std_out=6

!Real, fractionary and other constants
 real(dp), parameter :: zero =0._dp
 real(dp), parameter :: one  =1._dp
 real(dp), parameter :: two  =2._dp
 real(dp), parameter :: three=3._dp
 real(dp), parameter :: four =4._dp
 real(dp), parameter :: half =0.50_dp
 real(dp), parameter :: tol8= 0.00000001_dp
 real(dp), parameter :: tol10=0.0000000001_dp
 real(dp), parameter :: tol12=0.000000000001_dp
 real(dp), parameter :: tol14=0.00000000000001_dp
 character(len=1), parameter :: ch10 = char(10)

!Real physical constants
!Revised fundamental constants from http://physics.nist.gov/cuu/Constants/index.html
!(from 2002 least squares adjustment)
 real(dp), parameter :: Bohr_Ang=0.5291772108_dp    ! 1 Bohr, in Angstrom
 real(dp), parameter :: Ha_eV=27.2113845_dp ! 1 Hartree, in eV

!!=================================================================
!! DEFAULT PARAMETERS
!!=================================================================

!Version number
 character*10 :: atompaw2xmlver='2.0.1', verdate='june 2013'

!Default unit for XML file
 integer, parameter :: unit_xml=22

!Options
 logical :: usexcnhat_def=.false.

!Real Space Optimization (default parameters)
 logical  :: userso_def=.false.
 real(dp) :: ecut_rso_def=10._dp
 real(dp) :: gfact_rso_def=two
 real(dp) :: werror_rso_def=0.0001_dp

 real(dp),parameter :: rmax_vloc=10._dp !We decide to cut at r=10 u.a because of numeric noise...

!!=================================================================
!! STRUCTURED DATATYPES
!!=================================================================

!Definition of various radial meshes
 type mesh_data_type
  integer :: mesh_type           ! Default type of meshes (lin. or log.)
  integer :: nmesh               ! Number of meshes
  integer :: iwavmesh            ! Index of mesh for partial waves
  integer :: iprjmesh            ! Index of mesh for projectors
  integer :: icoremesh           ! Index of mesh for core density
  integer :: ivionmesh           ! Index of mesh for vion potential
  integer :: ivbaremesh          ! Index of mesh for vbare potential
  integer :: ivalemesh           ! Index of mesh for valence density
  integer :: wav_meshsz          ! Size of mesh for partial waves
  integer :: sph_meshsz          ! Size of mesh for partial waves
  integer :: prj_meshsz          ! Size of mesh for projectors
  integer :: core_meshsz         ! Size of mesh for core density
  integer :: vion_meshsz         ! Size of mesh for vion potential
  integer :: vbare_meshsz        ! Size of mesh for vbare potential
  integer :: vale_meshsz         ! Size of mesh for valence density
  integer :: prj_msz_max         ! Maximum size for projector (used for RSO)
  real(dp) :: rad_step           ! Default value for radial step
  real(dp) :: log_step           ! Default value for log step
  integer,pointer :: meshtp(:)   ! Array storing mesh type for all meshes
  integer,pointer :: meshsz(:)   ! Array storing mesh size for all meshes
  real(dp),pointer :: radstp(:)  ! Array storing radial step for all meshes
  real(dp),pointer :: logstp(:)  ! Array storing log step for all meshes
 end type mesh_data_type

!Options for REAL SPACE OPTIMIZATION of projectors
 type pawrso_type
  logical :: userso           ! TRUE if Real Space Optimization is required
  real(dp) :: ecut            ! Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
  real(dp) :: gfact           ! Real Space Optimization parameter: Gamma/Gmax ratio
  real(dp) :: werror          ! Real Space Optimization parameter: max. error W_l allowed
 end type pawrso_type


 CONTAINS

!!=================================================================
!! NAME
!! atompaw2xml
!!
!! FUNCTION
!! Main routine for building a XML PAW dataset file.
!!
!!       Units: Energies=Hartree, Lengths=Bohr
!!
!! PARENTS
!! atompaw
!!
!! CHILDREN
!! rdinputxml,build_mesh_data,opt_proj_rso,xmloutput
!!=================================================================

 subroutine Atompaw2XML(AEOrbit,AEPot,AESCF,PAW,FC,Grid,ifinput)

 integer,intent(in) :: ifinput
 type(OrbitInfo), intent(in)     :: AEOrbit
 type(PotentialInfo), intent(in) :: AEPot
 type (SCFInfo),  intent(in)     :: AESCF
 type(Pseudoinfo), intent(in)    :: PAW
 type(FCInfo), intent(in)        :: FC
 type(GridInfo), intent(in)      :: Grid

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: id,vlocopt
 type(mesh_data_type) :: mesh_data
 type(pawrso_type) :: pawrso
 character*(fnlen) :: author,file_xml,xcname
 character*(5000) :: input_string
 real(dp),allocatable :: tproj(:,:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 write(std_out,'(2a)') ch10,'========================================================'
 write(std_out,'(3a)')      '==   Atompaw2XML - v',trim(atompaw2xmlver),&
&                                                       ':                         =='
 write(std_out,'(a)')       '==   Writes a PAW dataset generated by "Atompaw"      =='
 write(std_out,'(a)')       '==   into a XML file...                                  =='
 write(std_out,'(a)')       '========================================================'

 call read_inputstring(ifinput,input_string)

!------------------------------------------------------------------
!---- Read choices from input file

 call rdinputxml(vlocopt,pawrso,author,input_string)

!------------------------------------------------------------------
!---- Build radial meshes definitions
 call build_mesh_data(mesh_data,Grid,PAW%irc,PAW%ivion,PAW%ivale,PAW%coretailpoints)

!------------------------------------------------------------------
!--Check that the pseudo valence density is zero at grid end
 if (mesh_data%vale_meshsz>0) then
  write(std_out,'(/,2x,a,a,f5.2,a,a,g11.4)') 'Atompaw2XML info:',&
&   '  At r=',Grid%r(mesh_data%vale_meshsz),' a.u.,',&
&   ' Pseudo_valence_density= ', PAW%tden(mesh_data%vale_meshsz) &
&                               /Grid%r(mesh_data%vale_meshsz)**2/(4*pi)
  write(std_out,'(2x,a)') '  This quantity must be as small as possible.'
 endif

!------------------------------------------------------------------
!--Check that the potential is zero at grid end
 if (vlocopt==1) then
  write(std_out,'(/,2x,a,a,f5.2,a,a,g11.4)') 'Atompaw2XML info:',&
&   '  At r_vloc=',Grid%r(mesh_data%vion_meshsz),' a.u.,',&
&   ' VHartree(ntild(Zv+Zc))= -Zv/r + ', half*PAW%abinitvloc(mesh_data%vion_meshsz) &
&                                    +(AEPot%nz-FC%zcore)/Grid%r(mesh_data%vion_meshsz)
  write(std_out,'(2x,a)') '  This quantity must be as small as possible.'
 endif
 if (vlocopt==2) then
  write(std_out,'(/,2x,a,a,f5.2,a,a,g11.4)') 'Atompaw2XML info:',&
&   '  At r_vloc=',Grid%r(mesh_data%vion_meshsz),' a.u.,',&
&   ' VHartree(ntild(Zv+Zc))= -Zv/r + ', half*PAW%abinitnohat(mesh_data%vion_meshsz) &
&                                    +(AEPot%nz-FC%zcore)/Grid%r(mesh_data%vion_meshsz)
  write(std_out,'(2x,a)') '  This quantity must be as small as possible.'
 endif

!------------------------------------------------------------------
!---- Eventually optimize projectors
 allocate(tproj(mesh_data%prj_msz_max,PAW%nbase))
 tproj=zero
 do id=1,PAW%nbase
  tproj(1:mesh_data%wav_meshsz,id)=PAW%otp(1:mesh_data%wav_meshsz,id)
 end do
 call opt_proj_rso(tproj,mesh_data,pawrso,Grid,PAW)

!------------------------------------------------------------------
!---- Write PAW dataset file in XML format
 xcname=exctype;if (have_libxc) call libxc_getshortname(exctype,xcname)
 file_xml=TRIM(AEpot%sym)//'.'//TRIM(xcname)//'-paw'

 call xmloutput(file_xml,Grid,AESCF,AEPot,FC,PAW,mesh_data,tproj,vlocopt,&
&               input_string,author)

!------------------------------------------------------------------
!---- End
 deallocate(tproj)
 call destroy_mesh_data(mesh_data)

 write(std_out,'(2x,a,a)') 'Atompaw2XML ended.',ch10

 end subroutine Atompaw2XML


!!=================================================================
!! NAME
!! rdinputxml
!!
!! FUNCTION
!! Read the input file in order to get XML dataset options
!!
!! INPUTS
!!
!! OUTPUT
!!  vlocopt= option for Vloc (use of nhat in XC or not)
!!  pawrso
!!    %ecut=Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
!!    %gfact=Real Space Optimization parameter: Gamma/Gmax ratio
!!    %userso=TRUE if REAL Space Optimization is required
!!    %werror=Real Space Optimization parameter: max. error W_l allowed
!!
!! SIDE EFFECTS
!!  input_string= string containing a copy of atompaw input file
!!                appended here
!!
!! PARENTS
!! atompaw2xml
!!
!! CHILDREN
!! extractword,uppercase
!!=================================================================

 subroutine rdinputxml(vlocopt,pawrso,author,input_string)

 integer,intent(out) :: vlocopt
 type(pawrso_type),intent(out) :: pawrso
 character(len=*),intent(out) :: author
 character(len=*),intent(inout) :: input_string

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: i_author,i_usexcnhat,i_rsoptim,iend,ok,nn
 character*(fnlen) :: readline,readline_u,inputline,inputword

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 read(std_in,'(a)',advance='no',iostat=ok) readline

 readline_u=readline
 call Uppercase(readline_u)

 i_usexcnhat=index(readline_u,'USEXCNHAT')
 i_rsoptim  =index(readline_u,'RSOPTIM')
 i_author   =index(readline_u,'AUTHOR')

!Option for use of NHAT in XC
 if (i_usexcnhat>0) then
  vlocopt=1
 else
  if (     usexcnhat_def) vlocopt=1
  if (.not.usexcnhat_def) vlocopt=2
 end if
 if (vlocopt==1) &
&  write(std_out,'(a,2x,a)') ch10,'Zero potential and Blochl local ionic potential will be output in XML file'
 if (vlocopt/=1) &
&  write(std_out,'(a,2x,a)') ch10,'Zero potential and Kresse local ionic potential will be output in XML file'

!Option for use of REAL SPACE OPTIMIZATION
 pawrso%userso=userso_def
 pawrso%ecut=ecut_rso_def
 pawrso%gfact=gfact_rso_def
 pawrso%werror=werror_rso_def
 if (i_rsoptim>0) then
  pawrso%userso=.true.
  iend=128
  if (i_usexcnhat>i_rsoptim.and.i_usexcnhat-1<iend) iend=i_usexcnhat-1
  inputline=""
  if (iend>i_rsoptim+7) inputline=trim(readline(i_rsoptim+7:iend))
  if (inputline/="") then
   call extractword(1,inputline,inputword);inputword=trim(inputword)
   if (inputword/="") then
    read(inputword,*) pawrso%ecut
    call extractword(2,inputline,inputword);inputword=trim(inputword)
    if (inputword/="") then
     read(inputword,*) pawrso%gfact
     call extractword(3,inputline,inputword);inputword=trim(inputword)
     if (inputword/="") read(inputword,*) pawrso%werror
    end if
   end if
  end if
 end if
 if (pawrso%userso) then
  write(std_out,'(a,2x,a,f4.1,a,f3.1,a,g6.1,a)') ch10,&
    'Real Space optim.: Ecut, Gamma/Gmax, Wl(error) [',&
    pawrso%ecut,', ',pawrso%gfact,', ',pawrso%werror,']'
 else
  write(std_out,'(a,2x,a)') ch10,'No Real Space Optimization of projectors'
 end if

 if(i_author>0) then
   inputline=trim(readline(i_author+6:))
   read(unit=inputline,fmt=*) author
   author=trim(author) ; nn=len(trim(author))
   write(unit=input_string,fmt='(6a)') trim(input_string),char(10),&
&   "XMLOUT",char(10),readline(1:i_author-1),trim(readline(i_author+nn+10:))
 else
   author=""
   write(unit=input_string,fmt='(5a)') trim(input_string),char(10),&
&   "XMLOUT",char(10),trim(readline)
 end if

 end subroutine rdinputxml


!!=================================================================
!! NAME
!! opt_proj_rso
!!
!! FUNCTION
!! Apply Real Space Optimization (RSO) on non-local projectors in order
!! to smooth them and cut large reciprocal vectors contribution.
!! Directly written from:
!!  RD King-smith, MC Payne and JS Lin, Phys. Rev. B, 44 (1991), 13063
!!
!! INPUTS
!!  mesh_data=datastructure containing mesh size definitions
!!  pawrso
!!    %ecut_rso=Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
!!    %gfact_rso=Real Space Optimization parameter: Gamma/Gmax ratio
!!    %userso=TRUE if REAL Space Optimization is required
!!    %werror_rso=Real Space Optimization parameter: max. error W_l allowed
!!  Grid= grid info from atompaw
!!  PAW= PAW datastructure from atompaw
!!
!! SIDE EFFECTS
!!  mesh_data
!!    %prj_meshsz= Dimension of radial mesh for tproj
!!  tproj(prj_msz_max,nbase)= projectors on partial waves
!!
!! PARENTS
!! atompaw2xml
!!
!! CHILDREN
!! aamat,simpson_int,gauleg,jbessel,DGETRF,DGETRS
!!=================================================================

 subroutine opt_proj_rso(tproj,mesh_data,pawrso,Grid,PAW)

 type(mesh_data_type),intent(inout) :: mesh_data
 type(pawrso_type),intent(in)    :: pawrso
 type(Gridinfo),intent(in) :: Grid
 type(Pseudoinfo),intent(in) ::  PAW
 real(dp),intent(inout) :: tproj(mesh_data%prj_msz_max,PAW%nbase)

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer,parameter :: ngaussmax=50
 integer, parameter :: nqmax=500

 integer :: ib,info,iq,iqp,ir,ll,nqgauss1,nqgauss2,r0_meshsz,rp_meshsz
 real(dp) :: amat,bess,bessp,dum,dq,gamma,gmax,lstep,r0,rc_prj,rstep,wwmax,wwl,xx

 integer, allocatable :: iwork(:)
 real(dp),allocatable :: am(:,:),bm(:),chi1g(:,:),chi2g(:),chireg(:,:),ff(:),&
&                        gg(:),qgauss1(:),qgauss2(:),wgauss1(:),wgauss2(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (.not.pawrso%userso) return

 write(std_out,'(3(/,2x,a),/,2x,a,f7.2,/,2x,a,f6.2,/,2x,a,g11.3)') &
&     'Atompaw2XML info:',&
&     '  Optimizing non-local projectors',&
&     '  using Real Space Optimization...',&
&     '  Parameters: Ecut (Hartree)=',pawrso%ecut,&
&     '              Gamma/Gmax    =',pawrso%gfact,&
&     '              Wl max (error)=',pawrso%werror

!Initialize data for optimization
 rstep=mesh_data%radstp(mesh_data%iprjmesh)
 lstep=mesh_data%logstp(mesh_data%iprjmesh)
 if (usingloggrid(Grid)) then
  rc_prj=rstep*(exp(dble(mesh_data%prj_meshsz-1)*lstep)-1.d0)
 else
  rc_prj=rstep*dble(mesh_data%prj_meshsz-1)
 end if
 r0=rc_prj/1.035_dp
 gmax=sqrt(two*pawrso%ecut)
 gamma=pawrso%gfact*gmax
 wwmax=10000000000._dp

 rp_meshsz=mesh_data%prj_meshsz;allocate(ff(rp_meshsz))

!Define q mesh for reciprocal space
 nqgauss1=ngaussmax
 allocate(qgauss1(nqgauss1),wgauss1(nqgauss1))
 call gauleg(zero,gmax,qgauss1,wgauss1,nqgauss1)
 nqgauss2=ngaussmax
 allocate(qgauss2(nqgauss2),wgauss2(nqgauss2))
 call gauleg(gmax,gamma,qgauss2,wgauss2,nqgauss2)
 allocate(chi1g(nqgauss1,PAW%nbase),chireg(nqmax,PAW%nbase))
 allocate(am(nqgauss2,nqgauss2),bm(nqgauss2),chi2g(nqgauss2))

!Transfer tproj(r) into chi(q) for 0<q<gmax
!-- On a Gaussian mesh
 do ib=1,PAW%nbase
  ll=PAW%l(ib)
  do iq=1,nqgauss1
   do ir=1,rp_meshsz
    call jbessel(bess,bessp,dum,ll,0,qgauss1(iq)*Grid%r(ir))
    ff(ir)=bess*Grid%r(ir)*tproj(ir,ib)
   enddo
   call simpson_int(ff,Grid,rp_meshsz,chi1g(iq,ib))
  enddo
!-- On a regular mesh
  dq=gmax/dble(nqmax-1)
  do iq=1,nqmax
   do ir=1,rp_meshsz
    call jbessel(bess,bessp,dum,ll,0,dble(iq-1)*dq*Grid%r(ir))
    ff(ir)=bess*Grid%r(ir)*tproj(ir,ib)
   enddo
   call simpson_int(ff,Grid,rp_meshsz,chireg(iq,ib))
  enddo
 enddo

!Loop on error Wl
 do while (wwmax>pawrso%werror)
  wwmax=-one
  if (usingloggrid(Grid)) then
   r0_meshsz=max(int(log(one+(r0*1.035_dp)/rstep)/lstep)+1,mesh_data%prj_meshsz+1)
  else
   r0_meshsz=max(int(r0*1.035_dp/rstep)+1,mesh_data%prj_meshsz+1)
   if (mod(r0_meshsz,2)==0) r0_meshsz=r0_meshsz-1
  end if
  if (r0_meshsz>mesh_data%prj_msz_max) then
   write(std_out,'(/,2x,a)') 'Error in Atompaw2XML(opt_proj_rsp): ro_meshsz too big !'
   stop
  endif
  r0=Grid%r(r0_meshsz)
  allocate(gg(r0_meshsz))

!Loop on (l,n) basis
 do ib=1,PAW%nbase
  ll=PAW%l(ib)

!Compute chi(q) for gmax<q<gamma on Gauss mesh
! --Loop on q
  do iq=1,nqgauss2

!  --Compute bm(q)
   bm(iq)=zero
   do iqp=1,nqgauss1
    call aamat(amat,qgauss2(iq),qgauss1(iqp),r0,ll)
    bm(iq)=bm(iq)+amat*chi1g(iqp,ib)*wgauss1(iqp)
   enddo

!  --Compute am(q,qp)
   do iqp=1,iq
    call aamat(amat,qgauss2(iq),qgauss2(iqp),r0,ll)
    am(iq,iqp)=-amat*wgauss2(iqp)
    am(iqp,iq)=-amat*wgauss2(iq)
   enddo
   am(iq,iq)=am(iq,iq)+(half*pi)*qgauss2(iq)**2

! --End loop on q
  enddo


! Solve Am(q,qp).X(q)=Bm(q)
  allocate(iwork(nqgauss2))
  call DGETRF(nqgauss2,nqgauss2,am,nqgauss2,iwork,info)
  call DGETRS('N',nqgauss2,1,am,nqgauss2,iwork,bm,nqgauss2,info)
  deallocate(iwork)

  chi2g=bm

! Transfer back chi(q) into tproj(r)
  do ir=1,r0_meshsz
   xx=zero
   do iq=1,nqgauss1
    call jbessel(bess,bessp,dum,ll,0,qgauss1(iq)*Grid%r(ir))
    xx=xx+wgauss1(iq)*bess*chi1g(iq,ib)*qgauss1(iq)**2
   enddo
   do iq=1,nqgauss2
    call jbessel(bess,bessp,dum,ll,0,qgauss2(iq)*Grid%r(ir))
    xx=xx+wgauss2(iq)*bess*chi2g(iq)*qgauss2(iq)**2
   enddo
   tproj(ir,ib)=two*Grid%r(ir)*xx/pi
  enddo

!Estimate the error W_l(q)
! --Compute Int(0,R0) [r**2.chi(r).jl(qr)] (and Wl)
! --for each q of the regular mesh
  wwl=-one
  do iq=1,nqmax
   do ir=1,r0_meshsz
    call jbessel(bess,bessp,dum,ll,0,dble(iq-1)*dq*Grid%r(ir))
    gg(ir)=bess*Grid%r(ir)*tproj(ir,ib)
   enddo
   call simpson_int(gg,Grid,r0_meshsz,xx)
   wwl=max(abs(chireg(iq,ib)-xx),wwl)
  enddo
  wwl=wwl/maxval(abs(chireg(:,ib)))
  wwmax=max(wwl,wwmax)

!End loop on ib
 enddo

!End loop on error
 deallocate(gg)
 enddo

 deallocate(am,bm,chi1g,chi2g,chireg,ff,wgauss1,qgauss1,wgauss2,qgauss2)

 mesh_data%prj_meshsz=r0_meshsz
 write(std_out,'(2(a,f7.4),a)') '  New radius R0 for nl projectors (Bohr)=',&
&                         r0,' (=',r0/rc_prj,'*Rc(proj))'
 if (r0>1.55_dp*rc_prj) &
&  write(std_out,'(4(/,a))') &
&                'Warning:',&
&                '  Radius for nl projectors (R0) seems to be high !',&
&                '  You should change parameters of Real Space Optimization',&
&                '  (increase Ecut, Gamma/Gmax or Wl).'

 end subroutine opt_proj_rso

!------------------------------------------------------------------
!---- aamat: useful matrix for Real Space Optimization of projectors
!------------------------------------------------------------------

 subroutine aamat(amat,qq,qqp,r0,ll)
 integer,intent(in) :: ll
 real(dp),intent(in) :: qq,qqp,r0
 real(dp),intent(out) :: amat

 real(dp) :: bess,bessp,bess_,bessp_,dum

 call jbessel(bess,bessp,dum,ll,1,qqp*r0)
 if (dabs(qq-qqp)<tol10) then
  amat=qq**3*half*r0**2*(bess*bessp &
&      +qq*r0*bessp**2+qq*r0*bess**2 &
&      -dble(ll*(ll+1))/(qq*r0)*bess**2)
 else
  call jbessel(bess_,bessp_,dum,ll,1,qq*r0)
  amat=qq**2*qqp**2/(qqp**2-qq**2)*r0**2 &
&     *(bess*bessp_*qq-bess_*bessp*qqp)
 endif

 end subroutine aamat


!!=================================================================
!! NAME
!! xmloutput
!!
!! FUNCTION
!! write the paw data file in XML format
!!
!! INPUTS
!! fname=file name (without .xml suffixe)
!! Grid= Grid datastructure from atompaw
!! AESCF= AESCF datastructure from atompaw
!! AEPOT= AEPOT datastructure from atompaw
!! PAW= PAW datastructure from atompaw
!! FC= FC datastructure from atompaw
!! mesh_data= datatructure containing the definition of
!!            all the radial meshes used in the XML file
!! tproj(prj_msz_max,nbase)= PAW projectors
!!       (might be modified by real space optimization)
!! vlocopt= option for local potential (1=Blochl, 2=Kresse)
!! input_string= string containing a copy of atompaw input file
!!
!! PARENTS
!! atompaw2xml
!!
!! CHILDREN
!!
!!=================================================================

 SUBROUTINE xmloutput(fname,Grid,AESCF,AEPot,FC,PAW,mesh_data,tproj,&
&                     vlocopt,input_string,author)

 integer,intent(in) :: vlocopt
 character(len=*),intent(in) :: input_string,author,fname
 TYPE(Gridinfo),intent(in) :: Grid
 TYPE (SCFInfo),intent(in) :: AESCF
 TYPE(Potentialinfo),intent(in) :: AEPot
 TYPE (FCInfo),intent(in) :: FC
 TYPE(Pseudoinfo),intent(in) :: PAW
 TYPE(mesh_data_type),intent(in) :: mesh_data
 real(dp) :: tproj(mesh_data%prj_msz_max,PAW%nbase)

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ib,ic,ii,n,irc,OK
 character(len=3) :: gridtype
 character(len=4) :: char4
 character(len=5) :: char5a,char5b,xc_type
 character(len=20) :: char20
 character(len=132) :: xc_name
 real(dp) :: radstp0,logstp0,sqr4pi
 character(len=3) :: gridt(mesh_data%nmesh)
 real(dp),allocatable :: dum(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!Some defs
 gridtype="log"
 sqr4pi=sqrt(4*pi)
 n=Grid%n;irc=PAW%irc

!Open file for writing
 OPEN(unit_xml,file=TRIM(fname)//'.xml',form='formatted')

!Write header
 WRITE(unit_xml,'("<?xml  version=""1.0""?>")')
 WRITE(unit_xml,'("<paw_setup version=""0.5"">")')
 write(unit=char5a,fmt='(f5.2)') AEPot%zz
 WRITE(unit_xml,'("<atom symbol=""",a,""" Z=""",a,$)') &
&   trim(ADJUSTL(AEPot%sym)),trim(ADJUSTL(char5a))
 write(unit=char5a,fmt='(f5.2)') FC%zcore
 write(unit=char5b,fmt='(f5.2)') AEPot%nz-FC%zcore
 WRITE(unit_xml,'(""" core=""",a,""" valence=""",a,"""/>")') &
&   trim(ADJUSTL(char5a)),trim(ADJUSTL(char5b))

!Write XC definition
 if (trim(ADJUSTL(exctype))=="LDA-PW") then
   xc_type="LDA"
   xc_name="PW"
 elseif (trim(ADJUSTL(exctype))=="GGA-PBE") then
   xc_type="GGA"
   xc_name="PBE"
 elseif (trim(ADJUSTL(exctype))=="GGA-PBESOL") then
   xc_type="GGA"
   xc_name="PBESOL"
 else if (have_libxc) then
   if (libxc_isgga()) then
     xc_type="GGA"
   else
     xc_type="LDA"
   end if
   xc_name=trim(exctype)
!else ! old coding
!  xc_type="LIBXC"
!  write(xc_name,'(i6)') -pshead%pspxc_abinit
 else
   stop "Error in Atompaw2XML: unknown XC type !"
 end if
 WRITE(unit_xml,'("<xc_functional type=""",a,""" name=""",a,"""/>")') &
&      TRIM(xc_type),TRIM(xc_name)

!Generator data
 if (scalarrelativistic) then
   WRITE(unit_xml,'("<generator type=""scalar-relativistic"" name=""atompaw"">")')
 else
   WRITE(unit_xml,'("<generator type=""non-relativistic"" name=""atompaw"">")')
 endif
 WRITE(unit_xml,'("</generator>)")')

!Echo input file
!WRITE(unit_xml,'("<!-- Contact info: email: natalie@wfu.edu")')
!WRITE(unit_xml,'("               web: pwpaw.wfu.edu")')
!WRITE(unit_xml,'(" Energy units=Hartree, length units=bohr")')
!WRITE(unit_xml,'(" Note: consistent with 06-14-06 standard")')
!WRITE(unit_xml,'(" As discussed at CECAM PAW Workshop     ")')
!Call PrintDate(unit_xml, ' PAW functions generated on ')
!WRITE(unit_xml,'("  ",a)') Trim(PAW%Vloc_description)
!WRITE(unit_xml,'("  ",a)') Trim(PAW%Proj_description)
!WRITE(unit_xml,'("  ",a)') Trim(PAW%Comp_description)
 WRITE(unit_xml,'("<!-- Atompaw ",a)') atp_version
 WRITE(unit_xml,'(" Contact info: email: natalie@wfu.edu")')
 WRITE(unit_xml,'("               web:   pwpaw.wfu.edu")')
 WRITE(unit_xml,'(" Energy units=Hartree, length units=bohr")')
 Call PrintDate(unit_xml, ' PAW functions generated on ')
 if (trim(author)/="") WRITE(unit_xml,'(a,a)') ' by ',trim(author)
 WRITE(unit_xml,'(" Program:  atompaw - input data follows: ")')
 WRITE(unit_xml,'(a)') trim(input_string)
 write(unit_xml,'(a)') "END"
 WRITE(unit_xml,'(" Program:  atompaw - input end -->")')

!Echo input file
 WRITE(unit_xml,'("  Program:  atompaw - input data follows: ")')
 WRITE(unit_xml,'(a)') trim(input_string)
 write(unit_xml,'(a)') "0"
 WRITE(unit_xml,'("  Program:  atompaw - input end -->")')

!Energies
 WRITE(unit_xml,'("<ae_energy kinetic=""",1pe25.17,""" xc=""",1pe25.17,"""")') &
&      AESCF%ekin/2,AESCF%eexc/2
 WRITE(unit_xml,'("  electrostatic=""",1pe25.17,""" total=""",1pe25.17,"""/>")') &
&      AESCF%estatic/2,AESCF%etot/2
 WRITE(unit_xml,'("<core_energy kinetic=""",1pe25.17,"""/>")') AESCF%corekin*0.5d0

!PAW radius
 WRITE(unit_xml,'("<PAW_radius rpaw=""",f13.10,"""/>")') PAW%rc

!Electronic configuration
 WRITE(unit_xml,'("<valence_states>")')
 do ib=1,PAW%nbase
   call mkname(ib,char4)
   char20=stripchar('"'//AEPot%sym//char4//'"')
   ii=min(ABS(PAW%np(ib)),20)
   WRITE(unit_xml,'("  <state n=""",i2,""" l=""",i1,""" f=""",1pe14.7,$)')&
&        ii,PAW%l(ib),PAW%occ(ib)
   WRITE(unit_xml,'(""" rc=""",f13.10,""" e=""",1pe14.7,""" id=",a6,"/>")')&
&        PAW%rcio(ib),PAW%eig(ib)*0.5d0,TRIM(char20)
 enddo
 WRITE(unit_xml,'("</valence_states>")')

!Radial meshes definitions
 do ii=1,mesh_data%nmesh
  select case(mesh_data%meshtp(ii))
   case(1)
    char20='r=d*i'
    gridt(ii)="lin"
    radstp0=zero
    logstp0=mesh_data%radstp(ii)
   case(2)
    char20='r=a*(exp(d*i)-1)'
    gridt(ii)="log"
    radstp0=mesh_data%radstp(ii)
    logstp0=mesh_data%logstp(ii)
   case default
    write(std_out, '(a,a,a)' )&
&     '  This mesh type ', ch10,&
&     '  is not implemented in Atompaw.'
    stop
  end select
  WRITE(unit_xml,'("<radial_grid eq=""",a,""" a=""",es22.16,$)') trim(char20),radstp0
  WRITE(unit_xml,'(""" d=""",es22.16,""" istart=""0"" iend=""",i5,$)') &
&  logstp0,mesh_data%meshsz(ii)-1
  WRITE(unit_xml,'(""" id=""",a,i1,"""/>")') gridt(ii),ii
 end do

!Compensation charge shape function
 if (gaussianshapefunction) then
   WRITE(unit_xml,'("<shape_function type=""gauss"" rc=""",f19.16"""/>")') PAW%gausslength
 else if (besselshapefunction) then
   WRITE(unit_xml,'("<shape_function type=""bessel"" rc=""",f19.16"""/>")') PAW%rc_shap
 else
   WRITE(unit_xml,'("<shape_function type=""sinc"" rc=""",f19.16"""/>")') PAW%rc_shap
 endif

!Core densities
 allocate(dum(mesh_data%meshsz(mesh_data%icoremesh)),stat=OK)
 dum=zero
 dum(2:mesh_data%meshsz(mesh_data%icoremesh))= &
&              sqr4pi*FC%coreden(2:mesh_data%meshsz(mesh_data%icoremesh))&
&              /(4*pi*(Grid%r(2:mesh_data%meshsz(mesh_data%icoremesh)))**2)
 call extrapolate(Grid,dum)
 WRITE(unit_xml,'("<ae_core_density grid=""",a,i1""">")') &
& gridt(mesh_data%icoremesh),mesh_data%icoremesh
 WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%icoremesh))
 WRITE(unit_xml,'("</ae_core_density>")')
 dum=zero
 dum(2:mesh_data%meshsz(mesh_data%icoremesh))= &
&              sqr4pi*PAW%tcore(2:mesh_data%meshsz(mesh_data%icoremesh))&
&               /(4*pi*(Grid%r(2:mesh_data%meshsz(mesh_data%icoremesh)))**2)
 call extrapolate(Grid,dum)
 WRITE(unit_xml,'("<pseudo_core_density grid=""",a,i1,""">")') &
& gridt(mesh_data%icoremesh),mesh_data%icoremesh
 WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%icoremesh))
 WRITE(unit_xml,'("</pseudo_core_density>")')
 DEALLOCATE(dum)

!Valence density
 ALLOCATE(dum(mesh_data%meshsz(mesh_data%ivalemesh)),stat=OK)
 dum=zero
 dum(2:mesh_data%meshsz(mesh_data%ivalemesh))= &
&              sqr4pi*PAW%tden(2:mesh_data%meshsz(mesh_data%ivalemesh))&
&               /(4*pi*(Grid%r(2:mesh_data%meshsz(mesh_data%ivalemesh)))**2)
 call extrapolate(Grid,dum)
 WRITE(unit_xml,'("<pseudo_valence_density grid=""",a,i1,""">")') &
& gridt(mesh_data%ivalemesh),mesh_data%ivalemesh
 WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%ivalemesh))
 WRITE(unit_xml,'("</pseudo_valence_density>")')
 DEALLOCATE(dum)

!Vbare potential
 ALLOCATE(dum(mesh_data%meshsz(mesh_data%ivbaremesh)),stat=OK)
 dum=zero
 dum(1:mesh_data%meshsz(mesh_data%ivbaremesh))= &
&              sqr4pi*half*PAW%vloc(1:mesh_data%meshsz(mesh_data%ivbaremesh))
 WRITE(unit_xml,'("<zero_potential grid=""",a,i1,""">")') &
& gridt(mesh_data%ivbaremesh),mesh_data%ivbaremesh
 WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%ivbaremesh))
 WRITE(unit_xml,'("</zero_potential>")')
 DEALLOCATE(dum)

!Local ionic potential
 if (vlocopt==1) then
  ALLOCATE(dum(mesh_data%meshsz(mesh_data%ivionmesh)),stat=OK)
  dum=zero
  dum(1:mesh_data%meshsz(mesh_data%ivionmesh))= &
&            sqr4pi*half*PAW%abinitvloc(1:mesh_data%meshsz(mesh_data%ivionmesh))
   WRITE(unit_xml,'("<kresse_joubert_local_ionic_potential grid=""",a,i1,""">")') &
&   gridt(mesh_data%ivionmesh),mesh_data%ivionmesh
   WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%ivionmesh))
   WRITE(unit_xml,'("</kresse_joubert_local_ionic_potential>")')
   DEALLOCATE(dum)
  end if

!Local Blochl''s potential
 if(vlocopt==2) then
  ALLOCATE(dum(mesh_data%meshsz(mesh_data%ivionmesh)),stat=OK)
  dum=zero
  dum(1:mesh_data%meshsz(mesh_data%ivionmesh))= &
&               sqr4pi*half*PAW%abinitnohat(1:mesh_data%meshsz(mesh_data%ivionmesh))
  WRITE(unit_xml,'("<blochl_local_ionic_potential grid=""",a,i1,""">")') &
&  gridt(mesh_data%ivionmesh),mesh_data%ivionmesh
  WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%ivionmesh))
  WRITE(unit_xml,'("</blochl_local_ionic_potential>")')
  DEALLOCATE(dum)
 endif

!Partial waves and projectors
 Do ib=1,PAW%nbase
  call mkname(ib,char4)
  char20=stripchar('"'//AEPot%sym//char4//'"')
  ALLOCATE(dum(mesh_data%meshsz(mesh_data%iwavmesh)),stat=OK)
  dum=zero
  dum(2:mesh_data%meshsz(mesh_data%iwavmesh))= &
&               PAW%ophi(2:mesh_data%meshsz(mesh_data%iwavmesh),ib) &
&              /Grid%r(2:mesh_data%meshsz(mesh_data%iwavmesh))
  call extrapolate(Grid,dum)
  WRITE(unit_xml,'("<ae_partial_wave state=",a6," grid=""",a,i1,""">")') &
&  TRIM(char20),gridt(mesh_data%iwavmesh),mesh_data%iwavmesh
  WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%iwavmesh))
  WRITE(unit_xml,'("</ae_partial_wave>")')
  dum=zero
  dum(2:mesh_data%meshsz(mesh_data%iwavmesh))= &
&               PAW%otphi(2:mesh_data%meshsz(mesh_data%iwavmesh),ib) &
&              /Grid%r(2:mesh_data%meshsz(mesh_data%iwavmesh))
  call extrapolate(Grid,dum)
  WRITE(unit_xml,'("<pseudo_partial_wave state=",a6," grid=""",a,i1,""">")')&
&  TRIM(char20),gridt(mesh_data%iwavmesh),mesh_data%iwavmesh
  WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%iwavmesh))
  WRITE(unit_xml,'("</pseudo_partial_wave>")')
  DEALLOCATE (dum)
  ALLOCATE(dum(mesh_data%meshsz(mesh_data%iprjmesh)),stat=OK)
  dum=zero
  dum(2:mesh_data%meshsz(mesh_data%iprjmesh))= &
&               tproj(2:mesh_data%meshsz(mesh_data%iprjmesh),ib) &
&              /Grid%r(2:mesh_data%meshsz(mesh_data%iprjmesh))
  call extrapolate(Grid,dum)
  WRITE(unit_xml,'("<projector_function state=",a6," grid=""",a,i1,""">")') &
&  TRIM(char20),gridt(mesh_data%iprjmesh),mesh_data%iprjmesh
  WRITE(unit_xml,'(3(1x,es23.16))') (dum(ii),ii=1,mesh_data%meshsz(mesh_data%iprjmesh))
  WRITE(unit_xml,'("</projector_function>")')
  DEALLOCATE (dum)
 Enddo

!Kinetic terms
 WRITE(unit_xml,'("<kinetic_energy_differences>")')
 WRITE(unit_xml,'(3(1x,es23.16))') ((PAW%kij(ib,ic)/2,ic=1,PAW%nbase),ib=1,PAW%nbase)
 WRITE(unit_xml,'("</kinetic_energy_differences>")')
 WRITE(unit_xml,'("</paw_setup>")')

!Cloe file
 CLOSE(unit_xml)

 END SUBROUTINE xmloutput


!!=================================================================
!! NAME
!! read_inputstring
!!
!! FUNCTION
!! Read the file echoing the atompaw input file
!! and transfer it into a character string
!!
!! INPUTS
!!  ifinput= unit number of the file to read
!!
!! OUTPUT
!!  input_string=character string containing the file
!!
!! PARENTS
!  xml2abinit
!!=================================================================

 subroutine read_inputstring(ifinput,input_string)

 integer,intent(in) :: ifinput
 character(len=*) :: input_string

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: OK
 character(len=132) :: inputline

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 open(ifinput,file='dummy',form='formatted')
 read(ifinput,'(a)',iostat=OK,end=10) inputline
 if (OK/=0) return
 input_string=trim(inputline)
 do
   read(ifinput,'(a)',iostat=OK,end=10) inputline
   if (OK/=0) exit
   write(unit=input_string,fmt='(3a)') trim(input_string),char(10),trim(inputline)
 enddo
 return
10 continue
 close(ifinput)

 end subroutine read_inputstring


!!=================================================================
!! NAME
!! simpson_int
!!
!! FUNCTION
!! Do integral using corrected Simpson rule (on a linear or logarithmic grid)
!!
!! INPUTS
!!  ff(meshsz)=integrand values
!!  meshsz=size of radial mesh for integration
!!  Grid
!!    %h= radial step used for integration
!!    %drdu(max_meshsz)= Factor used to compute radial integrals on generalized grid
!!
!! OUTPUT
!!  simp=resulting integral by corrected Simpson rule
!!
!! PARENTS
!! calc_dij0,calc_valden,calc_shapef
!!=================================================================

 subroutine simpson_int(ff,Grid,meshsz,simp)

 TYPE(Gridinfo),intent(in)    :: Grid
 integer,intent(in)           :: meshsz
 real(dp),intent(out)         :: simp
 real(dp),intent(in)          :: ff(meshsz)

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 real(dp),parameter :: eps=tol14**4
 integer :: ii,nn
 real(dp) :: simp1,simp2,simp4

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 ii=meshsz;do while(abs(ff(ii))<eps.and.ii>0);ii=ii-1;enddo
 if (ii<=0) then
   simp=zero;return
 end if
 nn=min(ii+1,meshsz)

 if (nn>=5) then

  if (mod(nn,2)==1) then
   simp1=ff(1)*Grid%drdu(1)+ff(nn)*Grid%drdu(nn)
   simp2=zero;simp4=ff(2)*Grid%drdu(2)
   do ii=3,nn-1,2
    simp2=simp2+ff(ii)  *Grid%drdu(ii)
    simp4=simp4+ff(ii+1)*Grid%drdu(ii+1)
   enddo
  else
   simp1=1.25_dp*ff(1)*Grid%drdu(1)+three*ff(2)*Grid%drdu(2) &
&       -0.25_dp*ff(3)*Grid%drdu(3)+ff(nn)*Grid%drdu(nn)
   simp2=zero;simp4=ff(3)*Grid%drdu(3)
   do ii=4,nn-1,2
    simp2=simp2+ff(ii)  *Grid%drdu(ii)
    simp4=simp4+ff(ii+1)*Grid%drdu(ii+1)
   enddo
  endif
  simp=Grid%h/three*(simp1+two*simp2+four*simp4)

 else if (nn==4) then
  simp=Grid%h*0.375_dp*(three*(ff(2)*Grid%drdu(2)+ff(3)*Grid%drdu(3)) &
&                   +     (ff(1)*Grid%drdu(1)+ff(4)*Grid%drdu(4)))
 else if (nn==3) then
  simp=Grid%h/three*(ff(1)*Grid%drdu(1)+four*ff(2)*Grid%drdu(2)+ff(3)*Grid%drdu(3))
 else if (nn==2) then
  simp=Grid%h*half*(ff(1)*Grid%drdu(1)+ff(2)*Grid%drdu(2))
 else
  simp=zero
 end if

 end subroutine simpson_int


!!=================================================================
!! NAME
!! gauleg
!!
!! FUNCTION
!! Compute supports and weights for Gauss-Legendre integration
!!
!! INPUTS
!!  x1=lower bound of integration
!!  x2=upper bound of integration
!!  n=order of integration
!!
!! OUTPUT
!!  x(n)=array of support points
!!  w(n)=array of integration weights
!!
!! PARENTS
!! opt_proj
!!
!! NOTES
!!   Code follows [W.H. Press et al., Numerical
!!   Recipes (Cambridge University Press, New York, 1986].
!!=================================================================

 subroutine gauleg(x1,x2,x,w,n)

 integer,intent(in)   :: n
 real(dp),intent(in)  :: x1,x2
 real(dp),intent(out) :: x(n),w(n)

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 real(dp), parameter :: eps=0.0000000000001_dp
 integer :: i,j,loop_root,m
 real(dp) :: p1,p2,p3,pp,xm,xl,z,z1

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 m=(n+1)/2
 xm=half*(x2+x1)
 xl=half*(x2-x1)

 do i=1,m
  z=cos(pi*(i-0.25_dp)/(n+half))

  do loop_root=1,100

   p1=one
   p2=zero

   do j=1,n
    p3=p2
    p2=p1
    p1=((two*j-one)*z*p2 - (j-one)*p3)/j
   enddo

   pp=n*(z*p1-p2)/(z*z-one)
   z1=z
   z=z1-p1/pp

   if (abs(z-z1) < eps) exit

  enddo
  x(i)=xm-xl*z
  x(n+1-i)=xm+xl*z
  w(i)=two*xl/((one-z*z)*pp*pp)
  w(n+1-i)=w(i)

 enddo

 end subroutine gauleg


!!=================================================================
!! NAME
!! build_mesh_data
!!
!! FUNCTION
!! Determine meshes definitions
!! (if necessary define a logarithmic radial grid)
!!
!! INPUTS
!!  Grid=grid datastructure in AtomPAW format
!!  irc=index of rc in Grid
!!
!! OUTPUT
!!  mesh_data
!!   Data defining various meshes
!!
!! PARENTS
!!
!!=================================================================

subroutine build_mesh_data(mesh_data,Grid,irc,ivion,ivale,coretailpoints)

 type(mesh_data_type),intent(out) :: mesh_data
 type(GridInfo),intent(in) :: Grid
 integer,intent(in) :: irc,ivion,ivale,coretailpoints

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer, parameter :: nmesh_max=10
 integer :: ii1

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 allocate(mesh_data%meshtp(nmesh_max))
 allocate(mesh_data%meshsz(nmesh_max))
 allocate(mesh_data%radstp(nmesh_max))
 allocate(mesh_data%logstp(nmesh_max))

!Mesh definition
 if (usingloggrid(Grid)) then
  mesh_data%mesh_type=2
  mesh_data%rad_step=Grid%drdu(1)
  mesh_data%log_step=Grid%h
 else
  mesh_data%mesh_type=1
  mesh_data%rad_step=Grid%h
  mesh_data%log_step=zero
 end if

!Various mesh sizes
 mesh_data%wav_meshsz=irc+Grid%ishift
 if (usingloggrid(Grid)) then
  mesh_data%sph_meshsz=min(1+nint(log(one+Grid%r(irc)/mesh_data%rad_step)/mesh_data%log_step),mesh_data%wav_meshsz)
 else
  mesh_data%sph_meshsz=min(1+nint(Grid%r(irc)/mesh_data%rad_step),mesh_data%wav_meshsz)
 end if
 mesh_data%prj_meshsz=mesh_data%sph_meshsz  ! To be modified by RSO

 mesh_data%core_meshsz=coretailpoints
 mesh_data%vale_meshsz=ivale

 mesh_data%vion_meshsz=ivion
 if (mesh_data%vion_meshsz<=0) then
  if (mesh_data%mesh_type==1) then
   ii1=int(one+rmax_vloc/mesh_data%rad_step)
  else
   ii1=int(log(one+rmax_vloc/mesh_data%rad_step)/mesh_data%log_step)+1
  endif
  mesh_data%vion_meshsz=max(mesh_data%sph_meshsz,mesh_data%core_meshsz,ii1)
 endif

!Mesh for vbare should be inside augmentation region
!For compatibility with other codes, could put it to the same mesh
!as ionic potential
!mesh_data%vbare_meshsz=mesh_data%vion_meshsz
 mesh_data%vbare_meshsz=mesh_data%sph_meshsz

 mesh_data%prj_msz_max=four*mesh_data%sph_meshsz

!=== Build mesh definitions ===

 mesh_data%meshtp=-1;mesh_data%meshsz=0
 mesh_data%radstp=zero;mesh_data%logstp=zero

!Partial waves
 mesh_data%nmesh=1
 mesh_data%iwavmesh=1
 mesh_data%meshtp(1)=mesh_data%mesh_type
 mesh_data%meshsz(1)=mesh_data%wav_meshsz
 mesh_data%radstp(1)=mesh_data%rad_step
 mesh_data%logstp(1)=mesh_data%log_step
!Projectors
 if (mesh_data%wav_meshsz/=mesh_data%prj_meshsz) then
  mesh_data%nmesh=mesh_data%nmesh+1
  mesh_data%iprjmesh=mesh_data%nmesh
  mesh_data%meshtp(mesh_data%nmesh)=mesh_data%mesh_type
  mesh_data%meshsz(mesh_data%nmesh)=mesh_data%prj_meshsz
  mesh_data%radstp(mesh_data%nmesh)=mesh_data%rad_step
  mesh_data%logstp(mesh_data%nmesh)=mesh_data%log_step
 else
  mesh_data%iprjmesh=mesh_data%iwavmesh
 endif
!Core density
 if (mesh_data%wav_meshsz/=mesh_data%core_meshsz) then
  if (mesh_data%prj_meshsz/=mesh_data%core_meshsz) then
   mesh_data%nmesh=mesh_data%nmesh+1
   mesh_data%icoremesh=mesh_data%nmesh
   mesh_data%meshtp(mesh_data%nmesh)=mesh_data%mesh_type
   mesh_data%meshsz(mesh_data%nmesh)=mesh_data%core_meshsz
   mesh_data%radstp(mesh_data%nmesh)=mesh_data%rad_step
   mesh_data%logstp(mesh_data%nmesh)=mesh_data%log_step
  else
   mesh_data%icoremesh=mesh_data%iprjmesh
  endif
 else
  mesh_data%icoremesh=mesh_data%iwavmesh
 endif
!Local ionic potential
 if (mesh_data%wav_meshsz/=mesh_data%vion_meshsz) then
  if(mesh_data%prj_meshsz/=mesh_data%vion_meshsz) then
   if(mesh_data%core_meshsz/=mesh_data%vion_meshsz) then
    mesh_data%nmesh=mesh_data%nmesh+1
    mesh_data%ivionmesh=mesh_data%nmesh
    mesh_data%meshtp(mesh_data%nmesh)=mesh_data%mesh_type
    mesh_data%meshsz(mesh_data%nmesh)=mesh_data%vion_meshsz
    mesh_data%radstp(mesh_data%nmesh)=mesh_data%rad_step
    mesh_data%logstp(mesh_data%nmesh)=mesh_data%log_step
   else
    mesh_data%ivionmesh=mesh_data%icoremesh
   endif
  else
   mesh_data%ivionmesh=mesh_data%iprjmesh
  endif
 else
  mesh_data%ivionmesh=mesh_data%iwavmesh
 endif
!Local vbare potential
 if (mesh_data%wav_meshsz/=mesh_data%vbare_meshsz) then
  if(mesh_data%prj_meshsz/=mesh_data%vbare_meshsz) then
   if(mesh_data%core_meshsz/=mesh_data%vbare_meshsz) then
    if(mesh_data%vion_meshsz/=mesh_data%vbare_meshsz) then
     mesh_data%nmesh=mesh_data%nmesh+1
     mesh_data%ivbaremesh=mesh_data%nmesh
     mesh_data%meshtp(mesh_data%nmesh)=mesh_data%mesh_type
     mesh_data%meshsz(mesh_data%nmesh)=mesh_data%vbare_meshsz
     mesh_data%radstp(mesh_data%nmesh)=mesh_data%rad_step
     mesh_data%logstp(mesh_data%nmesh)=mesh_data%log_step
    else
     mesh_data%ivbaremesh=mesh_data%ivionmesh
    endif
   else
    mesh_data%ivbaremesh=mesh_data%icoremesh
   endif
  else
   mesh_data%ivbaremesh=mesh_data%iprjmesh
  endif
 else
  mesh_data%ivbaremesh=mesh_data%iwavmesh
 endif
!Valence density
 if (mesh_data%wav_meshsz/=mesh_data%vale_meshsz) then
  if(mesh_data%prj_meshsz/=mesh_data%vale_meshsz) then
   if(mesh_data%core_meshsz/=mesh_data%vale_meshsz) then
    if(mesh_data%vion_meshsz/=mesh_data%vale_meshsz) then
     if(mesh_data%vbare_meshsz/=mesh_data%vale_meshsz) then
      mesh_data%nmesh=mesh_data%nmesh+1
      mesh_data%ivalemesh=mesh_data%nmesh
      mesh_data%meshtp(mesh_data%nmesh)=mesh_data%mesh_type
      mesh_data%meshsz(mesh_data%nmesh)=mesh_data%vale_meshsz
      mesh_data%radstp(mesh_data%nmesh)=mesh_data%rad_step
      mesh_data%logstp(mesh_data%nmesh)=mesh_data%log_step
     else
      mesh_data%ivalemesh=mesh_data%ivbaremesh
     endif
    else
     mesh_data%ivalemesh=mesh_data%ivionmesh
    endif
   else
    mesh_data%ivalemesh=mesh_data%icoremesh
   endif
  else
   mesh_data%ivalemesh=mesh_data%iprjmesh
  endif
 else
  mesh_data%ivalemesh=mesh_data%iwavmesh
 endif

 end subroutine build_mesh_data

!!=================================================================
!! NAME
!! destroy_mesh_data
!!
!! FUNCTION
!! Free all the memory in a mesh_data datastructure
!!
!! ISIDE EFFECTS
!!  mesh_data
!!   Data defining various meshes
!!
!! PARENTS
!!
!!=================================================================

subroutine destroy_mesh_data(mesh_data)

 type(mesh_data_type),intent(inout) :: mesh_data

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (associated(mesh_data%meshtp)) then
   deallocate(mesh_data%meshtp)
 end if
 if (associated(mesh_data%meshsz)) then
   deallocate(mesh_data%meshsz)
 end if
 if (associated(mesh_data%radstp)) then
   deallocate(mesh_data%radstp)
 end if
 if (associated(mesh_data%logstp)) then
   deallocate(mesh_data%logstp)
 end if

 end subroutine destroy_mesh_data

!!=================================================================
END Module XMLInterface
