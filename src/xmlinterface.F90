!! NAME
!!  XMLInterface
!!
!! FUNCTION
!!  This module contains routines used to convert PAW atomic dataset into a XML file
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ATOMPAW group (NHolzwarth, MTorrent, FJollet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains the following active subroutines:
!     Atompaw2XML, rdinputxml, opt_proj_rso, aamat, xmloutput, xmlprtcore,
!        read_inputstring, simpson_int, gauleg, build_mesh_data,
!        destroy_mesh_data, get_xc_data, get_xc_alias
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

Module XMLInterface

 use io_tools
 use Tools
 use GlobalMath
 use atomdata
 use gridmod
 use excor
 use pseudo
 use interpolation_mod
 use pkginfo
 use libxc_mod
 use input_dataset_mod

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
 character*10 :: atompaw2xmlver='3.3.0', verdate='dec. 2020'

!Default unit for XML file(s)
 integer, parameter :: unit_xml=22,unit_xml_core=23

!We decide to cut at r=10 u.a because of numeric noise...
 real(dp),parameter :: rmax_vloc=10._dp

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
  integer :: itaumesh            ! Index of mesh for kinetic energy core density
  integer :: ivionmesh           ! Index of mesh for vion potential
  integer :: ivbaremesh          ! Index of mesh for vbare potential
  integer :: ivlda12mesh         ! Index of mesh for LDA-1/2 potential
  integer :: ivalemesh           ! Index of mesh for valence density
  integer :: wav_meshsz          ! Size of mesh for partial waves
  integer :: sph_meshsz          ! Size of mesh for partial waves
  integer :: prj_meshsz          ! Size of mesh for projectors
  integer :: core_meshsz         ! Size of mesh for core density
  integer :: tau_meshsz          ! Size of mesh for kinetic energy core density
  integer :: vion_meshsz         ! Size of mesh for vion potential
  integer :: vbare_meshsz        ! Size of mesh for vbare potential
  integer :: vlda12_meshsz       ! Size of mesh for LDA-1/2 potential
  integer :: vale_meshsz         ! Size of mesh for valence density
  integer :: prj_msz_max         ! Maximum size for projector (used for RSO)
  real(dp) :: rad_step           ! Default value for radial step
  real(dp) :: log_step           ! Default value for log step
  integer,allocatable :: meshtp(:)  ! Array storing mesh type for all meshes
  integer,allocatable :: meshsz(:)  ! Array storing mesh size for all meshes
  real(dp),allocatable :: radstp(:) ! Array storing radial step for all meshes
  real(dp),allocatable :: logstp(:) ! Array storing log step for all meshes
 end type mesh_data_type

!Options for REAL SPACE OPTIMIZATION of projectors
 type pawrso_type
  logical :: userso           ! TRUE if Real Space Optimization is required
  real(dp) :: ecut            ! Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
  real(dp) :: gfact           ! Real Space Optimization parameter: Gamma/Gmax ratio
  real(dp) :: werror          ! Real Space Optimization parameter: max. error W_l allowed
 end type pawrso_type

!Options for LDA-1/2 METHOD
 type pawlda12_type
  logical :: uselda12            ! TRUE if LDA-1/2 potential calculation is required
  integer :: orb_l               ! LDA-1/2 parameter: l quantum number of the ionized orbital
  integer :: orb_n               ! LDA-1/2 parameter: n quantum number of the ionized orbital
  real(dp) :: ion                ! LDA-1/2 parameter: amount of charge removed from the ionized orbital
  real(dp) :: rcut               ! LDA-1/2 parameter: cut-off radius (in bohr)
  character(15) :: logfile       ! LDA-1/2 parameter: name of the log file
  real(dp),allocatable :: pot(:) ! LDA-1/2 parameter: local potential used to apply LDA-1/2 method
 end type pawlda12_type


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
!! rdinputxml,build_mesh_data,opt_proj_rso,xmloutput,xmlprtcore
!!=================================================================

 subroutine Atompaw2XML(AEOrbit,AEPot,AESCF,PAW,FC,Grid)

 type(OrbitInfo), intent(in)     :: AEOrbit
 type(PotentialInfo), intent(in) :: AEPot
 type (SCFInfo),  intent(in)     :: AESCF
 type(Pseudoinfo), intent(in)    :: PAW
 type(FCInfo), intent(in)        :: FC
 type(GridInfo), intent(in)      :: Grid

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: id,vlocopt,nsplgrid,err
 type(mesh_data_type) :: mesh_data
 type(pawrso_type) :: pawrso
 type(pawlda12_type) :: pawlda12
 logical :: print_corewf,print_xmlfile
 character*(fnlen) :: author,comment,file_xml,file_xml_core,xcname
 character*(10000) :: input_string
 real(dp),allocatable :: tproj(:,:)
 real(dp),pointer :: vion(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 write(std_out,'(/,a)') '========================================================'
 write(std_out,'(3a)')  '==   Atompaw2XML - v',trim(atompaw2xmlver),&
&                                                '                             =='
 write(std_out,'(a)')   '==   Writes a PAW dataset generated by "Atompaw"      =='
 write(std_out,'(a)')   '==   into a XML file...                               =='
 write(std_out,'(a)')   '========================================================'

 print_xmlfile=.true.
 print_corewf=.true.

 call read_inputstring(input_string)

!------------------------------------------------------------------
!---- Read choices from input file

 call rdinputxml(vlocopt,print_corewf,pawrso,pawlda12,author,comment,nsplgrid,input_string)

!Do not try to print XML dataset when it is not possible
 print_xmlfile=.true.
 if (diracrelativistic) then
  write(std_out,'(/,2x,a)') 'XML dataset output not compatible with Dirac relativistic approach!'
  write(std_out,'(2x,a)')   'XML dataset file will not be created.'
  if (print_corewf) write(std_out,'(2x,a)')   'Only core orbitals file will be created.'
  print_xmlfile=.false.
 end if

!------------------------------------------------------------------
!---- Build radial meshes definitions

 call build_mesh_data(mesh_data,Grid,PAW%irc,PAW%ivion,PAW%ivale,PAW%coretailpoints,PAW%itau)

!------------------------------------------------------------------
!---- Several checks on then densities/potentials

 if (print_xmlfile) then

! Check that the pseudo valence density is zero at grid end
  if (mesh_data%vale_meshsz>0) then
   write(std_out,'(/,2x,a,/2x,a,f5.2,a,a,g11.4)') 'Atompaw2XML info:',&
&    '  At r=',Grid%r(mesh_data%vale_meshsz),' a.u.,',&
&    ' Pseudo_valence_density= ', PAW%tden(mesh_data%vale_meshsz) &
&                                /Grid%r(mesh_data%vale_meshsz)**2/(4*pi)
   write(std_out,'(2x,a)') '  This quantity must be as small as possible.'
  endif

! Check that the local potential is zero at grid end
  if (vlocopt==1.or.vlocopt==2) then
   if (vlocopt==1) vion => PAW%abinitvloc
   if (vlocopt==2) vion => PAW%abinitnohat
   write(std_out,'(/,2x,a,/,2x,a,f5.2,a,a,g11.4)') 'Atompaw2XML info:',&
&     '  At r_vloc=',Grid%r(mesh_data%vion_meshsz),' a.u.,',&
&   ' VHartree(ntild(Zv+Zc))= -Zv/r + ', half*vion(mesh_data%vion_meshsz) &
&                      +(AEPot%nz-FC%zcore)/Grid%r(mesh_data%vion_meshsz)
   write(std_out,'(2x,a)') '  This quantity must be as small as possible.'
  endif

! Check that (pseudo-core + nhat) density is positive 
  if (vlocopt==1.and..not.PAW%poscorenhat) then 
   write(std_out,'(/,2x,a)') 'Detected negative values for pseudo core + nhat'
   write(std_out,'(2x,a)')   '  which is incompatible with usexcnhat.'
   write(std_out,'(2x,a)')   'Please try reducing rc_core.'
   write(std_out,'(2x,a)')   'XML file not created!'
   print_xmlfile=.false.
  endif

 endif

!------------------------------------------------------------------
!---- Possibly optimize projectors

 if (print_xmlfile) then
  allocate(tproj(mesh_data%prj_msz_max,PAW%nbase))
  tproj=zero
  do id=1,PAW%nbase
   tproj(1:mesh_data%wav_meshsz,id)=PAW%otp(1:mesh_data%wav_meshsz,id)
  end do
  call opt_proj_rso(tproj,mesh_data,pawrso,Grid,PAW,err)
  print_xmlfile=(print_xmlfile.and.err==0)
 endif

!------------------------------------------------------------------
!---- Possibly compute LDA-1/2 potential

 if (print_xmlfile) then
  if (pawlda12%uselda12) call lda12(pawlda12,mesh_data%vlda12_meshsz,Grid,PAW)
 endif

!------------------------------------------------------------------
!---- Write PAW dataset file in XML format

 if (print_xmlfile) then
  xcname=exctype;if (have_libxc) call libxc_getshortname(exctype,xcname)
  file_xml=TRIM(AEpot%sym)//'.'//TRIM(xcname)//'-paw.xml'
  call xmloutput(file_xml,Grid,AESCF,AEPot,FC,PAW,mesh_data,tproj,vlocopt,&
 &               input_string,author,comment,nsplgrid,pawlda12)
 endif

!------------------------------------------------------------------
!---- Write core wave functions in XML format

 if (print_corewf) then
  xcname=exctype;if (have_libxc) call libxc_getshortname(exctype,xcname)
  file_xml_core=TRIM(AEpot%sym)//'.'//TRIM(xcname)//'-paw.corewf.xml'
  call xmlprtcore(file_xml_core,Grid,AEOrbit,AEPot,FC,mesh_data,&
 &                input_string,author,comment)
 endif

!------------------------------------------------------------------
!---- End
 if (allocated(tproj)) deallocate(tproj)
 if (allocated(pawlda12%pot)) deallocate(pawlda12%pot)
 call destroy_mesh_data(mesh_data)

 write(std_out,'(/,2x,a)') 'Atompaw2XML ended.'

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
!!  prtcorewf= option for the printing of core wave functions
!!  pawrso
!!    %userso=TRUE if REAL Space Optimization is required
!!    %ecut=Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
!!    %gfact=Real Space Optimization parameter: Gamma/Gmax ratio
!!    %werror=Real Space Optimization parameter: max. error W_l allowed
!!  pawlda12
!!    %uselda12=TRUE if LDA-1/2 potential calculation is required
!!    %orb_l=LDA-1/2 parameter: l quantum number of the ionized orbital
!!    %orb_n=LDA-1/2 parameter: n quantum number of the ionized orbital
!!    %ion=LDA-1/2 parameter: amount of charge removed from the ionized orbital
!!    %rcut=LDA-1/2 parameter: cut-off radius (in bohr)
!!    %logfile=LDA-1/2 parameter: name of the log file
!!  author= name of the author(s)
!!  comment= additional comment line (usually table version)
!!  nsplgrid=if >0, size of a (reduced) grid on which interpolate all data in XML file
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

 subroutine rdinputxml(vlocopt,prtcorewf,pawrso,pawlda12,author,comment,nsplgrid,input_string)

 integer,intent(out)            :: vlocopt,nsplgrid
 logical,intent(out)            :: prtcorewf
 type(pawrso_type),intent(out)  :: pawrso
 type(pawlda12_type),intent(out):: pawlda12
 character(len=*),intent(out)   :: author,comment
 character(len=*),intent(inout) :: input_string

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 logical,parameter :: print_options=.false. !Already done in input_dataset_mod
 integer :: i_author,i_comment,len_author,len_comment
 logical :: usexcnhat
 character(200) :: xmlstr,xmlstr_tmp

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!Read XML options from standard input
 CALL input_dataset_read_xml(xml_string=xmlstr)

!Option for the PRinTing of CORE Wave Functions
 prtcorewf=input_dataset%xml_prtcorewf
 if (print_options) then
  if (prtcorewf) then
   write(std_out,'(a,2x,a)') ch10,'Printing of core wave functions file'
  else
   write(std_out,'(a,2x,a)') ch10,'No printing of core wave functions file'
  end if
 end if

!Option for use of NHAT in XC
 usexcnhat=input_dataset%xml_usexcnhat
 if (usexcnhat) then
  vlocopt=1
  if (print_options) then
   write(std_out,'(a,2x,a)') ch10,'Use of compensation charge in XC terms'
   write(std_out,'(a,2x,a)') ch10,'Zero potential and Blochl local ionic potential will be output in XML file'
  end if
 else
  vlocopt=2
  if (print_options) then
   write(std_out,'(a,2x,a)') ch10,'No use of compensation charge in XC terms'
   write(std_out,'(a,2x,a)') ch10,'Zero potential and Kresse local ionic potential will be output in XML file'
  end if
 end if

!Options for use of REAL SPACE OPTIMIZATION
 pawrso%userso=input_dataset%xml_userso
 pawrso%ecut=input_dataset%xml_rso_ecut
 pawrso%gfact=input_dataset%xml_rso_gfact
 pawrso%werror=input_dataset%xml_rso_werror
 if (print_options) then
  if (pawrso%userso) then
   write(std_out,'(a,2x,a,f6.2,a,f6.2,a,g8.2,a)') ch10,&
&    'Real Space optim.: Ecut, Gamma/Gmax, Wl(error) = [',&
&    pawrso%ecut,', ',pawrso%gfact,', ',pawrso%werror,']'
  else
   write(std_out,'(a,2x,a)') ch10,'No Real Space Optimization of projectors'
  end if
 end if

!Options for use of LDA-1/2 METHOD
 pawlda12%uselda12=input_dataset%xml_uselda12
 pawlda12%orb_l=input_dataset%xml_lda12_orb_l
 pawlda12%orb_n=input_dataset%xml_lda12_orb_n
 pawlda12%ion=input_dataset%xml_lda12_ion
 pawlda12%rcut=input_dataset%xml_lda12_rcut
 pawlda12%logfile=trim(input_dataset%xml_lda12_logfile)
 if (print_options) then
  if (pawlda12%uselda12) then
   write(std_out,'(a,2x,2(a,i1),2(a,f5.2),a)') ch10,&
&    'LDA-1/2 pot. calculation: orbital(n,l), ionization, rcut = [(',&
&    pawlda12%orb_n,',',pawlda12%orb_l,'), ',pawlda12%ion,', ',pawlda12%rcut,']'
   write(std_out,'(3a)') &
&    '  See ''',trim(input_dataset%xml_lda12_logfile),''' file to check convergence.'
  else
   write(std_out,'(a,2x,a)') ch10,'No computation of LDA-1/2 potential'
  end if
 end if

!Option for spline on a reduced log. grid
 nsplgrid=-1;if (input_dataset%xml_usespl) nsplgrid=input_dataset%xml_spl_meshsz
 if (print_options) then
  if (nsplgrid>0) then
   write(std_out,'(a,2x,a,i5,a)') ch10,&
&    'Spline on a reduced log grid with ',nsplgrid, ' points'
  else
   write(std_out,'(a,2x,a)') ch10,'No spline on a reduced log. grid'
  end if
 end if

!Author to be mentioned in XML file header
 author=trim(input_dataset%xml_author)

!Comment line to be added in XMLfile header
 comment=trim(input_dataset%xml_comment)

!Update input string (without author name and comment)
 i_author=index(xmlstr,'AUTHOR')
 i_comment=index(xmlstr,'COMMENT')
 xmlstr_tmp=xmlstr
 if (i_author>0) then
  len_author=len(trim(author))
  if (i_author==1) xmlstr_tmp=trim(xmlstr_tmp(i_author+len_author+10:))
  if (i_author >1) xmlstr_tmp=xmlstr_tmp(1:i_author-1)//trim(xmlstr_tmp(i_author+len_author+10:))
 end if
 if (i_comment>0) then
  len_comment=len(trim(comment))
  if (i_comment==1) xmlstr_tmp=trim(xmlstr_tmp(i_author+len_comment+10:))
  if (i_comment >1) xmlstr_tmp=xmlstr_tmp(1:i_comment-1)//trim(xmlstr_tmp(i_author+len_comment+10:))
 end if
 write(unit=input_string,fmt='(6a)') trim(input_string),char(10),"XMLOUT",char(10),trim(xmlstr_tmp)

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
!! OUTPUT
!!  err=error code (0=OK)
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

 subroutine opt_proj_rso(tproj,mesh_data,pawrso,Grid,PAW,err)

 type(mesh_data_type),intent(inout) :: mesh_data
 type(pawrso_type),intent(in)    :: pawrso
 type(Gridinfo),intent(in) :: Grid
 type(Pseudoinfo),intent(in) ::  PAW
 real(dp),intent(inout) :: tproj(mesh_data%prj_msz_max,PAW%nbase)
 integer,intent(out) :: err

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

 err=0
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
   err=1;exit
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

 if (err==0) then
  mesh_data%prj_meshsz=r0_meshsz
  write(std_out,'(4x,2(a,f7.4),a)') 'New radius R0 for projectors (a.u.)=',&
&                          r0,' (=',r0/rc_prj,'*Rc(proj))'
  if (r0>1.55_dp*rc_prj) &
&   write(std_out,'(4x,a,3(/,4x,a))') &
&                 'Warning:',&
&                 '  Radius for projectors (R0) seems to be high !',&
&                 '  You should change parameters of Real Space Optimization',&
&                 '  (increase Ecut, Gamma/Gmax or Wl).'
 end if

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
!! lda12
!!
!! FUNCTION
!! compute the potential for LDA-1/2 method
!!
!! INPUTS
!!  pawlda12
!!    %uselda12=TRUE if LDA-1/2 potential calculation is required
!!    %orb_l=LDA-1/2 parameter: l quantum number of the ionized orbital
!!    %orb_n=LDA-1/2 parameter: n quantum number of the ionized orbital
!!    %ion=LDA-1/2 parameter: amount of charge removed from the ionized orbital
!!    %rcut=LDA-1/2 parameter: cut-off radius (in bohr)
!!    %logfile=LDA-1/2 parameter: name of the log file
!!  vlda12_meshsz=size of output pawlda12%pot potential
!!  Grid=grid info from atompaw
!!  PAW=PAW datastructure from atompaw
!!
!! OUTPUT
!!  pawlda12
!!    %pot(:)=LDA-1/2 parameter: local potential used to apply LDA-1/2 method
!!
!! PARENTS
!! atompaw2xml
!!
!!=================================================================
 subroutine lda12(pawlda12,vlda12_meshsz,Grid,PAW)

 integer,intent(in) :: vlda12_meshsz
 type(pawlda12_type),intent(inout)  :: pawlda12
 type(Gridinfo),intent(in) :: Grid
 type(Pseudoinfo),intent(in) ::  PAW

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer,parameter :: unlog=55
 integer :: io,irc,n,unstdout
 logical :: found
 type(OrbitInfo)     :: wkOrbit
 type(PotentialInfo) :: wkPot
 type(SCFInfo)       :: wkSCF
 type(FCInfo)        :: wkFC
 real(dp),allocatable :: cut_off(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (.not.pawlda12%uselda12) return

!Redirect standard output
 open(unlog,file=trim(pawlda12%logfile),form='formatted')
 unstdout=STD_OUT
 STD_OUT=unlog

!Save current status of atomic AE data
 call CopyOrbit(AEOrbit,wkOrbit)
 call CopyPot(AEPot,wkPot)
 call CopySCF(AESCF,wkSCF)
 call CopyFC(FC,wkFC)
 
!Change occupations (ionize one orbital)
 found=.false.
 do io=1,AEOrbit%norbit
   if (AEOrbit%l(io)==pawlda12%orb_l.and.AEOrbit%np(io)==pawlda12%orb_n) then
     found=.true.
     AEOrbit%occ(io)=AEOrbit%occ(io)-pawlda12%ion
   end if
 end do
 if (.not.found) then
   write(unstdout,'(/,2x,a)') &
&   'Error in LDA-1/2 configuration: LDA-1/2 potential computation aborted!'
   pawlda12%uselda12=.false.
 else

!  Compute new AE atomic data
   call SCFatom('NC',.false.,skip_reading=.true.)

!  Compute LDA-1/2 potential
   n=min(vlda12_meshsz,size(AEPot%rv))
   irc=FindGridIndex(Grid,pawlda12%rcut)
   allocate(pawlda12%pot(n))
   allocate(cut_off(n))
   pawlda12%pot(:)=zero
   cut_off(1:n)=(one-(Grid%r(1:n)/pawlda12%rcut)**8)**3
   cut_off(1:n)=half*(cut_off(1:n)+abs(cut_off(1:n)))
   if (irc<n) cut_off(irc+1:n)=zero
   pawlda12%pot(1:n)=(AEPot%rv(1:n)-PAW%AErefrv(1:n))*cut_off(1:n)*sqrt(four*pi)*half
   pawlda12%pot(2:n)=pawlda12%pot(2:n)/Grid%r(2:n)
   call extrapolate(Grid,pawlda12%pot(1:n))
   deallocate(cut_off)

!  Restore initial AE atomic data
   call CopyOrbit(wkOrbit,AEOrbit)
   call CopyPot(wkPot,AEPot)
   call CopySCF(wkSCF,AESCF)
   call CopyFC(wkFC,FC)
 end if

!Release temporary memory
 call DestroyOrbit(wkOrbit)
 call DestroyPot(wkPot)
 call DestroyFC(wkFC)

!Restore stdout
 close(unlog)
 STD_OUT=unstdout

 end subroutine lda12


!!=================================================================
!! NAME
!! xmloutput
!!
!! FUNCTION
!! Write the PAW data file in XML format
!!
!! INPUTS
!! fname=file name (with .xml suffixe)
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
!! author= string containing the author(s) name
!! comment= additional comment line to be printed in the header (usually table version)
!! nsplgrid=if >0, size of a (reduced) grid on which interpolate all data in XML file
!! pawlda12
!!    %uselda12=TRUE if LDA-1/2 potential calculation is required
!!    %rcut=LDA-1/2 parameter: cut-off radius (in bohr)
!!    %pot(:)=LDA-1/2 parameter: local potential used to apply LDA-1/2 method
!!
!! PARENTS
!! atompaw2xml
!!
!! CHILDREN
!!
!!=================================================================

 subroutine xmloutput(fname,Grid,AESCF,AEPot,FC,PAW,mesh_data,tproj,&
&                     vlocopt,input_string,author,comment,nsplgrid,pawlda12)

 integer,intent(in) :: vlocopt,nsplgrid
 character(len=*),intent(in) :: input_string,author,comment,fname
 TYPE(Gridinfo),intent(in) :: Grid
 TYPE (SCFInfo),intent(in) :: AESCF
 TYPE(Potentialinfo),intent(in) :: AEPot
 TYPE (FCInfo),intent(in) :: FC
 TYPE(Pseudoinfo),intent(in) :: PAW
 TYPE(pawlda12_type),intent(in)  :: pawlda12
 TYPE(mesh_data_type),intent(inout) :: mesh_data
 real(dp),intent(in) :: tproj(mesh_data%prj_msz_max,PAW%nbase)

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 real(dp), parameter :: tol_zero=1.d-50 ! Threshold below which quantities are zero
 
 integer :: ib,ic,ii,n,n_aux,ir,irc_aux,meshsz,meshsz_aux,meshst_aux,nmesh,OK
 integer :: mesh_start(mesh_data%nmesh),mesh_size(mesh_data%nmesh)
 logical :: extra1
 character(len=4) :: char4
 character(len=5) :: char5a,char5b,xc_type
 character(len=20) :: char20
 character(len=132) :: xc_name,xcname_short,code_name
 real(dp) :: sqr4pi,rad,radstp0,logstp0,radstp_spl,logstp_spl
 character(len=3) :: gridt(mesh_data%nmesh)
 real(dp),allocatable :: dum(:),dum_aux(:),rad_aux(:),dudr(:),rveff_aux(:)
 real(dp),allocatable :: phi_aux(:,:),tphi_aux(:,:),proj_aux(:,:)
 TYPE(Gridinfo) :: Grid1

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!Some defs
 sqr4pi=sqrt(4*pi)
 n=Grid%n ; allocate(dum(n))
 extra1=.false.

!In a change of grid has been requested, determine new grid data
 radstp_spl=-1.d0 ; logstp_spl=-1.d0
 if (nsplgrid>0) then
   write(std_out,'(/,2x,a,/,2x,a,i5,a)') 'Atompaw2XML info:',&
&   '  All quantities will be interpolated on a ',nsplgrid,'-point log. grid.'
   n_aux=nsplgrid
   logstp_spl=0.02d0
   call findh(AEPot%zz,Grid%r(mesh_data%meshsz(1)-1),nsplgrid,logstp_spl,radstp_spl)
   irc_aux=int(tol8+log(1.d0+PAW%rc/radstp_spl)/logstp_spl)+1
   radstp_spl= PAW%rc/(exp(logstp_spl*(irc_aux-1))-1.d0)
   extra1=(Grid%r(mesh_data%meshsz(1))<radstp_spl*(exp(logstp_spl*(nsplgrid-1))-1.d0))
   if (extra1) then
     n_aux=n_aux+1
     logstp_spl=0.02d0
     call findh(AEPot%zz,Grid%r(mesh_data%meshsz(1)-1),nsplgrid+1,logstp_spl,radstp_spl)
     irc_aux=int(tol8+log(1.d0+PAW%rc/radstp_spl)/logstp_spl)+1
     radstp_spl= PAW%rc/(exp(logstp_spl*(irc_aux-1))-1.d0)
   end if
   allocate(rad_aux(n_aux),dum_aux(n_aux))
   call InitGrid(Grid1,logstp_spl,Grid%range,r0=radstp_spl,do_not_print=.true.)
   rad_aux(1:nsplgrid)=Grid1%r(1:nsplgrid)
 else
   allocate(rad_aux(n),dum_aux(n))
   rad_aux(1:n)=Grid%r(1:n)
   irc_aux=PAW%irc
 end if

!Open file for writing
 OPEN(unit_xml,file=TRIM(fname),form='formatted')

!Write XML header
 WRITE(unit_xml,'("<?xml  version=""1.0""?>")')
 WRITE(unit_xml,'("<paw_dataset version=""0.7"">")')
 WRITE(unit_xml,'("<!-- PAW-XML specification: http://esl.cecam.org/Paw-xml -->")')

!Write title
 WRITE(unit_xml,'(/,"<!-- PAW atomic dataset for ",a," -->")') trim(ADJUSTL(AEPot%sym))

!Write additional comment line (usually table version)
 if (trim(comment)/="") WRITE(unit_xml,'("<!-- ",a," -->")') trim(comment)

!Write Atompaw information
 WRITE(unit_xml,'(/,"<!-- Atompaw ",a)') atp_version
 WRITE(unit_xml,'("  Contact info: Natalie Holzwarth")')
 WRITE(unit_xml,'("  email: natalie@wfu.edu, web: pwpaw.wfu.edu")')
 WRITE(unit_xml,'("  Energy units=Hartree, length units=bohr")')
 CALL PrintDate(unit_xml, '  PAW functions generated on ')
 if (trim(author)/="") WRITE(unit_xml,'(a,a)') '  by ',trim(author)
 WRITE(unit_xml,'("  The input file is available at the end of this file")')
 WRITE(unit_xml,'("-->",/)')

!Write atom definition
 WRITE(unit=char5a,fmt='(f5.2)') AEPot%zz
 WRITE(unit_xml,'("<atom symbol=""",a,""" Z=""",a,$)') &
&   trim(ADJUSTL(AEPot%sym)),trim(ADJUSTL(char5a))
 WRITE(unit=char5a,fmt='(f5.2)') FC%zcore
 WRITE(unit=char5b,fmt='(f5.2)') AEPot%nz-FC%zcore
 WRITE(unit_xml,'(""" core=""",a,""" valence=""",a,"""/>")') &
&   trim(ADJUSTL(char5a)),trim(ADJUSTL(char5b))

!Write XC definition
 call get_xc_data(xc_type,xc_name)
 if (have_libxc.and.xc_type/="UNKNOWN") then
   call libxc_getshortname(xc_name,xcname_short)
   call get_xc_alias(xcname_short,xc_name)
 endif
 WRITE(unit_xml,'("<xc_functional type=""",a,""" name=""",a,"""/>")') &
&      TRIM(xc_type),TRIM(xc_name)
 code_name="atompaw-"//atp_version

!Generator data
 if (scalarrelativistic) then
   WRITE(unit_xml,'("<generator type=""scalar-relativistic"" name=""",a,""" orthogonalisation=""", a,"""/>")')&
&               TRIM(code_name),trim(PAW%orthogonalization_scheme)
 else if (diracrelativistic) then
   WRITE(unit_xml,'("<generator type=""dirac-relativistic"" name=""",a,""" orthogonalisation=""", a,"""/>")')&
&               TRIM(code_name),trim(PAW%orthogonalization_scheme)
 else
   WRITE(unit_xml,'("<generator type=""non-relativistic"" name=""",a,""" orthogonalisation=""", a,"""/>")')&
&               TRIM(code_name),trim(PAW%orthogonalization_scheme)
 endif

!Energies
 WRITE(unit_xml,'("<ae_energy kinetic=""",1pe25.17,""" xc=""",1pe25.17,"""")') &
&      AESCF%ekin/2,AESCF%eexc/2
 WRITE(unit_xml,'("  electrostatic=""",1pe25.17,""" total=""",1pe25.17,"""/>")') &
&      AESCF%estatic/2,AESCF%etot/2
 WRITE(unit_xml,'("<core_energy kinetic=""",1pe25.17,"""/>")') AESCF%corekin*0.5d0

!PAW radius
 WRITE(unit_xml,'("<paw_radius rc=""",f17.14,"""/>")') match_on_splgrid(PAW%rc)

!Electronic configuration
 WRITE(unit_xml,'("<valence_states>")')
 do ib=1,PAW%nbase
   call mkname(ib,char4)
   char20=stripchar('"'//AEPot%sym//char4//'"')
   ii=min(ABS(PAW%np(ib)),100)
   if (ii<100) then
     if(diracrelativistic) then
       WRITE(unit_xml,'("  <state n=""",i2,""" l=""",i1,""" kappa=""",i2,""" f=""",1pe14.7,$)')&
&          ii,PAW%l(ib),PAW%kappa(ib),AEOrbit%occ(ib)
     else
       WRITE(unit_xml,'("  <state n=""",i2,""" l=""",i1,""" f=""",1pe14.7,$)')&
&          ii,PAW%l(ib),PAW%occ(ib)
     end if
     WRITE(unit_xml,'(""" rc=""",f13.10,""" e=""",1pe14.7,""" id=",a6,"/>")')&
&        match_on_splgrid(PAW%rcio(ib)),PAW%eig(ib)*0.5d0,TRIM(char20)
   else
     if(diracrelativistic) then
       WRITE(unit_xml,'("  <state        l=""",i1,""" kappa=""",i2,$)') PAW%l(ib),PAW%kappa(ib)
     else
       WRITE(unit_xml,'("  <state        l=""",i1,$)') PAW%l(ib)
     end if
     WRITE(unit_xml,'("""                    rc=""",f13.10,""" e=""",1pe14.7,""" id=",a6,"/>")')&
&        match_on_splgrid(PAW%rcio(ib)),PAW%eig(ib)*0.5d0,TRIM(char20)
   end if
 enddo
 WRITE(unit_xml,'("</valence_states>")')

!Radial meshes definitions
 nmesh=mesh_data%nmesh
 if(maxval(mesh_data%meshtp(1:mesh_data%nmesh))==minval(mesh_data%meshtp(1:mesh_data%nmesh))) then
   mesh_data%meshsz(1:mesh_data%nmesh)=maxval(mesh_data%meshsz(1:mesh_data%nmesh))
   nmesh=1
   mesh_data%icoremesh=1
   mesh_data%itaumesh=1
   mesh_data%iprjmesh=1
   mesh_data%iwavmesh=1
   mesh_data%ivionmesh=1
   mesh_data%ivalemesh=1
   mesh_data%ivbaremesh=1
   mesh_data%ivlda12mesh=1
 endif
 if (nmesh>1.and.nsplgrid>0) stop '  Bug (1) in xmlinterface: nmesh>1 and nsplgrid>0!'

 do ii=1,nmesh
  allocate(dudr(mesh_data%meshsz(ii)),stat=OK)

  select case(mesh_data%meshtp(ii))

   case(1)
    char20='r=d*i'
    gridt(ii)="lin"
    mesh_start(ii)=1
    mesh_size(ii)=mesh_data%meshsz(ii)
    radstp0=zero
    logstp0=mesh_data%radstp(ii)
    dudr(1:mesh_data%meshsz(ii))=logstp0

   case(2)
    char20='r=a*(exp(d*i)-1)'
    gridt(ii)="log"
    mesh_start(ii)=1
    if (nsplgrid<=0) then
      mesh_size(ii)=mesh_data%meshsz(ii)
      radstp0=mesh_data%radstp(ii)
      logstp0=mesh_data%logstp(ii)
    else
      mesh_size(ii)=nsplgrid
      radstp0=radstp_spl
      logstp0=logstp_spl
    end if
    dudr(1:mesh_size(ii))=logstp0*(radstp0+rad_aux(1:mesh_size(ii)))

   case default
    stop '  Bug (2) in xmlinterface: mesh type not implemented in Atompaw!'
  end select

  WRITE(unit_xml,'("<radial_grid eq=""",a,""" a=""",es23.16,$)') trim(char20),radstp0
  WRITE(unit_xml,'(""" d=""",es23.16,""" istart=""0"" iend=""",i5,$)') logstp0,mesh_size(ii)-mesh_start(ii)
  WRITE(unit_xml,'(""" id=""",a,i1,""">")') gridt(ii),ii
  WRITE(unit_xml,'("  <values>")')
  WRITE(unit_xml,'(3(1x,es23.16))') (rad_aux(ir),ir=mesh_start(ii),mesh_size(ii))
  WRITE(unit_xml,'("  </values>")')
  WRITE(unit_xml,'("  <derivatives>")')
  WRITE(unit_xml,'(3(1x,es23.16))') (dudr(ir),ir=mesh_start(ii),mesh_size(ii))
  WRITE(unit_xml,'("  </derivatives>")')
  WRITE(unit_xml,'("</radial_grid>")')
  deallocate(dudr)

 end do

!Compensation charge shape function
 if (gaussianshapefunction) then
   WRITE(unit_xml,'("<shape_function type=""gauss"" rc=""",f19.16,"""/>")') match_on_splgrid(PAW%gausslength)
 else if (besselshapefunction) then
   WRITE(unit_xml,'("<shape_function type=""bessel"" rc=""",f19.16,"""/>")') match_on_splgrid(PAW%rc_shap)
 else
   WRITE(unit_xml,'("<shape_function type=""sinc"" rc=""",f19.16,"""/>")') match_on_splgrid(PAW%rc_shap)
 endif

!Core densities
 meshsz=mesh_data%meshsz(mesh_data%icoremesh)
 meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%icoremesh),extra1)
 meshst_aux=mesh_start(mesh_data%icoremesh)
 dum(2:meshsz)=sqr4pi*FC%coreden(2:meshsz)/(4*pi*Grid%r(2:meshsz)**2)
 call extrapolate(Grid,dum(1:meshsz))
 call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
 WRITE(unit_xml,'("<ae_core_density grid=""",a,i1,""" rc=""",f19.16,""">")') &
& gridt(mesh_data%icoremesh),mesh_data%icoremesh,match_on_splgrid(PAW%rc_core)
 WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
 WRITE(unit_xml,'("</ae_core_density>")')
 dum(2:meshsz)=sqr4pi*PAW%tcore(2:meshsz)/(4*pi*Grid%r(2:meshsz)**2)
 call extrapolate(Grid,dum(1:meshsz))
 call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
 WRITE(unit_xml,'("<pseudo_core_density grid=""",a,i1,""" rc=""",f19.16,""">")') &
& gridt(mesh_data%icoremesh),mesh_data%icoremesh,match_on_splgrid(PAW%rc_core)
 WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
 WRITE(unit_xml,'("</pseudo_core_density>")')

!Valence density
 meshsz=mesh_data%meshsz(mesh_data%ivalemesh)
 meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%ivalemesh),extra1)
 meshst_aux=mesh_start(mesh_data%ivalemesh)
 dum(2:meshsz)=sqr4pi*PAW%tden(2:meshsz)/(4*pi*Grid%r(2:meshsz)**2)
 call extrapolate(Grid,dum(1:meshsz))
 call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
 rad=maxval(PAW%rcio(1:PAW%nbase))
 WRITE(unit_xml,'("<pseudo_valence_density grid=""",a,i1,""" rc=""",f19.16,""">")') &
& gridt(mesh_data%ivalemesh),mesh_data%ivalemesh,match_on_splgrid(rad)
 WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
 WRITE(unit_xml,'("</pseudo_valence_density>")')

!Kinetic energy core densities
!Temporarily not available in scalar relativistic
 if (.not.scalarrelativistic) then
  meshsz=mesh_data%meshsz(mesh_data%itaumesh)
  meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%itaumesh),extra1)
  meshst_aux=mesh_start(mesh_data%itaumesh)
  dum(2:meshsz)= sqr4pi*FC%coretau(2:meshsz)/(4*pi*Grid%r(2:meshsz)**2)
  call extrapolate(Grid,dum(1:meshsz))
  call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
  WRITE(unit_xml,'("<ae_core_kinetic_energy_density grid=""",a,i1,""" rc=""",f19.16,""">")') &
&  gridt(mesh_data%itaumesh),mesh_data%itaumesh,match_on_splgrid(PAW%rc_core)
  WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
  WRITE(unit_xml,'("</ae_core_kinetic_energy_density>")')
  dum(2:meshsz)= sqr4pi*PAW%tcoretau(2:meshsz)/(4*pi*Grid%r(2:meshsz)**2)
  call extrapolate(Grid,dum(1:meshsz))
  call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
  WRITE(unit_xml,'("<pseudo_core_kinetic_energy_density grid=""",a,i1,""" rc=""",f19.16,""">")') &
&  gridt(mesh_data%itaumesh),mesh_data%itaumesh,match_on_splgrid(PAW%rc_core)
  WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
  WRITE(unit_xml,'("</pseudo_core_kinetic_energy_density>")')
 else
  write(std_out,'(5(/,2x,a))') 'Atompaw2XML WARNING!!!!!',&
&   '  Kinetic energy core density is not available',&
&   '    within scalar relativistic scheme!',&
&   '  Will not be present in the XML dataset.',&
&   '  This is temporary, sorry!'
 end if

!Vbare potential
 meshsz=mesh_data%meshsz(mesh_data%ivbaremesh)
 meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%ivbaremesh),extra1)
 meshst_aux=mesh_start(mesh_data%ivbaremesh)
 dum(1:meshsz)=sqr4pi*half*PAW%vloc(1:meshsz)
 call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
 dum_aux(irc_aux:meshsz_aux)=0.d0 ! Vbare has to be zero at rc
 WRITE(unit_xml,'("<zero_potential grid=""",a,i1,""" rc=""",f19.16,""">")') &
& gridt(mesh_data%ivbaremesh),mesh_data%ivbaremesh,match_on_splgrid(PAW%rc)
 WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
 WRITE(unit_xml,'("</zero_potential>")')

!Local ionic potential
 if (vlocopt==1) then
  meshsz=mesh_data%meshsz(mesh_data%ivionmesh)
  meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%ivionmesh),extra1)
  meshst_aux=mesh_start(mesh_data%ivionmesh)
  dum(1:meshsz)=sqr4pi*half*PAW%abinitvloc(1:meshsz)
  call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
   WRITE(unit_xml,'("<kresse_joubert_local_ionic_potential grid=""",a,i1,""" rc=""",f19.16,""">")') &
&   gridt(mesh_data%ivionmesh),mesh_data%ivionmesh,match_on_splgrid(PAW%rc)
   WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
   WRITE(unit_xml,'("</kresse_joubert_local_ionic_potential>")')
  end if

!Local Blochl''s potential
 if(vlocopt==2) then
  meshsz=mesh_data%meshsz(mesh_data%ivionmesh)
  meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%ivionmesh),extra1)
  meshst_aux=mesh_start(mesh_data%ivionmesh)
  dum(1:meshsz)=sqr4pi*half*PAW%abinitnohat(1:meshsz)
  call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
  WRITE(unit_xml,'("<blochl_local_ionic_potential grid=""",a,i1,""" rc=""",f19.16,""">")') &
&  gridt(mesh_data%ivionmesh),mesh_data%ivionmesh,match_on_splgrid(PAW%rc)
  WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
  WRITE(unit_xml,'("</blochl_local_ionic_potential>")')
 endif

!Local LDA-1/2 potential
 if (pawlda12%uselda12) then
  meshsz=mesh_data%meshsz(mesh_data%ivlda12mesh)
  meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%ivlda12mesh),extra1)
  meshst_aux=mesh_start(mesh_data%ivlda12mesh)
  dum(1:meshsz)=pawlda12%pot(1:meshsz)
  call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
  ir=index_on_splgrid(pawlda12%rcut)
  dum_aux(ir+1:meshsz_aux)=0.d0 ! V_lda12 has to be zero for r>rc
    WRITE(unit_xml,'("<LDA_minus_half_potential grid=""",a,i1,""" rc=""",f19.16,""">")') &
&   gridt(mesh_data%ivlda12mesh),mesh_data%ivlda12mesh,match_on_splgrid(pawlda12%rcut)
  WRITE(unit_xml,'(3(1x,es23.16))') (dum_aux(ii),ii=meshst_aux,meshsz_aux)
  WRITE(unit_xml,'("</LDA_minus_half_potential>")')
 end if

!Partial waves and projectors
!-- Partial waves
 meshsz=mesh_data%meshsz(mesh_data%iwavmesh)
 meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%iwavmesh),extra1)
 meshst_aux=mesh_start(mesh_data%iwavmesh)
 allocate(phi_aux(meshsz_aux-meshst_aux+1,PAW%nbase))
 allocate(tphi_aux(meshsz_aux-meshst_aux+1,PAW%nbase))
 Do ib=1,PAW%nbase
   dum(2:meshsz)=PAW%ophi(2:meshsz,ib)/Grid%r(2:meshsz)
   call extrapolate(Grid,dum(1:meshsz))
   call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
   phi_aux(1:meshsz_aux,ib)=dum_aux(1:meshsz_aux)
   dum(2:meshsz)=PAW%otphi(2:meshsz,ib)/Grid%r(2:meshsz)
   call extrapolate(Grid,dum(1:meshsz))
   call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
   if(nsplgrid>0) dum_aux(irc_aux:meshsz_aux)=phi_aux(irc_aux:meshsz_aux,ib)
   tphi_aux(meshst_aux:meshsz_aux,ib)=dum_aux(meshst_aux:meshsz_aux)
 Enddo
!-- Projectors
 meshsz=mesh_data%meshsz(mesh_data%iprjmesh)
 meshsz_aux=merge(nsplgrid,mesh_size(mesh_data%iprjmesh),extra1)
 meshst_aux=mesh_start(mesh_data%iprjmesh)
 allocate(proj_aux(meshsz_aux-meshst_aux+1,PAW%nbase))
 Do ib=1,PAW%nbase
   dum(2:meshsz)=tproj(2:meshsz,ib)/Grid%r(2:meshsz)
   call extrapolate(Grid,dum(1:meshsz))
   call interp_and_filter(dum(1:meshsz),dum_aux(1:meshsz_aux))
   if(nsplgrid>0) dum_aux(irc_aux:meshsz_aux)=zero
   proj_aux(meshst_aux:meshsz_aux,ib)=dum_aux(meshst_aux:meshsz_aux)
 Enddo
!-- In case of a spline, re-orthogonalize projectors
 if (nsplgrid >0) then
   Do ib=1,PAW%nbase
     tphi_aux(1:meshsz_aux,ib)=tphi_aux(1:meshsz_aux,ib)*Grid1%r(1:meshsz_aux)
     proj_aux(1:meshsz_aux,ib)=proj_aux(1:meshsz_aux,ib)*Grid1%r(1:meshsz_aux)
   Enddo
   call vdborth(irc_aux,tphi_aux,proj_aux)
   Do ib=1,PAW%nbase
     tphi_aux(2:meshsz_aux,ib)=tphi_aux(2:meshsz_aux,ib)/Grid1%r(2:meshsz_aux)
     call extrapolate(Grid1,tphi_aux(1:meshsz_aux,ib))
     proj_aux(2:meshsz_aux,ib)=proj_aux(2:meshsz_aux,ib)/Grid1%r(2:meshsz_aux)
     call extrapolate(Grid1,proj_aux(1:meshsz_aux,ib))
   Enddo
 end if
!-- Writing
 Do ib=1,PAW%nbase
   call mkname(ib,char4)
   char20=stripchar('"'//AEPot%sym//char4//'"')
   WRITE(unit_xml,'("<ae_partial_wave state=",a6," grid=""",a,i1,""">")') &
&   TRIM(char20),gridt(mesh_data%iwavmesh),mesh_data%iwavmesh
   WRITE(unit_xml,'(3(1x,es23.16))') (phi_aux(ii,ib),ii=meshst_aux,meshsz_aux)
   WRITE(unit_xml,'("</ae_partial_wave>")')
   WRITE(unit_xml,'("<pseudo_partial_wave state=",a6," grid=""",a,i1,""">")')&
&   TRIM(char20),gridt(mesh_data%iwavmesh),mesh_data%iwavmesh
   WRITE(unit_xml,'(3(1x,es23.16))') (tphi_aux(ii,ib),ii=meshst_aux,meshsz_aux)
   WRITE(unit_xml,'("</pseudo_partial_wave>")')
   WRITE(unit_xml,'("<projector_function state=",a6," grid=""",a,i1,""">")') &
&   TRIM(char20),gridt(mesh_data%iprjmesh),mesh_data%iprjmesh
   WRITE(unit_xml,'(3(1x,es23.16))') (proj_aux(ii,ib),ii=meshst_aux,meshsz_aux)
   WRITE(unit_xml,'("</projector_function>")')
   !do ic=1,PAW%nbase
   !  write(std_out,*) "splined", ib,ic
   !  dum(meshst_aux:meshsz_aux)=tphi_aux(meshst_aux:meshsz_aux,ib)*proj_aux(meshst_aux:meshsz_aux,ic)*dudr(meshst_aux:meshsz_aux)
   !  dum(meshst_aux:meshsz_aux)=dum(meshst_aux:meshsz_aux)*rad_aux(meshst_aux:meshsz_aux)*rad_aux(meshst_aux:meshsz_aux)
   !  xx=overint(meshsz_aux-meshst_aux+1,logstp0,dum(meshst_aux:meshsz_aux),-1)
   !  write(std_out,*) "ORTHO", ib,xx/logstp0
   !end do
 Enddo
!-- Release memory
 deallocate(phi_aux)
 deallocate(proj_aux)
 deallocate(tphi_aux)

!Kinetic terms
 WRITE(unit_xml,'("<kinetic_energy_differences>")')
 WRITE(unit_xml,'(3(1x,es23.16))') ((PAW%kij(ib,ic)/2,ic=1,PAW%nbase),ib=1,PAW%nbase)
 WRITE(unit_xml,'("</kinetic_energy_differences>")')

!Core-valence exchange terms
 WRITE(unit_xml,'("<exact_exchange_X_matrix>")')
 WRITE(unit_xml,'(3(1x,es23.16))') ((PAW%TXVC(ib,ic)/2,ic=1,PAW%nbase),ib=1,PAW%nbase)
 WRITE(unit_xml,'("</exact_exchange_X_matrix>")')

!Core-core exchange terms
 WRITE(unit_xml,'("<exact_exchange core-core=""", 1x,es23.16,"""/>")') &
&     PAW%XCORECORE/2 

! Input file
 WRITE(unit_xml,'("<!-- Program:  atompaw - input data follows: ")')
 WRITE(unit_xml,'(a)') trim(input_string)
 WRITE(unit_xml,'(a)') "END"
 WRITE(unit_xml,'(" Program:  atompaw - input end -->")')
 WRITE(unit_xml,'("</paw_dataset>")')

!Close file and end
 CLOSE(unit_xml)
 WRITE(STD_OUT,'(/,2x,a)') 'XML atomic dataset created.'

 if(nsplgrid>0) call destroygrid(Grid1)
 deallocate(rad_aux)
 deallocate(dum_aux)
 deallocate(dum)
 
 CONTAINS
 !**************************************************
 ! If an interpolation on an auxiliary grid is
 !   requested, give index of input radius
 !  on this auxiliary grid
 !**************************************************
  integer function index_on_splgrid(input_radius)
    real(dp),intent(in) :: input_radius
    if (nsplgrid>0) then
      if (logstp_spl>0.d0) then
        index_on_splgrid=int(tol8+log(1.d0+input_radius/radstp_spl)/logstp_spl)+1
      else if (radstp_spl>0.d0) then
        index_on_splgrid=int(tol8+input_radius/radstp_spl)+1
      end if
    else
      index_on_splgrid=FindGridIndex(Grid,input_radius)
    end if
  end function index_on_splgrid

 !**************************************************
 ! If an interpolation on an auxiliary grid is
 !   requested, match input radius on this auxiliary
 !   grid (defined by radstp_spl, logstp_spl)
 !**************************************************
  real(dp) function match_on_splgrid(input_radius)
    real(dp),intent(in) :: input_radius
    integer :: indx
    match_on_splgrid=input_radius
    if (nsplgrid>0) then
      if (logstp_spl>0.d0) then
        indx=int(tol8+log(1.d0+input_radius/radstp_spl)/logstp_spl)+1
        match_on_splgrid=radstp_spl*(exp(logstp_spl*(indx-1))-1.d0)
      else if (radstp_spl>0.d0) then
        indx=int(tol8+input_radius/radstp_spl)+1
        match_on_splgrid=radstp_spl*(indx-1)
      end if
    end if
  end function match_on_splgrid

 !**************************************************
 ! If an interpolation on an auxiliary grid is
 !   requested, inerpolate an input function on
 !   this auxiliary grid.
 ! Also filter the input function (put zero below
 !   a given threshold)
 !**************************************************
  subroutine interp_and_filter(func_in,func_out)
  real(dp),intent(in) :: func_in(:)
  real(dp),intent(out) :: func_out(:)
  integer :: jj,msz_in,msz_out,msz_spl
  real(dp) :: l1,l2,l3,x1,x2,x3,xx
  logical :: extra
  msz_in=size(func_in) ; msz_out=size(func_out)
  func_out=zero
  if (nsplgrid>0) then
    if (msz_out/=nsplgrid) stop '  Bug (1) in interp_and_filter: msz_out/=nsplgrid!'
    extra=(Grid%r(msz_in)<Grid1%r(msz_out))
    msz_spl=merge(msz_out-1,msz_out,extra)
    call interpfunc(msz_in,Grid%r,func_in,msz_spl,Grid1%r,func_out)
  else
    if (msz_out>msz_in) stop '  Bug (2) in interp_and_filter: msz_out>msz_in!'
    func_out(1:msz_out)=func_in(1:msz_out)
  end if
  do jj=1,msz_out
    if (abs(func_out(jj))<tol_zero) func_out(jj)=0.d0
  end do

  end subroutine interp_and_filter 

!***********************************************************************
!* In case of interpolation on an auxiliary grid, a
!* reorthonomalisation of projector and pseudo wavefunctions
!* is necessary. This is done thanks to a Vanderbilt orthonormalisation
!************************************************************************
  subroutine vdborth(irc_aux,tphi_aux,proj_aux)

  real(dp),intent(in) :: tphi_aux(:,:)
  real(dp),intent(inout) :: proj_aux(:,:)
  integer :: irc_aux

  integer :: i,icount,io,irc,j,jo,l,lmax,nbase
  real(dp), allocatable :: aa(:,:),ai(:,:),omap(:),proj_aux1(:,:)

  allocate(proj_aux1(size(proj_aux,1),PAW%nbase))
  lmax=PAW%lmax
  nbase=PAW%nbase
  irc=irc_aux
     do l=0,lmax
       icount=0
       do io=1,nbase
        if (PAW%l(io)==l) icount=icount+1
       enddo
       if (icount==0) cycle
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       aa=0;icount=0
       do io=1,nbase
        if (PAW%l(io)==l) then
          icount=icount+1
          omap(icount)=io
        endif
       enddo
       do i=1,icount
         io=omap(i)
         do j=1,icount
           jo=omap(j)
           aa(i,j)=overlap(Grid1,tphi_aux(:,io),proj_aux(:,jo),1,irc)
         enddo
       enddo
       ai=aa;call minverse(ai,icount,icount,icount)

       do i=1,icount
         io=omap(i)
         proj_aux1(:,io)=0
         do j=1,icount
           jo=omap(j)
           proj_aux1(:,io)=proj_aux1(:,io)+proj_aux(:,jo)*ai(j,i)
         enddo
       enddo
       deallocate(aa,ai,omap)
     enddo
     proj_aux=proj_aux1
     deallocate(proj_aux1)
  end subroutine vdborth 

 END SUBROUTINE xmloutput

!!=================================================================
!! NAME
!! xmlprtcore
!!
!! FUNCTION
!! Write the core wave functions in XML format
!!
!! INPUTS
!! fname=file name (with .xml suffixe)
!! Grid= Grid datastructure from atompaw
!! AEOrbit= AEOrbit datastructure from atompaw
!! AEPot= AEPOT datastructure from atompaw
!! FC= FC datastructure from atompaw
!! mesh_data= datatructure containing the definition of
!!            all the radial meshes used in the XML file
!! input_string= string containing a copy of atompaw input file
!! author= string containing the author(s) name
!! comment= additional comment line to be printed in the header (usually table version)
!!
!! PARENTS
!! atompaw2xml
!!
!! CHILDREN
!!
!!=================================================================

 SUBROUTINE xmlprtcore(fname,Grid,AEOrbit,AEPot,FC,mesh_data,input_string,author,comment)

 character(len=fnlen),intent(in) :: fname
 TYPE (FCInfo),intent(in) :: FC
 TYPE(mesh_data_type),intent(in) :: mesh_data
 type(OrbitInfo),intent(in) :: AEOrbit
 TYPE(Potentialinfo),intent(in) :: AEPot
 TYPE(Gridinfo),intent(in) :: Grid
 character(len=*),intent(in) :: input_string,author,comment

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: core_size,corewf_meshsz,ib,icor,ii,ir,isppol,nmesh,nsppol,OK
 real(dp) :: radstp0,logstp0,sqr4pi
 character(len=3) :: spstrg,gridt(mesh_data%nmesh)
 character(len=4) :: char4
 character(len=5) :: char5a,xc_type
 character(len=132) :: xc_name,xcname_short
 character(len=20) :: char20,char21
 character(len=1) :: char_orb(4)
 integer,allocatable :: irwf(:)
 real(dp),allocatable :: dum(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!Hard-coded values, spin unrestricted
!Spinors or collinear magnetism not yet supported
 nsppol=1
 char_orb(1)="s";char_orb(2)="p";char_orb(3)="d";char_orb(4)="f"
!Open file for writing
 OPEN(unit_xml_core,file=TRIM(fname),form='formatted')

!Write XML header
 WRITE(unit_xml_core,'("<?xml  version=""1.0""?>")')
 WRITE(unit_xml_core,'("<paw_setup version=""0.7"">")')

!Write title
 WRITE(unit_xml_core,'(/,"<!-- All-electron core wavefunctions for ",a," -->")') trim(ADJUSTL(AEPot%sym))

!Write additional comment line (usually table version)
 if (trim(comment)/="") WRITE(unit_xml_core,'("<!-- ",a," -->")') trim(comment)

!Write Atompaw information
 WRITE(unit_xml_core,'(/,"<!-- Atompaw ",a)') atp_version
 CALL PrintDate(unit_xml_core, '  Core orbitals generated on ')
 if (trim(author)/="") WRITE(unit_xml_core,'(a,a)') '  by ',trim(author)
 WRITE(unit_xml_core,'("  Energy units=Hartree, length units=bohr")')
 WRITE(unit_xml_core,'("  The input file is available at the end of this file")')
 WRITE(unit_xml_core,'("-->",/)')

!Write atom definition
 WRITE(unit=char5a,fmt='(f5.2)') AEPot%zz
 WRITE(unit_xml_core,'("<atom symbol=""",a,""" Z=""",a,$)') &
&      trim(ADJUSTL(AEPot%sym)),trim(ADJUSTL(char5a))
 WRITE(unit=char5a,fmt='(f5.2)') FC%zcore
 WRITE(unit_xml_core,'(""" core=""",a,"""/>")')trim(ADJUSTL(char5a))

!Write XC definition
 call get_xc_data(xc_type,xc_name)
 if (have_libxc.and.xc_type/="UNKNOWN") then
   call libxc_getshortname(xc_name,xcname_short)
   call get_xc_alias(xcname_short,xc_name)
 endif
 WRITE(unit_xml_core,'("<xc_functional type=""",a,""" name=""",a,"""/>")') &
&      TRIM(xc_type),TRIM(xc_name)

!Generator data
 if (scalarrelativistic) then
   WRITE(unit_xml_core,'("<generator type=""scalar-relativistic"" name=""atompaw"">")')
 else if (diracrelativistic) then
   WRITE(unit_xml_core,'("<generator type=""dirac-relativistic"" name=""atompaw"">")')
 else
   WRITE(unit_xml_core,'("<generator type=""non-relativistic"" name=""atompaw"">")')
 endif
 WRITE(unit_xml_core,'("</generator>)")')

!Number of core orbitals (not needed)
 core_size=count(AEOrbit%iscore(1:AEOrbit%norbit))
!WRITE(unit_xml_core,'("<orbitals norbs=""",i2,"""/>")') core_size

!Read mesh size 
 ii=0;allocate(irwf(core_size));irwf(:)=mesh_data%meshsz(mesh_data%iwavmesh)
 do ib=1,AEOrbit%norbit
  if (AEOrbit%iscore(ib)) then
   ii=ii+1;ir=mesh_data%meshsz(mesh_data%iwavmesh)+1
   do while (ir>1)
    ir=ir-1
    if (abs(AEOrbit%wfn(ir,ib))>tol10) then
     irwf(ii)=min(ir+1,mesh_data%meshsz(mesh_data%iwavmesh));ir=1
    end if
   end do
  end if
 end do
 corewf_meshsz=maxval(irwf(1:core_size))
 deallocate(irwf)
 corewf_meshsz=maxval(mesh_data%meshsz(1:mesh_data%nmesh))
!Electronic configuration
 WRITE(unit_xml_core,'("<core_states>")')
 icor=0
 do ib=1,AEOrbit%norbit
   if (AEOrbit%iscore(ib)) then
     icor=icor+1;if (icor>core_size) stop '  Bug (1) in xmlprtcore!'
     call mkname(icor,char4)
     char20=stripchar('"'//AEPot%sym//'_core'//char4//'"')
     if(diracrelativistic) then
       WRITE(unit_xml_core,'("  <state n=""",i2,""" l=""",i1,""" kappa=""",i2,""" f=""",1pe14.7,$)')&
&          AEOrbit%np(ib),AEOrbit%l(ib),AEOrbit%kappa(ib),AEOrbit%occ(ib)
     else
       WRITE(unit_xml_core,'("  <state n=""",i2,""" l=""",i1,""" f=""",1pe14.7,$)')&
&          AEOrbit%np(ib),AEOrbit%l(ib),AEOrbit%occ(ib)
     endif
     WRITE(unit_xml_core,'("""  e=""",1pe14.7,""" id=",a11,"/>")')&
&        AEOrbit%eig(ib)*0.5d0,TRIM(char20)
   end if
 enddo
 WRITE(unit_xml_core,'("</core_states>")')

!Radial meshes definitions
 nmesh=1
 do ii=1,nmesh
  allocate(dum(corewf_meshsz))
  select case(mesh_data%meshtp(ii))
   case(1)
    char21='r=d*i'
    gridt(ii)="lin"
    radstp0=zero
    logstp0=mesh_data%radstp(ii)
    dum(:)=logstp0
   case(2)
    char21='r=a*(exp(d*i)-1)'
    gridt(ii)="log"
    radstp0=mesh_data%radstp(ii)
    logstp0=mesh_data%logstp(ii)
    dum(1:corewf_meshsz)=logstp0*(radstp0+Grid%r(1:corewf_meshsz)) 
   case default
    stop '  Bug (2) in xmlprtcore: mesh type not implemented in Atompaw!'
  end select
  WRITE(unit_xml_core,'("<radial_grid eq=""",a,""" a=""",es23.16,$)') trim(char21),radstp0
  WRITE(unit_xml_core,'(""" d=""",es23.16,""" istart=""0"" iend=""",i5,$)') &
&  logstp0,corewf_meshsz-1
  WRITE(unit_xml_core,'(""" id=""",a,i1,"""/>")') gridt(ii),ii
  WRITE(unit_xml_core,'("  <values>")')
  WRITE(unit_xml_core,'(3(1x,es23.16))') (Grid%r(ir),ir=1,corewf_meshsz)
  WRITE(unit_xml_core,'("  </values>")')
  WRITE(unit_xml_core,'("  <derivatives>")')
  WRITE(unit_xml_core,'(3(1x,es23.16))') (dum(ir),ir=1,corewf_meshsz)
  WRITE(unit_xml_core,'("  </derivatives>")')
  deallocate(dum)
 end do

!Write the core wave functions
 icor=0
 do ib=1,AEOrbit%norbit
  if (AEOrbit%iscore(ib)) then
    icor=icor+1;if (icor>core_size) stop 'Bug (3) in xmlprtcore!'
    call mkname(AEOrbit%np(ib),char4)
    char20=stripchar('"'//AEPot%sym//char4//char_orb(AEOrbit%l(ib)+1)//'"')
    WRITE(unit_xml_core,'("<ae_core_wavefunction state=",a6," grid=""",a,i1,""">")') &
&     TRIM(char20),gridt(mesh_data%iwavmesh),mesh_data%iwavmesh
    ALLOCATE(dum(mesh_data%meshsz(mesh_data%iwavmesh)),stat=OK)
    dum=zero
    dum(2:mesh_data%meshsz(mesh_data%iwavmesh))= &
&               AEOrbit%wfn(2:mesh_data%meshsz(mesh_data%iwavmesh),ib) &
&              /Grid%r(2:mesh_data%meshsz(mesh_data%iwavmesh))
    call extrapolate(Grid,dum)
    WRITE(unit_xml_core,'(3(1x,es23.16))') (dum(ir),ir=1,corewf_meshsz)
    DEALLOCATE(dum)
    WRITE(unit_xml_core,'("</ae_core_wavefunction>")')
  end if ! if icore
 end do   !ib

!Write the core lwave functions
 if(diracrelativistic) then
   icor=0
   do ib=1,AEOrbit%norbit
    if (AEOrbit%iscore(ib)) then
      icor=icor+1;if (icor>core_size) stop '  Bug (4) in xmlprtcore!'
      call mkname(AEOrbit%np(ib),char4)
      char20=stripchar('"'//AEPot%sym//char4//char_orb(AEOrbit%l(ib)+1)//'"')
      WRITE(unit_xml_core,'("<ae_core_lwavefunction state=",a6," grid=""",a,i1,""">")') &
&       TRIM(char20),gridt(mesh_data%iwavmesh),mesh_data%iwavmesh
      ALLOCATE(dum(mesh_data%meshsz(mesh_data%iwavmesh)),stat=OK)
      dum=zero
      dum(2:mesh_data%meshsz(mesh_data%iwavmesh))= &
&                 AEOrbit%lwfn(2:mesh_data%meshsz(mesh_data%iwavmesh),ib) &
&                /Grid%r(2:mesh_data%meshsz(mesh_data%iwavmesh))
      call extrapolate(Grid,dum)
      WRITE(unit_xml_core,'(3(1x,es23.16))') (dum(ir),ir=1,corewf_meshsz)
      DEALLOCATE(dum)
      WRITE(unit_xml_core,'("</ae_core_lwavefunction>")')
    end if ! if icore
   end do   !ib
 end if ! diracrelativistic

!Echi input file content
 WRITE(unit_xml_core,'("<!-- Program:  atompaw - input data follows: ")')
 WRITE(unit_xml_core,'(a)') trim(input_string)
 WRITE(unit_xml_core,'(a)') "END"
 WRITE(unit_xml_core,'(" Program:  atompaw - input end -->")')
 WRITE(unit_xml_core,'("</paw_setup>")')

!Close the file
 close(unit_xml_core)
 WRITE(STD_OUT,'(/,2x,a)') 'XML core orbitals file created.'

 end subroutine xmlprtcore


!!=================================================================
!! NAME
!! read_inputstring
!!
!! FUNCTION
!! Read the file echoing the atompaw input file
!! and transfer it into a character string
!!
!! OUTPUT
!!  input_string=character string containing the file
!!
!! PARENTS
!  xml2abinit
!!=================================================================

 subroutine read_inputstring(input_string)

 character(len=*) :: input_string

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: OK,input_unit
 character(len=132) :: inputline

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 open(input_unit,file='dummy',form='formatted')
 read(input_unit,'(a)',iostat=OK,end=10) inputline
 if (OK/=0) return
 input_string=trim(inputline)
 do
   read(input_unit,'(a)',iostat=OK,end=10) inputline
   if (OK/=0) exit
   write(unit=input_string,fmt='(3a)') trim(input_string),char(10),trim(inputline)
 enddo
 return
10 continue
 close(input_unit)

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

subroutine build_mesh_data(mesh_data,Grid,irc,ivion,ivale,coretailpoints,itau)

 type(mesh_data_type),intent(out) :: mesh_data
 type(GridInfo),intent(in) :: Grid
 integer,intent(in) :: irc,ivion,ivale,coretailpoints,itau

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
 if (usingloggrid(Grid)) then
  mesh_data%wav_meshsz=Grid%n
  mesh_data%sph_meshsz=min(1+nint(log(one+Grid%r(irc)/mesh_data%rad_step)/mesh_data%log_step),mesh_data%wav_meshsz)
 else
  mesh_data%wav_meshsz=irc+Grid%ishift
  mesh_data%sph_meshsz=min(1+nint(Grid%r(irc)/mesh_data%rad_step),mesh_data%wav_meshsz)
 end if
 mesh_data%prj_meshsz=mesh_data%sph_meshsz  ! To be modified by RSO

 mesh_data%core_meshsz=coretailpoints
 mesh_data%vale_meshsz=ivale
 mesh_data%tau_meshsz=itau

 mesh_data%vion_meshsz=ivion
 if (mesh_data%vion_meshsz<=0) then
  if (mesh_data%mesh_type==1) then
   ii1=int(one+rmax_vloc/mesh_data%rad_step)
  else
   ii1=int(log(one+rmax_vloc/mesh_data%rad_step)/mesh_data%log_step)+1
  endif
  mesh_data%vion_meshsz=max(mesh_data%sph_meshsz,mesh_data%core_meshsz,ii1)
 endif

 mesh_data%vlda12_meshsz=max(mesh_data%vale_meshsz,mesh_data%vion_meshsz)

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
!Kinetic energy density
 if (mesh_data%wav_meshsz/=mesh_data%tau_meshsz) then
  if(mesh_data%prj_meshsz/=mesh_data%tau_meshsz) then
   if(mesh_data%core_meshsz/=mesh_data%tau_meshsz) then
    if(mesh_data%vion_meshsz/=mesh_data%tau_meshsz) then
     if(mesh_data%vbare_meshsz/=mesh_data%tau_meshsz) then
      if(mesh_data%vale_meshsz/=mesh_data%tau_meshsz) then
       mesh_data%nmesh=mesh_data%nmesh+1
       mesh_data%itaumesh=mesh_data%nmesh
       mesh_data%meshtp(mesh_data%nmesh)=mesh_data%mesh_type
       mesh_data%meshsz(mesh_data%nmesh)=mesh_data%tau_meshsz
       mesh_data%radstp(mesh_data%nmesh)=mesh_data%rad_step
       mesh_data%logstp(mesh_data%nmesh)=mesh_data%log_step
      else
       mesh_data%itaumesh=mesh_data%ivalemesh
      endif
     else
      mesh_data%itaumesh=mesh_data%ivbaremesh
     endif
    else
     mesh_data%itaumesh=mesh_data%ivionmesh
    endif
   else
    mesh_data%itaumesh=mesh_data%icoremesh
   endif
  else
   mesh_data%itaumesh=mesh_data%iprjmesh
  endif
 else
  mesh_data%itaumesh=mesh_data%iwavmesh
 endif
!LDA-1/2 potential
 if (mesh_data%wav_meshsz/=mesh_data%vlda12_meshsz) then
  if(mesh_data%prj_meshsz/=mesh_data%vlda12_meshsz) then
   if(mesh_data%core_meshsz/=mesh_data%vlda12_meshsz) then
    if(mesh_data%vion_meshsz/=mesh_data%vlda12_meshsz) then
     if(mesh_data%vbare_meshsz/=mesh_data%vlda12_meshsz) then
      if(mesh_data%vale_meshsz/=mesh_data%vlda12_meshsz) then
       if(mesh_data%tau_meshsz/=mesh_data%vlda12_meshsz) then 
        mesh_data%nmesh=mesh_data%nmesh+1
        mesh_data%ivlda12mesh=mesh_data%nmesh
        mesh_data%meshtp(mesh_data%nmesh)=mesh_data%mesh_type
        mesh_data%meshsz(mesh_data%nmesh)=mesh_data%vlda12_meshsz
        mesh_data%radstp(mesh_data%nmesh)=mesh_data%rad_step
        mesh_data%logstp(mesh_data%nmesh)=mesh_data%log_step
       else
        mesh_data%ivlda12mesh=mesh_data%itaumesh
       endif
      else
       mesh_data%ivlda12mesh=mesh_data%ivalemesh
      endif
     else
      mesh_data%ivlda12mesh=mesh_data%ivbaremesh
     endif
    else
     mesh_data%ivlda12mesh=mesh_data%ivionmesh
    endif
   else
    mesh_data%ivlda12mesh=mesh_data%icoremesh
   endif
  else
   mesh_data%ivlda12mesh=mesh_data%iprjmesh
  endif
 else
  mesh_data%ivlda12mesh=mesh_data%iwavmesh
 endif

 end subroutine build_mesh_data

!!=================================================================
!! NAME
!! destroy_mesh_data
!!
!! FUNCTION
!! Free all the memory in a mesh_data datastructure
!!
!! SIDE EFFECTS
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

 if (allocated(mesh_data%meshtp)) deallocate(mesh_data%meshtp)
 if (allocated(mesh_data%meshsz)) deallocate(mesh_data%meshsz)
 if (allocated(mesh_data%radstp)) deallocate(mesh_data%radstp)
 if (allocated(mesh_data%logstp)) deallocate(mesh_data%logstp)

 end subroutine destroy_mesh_data


!!=================================================================
!! NAME
!! get_xc_data
!!
!! FUNCTION
!! Get XC data in a suitable form for XML printing
!!
!! INPUTS
!!  exctype= string containing XC type
!!
!! OUTPUT
!!  xc_name= name of XC functional
!!  xc_type= LDA or GGA
!!
!! PARENTS
!!  xmlout,xmlprtcore
!!
!!=================================================================

subroutine get_xc_data(xctype,xcname)

 character(len=*),intent(out) :: xctype,xcname

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (trim(ADJUSTL(exctype))=="LDA-PW") then
   xctype="LDA"
   xcname="PW"
 elseif (trim(ADJUSTL(exctype))=="GGA-PBE") then
   xctype="GGA"
   xcname="PBE"
 elseif (trim(ADJUSTL(exctype))=="GGA-PBESOL") then
   xctype="GGA"
   xcname="PBESOL"
 else if (have_libxc) then
   if (libxc_isgga()) then
     xctype="GGA"
   else
     xctype="LDA"
   end if
   xcname=trim(exctype)
!else ! old coding
!  xctype="LIBXC"
!  write(xcname,'(i6)') -pshead%pspxc_abinit
 else
   write(std_out,'(2x,a)') 'Warning in xmlinterface: unknown XC type!'
   xctype="UNKNOWN"
   xcname="UNKNOWN"

 end if

 end subroutine get_xc_data


!!=================================================================
!! NAME
!! get_xc_alias
!!
!! FUNCTION
!! Get XC name alias (following PAW-XML specification)
!!   from a libXC functional name
!!
!! INPUTS
!!  xc_name= string containing long XC name
!!
!! OUTPUT
!!  xc_alias= alias of XC functional
!!
!! PARENTS
!!  xmlout,xmlprtcore
!!
!!=================================================================

subroutine get_xc_alias(xc_name,xc_alias)

 character(len=*),intent(in) :: xc_name
 character(len=*),intent(out) :: xc_alias

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 select case(trim(xc_name))
   case('LDA_X+LDA_C_PW')
	 xc_alias='PW'
   case('GGA_X_PBE+GGA_C_PBE')
	 xc_alias='PBE'
   case('LDA_X+LDA_C_PZ')
	 xc_alias='PZ'
   case('LDA_X+LDA_C_WIGNER')
	 xc_alias='W'
   case('LDA_X+LDA_C_HL')
	 xc_alias='HL'
   case('LDA_X+LDA_C_GL')
	 xc_alias='GL'
   case('LDA_X+LDA_C_VWN')
	 xc_alias='VWN'
   case('GGA_X_PBE_R+GGA_C_PBE')
	 xc_alias='revPBE'
   case('GGA_X_RPBE+GGA_C_PBE')
	 xc_alias='RPBE'
   case('GGA_X_PW91+GGA_C_PW91')
	 xc_alias='PW91'
   case('GGA_X_B88+GGA_C_LYP')
	 xc_alias='BLYP'
   case DEFAULT
	 xc_alias=xc_name
 end select

end subroutine get_xc_alias


!!=================================================================
END Module XMLInterface
