!! NAME
!!  ABINITInterface
!!
!! FUNCTION
!!  This module contains routines used to convert PAW atomic dataset
!!  into a suitable formatted file for ABINIT code
!!
!! COPYRIGHT
!! Copyright (C) 2002-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .

Module ABINITInterface

 use atomdata
 use GlobalMath
 use gridmod
 use excor
 use pseudo
 use interpolation_mod
 use pkginfo
 use libxc_mod

 implicit none

 private

 public :: atompaw2abinit

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
 character*10 :: atompaw2abinitver='3.4.1', abinitver='6.1.0+', verdate='june 2013'

!Default name for Abinit file
 integer, parameter :: unit_abinit=22

!Options
 logical :: usexcnhat_def=.false.
 logical :: prtcorewf_def=.false.

!Reduced logarithmic grid (default parameters)
 logical  :: uselog_def=.false.
 integer  :: loggrd_size_def=350
 real(dp) :: loggrd_step_def=0.035_dp

!Real Space Optimization (default parameters)
 logical  :: userso_def=.false.
 real(dp) :: ecut_rso_def=10._dp
 real(dp) :: gfact_rso_def=two
 real(dp) :: werror_rso_def=0.0001_dp

!!=================================================================
!! STRUCTURED DATATYPES
!!=================================================================

!PAW dataset "header" for ABINIT
 type pshead_type
  integer :: basis_size      ! Number of elements for the paw nl basis
  integer :: core_size       ! Number of core states
  integer :: core_meshsz     ! Dimension of radial mesh for core density
                             ! Must be greater or equal than sph_meshsz
  integer :: corewf_meshsz   ! Dimension of radial mesh for core wave functions
  integer :: creatorid       ! ID of psp generator (here creatorID=1 !)
  integer :: hat_meshsz      ! Dimension of radial mesh for hat density
  integer :: l_size          ! Maximum value of l+1 leading to a non zero Gaunt coef
  integer :: lambda          ! Lambda in gaussian type g(r)
  integer :: lmax            ! Maximum value of l for valence states
  integer :: lmax_core       ! Maximum value of l for core states
  integer :: lmn_size        ! Number of elements for the paw basis
  integer :: lmn_size_core   ! Number of lmn elements for the core
  integer :: lmn2_size       ! lmn2_size=lmn_size*(lmn_size+1)/2
  integer :: max_meshsz      ! max. dimension for radial meshes
  integer :: mesh_type       ! Flag: 1 if all meshes are regular grids
                             !       2 if all meshes are log. grids
  integer :: prj_meshsz      ! Dimension of radial mesh for tproj
  integer :: prj_msz_max     ! Maximal dimension of radial mesh for tproj
  integer :: pspcod          ! Psp code number for Abinit (here PAW->pspcod=7 !)
  integer :: pspxc_abinit    ! Abinit code number for the XC functionnal
  integer :: shape_type      ! Shape function type
                             ! shape_type=1 ; g(r)=sinc2
                             ! shape_type=2 ; g(r)=exp
                             ! shape_type=3 ; g(r)=bessel
                             ! shape_type=-1; g(r)=tabulated
  integer :: sph_meshsz      ! Dimension of radial mesh corresponding to PAW radius
  integer :: vale_meshsz     ! Dimension of radial mesh for pseudo valence density
  integer :: vloc_meshsz     ! Dimension of radial mesh for vloc=vhtnzc (with hat in XC)
                             ! For r>r(vloc_meshsz), we assume vloz=-Z/r
  integer :: vlocopt         ! Option for the format of Vloc
                             !  0: Vloc is Vbare (Blochl formulation)
                             !  1: Vloc is VH(tnzc) (Kresse formulation)
                             !  2: Vloc is VH(tnzc) (Kresse formulation without compensation in XC)
  integer :: wav_meshsz      ! Dimension of radial mesh for wave functions (phi, tphi)
  real(dp) :: atomic_charge  ! Total atomic charge
  real(dp) :: core_charge    ! Core charge
  real(dp) :: log_step       ! Step corresponding to exponential term in logarith. radial mesh
  real(dp) :: rad_step       ! Step corresponding to radial mesh
  real(dp) :: rc_hat         ! Radius used in sinus shape function
  real(dp) :: rc_sph         ! Default PAW sphere radius
  real(dp) :: sigma          ! Sigma for gaussian shape function
  character*(fnlen) :: title ! Title for pseudopotential
  integer, pointer :: orbitals(:)      ! orbitals(basis_size) (l quantum number for each basis function)
  integer, pointer :: orbitals_core(:) ! orbitals_core(basis_size) (l quantum number for each core WF)
  real(dp), pointer :: occ(:)          ! occ(basis_size) (occupation for each basis function)
 end type pshead_type

!ABINIT PAW dataset (except header)
 type pawps_type
  real(dp), pointer :: coreden4pr2(:)  ! coreden4pr2(core_meshsz)
                                       ! Gives the core density multiplied by 4Pi.r2
  real(dp), pointer :: tcoreden4pr2(:) ! tcoreden4pr2(core_meshsz)
                                       ! Gives the pseudized core density multiplied by 4Pi.r2
  real(dp), pointer :: tvaleden4pr2(:) ! tvaleden4pr2(core_meshsz)
                                       ! Gives the pseudized core density multiplied by 4Pi.r2 (up to r(vale_meshsz))
  real(dp), pointer :: phi(:,:)        ! phi(sph_meshsz,basis_size)
                                       ! Gives, on the radial grid, the PAW atomic wavefunctions
  real(dp), pointer :: tphi(:,:)       ! tphi(sph_meshsz,basis_size)
                                       ! Gives, on the radial grid, the PAW atomic pseudowavefunctions
  real(dp), pointer :: tproj(:,:)      ! tproj(prj_msz_max,basis_size)
                                       ! Gives, on the radial grid, the PAW projectors for each basis function
  real(dp), pointer :: dij0(:)         ! dij0(lmn2_size)
                                       ! Part of the Dij term (non-local operator) frozen during SCF
  real(dp), pointer :: rhoij0(:)       ! rhoij0(lmn2_size)
                                       ! Atomic initialization of rhoij
  real(dp), pointer :: vbare(:)        ! vbare(sph_meshsz)
                                       ! Gives the "bare" local potential
  real(dp), pointer :: vhtnzc(:)       ! vhtnzc(core_meshsz)
                                       ! Gives the Hartree potential of the pseudo density
                                       ! of the nucleus + core electrons of the atom
 end type pawps_type

!Grid definitions
 type pawrad_type
  integer :: islog                ! 0 if the radial grid is a linear grid, 1 if it is a logarithmic grid
  integer :: meshsz               ! Dimension of max. radial mesh (max. of all mesh sizes)
  real(dp) :: rstep               ! Step corresponding to radial mesh
  real(dp) :: lstep               ! Step corresponding to exponential term (logarithmic mesh)
  real(dp), pointer :: rad(:)     ! rad(max_meshsz)  Coordinates of all the pts of the radial grid
  real(dp), pointer :: radfact(:) ! radfact(max_meshsz) used to compute radial integrals on generalized grid
 end type pawrad_type

!Various useful arrays used during ABINIT PAW dataset generation
 type pawarray_type
  integer, pointer :: indlmn(:,:)      ! indlmn(6,lmn_size) Gives l,m,n,lm,ln,s for i=lmn
  real(dp), pointer :: hatden4pr2(:)   ! hatden4pr2(sph_meshsz) Gives the compensation density
                                       ! multiplied by 4Pi.r2 (following Abinit definition)
  real(dp), pointer :: kij(:)          ! kij(lmn2_size) Kinetic overlap operator
  real(dp), pointer :: shapefunc(:)    ! shapefunc(sph_meshsz) Gives the normalized shape function
                                       ! of the atom used for the compensation density
  real(dp), pointer :: shpnrm(:)       ! shpnrm(l_size) Gives the moment of shapefunction for each l
  real(dp), pointer :: tvaleden4pr2(:) ! tvaleden4pr2(sph_meshsz) Gives the pseudized valence density
                                       ! multiplied by 4Pi.r2 (only part inside spheres)
 end type pawarray_type

!Options for REAL SPACE OPTIMIZATION of projectors
 type pawrso_type
  logical :: userso           ! TRUE if Real Space Optimization is required
  real(dp) :: ecut            ! Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
  real(dp) :: gfact           ! Real Space Optimization parameter: Gamma/Gmax ratio
  real(dp) :: werror          ! Real Space Optimization parameter: max. error W_l allowed
 end type pawrso_type

!Options for transfer of dataset onto a reduced log. grid
 type loggrd_type
  logical :: uselog           ! TRUE if data are transfered on a log. grid before being written
  integer :: meshsz           ! Mesh size for the logarithmic grid
  real(dp) :: log_step        ! Logarithmic step for the logarithmic grid
  real(dp) :: rad_step        ! Radial step for the logarithmic grid
 end type loggrd_type

 CONTAINS

!!=================================================================
!! NAME
!! atompaw2abinit
!!
!! FUNCTION
!! Main routine for building a psp file for PAW calculations
!! in Abinit, starting from the code atomPAW (N. Holzwarth).
!!
!! The present main routine drives the following operations :
!! 1) Read useful data from the output file of AtomPAW
!! 2) Calculate quantities that do not depend from self-consistent PAW
!!    iterations
!! 3) Write an output file to be read as an input PAW atomic data
!!    file in abinit
!!
!!       Units: Energies=Hartree, Lengths=Bohr
!!
!! PARENTS
!! atompaw
!!
!! CHILDREN
!! calc_dij0,calc_rhoij0,calc_shapef,calc_tcore,initmesh,initpawps
!! opt_proj,rdinputabinit,rdpawps1,rdpawps2,wrpawps,wrcorewf
!!=================================================================

 subroutine Atompaw2Abinit(AEOrbit,AEPot,AESCF,PAW,FC,Grid,ifinput)

 type(OrbitInfo), intent(in)     :: AEOrbit
 type(PotentialInfo), intent(in) :: AEPot
 type (SCFInfo),  intent(in)     :: AESCF
 type(Pseudoinfo), intent(in)    ::  PAW
 type(FCInfo), intent(in)        :: FC
 type(GridInfo), intent(in)      :: Grid
 integer, intent(in)             :: ifinput

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 type(pshead_type)   :: pshead
 type(pawps_type)    :: pawps
 type(pawarray_type) :: pawarray
 type(pawrad_type)   :: pawrad
 type(loggrd_type)   :: loggrd
 type(pawrso_type)   :: pawrso
 integer :: id
 logical :: prtcorewf
 character*(fnlen) :: author,file_abinit,xcname
 character*(5000) :: input_string

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 write(std_out,'(2a)') ch10,'========================================================'
 write(std_out,'(3a)')      '==   Atompaw2Abinit - v',trim(atompaw2abinitver),&
&                                                       ':                         =='
 write(std_out,'(a)')       '==   Converts a PAW dataset generated by "Atompaw"    =='
 write(std_out,'(a)')       '==   into a PAW atomic data file readable             =='
 write(std_out,'(3a)')      '==   by ABINIT (v',trim(abinitver),&
&                                                    ')                            =='
 write(std_out,'(a)')       '==                                                    =='
 write(std_out,'(a)')       '==   "ABINIT" can be found at:  http://www.abinit.org =='
 write(std_out,'(a)')       '========================================================'

 call read_inputstring(ifinput,input_string)

!------------------------------------------------------------------
!---- Load PAW dataset header into ABINIT datastructure (sizes of arrays)

 call rdpawps1(AEOrbit,AEPot,PAW,FC,Grid,pshead)

!------------------------------------------------------------------
!---- Read choices from input file

 call rdinputabinit(pshead,pawrso,loggrd,prtcorewf,author,input_string)

!------------------------------------------------------------------
!---- Allocations and initializations

!---- Initialize some useful data
 allocate(pawarray%indlmn(6,pshead%lmn_size))
 call initpawps(pshead,pawarray)

!---- Allocate memory
 allocate(pawps%phi(pshead%wav_meshsz,pshead%basis_size))
 allocate(pawps%tphi(pshead%wav_meshsz,pshead%basis_size))
 allocate(pawps%tproj(pshead%prj_msz_max,pshead%basis_size))
 allocate(pawps%coreden4pr2(pshead%core_meshsz))
 allocate(pawps%tcoreden4pr2(pshead%core_meshsz))
 allocate(pawps%tvaleden4pr2(pshead%vale_meshsz))
 allocate(pawps%dij0(pshead%lmn2_size))
 allocate(pawps%rhoij0(pshead%lmn2_size))
 allocate(pawarray%shapefunc(pshead%wav_meshsz))
 allocate(pawarray%shpnrm(pshead%l_size))
 allocate(pawarray%kij(pshead%lmn2_size))
 allocate(pawarray%tvaleden4pr2(pshead%wav_meshsz))
 allocate(pawarray%hatden4pr2(pshead%sph_meshsz))
 allocate(pawps%vbare(pshead%sph_meshsz))
 allocate(pawps%vhtnzc(pshead%vloc_meshsz))

!---- Initialize radial mesh
 allocate(pawrad%rad(pshead%max_meshsz),pawrad%radfact(pshead%max_meshsz))
 call initmesh(pawrad,pshead%max_meshsz,pshead%mesh_type,pshead%rad_step,pshead%log_step)

!------------------------------------------------------------------
!---- Load PAW dataset "body" into ABINIT datastructure
 call rdpawps2(PAW,FC,pshead,pawps,pawarray,pawrad)

!------------------------------------------------------------------
!---- Compute new data needed by Abinit
 call calc_shapef(pshead,pawarray,pawrad)
 call calc_valden(pshead,pawps,pawarray,pawrad)
 call calc_dij0(pshead,pawps,pawarray,pawrad)
 call calc_rhoij0(pshead,pawps)

!------------------------------------------------------------------
!---- Eventually optimize projectors
 call opt_proj(pshead,pawps,pawrad,pawrso)

!------------------------------------------------------------------
!---- Write PAW dataset file for ABINIT
 xcname=exctype;if (have_libxc) call libxc_getshortname(exctype,xcname)
 file_abinit=TRIM(AEpot%sym)//'.'//TRIM(xcname)//'-paw.abinit'

 call wrpawps(pshead,pawps,pawarray,pawrad,loggrd,file_abinit,&
&             unit_abinit,author)

!------------------------------------------------------------------
!---- Write core Wave Function file for ABINIT
 if (prtcorewf) then
  xcname=exctype;if (have_libxc) call libxc_getshortname(exctype,xcname)
  file_abinit=TRIM(AEpot%sym)//'.'//TRIM(xcname)//'-corewf.abinit'
  call wrcorewf(AEOrbit,FC,pshead,pawrad,loggrd,file_abinit,unit_abinit)
 end if

!------------------------------------------------------------------
!---- Deallocate memory

 deallocate(pshead%orbitals,pshead%orbitals_core,pshead%occ)
 deallocate(pawps%phi,pawps%tphi,pawps%tproj)
 deallocate(pawps%coreden4pr2,pawps%tcoreden4pr2,pawps%tvaleden4pr2)
 deallocate(pawps%dij0,pawps%rhoij0)
 deallocate(pawarray%shapefunc,pawarray%shpnrm,pawarray%kij)
 deallocate(pawarray%tvaleden4pr2,pawarray%hatden4pr2)
 deallocate(pawps%vbare,pawps%vhtnzc)
 deallocate(pawarray%indlmn,pawrad%rad,pawrad%radfact)

 write(std_out,'(2x,a,a)') 'Atompaw2Abinit ended.',ch10

 end subroutine Atompaw2Abinit


!!=================================================================
!! NAME
!! rdpawps1
!!
!! FUNCTION
!! Load the header of a PAW dataset from AtomPAW (only scalars and dimensions)
!!
!! OUTPUT
!!  pshead
!!    %atomic_charge= Total atomic charge
!!    %basis_size= Number of elements for the paw nl basis
!!    %core_size= Number of core states
!!    %core_charge= Core charge
!!    %core_meshsz= Dimension of radial mesh for core density
!!    %corewf_meshsz= Dimension of radial mesh for core wave functions
!!    %flag_vloc_nohat=non-zero if VH(tnZc) "without hat in XC" is present in AtomPAW file
!!    %lambda= Lambda in gaussian type g(r)
!!    %lmn_size= Number of elements for the paw basis
!!    %lmn_size_core= Number of lmn elements for the core
!!    %log_step= Log. step corresponding to radial mesh
!!    %mesh_type= type of radial grid (regular or log)
!!    %occ(basis_size)= occupation for each basis function
!!    %orbitals(basis_size)= Quantum number l for each basis function
!!    %orbitals_core(basis_size)= Quantum number l for each core WF
!!    %pspxc_abinit= ABINIT code number for the exchange-correlation
!!    %rad_step= Step corresponding to radial mesh
!!    %rc_hat= radius for shape function
!!    %rc_sph= Default PAW sphere radius
!!    %shape_type= Shape function type
!!    %sigma= Sigma for gaussian type g(r)
!!    %title= Title for pseudopotential
!!    %vale_meshsz= Dimension of radial mesh for pseudo valence density
!!    %vloc_meshsz= Dimension of radial mesh for vloc=vhtnzc (with hat in XC)
!!    %wav_meshsz= Dimension of radial mesh for phi, tphi ...
!!
!! PARENTS
!! atompaw2abinit
!!=================================================================

 subroutine rdpawps1(AEOrbit,AEPot,PAW,FC,Grid,pshead)

 type(GridInfo),intent(in)       :: Grid
 type(OrbitInfo),intent(in)      :: AEOrbit
 type(PotentialInfo),intent(in)  :: AEPot
 type(Pseudoinfo),intent(in)     :: PAW
 type(FCInfo),intent(in)         :: FC
 type(pshead_type),intent(inout) :: pshead

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ib,ii,ir
 integer :: id(2)=(/0,0/)
 integer,allocatable :: irwf(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!--Read TITLE
 write(unit=pshead%title,fmt='(8a)') &
&   "Paw atomic data for element ",trim(ADJUSTL(AEPot%sym)),&
&   " - Generated with ",trim(atp_package), " v",trim(atp_version)

!--Read BASIS SIZE
 pshead%basis_size=PAW%nbase

!--Read CORE SIZE
 pshead%core_size=0
 do ib=1,AEOrbit%norbit
  if (AEOrbit%iscore(ib)) pshead%core_size=pshead%core_size+1
 end do

!--Read ORBITALS
 pshead%lmn_size=0
 allocate(pshead%orbitals(pshead%basis_size))
 do ib=1,pshead%basis_size
  pshead%orbitals(ib)=PAW%l(ib)
  pshead%lmn_size=pshead%lmn_size+2*pshead%orbitals(ib)+1
 end do

!--Read CORE ORBITALS
 ii=0;pshead%lmn_size_core=0
 allocate(pshead%orbitals_core(pshead%core_size))
 do ib=1,AEOrbit%norbit
  if (AEOrbit%iscore(ib)) then
   ii=ii+1;pshead%orbitals_core(ii)=AEOrbit%l(ib)
   pshead%lmn_size_core=pshead%lmn_size_core+2*AEOrbit%l(ib)+1
  end if
 end do

!--Read OCCUPATIONS
 allocate(pshead%occ(pshead%basis_size))
 do ib=1,pshead%basis_size
  pshead%occ(ib)=PAW%occ(ib)
 end do

!--Read GRID parameters
 pshead%wav_meshsz=PAW%irc+Grid%ishift
 if (usingloggrid(Grid)) then
  pshead%mesh_type=2
  pshead%rad_step=Grid%drdu(1)
  pshead%log_step=Grid%h
 else
  pshead%mesh_type=1
  pshead%rad_step=Grid%h
  pshead%log_step=zero
 end if

!--Read RC
 pshead%rc_sph=PAW%rc

!--Read MESH SIZE FOR CORE DENSITIES
 pshead%core_meshsz=PAW%coretailpoints

!--Read MESH SIZE FOR PSEUDO_VALENCE_DENSITY
 pshead%vale_meshsz=PAW%ivale

!--Read MESH SIZE FOR VLOCION
 pshead%vloc_meshsz=PAW%ivion

!--Read MESH SIZE FOR CORE WF
 ii=0;allocate(irwf(pshead%core_size));irwf(:)=Grid%n
 do ib=1,AEOrbit%norbit
  if (AEOrbit%iscore(ib)) then
   ii=ii+1;ir=Grid%n+1
   do while (ir>1)
    ir=ir-1
    if (abs(AEOrbit%wfn(ir,ib))>tol10) then
     irwf(ii)=min(ir+1,Grid%n);ir=1
    end if
   end do
  end if
 end do
 pshead%corewf_meshsz=maxval(irwf(1:pshead%core_size))
 deallocate(irwf)

!--Read ATOMIC CHARGE
 pshead%atomic_charge=AEPot%nz

!--Read CORE CHARGE
 pshead%core_charge=FC%zcore

!--Read EXCHANGE-CORRELATION TYPE
 if (trim(ADJUSTL(exctype))=="LDA-PW") then
  pshead%pspxc_abinit=7
 elseif (trim(ADJUSTL(exctype))=="GGA-PBE") then
  pshead%pspxc_abinit=11
 else if (have_libxc) then
  call libxc_getid(id)
  if (id(1)>=0.and.id(2)>=0) then
   pshead%pspxc_abinit=-(1000*id(1)+id(2))
  else
   write(std_out,'(/,2x,a)') "Error in Atompaw2Abinit(rdpawps1): unknown XC type !"
   stop
  end if
 else
  write(std_out,'(/,2x,a)') "Error in Atompaw2Abinit(rdpawps1): unknown XC type !"
  stop
 endif

!--Read SHAPE TYPE AND SHAPE PARAM
 pshead%rc_hat=PAW%rc_shap
 if (gaussianshapefunction) then
  pshead%shape_type=1
  pshead%lambda=2
  pshead%sigma=PAW%gausslength
 else if (besselshapefunction) then
  pshead%shape_type=3
  pshead%lambda=0
  pshead%sigma=zero
 else
  pshead%shape_type=2
  pshead%lambda=0
  pshead%sigma=zero
 end if

!Tests of consistency
 if (pshead%core_meshsz<pshead%wav_meshsz) then
  write(std_out,'(/,2x,a)') "Error in Atompaw2Abinit(rdpawps1):"
  write(std_out,'(2x,a)') "   Mesh size for tcore density (CORETAIL_POINTS)"
  write(std_out,'(2x,a)') "   must be greater or equal than MESH_SIZE !"
  stop
 endif
 if (pshead%mesh_type==1) then
  if (pshead%rc_sph>pshead%rad_step*dble(pshead%wav_meshsz-1)+tol8) then
   write(std_out,'(/,2x,a)') "Error in Atompaw2Abinit(rdpawps1):"
   write(std_out,'(2x,a)') "   Radius for PAW spheres (RC)"
   write(std_out,'(2x,a)') "   must be less (or equal) than R(MESH_SIZE) !"
   stop
  endif
 else
  if (pshead%rc_sph>pshead%rad_step*(exp(pshead%log_step*dble(pshead%wav_meshsz-1))-one)+tol8) then
   write(std_out,'(/,2x,a)') "Error in :Atompaw2Abinit(rdpawps1)"
   write(std_out,'(2x,a)') "   Radius for PAW spheres (RC)"
   write(std_out,'(2x,a)') "   must be less (or equal) than R(MESH_SIZE) !"
   stop
  endif
 endif

 end subroutine rdpawps1


!!=================================================================
!! NAME
!! rdinputabinit
!!
!! FUNCTION
!! Read the input file in order to get ABINIT dataset options
!!
!! INPUTS
!!
!! OUTPUT
!!  pshead
!!    %vlocopt= option for Vloc (use of nhat in XC or not)
!!  pawrso
!!    %ecut=Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
!!    %gfact=Real Space Optimization parameter: Gamma/Gmax ratio
!!    %userso=TRUE if REAL Space Optimization is required
!!    %werror=Real Space Optimization parameter: max. error W_l allowed
!!  loggrd
!!    %meshsz=mesh size for the logarithmic grid
!!    %uselog=TRUE if data are transfered on a log. grid before being written
!!    %log_step=logarithmic step for the logarithmic grid
!!  prtcorewf= TRUE is printing of core WF is required
!!
!! SIDE EFFECTS
!!  input_string= string containing a copy of atompaw input file
!!                appended here
!!
!! PARENTS
!! atompaw2abinit
!!
!! CHILDREN
!! extractword,uppercase
!!=================================================================

 subroutine rdinputabinit(pshead,pawrso,loggrd,prtcorewf,author,input_string)

 type(pshead_type),intent(inout) :: pshead
 type(pawrso_type),intent(out)   :: pawrso
 type(loggrd_type),intent(out)   :: loggrd
 logical,intent(out)             :: prtcorewf
 character(len=*),intent(out)    :: author
 character(len=*),intent(inout) :: input_string

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: i_author,i_usexcnhat,i_prtcorewf,i_logspline,i_rsoptim,iend,ok,nn
 character*(fnlen) :: readline,readline_u,inputline,inputword

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 read(std_in,'(a)',advance='no',iostat=ok) readline
 readline_u=readline
 call Uppercase(readline_u)

 i_usexcnhat=index(readline_u,'USEXCNHAT')
 i_prtcorewf=index(readline_u,'PRTCOREWF')
 i_logspline=index(readline_u,'LOGSPLINE')
 i_rsoptim  =index(readline_u,'RSOPTIM')
 i_author   =index(readline_u,'AUTHOR')

!Option for core WF printing
 if (i_prtcorewf>0) then
  prtcorewf=.true.
 else
  prtcorewf=prtcorewf_def
 end if
 if (prtcorewf) then
  write(std_out,'(a,2x,a)') ch10,'Printing of core wave functions file'
 else
  write(std_out,'(a,2x,a)') ch10,'No printing of core wave functions file'
 end if

!Option for use of NHAT in XC
 if (i_usexcnhat>0) then
  pshead%vlocopt=1
 else
  if (     usexcnhat_def) pshead%vlocopt=1
  if (.not.usexcnhat_def) pshead%vlocopt=2
 end if
!To be activated later:
!if (i_vbare>0) pshead%vlocopt=0
 if (pshead%vlocopt==1) &
&  write(std_out,'(a,2x,a)') ch10,'Use of compensation charge in XC terms'
 if (pshead%vlocopt/=1) &
&  write(std_out,'(a,2x,a)') ch10,'No use of compensation charge in XC terms'

!Option for use of REAL SPACE OPTIMIZATION
 pawrso%userso=userso_def
 pawrso%ecut=ecut_rso_def
 pawrso%gfact=gfact_rso_def
 pawrso%werror=werror_rso_def
 if (i_rsoptim>0) then
  pawrso%userso=.true.
  iend=128
  if (i_usexcnhat>i_rsoptim.and.i_usexcnhat-1<iend) iend=i_usexcnhat-1
  if (i_prtcorewf>i_rsoptim.and.i_prtcorewf-1<iend) iend=i_prtcorewf-1
  if (i_logspline>i_rsoptim.and.i_logspline-1<iend) iend=i_logspline-1
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
&   'Real Space optim.: Ecut, Gamma/Gmax, Wl(error) [',&
&   pawrso%ecut,', ',pawrso%gfact,', ',pawrso%werror,']'
 else
  write(std_out,'(a,2x,a)') ch10,'No Real Space Optimization of projectors'
 end if

!Option for spline on a reduced log. grid
 loggrd%uselog=uselog_def
 loggrd%meshsz=loggrd_size_def
 loggrd%log_step=loggrd_step_def
 if (i_logspline>0) then
  loggrd%uselog=.true.
  iend=128
  if (i_usexcnhat>i_logspline.and.i_usexcnhat-1<iend) iend=i_usexcnhat-1
  if (i_prtcorewf>i_logspline.and.i_prtcorewf-1<iend) iend=i_prtcorewf-1
  if (i_rsoptim  >i_logspline.and.i_rsoptim  -1<iend) iend=i_rsoptim  -1
  inputline=""
  if (iend>i_logspline+9) inputline=trim(readline(i_logspline+9:iend))
  if (inputline/="") then
   call extractword(1,inputline,inputword);inputword=trim(inputword)
   if (inputword/="") then
    read(inputword,*) loggrd%meshsz
    call extractword(2,inputline,inputword);inputword=trim(inputword)
    if (inputword/="") read(inputword,*) loggrd%log_step
   end if
  end if
 end if
 if (loggrd%uselog) then
  write(std_out,'(a,2x,a,i3,a,f5.3,a)') ch10,&
&   'Logarithmic grid: Number of pts, logarithmic step [',&
&   loggrd_size_def,', ',loggrd_step_def,']'
 else
  write(std_out,'(a,2x,a)') ch10,'No spline on a reduced log. grid'
 end if

 if(i_author>0) then
   inputline=trim(readline(i_author+6:))
   read(unit=inputline,fmt=*) author
   author=trim(author) ; nn=len(trim(author))
   write(unit=input_string,fmt='(6a)') trim(input_string),char(10),&
&   "ABINITOUT",char(10),readline(1:i_author-1),trim(readline(i_author+nn+10:))
 else
   author="" 
   write(unit=input_string,fmt='(5a)') trim(input_string),char(10),&
&   "ABINITOUT",char(10),trim(readline)
 end if

 end subroutine rdinputabinit


!!=================================================================
!! NAME
!! initpawps
!!
!! FUNCTION
!! Initialize several quantities
!!
!! INPUTS
!!  pshead
!!    %basis_size= Number of elements for the paw nl basis
!!    %core_meshsz= Dimension of radial mesh for core density
!!    %corewf_meshsz= Dimension of radial mesh for core wave functions
!!    %flag_vloc_nohat=non-zero if VH(tnZc) "without hat in XC" is present in AtomPAW file
!!    %lmn_size= Number of elements for the paw basis
!!    %mesh_type= type of radial grid (regular or log)
!!    %orbitals(basis_size)= Quantum number l for each basis function
!!    %rc_hat= radius for shape function
!!    %sph_meshsz= Dimension of radial mesh corresponding to PAW spheres
!!    %vale_meshsz= Dimension of radial mesh for pseudo valence density
!!    %vlocopt= Option for Vloc
!!
!! OUTPUT
!!  pawarray
!!    %indlmn(6,lmn_size)= Gives l,m,n,lm,ln,s for i=lmn
!!  pshead
!!    %creatorid= ID of psp generator (here pspfmt=1 !)
!!    %hat_meshsz= Dimension of radial mesh for shape function
!!    %l_size= Max. value of l+1 leading to a non zero Gaunt coeffs
!!    %lmax= Maximum value of l for valence states
!!    %lmax_core= Maximum value of l for core states
!!    %lmn2_size= lmn_size*(lmn_size+1)/2
!!    %max_meshsz=  Max. dimension for radial meshes
!!    %prj_msz_max= Maximal dimension of radial mesh for tproj
!!    %pspcod= Psp code number for Abinit (here PAW->pspcod=7 !)
!!    %rc_hat= radius for shape function
!!    %rc_sph= Default PAW sphere radius
!!    %vloc_meshsz= Dimension of radial mesh for vloc=htnzc
!!
!! PARENTS
!! atompaw2abinit
!!=================================================================

 subroutine initpawps(pshead,pawarray)

 type(pshead_type),intent(inout)   :: pshead
 type(pawarray_type),intent(inout) :: pawarray

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------
 real(dp),parameter :: rmax_vloc=10._dp !We decide to cut at r=10 u.a because of numeric noise...
 integer :: ii1,ii2,ii3,ib,il,ilm,ilmn,iln,nproj

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!Format codes for Abinit pseudopotential
 read(unit=atompaw2abinitver(1:1),fmt='(i1)') ii1
 read(unit=atompaw2abinitver(3:3),fmt='(i1)') ii2
 read(unit=atompaw2abinitver(5:5),fmt='(i1)') ii3
 pshead%creatorid=1000+100*ii1+10*ii2+ii3
 pshead%pspcod=7

!Maximum value for l quantum number
 pshead%lmax=maxval(pshead%orbitals(1:pshead%basis_size))
 pshead%lmax_core=maxval(pshead%orbitals_core(1:pshead%core_size))
 pshead%l_size=2*pshead%lmax+1
 pshead%lmn2_size=pshead%lmn_size*(pshead%lmn_size+1)/2

!Radial mesh for PAW spheres
 if (pshead%mesh_type==1) then
  pshead%sph_meshsz=min(1+nint(pshead%rc_sph/pshead%rad_step),pshead%wav_meshsz)
  pshead%rc_sph=pshead%rad_step*dble(pshead%sph_meshsz-1)
 else
  pshead%sph_meshsz=min(1+nint(log(one+pshead%rc_sph/pshead%rad_step)/pshead%log_step),pshead%wav_meshsz)
  pshead%rc_sph=pshead%rad_step*(exp(pshead%log_step*dble(pshead%sph_meshsz-1))-one)
 endif

!Radial mesh for hat density
 if (pshead%rc_hat<0) pshead%rc_hat=pshead%rc_sph
 if (abs(pshead%rc_hat-pshead%rc_sph)<tol8) then
  pshead%hat_meshsz=pshead%sph_meshsz
  pshead%rc_hat=pshead%rc_sph
 else
  if (pshead%mesh_type==1) then
   pshead%hat_meshsz=int(one+pshead%rc_hat/pshead%rad_step)
  else
   pshead%hat_meshsz=int(log(one+pshead%rc_hat/pshead%rad_step)/pshead%log_step)+1
  endif
 endif

!Radial mesh for projectors
 pshead%prj_msz_max=four*pshead%sph_meshsz
 pshead%prj_meshsz=pshead%sph_meshsz

!Radial mesh for local potential
 if (pshead%vloc_meshsz<=0) then
  if (pshead%mesh_type==1) then
   ii1=int(one+rmax_vloc/pshead%rad_step)
  else
   ii1=int(log(one+rmax_vloc/pshead%rad_step)/pshead%log_step)+1
  endif
  pshead%vloc_meshsz=max(pshead%sph_meshsz,pshead%core_meshsz,ii1)
 endif

!Radial mesh for core wave functions
 pshead%corewf_meshsz=max(pshead%corewf_meshsz,pshead%sph_meshsz)

!Total radial mesh
 pshead%max_meshsz=max(pshead%sph_meshsz,pshead%hat_meshsz,pshead%core_meshsz,pshead%vale_meshsz,&
&                      pshead%vloc_meshsz,pshead%prj_msz_max,pshead%prj_meshsz,pshead%corewf_meshsz)

!Initialization of the orbital basis indexes (indlmn)
 ilmn=0;iln=0
 do il=0,pshead%lmax
  nproj=0
  do ib=1,pshead%basis_size
   if (pshead%orbitals(ib)==il) then
    nproj=nproj+1;iln=iln+1
    do ilm=1,2*il+1
     pawarray%indlmn(1,ilmn+ilm)=il
     pawarray%indlmn(2,ilmn+ilm)=ilm-(il+1)
     pawarray%indlmn(3,ilmn+ilm)=nproj
     pawarray%indlmn(4,ilmn+ilm)=il*il+ilm
     pawarray%indlmn(5,ilmn+ilm)=iln
     pawarray%indlmn(6,ilmn+ilm)=1
    enddo
    ilmn=ilmn+2*il+1
   endif
  enddo
 enddo

 end subroutine initpawps


!!=================================================================
!! NAME
!! initmesh
!!
!! FUNCTION
!! Initialize the radial mesh (here only regular mesh)
!!
!! INPUTS
!!  meshsz= Size of radial mesh
!!  meshtype= type of radial mesh (lin. or log.)
!!  radstep= Radial step of radial mesh
!!  logstep= Log. step of radial mesh
!!
!! OUTPUT
!!  pawrad
!!   %meshsz= Size of radial grid
!!   %rad(meshsz)= Coordinates of the radial grid
!!   %radfact(meshsz)= Factors used in Simpson integrals
!!   %lstep= Step corresponding to exponential in log. mesh
!!   %rstep= Step corresponding to radial mesh
!!
!! PARENTS
!! atompaw2abinit
!!=================================================================

 subroutine initmesh(pawrad,meshsz,meshtype,radstep,logstep)

 integer,intent(in)              :: meshsz,meshtype
 real(dp),intent(in)             :: radstep,logstep
 type(pawrad_type),intent(inout) :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ir

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 pawrad%meshsz=meshsz
 pawrad%rstep=radstep

 if (meshtype==1) then
  pawrad%islog=0
  pawrad%lstep=zero
  do ir=1,meshsz
   pawrad%rad(ir)=(ir-1)*radstep
   pawrad%radfact(ir)=one
  enddo
 else
  pawrad%islog=1
  pawrad%lstep=logstep
  do ir=1,meshsz
   pawrad%rad(ir)=radstep*(exp((ir-1)*logstep)-1.d0)
   pawrad%radfact(ir)=pawrad%rad(ir)+radstep
  enddo
 end if

 end subroutine initmesh


!!=================================================================
!! NAME
!! rdpawps2
!!
!! FUNCTION
!! Load the body of a PAW dataset from AtomPAW
!! (the header must have been loaded before)
!!
!! INPUTS
!!  pshead
!!    %basis_size= Number of elements for the paw nl basis
!!    %core_meshsz= Dimension of radial mesh for core density
!!    %lmn_size= Number of elements for the paw basis
!!    %orbitals(basis_size)= Quantum number l for each basis function
!!    %shape_type= Shape function type
!!    %sph_meshsz= Dimension of radial mesh corresponding to PAW spheres
!!    %vale_meshsz= Dimension of radial mesh for pseudo valence
!!    %vloc_meshsz= Dimension of radial mesh for vloc=vhtnzc
!!    %vlocopt= Option for Vloc
!!    %wav_meshsz= Dimension of radial mesh for phi, tphi ...
!!  pawarray
!!    %indlmn(6,lmn_size)= Gives l,m,n,lm,ln,s for i=lmn
!!  pawrad= radial grid definitions
!!
!! OUTPUT
!!  pawarray
!!    %kij(lmn2_size)= Kinetic overlap operator
!!  pawps
!!    %coreden4pr2(core_meshsz)= Core density multiplied by 4Pi.r2
!!    %tcoreden4pr2(core_meshsz)= Pseudized core density multiplied by 4Pi.r2
!!    %tvaleden4pr2(vale_meshsz)= Pseudized valence density multiplied by 4Pi.r2
!!    %phi(wav_meshsz,basis_size)= PAW atomic wavefunctions on the radial grid
!!    %tphi(wav_meshsz,basis_size)= PAW atomic pseudo-wavefunctions on the radial grid
!!    %tproj(prj_msz_max,basis_size)= PAW projectors on the radial grid
!!    %vbare(sph_meshsz)= bare local potential (part of VH(tnzc))
!!
!! SIDE EFFECTS
!!  pshead
!!    %lambda= Lambda in gaussian type g(r)
!!    %sigma= Sigma for gaussian type g(r)
!!
!! PARENTS
!! atompaw2abinit
!!=================================================================

 subroutine rdpawps2(PAW,FC,pshead,pawps,pawarray,pawrad)

 type(Pseudoinfo),intent(in)       :: PAW
 type(FCInfo),intent(in)           :: FC
 type(pshead_type),intent(inout)   :: pshead
 type(pawps_type),intent(inout)    :: pawps
 type(pawarray_type),intent(inout) :: pawarray
 type(pawrad_type),intent(in)      :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ib,ilm,ilmn,iln,jlm,jlmn,j0lmn,jln,klmn

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!--Read WAVE FUNCTIONS PHI, PSEUDO WAVE FUNCTIONS TPHI and PROJECTORS TPROJ
 do ib=1,pshead%basis_size
  pawps%phi  (1:pshead%wav_meshsz,ib)=PAW%ophi (1:pshead%wav_meshsz,ib)
  pawps%tphi (1:pshead%wav_meshsz,ib)=PAW%otphi(1:pshead%wav_meshsz,ib)
  pawps%tproj(1:pshead%wav_meshsz,ib)=PAW%otp  (1:pshead%wav_meshsz,ib)
 end do

!--Read core density CORE_DENSTY
 pawps%coreden4pr2(1:pshead%wav_meshsz)=FC%coreden(1:pshead%wav_meshsz)

!--Read pseudized core density CORETAIL_DENSITY
 pawps%tcoreden4pr2(1:pshead%core_meshsz)=PAW%tcore(1:pshead%core_meshsz)

!--Read pseudized valence density PSEUDO_VALENCE_DENSITY
 pawps%tvaleden4pr2(1:pshead%vale_meshsz)=PAW%tden(1:pshead%vale_meshsz)

!--Read Vbare potential VLOCFUN     --no longer needed
 pawps%vbare(1:pshead%sph_meshsz)=half*PAW%vloc(1:pshead%sph_meshsz) ! Conversion Ry->Ha

!--Read local ionic potential VLOCION -- ionic local potential with correct
!  unscreening of vxc   (vxc(tilde_core + tilde_n + hat_hat_n) or vxc(tilde_core + tilde_n))
 if (pshead%vlocopt==1) then
  pawps%vhtnzc(1:pshead%vloc_meshsz)=half*PAW%abinitvloc(1:pshead%vloc_meshsz)
 else
  pawps%vhtnzc(1:pshead%vloc_meshsz)=half*PAW%abinitnohat(1:pshead%vloc_meshsz)
 end if

!--Read Kinetic term KINETIC_ENERGY_MATRIX
 pawarray%kij=zero
 do jlmn=1,pshead%lmn_size
  j0lmn=jlmn*(jlmn-1)/2
  jlm=pawarray%indlmn(4,jlmn);jln=pawarray%indlmn(5,jlmn)
  do ilmn=1,jlmn
   klmn=j0lmn+ilmn
   ilm=pawarray%indlmn(4,ilmn);iln=pawarray%indlmn(5,ilmn)
   if (ilm==jlm) pawarray%kij(klmn)=half*PAW%Kij(iln,jln) ! Conversion Ry->Ha
  enddo
 enddo

!--Copy pseudized core density to core density outside spheres
 if (pshead%core_meshsz>pshead%wav_meshsz) &
&  pawps%coreden4pr2(pshead%wav_meshsz+1:pshead%core_meshsz) = &
&  pawps%tcoreden4pr2(pshead%wav_meshsz+1:pshead%core_meshsz)

!--Test pseudo valence density
 if (pshead%vale_meshsz>0) then
  write(std_out,'(/,2x,a,a,f5.2,a,a,g11.4)') 'Atompaw2Abinit info:',&
&   '  At r=',pawrad%rad(pshead%vale_meshsz),' a.u.,',&
&   ' Pseudo_valence_density= ', pawps%tvaleden4pr2(pshead%vale_meshsz) &
&                               /pawrad%rad(pshead%vale_meshsz)**2/(4*pi)
  write(std_out,'(2x,a)') '  This quantity must be as small as possible.'
 endif

!--Test potential
 if (pshead%vlocopt==1.or.pshead%vlocopt==2) then
  write(std_out,'(/,2x,a,a,f5.2,a,a,g11.4)') 'Atompaw2Abinit info:',&
&   '  At r_vloc=',pawrad%rad(pshead%vloc_meshsz),' a.u.,',&
&   ' VHartree(ntild(Zv+Zc))= -Zv/r + ', pawps%vhtnzc(pshead%vloc_meshsz) &
&    +(pshead%atomic_charge-pshead%core_charge)/pawrad%rad(pshead%vloc_meshsz)
  write(std_out,'(2x,a)') '  This quantity must be as small as possible.'
 endif

 end subroutine rdpawps2


!!=================================================================
!! NAME
!! calc_shapef
!!
!! FUNCTION
!! Compute normalized shape function and its moments (shpnrm)
!!
!! INPUTS
!!  pawrad= radial grid definitions
!!  pshead
!!    %hat_meshsz= Dimension of radial mesh for hat density
!!    %lambda= Lambda in gaussian type g(r)
!!    %l_size= Max. value of l+1 leading to a non zero Gaunt coeffs
!!    %rad_step= Step corresponding to radial mesh
!!    %rc_hat= radius for hat density
!!    %shape_type= Shape function type
!!    %sigma= Sigma for gaussian type g(r)
!!
!! OUTPUT
!!  pawarray
!!    %shapefunc(wav_meshsz)= Normalized shape function
!!    %shpnrm(l_size)= Moments of shape function for each l
!!
!! PARENTS
!! atompaw2abinit
!!
!! CHILDREN
!! csimp,jbessel,shapebes
!!=================================================================

 subroutine calc_shapef(pshead,pawarray,pawrad)

 type(pshead_type),intent(in)    :: pshead
 type(pawarray_type),intent(out) :: pawarray
 type(pawrad_type),intent(in)    :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: il,ir
 real(dp) :: arg,dum1,dum2,jbes1,jbes2,qr,shapefunc1,shapefunc2
 real(dp) :: al(2),ql(2)
 real(dp), allocatable :: ff(:)
!Statement functions
 shapefunc1(arg)= exp(-(arg/pshead%sigma)**pshead%lambda)
 shapefunc2(arg)= (sin(pi*arg/pshead%rc_hat)/(pi*arg/pshead%rc_hat))**2

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!Select the suitable shape function
 pawarray%shapefunc(:)=zero
 pawarray%shapefunc(1)=one
 if (pshead%shape_type==1) then
  do ir=2,pshead%hat_meshsz
   pawarray%shapefunc(ir)=shapefunc1(pawrad%rad(ir))
  enddo
 else if (pshead%shape_type==2) then
  do ir=2,pshead%hat_meshsz
   pawarray%shapefunc(ir)=shapefunc2(pawrad%rad(ir))
  enddo
 else if (pshead%shape_type==3) then
  call shapebes(al,ql,0,pshead%rc_hat)
  do ir=1,pshead%hat_meshsz
   qr=ql(1)*pawrad%rad(ir);CALL jbessel(jbes1,dum1,dum2,0,0,qr)
   qr=ql(2)*pawrad%rad(ir);CALL jbessel(jbes2,dum1,dum2,0,0,qr)
   pawarray%shapefunc(ir)=al(1)*jbes1+al(2)*jbes2
  enddo
 else if (pshead%shape_type/=-1) then
  write(std_out,'(/,2x,a)') 'Error in Atompaw2Abinit(calc_shapef): invalid value of pshead%shape_type !'
  stop
 endif

!Compute moments of shape function
 allocate(ff(pshead%hat_meshsz));pawarray%shpnrm=zero;ff(1)=zero
 do il=1,pshead%l_size
  ff(2:pshead%hat_meshsz)=pawarray%shapefunc(2:pshead%hat_meshsz)*pawrad%rad(2:pshead%hat_meshsz)**(2*il)
  call csimp(ff,pawrad,pshead%hat_meshsz,pawarray%shpnrm(il))
 enddo
 deallocate(ff)

!Normalize shape function
 pawarray%shapefunc(1:pshead%hat_meshsz)=pawarray%shapefunc(1:pshead%hat_meshsz)/pawarray%shpnrm(1)

end subroutine calc_shapef


!!=================================================================
!! NAME
!! calc_valden
!!
!! FUNCTION
!! Compute pseudized valence density (n_tild)
!! and compensation density (n_hat, following ABINIT def.)
!!
!! INPUTS
!!  pawarray
!!    %shapefunc(wav_meshsz)= Normalized shape function
!!  pawps
!!    %coreden4pr2(core_meshsz)= Core density multiplied by 4Pi.r2
!!    %phi(wav_meshsz,basis_size)= PAW atomic wavefunctions on the radial grid
!!    %tcoreden4pr2(core_meshsz)= Pseudized core density multiplied by 4Pi.r2
!!    %tphi(wav_meshsz,basis_size)= PAW atomic pseudo-wavefunctions on the radial grid
!!    %tvaleden4pr2(core_meshsz)= Pseudized valence density multiplied by 4Pi.r2 (if read in input file)
!!  pawrad= radial grid definitions
!!  pshead
!!    %basis_size= Number of elements for the paw nl basis
!!    %occ(basis_size)= occupation for each basis function
!!    %sph_meshsz= Dimension of radial mesh corresponding to PAW spheres
!!    %vale_meshsz= Dimension of radial mesh for pseudized valence density (0 if not in input file)
!!    %wav_meshsz= Dimension of radial mesh for wave functions
!!
!! OUTPUT
!!  pawarray
!!    %hatden4pr2(sph_meshsz)= Compensation density *4pi.r2 (following Abinit s definition)
!!    %tvaleden4pr2(wav_meshsz)= Pseudized valence density *4pi.r2 (only part inside spheres)
!!
!! PARENTS
!! atompaw2abinit
!!
!! CHILDREN
!! csimp
!!=================================================================

 subroutine calc_valden(pshead,pawps,pawarray,pawrad)

 type(pshead_type),intent(in)      :: pshead
 type(pawps_type),intent(in)       :: pawps
 type(pawarray_type),intent(inout) :: pawarray
 type(pawrad_type),intent(in)      :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ib,meshsz
 real(dp) :: qq
 real(dp),allocatable :: den4pr2(:),valeden4pr2(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

!Compute pseudized valence density (only inside PAW augmentation region)
 meshsz=pshead%wav_meshsz
 pawarray%tvaleden4pr2=zero
 if (pshead%vale_meshsz>0) then
  pawarray%tvaleden4pr2(1:meshsz)=pawps%tvaleden4pr2(1:meshsz)
 else
  do ib=1,pshead%basis_size
   pawarray%tvaleden4pr2(1:meshsz)=pawarray%tvaleden4pr2(1:meshsz) &
&                                 +pshead%occ(ib)*pawps%tphi(1:meshsz,ib)**2
  enddo
 end if

!Compute valence density (only inside augmentation region) - temporary needed here
 meshsz=pshead%sph_meshsz
 allocate(valeden4pr2(meshsz));valeden4pr2=zero
 do ib=1,pshead%basis_size
  valeden4pr2(1:meshsz)=valeden4pr2(1:meshsz) &
&                      +pshead%occ(ib)*pawps%phi(1:meshsz,ib)**2
 enddo

!Compute compensation density
 allocate(den4pr2(meshsz));pawarray%hatden4pr2=zero
 den4pr2(1:meshsz)=valeden4pr2(1:meshsz)-pawarray%tvaleden4pr2(1:meshsz)
 call csimp(den4pr2,pawrad,meshsz,qq)
 pawarray%hatden4pr2(1:meshsz)=qq*pawrad%rad(1:meshsz)**2*pawarray%shapefunc(1:meshsz)

 deallocate(den4pr2,valeden4pr2)

 end subroutine calc_valden


!!=================================================================
!! NAME
!! calc_dij0
!!
!! FUNCTION
!! Compute "frozen" values of Dij = Dij0
!!
!! INPUTS
!!  pawarray
!!    %indlmn(6,lmn_size)= Gives l,m,n,lm,ln,s for i=lmn
!!    %kij(lmn2_size)= Kinetic overlap operator
!!    %shapefunc(wav_meshsz)= Normalized shape function
!!  pawps
!!    %coreden4pr2(core_meshsz)= Core density multiplied by 4Pi.r2
!!    %phi(wav_meshsz,basis_size)= PAW atomic wavefunctions on the radial grid
!!    %tphi(wav_meshsz,basis_size)= PAW atomic pseudo-wavefunctions on the radial grid
!!    %vhtnzc(core_meshsz)= Hartree potential of the ps-density
!!                          of the nucleus + core electrons
!!  pawrad= radial grid definitions
!!  pshead
!!    %atomic_charge= Total atomic charge
!!    %core_meshsz= Dimension of radial mesh for core density
!!    %lmn_size= Number of elements for the paw basis
!!    %lmn2_size=pshead%lmn_size*(pshead%lmn_size+1)/2
!!    %sph_meshsz= Dimension of radial mesh corresponding to PAW spheres
!!    %vlocopt= Option for Vloc
!!
!! OUTPUT
!!  pawps
!!    %dij0(lmn2_size)= Frozen part of the Dij term
!!
!! PARENTS
!! atompaw2abinit
!!
!! CHILDREN
!! csimp,grid2pawrad
!!=================================================================

 subroutine calc_dij0(pshead,pawps,pawarray,pawrad)

 type(pshead_type),intent(in)   :: pshead
 type(pawps_type),intent(inout) :: pawps
 type(pawarray_type),intent(in) :: pawarray
 type(pawrad_type),intent(inout)   :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: il,ilm,iln,ilmn,j0lmn,jl,jlm,jln,jlmn,klmn,meshsz
 real(dp) :: intg,intvh,ecoul,qq
 real(dp),allocatable :: ff(:),vhnzc(:)
 type(GridInfo) :: Grid_core

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 meshsz=pshead%sph_meshsz;allocate(ff(meshsz))
 call NullifyGrid(Grid_core)
 call grid2pawrad(Grid_core,pawrad,pshead%core_meshsz,-1)

!Kinetic part of Dij0
 pawps%dij0(1:pshead%lmn2_size)=pawarray%kij(1:pshead%lmn2_size)

!Computation of <phi_i|Vh(nZc)|phi_j>
 allocate(vhnzc(pshead%core_meshsz))
 call poisson_marc(Grid_core,qq,pawps%coreden4pr2,vhnzc,ecoul)
 vhnzc=half*vhnzc  ! Ryd. -> Ha
 vhnzc(2:pshead%core_meshsz)=(vhnzc(2:pshead%core_meshsz)-pshead%atomic_charge)/pawrad%rad(2:pshead%core_meshsz)
 call extrapolate(Grid_core,vhnzc)
 
 do jlmn=1,pshead%lmn_size
  j0lmn=jlmn*(jlmn-1)/2
  jlm=pawarray%indlmn(4,jlmn);jln=pawarray%indlmn(5,jlmn)
  do ilmn=1,jlmn
   klmn=j0lmn+ilmn
   ilm=pawarray%indlmn(4,ilmn);iln=pawarray%indlmn(5,ilmn)
   if (jlm==ilm) then
    ff(1:meshsz)=pawps%phi(1:meshsz,iln)*pawps%phi(1:meshsz,jln)*vhnzc(1:meshsz)
    call csimp(ff,pawrad,meshsz,intg)

    pawps%dij0(klmn)=pawps%dij0(klmn)+intg
   endif
  enddo
 enddo
 deallocate(vhnzc)

!Computation of <tphi_i|Vh(tnZc)|tphi_j>
 if (pshead%vlocopt==1.or.pshead%vlocopt==2) then
  do jlmn=1,pshead%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jlm=pawarray%indlmn(4,jlmn);jln=pawarray%indlmn(5,jlmn)
   do ilmn=1,jlmn
    klmn=j0lmn+ilmn
    ilm=pawarray%indlmn(4,ilmn);iln=pawarray%indlmn(5,ilmn)
    if (jlm==ilm) then
     ff(1:meshsz)=pawps%tphi(1:meshsz,iln)*pawps%tphi(1:meshsz,jln)*pawps%vhtnzc(1:meshsz)
     call csimp(ff,pawrad,meshsz,intg)
     pawps%dij0(klmn)=pawps%dij0(klmn)-intg
    endif
   enddo
  enddo
 endif

!Computation of int[Vh(tnzc)*Qijhat(r)dr]
 if (pshead%vlocopt==1.or.pshead%vlocopt==2) then
  ff(1:meshsz)=pawps%vhtnzc(1:meshsz)*pawarray%shapefunc(1:meshsz)*pawrad%rad(1:meshsz)**2
  call csimp(ff,pawrad,meshsz,intvh)
  do jlmn=1,pshead%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jl=pawarray%indlmn(1,jlmn);jln=pawarray%indlmn(5,jlmn);jlm=pawarray%indlmn(4,jlmn)
   do ilmn=1,jlmn
    klmn=j0lmn+ilmn
    il=pawarray%indlmn(1,ilmn);iln=pawarray%indlmn(5,ilmn);ilm=pawarray%indlmn(4,ilmn)
    if (ilm==jlm) then
     ff(1:meshsz)=(pawps%phi (1:meshsz,iln)*pawps%phi (1:meshsz,jln)&
&                 -pawps%tphi(1:meshsz,iln)*pawps%tphi(1:meshsz,jln))
     call csimp(ff,pawrad,meshsz,intg)
     pawps%dij0(klmn)=pawps%dij0(klmn)-intvh*intg
    endif
   enddo
  enddo
 endif

 deallocate(ff)
 call DestroyGrid(Grid_core)

 end subroutine calc_dij0


!!=================================================================
!! NAME
!! calc_rhoij0
!!
!! FUNCTION
!! Compute initial values (atomic initialization) of rhoij0
!!
!! INPUTS
!!  pshead
!!    %basis_size= Number of elements for the paw nl basis
!!    %occ(basis_size)= occupation for each basis function
!!    %orbitals(basis_size)= Quantum number l for each basis function
!!
!! OUTPUT
!!  pawps
!!    %rhoij0= Atomic initialization of rhoij
!!
!! PARENTS
!! atompaw2abinit
!!=================================================================

 subroutine calc_rhoij0(pshead,pawps)

 type(pshead_type),intent(in)   :: pshead
 type(pawps_type),intent(inout) :: pawps

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ib,ilmn,ilmn0,ll2

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 pawps%rhoij0=zero

!Store occupations into Rhoij0
!(This is only true if bound states are included in the wave function basis
 ilmn0=0
 do ib=1,pshead%basis_size
  ll2=2*pshead%orbitals(ib)+1
  do ilmn=ilmn0+1,ilmn0+ll2
   pawps%rhoij0(ilmn*(ilmn+1)/2)=pshead%occ(ib)/dble(ll2)
  enddo
  ilmn0=ilmn0+ll2
 enddo

 end subroutine calc_rhoij0


!!=================================================================
!! NAME
!! calc_vloc
!!
!! FUNCTION
!! Compute local potential in the Bloechl and Kresse formulation
!!
!! INPUTS
!!  pshead
!!    %basis_size= Number of elements for the paw nl basis
!!    %occ(basis_size)= occupation for each basis function
!!    %orbitals(basis_size)= Quantum number l for each basis function
!!
!! OUTPUT
!!  pawps
!!    %vhtnzc= Atomic initialization of rhoij
!!
!! PARENTS
!! atompaw2abinit
!!=================================================================

 subroutine calc_vloc(FC,nz,PAW,pawps,pawrad,pshead)

 real(dp) ,intent(in):: nz
 type(Pseudoinfo),intent(in)       :: PAW
 type(FCInfo),intent(in)           :: FC
 type(pshead_type),intent(in)   :: pshead
 type(pawps_type),intent(inout) :: pawps
 type(pawrad_type),intent(inout)   :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii,irc
 real(dp):: etxc,eexc,qeff,q00,rat
 type(GridInfo) :: Grid,Grid1
 real(dp), allocatable :: dd(:),vv(:),vxc1(:),vxc2(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 call NullifyGrid(Grid)
 call grid2pawrad(Grid,pawrad,pshead%core_meshsz,-1)

 allocate(dd(pshead%core_meshsz))
 allocate(vv(pshead%core_meshsz))
 dd(1:pshead%core_meshsz)=FC%coreden(1:pshead%core_meshsz)-PAW%tcore(1:pshead%core_meshsz)
 qeff=-nz+integrator(Grid,dd,1,PAW%irc_core)
 dd(1:pshead%core_meshsz)=PAW%tcore(1:pshead%core_meshsz)+qeff*PAW%hatden(1:pshead%core_meshsz)
 call poisson_marc(Grid,q00,dd,vv,rat)
 pawps%vhtnzc(2:pshead%core_meshsz)=half*(PAW%vloc(2:pshead%core_meshsz)+vv(2:pshead%core_meshsz)/pawrad%rad(2:pshead%core_meshsz))
 call extrapolate(Grid,pawps%vhtnzc)

 call DestroyGrid(Grid)
 if(pshead%vloc_meshsz>pshead%core_meshsz)then
   pawps%vhtnzc(pshead%core_meshsz+1:pshead%vloc_meshsz)=half*(PAW%vloc(pshead%core_meshsz+1:pshead%vloc_meshsz)+&
&       vv(pshead%core_meshsz)/pawrad%rad(pshead%core_meshsz+1:pshead%vloc_meshsz))
 end if
 deallocate(dd,vv)

 if (pshead%vlocopt==1) then
   allocate(dd(pshead%vloc_meshsz))
   allocate(vxc1(pshead%vloc_meshsz))
   allocate(vxc2(pshead%vloc_meshsz))
   call grid2pawrad(Grid,pawrad,pshead%vloc_meshsz,-1)
   dd(1:pshead%vloc_meshsz)=PAW%tden(1:pshead%vloc_meshsz)+PAW%tcore(1:pshead%vloc_meshsz)
   call exch(Grid,dd,vxc1,etxc,eexc)   !store vxc(\tilde(n)+\tilde(n_c))
   dd(1:pshead%vloc_meshsz)=PAW%den(1:pshead%vloc_meshsz)-PAW%tden(1:pshead%vloc_meshsz)
   irc=max(PAW%irc,PAW%irc_shap,PAW%irc_vloc,PAW%irc_core)
   qeff=integrator(Grid,dd,1,irc)
   dd(1:pshead%vloc_meshsz)=PAW%tden(1:pshead%vloc_meshsz)+PAW%tcore(1:pshead%vloc_meshsz)+qeff*PAW%hatden(1:pshead%vloc_meshsz)
   call exch(Grid,dd,vxc2,etxc,eexc) 
   pawps%vhtnzc(2:pshead%vloc_meshsz)=pawps%vhtnzc(2:pshead%vloc_meshsz)+&
&                        half*(vxc1(2:pshead%vloc_meshsz)-vxc2(2:pshead%vloc_meshsz))/pawrad%rad(2:pshead%vloc_meshsz)
   call extrapolate(Grid,pawps%vhtnzc)
   deallocate(dd,vxc1,vxc2)
   call DestroyGrid(Grid)
 end if

 end subroutine calc_vloc

 
!!=================================================================
!! NAME
!! opt_proj
!!
!! FUNCTION
!! Apply Real Space Optimization (RSO) on non-local projectors in order
!! to smooth them and cut large reciprocal vectors contribution.
!! Directly written from:
!!  RD King-smith, MC Payne and JS Lin, Phys. Rev. B, 44 (1991), 13063
!!
!! INPUTS
!!  pawrad=datastructure containing mesh size definitions
!!  pawrso
!!    %ecut_rso=Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
!!    %gfact_rso=Real Space Optimization parameter: Gamma/Gmax ratio
!!    %userso=TRUE if REAL Space Optimization is required
!!    %werror_rso=Real Space Optimization parameter: max. error W_l allowed
!!  pshead
!!    %basis_size= Number of elements for the paw nl basis
!!    %orbitals(basis_size)= Quantum number l for each basis function
!!    %prj_meshsz=Dimension of radial mesh for tproj (initialized at input)
!!
!! SIDE EFFECTS
!!  pshead
!!    %prj_meshsz= Dimension of radial mesh for tproj
!!  pawps
!!    %tproj(prj_msz_max,basis_size)= projectors on partial waves
!!
!! PARENTS
!! atompaw2abinit
!!
!! CHILDREN
!! aamat,csimp,gauleg,jbessel,DGETRF,DGETRS
!!=================================================================

 subroutine opt_proj(pshead,pawps,pawrad,pawrso)

 type(pshead_type),intent(inout) :: pshead
 type(pawps_type),intent(inout)  :: pawps
 type(pawrad_type),intent(in)    :: pawrad
 type(pawrso_type),intent(in)    :: pawrso

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer,parameter :: ngaussmax=50
 integer, parameter :: nqmax=500

 integer :: ib,info,iq,iqp,ir,ll,nqgauss1,nqgauss2,r0_meshsz,rp_meshsz
 real(dp) :: amat,bess,bessp,dum,dq,gamma,gmax,r0,rc_prj,wwmax,wwl,xx

 integer, allocatable :: iwork(:)
 real(dp),allocatable :: am(:,:),bm(:),chi1g(:,:),chi2g(:),chireg(:,:),ff(:),&
&                        gg(:),qgauss1(:),qgauss2(:),wgauss1(:),wgauss2(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (.not.pawrso%userso) return

 write(std_out,'(3(/,2x,a),/,2x,a,f7.2,/,2x,a,f6.2,/,2x,a,g11.3)') &
&     'Atompaw2Abinit info:',&
&     '  Optimizing non-local projectors',&
&     '  using Real Space Optimization...',&
&     '  Parameters: Ecut (Hartree)=',pawrso%ecut,&
&     '              Gamma/Gmax    =',pawrso%gfact,&
&     '              Wl max (error)=',pawrso%werror

!Initialize data for optimization
 if (pawrad%islog==1) then
  rc_prj=pawrad%rstep*(exp(dble(pshead%prj_meshsz-1)*pawrad%lstep)-1.d0)
 else
  rc_prj=pawrad%rstep*dble(pshead%prj_meshsz-1)
 end if
 r0=rc_prj/1.035_dp
 gmax=sqrt(two*pawrso%ecut)
 gamma=pawrso%gfact*gmax
 wwmax=10000000000._dp

 rp_meshsz=pshead%prj_meshsz;allocate(ff(rp_meshsz))

!Define q mesh for reciprocal space
 nqgauss1=ngaussmax
 allocate(qgauss1(nqgauss1),wgauss1(nqgauss1))
 call gauleg(zero,gmax,qgauss1,wgauss1,nqgauss1)
 nqgauss2=ngaussmax
 allocate(qgauss2(nqgauss2),wgauss2(nqgauss2))
 call gauleg(gmax,gamma,qgauss2,wgauss2,nqgauss2)
 allocate(chi1g(nqgauss1,pshead%basis_size),chireg(nqmax,pshead%basis_size))
 allocate(am(nqgauss2,nqgauss2),bm(nqgauss2),chi2g(nqgauss2))

!Transfer tproj(r) into chi(q) for 0<q<gmax
!-- On a Gaussian mesh
 do ib=1,pshead%basis_size
  ll=pshead%orbitals(ib)
  do iq=1,nqgauss1
   do ir=1,rp_meshsz
    call jbessel(bess,bessp,dum,ll,0,qgauss1(iq)*pawrad%rad(ir))
    ff(ir)=bess*pawrad%rad(ir)*pawps%tproj(ir,ib)
   enddo
   call csimp(ff,pawrad,rp_meshsz,chi1g(iq,ib))
  enddo
!-- On a regular mesh
  dq=gmax/dble(nqmax-1)
  do iq=1,nqmax
   do ir=1,rp_meshsz
    call jbessel(bess,bessp,dum,ll,0,dble(iq-1)*dq*pawrad%rad(ir))
    ff(ir)=bess*pawrad%rad(ir)*pawps%tproj(ir,ib)
   enddo
   call csimp(ff,pawrad,rp_meshsz,chireg(iq,ib))
  enddo
 enddo

!Loop on error Wl
 do while (wwmax>pawrso%werror)
  wwmax=-one
  if (pawrad%islog==1) then
   r0_meshsz=max(int(log(one+(r0*1.035_dp)/pawrad%rstep)/pawrad%lstep)+1,pshead%prj_meshsz+1)
  else
   r0_meshsz=max(int(r0*1.035_dp/pawrad%rstep)+1,pshead%prj_meshsz+1)
   if (mod(r0_meshsz,2)==0) r0_meshsz=r0_meshsz-1
  end if
  if (r0_meshsz>pshead%prj_msz_max) then
   write(std_out,'(/,2x,a)') 'Error in Atompaw2Abinit(opt_proj): ro_meshsz too big !'
   stop
  endif
  r0=pawrad%rad(r0_meshsz)
  allocate(gg(r0_meshsz))

!Loop on (l,n) basis
 do ib=1,pshead%basis_size
  ll=pshead%orbitals(ib)

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
    call jbessel(bess,bessp,dum,ll,0,qgauss1(iq)*pawrad%rad(ir))
    xx=xx+wgauss1(iq)*bess*chi1g(iq,ib)*qgauss1(iq)**2
   enddo
   do iq=1,nqgauss2
    call jbessel(bess,bessp,dum,ll,0,qgauss2(iq)*pawrad%rad(ir))
    xx=xx+wgauss2(iq)*bess*chi2g(iq)*qgauss2(iq)**2
   enddo
   pawps%tproj(ir,ib)=two*pawrad%rad(ir)*xx/pi
  enddo

!Estimate the error W_l(q)
! --Compute Int(0,R0) [r**2.chi(r).jl(qr)] (and Wl)
! --for each q of the regular mesh
  wwl=-one
  do iq=1,nqmax
   do ir=1,r0_meshsz
    call jbessel(bess,bessp,dum,ll,0,dble(iq-1)*dq*pawrad%rad(ir))
    gg(ir)=bess*pawrad%rad(ir)*pawps%tproj(ir,ib)
   enddo
   call csimp(gg,pawrad,r0_meshsz,xx)
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

 pshead%prj_meshsz=r0_meshsz
 write(std_out,'(2(a,f7.4),a)') '  New radius R0 for nl projectors (Bohr)=',&
&                         r0,' (=',r0/rc_prj,'*Rc(proj))'
 if (r0>1.55_dp*rc_prj) &
&  write(std_out,'(4(/,a))') &
&                'Warning:',&
&                '  Radius for nl projectors (R0) seems to be high !',&
&                '  You should change parameters of Real Space Optimization',&
&                '  (increase Ecut, Gamma/Gmax or Wl).'

 end subroutine opt_proj

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
!! wrpawps
!!
!! FUNCTION
!! Write a PAW dataset file formatted for ABINIT
!!
!! INPUTS
!!  filename= output file name for Abinit
!!  funit= output unit number
!!  loggrd
!!    %meshsz=mesh size for the logarithmic grid
!!    %uselog=TRUE if data are transfered on a log. grid before being written
!!    %log_step=logarithmic step for the logarithmic grid
!!    %rad_step=radial step for the logarithmic grid
!!  pawarray
!!    %shapefunc(wav_meshsz)= Normalized shape function
!!    %shpnrm(l_size)= Moments of shape function for each l
!!  pawps
!!    %coreden4pr2(core_meshsz)= Core density multiplied by 4Pi.r2
!!    %tcoreden4pr2(core_meshsz)= Pseudized core density multiplied by 4Pi.r2
!!    %tvaleden4pr2(vale_meshsz)= Pseudized valence density multiplied by 4Pi.r2
!!    %dij0(lmn2_size)= Part of the Dij term calculated in the psp part
!!    %phi(wav_meshsz,basis_size)= PAW atomic wavefunctions
!!                                on the radial grid
!!    %tphi(wav_meshsz,basis_size)= PAW atomic pseudo-wavefunctions
!!                                 on the radial grid
!!    %tproj(prj_msz_max,basis_size)= PAW projectors
!!                                 on the radial grid
!!    %rhoij0= Atomic initialization of rhoij
!!    %vbare(sph_meshsz)= bare local potential (part of VH(tnzc))
!!    %vhtnzc(core_meshsz)= Hartree potential of the ps-density
!!                           of the nucleus + core electrons
!!  pawrad= radial grid definitions
!!  pawrso
!!    %userso=TRUE if REAL Space Optimization is required
!!  pshead
!!    %atomic_charge= Total atomic charge
!!    %basis_size= Number of elements for the paw nl basis
!!    %core_charge= Core charge
!!    %core_meshsz= Dimension of radial mesh for core density
!!    %creatorid= ID of psp generator (here creatorID=1 !)
!!    %hat_meshsz= Dimension of radial mesh for shape function
!!    %lambda= Lambda in gaussian type g(r)
!!    %lmax= Maximum value of l
!!    %lmn_size= Number of elements for the paw basis
!!    %l_size= Max. value of l+1 leading to a non zero Gaunt coeffs
!!    %mesh_type=  Flag defining ther radial grid type
!!    %orbitals(basis_size)= Quantum number l for each basis function
!!    %prj_meshsz= Dimension of radial mesh for tproj
!!    %pspcod= Psp code number for Abinit (here PAW->pspcod=7 !)
!!    %pspxc_abinit= Abinit s code number for the exchange-correlation
!!    %rad_step= Step corresponding to radial mesh
!!    %rc_proj= Sphere radius for tproj
!!    %rc_hat= radius for shape function
!!    %rc_sph= Default PAW sphere radius
!!    %shape_type= Shape function type
!!    %sigma= Sigma for gaussian type g(r)
!!    %sph_meshsz= Dimension of radial mesh corresponding to PAW spheres
!!    %title= Title for pseudopotential
!!    %vale_meshsz= Dimension of radial mesh for pseudo valence density (0if not present)
!!    %vloc_meshsz= Dimension of radial mesh for vloc=vhtnzc
!!    %sph_meshsz= Dimension of radial mesh for phi, tphi ...
!!    %vlocopt= Option for Vloc
!!    %wav_meshsz= Dimension of radial mesh for Phi, tPhi, ...
!!
!! NOTES
!!  File format of formatted PAW psp input for ABINIT:
!!  --------------------------------------------------
!!  (1) title (character) line
!!  (2) znucl, zion, pspdat
!!  (3) pspcod, pspxc, lmax, lloc, mmax, r2well
!!  (4) psp_version, creatorID
!!  (5) basis_size, lmn_size
!!  (6) orbitals (for l=1 to basis_size)
!!  (7) number_of_meshes
!!  For imsh=1 to number_of_meshes
!!      (8)  mesh_index, mesh_type ,mesh_size, rad_step[, log_step]
!!  (9) r_cut(SPH)
!!  (10) shape_type, r_shape[, shapefunction arguments]
!!  For iln=1 to basis_size
!!      (11) comment(character)
!!      (12) radial mesh index for phi
!!      (13) phi(r) (for ir=1 to phi_meshsz)
!!  For iln=1 to basis_size
!!      (14) comment(character)
!!      (15) radial mesh index for tphi
!!      (16) tphi(r) (for ir=1 to phi_mesh_size)
!!  For iln=1 to basis_size
!!      (17) comment(character)
!!      (18) radial mesh index for tproj
!!      (19) tproj(r) (for ir=1 to proj_mesh_size)
!!  (20) comment(character)
!!  (21) radial mesh index for core_density
!!  (22) core_density (for ir=1 to phi_mesh_size)
!!  (23) comment(character)
!!  (24) radial mesh index for tcore_density
!!  (25) tcore_density (for ir=1 to phi_mesh_size)
!!  (26) comment(character)
!!  (27) Dij0 (for ij=1 to lmn_size*(lmn_size+1)/2)
!!  (28) comment(character)
!!  (29) Rhoij0 (for ij=1 to lmn_size*(lmn_size+1)/2)
!!  (30) comment(character)
!!  (31) radial mesh index for Vloc, format of Vloc (0=Vbare, 1=VH(tnzc) with hat in XC, 2=VH(tnzc) without hat in XC)
!!  (32) Vloc(r) (for ir=1 to vloc_mesh_size)
!!  ===== Following lines only if shape_type=-1 =====
!!  For il=1 to 2*max(orbitals)+1
!!      (33) comment(character)
!!      (34) radial mesh index for shapefunc
!!      (35) shapefunc(r)*gnorm(l)*r**l (for ir=1 to phi_meshsz)
!!  --------------------------------------------------
!!
!! PARENTS
!! atompaw2abinit
!!
!! CHILDREN
!! date_and_time,extrapolate,interpfunc,grid2pawrad
!!=================================================================

 subroutine wrpawps(pshead,pawps,pawarray,pawrad,loggrd,fname,funit,author)

 type(pshead_type),intent(in)    :: pshead
 type(pawps_type),intent(in)     :: pawps
 type(pawarray_type),intent(in)  :: pawarray
 type(pawrad_type),intent(inout) :: pawrad
 type(loggrd_type),intent(inout) :: loggrd
 character*(fnlen),intent(in)    :: author,fname
 integer,intent(in)              :: funit

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer, parameter :: nmesh_max=5
 integer :: coremeshsz,ib,icoremesh,ii,il,ilmn,iprjmesh,ir
 integer :: ivalemesh,ivlocmesh,iwavmesh,jlmn,mtyp,nmesh
 integer :: prjmeshsz,valemeshsz,vlocmeshsz,wavmeshsz
 real(dp) :: rstep,lstep
 character*4 :: pspfmt
 character*8 :: strdate
 integer :: meshsz(nmesh_max),meshtp(nmesh_max)
 real(dp) :: radstp(nmesh_max),logstp(nmesh_max)
 real(dp), allocatable :: ffit(:),ntmp(:),rr_log(:),shpf(:)
 type(GridInfo) :: Grid_core,Grid_vale

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 pspfmt="paw5"
 if (pshead%vale_meshsz==0) pspfmt="paw3"
 if (pshead%vlocopt/=2)     pspfmt="paw4"

 call meshes_def(coremeshsz,icoremesh,iprjmesh,ivalemesh,ivlocmesh,&
 &    iwavmesh,nmesh,nmesh_max,meshtp,meshsz,loggrd,logstp,prjmeshsz,&
 &    rr_log,valemeshsz,vlocmeshsz,wavmeshsz,pshead,pawrad,radstp)

!Some translations AtomPAW -> ABINIT
 call NullifyGrid(Grid_core);call NullifyGrid(Grid_vale)
 if (pshead%core_meshsz>0) call grid2pawrad(Grid_core,pawrad,pshead%core_meshsz,-1)
 if (pshead%vale_meshsz>0) call grid2pawrad(Grid_vale,pawrad,pshead%vale_meshsz,-1)

!Open the file for writing
 open(unit=funit,file=trim(fname),form='formatted',status='unknown')

!Write the header
 call date_and_time(strdate)
 if (len(trim(author))>0) then
   write(funit,'(3a)') trim(pshead%title)," by ",trim(author)
 else
   write(funit,'(a)') trim(pshead%title)
 end if
 
 write(funit,'(1x,f7.3,1x,f7.3,1x,a,14x,a)') &
&      pshead%atomic_charge,&
&      pshead%atomic_charge-pshead%core_charge,&
&      trim(strdate),&
&      " : zatom,zion,pspdat"
 write (funit,'(1x,i2,1x,i7,1x,i2," 0 ",i5," 0.",14x,a)') &
&      pshead%pspcod,&
&      pshead%pspxc_abinit,&
&      pshead%lmax,&
&      wavmeshsz,&
&      " : pspcod,pspxc,lmax,lloc,mmax,r2well"
 write (funit,'(1x,a4,1x,i4,29x,a)') &
&       pspfmt,pshead%creatorid,&
&       " : pspfmt,creatorID"
 write (funit,'(2(1x,i2),33x,a)') &
&       pshead%basis_size,&
&       pshead%lmn_size,&
&       " : basis_size,lmn_size"
 do ib=1,pshead%basis_size
  write (funit,'(1x,i1)',ADVANCE='NO') pshead%orbitals(ib)
 enddo
 if (pshead%basis_size<20) then
  do ib=pshead%basis_size+1,20
   write (funit,'(a)',ADVANCE='NO') '  '
  enddo
 endif
 write (funit,'(a)') ": orbitals"
 write (funit,'(1x,i1,37x,a)') &
&       nmesh," : number_of_meshes"
 do ii=1,nmesh
  if (meshtp(ii)==1) then
   write (funit,'(1x,i1,1x,i1,1x,i4,1x,es22.16,15x,a,i1,a)') &
&         ii,meshtp(ii),meshsz(ii),radstp(ii),&
&         " : mesh ",ii,", type,size,rad_step[,log_step]"
  else
   write (funit,'(1x,i1,1x,i1,1x,i4,1x,es22.16,1x,es22.16,a,i1,a)') &
&         ii,meshtp(ii),meshsz(ii),radstp(ii),logstp(ii),&
&         " : mesh ",ii,", type,size,rad_step[,log_step]"
  endif
 enddo
 write (funit,'(1x,f13.10,25x,a)') &
&       pshead%rc_sph,&
&      " : r_cut(PAW)"
 if (pshead%shape_type==1) then
  write (funit,'(1x,i1,1x," 0.",i2,1x,es20.12,10x,a)') &
&         pshead%shape_type,pshead%lambda,pshead%sigma,&
&        " : shape_type,rshape,lambda,sigma"
 else
  if (pshead%hat_meshsz==pshead%sph_meshsz) then
   write (funit,'(1x,i1," 0.",34x,a)') &
&         pshead%shape_type,&
&         " : shape_type,rshape"
  else
   write (funit,'(1x,i1,1x,f13.10,20x,a)') &
&         pshead%shape_type,pshead%rc_hat,&
&         " : shape_type,rshape"
  endif
 endif

!Write basis functions
 do ib=1,pshead%basis_size
  write(funit,'(a,i1,a)') "===== PHI ",ib,&
&      " =====   [phi(r)=PHI(r)/r*Ylm(th,ph)]"
  write(funit,'(i2,a)') iwavmesh,"  : radial mesh index"
  if (loggrd%uselog) then
   allocate(ffit(wavmeshsz))
   call interpfunc(pshead%wav_meshsz,pawrad%rad(1:pshead%wav_meshsz),&
&                  pawps%phi(:,ib),wavmeshsz,rr_log(1:wavmeshsz),ffit)
   write(funit,'(3(1x,es23.16))') ffit
   deallocate(ffit)
  else
   write(funit,'(3(1x,es23.16))') (pawps%phi(ir,ib),ir=1,wavmeshsz)
  endif
 enddo

!Write pseudo bais functions
 do ib=1,pshead%basis_size
  write(funit,'(a,i1,a)') "===== TPHI ",ib,&
&      " =====   [tphi(r)=TPHI(r)/r*Ylm(th,ph)]"
  write(funit,'(i2,a)') iwavmesh,"  : radial mesh index"
  if (loggrd%uselog) then
   allocate(ffit(wavmeshsz))
   call interpfunc(pshead%wav_meshsz,pawrad%rad(1:pshead%wav_meshsz),&
&                  pawps%tphi(:,ib),wavmeshsz,rr_log(1:wavmeshsz),ffit)
   write(funit,'(3(1x,es23.16))') ffit
   deallocate(ffit)
  else
   write(funit,'(3(1x,es23.16))') (pawps%tphi(ir,ib),ir=1,wavmeshsz)
  endif
 enddo

!Write projectors
 do ib=1,pshead%basis_size
  write(funit,'(a,i1,a)') "===== TPROJECTOR ",ib,&
&      " =====   [tp(r)=TPROJECTOR(r)/r*Ylm(th,ph)]"
  write(funit,'(i2,a)') iprjmesh,"  : radial mesh index"
  write(funit,'(3(1x,es23.16))') (pawps%tproj(ir,ib),ir=1,pshead%prj_meshsz)
 enddo

!Write core density
 write(funit,'(a)') "===== CORE_DENSITY ====="
 write(funit,'(i2,a)') icoremesh,"  : radial mesh index"
 allocate(ntmp(pshead%core_meshsz))
 ntmp(2:pshead%core_meshsz)=pawps%coreden4pr2(2:pshead%core_meshsz) &
&                          /(four*pi*pawrad%rad(2:pshead%core_meshsz)**2)
 call extrapolate(Grid_core,ntmp)
 if (loggrd%uselog) then
  allocate(ffit(coremeshsz))
  call interpfunc(pshead%core_meshsz,pawrad%rad(1:pshead%core_meshsz),&
&                 ntmp,coremeshsz,rr_log(1:coremeshsz),ffit)
  write(funit,'(3(1x,es23.16))') ffit
  deallocate(ffit)
 else
  write(funit,'(3(1x,es23.16))') (ntmp(ir),ir=1,coremeshsz)
 endif
 deallocate(ntmp)

!Write pseudo core density
 write(funit,'(a)') "===== PSEUDO_CORE_DENSITY ====="
 write(funit,'(i2,a)') icoremesh,"  : radial mesh index"
 allocate(ntmp(pshead%core_meshsz))
 ntmp(2:pshead%core_meshsz)=pawps%tcoreden4pr2(2:pshead%core_meshsz) &
&                          /(four*pi*pawrad%rad(2:pshead%core_meshsz)**2)
 call extrapolate(Grid_core,ntmp)
 if (loggrd%uselog) then
  allocate(ffit(coremeshsz))
  call interpfunc(pshead%core_meshsz,pawrad%rad(1:pshead%core_meshsz),&
&                 ntmp,coremeshsz,rr_log(1:coremeshsz),ffit)
  write(funit,'(3(1x,es23.16))') ffit
  deallocate(ffit)
 else
  write(funit,'(3(1x,es23.16))') (ntmp(ir),ir=1,coremeshsz)
 endif
 deallocate(ntmp)

!Write Dij0 and Rhoij0
 write(funit,'(a)') "===== Dij0 ====="
 ii=0
 do jlmn=1,pshead%lmn_size
  write(funit,'(100(1x,es23.16))') (pawps%dij0(ii+ilmn),ilmn=1,jlmn)
  ii=ii+jlmn
 enddo
 write(funit,'(a)') "===== Rhoij0 ====="
 ii=0
 do jlmn=1,pshead%lmn_size
  write(funit,'(100(1x,es23.16))') (pawps%rhoij0(ii+ilmn),ilmn=1,jlmn)
  ii=ii+jlmn
 enddo

!Write Vloc
 if (pshead%vlocopt==0) write(funit,'(a)') "===== Vbare (Vloc(r)) ====="
 if (pshead%vlocopt==1) write(funit,'(a)') "===== VHntZC (Vloc(r)) with compensation charge in XC ====="
 if (pshead%vlocopt==2) write(funit,'(a)') "===== VHntZC (Vloc(r)) without compensation charge in XC ====="
 write(funit,'(i2,1x,i2,a)') ivlocmesh,pshead%vlocopt, &
&   "  : radial mesh index, Vloc format (0=Vbare, 1=VH(tnzc) with hat in XC, 2=VH(tnzc) w/o hat in XC)"
 if (loggrd%uselog) then
  allocate(ffit(vlocmeshsz))
  if (pshead%vlocopt==0) then
   call interpfunc(pshead%sph_meshsz,pawrad%rad(1:pshead%sph_meshsz),&
&                  pawps%vbare,vlocmeshsz,rr_log(1:vlocmeshsz),ffit)
  else
   call interpfunc(pshead%vloc_meshsz,pawrad%rad(1:pshead%vloc_meshsz),&
&                  pawps%vhtnzc,vlocmeshsz,rr_log(1:vlocmeshsz),ffit)
  endif
  write(funit,'(3(1x,es23.16))') ffit
  deallocate(ffit)
 else
  if (pshead%vlocopt==0) then
   write(funit,'(3(1x,es23.16))') (pawps%vbare(ir),ir=1,vlocmeshsz)
  else
   write(funit,'(3(1x,es23.16))') (pawps%vhtnzc(ir),ir=1,vlocmeshsz)
  endif
 endif

!Write (eventually) shape functions
 if (pshead%shape_type==-1) then
  allocate(shpf(wavmeshsz))
  do il=1,pshead%l_size
   write(funit,'(a,i1,a)') "===== SHAPEF (l=",il-1,") ====="
   write(funit,'(a)') " 1  : radial mesh index"
   if (loggrd%uselog) then
    allocate(ffit(wavmeshsz))
    call interpfunc(pshead%wav_meshsz,pawrad%rad(1:pshead%wav_meshsz),pawarray%shapefunc,&
&                   wavmeshsz,rr_log(1:wavmeshsz),ffit)
    if (il==1) then
     shpf(1:wavmeshsz)=ffit(1:wavmeshsz)
    else
     do ir=1,wavmeshsz
      shpf(ir)=ffit(ir)*pawarray%shpnrm(il)*rr_log(ir)**(il-1)
     enddo
    endif
    deallocate(ffit)
   else
    if (il==1) then
     shpf(1:wavmeshsz)=pawarray%shapefunc(1:wavmeshsz)
    else
     do ir=1,wavmeshsz
      shpf(ir)=pawarray%shapefunc(ir)*pawarray%shpnrm(il)*pawrad%rad(ir)**(il-1)
     enddo
    endif
   endif
   write(funit,'(3(1x,es23.16))') (shpf(ir),ir=1,wavmeshsz)
  enddo
  deallocate(shpf)
 endif

!Write pseudo valence density
 if (pshead%vale_meshsz>0) then
  write(funit,'(a)') "===== PSEUDO_VALENCE_DENSITY ====="
  write(funit,'(i2,a)') ivalemesh,"  : radial mesh index"
  allocate(ntmp(pshead%vale_meshsz))
  ntmp(2:pshead%vale_meshsz)=pawps%tvaleden4pr2(2:pshead%vale_meshsz) &
&                           /(four*pi*pawrad%rad(2:pshead%vale_meshsz)**2)
  call extrapolate(Grid_vale,ntmp)
  if (loggrd%uselog) then
   allocate(ffit(valemeshsz))
   call interpfunc(pshead%vale_meshsz,pawrad%rad(1:pshead%vale_meshsz),&
&                  ntmp,valemeshsz,rr_log(1:valemeshsz),ffit)
   write(funit,'(3(1x,es23.16))') ffit
   deallocate(ffit)
  else
   write(funit,'(3(1x,es23.16))') (ntmp(ir),ir=1,valemeshsz)
  endif
  deallocate(ntmp)
 end if

!Close the file and end
 close(funit)

 if (loggrd%uselog) deallocate(rr_log)
 if (pshead%core_meshsz>0) call DestroyGrid(Grid_core)
 if (pshead%vale_meshsz>0) call DestroyGrid(Grid_vale)

 end subroutine wrpawps


!!=================================================================
!! NAME
!! wrcorewf
!!
!! FUNCTION
!! Write a core wave-function file formatted for ABINIT
!!
!! INPUTS
!!  filename= core WF file name for Abinit
!!  funit= output unit number
!!  pshead
!!    %atomic_charge= Total atomic charge
!!    %core_charge= Core charge
!!    %core_size= Number of core states
!!    %corewf_meshsz= Dimension of radial mesh for core WFs
!!    %creatorid= ID of psp generator (here creatorID=1 !)
!!    %lmax_core= Maximum value of l for core states
!!    %lmn_size_core= Number of lmn elements for the core states
!!    %log_step= Logarithmic step corresponding to radial mesh
!!    %mesh_type=  Flag defining ther radial grid type
!!    %orbitals_core(core_size)= Quantum number l for each core wave function
!!    %pspcod= Psp code number for Abinit (here PAW->pspcod=7 !)
!!    %pspxc_abinit= Abinit s code number for the exchange-correlation
!!    %rad_step= Step corresponding to radial mesh
!!    %title= Title for pseudopotential
!!  loggrd
!!    %meshsz=mesh size for the logarithmic grid
!!    %uselog=TRUE if data are transfered on a log. grid before being written
!!    %log_step=logarithmic step for the logarithmic grid
!!    %rad_step=radial step for the logarithmic grid
!!  pawrad= radial grid definitions
!!  FC= definition for core states (from AtomPAW)
!!  AEOrbit= definition of AE orbitals (from AtomPAW)
!!
!! NOTES
!!  File format of formatted core orbitals for ABINIT:
!!  --------------------------------------------------
!!  (1) title (character) line
!!  (2) method, nspinor,nsppol
!!  (3) znucl, zion, pspdat
!!  (4) ixc, lmax
!!  (5) version, creatorID
!!  (6) ln_size_core, lmn_size_core
!!  (7) orbitals_core (for i=1 to ln_size)
!!  (8) number_of_meshes
!!  For imsh=1 to number_of_meshes
!!    (9) mesh_index, mesh_type ,mesh_size, rad_step[, log_step]
!!  (10) rcore (SPH)
!!  For iln=1 to ln_size_core
!!    (11) comment(character)
!!    (12) radial mesh index for phi
!!    (13) nn, ll, spin
!!    (14) phi_core(r) (for ir=1 to phi_core_meshsz)
!!
!!Comments:
!! * allowed values for method are:
!!    1 for restricted, compatible only with nsppol=1.
!!    2 for spin unrestricted, compatible only with nsppol=2.
!!* psp_version= ID of core WF file version
!!    4 characters string of the form 'coren' (with n varying)
!!* creatorID= ID of psp generator
!!  creatorid=1xyz : psp generated from Holzwarth AtomPAW generator version x.yz
!!* mesh_type= type of radial mesh
!!    mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!    mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!    mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!    mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!  --------------------------------------------------
!!
!! PARENTS
!! atompaw2abinit
!!
!! CHILDREN
!! date_and_time,interpfunc
!!=================================================================

 subroutine wrcorewf(AEOrbit,FC,pshead,pawrad,loggrd,fname,funit)

 type(OrbitInfo), intent(in)  :: AEOrbit
 type(FCInfo), intent(in)     :: FC
 type(pshead_type),intent(in) :: pshead
 type(pawrad_type),intent(in) :: pawrad
 type(loggrd_type),intent(in) :: loggrd
 character*(fnlen),intent(in)   :: fname
 integer,intent(in)             :: funit

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer, parameter :: nmesh_max=2
 integer :: ib,icor,ii,imesh,ir,isppol,method,nmesh,nspinor,nsppol
 integer :: meshsz(nmesh_max),meshtp(nmesh_max)
 real(dp) :: radstp(nmesh_max),logstp(nmesh_max)
 character*3 :: spstrg
 character*5 :: pspfmt
 character*8 :: strdate
 real(dp), allocatable :: ffit(:),rr_log(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 pspfmt="core1"

 call date_and_time(strdate)

!Hard-coded values, spin unrestricted, spinors or collinear magnetism not yet supported
 method=1
 nspinor=1
 nsppol=1

!Mesh(es) definitions
 nmesh=1;imesh=1
!--- Use of an auxilliary log grid
 if (loggrd%uselog) then
  meshtp(1)=3
  meshsz(1)=loggrd%meshsz
  radstp(1)=(pawrad%rad(pshead%corewf_meshsz)*(one-tol12))*exp(-loggrd%log_step*dble(loggrd%meshsz-2))
  logstp(1)=loggrd%log_step
  allocate(rr_log(meshsz(1)));rr_log(1)=zero
  do ir=2,meshsz(1)
   rr_log(ir)=radstp(1)*exp(logstp(1)*dble(ir-2))
  enddo
 else
!--- No use of an auxilliary log grid
  meshtp(1)=pshead%mesh_type
  meshsz(1)=pshead%corewf_meshsz
  radstp(1)=pshead%rad_step
  logstp(1)=pshead%log_step
 endif

!Open the file for writing
 open(unit=funit,file=trim(fname),form='formatted',status='unknown')

!Write the header
 write(funit,'(2a)') "All-electron core wavefunctions - ",trim(pshead%title)
 write(funit,'(3(1x,i2),30x,a)') method,nspinor,nsppol," : method,nspinor,nsppol"
 write(funit,'(1x,f7.3,1x,f7.3,1x,a,14x,a)') &
&      pshead%atomic_charge,pshead%core_charge,&
&      trim(strdate),&
&      " : zatom,zcore,pspdat"
 write (funit,fmt='(1x,i2,1x,i7,1x,i2,25x,a)') &
&      pshead%pspcod,&
&      pshead%pspxc_abinit,&
&      pshead%lmax_core,&
&      " : pspcod,pspxc,lmax"
 write (funit,'(1x,a5,1x,i4,28x,a)') &
&       pspfmt,pshead%creatorid,&
&       " : pspfmt,creatorID"
 write(funit,'(2(1x,i2),33x,a)') &
&       pshead%core_size,pshead%lmn_size_core,&
&       " : norb_core, lmn_size"
 do ib=1,pshead%core_size
  write (funit,'(1x,i1)',ADVANCE='NO') pshead%orbitals_core(ib)
 enddo
 if (pshead%core_size<20) then
  do ib=pshead%core_size+1,20
   write (funit,'(a)',ADVANCE='NO') '  '
  enddo
 endif
 write (funit,'(a)') ": core_orbitals"
 write (funit,'(2x,i1,36x,a)') &
&       nmesh," : number_of_meshes"
 do ii=1,nmesh
  if (meshtp(ii)==1) then
   write (funit,'(1x,i1,1x,i1,1x,i4,1x,es22.16,15x,a,i1,a)') &
&         ii,meshtp(ii),meshsz(ii),radstp(ii),&
&         " : mesh ",ii,", type,size,rad_step[,log_step]"
  else
   write (funit,'(1x,i1,1x,i1,1x,i4,1x,es22.16,1x,es22.16,a,i1,a)') &
&         ii,meshtp(ii),meshsz(ii),radstp(ii),logstp(ii),&
&         " : mesh ",ii,", type,size,rad_step[,log_step]"
  endif
 enddo
 write (funit,'(1x,f14.10,24x,a)') &
&       pawrad%rad(pshead%corewf_meshsz),&
&      " : r_max(CORE)"

!Write the core wave functions
 icor=0
 do ib=1,AEOrbit%norbit
  if (AEOrbit%iscore(ib)) then
   icor=icor+1;if (icor>pshead%core_size) stop "Atompaw : bug in wrcorewf !"
   do isppol=1,nsppol
    if (nsppol==2) then
     if (isppol==1) spstrg=" UP"
     if (isppol==2) spstrg=" DN"
     write(funit,'(a,i1,2a)') "===== Core wave functions PHI ",icor,&
&        spstrg," =====   [phi(r)=PHI(r)/r*Ylm(th,ph)] "
    else
     write(funit,'(a,i1,a)') "===== Core wave functions PHI ",icor,&
&        " =====   [phi(r)=PHI(r)/r*Ylm(th,ph)]"
    end if
    write(funit,'(i2,a)') imesh,"  : radial mesh index"
    write(funit,'(1x,3(i3,1x),26x,a)') &
&       AEOrbit%np(ib),AEOrbit%l(ib),isppol,&
&       " : n,l,spin"
    write(funit,'((1x,2(es16.8,1x),4x,a))') &
&       AEOrbit%eig(ib),AEOrbit%occ(ib),&
&       " : ene,occ"
    if (loggrd%uselog) then
     allocate(ffit(meshsz(1)))
     call interpfunc(pshead%corewf_meshsz,pawrad%rad(1:pshead%corewf_meshsz),&
&                    AEOrbit%wfn(1:pshead%corewf_meshsz,ib),&
&                    meshsz(1),rr_log(1:meshsz(1)),ffit)
     write(funit,'(3(1x,es23.16))') ffit
     deallocate(ffit)
    else
     write(funit,'(3(1x,es23.16))') (AEOrbit%wfn(ir,ib),ir=1,meshsz(1))
    endif
   end do ! isppol
  end if ! if icore
 end do   !ib

!Close the file and end
 close(funit)
 if (loggrd%uselog) deallocate(rr_log)

 end subroutine wrcorewf


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
!! grid2pawrad
!!
!! FUNCTION
!! Convert a grid datastructure from ABINIT pawrad format to AomPAW Grid format
!! or the contrary
!!
!! INPUTS
!!  nn= size of mesh to be copied (can be smaller as default mesh size)
!!      if negative or zero, copy entire grid
!!  idir= direction of translation: +1: Grid->pawrad, -1: pawrad->Grid
!!
!! SIDE EFFECTS
!!  Grid=  grid datastructure in AtomPAW format
!!  pawrad=grid datastructure in ABINIT  format
!!
!! PARENTS
!! calc_dij0,rdpawps1
!!=================================================================

 subroutine grid2pawrad(Grid,pawrad,nn,idir)

 integer,intent(in) :: nn,idir
 type(GridInfo),intent(inout) :: Grid
 type(pawrad_type),intent(inout) :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ir,meshsz

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (idir==1) then

  if (nn>Grid%n) then
   write(std_out,'(/,2x,a)') "Error in Atompaw2Abinit(grid2pawrad): mesh size too large !"
  end if
  if (nn> 0) meshsz=nn
  if (nn<=0) meshsz=Grid%n

  pawrad%meshsz=meshsz
  if (usingloggrid(Grid)) then
   pawrad%islog=1
   pawrad%rstep=Grid%drdu(1)
   pawrad%lstep=Grid%h
  else
   pawrad%islog=0
   pawrad%rstep=Grid%h
   pawrad%lstep=zero
  end if
  if (associated(pawrad%rad)) deallocate(pawrad%rad)
  if (associated(pawrad%radfact)) deallocate(pawrad%radfact)
  allocate(pawrad%rad(meshsz),pawrad%radfact(meshsz))
  pawrad%rad    (1:meshsz)=Grid%r   (1:meshsz)
  pawrad%radfact(1:meshsz)=Grid%drdu(1:meshsz)

 else if (idir==-1) then

  if (nn>pawrad%meshsz) then
   write(std_out,'(/,2x,a)') "Error in Atompaw2Abinit(grid2pawrad): mesh size too large !"
  end if
  if (nn> 0) meshsz=nn
  if (nn<=0) meshsz=pawrad%meshsz
  Grid%n=meshsz

  if (pawrad%islog==1) then
   Grid%type=2
   Grid%h=pawrad%lstep
   Grid%ishift=5
  else
   Grid%type=1
   Grid%h=pawrad%rstep
   Grid%ishift=25
  end if
  if (associated(Grid%r)) deallocate(Grid%r)
  if (associated(Grid%drdu)) deallocate(Grid%drdu)
  allocate(Grid%r(meshsz),Grid%drdu(meshsz))
  Grid%r   (1:meshsz)=pawrad%rad    (1:meshsz)
  Grid%drdu(1:meshsz)=pawrad%radfact(1:meshsz)
  if (pawrad%islog==1) then
   if (associated(Grid%pref)) deallocate(Grid%pref)
   if (associated(Grid%rr02)) deallocate(Grid%rr02)
   allocate(Grid%pref(meshsz),Grid%rr02(meshsz))
   do ir=1,meshsz
    Grid%pref(ir)=pawrad%rstep*exp(pawrad%lstep*(ir-1)*half)
    Grid%rr02(ir)=(pawrad%rad(ir)+pawrad%rstep)**2
   end do
  end if

 end if

 end subroutine grid2pawrad


!!=================================================================
!! NAME
!! csimp
!!
!! FUNCTION
!! Do integral using corrected Simpson rule (on a linear or logarithmic grid)
!! (exactly like in ABINIT)
!!
!! INPUTS
!!  ff(meshsz)=integrand values
!!  meshsz=size of radial mesh for integration
!!  pawrad
!!    %rstep= Logarithmic step corresponding to radial mesh
!!    %lstep= Logarithmic step corresponding to radial mesh
!!    %radfact(max_meshsz)= Factor used to compute radial integrals on generalized grid
!!
!! OUTPUT
!!  simp=resulting integral by corrected Simpson rule
!!
!! PARENTS
!! calc_dij0,calc_valden,calc_shapef
!!=================================================================

 subroutine csimp(ff,pawrad,meshsz,simp)

 integer,intent(in)           :: meshsz
 real(dp),intent(out)         :: simp
 real(dp),intent(in)          :: ff(meshsz)
 type(pawrad_type),intent(in) :: pawrad

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 real(dp),parameter :: eps=tol14**4
 integer :: ii,nn
 real(dp) :: hh,simp1,simp2,simp4

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (pawrad%islog==1) then
  hh=pawrad%lstep
 else
  hh=pawrad%rstep
 end if

 ii=meshsz;do while(abs(ff(ii))<eps.and.ii>0);ii=ii-1;enddo
 if (ii<=0) then
   simp=zero;return
 end if
 nn=min(ii+1,meshsz)

 if (nn>=5) then

  if (mod(nn,2)==1) then
   simp1=ff(1)*pawrad%radfact(1)+ff(nn)*pawrad%radfact(nn)
   simp2=zero;simp4=ff(2)*pawrad%radfact(2)
   do ii=3,nn-1,2
    simp2=simp2+ff(ii)  *pawrad%radfact(ii)
    simp4=simp4+ff(ii+1)*pawrad%radfact(ii+1)
   enddo
  else
   simp1=1.25_dp*ff(1)*pawrad%radfact(1)+three*ff(2)*pawrad%radfact(2) &
&       -0.25_dp*ff(3)*pawrad%radfact(3)+ff(nn)*pawrad%radfact(nn)
   simp2=zero;simp4=ff(3)*pawrad%radfact(3)
   do ii=4,nn-1,2
    simp2=simp2+ff(ii)  *pawrad%radfact(ii)
    simp4=simp4+ff(ii+1)*pawrad%radfact(ii+1)
   enddo
  endif
  simp=hh/three*(simp1+two*simp2+four*simp4)

 else if (nn==4) then
  simp=hh*0.375_dp*(three*(ff(2)*pawrad%radfact(2)+ff(3)*pawrad%radfact(3)) &
&                 +     (ff(1)*pawrad%radfact(1)+ff(4)*pawrad%radfact(4)))
 else if (nn==3) then
  simp=hh/three*(ff(1)*pawrad%radfact(1)+four*ff(2)*pawrad%radfact(2)+ff(3)*pawrad%radfact(3))
 else if (nn==2) then
  simp=hh*half*(ff(1)*pawrad%radfact(1)+ff(2)*pawrad%radfact(2))
 else
  simp=zero
 end if

 end subroutine csimp


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
!! meshes_defs
!!
!! FUNCTION
!! Determine meshes definitions
!! (if necessary define a reduced logarithmic radial grid)
!!
!! INPUTS
!!  pshead
!!  pawrad
!!  loggrd
!!  nmesh_max
!!
!! OUTPUT
!!  nmesh
!!  icoremesh,iprjmesh,ivalemesh,ivlocmesh,iwavmesh
!!  coremeshsz,prjmeshsz,valemeshsz,vlocmeshsz,wavmeshsz
!!  meshsz(nmesh_max),meshtp(nmesh_max)
!!  radstp(nmesh_max),logstp(nmesh_max)
!!
!! PARENTS
!!
!!=================================================================

subroutine meshes_def(coremeshsz,icoremesh,iprjmesh,ivalemesh,ivlocmesh,&
&          iwavmesh,nmesh,nmesh_max,meshtp,meshsz,loggrd,&
&          logstp,prjmeshsz,rr_log,valemeshsz,vlocmeshsz,&
&          wavmeshsz,pshead,pawrad,radstp)

 type(pshead_type),intent(in)   :: pshead
 type(pawrad_type),intent(inout)   :: pawrad
 type(loggrd_type),intent(inout)   :: loggrd
 integer, intent(in) :: nmesh_max
 integer,intent(inout) :: nmesh
 integer,intent(inout) :: icoremesh,iprjmesh,ivalemesh,ivlocmesh,iwavmesh
 integer,intent(inout) :: coremeshsz,prjmeshsz,valemeshsz,vlocmeshsz,wavmeshsz
 integer ,intent(inout):: meshsz(nmesh_max),meshtp(nmesh_max)
 real(dp),intent(inout) :: radstp(nmesh_max),logstp(nmesh_max)

!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ib,ii,il,ilmn,ir
 integer :: jlmn,mtyp
 real(dp) :: rstep,lstep
 type(GridInfo) :: Grid_core,Grid_vale
 real(dp), allocatable ,intent(inout):: rr_log(:)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

  radstp=zero;logstp=zero
 
!--- Use of an auxilliary log grid
 if (loggrd%uselog) then
  loggrd%rad_step=(pawrad%rad(pshead%wav_meshsz)*(one-tol12))*exp(-loggrd%log_step*dble(loggrd%meshsz-2))
  mtyp=3;rstep=loggrd%rad_step;lstep=loggrd%log_step
  wavmeshsz=loggrd%meshsz
  prjmeshsz=pshead%prj_meshsz
  coremeshsz=int((one+tol12)*log(pshead%rad_step*dble(pshead%core_meshsz-1)/loggrd%rad_step)/loggrd%log_step)+2
  if (pshead%vale_meshsz>0) then
   valemeshsz=int((one+tol12)*log(pshead%rad_step*dble(pshead%vale_meshsz-1)/loggrd%rad_step)/loggrd%log_step)+2
  else
   valemeshsz=0
  end if
  if (pshead%vlocopt==0) then
   vlocmeshsz=int((one+tol12)*log(pshead%rad_step*dble(pshead%sph_meshsz-1)/loggrd%rad_step)/loggrd%log_step)+2
  else
   vlocmeshsz=int((one+tol12)*log(pshead%rad_step*dble(pshead%vloc_meshsz-1)/loggrd%rad_step)/loggrd%log_step)+2
  endif
  allocate(rr_log(max(wavmeshsz,coremeshsz,valemeshsz,vlocmeshsz)));rr_log(1)=zero
  do ir=2,max(wavmeshsz,coremeshsz,valemeshsz,vlocmeshsz)
   rr_log(ir)=loggrd%rad_step*exp(loggrd%log_step*dble(ir-2))
  enddo
  write(std_out,'(3(/,2x,a),2(/,2x,a,g10.4),4(/,2x,a,i4))') &
&    'Atompaw2Abinit info:',&
&    '  All quantities (except nl projectors) are transfered',&
&    '  into a logarithmic grid (r(i)=A*exp[B(i-2)])...',&
&    '  Log. grid parameters: rad_step=',loggrd%rad_step,&
&    '                        log_step=',loggrd%log_step,&
&    '                        Size       =',wavmeshsz,&
&    '                        Size (core)=',coremeshsz,&
&    '                        Size (vale)=',valemeshsz,&
&    '                        Size (Vloc)=',vlocmeshsz
 else

!--- No use of an auxilliary log grid
  mtyp=pshead%mesh_type;rstep=pshead%rad_step;lstep=pshead%log_step
  wavmeshsz=min(pshead%sph_meshsz+5,pshead%wav_meshsz)
  prjmeshsz=pshead%prj_meshsz
  coremeshsz=pshead%core_meshsz
  valemeshsz=pshead%vale_meshsz
  if (pshead%vlocopt==0)then
   vlocmeshsz=pshead%sph_meshsz
  else
   vlocmeshsz=pshead%vloc_meshsz
  endif
 endif

!--- Build mesh definitions
 nmesh=1;iwavmesh=1
 meshtp(1)=mtyp;meshsz(1)=wavmeshsz;radstp(1)=rstep;logstp(1)=lstep
 if (loggrd%uselog.or.wavmeshsz/=prjmeshsz) then
  nmesh=nmesh+1;iprjmesh=nmesh;meshsz(nmesh)=pshead%prj_meshsz
  meshtp(nmesh)=pshead%mesh_type;radstp(nmesh)=pshead%rad_step;logstp(nmesh)=pshead%log_step
 else
  iprjmesh=iwavmesh
 endif
 if (wavmeshsz/=coremeshsz) then
  if (prjmeshsz/=coremeshsz) then
   nmesh=nmesh+1;icoremesh=nmesh;meshsz(nmesh)=coremeshsz
   meshtp(nmesh)=mtyp;radstp(nmesh)=rstep;logstp(nmesh)=lstep
  else
   icoremesh=iprjmesh
  endif
 else
  icoremesh=iwavmesh
 endif
 if (wavmeshsz/=vlocmeshsz) then
  if(prjmeshsz/=vlocmeshsz) then
   if(coremeshsz/=vlocmeshsz) then
    nmesh=nmesh+1;ivlocmesh=nmesh;meshsz(nmesh)=vlocmeshsz
    meshtp(nmesh)=mtyp;radstp(nmesh)=rstep;logstp(nmesh)=lstep
   else
    ivlocmesh=icoremesh
   endif
  else
   ivlocmesh=iprjmesh
  endif
 else
  ivlocmesh=iwavmesh
 endif
 if (wavmeshsz/=valemeshsz) then
  if(prjmeshsz/=valemeshsz) then
   if(coremeshsz/=valemeshsz) then
    if(vlocmeshsz/=valemeshsz) then
     nmesh=nmesh+1;ivalemesh=nmesh;meshsz(nmesh)=valemeshsz
     meshtp(nmesh)=mtyp;radstp(nmesh)=rstep;logstp(nmesh)=lstep
    else
     ivalemesh=ivlocmesh
    endif
   else
    ivalemesh=icoremesh
   endif
  else
   ivalemesh=iprjmesh
  endif
 else
  ivalemesh=iwavmesh
 endif

 end subroutine meshes_def

!!=================================================================
END Module ABINITInterface
