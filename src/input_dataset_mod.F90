!! NAME
!!  input_dataset_mod
!!
!! FUNCTION
!!  This module contains objects and routines used to parse the atompaw
!!  input file. The complete input dataset is stored in the `input_dataset`
!!  Fortran data-structure.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Contains the following subroutines
!
! input_dataset_read, input_dataset_free, input_dataset_copy, input_dataset_read_occ,
!  input_dataset_read_abinit, input_dataset_read_xml, input_dataset_read_upf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE input_dataset_mod

 USE io_tools
 USE tools

 IMPLICIT NONE
  
 PRIVATE

!Public functions
 PUBLIC :: input_dataset_read        ! Initialize input_dataset by reading it from file
 PUBLIC :: input_dataset_free        ! Destroy input_dataset data-structure
 PUBLIC :: input_dataset_copy        ! Copy a input data-structure into another
 PUBLIC :: input_dataset_read_occ    ! Read and change modified occupations in input_dataset
 PUBLIC :: input_dataset_read_abinit ! Read ABINIT options and store them in input_dataset
 PUBLIC :: input_dataset_read_xml    ! Read XML options and store them in input_dataset
 PUBLIC :: input_dataset_read_upf    ! Read UPF options and store them in input_dataset

!Public structured datatype
!All data from input file
 TYPE,PUBLIC :: input_dataset_t
   CHARACTER(2) :: atomic_symbol    ! Atomic symbol
   INTEGER      :: atomic_charge    ! Atomic charge
   LOGICAL :: scalarrelativistic    ! Flag activating the scalar relativistic scheme
   LOGICAL :: diracrelativistic     ! Flag activating the Dirac relativistic scheme
   LOGICAL :: usespline             ! Flag for using splinesolver for lda/gga
   INTEGER :: splns=400             ! Spline interpolation grid length
   REAL(8) :: splr0=0.1d0           ! Spline interpolation r0 value
   LOGICAL :: finitenucleus         ! Flag activating finite nucleus model
   INTEGER :: finitenucleusmodel    ! Option for the finite nucleus model
   LOGICAL :: HFpostprocess         ! Option for the post-processing of Hartree-Fock
   LOGICAL :: localizedcoreexchange ! Flag activating the localization of the core exchange (HF)
   CHARACTER(10) :: gridkey         ! Type of radial grid
   REAL(8) :: gridrange             ! Range of the radial grid
   REAL(8) :: gridmatch             ! A matching radius in the radial grid
   INTEGER :: gridpoints            ! Number of points of the radial grid
   INTEGER :: nlogderiv             ! Number of points for the logarithmic derivative calculation
   REAL(8) :: minlogderiv           ! Minimum energy value for the logarithmic derivative calculation
   REAL(8) :: maxlogderiv           ! Maximum energy value for the logarithmic derivative calculation
   LOGICAL :: BDsolve               ! Flag activating the Block-Davidson solver
   LOGICAL :: fixed_zero            ! Flag activating the "fixed zero" exact exchange potential calculation
   INTEGER :: fixed_zero_index      ! Option for "fixed zero" calculation in the exact exchange potential
   CHARACTER(132) :: exctype        ! Exchange-correlation type (string)
   LOGICAL :: needvtau              ! TRUE if Calculation is performed with full kinetic energy functional
   INTEGER :: np(5)                 ! Electronic configuration: number of s,p,d,f,g shells
   INTEGER :: norbit                ! Electronic configuration: number of orbitals
   LOGICAL,ALLOCATABLE :: orbit_iscore(:)  ! Electronic configuration: TRUE for the core orbitals
   INTEGER :: norbit_mod            ! Electronic configuration: number of orbitals with modified occupations
   INTEGER,ALLOCATABLE :: orbit_mod_n(:)   ! Electronic config.: n number of the modified orbital
   INTEGER,ALLOCATABLE :: orbit_mod_l(:)   ! Electronic config.: l number of the modified orbital
   INTEGER,ALLOCATABLE :: orbit_mod_k(:)   ! Electronic config.: kappa number of the modified orbital
   REAL(8),ALLOCATABLE :: orbit_mod_occ(:) ! Electronic config.: occupation of the modified orbital
   INTEGER :: norbit_val            ! Electronic configuration: number of valence orbitals 
   INTEGER,ALLOCATABLE :: orbit_val_n(:)   ! Electronic config.: n number of the valence orbital
   INTEGER,ALLOCATABLE :: orbit_val_l(:)   ! Electronic config.: l number of the valence orbital
   INTEGER,ALLOCATABLE :: orbit_val_k(:)   ! Electronic config.: kappa number of the valence orbital 
   INTEGER :: lmax=-1               ! PAW Basis: maximum l value
   REAL(8) :: rc=0.d0               ! PAW basis: cut-off radius for the augmentation regions
   REAL(8) :: rc_shap=0.d0          ! PAW basis: cut-off radius of the compensation charge shape function
   REAL(8) :: rc_vloc=0.d0          ! PAW basis: matching radius for the local potential
   REAL(8) :: rc_core=0.d0          ! PAW basis: matching radius for the pseudo-core density
   LOGICAL :: coreshapemod=.false.  ! Optional alternative shape for tcoreden and tcoretau
   INTEGER :: nbasis                ! PAW basis : number of basis functions
   INTEGER :: nbasis_add            ! PAW basis: number of additional basis functions (unbound states)
   INTEGER,ALLOCATABLE :: basis_add_l(:)      ! PAW basis: l number for the additional basis func.
   INTEGER,ALLOCATABLE :: basis_add_k(:)      ! PAW basis: kappa number for the additional basis func.
   REAL(8),ALLOCATABLE :: basis_add_energy(:) ! PAW basis: ref. energy for the additional basis func.
   REAL(8),ALLOCATABLE :: basis_func_rc(:)    ! PAW basis: matching radii for the pseudo partial waves
   INTEGER :: projector_type        ! Type of projectors (Bloechl, Vanderbilt, ...)
   INTEGER :: pseudo_type           ! Type of pseudization scheme (Bessel, polynom, ...)
   INTEGER :: ortho_type            ! Type of orthogonalization scheme (Gram-Schmidt, ...)
   INTEGER :: pseudo_polynom2_pdeg  ! Polynom2 projectors: degree of the polynom
   REAL(8) :: pseudo_polynom2_qcut  ! Polynom2 projectors: q-value for Fourier filtering
   INTEGER :: shapefunc_type           ! Compensation shape function type (sinc2, gaussian, ...)
   REAL(8) :: shapefunc_gaussian_param ! Compensation shape function: parameter for gaussian type
   LOGICAL :: shapetcore            ! Flag activating building of tcore cancelling a negative compensation charge
   REAL(8) :: hf_coretol            ! Tolerance for core density (Hartree-Fock only)
   INTEGER :: vloc_type             ! Type of local potential pseudization
   INTEGER :: vloc_l                ! Local potential: l quantum number (MTrouillier, Ultrasoft)
   REAL(8) :: vloc_ene              ! Local potential: reference energy (MTrouillier, Ultrasoft)
   REAL(8) :: vloc_setvloc_coef     ! "SetVloc" local potential: coefficient
   REAL(8) :: vloc_setvloc_rad      ! "SetVloc" local potential: radius
   INTEGER :: vloc_kerker_power(4)  ! "Kerker" locazl potential: polynomial powers
   LOGICAL :: abinit_usexcnhat      ! ABINIT option: TRUE if nhat density (compensation) is included in XC
   LOGICAL :: abinit_prtcorewf      ! ABINIT option: TRUE is printing of core WF is required
   CHARACTER(132) :: abinit_author  ! ABINIT option: author to be printed in the PAW dataset file
   LOGICAL :: abinit_uselog         ! ABINIT option: TRUE if data are transfered on a log. grid
   INTEGER :: abinit_log_meshsz     ! ABINIT option: mesh size for the logarithmic grid
   REAL(8) :: abinit_log_step       ! ABINIT option: logarithmic step for the logarithmic grid
   LOGICAL :: abinit_userso         ! ABINIT option: TRUE if REAL Space Optimization is required
   REAL(8) :: abinit_rso_ecut       ! ABINIT option: Real Space Optimization parameter - plane wave cutoff
   REAL(8) :: abinit_rso_gfact      ! ABINIT option: Real Space Optimization parameter - Gamma/Gmax ratio
   REAL(8) :: abinit_rso_werror     ! ABINIT option: Real Space Optimization parameter: max. error
   LOGICAL :: xml_usexcnhat         ! XML option: TRUE if nhat density (compensation) is included in XC
   LOGICAL :: xml_prtcorewf         ! XML option: TRUE is printing of core WF is required
   CHARACTER(132) :: xml_author     ! XML option: author to be printed in the PAW dataset file
   CHARACTER(132) :: xml_comment    ! XML option: comment to be printed in the header of PAW dataset file
   LOGICAL :: xml_usespl            ! XML option: TRUE if data are interpolated on a log. grid
   INTEGER :: xml_spl_meshsz        ! XML option: mesh size for the reduced grid
   LOGICAL :: xml_userso            ! XML option: TRUE if REAL Space Optimization is required
   REAL(8) :: xml_rso_ecut          ! XML option: Real Space Optimization parameter - plane wave cutoff
   REAL(8) :: xml_rso_gfact         ! XML option: Real Space Optimization parameter - Gamma/Gmax ratio
   REAL(8) :: xml_rso_werror        ! XML option: Real Space Optimization parameter: max. error
   LOGICAL :: xml_uselda12          ! XML option: TRUE if LDA-1/2 potential calculation is required
   INTEGER :: xml_lda12_orb_l       ! XML option: LDA-1/2 parameter - l quantum number of the ionized orbital
   INTEGER :: xml_lda12_orb_n       ! XML option: LDA-1/2 parameter - n quantum number of the ionized orbital
   REAL(8) :: xml_lda12_ion         ! XML option: LDA-1/2 parameter - amount of charge removed from the ionized orbital
   REAL(8) :: xml_lda12_rcut        ! XML option: LDA-1/2 parameter - cut-off radius (in bohr)
   CHARACTER(15) :: xml_lda12_logfile ! XML option: LDA-1/2 parameter - name of the log file
   REAL(8) :: upf_grid_xmin         ! UPF option: minimum radius given by the grid
   REAL(8) :: upf_grid_zmesh        ! UPF option: inverse of the rdial step of the grid
   REAL(8) :: upf_grid_dx           ! UPF option: logarithmic step of the grid
   REAL(8) :: upf_grid_range        ! UPF option: range of the grid
 END TYPE input_dataset_t

!Public variable containing the complete input dataset
 TYPE(input_dataset_t),PUBLIC,SAVE,TARGET :: input_dataset

!Public parameters
 INTEGER,PARAMETER,PUBLIC :: PROJECTOR_TYPE_BLOECHL   = 1
 INTEGER,PARAMETER,PUBLIC :: PROJECTOR_TYPE_VANDERBILT= 2
 INTEGER,PARAMETER,PUBLIC :: PROJECTOR_TYPE_CUSTOM    = 3
 INTEGER,PARAMETER,PUBLIC :: PROJECTOR_TYPE_MODRRKJ   = 4
 INTEGER,PARAMETER,PUBLIC :: PROJECTOR_TYPE_HF        = 5
 INTEGER,PARAMETER,PUBLIC :: PSEUDO_TYPE_BLOECHL      = 1
 INTEGER,PARAMETER,PUBLIC :: PSEUDO_TYPE_POLYNOM      = 2
 INTEGER,PARAMETER,PUBLIC :: PSEUDO_TYPE_POLYNOM2     = 3
 INTEGER,PARAMETER,PUBLIC :: PSEUDO_TYPE_RRKJ         = 4
 INTEGER,PARAMETER,PUBLIC :: PSEUDO_TYPE_BLOECHL_K    = 5
 INTEGER,PARAMETER,PUBLIC :: PSEUDO_TYPE_HF           = 6
 INTEGER,PARAMETER,PUBLIC :: ORTHO_TYPE_GRAMSCHMIDT   = 1
 INTEGER,PARAMETER,PUBLIC :: ORTHO_TYPE_VANDERBILT    = 2
 INTEGER,PARAMETER,PUBLIC :: ORTHO_TYPE_SVD           = 3
 INTEGER,PARAMETER,PUBLIC :: ORTHO_TYPE_HF            = 4
 INTEGER,PARAMETER,PUBLIC :: SHAPEFUNC_TYPE_GAUSSIAN  = 1
 INTEGER,PARAMETER,PUBLIC :: SHAPEFUNC_TYPE_SINC      = 2
 INTEGER,PARAMETER,PUBLIC :: SHAPEFUNC_TYPE_BESSEL    = 3
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_MTROULLIER     = 1
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_ULTRASOFT      = 2
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_BESSEL         = 3
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_SETVLOC        = 4
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_KERKER_EXPF    = 5
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_KERKER_POLY    = 6
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_VPSMATCHNC     = 7
 INTEGER,PARAMETER,PUBLIC :: VLOC_TYPE_VPSMATCHNNC    = 8
 INTEGER,PARAMETER,PUBLIC :: UNKNOWN_TYPE             =-1

 !Default values
 LOGICAL,PARAMETER,PRIVATE :: PRTCOREWF_DEF      =.false.
 LOGICAL,PARAMETER,PRIVATE :: USEXCNHAT_DEF      =.false.
 LOGICAL,PARAMETER,PRIVATE :: USELOG_DEF         =.false.
 LOGICAL,PARAMETER,PRIVATE :: USERSO_DEF         =.false.
 LOGICAL,PARAMETER,PRIVATE :: USELDA12_DEF       =.false.
 INTEGER,PARAMETER,PRIVATE :: LOGGRD_SIZE_DEF    =350
 REAL(8),PARAMETER,PRIVATE :: LOGGRD_STEP_DEF    =0.035d0
 REAL(8),PARAMETER,PRIVATE :: RSO_ECUT_DEF       =10.0d0
 REAL(8),PARAMETER,PRIVATE :: RSO_GFACT_DEF      =2.d0
 REAL(8),PARAMETER,PRIVATE :: RSO_WERROR_DEF     =0.0001d0
 INTEGER,PARAMETER,PRIVATE :: LDA12_ORB_L_DEF    =-1
 INTEGER,PARAMETER,PRIVATE :: LDA12_ORB_N_DEF    =-1
 REAL(8),PARAMETER,PRIVATE :: LDA12_ION_DEF      =0.d0
 REAL(8),PARAMETER,PRIVATE :: LDA12_RCUT_DEF     =0.d0
 CHARACTER(15),PARAMETER,PRIVATE :: lda12_logfile='lda-12.log'
 REAL(8),PARAMETER,PRIVATE :: UPF_DX_DEF         =0.005d0
 REAL(8),PARAMETER,PRIVATE :: UPF_XMIN_DEF       =-9.d0
 REAL(8),PARAMETER,PRIVATE :: UPF_ZMESH_DEF      =1.d0
 REAL(8),PARAMETER,PRIVATE :: UPF_RANGE_DEF      =15.d0
 CHARACTER(132),PARAMETER,PRIVATE :: AUTHOR_DEF  ='', COMMENT_DEF=''

!Private parameters
 INTEGER,PARAMETER,PRIVATE :: norbit_max=50,nbasis_add_max=25
 REAL(8),PARAMETER,PRIVATE :: linrange=50.d0,mxgridlin=20001
 REAL(8),PARAMETER,PRIVATE :: logrange=80.d0,v4logrange=100.d0,mxgridlog=2001
 REAL(8),PARAMETER,PRIVATE :: logder_min=-5.d0,logder_max=4.95d0,logder_pts=200
 REAL(8),PARAMETER,PRIVATE :: polynom2_pdeg_def=4
 REAL(8),PARAMETER,PRIVATE :: polynom2_qcut_def=10.d0
 REAL(8),PARAMETER,PRIVATE :: gausstol_def=1.d-4
 REAL(8),PARAMETER,PRIVATE :: hf_coretol_def=1.d-4

!Private variables
 LOGICAL,PRIVATE,SAVE :: has_to_print=.true.
 LOGICAL,PRIVATE,SAVE :: has_to_ask=.true.

CONTAINS

!!=================================================================
!! NAME
!!  input_dataset_read
!!
!! FUNCTION
!!  Initialize an input_dataset datastructure by reading it from
!!  a file. If file is omitted, then read from standard input.
!!  Note: we only read here data used to generate the PAW dataset,
!!    not data used for the post-processing (output, explore, scfpaw, ...)
!!
!! INPUTS (all optionals)
!!  [inputfile]= name of input file to be read
!!  [echofile]= name of a file to echo input file content
!!  [read_global_data]= if TRUE, read global data (atom, XC, grid, ...) - Default TRUE
!!  [read_elec_data]= if TRUE, read electronic configuration (orbital & occupations) - Default TRUE
!!  [read_coreval_data]= if TRUE, read electronic config (core and valence) - Default TRUE
!!  [read_basis_data]= if TRUE, read basis data (radii, pseudo scheme, ...) - Default TRUE
!!
!! OUTPUT
!!  [input_dt]= datastructure containing the complete input file.
!!              If omitted, then the global public `input_dataset`
!!              is used.
!!
!!=================================================================
 SUBROUTINE input_dataset_read(input_dt,inputfile,echofile,&
&           read_global_data,read_elec_data,read_coreval_data,read_basis_data)

!---- Arguments
 CHARACTER*(*),INTENT(IN),OPTIONAL :: inputfile,echofile
 LOGICAL,INTENT(IN),OPTIONAL :: read_global_data,read_elec_data,&
&                               read_coreval_data,read_basis_data
 TYPE(input_dataset_t),INTENT(INOUT),OPTIONAL,TARGET :: input_dt
 
!---- Local variables
 INTEGER,PARAMETER :: ifunit=111,ecunit=222
 INTEGER,PARAMETER :: nkappa(5)=(/1,2,2,2,2/)
 INTEGER :: input_unit
 INTEGER :: ii,io,nadd,norb,nval,nbl,nn,ik,kk
 INTEGER :: ilin,ilog,inrl,iscl,ipnt,ifin,iend,ihfpp,ilcex,itau
 INTEGER :: igrid,irelat,ilogder,ilogv4,ibd,idirac,ifixz,ll,nstart
 INTEGER :: ispline,isplr0,isplns
 LOGICAL :: has_to_echo
 LOGICAL :: read_global_data_,read_elec_data_,read_coreval_data_,read_basis_data_
 CHARACTER(200) :: inputline,inputword
 !CHARACTER(128) :: exchangecorrelationandgridline
 CHARACTER(256) :: exchangecorrelationandgridline
 CHARACTER(1) :: CHR
 REAL(8) :: x1,x2
 INTEGER :: basis_add_l(nbasis_add_max)
 INTEGER :: basis_add_k(nbasis_add_max)
 REAL(8) :: basis_add_energy(nbasis_add_max)
 TYPE(input_dataset_t),POINTER :: dataset

!------------------------------------------------------------------

!Select datastruture to be read
 IF (PRESENT(input_dt)) THEN
   dataset => input_dt
 ELSE
   dataset => input_dataset
 ENDIF

!Select input file logical unit
 IF (PRESENT(inputfile)) THEN
   input_unit=ifunit
   OPEN(ifunit,file=trim(inputfile),form='formatted')
 ELSE
   input_unit=STD_IN
 END IF

!Check if we read a TTY or a FILE
 has_to_ask=unit_isatty(input_unit)

!Do we echo input file content?
 has_to_echo=PRESENT(echofile)
 IF (has_to_echo) THEN
   OPEN(ecunit,file=trim(echofile),form='formatted')
 END IF
 
!Select which components have to be read
 read_global_data_=.true.;if (PRESENT(read_global_data)) read_global_data_=read_global_data
 read_elec_data_=.true.;if (PRESENT(read_elec_data)) read_elec_data_=read_elec_data
 read_coreval_data_=.true.;if (PRESENT(read_coreval_data)) read_coreval_data_=read_coreval_data
 read_basis_data_=.true.;if (PRESENT(read_basis_data)) read_basis_data_=read_basis_data

!Print a title
 IF (read_global_data_.OR.read_elec_data_.OR.read_coreval_data_.OR.read_basis_data_) THEN
   IF (has_to_ask) THEN
     WRITE(STD_OUT,'(/,1x,a)') "===== READING OF INPUT DATA ====="
   ELSE
     WRITE(STD_OUT,'(/,3x,a)') "===== READING OF INPUT FILE ====="
   END IF
 END IF
 
 
!------------------------------------------------------------------
!Start reading of AE data
 IF (read_global_data_) THEN

!------------------------------------------------------------------
!=== 1st line: read atomic symbol, atomic number

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'Enter atomic symbol and atomic number'
 END IF

 READ(input_unit,'(a)') inputline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
 CALL eliminate_comment(inputline)

 READ(inputline,*) dataset%atomic_symbol,dataset%atomic_charge

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a,a2)') "Atomic symbol : ",dataset%atomic_symbol
   WRITE(STD_OUT,'(3x,a,i0)') "Atomic charge : ",dataset%atomic_charge
 END IF

!------------------------------------------------------------------
!=== 2nd line: read XC type, grid data, relativistic,point-nucleus,
!              logderiv data, HF data, Block-Davidson keyword 

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'Enter exchange-correlation type, among the following:'
   WRITE(STD_OUT,*) '    * LDA-PW (default), GGA-PBE, GGA-PBESOL'
   WRITE(STD_OUT,*) '    * MGGA-R2SCAN-001, MGGA-R2SCAN-01'
   WRITE(STD_OUT,*) '    * LibXC keyword beginning with XC_'
   WRITE(STD_OUT,*) '    * WTAU followed by LibXC keyword beginning with XC_ '
   WRITE(STD_OUT,*) '       for mgga with altered Kohn-Sham equations'
   WRITE(STD_OUT,*) '    * EXX, EXXOCC, EXXKLI, EXXCS'
   WRITE(STD_OUT,*) '    * HF, HFV'
   WRITE(STD_OUT,*) ' further optionally (space) "nonrelativistic/scalarrelativistic/diracrelativistic" keyword'
   WRITE(STD_OUT,*) ' further optionally (space) "point-nucleus/finite-nucleus" keyword'
   WRITE(STD_OUT,*) ' further optionally (space) "loggrid/lineargrid" keyword if appropriate'
   WRITE(STD_OUT,*) ' further optionally (space) "splineinterp" keyword to use splinesolver for lda/gga'
   WRITE(STD_OUT,*) '   Note: "loggridv4" puts more points near origin'
   WRITE(STD_OUT,*) ' further optionally n (number of grid points)'
   WRITE(STD_OUT,*) '                       r_max (max. grid radius)'
   WRITE(STD_OUT,*) '                       r_match (exact value of r(n))'
   WRITE(STD_OUT,*) ' further optionally (space) "logderivrange" keyword'
   WRITE(STD_OUT,*) ' further optionally (space) "splr0xxx" to change default initial spline pt'
   WRITE(STD_OUT,*) ' further optionally (space) "splnsxxx" to change default spline grid'
   WRITE(STD_OUT,*) '  further optionally emin (minimum energy for log. deriv. plot)'
   WRITE(STD_OUT,*) '                     emax (maximum energy for log. deriv. plot)'
   WRITE(STD_OUT,*) '                     ne   (#  of energies for log. deriv. plot)'
   WRITE(STD_OUT,*) ' and optionally "BDsolve" keyword for Block-Davidson solver'
 END IF

!Read full line
 READ(input_unit,'(a)') exchangecorrelationandgridline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(exchangecorrelationandgridline)
 CALL eliminate_comment(exchangecorrelationandgridline)

 CALL Uppercase(exchangecorrelationandgridline)
 exchangecorrelationandgridline=trim(exchangecorrelationandgridline)


!Retrieve keyword indexes
 ilin=0;ilin=0;ilog=0;ilogv4=0;inrl=0;iscl=0;ipnt=0;ifin=0
 ihfpp=0;ilcex=0;igrid=0;irelat=0;ilogder=0;ibd=0;idirac=0
 ispline=0;isplr0=0;isplns=0
 ilin=INDEX(exchangecorrelationandgridline,'LINEARGRID')
 ilog=INDEX(exchangecorrelationandgridline,'LOGGRID')
 ilogv4=INDEX(exchangecorrelationandgridline,'LOGGRIDV4')
 ibd=INDEX(exchangecorrelationandgridline,'BDSOLVE')
 inrl=INDEX(exchangecorrelationandgridline,'NONRELATIVISTIC')
 iscl=INDEX(exchangecorrelationandgridline,'SCALARRELATIVISTIC')
 idirac=INDEX(exchangecorrelationandgridline,'DIRACRELATIVISTIC')
 ipnt=INDEX(exchangecorrelationandgridline,'POINT-NUCLEUS')
 ifin=INDEX(exchangecorrelationandgridline,'FINITE-NUCLEUS')
 ilogder=INDEX(exchangecorrelationandgridline,'LOGDERIVRANGE')
 ihfpp=INDEX(exchangecorrelationandgridline,'HFPOSTPROCESS')
 ilcex=INDEX(exchangecorrelationandgridline,'LOCALIZEDCOREEXCHANGE')
 ifixz=INDEX(exchangecorrelationandgridline,'FIXED_ZERO')
 itau=INDEX(exchangecorrelationandgridline,'WTAU')
 ispline=INDEX(exchangecorrelationandgridline,'SPLINEINTERP')
 isplr0=INDEX(exchangecorrelationandgridline,'SPLR0')
 isplns=INDEX(exchangecorrelationandgridline,'SPLNS')
 igrid=max(ilin,ilog)  !This line may need attention....
 irelat=max(inrl,iscl) !This line may need attention....

!Treat simple logical variables
 dataset%scalarrelativistic=(iscl>0.and.inrl==0)
 dataset%diracrelativistic=(idirac>0.and.inrl==0)
 dataset%usespline=(itau>0.or.ispline>0.and.inrl==0)
 dataset%finitenucleus=(ifin>0.and.ipnt==0)
 dataset%BDsolve=(ibd>0)
 dataset%HFpostprocess=(ihfpp>0)

!Treat finite nucleus option
 dataset%finitenucleusmodel=-1 
 IF (dataset%finitenucleus) THEN
   READ(exchangecorrelationandgridline(ifin+14:ifin+14),'(a)') CHR
   IF (CHR=="2") dataset%finitenucleusmodel=2      
   IF (CHR=="3") dataset%finitenucleusmodel=3      
   IF (CHR=="4") dataset%finitenucleusmodel=4      
   IF (CHR=="5") dataset%finitenucleusmodel=5      
 END IF  

!Treat possible changes to spline grid
  if (isplr0>0) then
    READ(exchangecorrelationandgridline(isplr0+5:),*) dataset%splr0
  end if
  if (isplns>0) then
    READ(exchangecorrelationandgridline(isplns+5:),*) dataset%splns
  end if

!Treat grid data
 dataset%gridkey='LINEAR'
 dataset%gridpoints=mxgridlin
 dataset%gridrange=linrange
 dataset%gridmatch=linrange
 IF (ilog>0.and.ilin==0.and.ilogv4==0) THEN
   dataset%gridkey='LOGGRID'
   dataset%gridpoints=mxgridlog;
   dataset%gridrange=logrange
   dataset%gridmatch=logrange
 END IF
 IF (ilog>0.and.ilin==0.and.ilogv4>0) THEN
   dataset%gridkey='LOGGRID4'
   dataset%gridpoints=mxgridlog;
   dataset%gridrange=v4logrange
   dataset%gridmatch=v4logrange
 END IF
 IF (igrid>0) THEN
   iend=256
   IF (irelat >igrid.and.irelat-1 <iend) iend=irelat -1
   IF (ilogder>igrid.and.ilogder-1<iend) iend=ilogder-1
   IF (ibd>igrid.and.ibd-1<iend) iend=ibd-1
   inputline=""
   IF (ilog>0.and.ilogv4==0.and.iend>igrid+7) &
&    inputline=TRIM(exchangecorrelationandgridline(igrid+7:iend))
   IF (ilog>0.and.ilogv4>0.and.iend>igrid+9) &
&    inputline=TRIM(exchangecorrelationandgridline(igrid+9:iend))
   IF (ilin>0.and.iend>igrid+10) &
&    inputline=TRIM(exchangecorrelationandgridline(igrid+10:iend))
   IF (inputline/="") THEN
     CALL extractword(1,inputline,inputword);inputword=trim(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) dataset%gridpoints
       CALL extractword(2,inputline,inputword);inputword=trim(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) dataset%gridrange
         dataset%gridmatch=dataset%gridrange
         CALL extractword(3,inputline,inputword);inputword=trim(inputword)
         IF (inputword/="") read(inputword,*) dataset%gridmatch
       END IF
     END IF
   END IF
   IF (dataset%gridpoints<=0) STOP "input_dataset: error -- number of grid points should be >0!"
 END IF

!Treat logderiv data
 dataset%minlogderiv=logder_min
 dataset%maxlogderiv=logder_max
 dataset%nlogderiv=logder_pts
 IF (ilogder>0) THEN
   iend=256
   IF (igrid >ilogder.and.igrid-1 <iend) iend=igrid -1
   IF (irelat>ilogder.and.irelat-1<iend) iend=irelat-1
   inputline=""
   IF (iend>ilogder+13) inputline=trim(exchangecorrelationandgridline(ilogder+13:iend))
   IF (inputline/="") THEN
     CALL extractword(1,inputline,inputword);inputword=trim(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) dataset%minlogderiv
       CALL extractword(2,inputline,inputword);inputword=trim(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) dataset%maxlogderiv
         CALL extractword(3,inputline,inputword);inputword=trim(inputword)
         IF (inputword/="") READ(inputword,*) dataset%nlogderiv
       END IF
     END IF
   END IF
 END IF

!Treat XC/HF
 if (itau>0) then     
   READ(unit=exchangecorrelationandgridline(itau+5:),fmt=*) dataset%exctype
 else  
   READ(unit=exchangecorrelationandgridline(1:),fmt=*) dataset%exctype
 endif  
 dataset%needvtau=(itau>0.or.TRIM(dataset%exctype)=='MGGA-R2SCAN-001'.or.TRIM(dataset%exctype)=='MGGA-R2SCAN-01')
 dataset%localizedcoreexchange=(ilcex>0)
 dataset%fixed_zero=(ifixz>0) ; dataset%fixed_zero_index=-1
 IF (dataset%fixed_zero) &
&   READ(unit=exchangecorrelationandgridline(ifixz+10:),fmt=*) dataset%fixed_zero_index

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,2a)')     "Scalar-relativistic calculation : ",MERGE("YES"," NO",dataset%scalarrelativistic)
   WRITE(STD_OUT,'(3x,2a)')     "Dirac-relativistic calculation  : ",MERGE("YES"," NO",dataset%diracrelativistic)
   IF (dataset%usespline) THEN
     WRITE(STD_OUT,'(3x,a)')    "    - Use a spline solver"
   END IF
   WRITE(STD_OUT,'(3x,2a)')     "Exchange-correlation functional : ",TRIM(dataset%exctype)
   WRITE(STD_OUT,'(3x,3a)')     " (mGGA kinetic energy functional : ",MERGE("YES"," NO",dataset%needvtau),")"
   WRITE(STD_OUT,'(3x,2a)')     "Finite-nucleus calculation      : ",MERGE("YES"," NO",dataset%finitenucleus)
   IF (dataset%finitenucleus) THEN
     WRITE(STD_OUT,'(3x,a,i0)') "    - Finite-nucleus model      : ",dataset%finitenucleusmodel
   END IF
   WRITE(STD_OUT,'(3x,2a)')     "Block-Davidson calculation      : ",MERGE("YES"," NO",dataset%BDsolve)
   WRITE(STD_OUT,'(3x,2a)')     "Grid type                       : ",TRIM(dataset%gridkey)
   WRITE(STD_OUT,'(3x,a,i0)')   "Grid size                       : ",dataset%gridpoints
   WRITE(STD_OUT,'(3x,a,f7.3)') "Grid maximum value              : ",dataset%gridrange
   WRITE(STD_OUT,'(3x,a,f7.3)') "Grid imposed value              : ",dataset%gridmatch
   if(dataset%usespline) then
   WRITE(STD_OUT,'(3x,a,f7.3,2x,i0)') "Spline grid r0, ns              :",&
&      dataset%splr0,dataset%splns
   endif
   WRITE(STD_OUT,'(3x,a,i0)')   "Log. derivative, number of pts  : ",dataset%nlogderiv
   WRITE(STD_OUT,'(3x,a,f7.3)') "Log. derivative, min. energy    : ",dataset%minlogderiv
   WRITE(STD_OUT,'(3x,a,f7.3)') "Log. derivative, max. energy    : ",dataset%maxlogderiv
   WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, post-processing   : ",MERGE("YES"," NO",dataset%HFpostprocess)
   WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, localized core ex.: ",MERGE("YES"," NO",dataset%localizedcoreexchange)
   WRITE(STD_OUT,'(3x,2a)')     "Hartree-Fock, fixed zero        : ",MERGE("YES"," NO",dataset%fixed_zero)
   IF (dataset%fixed_zero) THEN
     WRITE(STD_OUT,'(3x,a,i0)') "    - HF fixed zero index       : ",dataset%fixed_zero_index
   END IF
   IF (dataset%BDsolve.and.dataset%gridkey=='LINEAR') THEN
     WRITE(STD_OUT,'(/,3x,a)') "WARNING: BlockDavidson solver works very slowly with linear grid!"
   END IF
END IF


!------------------------------------------------------------------
!End reading of global data. Start reading of electronic configuration data
 ENDIF
 IF (read_elec_data_) THEN


!------------------------------------------------------------------
!=== 3rd line and following: electronic configuration of atom

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'Enter maximum principal quantum numbers for s,p,d,f,g'
 END IF

 READ(input_unit,'(a)') inputline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
 CALL eliminate_comment(inputline)

 READ(inputline,*) dataset%np(1:5)

 DO ll=1,5
   IF(dataset%np(ll)<0) dataset%np(ll)=0
 END DO
 dataset%norbit=dataset%np(1)+max(dataset%np(2)-1,0)+max(dataset%np(3)-2,0) &
&                            +max(dataset%np(4)-3,0)+max(dataset%np(5)-4,0)
 IF (dataset%diracrelativistic) dataset%norbit=2*dataset%norbit-dataset%np(1)

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a,5(1x,i0))') "Max. quantum numbers (s,p,d,f,g): ",dataset%np(1:5)
   WRITE(STD_OUT,'(3x,a,i0)') "Total number of orbitals: ",dataset%norbit
 END IF
 CALL input_dataset_read_occ(dataset%norbit_mod,dataset%orbit_mod_l,&
&                   dataset%orbit_mod_n,dataset%orbit_mod_k,dataset%orbit_mod_occ,&
&                   dataset%np,dataset%diracrelativistic,&
&                   inputfile_unit=input_unit,echofile_unit=ecunit)


!------------------------------------------------------------------
!End reading of electronic data. Start reading of core/valence data
 ENDIF
 IF (read_coreval_data_) THEN


!------------------------------------------------------------------
!=== Core and valence states

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'For each state enter c for core or v for valence'
 END IF

!Read core and valence states
 IF (ALLOCATED(dataset%orbit_iscore)) DEALLOCATE(dataset%orbit_iscore)
 ALLOCATE(dataset%orbit_iscore(dataset%norbit))
 DO io=1,dataset%norbit
   DO
     READ(input_unit,'(a)') inputline
     CALL eliminate_comment(inputline)
     READ(inputline,*) CHR
     IF (CHR=='c'.OR.CHR=='C'.OR.&
&        CHR=='v'.OR.CHR=='V') THEN
       IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
       EXIT
     ELSE
       WRITE(STD_OUT,*) ' >> Please input c or v!'
     END IF
   END DO
   dataset%orbit_iscore(io)=(CHR=='c'.OR.CHR=='C')
 END DO

!Store valence states
 dataset%norbit_val=dataset%norbit-COUNT(dataset%orbit_iscore(:))
 IF (ALLOCATED(dataset%orbit_val_n)) DEALLOCATE(dataset%orbit_val_n)
 IF (ALLOCATED(dataset%orbit_val_l)) DEALLOCATE(dataset%orbit_val_l)
 IF (ALLOCATED(dataset%orbit_val_k)) DEALLOCATE(dataset%orbit_val_k)
 ALLOCATE(dataset%orbit_val_n(dataset%norbit_val))
 ALLOCATE(dataset%orbit_val_l(dataset%norbit_val))
 ALLOCATE(dataset%orbit_val_k(dataset%norbit_val))
 kk=0
 io=0;nval=0
 DO ll=0,4
   nn=dataset%np(ll+1)
   IF (nn>0) THEN
     DO ik=1,MERGE(nkappa(ll+1),1,dataset%diracrelativistic)
       kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
       IF (.NOT.dataset%diracrelativistic) kk=0
       DO ii=1+ll,nn
         io=io+1
         IF (.NOT.dataset%orbit_iscore(io)) THEN
           nval=nval+1
           dataset%orbit_val_n(nval)=ii
           dataset%orbit_val_l(nval)=ll
           dataset%orbit_val_k(nval)=kk    
         END IF  
       END DO
     END DO
   END IF
 END DO
 IF (dataset%norbit_val/=nval) STOP 'input_dataset: bug -- wrong nval!'

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a)') "Core and valence orbitals:"
   IF (.NOT.dataset%diracrelativistic) WRITE(STD_OUT,'(7x,a)') "n l : type"
   IF (dataset%diracrelativistic)      WRITE(STD_OUT,'(7x,a)') "n l kappa : type"
   io=0
   DO ll=0,4
     nn=dataset%np(ll+1)
     IF (nn>0) THEN
       IF (.NOT.dataset%diracrelativistic) THEN
         DO ii=1+ll,nn
           io=io+1
           WRITE(STD_OUT,'(7x,i1,1x,i1,2a)') ii,ll," : ", &
&            MERGE("CORE   ","VALENCE",dataset%orbit_iscore(io))
         END DO
       ELSE
         DO ik=1,nkappa(ll+1) 
           kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
           DO ii=1+ll,nn
             io=io+1
             WRITE(STD_OUT,'(7x,i1,1x,i1,2x,i2,2x,2a)') ii,ll,kk," : ", &
&              MERGE("CORE   ","VALENCE",dataset%orbit_iscore(io))
           END DO
         END DO
       END IF
     END IF
   END DO
 END IF


!------------------------------------------------------------------
!End reading of AE data. Start reading of basis data
 ENDIF
 IF (read_basis_data_) THEN


!------------------------------------------------------------------
!=== Maximum L for basis functions

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'Enter maximum L for basis and projector functions'
 END IF

 READ(STD_IN,'(a)') inputline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
 CALL eliminate_comment(inputline)

 READ(inputline,*) dataset%lmax

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a,i0)') "Basis, maximum L : ",dataset%lmax
 END IF


!------------------------------------------------------------------
!=== Cut-off radii

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'Enter rc [and possibly: rc_shape, rc_vloc, rc_core and coreshapemod]'
 END IF

 READ(STD_IN,'(a)') inputline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
 CALL eliminate_comment(inputline)

 CALL extractword(1,inputline,inputword);inputword=trim(inputword)
 IF (inputword/="") READ(inputword,*) dataset%rc
 IF (dataset%rc<=1.d-12) THEN
   WRITE(STD_OUT,*) 'input_dataset: error -- rc too small ',dataset%rc,'!'
   STOP
 END IF

 dataset%rc_shap=dataset%rc
 dataset%rc_vloc=dataset%rc
 dataset%rc_core=dataset%rc

 CALL extractword(2,inputline,inputword);inputword=trim(inputword)
 IF (inputword/="") THEN
   READ(inputword,*) dataset%rc_shap
   CALL extractword(3,inputline,inputword);inputword=trim(inputword)
   IF (inputword/="") THEN
     READ(inputword,*) dataset%rc_vloc
     CALL extractword(4,inputline,inputword);inputword=trim(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) dataset%rc_core
     ELSE
       WRITE(STD_OUT,*) 'input_dataset: error -- rc(core) is missing!'
       STOP
     END IF
   ELSE
     WRITE(STD_OUT,*) 'input_dataset: error -- rc(Vloc) is missing!'
     STOP
   END IF
   IF (dataset%rc_shap<=1.d-12.OR.dataset%rc_vloc<=1.d-12.OR.&
&      dataset%rc_core<=1.d-12) THEN
     WRITE(STD_OUT,*) 'input_dataset: error -- one rc is too small!'
     STOP
   END IF
   IF (dataset%rc_shap>dataset%rc.OR.dataset%rc_vloc>dataset%rc.OR.&
&      dataset%rc_core>dataset%rc) THEN
     WRITE(STD_OUT,*) 'input_dataset: error -- rc_shape, rc_vloc and rc_core must be <rc!'
     STOP
   END IF
 ENDIF

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a,f7.4)') "Augmentation region radius : ",dataset%rc
   WRITE(STD_OUT,'(3x,a,f7.4)') "Core dens. matching radius : ",dataset%rc_core
   WRITE(STD_OUT,'(3x,a,f7.4)') "Local pot. matching radius : ",dataset%rc_vloc
   WRITE(STD_OUT,'(3x,a,f7.4)') "Compens. shape func radius : ",dataset%rc_shap
 END IF

 !Optional flag for coreshapemod (perhaps useful for r2SCAN functional)
   CALL Uppercase(inputline)
   if (INDEX(inputline,'CORESHAPEMOD')>0) then
           dataset%coreshapemod=.true.
           IF (has_to_print) THEN
             WRITE(STD_OUT,'(3x,a,f7.4)') "Using alt core shape form    "
           ENDIF  
   endif

!------------------------------------------------------------------
!=== Additional basis functions

 nstart=0 ; dataset%nbasis_add=0 ; basis_add_k(:)=0
 DO ll=0,dataset%lmax
   nbl=0
   nadd = MERGE(nkappa(ll+1),1,dataset%diracrelativistic)
   IF (dataset%np(ll+1)>0) THEN
     nbl=COUNT(.NOT.dataset%orbit_iscore(nstart+1:nstart+dataset%np(ll+1)-ll))
     nstart=nstart+dataset%np(ll+1)-ll
   END IF
   DO
     IF (has_to_ask) THEN
       WRITE(STD_OUT,*) 'For l = ',ll,' there are currently ',nbl,' basis functions'
       IF (.NOT.dataset%diracrelativistic) THEN
         WRITE(STD_OUT,*) 'Enter y to add an additional function or n to go to next l'
       ELSE
         WRITE(STD_OUT,*) 'Enter y to add additional functions (one per kappa) or n to go to next l'
       END IF
     END IF
     READ(STD_IN,'(a)') inputline
     IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
     CALL eliminate_comment(inputline)
     READ(inputline,*) CHR
     IF (CHR/='y'.AND.CHR/='Y') THEN
       IF (CHR/='n'.AND.CHR/='N') STOP 'input_dataset: error -- Please enter Y or N!'
       EXIT
     END IF
     dataset%nbasis_add=dataset%nbasis_add+nadd
     IF (dataset%nbasis_add>nbasis_add_max) STOP 'Too many additional basis functions!'
     basis_add_l(dataset%nbasis_add-nadd+1:dataset%nbasis_add)=ll
     IF (dataset%diracrelativistic) THEN
       basis_add_k(dataset%nbasis_add)=-1
       IF (ll/=0) THEN
         basis_add_k(dataset%nbasis_add-1)=ll
         basis_add_k(dataset%nbasis_add)=-(ll+1)
       END IF
     END IF
     IF (has_to_ask) THEN
       IF (.NOT.dataset%diracrelativistic) THEN
         WRITE(STD_OUT,*) 'Enter energy for generalized function'
       ELSE
         IF (ll==0) WRITE(STD_OUT,*) 'Enter 1 energy for generalized function (kappa=-1)'
         IF (ll/=0) WRITE(STD_OUT,*) 'Enter 2 energies for generalized functions (one energy per kappa)'
       END IF
     END IF
     READ(STD_IN,'(a)') inputline
     IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
     CALL eliminate_comment(inputline)
     READ(inputline,*) basis_add_energy(dataset%nbasis_add-nadd+1:dataset%nbasis_add)
   END DO
 END DO

 IF (ALLOCATED(dataset%basis_add_l)) DEALLOCATE(dataset%basis_add_l)
 IF (ALLOCATED(dataset%basis_add_k)) DEALLOCATE(dataset%basis_add_k)
 IF (ALLOCATED(dataset%basis_add_energy)) DEALLOCATE(dataset%basis_add_energy)
 ALLOCATE(dataset%basis_add_l(dataset%nbasis_add))
 ALLOCATE(dataset%basis_add_k(dataset%nbasis_add))
 ALLOCATE(dataset%basis_add_energy(dataset%nbasis_add))
 IF (dataset%nbasis_add>0) THEN
   dataset%basis_add_l(1:dataset%nbasis_add)=basis_add_l(1:dataset%nbasis_add)
   dataset%basis_add_k(1:dataset%nbasis_add)=basis_add_k(1:dataset%nbasis_add)
   dataset%basis_add_energy(1:dataset%nbasis_add)=basis_add_energy(1:dataset%nbasis_add)
 END IF

 dataset%nbasis=COUNT(.NOT.dataset%orbit_iscore(:))+dataset%nbasis_add

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a,i0)') "Initial number of basis functions    : ",dataset%nbasis-dataset%nbasis_add
   WRITE(STD_OUT,'(3x,a,i0)') "Number of additional basis functions : ",dataset%nbasis_add
   WRITE(STD_OUT,'(3x,a,i0)') "Total number of basis functions      : ",dataset%nbasis
   WRITE(STD_OUT,'(3x,a)') "Additional basis functions:"
   IF (.NOT.dataset%diracrelativistic) THEN
     WRITE(STD_OUT,'(7x,a)') "l : energy"
     DO io=1,dataset%nbasis_add
       WRITE(STD_OUT,'(7x,i1,a,f7.4)') dataset%basis_add_l(io)," : " ,dataset%basis_add_energy(io)
     END DO
   ELSE
     WRITE(STD_OUT,'(7x,a)') "l kappa : energy"
     DO io=1,dataset%nbasis_add
       WRITE(STD_OUT,'(7x,i1,2x,i2,2x,a,f7.4)') dataset%basis_add_l(io), &
&          dataset%basis_add_k(io)," : " ,dataset%basis_add_energy(io)
     END DO
   END IF
 END IF


!------------------------------------------------------------------
!=== Projectors, compensation charge shape function, core tolerance

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'Enter "Bloechl", "Vanderbilt", "modrrkj" or "custom" keywords',&
&             ' for projector generation method.'
   WRITE(STD_OUT,*) ' In case of "custom" choice, enter additional (optional) keywords:'
   WRITE(STD_OUT,*) ' - For partial waves pseudization scheme:'
   WRITE(STD_OUT,*) '     "bloechlps", "polynom", "polynom2 p qcut" or "RRKJ"'
   WRITE(STD_OUT,*) ' - For orthogonalization scheme:'
   WRITE(STD_OUT,*) '     "gramschmidtortho" or "vanderbiltortho" or "svdortho"'
   WRITE(STD_OUT,*) ' - Sinc^2 compensation charge shape is set by default,'
   WRITE(STD_OUT,*) '   Gaussian shape can be specified by adding "Gaussian" keyword and tol,'
   WRITE(STD_OUT,*) '   Bessel shape can be specified by adding "Besselshape" keyword'
   WRITE(STD_OUT,*) ' - Optional shapetcore keyword to modify smooth core shape '
 END IF

 READ(STD_IN,'(a)') inputline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
 CALL eliminate_comment(inputline)

 CALL Uppercase(inputline)
 inputline=TRIM(inputline)

 dataset%pseudo_type=PSEUDO_TYPE_BLOECHL
 dataset%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
 dataset%pseudo_polynom2_pdeg=polynom2_pdeg_def
 dataset%pseudo_polynom2_qcut=polynom2_qcut_def
 dataset%shapefunc_type=SHAPEFUNC_TYPE_SINC
 dataset%shapefunc_gaussian_param=gausstol_def
 dataset%hf_coretol=hf_coretol_def

 READ(unit=inputline,fmt=*) inputword

 IF (TRIM(inputword)=='BLOECHL'.OR.TRIM(inputword)=='VNCT') THEN
   dataset%projector_type=PROJECTOR_TYPE_BLOECHL
   dataset%pseudo_type=PSEUDO_TYPE_BLOECHL
   dataset%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
 ELSE IF (TRIM(inputword)=='VNCK') THEN
   dataset%projector_type=PROJECTOR_TYPE_BLOECHL
   dataset%pseudo_type=PSEUDO_TYPE_BLOECHL_K
   dataset%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
 ELSE IF (TRIM(inputword)=='VANDERBILT'.OR.TRIM(inputword)=='VNCTV') THEN
   dataset%projector_type=PROJECTOR_TYPE_VANDERBILT
   dataset%pseudo_type=PSEUDO_TYPE_POLYNOM
   dataset%ortho_type=ORTHO_TYPE_VANDERBILT
 ELSE IF(TRIM(inputword)=='MODRRKJ') THEN
   dataset%projector_type=PROJECTOR_TYPE_MODRRKJ
   dataset%pseudo_type=PSEUDO_TYPE_RRKJ
   dataset%ortho_type=ORTHO_TYPE_VANDERBILT
   IF (INDEX(inputline,'VANDERBILTORTHO')>0) dataset%ortho_type=ORTHO_TYPE_VANDERBILT
   IF (INDEX(inputline,'GRAMSCHMIDTORTHO')>0) dataset%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
   IF (INDEX(inputline,'SVDORTHO')>0) dataset%ortho_type=ORTHO_TYPE_SVD
 ELSE IF (TRIM(inputword)=='CUSTOM') THEN
   dataset%projector_type=PROJECTOR_TYPE_CUSTOM
   IF (INDEX(inputline,'BLOECHLPS')>0) THEN
     dataset%pseudo_type=PSEUDO_TYPE_BLOECHL
     dataset%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
   ELSE IF (INDEX(inputline,'POLYNOM2')>0) THEN
     dataset%pseudo_type=PSEUDO_TYPE_POLYNOM2
     nstart=INDEX(inputline,'POLYNOM2')
     READ(unit=inputline(nstart+8:),fmt=*,err=111,end=111,iostat=io) &
&         dataset%pseudo_polynom2_pdeg,dataset%pseudo_polynom2_qcut
111  CONTINUE
   ELSE IF (INDEX(inputline,'POLYNOM')>0) THEN
     dataset%pseudo_type=PSEUDO_TYPE_POLYNOM
   ELSE IF (INDEX(inputline,'RRKJ')>0) THEN
     dataset%pseudo_type=PSEUDO_TYPE_RRKJ
   END IF
   IF (INDEX(inputline,'VANDERBILTORTHO')>0) dataset%ortho_type=ORTHO_TYPE_VANDERBILT
   IF (INDEX(inputline,'GRAMSCHMIDTORTHO')>0) dataset%ortho_type=ORTHO_TYPE_GRAMSCHMIDT
 END IF

 IF (TRIM(dataset%exctype)=='HF') THEN
   dataset%projector_type=PROJECTOR_TYPE_HF
   dataset%pseudo_type=PSEUDO_TYPE_HF
   dataset%ortho_type=ORTHO_TYPE_HF
   WRITE(STD_OUT,'(3x,a)') '>> You are using HF XC type: pseudo and orthogonalization line will be ignored!'
 END IF

 IF ((dataset%pseudo_type==PSEUDO_TYPE_BLOECHL.OR. &
&     dataset%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&   .AND.dataset%ortho_type==ORTHO_TYPE_VANDERBILT) STOP &
&  'input_dataset: error -- Vanderbilt orthogonalization not compatible with Bloechls projector scheme!'

 IF ((dataset%pseudo_type==PSEUDO_TYPE_BLOECHL.OR. &
&     dataset%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&   .AND.dataset%ortho_type==ORTHO_TYPE_VANDERBILT) STOP &
&  'input_dataset: error -- Vanderbilt orthogonalization not compatible with Bloechls projector scheme!'

 IF ((dataset%projector_type==PROJECTOR_TYPE_BLOECHL) &
&   .AND.dataset%needvtau) STOP &
&   'input_dataset: error -- mGGA not compatible the Bloechl projector scheme!'

 !!!! Hopefully this will never happen             
 IF ((dataset%projector_type==PROJECTOR_TYPE_HF) &
&   .AND.dataset%needvtau) STOP &
&   'input_dataset: error -- mGGA and Hartree-Fock are not compatible!'

 IF ((dataset%pseudo_type==PSEUDO_TYPE_BLOECHL.OR. &
&     dataset%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&   .AND.dataset%needvtau) STOP &
&   'input_dataset: error -- mGGA not compatible the Bloechl pseudization scheme!'

 IF (INDEX(inputline,'SINC2')>0) THEN
   dataset%shapefunc_type=SHAPEFUNC_TYPE_SINC
 ELSE IF (INDEX(inputline,'GAUSSIAN')>0) THEN
   dataset%shapefunc_type=SHAPEFUNC_TYPE_GAUSSIAN
   nstart=INDEX(inputline,'GAUSSIAN')
   READ(unit=inputline(nstart+8:),fmt=*,err=222,end=222,iostat=io) &
&       dataset%shapefunc_gaussian_param
222 CONTINUE
 ELSE IF (INDEX(inputline,'BESSELSHAPE')>0) THEN
   dataset%shapefunc_type=SHAPEFUNC_TYPE_BESSEL
 END IF
 
 nstart=INDEX(inputline,'CORETOL')
 IF (nstart>0) THEN
   READ(unit=inputline(nstart+7:),fmt=*) dataset%hf_coretol
 END IF

 dataset%shapetcore=(INDEX(inputline,'SHAPETCORE')>0)

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a)') "Projectors description:"
   IF (dataset%projector_type==PROJECTOR_TYPE_BLOECHL) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : BLOECHL"
   IF (dataset%projector_type==PROJECTOR_TYPE_VANDERBILT) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : VANDERBILT"
   IF (dataset%projector_type==PROJECTOR_TYPE_MODRRKJ) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : MODRRKJ"
   IF (dataset%projector_type==PROJECTOR_TYPE_CUSTOM) &
&    WRITE(STD_OUT,'(7x,a)') "Type              : CUSTOM"
   IF (dataset%projector_type==PROJECTOR_TYPE_HF) &
&    WRITE(STD_OUT,'(7x,a)') "Type : HARTREE-FOCK"
   IF (dataset%projector_type/=PROJECTOR_TYPE_HF) THEN
     IF (dataset%pseudo_type==PSEUDO_TYPE_BLOECHL) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : BLOECHL"
     IF (dataset%pseudo_type==PSEUDO_TYPE_POLYNOM) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : POLYNOM"
     IF (dataset%pseudo_type==PSEUDO_TYPE_RRKJ) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : RRKJ"
     IF (dataset%pseudo_type==PSEUDO_TYPE_BLOECHL_K) &
&      WRITE(STD_OUT,'(7x,a)') "Pseudization      : BLOECHL KERKER"
     IF (dataset%pseudo_type==PSEUDO_TYPE_POLYNOM2) &
&      WRITE(STD_OUT,'(7x,a,i0,a,es9.3)') "Pseudization      : POLYNOM2, pdeg=", &
&       dataset%pseudo_polynom2_pdeg,", qcut=",dataset%pseudo_polynom2_qcut
     IF (dataset%ortho_type==ORTHO_TYPE_GRAMSCHMIDT) &
&      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : GRAM-SCHMIDT"
     IF (dataset%ortho_type==ORTHO_TYPE_VANDERBILT) &
&      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : VANDERBILT"
     IF (dataset%ortho_type==ORTHO_TYPE_SVD) &
&      WRITE(STD_OUT,'(7x,a)') "Orthogonalisation : SVD"
   END IF
   IF (dataset%shapefunc_type==SHAPEFUNC_TYPE_GAUSSIAN) &
&    WRITE(STD_OUT,'(3x,a,es9.3)') "Compensation charge shape function : GAUSSIAN, tol=",&
&    dataset%shapefunc_gaussian_param
   IF (dataset%shapefunc_type==SHAPEFUNC_TYPE_SINC) &
&    WRITE(STD_OUT,'(3x,a)') "Compensation charge shape function : SINC2"
   IF (dataset%shapefunc_type==SHAPEFUNC_TYPE_BESSEL) &
&    WRITE(STD_OUT,'(3x,a)') "Compensation charge shape function : BESSEL"
   IF (INDEX(inputline,'CORETOL')>0) &
&    WRITE(STD_OUT,'(3x,a,es9.3)') "Core tolerance for Hartree-Fock : ",dataset%hf_coretol
   WRITE(STD_OUT,'(3x,2a)') "Smooth tcore shape (no negative nhat) : ",MERGE("YES"," NO",dataset%shapetcore)
 END IF


!------------------------------------------------------------------
!=== Local pseudopotential

 IF (has_to_ask) THEN
   WRITE(STD_OUT,*) 'To generate the local pseudopotential, this code can use:'
   WRITE(STD_OUT,*) '  1- a Troullier-Martins scheme for specified l value and energy'
   WRITE(STD_OUT,*) '  2- a non norm-conserving pseudopotential scheme for specified l value and energy'
   WRITE(STD_OUT,*) '  3- a simple pseudization scheme using a single spherical Bessel function'
   WRITE(STD_OUT,*) '  4- VPS match with norm-conservation'
   WRITE(STD_OUT,*) '  5- VPS match without norm-conservation'
   WRITE(STD_OUT,*) '  6- Vloc =  VlocCoef*Shapefunc'
   WRITE(STD_OUT,*) 'For choice 1, enter (high) l value and energy e'
   WRITE(STD_OUT,*) 'For choice 2, enter (high) l value, energy e and "ultrasoft"'
   WRITE(STD_OUT,*) 'For choice 3, enter "bessel"'
   WRITE(STD_OUT,*) 'For choice 4, enter l value, energy e and "vpsmatchnc"'
   WRITE(STD_OUT,*) 'For choice 5, enter l value, energy e and "vpsmatchnnc"'
   WRITE(STD_OUT,*) 'For choice 6, enter "setvloc x y" - x is VlocCoef, y is VlocRad'
 END IF

 READ(STD_IN,'(a)') inputline
 IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
 CALL eliminate_comment(inputline)

 call Uppercase(inputline)
 inputline=TRIM(inputline)
 
 dataset%vloc_type=VLOC_TYPE_MTROULLIER
 dataset%vloc_l=-1
 dataset%vloc_ene=0.d0
 dataset%vloc_setvloc_coef=0.d0
 dataset%vloc_setvloc_rad=dataset%rc
 dataset%vloc_kerker_power(:)=0

 IF (INDEX(inputline,'MTROULLIER')>0) THEN
   dataset%vloc_type=VLOC_TYPE_MTROULLIER
 ELSE IF (INDEX(inputline,'ULTRASOFT')>0) THEN
   dataset%vloc_type=VLOC_TYPE_ULTRASOFT
 ELSE IF (INDEX(inputline,'BESSEL')>0) THEN
   dataset%vloc_type=VLOC_TYPE_BESSEL
 ELSE IF (INDEX(inputline,'VPSMATCHNC')>0) THEN
   dataset%vloc_type=VLOC_TYPE_VPSMATCHNC
 ELSE IF (INDEX(inputline,'VPSMATCHNNC')>0) THEN
   dataset%vloc_type=VLOC_TYPE_VPSMATCHNNC
 ELSE IF (INDEX(inputline,'SETVLOC')>0) THEN
   dataset%vloc_type=VLOC_TYPE_SETVLOC
   nstart=INDEX(inputline,'SETVLOC')
   READ(unit=inputline(nstart+8:),fmt=*,err=333,end=333,iostat=io) x1,x2
   IF (x1<1.d3.AND.x1>-1.d3) dataset%vloc_setvloc_coef=x1
   IF (x2>1.d-8.AND.x2<dataset%rc) dataset%vloc_setvloc_rad=x2
333  CONTINUE
 ELSE IF (INDEX(inputline,'KERKER')>0.OR.dataset%pseudo_type==PSEUDO_TYPE_BLOECHL_K) THEN
   IF (INDEX(inputline,'EXPF')>0) THEN
     dataset%vloc_type=VLOC_TYPE_KERKER_EXPF
     nstart=INDEX(inputline,'EXPF')
   ELSE IF (INDEX(inputline,'POLY')>0) THEN
     dataset%vloc_type=VLOC_TYPE_KERKER_POLY
     nstart=INDEX(inputline,'POLY')
   ELSE
     STOP "EXPF or POLY keyword missing!"
   END IF   
   READ(unit=inputline(nstart+5:),fmt=*,err=334,end=334,iostat=io) &
&    dataset%vloc_kerker_power(1:4)
334  CONTINUE
 END IF

 IF ((dataset%vloc_type==VLOC_TYPE_SETVLOC.OR. &
&     dataset%vloc_type==VLOC_TYPE_KERKER_EXPF.OR. &
&     dataset%vloc_type==VLOC_TYPE_KERKER_POLY) &
&   .AND.dataset%needvtau) STOP &
&   'input_dataset: error -- mGGA not compatible the chosen Vloc scheme!'

 IF (dataset%vloc_type==VLOC_TYPE_MTROULLIER.OR. &
&    dataset%vloc_type==VLOC_TYPE_VPSMATCHNC.OR. &
&    dataset%vloc_type==VLOC_TYPE_VPSMATCHNNC.OR. &
&    dataset%vloc_type==VLOC_TYPE_ULTRASOFT) THEN
   READ(unit=inputline,fmt=*,err=444,end=444,iostat=io) dataset%vloc_l,dataset%vloc_ene
444  CONTINUE
   IF (dataset%vloc_l<0.or.dataset%vloc_l>10) STOP 'input_dataset: error while reading Vloc parameters!'
 END IF

!Print read data
 IF (has_to_print) THEN
   IF (dataset%vloc_type==VLOC_TYPE_MTROULLIER) &
&    WRITE(STD_OUT,'(7x,a,i0,a,f7.4)') "Local pseudopotential type : MTROULLIER, l=",&
&          dataset%vloc_l,", energy=",dataset%vloc_ene
   IF (dataset%vloc_type==VLOC_TYPE_ULTRASOFT) &
&    WRITE(STD_OUT,'(7x,a,i0,a,f7.4)') "Local pseudopotential type : ULTRASOFT, l=",&
&          dataset%vloc_l,", energy=",dataset%vloc_ene
   IF (dataset%vloc_type==VLOC_TYPE_BESSEL) &
&    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : BESSEL"
   IF (dataset%vloc_type==VLOC_TYPE_VPSMATCHNC) &
&    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNC"
   IF (dataset%vloc_type==VLOC_TYPE_VPSMATCHNNC) &
&    WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNNC"
   IF (dataset%vloc_type==VLOC_TYPE_SETVLOC) THEN
     WRITE(STD_OUT,'(7x,a,es9.4,a,es9.4)') "Local pseudopotential type : SETVLOC, coef=",&
&          dataset%vloc_setvloc_coef,", rad=",dataset%vloc_setvloc_rad
     IF (dataset%needvtau) THEN
       write(STD_OUT,*) 'SETVLOC  option not available for MGGA'
         stop
     ENDIF     
   ENDIF  
   IF (dataset%vloc_type==VLOC_TYPE_KERKER_EXPF) &
&    WRITE(STD_OUT,'(7x,a,4(1x,i0))') "Local pseudopotential type : KERKER EXPF, powers=",&
&          dataset%vloc_kerker_power(1:4)
   IF (dataset%vloc_type==VLOC_TYPE_KERKER_POLY) &
&    WRITE(STD_OUT,'(7x,a,4(1x,i0))') "Local pseudopotential type : KERKER POLY, powers=",&
&          dataset%vloc_kerker_power(1:4)
!!!!!!!! removed by NAWH 10-09-2023
!!!   IF (dataset%vloc_type==VLOC_TYPE_MTROULLIER.AND.dataset%needvtau) THEN
!!!     WRITE(STD_OUT,'(7x,a)') 'NOTE: MTROULLIER Vloc not available for mGGA!'
!!!     WRITE(STD_OUT,'(7x,a)') '      Calling VPSmatch with norm conservation instead.'
!!!     dataset%vloc_type=VLOC_TYPE_VPSMATCHNC
!!!     WRITE(STD_OUT,'(7x,a)') "Local pseudopotential type : VPS MATCHNC"
!!!   END IF
 END IF


!------------------------------------------------------------------
!=== Matching radii for the basis functions

!Not for all choice of projectors
 IF (dataset%projector_type==PROJECTOR_TYPE_CUSTOM.OR.&
&    dataset%projector_type==PROJECTOR_TYPE_VANDERBILT.OR.&
&    dataset%projector_type==PROJECTOR_TYPE_MODRRKJ.OR.&
&    dataset%projector_type==PROJECTOR_TYPE_HF.AND.dataset%vloc_type==VLOC_TYPE_MTROULLIER) THEN

   IF (ALLOCATED(dataset%basis_func_rc)) DEALLOCATE(dataset%basis_func_rc)
   ALLOCATE(dataset%basis_func_rc(dataset%nbasis))

   IF (has_to_ask) WRITE(STD_OUT,*) 'For each of the following basis functions enter rc'

   norb=0
   DO ll=0,dataset%lmax
     DO ik=1,MERGE(nkappa(ll+1),1,dataset%diracrelativistic)
       kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
       IF (.NOT.dataset%diracrelativistic) kk=0
       DO io=1,dataset%norbit_val
         IF (dataset%orbit_val_l(io)==ll.AND. &
&           ((.NOT.dataset%diracrelativistic).OR.dataset%orbit_val_k(io)==kk)) THEN     
           norb=norb+1     
           IF (has_to_ask.AND.(.NOT.dataset%diracrelativistic)) &
&            WRITE(STD_OUT,'(a,i2,a,2i4)') '  rc for basis function ',norb,' - n,l= ',norb,ll
           IF (has_to_ask.AND.dataset%diracrelativistic) &
&            WRITE(STD_OUT,'(a,i2,a,3i4)') '  rc for basis function ',norb,' - n,l,kappa= ',norb,ll,kk
           READ(STD_IN,'(a)') inputline
           IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
           CALL eliminate_comment(inputline)
           READ(inputline,*) dataset%basis_func_rc(norb)
         END IF
       END DO
       IF (dataset%nbasis_add>0) THEN
         DO io=1,dataset%nbasis_add
           IF (dataset%basis_add_l(io)==ll.AND. &
&             ((.NOT.dataset%diracrelativistic).OR.dataset%basis_add_k(io)==kk)) THEN     
             norb=norb+1     
             IF (has_to_ask.AND.(.NOT.dataset%diracrelativistic)) &
&              WRITE(STD_OUT,'(a,i2,a,2i4)') '  rc for basis function ',norb,' - n,l= ',999,ll
             IF (has_to_ask.AND.dataset%diracrelativistic) &
&              WRITE(STD_OUT,'(a,i2,a,3i4)') '  rc for basis function ',norb,' - n,l,kappa= ',999,ll,kk
             READ(STD_IN,'(a)') inputline
             IF (has_to_echo) WRITE(ecunit,'(a)') TRIM(inputline)
             CALL eliminate_comment(inputline)
             READ(inputline,*) dataset%basis_func_rc(norb)
           END IF
         END DO
       END IF
     END DO
   END DO

   IF (dataset%nbasis/=norb) STOP 'input_dataset: error -- inconsistency in the number of basis functions!'

!  Print read data
   IF (has_to_print) THEN
     WRITE(STD_OUT,'(3x,a)') "Matching radius for basis functions:"
     IF (.NOT.dataset%diracrelativistic) WRITE(STD_OUT,'(7x,a)') " # - n l : radius"
     IF (dataset%diracrelativistic) WRITE(STD_OUT,'(7x,a)') " # - n l kappa : radius"
     norb=0
     DO ll=0,dataset%lmax
       DO ik=1,MERGE(nkappa(ll+1),1,dataset%diracrelativistic)
         kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
         IF (.NOT.dataset%diracrelativistic) kk=0
         DO io=1,dataset%norbit_val
           IF (dataset%orbit_val_l(io)==ll.AND. &
&             ((.NOT.dataset%diracrelativistic).OR.dataset%orbit_val_k(io)==kk)) THEN     
             norb=norb+1
             IF (.NOT.dataset%diracrelativistic) &
&              WRITE(STD_OUT,'(7x,i2,a,i1,1x,i1,a,f7.4)') &
&              norb," - ",dataset%orbit_val_n(io),ll," : ",dataset%basis_func_rc(norb)
             IF (dataset%diracrelativistic) &
&              WRITE(STD_OUT,'(7x,i2,a,i1,1x,i1,2x,i2,2x,a,f7.4)') &
&              norb," - ",dataset%orbit_val_n(io),ll,kk," : ",dataset%basis_func_rc(norb)
           END IF
         END DO
         IF (dataset%nbasis_add>0) THEN
           DO io=1,dataset%nbasis_add
             IF (dataset%basis_add_l(io)==ll.AND. &
&             ((.NOT.dataset%diracrelativistic).OR.dataset%basis_add_k(io)==kk)) THEN     
               norb=norb+1
               IF (.NOT.dataset%diracrelativistic) &
&                WRITE(STD_OUT,'(7x,i2,a,a1,1x,i1,a,f7.4)') &
&                norb," - ",".",ll," : ",dataset%basis_func_rc(norb)
               IF (dataset%diracrelativistic) &
&                WRITE(STD_OUT,'(7x,i2,a,a1,1x,i1,2x,i2,2x,a,f7.4)') &
&                norb," - ",".",ll,kk," : ",dataset%basis_func_rc(norb)
             END IF
           END DO
         END IF
       END DO
     END DO
   END IF
      
 ELSE ! Other projectors
   IF (ALLOCATED(dataset%basis_func_rc)) DEALLOCATE(dataset%basis_func_rc)
   ALLOCATE(dataset%basis_func_rc(0))
 END IF


!------------------------------------------------------------------
!End reading of basis data
 ENDIF

!Final message
 IF (read_global_data_.OR.read_elec_data_.OR.read_coreval_data_.OR.read_basis_data_) THEN
   IF (has_to_ask) THEN
     WRITE(STD_OUT,'(1x,a)') "===== END READING OF INPUT DATA ====="
   ELSE
     WRITE(STD_OUT,'(3x,a)') "===== END READING OF INPUT FILE ====="
   END IF
 END IF
 WRITE(STD_OUT,'(2/)')


!------------------------------------------------------------------
!Close files
 IF (PRESENT(inputfile)) THEN
   CLOSE(ifunit)
 END IF
 IF (has_to_echo) THEN
   CLOSE(ecunit)
 END IF
 
 END SUBROUTINE input_dataset_read


!!=================================================================
!! NAME
!!  input_dataset_free
!!
!! FUNCTION
!!  Free and destroy a input_dataset datastructure
!!
!! SIDE EFFECT
!!  [input_dt]= datastructure to be destroyed (optional)
!!              If omitted, then the global public `input_dataset`
!!              is used.
!!
!!=================================================================
 SUBROUTINE input_dataset_free(input_dt)

!---- Arguments
 TYPE(input_dataset_t),INTENT(inout),OPTIONAL,TARGET :: input_dt

!---- Local variables
 TYPE(input_dataset_t),POINTER :: dataset

!------------------------------------------------------------------

!Select datastruture to be read
 IF (PRESENT(input_dt)) THEN
   dataset => input_dt
 ELSE
   dataset => input_dataset
 ENDIF

!Integers
 dataset%atomic_charge            =-1
 dataset%finitenucleusmodel       =-1
 dataset%gridpoints               = 0
 dataset%nlogderiv                = 0
 dataset%Fixed_zero_index         = 0
 dataset%np(1:5)                  = 0
 dataset%norbit                   = 0
 dataset%norbit_mod               = 0
 dataset%norbit_val               = 0
 dataset%lmax                     =-1
 dataset%nbasis                   = 0
 dataset%nbasis_add               = 0
 dataset%projector_type           = UNKNOWN_TYPE
 dataset%pseudo_type              = UNKNOWN_TYPE
 dataset%ortho_type               = UNKNOWN_TYPE
 dataset%vloc_type                = UNKNOWN_TYPE
 dataset%vloc_l                   = 0
 dataset%pseudo_polynom2_pdeg     = 0
 dataset%shapefunc_type           = UNKNOWN_TYPE
 dataset%abinit_log_meshsz        = 0
 dataset%xml_spl_meshsz           = 0
 dataset%xml_lda12_orb_l          = -1
 dataset%xml_lda12_orb_n          = -1
 dataset%vloc_kerker_power(1:4)   = 0

!Logicals
 dataset%scalarrelativistic       = .false.
 dataset%diracrelativistic        = .false.
 dataset%finitenucleus            = .false.
 dataset%HFpostprocess            = .false.
 dataset%localizedcoreexchange    = .false.
 dataset%BDsolve                  = .false.
 dataset%Fixed_zero               = .false.
 dataset%needvtau                 = .false.
 dataset%shapetcore               = .false.
 dataset%abinit_usexcnhat         = .false.
 dataset%abinit_prtcorewf         = .false.
 dataset%abinit_uselog            = .false.
 dataset%abinit_userso            = .false.
 dataset%xml_usexcnhat            = .false.
 dataset%xml_prtcorewf            = .false.
 dataset%xml_usespl               = .false.
 dataset%xml_userso               = .false.
 dataset%xml_uselda12             = .false.

!Reals
 dataset%gridrange                = 0.d0
 dataset%gridmatch                = 0.d0
 dataset%minlogderiv              = 0.d0
 dataset%maxlogderiv              = 0.d0
 dataset%rc                       = 0.d0
 dataset%rc_shap                  = 0.d0
 dataset%rc_vloc                  = 0.d0
 dataset%rc_core                  = 0.d0
 dataset%pseudo_polynom2_qcut     = 0.d0
 dataset%shapefunc_gaussian_param = 0.d0
 dataset%hf_coretol               = 0.d0
 dataset%vloc_ene                 = 0.d0
 dataset%vloc_setvloc_coef        = 0.d0
 dataset%vloc_setvloc_rad         = 0.d0
 dataset%abinit_log_step          = 0.d0
 dataset%abinit_rso_ecut          = 0.d0
 dataset%abinit_rso_gfact         = 0.d0
 dataset%abinit_rso_werror        = 0.d0
 dataset%xml_rso_ecut             = 0.d0
 dataset%xml_rso_gfact            = 0.d0
 dataset%xml_rso_werror           = 0.d0
 dataset%xml_lda12_ion            = 0.d0
 dataset%xml_lda12_rcut           = 0.d0
 dataset%upf_grid_xmin            = 0.d0
 dataset%upf_grid_zmesh           = 0.d0
 dataset%upf_grid_dx              = 0.d0
 dataset%upf_grid_range           = 0.d0

!Character strings
 dataset%atomic_symbol            = ''
 dataset%gridkey                  = ''
 dataset%exctype                  = ''
 dataset%abinit_author            = ''
 dataset%xml_author               = ''
 dataset%xml_comment              = ''
 dataset%xml_lda12_logfile        = ''

!Integer allocatable arrays
 IF (ALLOCATED(dataset%orbit_mod_l)) DEALLOCATE(dataset%orbit_mod_l)
 IF (ALLOCATED(dataset%orbit_mod_n)) DEALLOCATE(dataset%orbit_mod_n)
 IF (ALLOCATED(dataset%orbit_mod_k)) DEALLOCATE(dataset%orbit_mod_k)
 IF (ALLOCATED(dataset%orbit_val_l)) DEALLOCATE(dataset%orbit_val_l)
 IF (ALLOCATED(dataset%orbit_val_n)) DEALLOCATE(dataset%orbit_val_n)
 IF (ALLOCATED(dataset%orbit_val_k)) DEALLOCATE(dataset%orbit_val_k)
 IF (ALLOCATED(dataset%basis_add_l)) DEALLOCATE(dataset%basis_add_l)
 IF (ALLOCATED(dataset%basis_add_k)) DEALLOCATE(dataset%basis_add_k)

!Logical allocatable arrays
 IF (ALLOCATED(dataset%orbit_iscore)) DEALLOCATE(dataset%orbit_iscore)

!Real allocatable arrays
 IF (ALLOCATED(dataset%orbit_mod_occ))    DEALLOCATE(dataset%orbit_mod_occ)
 IF (ALLOCATED(dataset%basis_add_energy)) DEALLOCATE(dataset%basis_add_energy)
 IF (ALLOCATED(dataset%basis_func_rc))    DEALLOCATE(dataset%basis_func_rc)

 END SUBROUTINE input_dataset_free


!!=================================================================
!! NAME
!!  input_dataset_copy
!!
!! FUNCTION
!!  Copy a input_dataset datastructure into another
!!
!! INPUTS
!!  input_dt= datastructure to copy
!!
!! OUTPUT
!!  output_dt= copied datastructure
!!
!!=================================================================
 SUBROUTINE input_dataset_copy(input_dt,output_dt)

!---- Arguments
 TYPE(input_dataset_t),INTENT(in) :: input_dt
 TYPE(input_dataset_t),INTENT(inout) :: output_dt

!------------------------------------------------------------------

!Integers
 output_dt%atomic_charge           =input_dt%atomic_charge
 output_dt%finitenucleusmodel      =input_dt%finitenucleusmodel
 output_dt%gridpoints              =input_dt%gridpoints
 output_dt%nlogderiv               =input_dt%nlogderiv
 output_dt%Fixed_zero_index        =input_dt%Fixed_zero_index
 output_dt%np(1:5)                 =input_dt%np(1:5)
 output_dt%norbit                  =input_dt%norbit
 output_dt%norbit_mod              =input_dt%norbit_mod
 output_dt%norbit_val              =input_dt%norbit_val
 output_dt%lmax                    =input_dt%lmax
 output_dt%nbasis                  =input_dt%nbasis
 output_dt%nbasis_add              =input_dt%nbasis_add
 output_dt%projector_type          =input_dt%projector_type
 output_dt%pseudo_type             =input_dt%pseudo_type
 output_dt%ortho_type              =input_dt%ortho_type
 output_dt%pseudo_polynom2_pdeg    =input_dt%pseudo_polynom2_pdeg
 output_dt%shapefunc_type          =input_dt%shapefunc_type
 output_dt%vloc_type               =input_dt%vloc_type
 output_dt%vloc_l                  =input_dt%vloc_l
 output_dt%abinit_log_meshsz       =input_dt%abinit_log_meshsz
 output_dt%xml_spl_meshsz          =input_dt%xml_spl_meshsz
 output_dt%xml_lda12_orb_l         =input_dt%xml_lda12_orb_l
 output_dt%xml_lda12_orb_n         =input_dt%xml_lda12_orb_n
 output_dt%vloc_kerker_power(1:4)  =input_dt%vloc_kerker_power(1:4)

!Logicals
 output_dt%scalarrelativistic      =input_dt%scalarrelativistic
 output_dt%diracrelativistic       =input_dt%diracrelativistic
 output_dt%finitenucleus           =input_dt%finitenucleus
 output_dt%HFpostprocess           =input_dt%HFpostprocess
 output_dt%localizedcoreexchange   =input_dt%localizedcoreexchange
 output_dt%BDsolve                 =input_dt%BDsolve
 output_dt%Fixed_zero              =input_dt%Fixed_zero
 output_dt%needvtau                =input_dt%needvtau
 output_dt%shapetcore              =input_dt%shapetcore
 output_dt%abinit_usexcnhat        =input_dt%abinit_usexcnhat
 output_dt%abinit_prtcorewf        =input_dt%abinit_prtcorewf
 output_dt%abinit_uselog           =input_dt%abinit_uselog
 output_dt%abinit_userso           =input_dt%abinit_userso
 output_dt%xml_usexcnhat           =input_dt%xml_usexcnhat
 output_dt%xml_prtcorewf           =input_dt%xml_prtcorewf
 output_dt%xml_usespl              =input_dt%xml_usespl
 output_dt%xml_userso              =input_dt%xml_userso
 output_dt%xml_uselda12            =input_dt%xml_uselda12

!Reals
 output_dt%gridrange               =input_dt%gridrange
 output_dt%gridmatch               =input_dt%gridmatch
 output_dt%minlogderiv             =input_dt%minlogderiv
 output_dt%maxlogderiv             =input_dt%maxlogderiv
 output_dt%rc                      =input_dt%rc
 output_dt%rc_shap                 =input_dt%rc_shap
 output_dt%rc_vloc                 =input_dt%rc_vloc
 output_dt%rc_core                 =input_dt%rc_core
 output_dt%pseudo_polynom2_qcut    =input_dt%pseudo_polynom2_qcut
 output_dt%shapefunc_gaussian_param=input_dt%shapefunc_gaussian_param
 output_dt%hf_coretol              =input_dt%hf_coretol
 output_dt%vloc_ene                =input_dt%vloc_ene
 output_dt%vloc_setvloc_coef       =input_dt%vloc_setvloc_coef
 output_dt%vloc_setvloc_rad        =input_dt%vloc_setvloc_rad
 output_dt%abinit_log_step         =input_dt%abinit_log_step
 output_dt%abinit_rso_ecut         =input_dt%abinit_rso_ecut
 output_dt%abinit_rso_gfact        =input_dt%abinit_rso_gfact
 output_dt%abinit_rso_werror       =input_dt%abinit_rso_werror
 output_dt%xml_rso_ecut            =input_dt%xml_rso_ecut
 output_dt%xml_rso_gfact           =input_dt%xml_rso_gfact
 output_dt%xml_rso_werror          =input_dt%xml_rso_werror
 output_dt%xml_lda12_ion           =input_dt%xml_lda12_ion
 output_dt%xml_lda12_rcut          =input_dt%xml_lda12_rcut
 output_dt%upf_grid_xmin           =input_dt%upf_grid_xmin
 output_dt%upf_grid_zmesh          =input_dt%upf_grid_zmesh
 output_dt%upf_grid_dx             =input_dt%upf_grid_dx
 output_dt%upf_grid_range          =input_dt%upf_grid_range

!Character strings
 output_dt%atomic_symbol           =TRIM(input_dt%atomic_symbol)
 output_dt%gridkey                 =TRIM(input_dt%gridkey)
 output_dt%exctype                 =TRIM(input_dt%exctype)
 output_dt%abinit_author           =TRIM(input_dt%abinit_author)
 output_dt%xml_author              =TRIM(input_dt%xml_author)
 output_dt%xml_comment             =TRIM(input_dt%xml_comment)
 output_dt%xml_lda12_logfile       =TRIM(input_dt%xml_lda12_logfile)

!Integer allocatable arrays
 IF (ALLOCATED(output_dt%orbit_mod_l)) DEALLOCATE(output_dt%orbit_mod_l)
 IF (ALLOCATED(output_dt%orbit_mod_n)) DEALLOCATE(output_dt%orbit_mod_n)
 IF (ALLOCATED(output_dt%orbit_mod_k)) DEALLOCATE(output_dt%orbit_mod_k)
 IF (ALLOCATED(output_dt%orbit_val_l)) DEALLOCATE(output_dt%orbit_val_l)
 IF (ALLOCATED(output_dt%orbit_val_n)) DEALLOCATE(output_dt%orbit_val_n)
 IF (ALLOCATED(output_dt%orbit_val_k)) DEALLOCATE(output_dt%orbit_val_k)
 IF (ALLOCATED(output_dt%basis_add_l)) DEALLOCATE(output_dt%basis_add_l)
 IF (ALLOCATED(output_dt%basis_add_k)) DEALLOCATE(output_dt%basis_add_k)
 ALLOCATE(output_dt%orbit_mod_l(input_dt%norbit_mod))
 ALLOCATE(output_dt%orbit_mod_n(input_dt%norbit_mod))
 ALLOCATE(output_dt%orbit_mod_k(input_dt%norbit_mod))
 ALLOCATE(output_dt%orbit_val_l(input_dt%norbit_val))
 ALLOCATE(output_dt%orbit_val_n(input_dt%norbit_val))
 ALLOCATE(output_dt%orbit_val_k(input_dt%norbit_val))
 ALLOCATE(output_dt%basis_add_l(input_dt%nbasis_add))
 ALLOCATE(output_dt%basis_add_k(input_dt%nbasis_add))
 output_dt%orbit_mod_l(1:input_dt%norbit_mod)=input_dt%orbit_mod_l(1:input_dt%norbit_mod)
 output_dt%orbit_mod_n(1:input_dt%norbit_mod)=input_dt%orbit_mod_n(1:input_dt%norbit_mod)
 output_dt%orbit_mod_k(1:input_dt%norbit_mod)=input_dt%orbit_mod_k(1:input_dt%norbit_mod)
 output_dt%orbit_val_l(1:input_dt%norbit_val)=input_dt%orbit_val_l(1:input_dt%norbit_val)
 output_dt%orbit_val_n(1:input_dt%norbit_val)=input_dt%orbit_val_n(1:input_dt%norbit_val)
 output_dt%orbit_val_k(1:input_dt%norbit_val)=input_dt%orbit_val_k(1:input_dt%norbit_val)
 output_dt%basis_add_l(1:input_dt%nbasis_add)=input_dt%basis_add_l(1:input_dt%nbasis_add)
 output_dt%basis_add_k(1:input_dt%nbasis_add)=input_dt%basis_add_k(1:input_dt%nbasis_add)

!Logical allocatable arrays
 IF (ALLOCATED(output_dt%orbit_iscore)) DEALLOCATE(output_dt%orbit_iscore)
 ALLOCATE(output_dt%orbit_iscore(input_dt%norbit))
 output_dt%orbit_iscore(1:input_dt%norbit)=input_dt%orbit_iscore(1:input_dt%norbit)

!Real allocatable arrays
 IF (ALLOCATED(output_dt%orbit_mod_occ))    DEALLOCATE(output_dt%orbit_mod_occ)
 IF (ALLOCATED(output_dt%basis_add_energy)) DEALLOCATE(output_dt%basis_add_energy)
 IF (ALLOCATED(output_dt%basis_func_rc))    DEALLOCATE(output_dt%basis_func_rc)
 ALLOCATE(output_dt%orbit_mod_occ(input_dt%norbit_mod))
 ALLOCATE(output_dt%basis_add_energy(input_dt%nbasis_add))
 ALLOCATE(output_dt%basis_func_rc(input_dt%nbasis))
 output_dt%orbit_mod_occ(1:input_dt%norbit_mod)=input_dt%orbit_mod_occ(1:input_dt%norbit_mod)
 output_dt%basis_add_energy(1:input_dt%nbasis_add)=input_dt%basis_add_energy(1:input_dt%nbasis_add)
 output_dt%basis_func_rc(1:input_dt%nbasis)=input_dt%basis_func_rc(1:input_dt%nbasis)

 END SUBROUTINE input_dataset_copy


!!=================================================================
!! NAME
!!  input_dataset_read_occ
!!
!! FUNCTION
!!  Initialize only modified occupations in a data-structure by reading them
!!  from a file. If file is omitted, then read from standard input.
!!  Note: file is supposed to be opened.
!!
!! INPUTS
!!  np(5)= electronic configuration: number of s,p,d,f,g shells
!!  diracrel= TRUE if we perform a Dirac relativistic calculation
!!  [inputfile_unit]= logical unit of input file to be read (optional)
!!  [echofile_unit]= logical unit of a file to echo input file content (optional)
!!
!! OUTPUT
!!  norbit_mod= number of orbital with modified occupations
!!  orbit_mod_l(norbit_mod)= l quantum numbers for the modifed orbitals
!!  orbit_mod_n(norbit_mod)= n quantum numbers for the modifed orbitals
!!  orbit_mod_k(norbit_mod)= kappa quantum numbers for the modifed orbitals
!!  orbit_mod_occ(norbit_mod)= occupations for the modifed orbitals
!!  [input_dt]= data-structure containing the complete input file (optional).
!!              If omitted, then the global public `input_dataset`
!!              is used.
!!
!!=================================================================
 SUBROUTINE input_dataset_read_occ(norbit_mod,orbit_mod_l,orbit_mod_n,orbit_mod_k,&
 &                          orbit_mod_occ,np,diracrel,inputfile_unit,echofile_unit)

!---- Arguments
 INTEGER,INTENT(IN) :: np(5)
 LOGICAL,INTENT(IN) :: diracrel
 INTEGER,INTENT(IN),OPTIONAL :: inputfile_unit,echofile_unit
 INTEGER,INTENT(OUT) :: norbit_mod
 INTEGER,ALLOCATABLE,INTENT(INOUT) :: orbit_mod_l(:),orbit_mod_n(:),orbit_mod_k(:)
 REAL(8),ALLOCATABLE,INTENT(INOUT) :: orbit_mod_occ(:)

!---- Local variables
 INTEGER,PARAMETER :: nkappa(5)=(/1,2,2,2,2/)
 INTEGER :: input_unit,echo_unit,ll,nn,io,ii,kk,ik
 LOGICAL :: has_to_echo
 CHARACTER(200) :: inputline
 REAL(8) :: xocc
 INTEGER :: tmp_n(norbit_max),tmp_l(norbit_max),tmp_k(norbit_max)
 REAL(8) :: tmp_occ(norbit_max),basis_add_energy(nbasis_add_max)

!------------------------------------------------------------------

!Select input file logical unit
 input_unit=STD_IN;IF (PRESENT(inputfile_unit)) input_unit=inputfile_unit

!Do we echo input file content?
 has_to_echo=PRESENT(echofile_unit)
 echo_unit=-1;IF (has_to_echo) echo_unit=echofile_unit

 IF (has_to_ask) THEN
   IF (.NOT.diracrel) THEN
     WRITE(STD_OUT,*)'Enter np l occ for all occupations for all revisions'
     WRITE(STD_OUT,*)'Enter 0 0 0 to end'
   ELSE
     WRITE(STD_OUT,*)'Enter np l kappa occ for all occupations for all revisions'
     WRITE(STD_OUT,*)'Enter 0 0 0 0 to end'
   END IF
 END IF

 norbit_mod=0 ; kk=0
 DO
   READ(input_unit,'(a)') inputline
   IF (has_to_echo) WRITE(echo_unit,'(a)') TRIM(inputline)
   if (.not.diracrel) READ(inputline,*) nn,ll,xocc
   if (     diracrel) READ(inputline,*) nn,ll,kk,xocc
   IF (nn<=0) EXIT
   IF (xocc<0.d0.OR.&
&    ((.NOT.diracrel).AND.(xocc>2.d0*(2*ll+1))).OR.&
&    ((     diracrel).AND.(xocc>2.d0*ABS(kk)))) THEN
     WRITE(STD_OUT,*) 'input_dataset: error in occupations -- ip,l,kap,xocc: ',nn,ll,kk,xocc,'!'
     STOP
   END IF
   norbit_mod=norbit_mod+1
   if (norbit_mod>norbit_max) stop 'input_dataset_occ: error -- too many occupation lines!'
   tmp_l(norbit_mod)=ll
   tmp_n(norbit_mod)=nn
   tmp_k(norbit_mod)=kk
   tmp_occ(norbit_mod)=xocc
 END DO

 IF (ALLOCATED(orbit_mod_l)) DEALLOCATE(orbit_mod_l)
 IF (ALLOCATED(orbit_mod_n)) DEALLOCATE(orbit_mod_n)
 IF (ALLOCATED(orbit_mod_k)) DEALLOCATE(orbit_mod_k)
 IF (ALLOCATED(orbit_mod_occ)) DEALLOCATE(orbit_mod_occ)
 ALLOCATE(orbit_mod_l(norbit_mod))
 ALLOCATE(orbit_mod_n(norbit_mod))
 ALLOCATE(orbit_mod_k(norbit_mod))
 ALLOCATE(orbit_mod_occ(norbit_mod))
 orbit_mod_l(1:norbit_mod)=tmp_l(1:norbit_mod)
 orbit_mod_n(1:norbit_mod)=tmp_n(1:norbit_mod)
 orbit_mod_k(1:norbit_mod)=tmp_k(1:norbit_mod)
 orbit_mod_occ(1:norbit_mod)=tmp_occ(1:norbit_mod)

!Print read data
IF (has_to_print) THEN
   WRITE(STD_OUT,'(3x,a)') "Occupations of orbitals:"
   IF (.NOT.diracrel) THEN
     WRITE(STD_OUT,'(7x,a)') "n l :  occ"
     DO ll=0,4
       IF (np(ll+1)>0) THEN
         DO ii=1+ll,np(ll+1)
           nn=-1
           DO io=1,norbit_mod
             IF (orbit_mod_l(io)==ll.AND.orbit_mod_n(io)==ii) THEN
               nn=io ; EXIT
             END IF
           END DO
           IF (nn<=0) WRITE(STD_OUT,'(7x,i1,1x,i1,a,f4.1)') ii,ll," : ",real(2*(2*ll+1))
           IF (nn >0) WRITE(STD_OUT,'(7x,i1,1x,i1,a,f4.1)') ii,ll," : ",orbit_mod_occ(nn)
         END DO
       END IF
     END DO
   ELSE
     WRITE(STD_OUT,'(7x,a)') "n l kappa :  occ"
     DO ll=0,4
       IF (np(ll+1)>0) THEN
         DO ik=1,nkappa(ll+1)
           kk=MERGE(ll,-(ll+1),ik==1);IF (ll==0) kk=-1
           DO ii=1+ll,np(ll+1)
             nn=-1
             DO io=1,norbit_mod
               IF (orbit_mod_l(io)==ll.and.orbit_mod_n(io)==ii.AND.orbit_mod_k(io)==kk) THEN
                 nn=io ; EXIT
               END IF
             END DO
             IF (nn<=0) WRITE(STD_OUT,'(7x,i1,1x,i1,2x,i2,2x,a,f4.1)') ii,ll,kk," : ",real(2*abs(kk))
             IF (nn >0) WRITE(STD_OUT,'(7x,i1,1x,i1,2x,i2,2x,a,f4.1)') ii,ll,kk," : ",orbit_mod_occ(nn)
           END DO
         END DO
       END IF
     END DO
   END IF
END IF

 END SUBROUTINE input_dataset_read_occ


!!=================================================================
!! NAME
!! input_dataset_read_abinit
!!
!! FUNCTION
!! Read the input file (or standard input) in order to get ABINIT options.
!! Put them in a input_dataset datastructure.
!! If file is omitted, then read from standard input.
!! File is supposed to be opened.
!!
!! INPUTS
!!  [inputfile_unit]= logical unit of input file to be read (optional)
!!  [echofile_unit]= logical unit of a file to echo input file content (optional)
!!
!! OUTPUT
!!  [input_dt]= data-structure containing the complete input file (optional).
!!              If omitted, then the global public `input_dataset`
!!    %abinit_usexcnhat=TRUE if nhat density (compensation) is included in XC potential
!!    %abinit_prtcorewf= TRUE is printing of core WF is required
!!    %abinit_author= author to be printed in the PAW dataset file
!!    %abinit_uselog=TRUE if data are transfered on a log. grid before being written
!!    %abinit_log_meshsz=mesh size for the logarithmic grid
!!    %abinit_log_step=logarithmic step for the logarithmic grid
!!    %abinit_userso=TRUE if REAL Space Optimization is required
!!    %abinit_rso_ecut=Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
!!    %abinit_rso_gfact=Real Space Optimization parameter: Gamma/Gmax ratio
!!    %abinit_rso_werror=Real Space Optimization parameter: max. error W_l allowed
!!  [abinit_string]= character string containing the ABINIT options line from input file
!!
!!=================================================================
 SUBROUTINE input_dataset_read_abinit(input_dt,inputfile_unit,echofile_unit,abinit_string)

!---- Arguments
 INTEGER,INTENT(IN),OPTIONAL :: inputfile_unit,echofile_unit
 CHARACTER(*),INTENT(OUT),OPTIONAL :: abinit_string
 TYPE(input_dataset_t),INTENT(INOUT),OPTIONAL,TARGET :: input_dt

!---- Local variables
 INTEGER :: input_unit,echo_unit
 INTEGER :: i_author,i_usexcnhat,i_prtcorewf,i_logspline,i_rsoptim,iend
 LOGICAL :: has_to_echo
 CHARACTER(200) :: inputline,inputstring,inputword
 TYPE(input_dataset_t),POINTER :: dataset

!------------------------------------------------------------------

!Select datastruture to be read
 IF (PRESENT(input_dt)) THEN
   dataset => input_dt
 ELSE
   dataset => input_dataset
 ENDIF

!Select input file logical unit
 input_unit=STD_IN;IF (PRESENT(inputfile_unit)) input_unit=inputfile_unit

!Do we echo input file content?
 has_to_echo=PRESENT(echofile_unit)
 echo_unit=-1;IF (has_to_echo) echo_unit=echofile_unit

!Reading of keywords
 READ(input_unit,'(a)') inputline
 IF (has_to_echo) WRITE(echo_unit,'(a)') TRIM(inputline)
 IF (PRESENT(abinit_string)) abinit_string=TRIM(inputline)
 CALL Uppercase(inputline)
 i_usexcnhat=INDEX(inputline,'USEXCNHAT')
 i_prtcorewf=INDEX(inputline,'PRTCOREWF')
 i_logspline=INDEX(inputline,'LOGSPLINE')
 i_rsoptim  =INDEX(inputline,'RSOPTIM')
 i_author   =INDEX(inputline,'AUTHOR')

!Option for core WF printing
 dataset%abinit_prtcorewf=MERGE(.true.,PRTCOREWF_DEF,i_prtcorewf>0)

!Option for use of NHAT in XC
 dataset%abinit_usexcnhat=MERGE(.true.,USEXCNHAT_DEF,i_usexcnhat>0)

!To be activated later:
!dataset%abinit_vbare=merge(.true.,vbare_def,i_vbare>0)

!Options related to the use of REAL SPACE OPTIMIZATION
 dataset%abinit_userso=MERGE(.true.,USERSO_DEF,i_rsoptim>0)
 dataset%abinit_rso_ecut=RSO_ECUT_DEF
 dataset%abinit_rso_gfact=RSO_GFACT_DEF
 dataset%abinit_rso_werror=RSO_WERROR_DEF
 IF (dataset%abinit_userso) THEN
   iend=200
   IF (i_usexcnhat>i_rsoptim.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_rsoptim.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_logspline>i_rsoptim.AND.i_logspline-1<iend) iend=i_logspline-1
   IF (i_author   >i_rsoptim.AND.i_author   -1<iend) iend=i_author   -1
   inputstring="";IF (iend>i_rsoptim+7) inputstring=TRIM(inputline(i_rsoptim+7:iend))
   IF (inputstring/="") THEN
     CALL extractword(1,inputstring,inputword);inputword=TRIM(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) dataset%abinit_rso_ecut
       CALL extractword(2,inputstring,inputword);inputword=TRIM(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) dataset%abinit_rso_gfact
         CALL extractword(3,inputstring,inputword);inputword=TRIM(inputword)
         IF (inputword/="") READ(inputword,*) dataset%abinit_rso_werror
       END IF
     END IF
   END IF
 END IF

!Options for related to the transfer to a reduced logarihmic grid
 dataset%abinit_uselog=MERGE(.true.,USELOG_DEF,i_logspline>0)
 dataset%abinit_log_meshsz=LOGGRD_SIZE_DEF
 dataset%abinit_log_step=LOGGRD_STEP_DEF
 IF (dataset%abinit_uselog) THEN
   iend=200
   IF (i_usexcnhat>i_logspline.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_logspline.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_rsoptim  >i_logspline.AND.i_rsoptim  -1<iend) iend=i_rsoptim  -1
   IF (i_author   >i_logspline.AND.i_author   -1<iend) iend=i_author   -1
   inputstring="";IF (iend>i_logspline+9) inputstring=TRIM(inputline(i_logspline+9:iend))
   IF (inputstring/="") THEN
     CALL extractword(1,inputstring,inputword);inputword=TRIM(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) dataset%abinit_log_meshsz
       CALL extractword(2,inputstring,inputword);inputword=TRIM(inputword)
       IF (inputword/="") READ(inputword,*) dataset%abinit_log_step
     END IF
   END IF
 END IF

!Author to be mentioned in the ABINIT file
 dataset%abinit_author=TRIM(AUTHOR_DEF)
 IF (i_author>0) then
   iend=200
   IF (i_usexcnhat>i_author.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_author.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_logspline>i_author.AND.i_logspline-1<iend) iend=i_logspline-1
   IF (i_rsoptim  >i_author.AND.i_rsoptim  -1<iend) iend=i_rsoptim  -1
   inputstring=TRIM(inputline(i_author+6:iend))
   READ(unit=inputstring,fmt=*) dataset%abinit_author
   dataset%abinit_author=TRIM(dataset%abinit_author)
 END IF

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(/,2x,a)')'Options for ABINIT file:'
   WRITE(STD_OUT,'(2x,2a)') 'ABINIT option : output of core wave funcs =',MERGE("YES"," NO",dataset%abinit_prtcorewf)
   WRITE(STD_OUT,'(2x,2a)') 'ABINIT option : transfer to a log grid    =',MERGE("YES"," NO",dataset%abinit_uselog)
   IF (dataset%abinit_uselog) THEN
     WRITE(STD_OUT,'(5x,a,i5)') '- mesh size of the log grid =',dataset%abinit_log_meshsz
     WRITE(STD_OUT,'(5x,a,g8.2)') '- step of the log grid      = ',dataset%abinit_log_step
   END IF
   WRITE(STD_OUT,'(2x,2a)') 'ABINIT option : Real Space Optimization   =',MERGE("YES"," NO",dataset%abinit_userso)
   IF (dataset%abinit_userso) THEN
     WRITE(STD_OUT,'(5x,a,f7.2)') '- RSO plane wave cutoff =',dataset%abinit_rso_ecut
     WRITE(STD_OUT,'(5x,a,1x,f6.2)') '- RSO Gamma/Gmax ratio  =',dataset%abinit_rso_gfact
     WRITE(STD_OUT,'(5x,a,g11.3)') '- RSO maximum error     = ',dataset%abinit_rso_werror
   END IF
   WRITE(STD_OUT,'(2x,2a)') 'ABINIT option : use of compensation in XC =',MERGE("YES"," NO",dataset%abinit_usexcnhat)
   IF (dataset%abinit_usexcnhat) THEN
     WRITE(STD_OUT,'(5x,a)') '- Kresse local ionic potential output in XML file'
   ELSE
     WRITE(STD_OUT,'(5x,a)') '- Blochl local ionic potential output in XML file'
   END IF
 END IF

 END SUBROUTINE input_dataset_read_abinit


!!=================================================================
!! NAME
!! input_dataset_read_xml
!!
!! FUNCTION
!! Read the input file (or standard input) in order to get UPF options.
!! Put them in a input_dataset datastructure.
!! If file is omitted, then read from standard input.
!! File is supposed to be opened.
!!
!! INPUTS
!!  [inputfile_unit]= logical unit of input file to be read (optional)
!!  [echofile_unit]= logical unit of a file to echo input file content (optional)
!!
!! OUTPUT
!!  [input_dt]= data-structure containing the complete input file (optional).
!!              If omitted, then the global public `input_dataset`
!!    %_usexcnhat=TRUE if nhat density (compensation) is included in XC potential
!!    %xml_prtcorewf= TRUE is printing of core WF is required
!!    %xml_author= author to be printed in the PAW dataset file
!!    %xml_comment= comment line to be printed in the heder of the PAW dataset file
!!    %xml_uselog=TRUE if data are transfered on a log. grid before being written
!!    %xml_log_meshsz=mesh size for the logarithmic grid
!!    %xml_log_step=logarithmic step for the logarithmic grid
!!    %xml_userso=TRUE if REAL Space Optimization is required
!!    %xml_rso_ecut=Real Space Optimization parameter: plane wave cutoff = 1/2 Gmax**2
!!    %xml_rso_gfact=Real Space Optimization parameter: Gamma/Gmax ratio
!!    %xml_rso_werror=Real Space Optimization parameter: max. error W_l allowed
!!    %xml_uselda12= TRUE if LDA-1/2 potential calculation is required
!!    %xml_lda12_orb_l=LDA-1/2 parameter: l quantum number of the orbital to be ionized
!!    %xml_lda12_orb_n=LDA-1/2 parameter: n quantum number of the orbital to be ionized
!!    %xml_lda12_ion=LDA-1/2 parameter: amount of charge to be removed from the ionized orbital
!!    %xml_lda12_rcut=LDA-1/2 parameter: cut-off radius (in bohr)
!!    %xml_lda12_logfile=LDA-1/2 parameter: name of the log file
!!  [xml_string]= character string containing the ABINIT options line from input file
!!
!!=================================================================
 SUBROUTINE input_dataset_read_xml(input_dt,inputfile_unit,echofile_unit,xml_string)

!---- Arguments
 INTEGER,INTENT(IN),OPTIONAL :: inputfile_unit,echofile_unit
 CHARACTER(*),INTENT(OUT),OPTIONAL :: xml_string
 TYPE(input_dataset_t),INTENT(INOUT),OPTIONAL,TARGET :: input_dt

!---- Local variables
 INTEGER :: input_unit,echo_unit
 INTEGER :: i_author,i_comment,i_usexcnhat,i_prtcorewf,i_logspline,i_rsoptim,i_lda12,iend
 LOGICAL :: has_to_echo
 CHARACTER(200) :: inputline,inputstring,inputword
 TYPE(input_dataset_t),POINTER :: dataset

!------------------------------------------------------------------

!Select datastruture to be read
 IF (PRESENT(input_dt)) THEN
   dataset => input_dt
 ELSE
   dataset => input_dataset
 ENDIF

!Select input file logical unit
 input_unit=STD_IN;IF (PRESENT(inputfile_unit)) input_unit=inputfile_unit

!Do we echo input file content?
 has_to_echo=PRESENT(echofile_unit)
 echo_unit=-1;IF (has_to_echo) echo_unit=echofile_unit

!Reading of keywords
 READ(input_unit,'(a)') inputline
 IF (has_to_echo) WRITE(echo_unit,'(a)') TRIM(inputline)
 IF (PRESENT(xml_string)) xml_string=TRIM(inputline)
 CALL Uppercase(inputline)
 i_usexcnhat=INDEX(inputline,'USEXCNHAT')
 i_prtcorewf=INDEX(inputline,'PRTCOREWF')
 i_logspline=INDEX(inputline,'WITHSPLGRID')
 i_rsoptim  =INDEX(inputline,'RSOPTIM')
 i_lda12    =INDEX(inputline,'LDA12')
 i_author   =INDEX(inputline,'AUTHOR')
 i_comment  =INDEX(inputline,'COMMENT')

!Option for core WF printing
 dataset%xml_prtcorewf=MERGE(.true.,PRTCOREWF_DEF,i_prtcorewf>0)

!Option for use of NHAT in XC
 dataset%xml_usexcnhat=MERGE(.true.,USEXCNHAT_DEF,i_usexcnhat>0)

!To be activated later:
!dataset%xml_vbare=merge(.true.,vbare_def,i_vbare>0)

!Options related to the use of REAL SPACE OPTIMIZATION
 dataset%xml_userso=MERGE(.true.,USERSO_DEF,i_rsoptim>0)
 dataset%xml_rso_ecut=RSO_ECUT_DEF
 dataset%xml_rso_gfact=RSO_GFACT_DEF
 dataset%xml_rso_werror=RSO_WERROR_DEF
 IF (dataset%xml_userso) THEN
   iend=200
   IF (i_usexcnhat>i_rsoptim.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_rsoptim.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_logspline>i_rsoptim.AND.i_logspline-1<iend) iend=i_logspline-1
   IF (i_lda12    >i_rsoptim.AND.i_lda12    -1<iend) iend=i_lda12    -1
   IF (i_author   >i_rsoptim.AND.i_author   -1<iend) iend=i_author   -1
   IF (i_comment  >i_rsoptim.AND.i_comment  -1<iend) iend=i_comment  -1
   inputstring="";IF (iend>i_rsoptim+7) inputstring=TRIM(inputline(i_rsoptim+7:iend))
   IF (inputstring/="") THEN
     CALL extractword(1,inputstring,inputword);inputword=TRIM(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) dataset%xml_rso_ecut
       CALL extractword(2,inputstring,inputword);inputword=TRIM(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) dataset%xml_rso_gfact
         CALL extractword(3,inputstring,inputword);inputword=TRIM(inputword)
         IF (inputword/="") READ(inputword,*) dataset%xml_rso_werror
       END IF
     END IF
   END IF
 END IF

!Options related to the use of LDA-1/2 technique
 dataset%xml_uselda12=MERGE(.true.,USELDA12_DEF,i_lda12>0)
 dataset%xml_lda12_orb_l=LDA12_ORB_L_DEF
 dataset%xml_lda12_orb_n=LDA12_ORB_N_DEF
 dataset%xml_lda12_ion=LDA12_ION_DEF
 dataset%xml_lda12_rcut=LDA12_RCUT_DEF
 dataset%xml_lda12_logfile=TRIM(LDA12_LOGFILE)
 IF (dataset%xml_uselda12) THEN
   iend=200
   IF (i_usexcnhat>i_lda12.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_lda12.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_logspline>i_lda12.AND.i_logspline-1<iend) iend=i_logspline-1
   IF (i_rsoptim  >i_lda12.AND.i_rsoptim  -1<iend) iend=i_rsoptim  -1
   IF (i_author   >i_lda12.AND.i_author   -1<iend) iend=i_author   -1
   IF (i_comment  >i_lda12.AND.i_comment  -1<iend) iend=i_comment  -1
   inputstring="";IF (iend>i_lda12+5) inputstring=TRIM(inputline(i_lda12+5:iend))
   IF (inputstring/="") THEN
     CALL extractword(1,inputstring,inputword);inputword=TRIM(inputword)
     IF (inputword/="") THEN
       READ(inputword,*) dataset%xml_lda12_orb_n
       CALL extractword(2,inputstring,inputword);inputword=TRIM(inputword)
       IF (inputword/="") THEN
         READ(inputword,*) dataset%xml_lda12_orb_l
         CALL extractword(3,inputstring,inputword);inputword=TRIM(inputword)
         IF (inputword/="") THEN
           READ(inputword,*) dataset%xml_lda12_ion
           CALL extractword(4,inputstring,inputword);inputword=TRIM(inputword)
           IF (inputword/="") READ(inputword,*) dataset%xml_lda12_rcut
         END IF
       END IF
     END IF
   END IF
 END IF

!Options for related to the transfer to a reduced logarihmic grid
 dataset%xml_usespl=MERGE(.true.,USELOG_DEF,i_logspline>0)
 dataset%xml_spl_meshsz=LOGGRD_SIZE_DEF
 IF (dataset%xml_usespl) THEN
   iend=200
   IF (i_usexcnhat>i_logspline.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_logspline.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_rsoptim  >i_logspline.AND.i_rsoptim  -1<iend) iend=i_rsoptim  -1
   IF (i_lda12    >i_logspline.AND.i_lda12    -1<iend) iend=i_lda12    -1
   IF (i_author   >i_logspline.AND.i_author   -1<iend) iend=i_author   -1
   IF (i_comment  >i_logspline.AND.i_comment  -1<iend) iend=i_comment  -1
   inputstring="";IF (iend>i_logspline+11) inputstring=TRIM(inputline(i_logspline+11:iend))
   IF (inputstring/="") THEN
     CALL extractword(1,inputstring,inputword);inputword=TRIM(inputword)
     IF (inputword/="") READ(inputword,*) dataset%xml_spl_meshsz
   END IF
 END IF

!Author to be mentioned in the XML file
 dataset%xml_author=TRIM(AUTHOR_DEF)
 IF (i_author>0) then
   iend=200
   IF (i_usexcnhat>i_author.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_author.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_logspline>i_author.AND.i_logspline-1<iend) iend=i_logspline-1
   IF (i_rsoptim  >i_author.AND.i_rsoptim  -1<iend) iend=i_rsoptim  -1
   IF (i_lda12    >i_author.AND.i_lda12    -1<iend) iend=i_lda12    -1
   IF (i_comment  >i_author.AND.i_comment  -1<iend) iend=i_comment  -1
   inputstring=TRIM(inputline(i_author+6:))
   READ(unit=inputstring,fmt=*) dataset%xml_author
   dataset%xml_author=TRIM(dataset%xml_author)
 END IF

!Comment to be added in the XML file
 dataset%xml_comment=TRIM(COMMENT_DEF)
 IF (i_comment>0) then
   iend=200
   IF (i_usexcnhat>i_comment.AND.i_usexcnhat-1<iend) iend=i_usexcnhat-1
   IF (i_prtcorewf>i_comment.AND.i_prtcorewf-1<iend) iend=i_prtcorewf-1
   IF (i_logspline>i_comment.AND.i_logspline-1<iend) iend=i_logspline-1
   IF (i_rsoptim  >i_comment.AND.i_rsoptim  -1<iend) iend=i_rsoptim  -1
   IF (i_lda12    >i_comment.AND.i_lda12    -1<iend) iend=i_lda12    -1
   IF (i_author   >i_comment.AND.i_author   -1<iend) iend=i_author   -1
   inputstring=TRIM(inputline(i_comment+7:))
   READ(unit=inputstring,fmt=*) dataset%xml_comment
   dataset%xml_comment=TRIM(dataset%xml_comment)
 END IF

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(/,2x,a)')'Options for XML file:'
   WRITE(STD_OUT,'(2x,2a)') 'XML option : output of core wave funcs =',MERGE("YES"," NO",dataset%xml_prtcorewf)
   WRITE(STD_OUT,'(2x,2a)') 'XML option : spline to a log grid      =',MERGE("YES"," NO",dataset%xml_usespl)
   IF (dataset%xml_usespl) THEN
     WRITE(STD_OUT,'(5x,a,i5)') '- mesh size of the log grid =',dataset%xml_spl_meshsz
   END IF
   WRITE(STD_OUT,'(2x,2a)') 'XML option : Real Space Optimization   =',MERGE("YES"," NO",dataset%xml_userso)
   IF (dataset%xml_userso) THEN
     WRITE(STD_OUT,'(5x,a,f7.2)') '- RSO plane wave cutoff =',dataset%xml_rso_ecut
     WRITE(STD_OUT,'(5x,a,1x,f6.2)') '- RSO Gamma/Gmax ratio  =',dataset%xml_rso_gfact
     WRITE(STD_OUT,'(5x,a,g11.3)') '- RSO maximum error     = ',dataset%xml_rso_werror
   END IF
   WRITE(STD_OUT,'(2x,2a)') 'XML option : output of LDA-1/2 pot.    =',MERGE("YES"," NO",dataset%xml_uselda12)
   IF (dataset%xml_uselda12) THEN
     WRITE(STD_OUT,'(5x,a,i2)') '- LDA-1/2 ionized orbital n =',dataset%xml_lda12_orb_n
     WRITE(STD_OUT,'(5x,a,i2)') '- LDA-1/2 ionized orbital l =',dataset%xml_lda12_orb_l
     WRITE(STD_OUT,'(5x,a,f5.2)') '- LDA-1/2 ionization        =',dataset%xml_lda12_ion
     WRITE(STD_OUT,'(5x,a,f7.4)') '- LDA-1/2 cutoff radius (au)=',dataset%xml_lda12_rcut
     WRITE(STD_OUT,'(5x,3a)') '- LDA-1/2: see ''',trim(dataset%xml_lda12_logfile),&
&                             ''' file to check convergence'
   END IF
   WRITE(STD_OUT,'(2x,2a)') 'XML option : use of compensation in XC =',MERGE("YES"," NO",dataset%xml_usexcnhat)
   IF (dataset%xml_usexcnhat) THEN
     WRITE(STD_OUT,'(5x,a)') '- Zero potential and Kresse local ionic potential output in XML file'
   ELSE
     WRITE(STD_OUT,'(5x,a)') '- Zero potential and Blochl local ionic potential output in XML file'
   END IF
 END IF

 END SUBROUTINE input_dataset_read_xml


!!=================================================================
!! NAME
!! input_dataset_read_upf
!!
!! FUNCTION
!! Read the input file (or standard input) in order to get XML options.
!! Put them in a input_dataset datastructure.
!! If file is omitted, then read from standard input.
!! File is supposed to be opened.
!!
!! INPUTS
!!  [inputfile_unit]= logical unit of input file to be read (optional)
!!  [echofile_unit]= logical unit of a file to echo input file content (optional)
!!
!! OUTPUT
!!  [input_dt]= data-structure containing the complete input file (optional).
!!              If omitted, then the global public `input_dataset`
!!    %upf_grid_xmin= minimum radius given by the grid
!!    %upf_grid_zmesh= inverse of the radial step of the grid
!!    %upf_grid_dx= logarithmic step of the grid
!!    %upf_grid_range= range of the grid
!!  [upf_string]= character string containing the ABINIT options line from input file
!!
!!=================================================================
 SUBROUTINE input_dataset_read_upf(input_dt,inputfile_unit,echofile_unit,upf_string)

!---- Arguments
 INTEGER,INTENT(IN),OPTIONAL :: inputfile_unit,echofile_unit
 CHARACTER(*),INTENT(OUT),OPTIONAL :: upf_string
 TYPE(input_dataset_t),INTENT(INOUT),OPTIONAL,TARGET :: input_dt

!---- Local variables
 INTEGER :: input_unit,echo_unit
 INTEGER :: i_all,i_dx,i_xmin,i_zmesh,i_range
 LOGICAL :: has_to_echo
 CHARACTER(200) :: inputline
 TYPE(input_dataset_t),POINTER :: dataset

!------------------------------------------------------------------

!Select datastruture to be read
 IF (PRESENT(input_dt)) THEN
   dataset => input_dt
 ELSE
   dataset => input_dataset
 ENDIF

!Select input file logical unit
 input_unit=STD_IN;IF (PRESENT(inputfile_unit)) input_unit=inputfile_unit

!Do we echo input file content?
 has_to_echo=PRESENT(echofile_unit)
 echo_unit=-1;IF (has_to_echo) echo_unit=echofile_unit

!Reading of keywords
 READ(input_unit,'(a)') inputline
 IF (has_to_echo) WRITE(echo_unit,'(a)') TRIM(inputline)
 IF (PRESENT(upf_string)) upf_string=TRIM(inputline)
 CALL Uppercase(inputline)

!Default values
 dataset%upf_grid_dx=UPF_DX_DEF
 dataset%upf_grid_xmin=UPF_XMIN_DEF
 dataset%upf_grid_zmesh=UPF_ZMESH_DEF
 dataset%upf_grid_range=UPF_RANGE_DEF

 i_all=INDEX(inputline,'UPFDX')+INDEX(inputline,'UPFXMIN')+&
&      INDEX(inputline,'UPFZMESH')+INDEX(inputline,'UPFRANGE')
 IF (i_all>0) then

!  Logarithmic step
   i_dx=INDEX(inputline,'UPFDX')
   IF (i_dx>0) READ(inputline(i_dx+5:256),*) dataset%upf_grid_dx

!  Minimum radius
   i_xmin=INDEX(inputline,'UPFXMIN')
   IF (i_xmin>0) READ(inputline(i_xmin+7:256),*) dataset%upf_grid_xmin

!  Inverse of radial step
   i_zmesh=INDEX(inputline,'UPFZMESH')
   IF (i_zmesh>0) READ(inputline(i_zmesh+8:256),*) dataset%upf_grid_zmesh

!  Range
   i_range=INDEX(inputline,'UPFRANGE')
   IF (i_range>0) READ(inputline(i_range+8:256),*) dataset%upf_grid_range

 END IF

!Print read data
 IF (has_to_print) THEN
   WRITE(STD_OUT,'(/,2x,a)')      'Options for UPF file:'
   WRITE(STD_OUT,'(2x,a,es10.3)') 'UPF grid option : logarithmic step (upfdx)      =',dataset%upf_grid_dx
   WRITE(STD_OUT,'(2x,a,es10.3)') 'UPF grid option : minimum grid x (upfxmin)      =',dataset%upf_grid_xmin
   WRITE(STD_OUT,'(2x,a,es10.3)') 'UPF grid option : inv. of radial step (upfzmesh)=',dataset%upf_grid_zmesh
   WRITE(STD_OUT,'(2x,a,es10.3)') 'UPF grid option : grid range (upf_gridrange)    =',dataset%upf_grid_range
 END IF

 END SUBROUTINE input_dataset_read_upf

!!=================================================================

END MODULE input_dataset_mod
