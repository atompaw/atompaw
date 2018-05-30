!! NAME
!!  libxc_mod
!!
!! FUNCTION
!!  This module contains routines using features from libXC library
!!  (http://www.tddft.org/programs/octopus/wiki/index.php/Libxc)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2010 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module libxc_mod

 use Tools
 use globalmath

!ISO C bindings are mandatory
#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#endif

 implicit none

!!=================================================================
!! CONSTANTS
!!=================================================================

#if defined HAVE_LIBXC
 logical,parameter,public :: have_libxc=.true.
#else
 logical,parameter,public :: have_libxc=.false.
#endif

!Public constants (use libxc_functionals_constants_load to init them)
 integer,public,save :: XC_FAMILY_UNKNOWN       = -1
 integer,public,save :: XC_FAMILY_LDA           =  1
 integer,public,save :: XC_FAMILY_GGA           =  2
 integer,public,save :: XC_FAMILY_MGGA          =  4
 integer,public,save :: XC_FAMILY_LCA           =  8
 integer,public,save :: XC_FAMILY_OEP           = 16
 integer,public,save :: XC_FAMILY_HYB_GGA       = 32
 integer,public,save :: XC_FAMILY_HYB_MGGA      = 64
 integer,public,save :: XC_FLAGS_HAVE_EXC       =  1
 integer,public,save :: XC_FLAGS_HAVE_VXC       =  2
 integer,public,save :: XC_FLAGS_HAVE_FXC       =  4
 integer,public,save :: XC_FLAGS_HAVE_KXC       =  8
 integer,public,save :: XC_FLAGS_HAVE_LXC       = 16
 integer,public,save :: XC_FLAGS_NEEDS_LAPLACIAN= 32768
 integer,public,save :: XC_EXCHANGE             =  0
 integer,public,save :: XC_CORRELATION          =  1
 integer,public,save :: XC_EXCHANGE_CORRELATION =  2
 integer,public,save :: XC_KINETIC              =  3
 integer,public,save :: XC_SINGLE_PRECISION     =  0
 logical,private,save :: libxc_constants_initialized=.false.

!!=================================================================
!! FUNCTIONS
!!=================================================================

!Public functions
 public :: libxc_init_func       ! Initialize libXC functional(s)
 public :: libxc_end_func        ! Destroy libXC functional(s)
 public :: libxc_print_func      ! Print libXC functionnal(s) details
 public :: libxc_getid_fromInput ! From XC input string, gives the libXC id(s)
 public :: libxc_getid_fromName  ! From a char. string, gives the libXC id
 public :: libxc_getshortname    ! Gives a short version of the libXC name (w/o XC_)
 public :: libxc_getid           ! From LibXC datastructure, gives the libXC id(s)
 public :: libxc_family_from_id  ! Retrieve family of a XC functional from its libXC id
 public :: libxc_islda           ! Return TRUE if the XC functional is LDA
 public :: libxc_isgga           ! Return TRUE if the XC functional is GGA
 public :: libxc_ismgga          ! Return TRUE if the XC functional is MGGA
 public :: libxc_needlap       ! Return TRUE if functional uses LAPLACIAN 
 public :: libxc_ishybrid        ! Return TRUE if the XC functional is Hybrid
 public :: libxc_getvxc          ! Return XC potential and energy, from input density

!Private functions
 private :: libxc_check          ! Check if the code has been compiled with libXC
 private :: libxc_constants_load ! Load libXC constants from C headers
#if defined HAVE_FC_ISO_C_BINDING
 private :: xc_char_to_c         ! Convert a string from Fortran to C
 private :: xc_char_to_f         ! Convert a string from C to Fortran
#endif

!!=================================================================
!! GLOBAL VARIABLE FOR XC FUNCTIONAL
!!=================================================================

!XC functional public type
 type,public :: libxc_functional_t
   integer  :: id              ! identifier
   integer  :: family          ! LDA, GGA, etc.
   integer  :: xckind          ! EXCHANGE, CORRELATION, etc.
   integer  :: nspin           ! # of spin components
   integer  :: abi_ixc         ! Abinit IXC id for this functional
   logical  :: has_exc         ! TRUE is exc is available for the functional
   logical  :: has_vxc         ! TRUE is vxc is available for the functional
   logical  :: has_fxc         ! TRUE is fxc is available for the functional
   logical  :: has_kxc         ! TRUE is kxc is available for the functional
   logical  :: needs_lapl      ! TRUE is functional needs laplacian of n
   real(8) :: hyb_mixing      ! Hybrid functional: mixing factor of Fock contribution (default=0)
   real(8) :: hyb_mixing_sr   ! Hybrid functional: mixing factor of SR Fock contribution (default=0)
   real(8) :: hyb_range       ! Range (for separation) for a hybrid functional (default=0)
#ifdef HAVE_FC_ISO_C_BINDING
   type(C_PTR),pointer :: conf => null() ! C pointer to the functional itself
#endif
 end type libxc_functional_t

!Private global XC functional
 type(libxc_functional_t),target,save,private :: libxc_funcs(2)

!!=================================================================
!! INTERFACES for C BINDINGS
!!=================================================================
#ifdef HAVE_FC_ISO_C_BINDING

 interface
   integer(C_INT) function xc_func_init(xc_func,functional,nspin) bind(C)
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: functional,nspin
     type(C_PTR) :: xc_func
   end function xc_func_init
 end interface

 interface
   subroutine xc_func_end(xc_func) bind(C)
     use iso_c_binding, only : C_PTR
     type(C_PTR) :: xc_func
   end subroutine xc_func_end
 end interface

!
interface
  integer(C_INT) function xc_functional_get_number(name) bind(C)
  use iso_c_binding, only : C_INT,C_PTR
  type(C_PTR),value :: name
  end function xc_functional_get_number
end interface
!
interface
  type(C_PTR) function xc_functional_get_name(number) bind(C)
  use iso_c_binding, only : C_INT,C_PTR
  integer(C_INT),value :: number
  end function xc_functional_get_name
end interface
!
interface
  integer(C_INT) function xc_family_from_id(id,family,number) bind(C)
  use iso_c_binding, only : C_INT,C_PTR
  integer(C_INT),value :: id
  type(C_PTR),value :: family,number
  end function xc_family_from_id
end interface
!

!added by NAWH
interface
  subroutine xc_func_set_dens_threshold(xc_func, dens_threshold) bind(C)
    use iso_c_binding, only : C_DOUBLE,C_PTR
    type(C_PTR) :: xc_func
    real(C_double) :: dens_threshold
    end subroutine xc_func_set_dens_threshold
end interface

interface
  subroutine xc_hyb_cam_coef(xc_func,omega,alpha,beta) bind(C)
  use iso_c_binding, only : C_DOUBLE,C_PTR
  real(C_DOUBLE) :: omega,alpha,beta
  type(C_PTR) :: xc_func
  end subroutine xc_hyb_cam_coef
end interface
!
 interface
   subroutine xc_lda(xc_func,np,rho,zk,vrho,v2rho2,v3rho3) bind(C)
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: np
     type(C_PTR),value :: rho,zk,vrho,v2rho2,v3rho3
     type(C_PTR) :: xc_func
   end subroutine xc_lda
 end interface

 interface
   subroutine xc_gga(xc_func,np,rho,sigma,zk,vrho,vsigma,v2rho2,v2rhosigma,v2sigma2, &
&                    v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3) bind(C)
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: np
     type(C_PTR),value :: rho,sigma,zk,vrho,vsigma,v2rho2,v2rhosigma,v2sigma2, &
&                         v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3
     type(C_PTR) :: xc_func
   end subroutine xc_gga
 end interface

interface
  subroutine xc_mgga(xc_func,np,rho,sigma,lapl,tau,zk,vrho,vsigma,vlapl,vtau, &
& v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,v2rhotau, &
&                     v2sigmalapl,v2sigmatau,v2lapltau) bind(C)
  use iso_c_binding, only : C_INT,C_PTR
  integer(C_INT),value :: np
  type(C_PTR),value :: rho,sigma,lapl,tau,zk,vrho,vsigma,vlapl,vtau, &
& v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,v2rhotau, &
&                         v2sigmalapl,v2sigmatau,v2lapltau
  type(C_PTR) :: xc_func
  end subroutine xc_mgga
end interface

 interface
   subroutine xc_func_set_params(xc_func,params,n_params) bind(C)
     use iso_c_binding, only : C_INT,C_DOUBLE,C_PTR
     integer(C_INT),value :: n_params
     real(C_DOUBLE) :: params(*)
     type(C_PTR) :: xc_func
   end subroutine xc_func_set_params
 end interface

interface
  subroutine xc_mgga_x_tb09_set_params(xc_func,c) bind(C)
  use iso_c_binding, only : C_DOUBLE,C_PTR
  real(C_DOUBLE),value :: c
  type(C_PTR) :: xc_func
  end subroutine xc_mgga_x_tb09_set_params
end interface



 interface
   subroutine xc_get_singleprecision_constant(xc_cst_singleprecision) bind(C)
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_singleprecision
   end subroutine xc_get_singleprecision_constant
 end interface

 interface
   subroutine xc_get_family_constants(xc_cst_unknown,xc_cst_lda,xc_cst_gga,xc_cst_mgga, &
&                                     xc_cst_lca,xc_cst_oep,xc_cst_hyb_gga,xc_cst_hyb_mgga) &
&                                     bind(C)
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_unknown,xc_cst_lda,xc_cst_gga,xc_cst_mgga, &
&                      xc_cst_lca,xc_cst_oep,xc_cst_hyb_gga,xc_cst_hyb_mgga
   end subroutine xc_get_family_constants
 end interface

 interface
   subroutine xc_get_flags_constants(xc_cst_flags_have_exc,xc_cst_flags_have_vxc, &
&   xc_cst_flags_have_fxc,xc_cst_flags_have_kxc, &
&   xc_cst_flags_have_lxc,xc_cst_flags_needs_laplacian) bind(C)
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_flags_have_exc,xc_cst_flags_have_vxc,xc_cst_flags_have_fxc, &
&  xc_cst_flags_have_kxc,xc_cst_flags_have_lxc,xc_cst_flags_needs_laplacian
   end subroutine xc_get_flags_constants
 end interface

 interface
   subroutine xc_get_kind_constants(xc_cst_exchange,xc_cst_correlation, &
&                                   xc_cst_exchange_correlation,xc_cst_kinetic) bind(C)
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_exchange,xc_cst_correlation, &
&                      xc_cst_exchange_correlation,xc_cst_kinetic
   end subroutine xc_get_kind_constants
 end interface

 interface
   type(C_PTR) function xc_func_type_malloc() bind(C)
     use iso_c_binding, only : C_PTR
   end function xc_func_type_malloc
 end interface

 interface
   subroutine xc_func_type_free(xc_func) bind(C)
     use iso_c_binding, only : C_PTR
     type(C_PTR) :: xc_func
   end subroutine xc_func_type_free
 end interface

 interface
   type(C_PTR) function xc_get_info_name(xc_func) bind(C)
     use iso_c_binding, only : C_PTR
     type(C_PTR) :: xc_func
   end function xc_get_info_name
 end interface

 interface
   type(C_PTR) function xc_get_info_refs(xc_func,iref) bind(C)
     use iso_c_binding, only : C_INT,C_PTR
     type(C_PTR) :: xc_func
     integer(C_INT) :: iref
   end function xc_get_info_refs
 end interface

 interface
   integer(C_INT) function xc_get_info_flags(xc_func) bind(C)
     use iso_c_binding, only : C_INT,C_PTR
     type(C_PTR) :: xc_func
   end function xc_get_info_flags
 end interface

 interface
   integer(C_INT) function xc_get_info_kind(xc_func) bind(C)
     use iso_c_binding, only : C_INT,C_PTR
     type(C_PTR) :: xc_func
   end function xc_get_info_kind
 end interface

#endif

 private

 CONTAINS


!!=================================================================
!! NAME
!!  libxc_getid_fromInput
!!
!! FUNCTION
!!  From a character string (as given in ATOMPAW input file),
!!   gives the libXC id(s)
!!
!! INPUTS
!!  xcname= string containing the name of a XC functional
!!
!! OUTPUT
!!  id(2)= id(s) of the libXC functional
!!  [xcname_short]= short version of the libXC name (optional)
!!
!!=================================================================
 subroutine libxc_getid_fromInput(xcname,id,xcname_short)

 implicit none
 integer,intent(inout) :: id(2)
 character*(*),intent(in) :: xcname
 character*(*),intent(out),optional :: xcname_short

#if defined HAVE_LIBXC

!---- Local variables
 integer :: ii,i_plus
 character*50 :: xcstrg(2)

!------------------------------------------------------------------
!---- Executable code

 i_plus=index(xcname,'+')
 if (i_plus<=0) then
  xcstrg(1)=trim(xcname)
  xcstrg(2)=""
 else
  xcstrg(1)=trim(xcname(1:i_plus-1))
  xcstrg(2)=trim(xcname(i_plus+1:))
 end if

 do ii=1,2
  id(ii)=-1
  call uppercase(xcstrg(ii))

  id(ii)=libxc_getid_fromName(xcstrg(ii))

  if (id(ii)==-1.and.ii==2) exit

  if (id(ii)==-1.and.xcstrg(ii)(1:6)=="LIBXC_") then
   read(unit=xcstrg(ii)(7:),fmt=*,err=333,end=333) id(ii)
333 continue
  end if

  if (id(ii)==-1) then
   write(6,'(/,2x,a)') "Error in get_id_from_name:"
   write(6,'(2x,3a)')  " Unknown X, C or XC functionnal (", &
&                     trim(xcstrg(ii)),") !"
   stop
  end if
 end do

 if (present(xcname_short)) then
  xcname_short=""
  if (id(1)>0.and.xcstrg(1)(1:3)=="XC_")    xcname_short=trim(xcname_short)     //trim(xcstrg(1)(4:))
  if (id(1)>0.and.xcstrg(1)(1:6)=="LIBXC_") xcname_short=trim(xcname_short)     //trim(xcstrg(1))
  if (id(2)>0.and.xcstrg(2)(1:3)=="XC_")    xcname_short=trim(xcname_short)//"+"//trim(xcstrg(2)(4:))
  if (id(2)>0.and.xcstrg(2)(1:6)=="LIBXC_") xcname_short=trim(xcname_short)//"+"//trim(xcstrg(2))
  if (trim(xcname_short)=="") xcname_short=xcname
 end if

#else
 id(1:2)=-2
 if (present(xcname_short)) xcname_short=xcname
#endif

 end subroutine libxc_getid_fromInput


!!=================================================================
!! NAME
!!  libxc_getid_fromName
!!
!! FUNCTION
!!  Return identifer of a XC functional from its name
!!  Return -1 if undefined
!!
!! INPUTS
!!  xcname= string containing the name of a XC functional
!!
!!=================================================================
 function libxc_getid_fromName(xcname)

 implicit none
 integer :: libxc_getid_fromName
 character(len=*),intent(in) :: xcname

#if defined HAVE_LIBXC

!---- Local variables
#if defined HAVE_FC_ISO_C_BINDING
 character(len=256) :: str
 character(kind=C_CHAR,len=1),target :: name_c(len_trim(xcname)+1)
 character(kind=C_CHAR,len=1),target :: name_c_xc(len_trim(xcname)-2)
 type(C_PTR) :: name_c_ptr
#endif

!------------------------------------------------------------------
!---- Executable code

#if defined HAVE_FC_ISO_C_BINDING
 str=trim(xcname)
 if (xcname(1:3)=="XC_".or.xcname(1:3)=="xc_") then
   str=xcname(4:);name_c_xc=xc_char_to_c(str)
   name_c_ptr=c_loc(name_c_xc)
 else
   name_c=xc_char_to_c(str)
   name_c_ptr=c_loc(name_c)
 end if
 libxc_getid_fromName=int(xc_functional_get_number(name_c_ptr))
#endif

#else
 libxc_getid_fromName=-1
#endif

end function libxc_getid_fromName


!!=================================================================
!! NAME
!!  libxc_getid
!!
!! FUNCTION
!!  From LibXC datastructure, gives the libXC id(s)
!!
!! OUTPUT
!!  id(2)= id(s) of the XC functional
!!
!!=================================================================
 subroutine libxc_getid(id)

 implicit none
 integer :: id(2)

#if defined HAVE_LIBXC

!---- Local variables
 integer :: ii

!------------------------------------------------------------------
!---- Executable code

 do ii=1,2
  id(ii)=libxc_funcs(ii)%id
 end do

#else
 id(1:2)=-2
#endif

 end subroutine libxc_getid


!!=================================================================
!! NAME
!!  libxc_getshortname
!!
!! FUNCTION
!!  From a character string (given in input file), gives a short
!!  version of the libXC name (without XC_)
!!
!! INPUTS
!!  xcname= string containing the name of a XC functional
!!
!! OUTPUT
!!  xcname_short= short version of the libXC name
!!
!!=================================================================
 subroutine libxc_getshortname(xcname,xcname_short)

 implicit none
 character*(*),intent(in) :: xcname
 character*(*),intent(out) :: xcname_short

#if defined HAVE_LIBXC
!---- Local variables
 integer :: i_plus
 character*50 :: xcstrg(2)

!------------------------------------------------------------------
!---- Executable code

 i_plus=index(xcname,'+')
 if (i_plus<=0) then
  xcstrg(1)=trim(xcname)
  xcstrg(2)=""
 else
  xcstrg(1)=trim(xcname(1:i_plus-1))
  xcstrg(2)=trim(xcname(i_plus+1:))
 end if
 call uppercase(xcstrg(1))
 call uppercase(xcstrg(2))

 xcname_short=""
 if (xcstrg(1)(1:3)=="XC_")    xcname_short=trim(xcname_short)     //trim(xcstrg(1)(4:))
 if (xcstrg(1)(1:6)=="LIBXC_") xcname_short=trim(xcname_short)     //trim(xcstrg(1))
 if (xcstrg(2)(1:3)=="XC_")    xcname_short=trim(xcname_short)//"+"//trim(xcstrg(2)(4:))
 if (xcstrg(2)(1:6)=="LIBXC_") xcname_short=trim(xcname_short)//"+"//trim(xcstrg(2))
 if (trim(xcname_short)=="") xcname_short=xcname

#else
 xcname_short=xcname
#endif

 end subroutine libxc_getshortname


!!=================================================================
!! NAME
!! libxc_init_func
!!
!! FUNCTION
!! Initialize libXC functional(s)
!!
!! INPUTS
!!  id(2)= libXC id(s) a XC functional
!!  nsp=number of psin component
!!
!!=================================================================
 subroutine libxc_init_func(id,nsp)

 implicit none
 integer,intent(in)  :: id(2)
 integer,intent(in) :: nsp

!---- Local variables
 integer :: ii
 type(libxc_functional_t),pointer :: xc_func
#if defined HAVE_LIBXC && defined HAVE_FC_ISO_C_BINDING
 integer :: flags
 integer(C_INT) :: func_id_c,npar_c,nspin_c,success_c
 real(C_DOUBLE) :: alpha_c,beta_c,omega_c
 real(C_DOUBLE) :: param_c(3)
 character(kind=C_CHAR,len=1),pointer :: strg_c
 type(C_PTR) :: func_ptr_c
#endif

!------------------------------------------------------------------
!---- Executable code

!Check libXC
 if (.not.libxc_check(stop_if_error=.true.)) return
 if (.not.libxc_constants_initialized) call libxc_constants_load()

 libxc_funcs(1:2)%id=id(1:2)

 do ii=1,2

   xc_func => libxc_funcs(ii)

   xc_func%family=XC_FAMILY_UNKNOWN
   xc_func%xckind=-1
   xc_func%nspin=nsp
   xc_func%has_exc=.false.
   xc_func%has_vxc=.false.
   xc_func%has_fxc=.false.
   xc_func%has_kxc=.false.
   xc_func%needs_lapl=.false.
   xc_func%hyb_mixing=0.d0
   xc_func%hyb_mixing_sr=0.d0
   xc_func%hyb_range=0.d0

   if (xc_func%id<=0) cycle

!  Get XC functional family
   libxc_funcs%family=libxc_family_from_id(xc_func%id)
!   live dangerously
   write(6,*) 'The LibXC functional family ',xc_func%family
!   if (xc_func%family/=XC_FAMILY_LDA.and.xc_func%family/=XC_FAMILY_GGA) then
 !    write(6,'(a,i4,a)') 'The LibXC functional family ',xc_func%family, &
!&                        ' is currently unsupported by ATOMPAW!'
!     write(6,'(a)') '(-1 means the family is unknown to the LibXC itself)'
!     write(6,'(a)') 'Please consult the LibXC documentation.'
!     stop
 !  end if

#if defined HAVE_LIBXC && defined HAVE_FC_ISO_C_BINDING

!  Allocate functional
   func_ptr_c=xc_func_type_malloc()
   call c_f_pointer(func_ptr_c,xc_func%conf)

!  Initialize functional
   func_id_c=int(xc_func%id,kind=C_INT)
   nspin_c=int(nsp,kind=C_INT)
   success_c=xc_func_init(xc_func%conf,func_id_c,nspin_c)
   if (success_c/=0) then
     write(6,'(a)') 'Error in libXC functional initialization!'
     stop
   end if

!  Special treatment for LDA_C_XALPHA functional
   if (xc_func%id==libxc_getid_fromName('XC_LDA_C_XALPHA')) then
     param_c(1)=real(0.d0,kind=C_DOUBLE);npar_c=int(1,kind=C_INT)
     call xc_func_set_params(xc_func%conf,param_c,npar_c)
   end if

!  Get functional kind
   xc_func%xckind=int(xc_get_info_kind(xc_func%conf))

!  Get functional flags
   flags=int(xc_get_info_flags(xc_func%conf))
   xc_func%has_exc=(iand(flags,XC_FLAGS_HAVE_EXC)>0)
   xc_func%has_vxc=(iand(flags,XC_FLAGS_HAVE_VXC)>0)
   xc_func%has_fxc=(iand(flags,XC_FLAGS_HAVE_FXC)>0)
   xc_func%has_kxc=(iand(flags,XC_FLAGS_HAVE_KXC)>0)
   xc_func%needs_lapl=(iand(flags,XC_FLAGS_NEEDS_LAPLACIAN)>0)

!  Retrieve parameters for hybrid functionals
   if (xc_func%family==XC_FAMILY_HYB_GGA) then
     call xc_hyb_cam_coef(xc_func%conf,omega_c,alpha_c,beta_c)
     xc_func%hyb_mixing=real(alpha_c,kind=8)
     xc_func%hyb_mixing_sr=real(beta_c,kind=8)
     xc_func%hyb_range=real(omega_c,kind=8)
   endif

#endif

 end do ! loop on functionals

!Print functional(s) information
 call libxc_print_func(6)

 end subroutine libxc_init_func


!!=================================================================
!! NAME
!! libxc_end_func
!!
!! FUNCTION
!!  Destroy libXC functional(s)
!!
!!=================================================================
 subroutine libxc_end_func()

  implicit none

#if defined HAVE_LIBXC

!---- Local variables
 integer :: ii
 type(libxc_functional_t),pointer :: xc_func

!------------------------------------------------------------------
!---- Executable code

 do ii = 1,2

   xc_func => libxc_funcs(ii)

   if (xc_func%id <= 0) cycle
   xc_func%id=-1
   xc_func%family=-1
   xc_func%xckind=-1
   xc_func%nspin=1
   xc_func%abi_ixc=huge(0)
   xc_func%has_exc=.false.
   xc_func%has_vxc=.false.
   xc_func%has_fxc=.false.
   xc_func%has_kxc=.false.
   xc_func%hyb_mixing=0.d0
   xc_func%hyb_mixing_sr=0.d0
   xc_func%hyb_range=0.d0
#if defined HAVE_FC_ISO_C_BINDING
   if (associated(xc_func%conf)) then
     call xc_func_end(xc_func%conf)
     call xc_func_type_free(c_loc(xc_func%conf))
   end if
#endif

 end do

#endif

 end subroutine libxc_end_func


!!=================================================================
!! NAME
!! libxc_print_func
!!
!! FUNCTION
!!  Print libXC functionnal(s) details
!!
!! INPUTS
!!  unt= logical unit to print
!!
!!=================================================================
 subroutine libxc_print_func(unt)

 implicit none
 integer :: unt

#if defined HAVE_LIBXC

!---- Local variables
 integer :: ii
 character(len=500) :: msg
 type(libxc_functional_t),pointer :: xc_func
#if defined HAVE_FC_ISO_C_BINDING
 integer(C_INT) :: iref_c
 character(kind=C_CHAR,len=1),pointer :: strg_c
#endif

!------------------------------------------------------------------
!---- Executable code

 do ii=1,2

   xc_func => libxc_funcs(ii)
   if (xc_func%id<=0) cycle

   if (xc_func%xckind==XC_EXCHANGE) then
     write(unt,'(a)') 'Exchange functional (LibXC):'
   else if (xc_func%xckind==XC_CORRELATION) then
     write(unt,'(a)') 'Correlation functional (LibXC):'
   else if (xc_func%xckind==XC_EXCHANGE_CORRELATION) then
     write(unt,'(a)') 'Exchange-Correlation functional (LibXC):'
   end if

#if defined HAVE_FC_ISO_C_BINDING
   call c_f_pointer(xc_get_info_name(xc_func%conf),strg_c)
   call xc_char_to_f(strg_c,msg)

   iref_c=0
   do while (iref_c>=0)
     call c_f_pointer(xc_get_info_refs(xc_func%conf,iref_c),strg_c)
     if (associated(strg_c)) then
       call xc_char_to_f(strg_c,msg)
       write(unt,'(2x,a)') trim(msg)
       iref_c=iref_c+1
     else
       iref_c=-1
     end if
   end do
#endif

 end do

#endif

 end subroutine libxc_print_func


!!=================================================================
!!  libxc_family_from_id
!!
!! FUNCTION
!!  Return family of a XC functional from its id
!!
!! INPUTS
!!  xcid= id of a LibXC functional
!!
!!=================================================================
 function libxc_family_from_id(xcid)

 implicit none
 integer :: libxc_family_from_id
 integer,intent(in) :: xcid

#if defined HAVE_LIBXC

!---- Local variables
#if defined HAVE_FC_ISO_C_BINDING
 integer(C_INT) :: xcid_c
#endif

!------------------------------------------------------------------
!---- Executable code

#if defined HAVE_FC_ISO_C_BINDING
 xcid_c=int(xcid,kind=C_INT)
 libxc_family_from_id=int(xc_family_from_id(xcid_c,C_NULL_PTR,C_NULL_PTR))
#endif

#else
 libxc_family_from_id=-1
#endif

end function libxc_family_from_id


!!=================================================================
!! NAME
!!  libxc_isgga
!!
!! FUNCTION
!!  Return TRUE is LibXC functional is GGA
!!
!!=================================================================
 function libxc_islda()

 implicit none
 logical :: libxc_islda

!------------------------------------------------------------------
!---- Executable code

 libxc_islda = .false.
 if (.not.libxc_constants_initialized) call libxc_constants_load()
 
 libxc_islda = (any(libxc_funcs(:)%family == XC_FAMILY_LDA))

 end function libxc_islda

!!=================================================================
!!=================================================================
!! NAME
!!  libxc_isgga
!!
!! FUNCTION
!!  Return TRUE is LibXC functional is GGA
!!
!!=================================================================
 function libxc_isgga()

 implicit none
 logical :: libxc_isgga

!------------------------------------------------------------------

 libxc_isgga = .false.
 if (.not.libxc_constants_initialized) call libxc_constants_load()

 libxc_isgga = (any(libxc_funcs(:)%family == XC_FAMILY_GGA))

 end function libxc_isgga


!!=================================================================
 function libxc_ismgga()

 implicit none
 logical :: libxc_ismgga

!------------------------------------------------------------------
!---- Executable code

 libxc_ismgga = .false.
 if (.not.libxc_constants_initialized) call libxc_constants_load()

 libxc_ismgga = (any(libxc_funcs(:)%family == XC_FAMILY_MGGA))

 end function libxc_ismgga


!!=================================================================
 function libxc_needlap()

 implicit none
 logical :: libxc_needlap

!------------------------------------------------------------------
!---- Executable code

 libxc_needlap = .false.
 if (.not.libxc_constants_initialized) call libxc_constants_load()

 libxc_needlap = (any(libxc_funcs(:)%needs_lapl))

 end function libxc_needlap

!!=================================================================
!! NAME
!!  libxc_ishybrid
!!
!! FUNCTION
!!  Return TRUE is LibXC functional is Hybrid
!!
!!=================================================================
 function libxc_ishybrid()

 implicit none
 logical :: libxc_ishybrid

!------------------------------------------------------------------
!---- Executable code

 libxc_ishybrid = .false.
 if (.not.libxc_constants_initialized) call libxc_constants_load()

 libxc_ishybrid = (any(libxc_funcs(:)%family == XC_FAMILY_HYB_GGA))

 end function libxc_ishybrid


!!=================================================================
!! NAME
!!  libxc_getvxc
!!
!! FUNCTION
!!  Return XC potential and energy, from input density (gradient etc...)
!!
!! INPUTS
!!  npts= number of of points for the density
!!  nsp= number of spin-density components
!!  rho(npts,nsp)= electronic density
!!  [grho(npts,2*nsp-1)]= gradient of the density (optional)
!!  [lrho(npts,2*nsp-1)]= laplacian of the density (optional)
!!  [tau(npts,2*nsp-1)]= sum of squared gradient of occ wf's (optional)
!!
!! OUTPUT
!!  exc(npts)=XC energy density
!!  vxc(npts,nsp)=derivative of the energy density wrt to the density
!!  [vxcgr(npts,2*nsp-1)]=2nd der. of the energy density wrt to the gradient
!!                        2nd der. of the energy density wrt to the density and the gradient
!!                        (optional)
!!
!! NOTES
!!  nsp=1 : rho is total density (not half)
!!          grho is abs(grad(rho))
!!  nsp=2 : rho is [rho^up,rho^dn]
!!          grho is [abs(grad(rho^up)),abs(grad(rho^dn)),abs(grad(rho^tot))]
!!
!!=================================================================
 subroutine &
&   libxc_getvxc(npts,exc,vxc,nsp,rho,   &
&   grho,lrho,tau,vxcgr,dvxc,d2vxc,vxclrho,vxctau)

 implicit none
 integer, intent(in)            :: npts,nsp
 real(8),intent(inout)          :: exc(npts),vxc(npts,nsp)
 real(8),intent(in)             :: rho(npts,nsp)
 real(8),intent(in),optional    :: grho(npts,2*nsp-1)
 real(8),intent(in),optional    :: lrho(npts,2*nsp-1)
 real(8),intent(in),optional    :: tau(npts,2*nsp-1)
 real(8),intent(inout),optional :: vxcgr(npts,2*nsp-1)
 real(8),intent(inout), optional :: dvxc(npts,2*nsp-1)
 real(8),intent(inout), optional :: d2vxc(npts,2*nsp-1)
 real(8),intent(inout), optional :: vxclrho(npts,2*nsp-1)
 real(8),intent(inout), optional :: vxctau(npts,2*nsp-1)

#if defined HAVE_LIBXC

!---- Local variables
 real(8),parameter :: tol=1.d-14
 
 integer :: ii,ipts,izero
 logical :: is_lda,is_gga,is_mgga,needs_lapl
 real(8),target :: exctmp
 real(8),target :: rhotmp(nsp),vxctmp(nsp),sigma(3),vsigma(3)
 real(8),target :: v2rho2(3),v2rhosigma(6),v2sigma2(6),v3rho3(4)
 real(8),target :: lrhotmp(nsp),tautmp(nsp),vlrho(nsp),vtau(nsp)

#if defined HAVE_LIBXC && defined HAVE_FC_ISO_C_BINDING
 type(C_PTR) :: rho_c,sigma_c,lrho_c,tau_c
 type(C_PTR) :: exc_c(2),vxc_c(2),vsigma_c(2)
 type(C_PTR) :: v2rho2_c(2),v2rhosigma_c(2),v2sigma2_c(2)
 type(C_PTR) :: v3rho3_c(2),vlrho_c(2),vtau_c(2)
#endif

!------------------------------------------------------------------
!---- Executable code

 if (.not.libxc_constants_initialized) call libxc_constants_load()

 is_lda=libxc_islda()
 is_gga=libxc_isgga()
 is_mgga=libxc_ismgga()
 needs_lapl=libxc_needlap()
 if (is_gga.and.((.not.present(vxcgr).or.(.not.present(grho))))) then
   write(6,'(/,2x,a)') "Bug in libxc_getvxc:"
   write(6,'(2x,3a)')  " GGA called without grho or vxcgr!"
   stop
 end if
 if(needs_lapl.and.(.not.present(lrho))) then
    write(6,'((a))') 'getvxc need Laplacian of density error'
    stop
 endif    
! Need to add more consistency tests

!Initialize all output arrays to zero
 exc=0.d0 ; vxc=0.d0
 if (is_gga.or.is_mgga.and.present(vxcgr)) vxcgr=0.d0
 if (present(dvxc)) dvxc=0.d0
 if (present(d2vxc)) d2vxc=0.d0
 if (is_mgga.and.present(vxclrho)) vxclrho=0.d0
 if (is_mgga.and.present(vxctau)) vxctau=0.d0


!Filter density/gradient when density goes to zero
!This is useless ; libxc has its own filtering process
  izero=npts
! do ipts=npts,2,-1
!  if (all(rho(ipts,:)<tol)) izero=ipts-1
! end do

!Define C pointers to libXC routine arguments
!#if defined HAVE_FC_ISO_C_BINDING
! do ii = 1,2
!   if (libxc_funcs(ii)%has_exc) then
!     exc_c(ii)=c_loc(exctmp)
!   else
!     exc_c(ii)=C_NULL_PTR
!   end if
!   if (libxc_funcs(ii)%has_vxc) then
!     vxc_c(ii)=c_loc(vxctmp)
!     vsigma_c(ii)=c_loc(vsigma)
!   else
!     vxc_c(ii)=C_NULL_PTR
!     vsigma_c(ii)=c_NULL_PTR
!   end if
! end do
! rhotmp=0.d0 ; rho_c=c_loc(rhotmp)
! if (is_gga.or.is_mgga) then
!   sigma=0.d0 ; sigma_c=c_loc(sigma)
! endif  
! if (is_mgga) then  
!   tautmp=0.d0; tau_c=c_loc(tautmp)      
!   lrhotmp=0.d0; lrho_c=c_loc(lrhotmp)      
! end if
!#endif

!Define C pointers to libXC routine arguments
#if defined HAVE_FC_ISO_C_BINDING
 do ii = 1,2
     call xc_func_set_dens_threshold(libxc_funcs(ii)%conf, tol)
     exc_c(ii)=c_NULL_PTR
     vxc_c(ii)=c_NULL_PTR
     vsigma_c(ii)=c_NULL_PTR
     v2rho2_c(ii)=c_NULL_PTR
     v2rhosigma_c(ii)=c_NULL_PTR
     v2sigma2_c(ii)=c_NULL_PTR
     v3rho3_c(ii)=c_NULL_PTR
     vlrho_c(ii)=c_NULL_PTR
     vtau_c(ii)=c_NULL_PTR
   if (libxc_funcs(ii)%has_exc) then
       exc_c(ii)=c_loc(exctmp)
       vxc_c(ii)=c_loc(vxctmp)
       vsigma_c(ii)=c_loc(vsigma)
       vlrho_c(ii)=c_loc(vlrho)
       vtau_c(ii)=c_loc(vtau)
   endif    
 enddo  

 rhotmp=0.d0 ; rho_c=c_loc(rhotmp)
 if (is_gga.or.is_mgga) then
   sigma=0.d0 ; sigma_c=c_loc(sigma)
 endif  
 if (is_mgga) then  
   tautmp=0.d0; tau_c=c_loc(tautmp)      
   lrhotmp=0.d0; lrho_c=c_loc(lrhotmp)      
 end if
#endif

!Loop over points
 do ipts=1,npts

   !Load density (and gradient) for this point
   vxctmp=0.d0;exctmp=0.d0
   if (ipts<=izero) then
     rhotmp(1:nsp)=rho(ipts,1:nsp)
   else
     rhotmp=tol
   end if
   if (is_gga.or.is_mgga) then
     if (ipts<=izero) then
       if (nsp==1) then
         !AtomPAW passes |grho| while LibXC needs |grho|^2
         sigma(1)=grho(ipts,1)**2
       else
         !AtomPAW passes |grho_up|, |grho_dn|, and |grho_tot|
         !while Libxc needs |grho_up|^2, grho_up.grho_dn, and |grho_dn|^2
         sigma(1)= grho(ipts,1)**2
         sigma(3)= grho(ipts,2)**2
         sigma(2)=(grho(ipts,3)**2-sigma(1)-sigma(3))*0.5d0
       end if
     else
       sigma=0.d0
     end if
   end if
   if (is_mgga) then
     if (ipts<=izero) then
       if (nsp==1) then
         !AtomPAW passes tau    while LibXC needs tau     
         tautmp(1)=tau(ipts,1)
         if (needs_lapl) lrhotmp(1)=lrho(ipts,1)
       else
         !AtomPAW passes tau_up, tau_dn, and tau_tot
         !while Libxc needs tau_up, tau_dn ??????(wild guess)
         !not at all clear!
         tautmp(1)= tau(ipts,1)
         tautmp(2)= tau(ipts,2)
         if (needs_lapl) then
            lrhotmp(1)= lrho(ipts,1)
            lrhotmp(2)= lrho(ipts,2)
         endif   
       end if
     else
       tautmp=0.d0
       if (needs_lapl) lrhotmp=0.d0
     end if
   end if


! Loop over functionals
  do ii=1,2
    if (libxc_funcs(ii)%id<=0) cycle

!   Get the potential (and possibly the energy)
#if defined HAVE_FC_ISO_C_BINDING
!   ===== LDA =====
    if (libxc_funcs(ii)%family==XC_FAMILY_LDA) then
      exctmp=0.d0 ; vxctmp=0.d0
      call xc_lda(libxc_funcs(ii)%conf,1,rho_c, &
&                 exc_c(ii),vxc_c(ii), &
&                 C_NULL_PTR,C_NULL_PTR)
!   ===== GGA =====
    else if (libxc_funcs(ii)%family==XC_FAMILY_GGA.or. &
&            libxc_funcs(ii)%family==XC_FAMILY_HYB_GGA) then
      exctmp=0.d0 ; vxctmp=0.d0 ; vsigma=0.d0
      call xc_gga(libxc_funcs(ii)%conf,1,rho_c,sigma_c, &
&                 exc_c(ii),vxc_c(ii),vsigma_c(ii), &
&                 C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR, &
&                 C_NULL_PTR,C_NULL_PTR,C_NULL_PTR)
!    ===== mGGA =====
     else if (libxc_funcs(ii)%family==XC_FAMILY_MGGA) then
       exctmp=0.d0 ; vxctmp=0.d0 ; vsigma=0.d0 ; vlrho=0.d0 ; vtau=0.d0
       call xc_mgga(libxc_funcs(ii)%conf,1,rho_c,sigma_c,lrho_c,tau_c, &
&         exc_c(ii),vxc_c(ii),vsigma_c(ii),vlrho_c(ii),vtau_c(ii), &
&         C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,&
&         C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR)
       write(6,'(2i5,2x,1P,20e15.7)') ii,ipts,rhotmp(1),tautmp(1),&
&          lrhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1),&
&          vlrho(1),vtau(1)                   
    end if
#endif

    exc(ipts)=exc(ipts)+2.d0*exctmp ! From Ha to Ry
    vxc(ipts,1:nsp)=vxc(ipts,1:nsp)+2.d0*vxctmp(1:nsp) ! From Ha to Ry
    if (is_gga.or.is_mgga) then
      if (nsp==1) then
        !Note: for nsp=1, vsigma(1) contains 1/2 (up part only)
        !Natalie thinks the extra 2 comes from sigma=grad**2
        vxcgr(ipts,1)=vxcgr(ipts,1)+4.d0*vsigma(1) ! From Ha to Ry
      else
        !Natalie thinks this may be incorrect      
        vxcgr(ipts,2*nsp-1)=vxcgr(ipts,2*nsp-1)+2.d0*vsigma(2*nsp-1) ! From Ha to Ry
      end if
!     if (nsp==1) then
!       vxcgr(ipts,3)=vxcgr(ipts,3)+vsigma(1)*2.d0
!     else
!       vxcgr(ipts,1)=vxcgr(ipts,1)+2.d0*vsigma(1)-vsigma(2)
!       vxcgr(ipts,2)=vxcgr(ipts,2)+2.d0*vsigma(3)-vsigma(2)
!       vxcgr(ipts,3)=vxcgr(ipts,3)+     vsigma(2)
!     end if
    end if
   if (is_mgga) then
       vxctau(ipts,1:nsp)=vxctau(ipts,1:nsp)+2*vtau(1:nsp)  !H to Ry    
       if(needs_lapl) &
&        vxclrho(ipts,1:nsp)=vxclrho(ipts,1:nsp)+2*vlrho(1:nsp)
   end if
  end do ! loop over functional(s)
 end do  ! loop over points

#endif
 end subroutine libxc_getvxc


!!=================================================================
!! NAME
!!  libxc_constants_load
!!
!! FUNCTION
!!  Load libXC constants from C headers
!!
!!=================================================================
 subroutine libxc_constants_load()

 implicit none

#if defined HAVE_LIBXC

!---- Local variables
#if defined HAVE_FC_ISO_C_BINDING
 integer(C_INT) :: i1,i2,i3,i4,i5,i6,i7,i8
#endif

!------------------------------------------------------------------
!---- Executable code

#if defined HAVE_FC_ISO_C_BINDING
  call xc_get_singleprecision_constant(i1)
  XC_SINGLE_PRECISION     = int(i1)
  call xc_get_family_constants(i1,i2,i3,i4,i5,i6,i7,i8)
  XC_FAMILY_UNKNOWN       = int(i1)
  XC_FAMILY_LDA           = int(i2)
  XC_FAMILY_GGA           = int(i3)
  XC_FAMILY_MGGA          = int(i4)
  XC_FAMILY_LCA           = int(i5)
  XC_FAMILY_OEP           = int(i6)
  XC_FAMILY_HYB_GGA       = int(i7)
  XC_FAMILY_HYB_MGGA      = int(i8)
  call xc_get_flags_constants(i1,i2,i3,i4,i5,i6)
  XC_FLAGS_HAVE_EXC       = int(i1)
  XC_FLAGS_HAVE_VXC       = int(i2)
  XC_FLAGS_HAVE_FXC       = int(i3)
  XC_FLAGS_HAVE_KXC       = int(i4)
  XC_FLAGS_HAVE_LXC       = int(i5)
  XC_FLAGS_NEEDS_LAPLACIAN= int(i6)
  call xc_get_kind_constants(i1,i2,i3,i4)
  XC_EXCHANGE             = int(i1)
  XC_CORRELATION          = int(i2)
  XC_EXCHANGE_CORRELATION = int(i3)
  XC_KINETIC              = int(i4)
 libxc_constants_initialized=.true.
#endif

#else
 libxc_constants_initialized=.false.
#endif

 end subroutine libxc_constants_load


!!=================================================================
!! NAME
!!  libxc_check
!!
!! FUNCTION
!!  Check if the code has been compiled with libXC
!!
!! INPUTS
!!  [stop_if_error]=optional flag; if TRUE the code stops
!!                  if libXC is not correctly used
!!
!!=================================================================
 function libxc_check(stop_if_error)

 implicit none
 logical :: libxc_check
 logical,intent(in),optional :: stop_if_error

!---- Local variables
 character(len=100) :: msg

!------------------------------------------------------------------
!---- Executable code

 libxc_check=.true. ; msg=""

#if defined HAVE_LIBXC
#if defined FC_G95
 libxc_check=.false.
 msg='LibXC cannot be used with G95 Fortran compiler!'
#endif
#if defined HAVE_FC_ISO_C_BINDING
 if (.not.libxc_constants_initialized) call libxc_constants_load()
 if (XC_SINGLE_PRECISION==1) then
   libxc_check=.false.
   msg='LibXC should be compiled with double precision!'
 end if
#else
 libxc_check=.false.
 msg='LibXC cannot be used without ISO_C_BINDING support by the Fortran compiler!'
#endif
#else
 libxc_check=.false.
 msg='ATOMPAW was not compiled with LibXC support.'
#endif

 if (present(stop_if_error)) then
   if (stop_if_error.and.trim(msg)/="") then
     write(6,'(a)') trim(msg) ; stop
   end if
 end if

 end function libxc_check


!!=================================================================
!! NAME
!!  xc_char_to_c
!!
!! FUNCTION
!!  Helper function to convert a Fortran string to a C string
!!  Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  f_string=Fortran string
!!
!! OUTPUT
!!  c_string=C string
!!
!!=================================================================
#if defined HAVE_FC_ISO_C_BINDING
function xc_char_to_c(f_string) result(c_string)

 character(len=*),intent(in) :: f_string
 character(kind=C_CHAR,len=1) :: c_string(len_trim(f_string)+1)

!---- Local variables
 integer :: ii,strlen

!------------------------------------------------------------------
!---- Executable code

 strlen=len_trim(f_string)
 forall(ii=1:strlen)
   c_string(ii)=f_string(ii:ii)
 end forall
 c_string(strlen+1)=C_NULL_CHAR

 end function xc_char_to_c
#endif


!!=================================================================
!! NAME
!!  xc_char_to_f
!!
!! FUNCTION
!!  Helper function to convert a C string to a Fortran string
!!  Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  c_string=C string
!!
!! OUTPUT
!!  f_string=Fortran string
!!
!!=================================================================
#if defined HAVE_FC_ISO_C_BINDING
subroutine xc_char_to_f(c_string,f_string)

 character(kind=C_CHAR,len=1),intent(in) :: c_string(*)
 character(len=*),intent(out) :: f_string

!---- Local variables
 integer :: ii

!------------------------------------------------------------------
!---- Executable code

 ii=1
 do while(c_string(ii)/=C_NULL_CHAR.and.ii<=len(f_string))
   f_string(ii:ii)=c_string(ii) ; ii=ii+1
 end do
 if (ii<len(f_string)) f_string(ii:)=' '

 end subroutine xc_char_to_f
#endif


!!=================================================================

end  Module libxc_mod
