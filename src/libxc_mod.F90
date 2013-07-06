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

 use globalmath
#if defined HAVE_LIBXC
 use xc_f90_types_m
 use libxc_funcs_m
 use xc_f90_lib_m
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

!!=================================================================
!! STRUCTURED DATATYPES
!!=================================================================

#if defined HAVE_LIBXC
 type libxc_functional
  private
  integer :: family ! LDA, GGA, etc.
  integer :: id     ! identifier
  type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
  type(xc_f90_pointer_t) :: info ! information about the functional
 end type libxc_functional
 type(libxc_functional) :: libxc_funcs(2)
#endif

 private

 public :: libxc_init_func,&
&          libxc_print_func,&
&          libxc_getid_fromname,&
&          libxc_getshortname,&
&          libxc_getid,&
&          libxc_getvxc,&
&          libxc_isgga,&
&          libxc_end

 CONTAINS

#if defined HAVE_LIBXC
#include "libxc_id.in"
#endif

!!=================================================================
!! NAME
!! libxc_getid_fromname
!!
!! FUNCTION
!! From a character string (given in input file), gives the libXC id(s)
!!
!! PARENTS
!! excor
!!
!!=================================================================
 subroutine libxc_getid_fromname(xcname,id,xcname_short)

 implicit none
 integer,intent(inout) :: id(2)
 character*(*),intent(in) :: xcname
 character*(*),intent(out),optional :: xcname_short

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii,i_plus
 character*50 :: xcstrg(2)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 i_plus=index(xcname,'+')
 if (i_plus<=0) then
  xcstrg(1)=trim(xcname)
  xcstrg(2)=""
 else
  xcstrg(1)=trim(xcname(1:i_plus-1))
  xcstrg(2)=trim(xcname(i_plus+1:))
 end if

 do ii=1,2
  id(ii)=0
  call uppercase(xcstrg(ii))

  id(ii)=libxc_id(xcstrg(ii))

  if (id(ii)<0.and.ii==2) then
    id(ii)=0
    exit
  end if

  if (id(ii)==0.and.xcstrg(ii)(1:6)=="LIBXC_") then
   read(unit=xcstrg(ii)(7:),fmt=*,err=333,end=333) id(ii)
333 continue
  end if

  if (id(ii)==0) then
   id(ii)=-1
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

 end subroutine libxc_getid_fromname


!!=================================================================
!! NAME
!! libxc_getid
!!
!! FUNCTION
!! From LibXC datastructure, gives the libXC id(s)
!!
!! PARENTS
!! rdpawps1
!!
!!=================================================================
 subroutine libxc_getid(id)

 implicit none
 integer :: id(2)

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 do ii=1,2
  id(ii)=libxc_funcs(ii)%id
 end do

#else
 id(1:2)=-2
#endif

 end subroutine libxc_getid


!!=================================================================
!! NAME
!! libxc_getshortname
!!
!! FUNCTION
!! From a character string (given in input file), gives a short
!! version of the libXC name (without XC_)
!!
!! PARENTS
!! abinitinterface,xmlinterface
!!
!!=================================================================
 subroutine libxc_getshortname(xcname,xcname_short)

 implicit none
 character*(*),intent(in) :: xcname
 character*(*),intent(out) :: xcname_short

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: i_plus
 character*50 :: xcstrg(2)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

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
!! PARENTS
!! initexch
!!
!!=================================================================
 subroutine libxc_init_func(id,nsp)

 implicit none
 integer,intent(in)  :: id(2)
 integer, intent(in) :: nsp

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 libxc_funcs(1)%id=id(1)
 libxc_funcs(2)%id=id(2)

 do ii=1,2

  if (libxc_funcs(ii)%id==0) then
   libxc_funcs(ii)%family = 0
   cycle
  end if

! Get XC functional family
  libxc_funcs(ii)%family = xc_f90_family_from_id(libxc_funcs(ii)%id)
  select case (libxc_funcs(ii)%family)
   case (XC_FAMILY_LDA, XC_FAMILY_GGA)
    call xc_f90_func_init(libxc_funcs(ii)%conf,libxc_funcs(ii)%info,libxc_funcs(ii)%id,nsp)
   case default
    write(6,'(4a,i3,4a)' ) char(10),&
&    ' Error in libxc_functionals_init:',char(10),&
&    '  The LibXC functional family ',libxc_funcs(ii)%family,char(10),&
&    '  is currently unsupported by AtomPAW !',char(10),&
&    '  (at present only LGA or GGA are supported)'
    stop
  end select

  if (libxc_funcs(ii)%id == XC_LDA_C_XALPHA) then
   call xc_f90_lda_c_xalpha_set_par(libxc_funcs(ii)%conf,0.d0)
  end if

 end do ! loop on functionals

!Print functional(s) information
 call libxc_print_func(6)

#endif
 end subroutine libxc_init_func


!!=================================================================
!! NAME
!! libxc_end
!!
!! FUNCTION
!! Free libXC functional(s)
!!
!! PARENTS
!!
!!=================================================================
 subroutine libxc_end()

  implicit none

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------
 do ii=1,2
  if (libxc_funcs(ii)%id==0) cycle
  call xc_f90_func_end(libxc_funcs(ii)%conf)
 end do

#endif
 end subroutine libxc_end


!!=================================================================
!! NAME
!! libxc_print_func
!!
!! FUNCTION
!! Print libXC functionnal(s) details
!!
!! PARENTS
!! atompaw,libxc_init_func
!!
!!=================================================================
 subroutine libxc_print_func(unt)

 implicit none
 integer :: unt

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 integer :: ii,jj
 character(len=500) :: msg
 type(xc_f90_pointer_t) :: libxc_str

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 do ii=1,2

  if (libxc_funcs(ii)%id==0) cycle

  select case (xc_f90_info_kind(libxc_funcs(ii)%info))
   case (XC_EXCHANGE)
    write(unt,'(a)') 'Exchange functional (LibXC):'
   case (XC_CORRELATION)
    write(unt,'(a)') 'Correlation functional (LibXC):'
   case (XC_EXCHANGE_CORRELATION)
    write(unt,'(a)') 'Exchange-Correlation functional (LibXC):'
  end select

  call xc_f90_info_name(libxc_funcs(ii)%info,msg)
  write(unt,'(2x,a)') trim(msg)
  jj=0
  call xc_f90_info_refs(libxc_funcs(ii)%info,jj,libxc_str,msg)
  do while (jj>=0)
   write(unt,'(2x,a)') trim(msg)
   call xc_f90_info_refs(libxc_funcs(ii)%info,jj,libxc_str,msg)
  end do
 end do

#endif
 end subroutine libxc_print_func


!!=================================================================
!! NAME
!! libxc_isgga
!!
!! FUNCTION
!! Returns TRUE is LibXC functional is GGA
!!
!! PARENTS
!! exch
!!
!!=================================================================
 function libxc_isgga()

 implicit none
 logical :: libxc_isgga

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------
 libxc_isgga = .false.
#if defined HAVE_LIBXC
 libxc_isgga = (any(libxc_funcs(:)%family == XC_FAMILY_GGA))
#endif

 end function libxc_isgga


!!=================================================================
!! NAME
!! libxc_getvxc
!!
!! FUNCTION
!! Returns TRUE is LibXC functional is GGA
!!
!! NOTES
!!  nsp=1 : rho is total density (not half)
!!          grho is abs(grad(rho))
!!  nsp=2 : rho is [rho^up,rho^dn]
!!          grho is [abs(grad(rho^up)),abs(grad(rho^dn)),abs(grad(rho^tot))]
!! PARENTS
!! exch
!!
!!=================================================================
 subroutine libxc_getvxc(npts,exc,vxc,nsp,rho,grho,vxcgr)

 implicit none
 integer, intent(in)            :: npts,nsp
 real(8),intent(in)             :: rho(npts,nsp)
 real(8),intent(inout)          :: exc(npts),vxc(npts,nsp)
 real(8),intent(in),optional    :: grho(npts,2*nsp-1)
 real(8),intent(inout),optional :: vxcgr(npts,2*nsp-1)

#if defined HAVE_LIBXC
!------------------------------------------------------------------
!---- Local variables
!------------------------------------------------------------------

 real(8),parameter :: tol=1.d-14
 integer :: ii,ipts,izero
 real(8) :: rhotmp(nsp),exctmp
 real(8) :: sigma(3),vsigma(3),vxctmp(nsp)

!------------------------------------------------------------------
!---- Executable code
!------------------------------------------------------------------

 if (libxc_isgga().and.((.not.present(vxcgr).or.(.not.present(grho))))) then
   write(6,'(/,2x,a)') "Bug in libxc_getvxc:"
   write(6,'(2x,3a)')  " GGA called without grho or vxcgr !"
   stop
 end if

!Initializations
 vxc=0.d0;exc=0.d0
 if (libxc_isgga()) vxcgr=0.d0

!Filter density/gradient when density goes to zero
 izero=0
 do ipts=1,npts
  if (rho(ipts,1)>tol) izero=ipts
 end do
 if (nsp==2.and.izero<npts) then
  do ipts=izero+1,npts
   if (rho(ipts,2)>tol) izero=ipts
  end do
 end if

!Loop over points
 do ipts=1,npts

  vxctmp=0.d0;exctmp=0.d0
  if (ipts<=izero) then
   rhotmp(1:nsp)=rho(ipts,1:nsp)
  else
   rhotmp=tol
  end if

  if (libxc_isgga()) then
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

! Loop over functionals
  do ii=1,2
   if (libxc_funcs(ii)%id==0) cycle

!  Get the potential (and possibly the energy)
   if (iand(xc_f90_info_flags(libxc_funcs(ii)%info),XC_FLAGS_HAVE_EXC)/=0) then
    select case (libxc_funcs(ii)%family)
     case (XC_FAMILY_LDA)
      call xc_f90_lda_exc_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),exctmp,vxctmp(1))
     case (XC_FAMILY_GGA)
      call xc_f90_gga_exc_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),sigma(1),&
&                             exctmp,vxctmp(1),vsigma(1))
    end select
   else
    exctmp=0.d0
    select case (libxc_funcs(ii)%family)
     case (XC_FAMILY_LDA)
      call xc_f90_lda_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),vxctmp(1))
     case (XC_FAMILY_GGA)
      call xc_f90_gga_vxc(libxc_funcs(ii)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))
    end select
   end if

   exc(ipts)      =exc(ipts)               +2.d0*exctmp           ! From Ha to Ry
   vxc(ipts,1:nsp)=vxc(ipts,1:nsp)         +2.d0*vxctmp(1:nsp)    ! From Ha to Ry
   if (libxc_isgga()) then
    if (nsp==1) then
     vxcgr(ipts,2*nsp-1)=vxcgr(ipts,2*nsp-1)+4.d0*vsigma(2*nsp-1) ! From Ha to Ry
                         ! Note: for nsp=1, vsigma(1) contains 1/2 (up part only)
    else
     vxcgr(ipts,2*nsp-1)=vxcgr(ipts,2*nsp-1)+2.d0*vsigma(2*nsp-1) ! From Ha to Ry
    end if
!   if (nsp==1) then
!    vxcgr(ipts,3)=vxcgr(ipts,3)+vsigma(1)*2.d0
!   else
!    vxcgr(ipts,1)=vxcgr(ipts,1)+2.d0*vsigma(1)-vsigma(2)
!    vxcgr(ipts,2)=vxcgr(ipts,2)+2.d0*vsigma(3)-vsigma(2)
!    vxcgr(ipts,3)=vxcgr(ipts,3)+     vsigma(2)
!   end if
   end if

  end do ! loop over functional(s)
 end do  ! loop over points

#endif
 end subroutine libxc_getvxc

end  Module libxc_mod
