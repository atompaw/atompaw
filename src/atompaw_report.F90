!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains the following active subroutines:
!     Report_Atomres, Report_Pseudobasis, Report_Pseudopotential,
!        WRITE_ATOMDATA, Report_pseudo_energies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE  atompaw_report

  USE io_tools
  USE atomdata
  USE fock
  USE globalmath
  USE gridmod
  USE pseudo
  USE pseudodata
  USE libxc_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Report_Atomres(key,Grid,Orbit,Pot,SCF,ifen)
    CHARACTER(2),INTENT(IN) :: key
    TYPE(GridInfo),INTENT(IN) :: Grid
    TYPE(PotentialInfo),INTENT(IN) :: Pot
    TYPE(OrbitInfo),INTENT(IN) :: Orbit
    TYPE(SCFInfo),INTENT(IN) :: SCF
    INTEGER, INTENT(IN) :: ifen

    INTEGER :: i,j,io,many,l,istart
    CHARACTER (len=4) :: flnm
    CHARACTER (len=20) :: nm
    CHARACTER (len=2) :: sym
    CHARACTER (len=1) :: syml
    REAL(8), POINTER :: r(:),den(:),rv(:),wfn(:,:)
    INTEGER, POINTER :: n,nps,npp,npd,npf,npg

    WRITE(ifen,*) 'Completed calculations for ',TRIM(POT%sym)

    SELECT CASE(TRIM(Orbit%exctype))
    CASE default
      if (have_libxc) then
        WRITE(ifen,*) 'Exchange-correlation type:'
        call libxc_print_func(ifen)
        if(needvtau) then
           WRITE(ifen,*) 'Full generalized Kohn-Sham equations solved'
        else
           WRITE(ifen,*) 'Conventional Kohn-Sham equations solved'
        endif    
      else
        WRITE(ifen,*) 'Exchange-correlation type: LDA, Perdew-Wang correlation'
      end if
    CASE('LDA-PW')
       WRITE(ifen,*) 'Exchange-correlation type: LDA, Perdew-Wang correlation'
    CASE('GGA-PBE')
       WRITE(ifen,*) 'Exchange-correlation type: GGA, Perdew-Burke-Ernzerhof'
    CASE('GGA-PBESOL')
       WRITE(ifen,*) &
&     'Exchange-correlation type: GGA, Perdew-Burke-Ernzerhof modified (PBEsol)'
    CASE('HF')
       WRITE(ifen,*) 'Exchange-correlation type: Hartree-Fock (in devel)'
    CASE('HFV')
       WRITE(ifen,*) &
&    'Exchange-correlation type: Hartree-Fock, frozen-core potential (in devel)'
    CASE('EXX')
       WRITE(ifen,*) 'Exchange-correlation type: Exact-exchange (in devel)'
    CASE('EXXOCC')
       WRITE(ifen,*) 'Exchange-correlation type: Exact-exchange OCC (in devel)'
    CASE('EXXKLI')
       WRITE(ifen,*) 'Exchange-correlation type: Exact-exchange KLI (in devel)'
    CASE('EXXCS')
       WRITE(ifen,*) &
&     'Exchange-correlation type: Exact-exchange Colle-Salvetti (in devel)'
    END SELECT

    CALL reportgrid(Grid,ifen)
    IF (scalarrelativistic) THEN
       if (finitenucleus) then
         WRITE(ifen,*) 'Scalar relativistic calculation -- finite nucleus'
       else
         WRITE(ifen,*) 'Scalar relativistic calculation'
       endif
    ELSEIF (diracrelativistic) then   
         WRITE(ifen,*) 'Dirac-relativistic calculation'
    ELSE
       if (finitenucleus) then
         WRITE(ifen,*) 'Non-relativistic calculation -- finite nucleus'
       else
         WRITE(ifen,*) 'Non-relativistic calculation'
       endif
    ENDIF
    IF (key=='AE') &
         WRITE(ifen,*) '  AEatom converged in',SCF%iter,' iterations'
    IF (key=='SC') &
         WRITE(ifen,*) '  SCatom converged in',SCF%iter,' iterations'
    WRITE(ifen,'(a,f6.2)') '     for nz = ',Pot%zz
    WRITE(ifen,*) '    delta  = ', SCF%delta
    IF (key=='AE') THEN
      WRITE(ifen,*) ' All Electron Orbital energies:         '
      if(diracrelativistic) then
      WRITE(ifen,*) ' n   kappa  l   occupancy       energy'
       DO io=1,Orbit%norbit
          write(std_out,*) 'io',io, Orbit%norbit; call flush_unit(std_out)
          WRITE(ifen,'(i2,1x,i2,2x,i2,4x,1p,2e15.7)') &
&               Orbit%np(io),Orbit%kappa(io),Orbit%l(io),&
&               Orbit%occ(io),Orbit%eig(io)
       ENDDO
      else        
      WRITE(ifen,*) ' n  l     occupancy       energy'
       DO io=1,Orbit%norbit
          write(std_out,*) 'io',io, Orbit%norbit; call flush_unit(std_out)
          WRITE(ifen,'(i2,1x,i2,4x,1p,2e15.7)') &
               Orbit%np(io),Orbit%l(io),&
               Orbit%occ(io),Orbit%eig(io)
       ENDDO
       endif 
    ELSE IF (key=='SC') THEN
      WRITE(ifen,*) '  Valence Electron Orbital energies:         '
      if(diracrelativistic) then
      WRITE(ifen,*) ' n  kappa l   occupancy       energy'
       DO io=1,Orbit%norbit
          IF (.NOT.orbit%iscore(io)) &
&          WRITE(ifen,'(i2,1x,i2,2x,i2,4x,1p,2e15.7)') &
&               Orbit%np(io),Orbit%kappa(io),Orbit%l(io),&
&               Orbit%occ(io),Orbit%eig(io)
       ENDDO
      else        
      WRITE(ifen,*) ' n  l     occupancy       energy'
       DO io=1,Orbit%norbit
          IF (.NOT.orbit%iscore(io))WRITE(ifen,'(i2,1x,i2,4x,1p,2e15.7)') &
               Orbit%np(io),Orbit%l(io),&
               Orbit%occ(io),Orbit%eig(io)
       ENDDO
       endif 
    ENDIF
    WRITE(ifen,*)
    WRITE(ifen,*) ' Total energy'
    WRITE(ifen,*) '    Total                    :  ',SCF%etot
    IF (key=='SC') THEN
       WRITE(ifen,*) '    Valence                  :  ',SCF%evale
    ENDIF

  END SUBROUTINE Report_Atomres


  SUBROUTINE Report_Pseudobasis(Grid,PAW,ifen)
     Type(GridInfo), INTENT(IN)  :: Grid
     Type(PseudoInfo), INTENT(IN) :: PAW
     INTEGER, INTENT(IN) :: ifen

     INTEGER, parameter :: ifout=15
     INTEGER :: i,j,io,nbase,irc,icount,n
     INTEGER, ALLOCATABLE :: mapp(:)
     CHARACTER (len=4) :: flnm

     nbase=PAW%nbase;irc=PAW%irc;n=Grid%n
     WRITE(ifen,'(/"Number of basis functions ",i5)') nbase
     WRITE(ifen,*)'No.   n    l      Energy         Cp coeff         Occ'


     DO io=1,nbase
        WRITE(ifen,'(3i5,1p,3e15.7)') io,PAW%np(io),PAW%l(io),PAW%eig(io),&
          PAW%ck(io),PAW%occ(io)
        CALL mkname(io,flnm)
        OPEN(ifout,file='wfn'//TRIM(flnm),form='formatted')
        WRITE(ifout,*) '# l=',PAW%l(io),'basis function with energy  ',&
             PAW%eig(io)
          DO i=1,irc+50
             WRITE(ifout,'(1p,5e12.4)') Grid%r(i),PAW%ophi(i,io),&
                    PAW%otphi(i,io),PAW%otp(i,io)
          ENDDO
       CLOSE(ifout)
    ENDDO

    ! also write "raw" wavefunctions
     DO io=1,nbase
        CALL mkname(io,flnm)
        OPEN(ifout,file='wfn00'//TRIM(flnm),form='formatted')
        WRITE(ifout,*) '# l=',PAW%l(io),'basis function with energy  ',&
             PAW%eig(io)
          DO i=1,irc+50
             WRITE(ifout,'(1p,5e12.4)') Grid%r(i),PAW%phi(i,io),&
                    PAW%tphi(i,io),PAW%tp(i,io)
          ENDDO
       CLOSE(ifout)
    ENDDO

    allocate(mapp(PAW%OCCWFN%norbit))
    mapp=0
    icount=0
    do io=1,PAW%OCCWFN%norbit
       if(PAW%OCCWFN%iscore(io)) then
       else
         icount=icount+1
         mapp(icount)=io
       endif
    enddo
    OPEN(ifout,file='OCCWFN',form='formatted')
    WRITE(ifout,'("#            ",50i30)') (mapp(j),j=1,icount)
    do i=1,n
       write(ifout,'(1p,51e15.7)') Grid%r(i),(PAW%OCCWFN%wfn(i,mapp(j)),&
                PAW%TOCCWFN%wfn(i,mapp(j)),j=1,icount)
    enddo
    close(ifout)

    deallocate(mapp)

  END SUBROUTINE Report_Pseudobasis

  SUBROUTINE Report_Pseudopotential(Grid,PAW)
     Type(GridInfo), INTENT(IN)  :: Grid
     Type(PseudoInfo), INTENT(INOUT) :: PAW

     INTEGER, parameter :: ifout=15
     INTEGER :: i,n,irc
     REAL(8) :: sqr4pi

     n=Grid%n; irc=PAW%irc

     OPEN(ifout,file='density', form='formatted')
     WRITE(ifout,*) '#  r       coreden        valeden         tcoreden        tvaleden      nhatden'
     DO i=1,n
        IF (abs(PAW%core(i))<machine_zero) PAW%core(i)=0
        IF (abs(PAW%tcore(i))<machine_zero) PAW%tcore(i)=0
        IF (abs(PAW%den(i))<machine_zero) PAW%den(i)=0
        IF (abs(PAW%tden(i))<machine_zero) PAW%tden(i)=0
        IF (abs(PAW%nhatv(i))<machine_zero) PAW%nhatv(i)=0
        WRITE(ifout,'(1p,1e15.7,1p,5e25.17)') Grid%r(i),PAW%core(i),&
             PAW%den(i),PAW%tcore(i),PAW%tden(i),PAW%nhatv(i)
     ENDDO
     CLOSE(ifout)

     OPEN(ifout,file='tau', form='formatted')
     WRITE(ifout,*) '#  r       coretau        valetau         tcoretau        tvaletau '
     DO i=1,n
        IF (PAW%coretau(i)<machine_zero) PAW%coretau(i)=0
        IF (PAW%tcoretau(i)<machine_zero) PAW%tcoretau(i)=0
        IF (PAW%valetau(i)<machine_zero) PAW%valetau(i)=0
        IF (PAW%tvaletau(i)<machine_zero) PAW%tvaletau(i)=0
        WRITE(ifout,'(1p,1e15.7,1p,5e25.17)') Grid%r(i),PAW%coretau(i),&
&           PAW%valetau(i),PAW%tcoretau(i),PAW%tvaletau(i)
     ENDDO
     CLOSE(ifout)

     OPEN(ifout,file='potential', form='formatted')
     DO i=1,n
        IF (ABS(PAW%AErefrv(i))<machine_zero) PAW%AErefrv(i)=0
        IF (ABS(PAW%rveff(i))<machine_zero) PAW%rveff(i)=0
        WRITE(ifout,'(1p,1e15.7,1p,3e25.17)')Grid%r(i),PAW%AErefrv(i),PAW%rveff(i)
     ENDDO
     CLOSE(ifout)

     OPEN(ifout,file='vloc', form='formatted')
     DO i=1,irc+10
        WRITE(ifout,'(1p,1e15.7,1p,3e25.17)') Grid%r(i),PAW%vloc(i)
     ENDDO
     CLOSE(ifout)

     If (TRIM(PAW%exctype)/='HF') then
     OPEN(ifout,file='rVx', form='formatted')
     DO i=1,n
        WRITE(ifout,'(1p,4e15.7)') Grid%r(i),PAW%rvx(i),PAW%trvx(i)
     ENDDO
     CLOSE(ifout)
     EndIf

!   Find radii ensuring:
!     - Abs(density)<10e-10 (for the pseudo-valence density)
!     - r>=10 bohr (for the ionic potential)
    sqr4pi=sqrt(4*pi)*1.d-10;PAW%ivion=gridindex(Grid,10.d0)
    PAW%ivale=PAW%ivion
    do while (PAW%ivale<Grid%n.and. &
&        abs(PAW%tden(PAW%ivale))>sqr4pi*Grid%r(PAW%ivale)**2)
     PAW%ivale=PAW%ivale+1
    end do

  End SUBROUTINE Report_Pseudopotential

  SUBROUTINE WRITE_ATOMDATA(Grid,Pot,Orbit,FC,PAW)
    TYPE(GridInfo) , INTENT(IN):: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    TYPE(FCInfo), INTENT(IN) :: FC
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW

    integer, parameter :: ifatompaw=1001
    integer :: i,j,k,l,ib,ic,io,jo,li,lj,n,ishift,ivion,ivale,lcao_points,icount
    integer :: lmin,lmax,id,ie,lp
    integer, allocatable :: llist(:)
    real(8) :: sqr4pi, term,accum,occ
    real(8), allocatable :: f(:)
    logical :: even,print_pwpawfile

	write(std_out,'(/,a)') '========================================================'
	write(std_out,'(a)')   '==   Atompaw2PWPAW                                    =='
	write(std_out,'(a)')   '==   Writes a PAW dataset generated by "Atompaw"      =='
	write(std_out,'(a)')   '==   for PWPAW software...                            =='
	write(std_out,'(a)')   '========================================================'

	print_pwpawfile=.true.

   !Do not try to print UPF dataset when it is not possible
	if (diracrelativistic) then
	  write(std_out,'(/,2x,a)') 'PWPAW dataset output not compatible with Dirac relativistic approach!'
	  write(std_out,'(2x,a)')   'PWPAW dataset file will not be created.'
	  print_pwpawfile=.false.;return
	end if

    ishift=Grid%ishift

    PAW%mesh_size=PAW%irc+ishift
    PAW%coretailpoints=MAX(PAW%coretailpoints,PAW%mesh_size)

!   Find radii ensuring:
!     - Abs(density)<10e-10 (for the pseudo-valence density)
!     - r>=10 bohr (for the ionic potential)
    sqr4pi=sqrt(4*pi)*1.d-10;ivion=gridindex(Grid,10.d0);ivale=ivion
    do while (ivale<Grid%n.and.abs(PAW%tden(ivale))>sqr4pi*Grid%r(ivale)**2)
      ivale=ivale+1
    end do

    allocate(f(Grid%n))

    OPEN(ifatompaw,file=TRIM(POT%sym)//'.atomicdata',form='formatted')
    WRITE(ifatompaw,'("  ATOMTYPE     ",a2)') POT%sym
    WRITE(ifatompaw,'("  ATOMXCTYPE     ",a10)') TRIM(Orbit%exctype)
    WRITE(ifatompaw,'("  ATOMIC_CHARGE    ",f5.0)') POT%zz
    WRITE(ifatompaw,'("  MOMENTLESSHARTREE    ")')   ! new grouping of terms
    WRITE(ifatompaw,'("  CORE_CHARGE    ",1p,1e20.13)') FC%zcore
    WRITE(ifatompaw,'("  RC         ",1p,1e20.13)') PAW%rc
        if (gaussianshapefunction) then
     WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1p,1e20.13)') 'gaussian', &
           PAW%gausslength
    else if (besselshapefunction) then
     if (PAW%multi_rc) then
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1p,1e20.13)') 'bessel',PAW%rc_shap
     else
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20)') 'bessel'
     endif
    else
     if (PAW%multi_rc) then
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1p,1e20.13)') 'sinc2',PAW%rc_shap
     else
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20)') 'sinc2'
     endif
    endif
    WRITE(ifatompaw,'("  BASIS_SIZE    ",i5)') PAW%nbase
    WRITE(ifatompaw,'("  ORBITALS      ")')
    WRITE(ifatompaw,'(20i4)') (PAW%l(ib),ib=1,PAW%nbase)
    WRITE(ifatompaw,'("  END     ")')
    WRITE(ifatompaw,'("  INITOCC       ")')
    WRITE(ifatompaw,'(8f10.6)') (PAW%occ(ib),ib=1,PAW%nbase)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("  MESH_SIZE    ",i10)') PAW%mesh_size
    WRITE(ifatompaw,'("  MESH_STEP    ",1p,1e20.13)') Grid%h
    if (usingloggrid(Grid)) then
       WRITE(ifatompaw,'("  LOG_GRID    ",1p,1e20.13)') Grid%drdu(1)
    endif

    WRITE(ifatompaw,'("  CORETAIL_POINTS   ",i10)') PAW%coretailpoints
!   Find index for Grid%r(i)>10
    j=gridindex(Grid,10.d0); lcao_points=j
    WRITE(ifatompaw,'("  LCAO_SIZE  ",i10)') j     !lcao_points
    WRITE(ifatompaw,'("  LCAO_STEP   ",1p,1e20.13)') Grid%h     !hlcao

    WRITE(ifatompaw,'("   CORE_DENSITY   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (FC%coreden(i),i=1,PAW%mesh_size)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   CORETAIL_DENSITY   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%tcore(i),i=1,PAW%coretailpoints)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   PSEUDO_VALENCE_DENSITY   ",3x,i8)') ivale
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%tden(i),i=1,ivale)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   SHAPE_FUNC   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%hatshape(i),i=1,PAW%mesh_size)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   VLOCFUN      ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%vloc(i),i=1,PAW%mesh_size)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,&
     '("   VLOCION      ",3x,i8," #ionic vloc for abinit in Ryd units")') ivion
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%abinitvloc(i),i=1,ivion)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,&
     '("   VLOCION_NOHAT",3x,i8," #ionic vlocnohat for abinit in Ryd units")') ivion
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%abinitnohat(i),i=1,ivion)
    WRITE(ifatompaw,'("  END     ")')

    DO ib=1,PAW%nbase
       WRITE(ifatompaw,'("   TPROJECTOR",i4," #p(r), for p(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%otp(i,ib),i=1,PAW%mesh_size)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO


    DO ib=1,PAW%nbase
       WRITE(ifatompaw,'("   PHI",i5," #phi(r), for phi(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%ophi(i,ib),i=1,PAW%mesh_size)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO

    DO ib=1,PAW%nbase
       WRITE(ifatompaw,'("   TPHI",i5," #tphi(r), for tphi(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%otphi(i,ib),i=1,PAW%mesh_size)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO

    DO ib=1,PAW%nbase
       f=PAW%tphi(:,ib)
       CALL trunk(Grid,f(1:Grid%n),6.d0,10.d0)
       WRITE(ifatompaw,'("   TPHI_LCAO",i4," #tphi0(r) for tphi0(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (f(j),j=1,lcao_points)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO


!! optional l-dependent shape functions
    if (besselshapefunction) then
      li=MAXVAL(PAW%TOCCWFN%l(:)); li=MAX(li,PAW%lmax); li=2*li
      WRITE(ifatompaw,'("    SHAPE_L",i4,"  # for  l= 0 .. lmax")') li
      Do l=1,li+1
         f(2:PAW%mesh_size)=PAW%g(2:PAW%mesh_size,l)/(Grid%r(2:PAW%mesh_size)**2)
         call extrapolate(Grid,f)
         WRITE(ifatompaw,'(1p,3e25.17)') (f(i),i=1,PAW%mesh_size)
         WRITE(ifatompaw,'( "END")')
      ENddo
    endif



    ! spherical matrix elements
    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          IF (PAW%l(ib)==PAW%l(ic)) THEN
             icount=icount+1
             f(icount)=PAW%oij(ib,ic)
          ENDIF
       ENDDO
    ENDDO


    WRITE(ifatompaw,'("   OVERLAP_SIZE    ",i10)') icount
    WRITE(ifatompaw,'("   OVERLAP_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (f(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')


    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          IF (PAW%l(ib)==PAW%l(ic)) THEN
             icount=icount+1
             f(icount)=PAW%Kij(ib,ic)
          ENDIF
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   KINETIC_ENERGY_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (f(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')


    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          IF (PAW%l(ib)==PAW%l(ic)) THEN
             icount=icount+1
             f(icount)=PAW%VFij(ib,ic)
          ENDIF
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   V_ION_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (f(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')

    !
    !  angularly dependent matrix elements
    !
    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO l=lmin,lmax,2
             icount=icount+1
          ENDDO
       ENDDO
    ENDDO

    WRITE(ifatompaw,'("   DENVHAT_SIZE   ",i10)') icount
    WRITE(ifatompaw,'("   DENSITY   ")')
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO l=lmin,lmax,2
             WRITE(ifatompaw,'(3i10,1p,1e25.17)') ib,ic,l,PAW%mLij(ib,ic,l+1)
          ENDDO
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   END   ")')

    !
    !  New (momentless) form of Hartree matrix elements
    !

    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO l=lmin,lmax,2
             DO id=1,PAW%nbase
                DO ie=id,PAW%nbase
                   lp=lmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (l.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
                        l.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      icount=icount+1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    WRITE(ifatompaw,'("   HARTREE_SIZE   ",i10)') icount
    WRITE(ifatompaw,'("   V_HARTREE   ")')

    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO l=lmin,lmax,2
             DO id=1,PAW%nbase
                DO ie=id,PAW%nbase
                   lp=lmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (l.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
                        l.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                    WRITE(ifatompaw,'(5i5,1p,1e25.17)')ib,ic,id,ie,l,&
                           PAW%DR(ib,ic,id,ie,l+1)
                      icount=icount+1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    WRITE(ifatompaw,'("   END         ")')

    WRITE(ifatompaw,'("   EAION   ",1p,1e25.17)') PAW%Eaion
    WRITE(ifatompaw,'("   EAIONHAT   ",1p,1e25.17)') PAW%Eaionhat
    WRITE(ifatompaw,'("   ENERGY   ",1p,1e25.17)') PAW%Etotal

  ! Core-valence exchange energies now generally provided for possible
  !    use in exact exchange and hybrid evaluations
    Write(STD_OUT,'(3(/,2x,a))') 'Atompaw2PWPAW info:',&
&             '  Core-valence exchange integrals provided',&
&             '  with the assumption that all core states confined 0 < r < rc'

    WRITE(ifatompaw,'("   XCOREVAL  ")')
      do ib=1,PAW%nbase
         do ic=1,PAW%nbase
            WRITE(ifatompaw, '(2i10, 1p,1e25.17)' ) ib,ic, PAW%TXVC(ib,ic)
         enddo
      enddo   
   
    WRITE(ifatompaw,'("   XCORECORE   ",1p, 1e25.17)') PAW%XCORECORE   


    IF (PAW%OCCWFN%exctype=='HF'.or.PAW%OCCWFN%exctype=='EXXKLI') THEN
       Write(STD_OUT,*) 'For HF/KLI -- additional core information provided'
       Write(STD_OUT,*) 'Warning ... assume all core states confined 0 < r < rc'

       icount=0
       do io=1,PAW%OCCWFN%norbit
          if (PAW%OCCWFN%iscore(io)) icount=icount+1
       enddo

       WRITE(ifatompaw,'(" CORE_SIZE  ",i10)') icount
       if (icount>0) then
          allocate(llist(icount))
          i=0
          do io=1,PAW%OCCWFN%norbit
             if (PAW%OCCWFN%iscore(io)) then
                i=i+1
                llist(i)=PAW%OCCWFN%l(io)
             endif
          enddo
          WRITE(ifatompaw,'("  CORE_L      ")')
          WRITE(ifatompaw,'(20i4)') (llist(i),i=1,icount)
          WRITE(ifatompaw,'("  END     ")')

          deallocate(llist)

          i=0
          do io=1,PAW%OCCWFN%norbit
             if (PAW%OCCWFN%iscore(io)) then
                i=i+1
                WRITE(ifatompaw,'("  CORE_PSI  ",i5)') i
                WRITE(ifatompaw,'(1p,3e25.17)') &
                  (PAW%OCCWFN%wfn(j,io), j=1,PAW%coretailpoints)
                WRITE(ifatompaw,'("  END     ")')
             endif
          enddo

          f=0;i=0       !core-core terms
          do io=1,PAW%OCCWFN%norbit
             if (PAW%OCCWFN%iscore(io)) then
                i=i+1; li=PAW%OCCWFN%l(io)
                do  jo=1,PAW%OCCWFN%norbit
                   if (PAW%OCCWFN%iscore(jo)) then
                      occ=PAW%OCCWFN%occ(jo);lj=PAW%OCCWFN%l(jo)
                   do l=abs(li-lj),(li+lj),2
                      call EXXwgt(1.d0,occ,1,li,2,lj,l,accum)
                      call CondonShortley(Grid,l,PAW%OCCWFN%wfn(:,io), &
                         PAW%OCCWFN%wfn(:,jo),PAW%OCCWFN%wfn(:,io), &
                         PAW%OCCWFN%wfn(:,jo),term)
                      f(i)=f(i)+accum*term
                      write(std_out,*) 'core-core CondonShortley', i,accum,term
                   enddo
                   endif
                enddo
             endif
           enddo
           WRITE(ifatompaw,'("  CORECORE_R   ")')
           WRITE(ifatompaw,'(1p,3e25.17)') (f(i),i=1,icount)
           WRITE(ifatompaw,'("  END     ")')

!          k=0            ! core-valence terms
!          do io=1,PAW%OCCWFN%norbit
!             if (PAW%OCCWFN%iscore(io)) then
!                li=PAW%OCCWFN%l(io)
!                do ib=1,PAW%nbase
!                   do ic=1,PAW%nbase
!                      if (PAW%l(ib)==PAW%l(ic)) then
!                         k=k+1
!                      endif
!                   enddo
!                enddo
!              endif
!           enddo
!          WRITE(ifatompaw,'("  COREVAL_LIST   ",i10)') k
!          WRITE(ifatompaw,'("  COREVAL_R      ")')
!          i=0;f=0
!          do io=1,PAW%OCCWFN%norbit
!             if (PAW%OCCWFN%iscore(io)) then
!                i=i+1
!                li=PAW%OCCWFN%l(io)
!                do ib=1,PAW%nbase
!                   lj=PAW%l(ib)
!                   do ic=1,PAW%nbase
!                      if (PAW%l(ib)==PAW%l(ic)) then
!                         f(i)=0
!                         do l=abs(li-lj),(li+lj),2
!                           call EXXwgt(1.d0,1.d0,1,li,2,lj,l,accum)
!                           call CondonShortley(Grid,l,PAW%OCCWFN%wfn(:,io), &
!                              PAW%ophi(:,ib),PAW%OCCWFN%wfn(:,io), &
!                              PAW%ophi(:,ic),term)
!                              f(i)=f(i)+2*accum*term !EXXwgt returns 1/2*3J
!                              write(std_out,*) 'core-val CondonShortley',&
!                                    i,li,lj,l,2*accum,term
!                         enddo
!                        WRITE(ifatompaw,'(3i10,1p,1e25.17)') i, ib,ic,f(i)
!                      endif
!                   enddo
!                enddo
!              endif
!           enddo
          WRITE(ifatompaw,'("  COREVAL_LIST   ",i10)') PAW%ncoreshell
          If (PAW%ncoreshell>0) then
             WRITE(ifatompaw,'("  COREVAL_R      ")')
             Do k=1,PAW%ncoreshell
                do ib=1,PAW%nbase
                   do ic=1,PAW%nbase
                      if (ABS(PAW%DRVC(k,ib,ic))>1.d-8) then
                        WRITE(ifatompaw,'(3i10,1p,1e25.17)') k, ib,ic,&
                                   PAW%DRVC(k,ib,ic)
                      endif
                   enddo
                enddo
             Enddo
          WRITE(ifatompaw,'("  END     ")')
          ENdif

    ENDIF
    ENDIF
    WRITE(ifatompaw,'("  END     ")')

    CLOSE(ifatompaw)

    DEALLOCATE(f)
  END SUBROUTINE WRITE_ATOMDATA

  SUBROUTINE WRITE_SOCORROATOMDATA(Grid,Pot,Orbit,FC,PAW)
    TYPE(GridInfo) , INTENT(IN):: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(OrbitInfo), INTENT(IN) :: Orbit
    TYPE(FCInfo), INTENT(IN) :: FC
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW

    integer, parameter :: ifatompaw=1001
    integer :: i,j,k,l,ib,ic,io,jo,li,lj,n,ishift,ivion,ivale,lcao_points,icount
    integer :: lmin,lmax,id,ie,lp,lcount,jcount,ll,ok
    integer, allocatable :: llist(:)
    real(8) :: sqr4pi, term,accum,occ,ctctse,cthatse,rr,eself
    real(8), allocatable :: f(:),dum(:),dum1(:)
    real(8), allocatable :: wf(:),twf(:),rh(:),rth(:),shartree(:),sshartree(:,:,:)
    logical :: even,print_socorrofile

	write(std_out,'(/,a)') '========================================================'
	write(std_out,'(a)')   '==   Atompaw2Socorro                                  =='
	write(std_out,'(a)')   '==   Writes a PAW dataset generated by "Atompaw"      =='
	write(std_out,'(a)')   '==   for SOCORRO software...                          =='
	write(std_out,'(a)')   '========================================================'

	print_socorrofile=.true.

   !Do not try to print UPF dataset when it is not possible
	if (diracrelativistic) then
	  write(std_out,'(/,2x,a)') 'SOCORRO dataset output not compatible with Dirac relativistic approach!'
	  write(std_out,'(2x,a)')   'SOCORRO dataset file will not be created.'
	  print_socorrofile=.false.;return
	end if

    if (besselshapefunction) then
      write(std_out,'(/,2x,a)') 'SOCORRO does not yet support l-dependent shape func!'
	  write(std_out,'(2x,a)')   'SOCORRO dataset file will not be created.'
	  print_socorrofile=.false.;return
	end if

    ishift=Grid%ishift

    PAW%mesh_size=PAW%irc+ishift
    PAW%coretailpoints=MAX(PAW%coretailpoints,PAW%mesh_size)

!   Find radii ensuring:
!     - Abs(density)<10e-10 (for the pseudo-valence density)
!     - r>=10 bohr (for the ionic potential)
    sqr4pi=sqrt(4*pi)*1.d-10;ivion=gridindex(Grid,10.d0);ivale=ivion
    do while (ivale<Grid%n.and.abs(PAW%tden(ivale))>sqr4pi*Grid%r(ivale)**2)
      ivale=ivale+1
    end do

    allocate(f(Grid%n))

    OPEN(ifatompaw,file=TRIM(POT%sym)//'.SOCORRO.atomicdata',form='formatted')
    WRITE(ifatompaw,'("  ATOMTYPE     ",a2)') POT%sym
    WRITE(ifatompaw,'("  ATOMXCTYPE     ",a10)') TRIM(Orbit%exctype)
    WRITE(ifatompaw,'("  ATOMIC_CHARGE    ",i5)') POT%nz
    WRITE(ifatompaw,'("  CORE_CHARGE    ",1p,1e20.13)') FC%zcore
    WRITE(ifatompaw,'("  RC         ",1p,1e20.13)') PAW%rc
        if (gaussianshapefunction) then
     WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1p,1e20.13)') 'gaussian', &
           PAW%gausslength
    else if (besselshapefunction) then
     if (PAW%multi_rc) then
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1p,1e20.13)') 'bessel',PAW%rc_shap
     else
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20)') 'bessel'
     endif
    else
     if (PAW%multi_rc) then
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20,2x,1p,1e20.13)') 'sinc2',PAW%rc_shap
     else
      WRITE(ifatompaw,'("  SHAPE_TYPE  ",a20)') 'sinc2'
     endif
    endif
    WRITE(ifatompaw,'("  BASIS_SIZE    ",i5)') PAW%nbase
    WRITE(ifatompaw,'("  ORBITALS      ")')
    WRITE(ifatompaw,'(20i4)') (PAW%l(ib),ib=1,PAW%nbase)
    WRITE(ifatompaw,'("  END     ")')
    WRITE(ifatompaw,'("  INITOCC       ")')
    WRITE(ifatompaw,'(8f10.6)') (PAW%occ(ib),ib=1,PAW%nbase)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("  MESH_SIZE    ",i10)') PAW%mesh_size
    WRITE(ifatompaw,'("  MESH_STEP    ",1p,1e20.13)') Grid%h
    if (usingloggrid(Grid)) then
       WRITE(ifatompaw,'("  LOG_GRID    ",1p,1e20.13)') Grid%drdu(1)
    endif

    WRITE(ifatompaw,'("  CORETAIL_POINTS   ",i10)') PAW%coretailpoints
!   Find index for Grid%r(i)>10
    j=gridindex(Grid,10.d0); lcao_points=j
    WRITE(ifatompaw,'("  LCAO_SIZE  ",i10)') j     !lcao_points
    WRITE(ifatompaw,'("  LCAO_STEP   ",1p,1e20.13)') Grid%h     !hlcao

    WRITE(ifatompaw,'("   CORE_DENSITY   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (FC%coreden(i),i=1,PAW%mesh_size)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   CORETAIL_DENSITY   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%tcore(i),i=1,PAW%coretailpoints)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   PSEUDO_VALENCE_DENSITY   ",3x,i8)') ivale
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%tden(i),i=1,ivale)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   SHAPE_FUNC   ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%hatshape(i),i=1,PAW%mesh_size)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,'("   VLOCFUN      ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%vloc(i),i=1,PAW%mesh_size)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,&
     '("   VLOCION      ",3x,i8," #ionic vloc for abinit in Ryd units")') ivion
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%abinitvloc(i),i=1,ivion)
    WRITE(ifatompaw,'("  END     ")')

    WRITE(ifatompaw,&
     '("   VLOCION_NOHAT",3x,i8," #ionic vlocnohat for abinit in Ryd units")') ivion
    WRITE(ifatompaw,'(1p,3e25.17)') (PAW%abinitnohat(i),i=1,ivion)
    WRITE(ifatompaw,'("  END     ")')

    DO ib=1,PAW%nbase
       WRITE(ifatompaw,'("   TPROJECTOR",i4," #p(r), for p(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%otp(i,ib),i=1,PAW%mesh_size)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO


    DO ib=1,PAW%nbase
       WRITE(ifatompaw,'("   PHI",i5," #phi(r), for phi(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%ophi(i,ib),i=1,PAW%mesh_size)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO

    DO ib=1,PAW%nbase
       WRITE(ifatompaw,'("   TPHI",i5," #tphi(r), for tphi(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (PAW%otphi(i,ib),i=1,PAW%mesh_size)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO

    DO ib=1,PAW%nbase
       f=PAW%tphi(:,ib)
       CALL trunk(Grid,f(1:Grid%n),6.d0,10.d0)
       WRITE(ifatompaw,'("   TPHI_LCAO",i4," #tphi0(r) for tphi0(r)/r*Ylm)")') ib
       WRITE(ifatompaw,'(1p,3e25.17)') (f(j),j=1,lcao_points)
       WRITE(ifatompaw,'("  END     ")')
    ENDDO


!! optional l-dependent shape functions
!    if (besselshapefunction) then
!      li=MAXVAL(PAW%TOCCWFN%l(:)); li=MAX(li,PAW%lmax); li=2*li
!      WRITE(ifatompaw,'("    SHAPE_L",i4,"  # for  l= 0 .. lmax")') li
!      Do l=1,li+1
!         f(2:PAW%mesh_size)=PAW%g(2:PAW%mesh_size,l)/(Grid%r(2:PAW%mesh_size)**2)
!         call extrapolate(Grid,f)
!         WRITE(ifatompaw,'(1p,3e25.17)') (f(i),i=1,PAW%mesh_size)
!         WRITE(ifatompaw,'( "END")')
!      ENddo
!    endif

    ! spherical matrix elements
    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          IF (PAW%l(ib)==PAW%l(ic)) THEN
             icount=icount+1
             f(icount)=PAW%oij(ib,ic)
          ENDIF
       ENDDO
    ENDDO


    WRITE(ifatompaw,'("   OVERLAP_SIZE    ",i10)') icount
    WRITE(ifatompaw,'("   OVERLAP_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (f(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')


    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          IF (PAW%l(ib)==PAW%l(ic)) THEN
             icount=icount+1
             f(icount)=PAW%Kij(ib,ic)
          ENDIF
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   KINETIC_ENERGY_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (f(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')


    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          IF (PAW%l(ib)==PAW%l(ic)) THEN
             icount=icount+1
             call dvij(Grid,PAW,FC,Pot%zz,ib,ic,f(icount))
          ENDIF
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   V_ION_MATRIX  ")')
    WRITE(ifatompaw,'(1p,3e25.17)') (f(ic),ic=1,icount)
    WRITE(ifatompaw,'("  END     ")')

    !
    !  angularly dependent matrix elements
    !
    icount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO l=lmin,lmax,2
             icount=icount+1
          ENDDO
       ENDDO
    ENDDO

    WRITE(ifatompaw,'("   DENVHAT_SIZE   ",i10)') icount
    WRITE(ifatompaw,'("   DENSITY   ")')
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO l=lmin,lmax,2
             WRITE(ifatompaw,'(3i10,1p,1e25.17)') ib,ic,l,PAW%mLij(ib,ic,l+1)
          ENDDO
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   END   ")')

! Old VHAT matrix elements
    WRITE(ifatompaw,'("   V_HAT   ")')
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO ll=lmin,lmax,2
             CALL hatpotL(Grid,PAW,ll,f)
             DO i=1,PAW%irc
                f(i)=f(i)*PAW%otphi(i,ib)*PAW%otphi(i,ic)
             ENDDO
             WRITE(ifatompaw,'(3i10,1pe25.17)') ib,ic,ll,&
&             integrator(Grid,f,1,PAW%irc)
          ENDDO
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   END   ")')

    !
    !  Old form of Hartree matrix elements
    !

    lcount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO l=lmin,lmax,2
             DO id=1,PAW%nbase
                DO ie=id,PAW%nbase
                   lp=lmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (l.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
                        l.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    WRITE(ifatompaw,'("   HARTREE_SIZE   ",i10)') lcount
    WRITE(ifatompaw,'("   V_HARTREE   ")')

    icount=(PAW%nbase*(PAW%nbase+1))/2; n=Grid%n
    ALLOCATE(shartree(lcount),sshartree(icount,icount,2*PAW%lmax+2),&
&        wf(n),twf(n),rh(n),rth(n),dum(n),dum1(n),stat=ok)
    IF (ok/=0) THEN
      WRITE(STD_OUT,'(/,2x,a)') 'Atompaw2Socorro: error in hartree allocation!'
      WRITE(STD_OUT,'(2x,a,4i3)') '  icount,lcount,irc,ok=',icount,lcount,PAW%irc,ok
      GOTO 111
    ENDIF
    !
    !  Hartree matrix elements
    !
    !! Due to precision error in apoisson solver, some equivalent Hartree
    !!  matrix elements are not equal with an error of e-5; take average
    !!  these elements in order to remove assymmetry errors
    !!  Thanks to Francois Jollet for pointing out this problem
    sshartree=0
    shartree=0
    icount=0
    lcount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          icount=icount+1
          DO i=1,PAW%irc
             wf(i)=PAW%ophi(i,ib)*PAW%ophi(i,ic)
             twf(i)=PAW%otphi(i,ib)*PAW%otphi(i,ic)
          ENDDO
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO ll=lmin,lmax,2
             CALL apoisson(Grid,ll,PAW%irc,wf,rh)
             CALL apoisson(Grid,ll,PAW%irc,twf,rth)
             jcount=0
             DO id=1,PAW%nbase
                DO ie=id,PAW%nbase
                   jcount=jcount+1
                   lp=lmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (ll.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
&                       ll.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                      dum=0;dum1=0
                      DO i=2,PAW%irc
                         rr=Grid%r(i)
                         dum(i)=PAW%ophi(i,id)*PAW%ophi(i,ie)*rh(i)/rr     &
&                             -PAW%otphi(i,id)*PAW%otphi(i,ie)*rth(i)/rr
                         dum1(i)=PAW%ophi(i,id)*PAW%ophi(i,ie)*rh(i)/rr
                      ENDDO
                      shartree(lcount)=integrator(Grid,dum(1:PAW%irc),1,PAW%irc)
                      sshartree(icount,jcount,ll+1)=shartree(lcount)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    icount=0
    lcount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          icount=icount+1
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO ll=lmin,lmax,2
             jcount=0
             DO id=1,PAW%nbase
                DO ie=id,PAW%nbase
                   jcount=jcount+1
                   lp=lmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (ll.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
&                       ll.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                      shartree(lcount)=                      &
&                          0.5d0*(sshartree(icount,jcount,ll+1) + &
&                          sshartree(jcount,icount,ll+1))
                       !write(std_out,*) ll,lcount,shartree(lcount); call flush_unit(std_out)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    lcount=0
    DO ib=1,PAW%nbase
       DO ic=ib,PAW%nbase
          lmin=ABS(PAW%l(ib)-PAW%l(ic))
          lmax=PAW%l(ib)+PAW%l(ic)
          DO ll=lmin,lmax,2
             DO id=1,PAW%nbase
                DO ie=id,PAW%nbase
                   lp=lmax+PAW%l(id)+PAW%l(ie)
                   even=.false.
                   if (2*(lp/2)==lp) even=.true.
                   IF (ll.GE.ABS(PAW%l(id)-PAW%l(ie)).AND.           &
&                       ll.LE.PAW%l(id)+PAW%l(ie).AND.even) THEN
                      lcount=lcount+1
                      WRITE(ifatompaw,'(5i5,1pe25.17)')ib,ic,id,ie,ll, &
&                        shartree(lcount)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    WRITE(ifatompaw,'("   END         ")')

    WRITE(ifatompaw,'("   HAT_SELF-ENERGY    ",i10)') 2*PAW%lmax
    DO ll=0,2*PAW%lmax
       CALL selfhatpot(Grid,PAW,ll,eself)
       WRITE(ifatompaw,'(i10,1pe25.17)') ll,eself
    ENDDO
    WRITE(ifatompaw,'("   END         ")')

    call coretailselfenergy(Grid,PAW,ctctse,cthatse)
    WRITE(ifatompaw,'("   CORETAILSELFENERGY   ",1pe25.17,"  END")') ctctse
    WRITE(ifatompaw,'("   CORETAILHATENERGY   ",1pe25.17,"  END")') cthatse

    WRITE(ifatompaw,'("   ENERGY   ",1p,1e25.17)') PAW%Etotal

111 CLOSE(ifatompaw)
    DEALLOCATE(f)
    DEALLOCATE(shartree,sshartree,wf,twf,rh,rth,dum,dum1)

  END SUBROUTINE WRITE_SOCORROATOMDATA

  Subroutine Report_pseudo_energies(PAW,ien)
       Type(PseudoInfo), INTENT(IN) :: PAW
       Integer , INTENT(IN) :: ien

       write(ien,*)' Summary of PAW energies'
       write(ien,*)'       Total valence energy     ', PAW%Etotal
       write(ien,*)'         Smooth energy          ', PAW%tvale
       write(ien,*)'         One center             ', PAW%Ea
       write(ien,*)'         Smooth kinetic         ', PAW%tkin
       write(ien,*)'         Vloc energy            ', PAW%tion
       write(ien,*)'         Smooth exch-corr       ', PAW%txc
       write(ien,*)'         One-center xc          ', PAW%Eaxc

  End Subroutine Report_pseudo_energies

END MODULE  atompaw_report
