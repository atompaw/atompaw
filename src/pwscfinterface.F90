!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains the subroutines Atompaw2PWscf and libxc2UPF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

Module PWscfInterface

  USE Tools
  USE GlobalMath
  USE atomdata
  USE excor
  USE gridmod
  USE interpolation_mod
  USE pseudo
  USE libxc_mod

  Implicit none

  CONTAINS

  SUBROUTINE Atompaw2PWscf(Grid,Pot,FC,PAW,Orbit,ifinput)
     Type(GridInfo), INTENT(IN) :: Grid
     Type(PotentialInfo), INTENT(IN) :: Pot
     Type(FCInfo), INTENT(IN) :: FC
     Type(PseudoInfo), INTENT(IN) :: PAW
     Type (OrbitInfo), INTENT(IN) :: Orbit
     INTEGER, INTENT(IN) :: ifinput

     INTEGER :: i,j,k,l,n,m,io,jo,llmin,llmax,number_of_wfc,nn,ok, ncore
     CHARACTER(132) :: inputfileline,stuff
     CHARACTER(4) :: s1,s2,s3,s4
     CHARACTER(1), parameter :: label(4) = 'SDPF'
     REAL(8) :: x,y,z,q
     REAL(8) :: upfdx=0.005d0,upfxmin=-9.d0,upfzmesh=1.d0 ! defaults
     REAL(8) :: upfrange=15.d0       !default
     INTEGER :: upfmesh,upfirc
     REAL(8), POINTER :: vps(:)
     REAL(8), ALLOCATABLE :: dum(:),arg(:),rvion(:),rtvion(:),dij0(:,:)
     REAL(8), ALLOCATABLE :: upfr(:),upff(:)
     CHARACTER(256) :: inputline
     CHARACTER(50) :: UPFlabel,tmplabel
     LOGICAL :: testing
     ALLOCATE(dum(Grid%n),arg(Grid%n),rvion(Grid%n),rtvion(Grid%n),&
&            dij0(PAW%nbase,PAW%nbase))
   
    testing=.false.
    if (TRIM(exctype)=='LDA-PW')  then
           UPFlabel='SLA PW NOGX NOGC'    
           tmplabel=TRIM(exctype)
           testing=.true.
    elseif (TRIM(exctype)=='GGA-PBE')  then
           UPFlabel='SLA PW PBX PBC'
           tmplabel=TRIM(exctype)
           testing=.true.
    elseif (TRIM(exctype)=='GGA-PBESOL')  then
           UPFlabel='SLA PW PSX PSC'    
           tmplabel=TRIM(exctype)
           testing=.true.
    elseif (.not.testing.and.have_libxc) then
           call libxc2UPF(exctype,UPFlabel,tmplabel,testing)
           if (.not.testing) then
              write(6,'(/,2x,a)') "Error in Atompaw2PWscf:"
              write(6,'(2x,a)') " libxc2UPF incorrectly executed "
              write(6,'(2x,a,a)') " exctype = ", exctype
              write(6,'(2x,a,a)') " UPFlabel = ", UPFlabel
              write(6,'(2x,a,a)') " tmplabel = ", tmplabel
              write(6,'(2x,a)') "   Please contact natalie@wfu.edu"
              stop 
           endif   
     else      
           write(6,'(/,2x,a)') "Error in Atompaw2PWscf:"
           write(6,'(2x,a,a)') " exctype = ", exctype
           write(6,'(2x,a)') "   Limited functionals supported for UPF output "
           write(6,'(2x,a)') "   No UPDF file written !"
           write(6,'(2x,a)') "   Please contact natalie@wfu.edu"
           return
     end if

  ! set up grid
     read(5,'(a)',advance='no',iostat=ok) inputline
     if (ok<=0) then
        write(6,'(a)') inputline
        call Uppercase(inputline)
        i=0;i=INDEX(inputline,'UPFDX')+INDEX(inputline,'UPFXMIN')+&
&              INDEX(inputline,'UPFZMESH')+INDEX(inputline,'UPFRANGE')
        if (i>0) then
           k=0;k=INDEX(inputline,'UPFDX')
             if (k>0) then
                read(inputline(k+5:256),*) upfdx
                write(6,*) 'upfdx reset ', upfdx
             endif
           k=0;k=INDEX(inputline,'UPFXMIN')
             if (k>0) then
                read(inputline(k+7:256),*) upfxmin
                write(6,*) 'upfxmin reset ', upfxmin
             endif
           k=0;k=INDEX(inputline,'UPFZMESH')
             if (k>0) then
                read(inputline(k+8:256),*) upfzmesh
                write(6,*) 'upfzmesh reset ', upfzmesh
             endif
           k=0;k=INDEX(inputline,'UPFRANGE')
             if (k>0) then
                read(inputline(k+8:256),*) upfrange
                write(6,*) 'upfrange reset ', upfrange
             endif
         endif
    endif

         !upfmesh=1+(LOG(upfzmesh*Grid%r(Grid%n))-upfxmin)/upfdx
         upfmesh=1+(LOG(upfzmesh*upfrange)-upfxmin)/upfdx

         write(6,*) 'UPF mesh size = ', upfmesh
         write(6,*) 'UPF xmin = ', upfxmin
         write(6,*) 'UPF zmesh = ', upfzmesh
         write(6,*) 'UPF zmesh = ', upfrange

         ALLOCATE(upfr(upfmesh),upff(upfmesh))

         upfr=0;upff=0
         do i=1,upfmesh
            upfr(i)=EXP(upfxmin+upfdx*(i-1))/upfzmesh
         enddo

         upfirc=upfmesh
         do i=1,upfmesh
            if (upfr(i)>Grid%r(PAW%irc)) then
               upfirc=i+1
               exit
            endif
         enddo

         write(6,*) 'UPF irc = ',   upfirc

     vps=>PAW%abinitvloc       ! could have flag for using PAW%abinitnohat
     rtvion=vps*Grid%r
     n=Grid%n
   ! reconstruct AE ionic potential
      dum=-2*Pot%nz
      call poisson(Grid,x,FC%coreden,arg,y,z)
      dum=dum+arg
      rvion=dum
   ! reconstruct DIJ0 coefficients (Corresponds to Eq. 97 in atompawEqns.pdf)
      dij0=0
      Do io=1,PAW%nbase
         Do jo=1,PAW%nbase
            If (PAW%l(io)==PAW%l(jo)) then
               dij0(io,jo)=PAW%Kij(io,jo)
               arg=rvion*PAW%ophi(:,io)*PAW%ophi(:,jo) - &
&                       rtvion*(PAW%otphi(:,io)*PAW%otphi(:,jo) + &
&                             PAW%oij(io,jo)*PAW%hatden)
               arg(1)=0.d0; arg(2:Grid%n)=arg(2:Grid%n)/Grid%r(2:Grid%n)
               dij0(io,jo)=dij0(io,jo)+integrator(Grid,arg,1,PAW%irc)
            endif
         enddo
      enddo

      OPEN(1001,file=TRIM(Pot%sym)//'.'//TRIM(tmplabel)//'-paw.UPF',&
 &          form='formatted')
      write(1001,'("<UPF version=""2.0.1"">")')
      write(1001,'("   <PP_INFO>")')
      write(1001,'("         UPF file from ATOMPAW code with following input")')

      open(ifinput,file='dummy',form='formatted')
      do
        read(ifinput,'(a)',iostat=i) inputfileline
           if (i/=0) exit
           write(1001,'("           ",a)') TRIM(inputfileline)
      enddo
      close(ifinput)

      WRITE(1001,'("  </PP_INFO>")')
      WRITE(1001,'("  <!--                               -->")')
      WRITE(1001,'("  <!-- END OF HUMAN READABLE SECTION -->")')
      WRITE(1001,'("  <PP_HEADER generated=""Generated using ATOMPAW code""")')
      WRITE(1001,'("             author=""   """)')
      Call PrintDateStr(inputfileline)
      WRITE(1001,'("             date=""",(a),"""")') TRIM(inputfileline)
      WRITE(1001,'("             element=""",a2,"""")') Pot%sym
      WRITE(1001,'("             pseudo_type=""PAW""")')
      If(scalarrelativistic) then
      WRITE(1001,'("             relativistic=""scalar""")')
      else
      WRITE(1001,'("             relativistic=""no""")')
      endif
      WRITE(1001,'("             is_ultrasoft=""T""")')
      WRITE(1001,'("             is_paw=""T""")')
      WRITE(1001,'("             is_coulomb=""F""")')
      WRITE(1001,'("             has_so=""F""")')
      WRITE(1001,'("             has_wfc=""T""")')
      WRITE(1001,'("             has_gipaw=""T""")')
      WRITE(1001,'("             paw_as_gipaw=""T""")')
      WRITE(1001,'("             core_correction=""T""")')
      WRITE(1001,'("             functional=""",(a),"""")') TRIM(UPFlabel)
      WRITE(1001,'("             z_valence=""",f6.3,"""")')FC%zvale
      WRITE(1001,'("             l_max=""",i1,"""")') PAW%lmax
      WRITE(1001,'("             l_max_rho=""",i1,"""")') 2*PAW%lmax
      WRITE(1001,'("             mesh_size=""",i6,"""")') n-1
      WRITE(1001,'("             number_of_wfc=""",i6,"""")') PAW%nbase
      WRITE(1001,'("             number_of_proj=""",i2,"""")')PAW%nbase
      WRITE(1001,'("                                    /> ")')

      WRITE(1001,'("<PP_MESH dx=""",1p,1e20.13,""" mesh=""",i6,""" xmin=""",&
&      1p,1e20.13,""" zmesh=""",1p,1e20.13,""">")')upfdx,upfmesh,upfxmin,upfzmesh
      WRITE(1001,'(" <PP_R type=""real"" size=""",i6,""" columns=""3"">")')&
&                  upfmesh
      WRITE(1001,'(1p,3e25.17)') (upfr(i),i=1,upfmesh)
      WRITE(1001,'("                               </PP_R> ")')
      WRITE(1001,'(" <PP_RAB type=""real"" size=""",i6,""" columns=""3"">")')&
&                upfmesh
      WRITE(1001,'(1p,3e25.17)') (upfdx*upfr(i),i=1,upfmesh)
      WRITE(1001,'("                             </PP_RAB> ")')
      WRITE(1001,'("                                  </PP_MESH>")')

      WRITE(1001,'(" <PP_NLCC type=""real"" size=""",i6,""" columns=""3"">")')&
&                upfmesh
      dum=0
      dum(2:n)=PAW%tcore(2:n)/(4*pi*(Grid%r(2:n))**2)
      call extrapolate(Grid,dum)
      upff=0;call interpfunc(n,Grid%r,dum,upfmesh,upfr,upff)
      call filter(upfmesh,upff,machine_zero)
      WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
      WRITE(1001,'("                             </PP_NLCC> ")')


      WRITE(1001,'(" <PP_LOCAL type=""real"" size=""",i6,""" columns=""3"">")')&
&                    upfmesh
      upff=0;call interpfunc(n,Grid%r,vps,upfmesh,upfr,upff)
      call filter(upfmesh,upff,machine_zero)
      WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
      WRITE(1001,'("                             </PP_LOCAL> ")')


      WRITE(1001,'("<PP_NONLOCAL> ")')
      do io=1,PAW%nbase
        call mkname(io,s1)  ; inputfileline=' <PP_BETA.'//TRIM(s1)
        WRITE(1001,&
&       '((a)," type=""real"" size=""",i6,""" angular_momentum=""",i1,&
&       """ columns=""3"" cutoff_radius_index=""",i6,""" norm_conserving_radius=""",F12.6,""">")')&
&           TRIM(inputfileline),upfmesh,PAW%l(io),upfirc, Grid%r(PAW%irc)
            upff=0;call interpfunc(n,Grid%r,PAW%otp(:,io),upfmesh,upfr,upff)
            call filter(upfmesh,upff,machine_zero)
            upff(upfirc+1:upfmesh)=0.d0
            WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
         WRITE(1001,'(30x,(a))') '</PP_BETA.'//TRIM(s1)//'>'
      enddo

      WRITE(1001,'("   <PP_DIJ type=""real"" size=""",i5,""" columns=""3"">")')&
&                 PAW%nbase*PAW%nbase
         WRITE(1001,'(1p,3e25.17)') ((dij0(i,j),i=1,PAW%nbase),j=1,PAW%nbase)
      WRITE(1001,'("                                       </PP_DIJ>")')

      stuff='SINC**2'
      if (gaussianshapefunction) stuff='GAUSSIAN'
      if (besselshapefunction) stuff='BESSEL'
      WRITE(1001,'("  <PP_AUGMENTATION q_with_l=""T"" nqf=""0"" cutoff_r=""",&
&      1p,1e15.7,""" cutoff_r_index=""",i6,""" iraug=""",i6,&
&      """ nqlc=""",i3,""" lmax_aug=""",i3,""" shape=""",(a),&
&      """ augmentation_epsilon=""1.d-12"">")')  PAW%rc,upfirc,upfirc,&
&          2*PAW%lmax+1,2*PAW%lmax,TRIM(stuff)

      WRITE(1001,'("    <PP_Q type=""real"" size=""",i5,""" columns=""3"">")')&
&                 PAW%nbase*PAW%nbase
         WRITE(1001,'(1p,3e25.17)') ((PAW%oij(i,j),i=1,PAW%nbase),j=1,PAW%nbase)
      WRITE(1001,'("                                       </PP_Q>")')

      WRITE(1001,'("    <PP_MULTIPOLES type=""real"" size=""",&
&       i6,""" columns=""3"">")') PAW%nbase*PAW%nbase*(2*PAW%lmax+1)
      WRITE(1001,'(1p,3e25.17)') (((PAW%mLij(i,j,k),i=1,PAW%nbase),&
&                      j=1,PAW%nbase), k=1,(2*PAW%lmax+1))
      WRITE(1001,'("                                       </PP_MULTIPOLES>")')
      do io=1,PAW%nbase
         do jo=io,PAW%nbase
            llmin=ABS(PAW%l(io)-PAW%l(jo))
            llmax=PAW%l(io)+PAW%l(jo)
            do l=llmin,llmax,2
               call hatL(Grid,PAW,l,dum)
               dum=dum*PAW%mLij(io,jo,l+1)
               arg=dum*(Grid%r**l)
               write(6,*) 'Chk aug', io,jo,l,integrator(Grid,arg)
               call mkname(io,s1)
               call mkname(jo,s2)
               call mkname(l,s3)
               inputfileline=' <PP_QIJL.'//TRIM(s1)//'.'//TRIM(s2)//'.'//TRIM(s3)
               WRITE(1001,'((a)," type=""real"" size=""",i6,&
&                """ columns=""3"" first_index=""",i2,&
&                """ second_index=""",i2,""" composite_index=""",i2,&
&                """ angular_momentum=""",i1,""">")') &
&                TRIM(inputfileline),upfmesh,io,jo,(jo*(jo-1))/2+io,l
               upff=0;call interpfunc(n,Grid%r,dum,upfmesh,upfr,upff)
               call filter(upfmesh,upff,machine_zero)
               upff(upfirc+1:upfmesh)=0.d0
               WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
               inputfileline=' </PP_QIJL.'//TRIM(s1)//'.'//TRIM(s2)//'.'//TRIM(s3)//'>'
              WRITE(1001,'("                      ", (a))') TRIM(inputfileline)
            enddo
          enddo
       enddo
      WRITE(1001,'("                                   </PP_AUGMENTATION> ")')
      WRITE(1001,'("                                   </PP_NONLOCAL> ")')

      WRITE(1001,'("<PP_PSWFC> ")')
      j=0
      do io=1,PAW%nbase
            j=j+1
            call mkname(j,s1)  ; inputfileline=' <PP_CHI.'//TRIM(s1)
            WRITE(1001, '((a)," type=""real"" size=""",i6,""" l=""",i1,&
&             """ occupation=""",f7.4,""" columns=""3"">")')&
&           TRIM(inputfileline),upfmesh, PAW%l(io),PAW%occ(io)
            upff=0;call interpfunc(n,Grid%r,PAW%otphi(:,io),upfmesh,upfr,upff)
            call filter(upfmesh,upff,machine_zero)
            WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
            WRITE(1001,'(30x,(a))') '</PP_CHI.'//TRIM(s1)//'>'
      enddo
      WRITE(1001,'("                             </PP_PSWFC> ")')

      !!! Note number_of_wfc is ambiguous in PWscf code.   For now,
      !!!    set to the same as the number of basis functions
      WRITE(1001,'("<PP_FULL_WFC number_of_wfc=""",i2,""">")')PAW%nbase
      do io=1,PAW%nbase
        call mkname(io,s1)  ; inputfileline=' <PP_AEWFC.'//TRIM(s1)
        WRITE(1001, '((a)," type=""real"" size=""",i6,&
&       """ l=""",i1,""" columns=""3"">")')&
&           TRIM(inputfileline),upfmesh,PAW%l(io)
            upff=0;call interpfunc(n,Grid%r,PAW%ophi(:,io),upfmesh,upfr,upff)
            call filter(upfmesh,upff,machine_zero)
            WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
         WRITE(1001,'(30x,(a))') '</PP_AEWFC.'//TRIM(s1)//'>'
      enddo

      do io=1,PAW%nbase
        call mkname(io,s1)  ; inputfileline=' <PP_PSWFC.'//TRIM(s1)
        WRITE(1001, '((a)," type=""real"" size=""",i6,&
&       """ l=""",i1,""" columns=""3"">")')&
&           TRIM(inputfileline),upfmesh,PAW%l(io)
            upff=0;call interpfunc(n,Grid%r,PAW%otphi(:,io),upfmesh,upfr,upff)
            call filter(upfmesh,upff,machine_zero)
            WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
         WRITE(1001,'(30x,(a))') '</PP_PSWFC.'//TRIM(s1)//'>'
      enddo
      WRITE(1001,'("                             </PP_FULL_WFC> ")')

      arg=PAW%den-PAW%tden
      q=integrator(Grid,arg)
      write(6,*) 'check augmentation charge ', q
      dum=PAW%tden+q*PAW%hatden
      write(6,*) 'check valence charge ', integrator(Grid,dum)
      WRITE(1001,'(" <PP_RHOATOM type=""real"" size=""",i6,&
&            """ columns=""3"">")')upfmesh
            upff=0;call interpfunc(n,Grid%r,dum,upfmesh,upfr,upff)
            call filter(upfmesh,upff,machine_zero)
            WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
      WRITE(1001,'("                             </PP_RHOATOM> ")')

      WRITE(1001,'("<PP_PAW paw_data_format=""2""> ")')
      WRITE(1001,&
&        '("   <PP_OCCUPATIONS type=""real"" columns=""3"" size=""",i2,""">")')&
&         PAW%nbase
      WRITE(1001,'(1p,3e25.17)') (PAW%occ(io),io=1,PAW%nbase)
      WRITE(1001,'("                             </PP_OCCUPATIONS> ")')

      WRITE(1001,&
&       '(" <PP_AE_NLCC type=""real"" size=""",i6,""" columns=""3"">")')upfmesh
      dum=0
      dum(2:n)=FC%coreden(2:n)/(4*pi*(Grid%r(2:n))**2)
      call extrapolate(Grid,dum)
        upff=0;call interpfunc(n,Grid%r,dum,upfmesh,upfr,upff)
        call filter(upfmesh,upff,machine_zero)
        WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
      WRITE(1001,'("                             </PP_AE_NLCC> ")')

      upff=0;call interpfunc(n,Grid%r,rvion,upfmesh,upfr,upff)
      call filter(upfmesh,upff,machine_zero)
      upff=upff/upfr
      WRITE(1001,&
&      '(" <PP_AE_VLOC type=""real"" size=""",i6,""" columns=""3"">")')upfmesh
      WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
      WRITE(1001,'("                             </PP_AE_VLOC> ")')
      WRITE(1001,'("                                      </PP_PAW> ")')

      WRITE(1001,'(A)') ' <PP_GIPAW gipaw_data_format="2">'
      ncore = 0
      do io = 1, Orbit%norbit
        if (Orbit%iscore(io) .eqv. .true.) ncore = ncore + 1
      enddo
      WRITE(1001,'(A,I2,A)') '  <PP_GIPAW_CORE_ORBITALS number_of_core_orbitals="', ncore, '">'
      ncore = 0
      do io = 1, Orbit%norbit
        if (Orbit%iscore(io) .eqv. .false.) cycle
        ncore = ncore + 1
        if (ncore.lt.10) then
           WRITE(1001,'("   <PP_GIPAW_CORE_ORBITAL.",I1," type=""real"" size=""",i6,"""")',advance="no") ncore, upfmesh
           WRITE(1001,'(" columns=""3"" index=""",I1,""" label=""",I1,A,""" n=""",I1,""" l=""",I1,""">")') &
&          ncore, Orbit%np(io), label(Orbit%l(io)+1), Orbit%np(io), Orbit%l(io)
        else
           WRITE(1001,'("   <PP_GIPAW_CORE_ORBITAL.",I2," type=""real"" size=""",i6,"""")',advance="no") ncore, upfmesh
           WRITE(1001,'(" columns=""3"" index=""",I2,""" label=""",I1,A,""" n=""",I1,""" l=""",I1,""">")') &
&          ncore, Orbit%np(io), label(Orbit%l(io)+1), Orbit%np(io), Orbit%l(io)
        endif   
        upff=0;call interpfunc(n,Grid%r,Orbit%wfn(:,io),upfmesh,upfr,upff)
        call filter(upfmesh,upff,machine_zero)
        WRITE(1001,'(1p,3e25.17)') (upff(i),i=1,upfmesh)
        if (ncore.lt.10) then
           WRITE(1001,'("   </PP_GIPAW_CORE_ORBITAL.",I1,">")') ncore
        else   
           WRITE(1001,'(" </PP_GIPAW_CORE_ORBITAL.",I2,">")') ncore
        endif   
      enddo

      WRITE(1001,'(A)') '  </PP_GIPAW_CORE_ORBITALS>'
      WRITE(1001,'(A)') ' </PP_GIPAW>'

      WRITE(1001,'("</UPF> ")')

      CLOSE(1001)

    DEALLOCATE(dum,arg,rvion,rtvion,dij0,upfr,upff)
  END SUBROUTINE Atompaw2PWscf

   SUBROUTINE libxc2UPF(libxcform,qeform,shortqe,OK)
   IMPLICIT NONE
   CHARACTER(50), INTENT(IN) :: libxcform
   CHARACTER(50), INTENT(OUT) :: qeform,shortqe
   Logical, INTENT(OUT) :: OK

   OK=.false.

   if (TRIM(libxcform)=='XC_LDA_X+XC_LDA_C_PW') then
       qeform='SLA PW NOGX NOGC'
       shortqe='LDA-PW'
       OK=.true.
       return
   elseif (TRIM(libxcform)=='XC_LDA_X+XC_LDA_C_PZ') then
       qeform='SLA PZ NOGX NOGC'
       shortqe='LDA-PZ'
       OK=.true.
       return
   elseif (TRIM(libxcform)=='XC_GGA_X_PBE_SOL+XC_GGA_C_PBE_SOL') then
       qeform='SLA PW PSX PSC' 
       shortqe='GGA-PBESOL'
       OK=.true.
       return
   elseif (TRIM(libxcform)=='XC_GGA_X_PBE+XC_GGA_C_PBE') then
       qeform='SLA PW PBX PBC' 
       shortqe='GGA-PBE'
       OK=.true.
       return
   elseif (TRIM(libxcform)=='XC_GGA_X_WC+XC_GGA_C_PBE') then
       qeform='SLA PW WCX PBC' 
       shortqe='GGA-WC'
       OK=.true.
       return
   endif 

  END   SUBROUTINE libxc2UPF


END Module PWscfInterface
