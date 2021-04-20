!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!     LDAGGA_SCF, LDAGGASub, Get_EXC, Get_FCEXC, Report_LDAGGA_functions
!         fixdensity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE ldagga_mod

  USE io_tools
  USE atomdata
  USE anderson_driver
  USE general_mod
  USE excor
  USE globalmath
  USE gridmod
  USE report_mod

  IMPLICIT NONE


  TYPE(GridInfo), PRIVATE, POINTER :: Gridwk
  TYPE(OrbitInfo), PRIVATE, POINTER :: Orbitwk
  TYPE(PotentialInfo), PRIVATE, POINTER :: Potwk
  TYPE(FCInfo), PRIVATE, POINTER :: FCwk
  TYPE(SCFInfo), PRIVATE, POINTER :: SCFwk
  TYPE(Anderson_context), PRIVATE :: AC

CONTAINS

  SUBROUTINE LDAGGA_SCF(scftype,lotsofoutput,Gridin,Orbitin,Potin,FCin,SCFin)
    CHARACTER(2), INTENT(IN) :: scftype
    LOGICAL, INTENT(IN) :: lotsofoutput
    TYPE(GridInfo),TARGET :: Gridin
    TYPE(OrbitInfo),TARGET :: Orbitin
    TYPE(PotentialInfo),TARGET :: Potin
    TYPE(FCInfo),TARGET :: FCin
    TYPE(SCFInfo),TARGET :: SCFin

    INTEGER :: n,i,j
    REAL(8) :: en1,etxc,eex
    REAL(8), ALLOCATABLE :: arg(:)
    LOGICAL :: success
    INTEGER :: counter=0
    CHARACTER(132) :: exctypesave

    Gridwk=>Gridin
    Orbitwk=>Orbitin
    Potwk=>Potin
    FCwk=>FCin
    SCFwk=>SCFin

   
    n=Gridwk%n 
    ALLOCATE(arg(n))

    if(needvtau.and.counter==0) then     ! first converge LDA  works better
   !  if(needvtau.and.counter==0) then     ! first converge GGA-PBE
      needvtau=.false.      
      exctypesave=exctype
      exctype="LDA-PW"
      !exctype="GGA-PBE"
      call initexch
      CALL exch(Gridwk,Orbitwk%den,Potwk%rvx,etxc,eex,&
&       tau=Orbitwk%tau,vtau=Potwk%vtau)
      arg=Potwk%rvh+Potwk%rvx   ! iterating only on electronic part of pot
      Potwk%rv=Potwk%rvh+Potwk%rvx-Potwk%rvx(1)
      CALL zeropot(Gridwk,Potwk%rv,Potwk%v0,Potwk%v0p)
      Potwk%rv=Potwk%rv+Potwk%rvn+Potwk%rvx(1)
      CALL InitAnderson_dr(AC,6,5,n,0.5d0,1.d3,2000,1.d-11,1.d-16,.true.)
      CALL DoAndersonMix(AC,arg,en1,LDAGGAsub,success)

      WRITE(STD_OUT,*) 'Anderson Mix with LDA',AC%res ,' iter = ',AC%CurIter
      !WRITE(STD_OUT,*) 'Anderson Mix with GGA',AC%res ,' iter = ',AC%CurIter
      OPEN (unit=1001,file='potinitLDA',form='formatted')
      !OPEN (unit=1001,file='potinitGGA',form='formatted')
      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau  '      
      DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i)
      ENDDO
      CLOSE(1001)
      CALL FreeAnderson(AC)
      write(STD_OUT,*) 'Completed initial iteration ' ; call flush_unit(std_out)
      counter=1
      needvtau=.true.
      exctype=exctypesave
      call initexch
    Endif  

    CALL exch(Gridwk,Orbitwk%den,Potwk%rvx,etxc,eex,&
&       tau=Orbitwk%tau,vtau=Potwk%vtau)
     !!! testing
     if (needvtau) then
      OPEN (unit=1001,file='potinitRSCAN',form='formatted')
      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau   vtau '      
      DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i),Potwk%rvh(i),&
&         Potwk%rvx(i),Orbitwk%den(i),Orbitwk%tau(i),Potwk%vtau(i)
      ENDDO
      CLOSE(1001)
      endif 
    

    Potwk%rv=Potwk%rvh+Potwk%rvx-Potwk%rvx(1)
    CALL zeropot(Gridwk,Potwk%rv,Potwk%v0,Potwk%v0p)
    Potwk%rv=Potwk%rv+Potwk%rvn+Potwk%rvx(1)

      SCFwk%iter=0
      SCFwk%delta=0

    !arg=Potwk%rv   !no longer used

      arg=Potwk%rvh+Potwk%rvx   ! iterating only on electronic part of pot
      CALL InitAnderson_dr(AC,6,5,n,0.5d0,1.d3,2000,1.d-11,1.d-16,.true.)
      CALL DoAndersonMix(AC,arg,en1,LDAGGAsub,success)
      SCFwk%iter=SCFwk%iter+AC%CurIter
      SCFwk%delta=AC%res

    CALL Report_LDAGGA_functions(scftype)

    CALL FreeAnderson(AC)
    WRITE(STD_OUT,*) 'Finished Anderson Mix', en1 ,' success = ', success
    DEALLOCATE(arg)
    counter=counter+1
  END SUBROUTINE LDAGGA_SCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  LDAGGASub	-- w is electronic part of potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE LDAGGASub(w,energy,residue,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: residue(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io,nw
    REAL(8),ALLOCATABLE :: dum(:)
    REAL(8) :: x
    INTEGER:: fcount=0
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot

    n=Gridwk%n
    nw=SIZE(w)
    ALLOCATE(dum(nw),stat=i)
    IF (i/=0) THEN
        WRITE(6,*) 'Error in LDAGGAsub allocation' ,nw
        STOP
    ENDIF


    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,w,success)
    !tmpPot%rv=w         ! no longer used
      If (.not.success) then   !  attempt to stablize solution
              !   probably needs rethinking....
          write(std_out,*) 'Current eigs', (Orbitwk%eig(io),io=1,Orbitwk%norbit)
          j=n
          x=Orbitwk%eig(1)
          if (Orbitwk%norbit>1) then
             do io = 2, Orbitwk%norbit
                 if (Orbitwk%eig(io)<0.d0.and.x<Orbitwk%eig(io)) &
&                      x=Orbitwk%eig(io)
             enddo
          endif
          write(std_out,*) x,1.d0/sqrt(abs(x))
          x=1.d0/sqrt(abs(x))
          j=FindGridIndex(Gridwk,x)
          write(std_out,*) 'index', x,j,Gridwk%r(j)
          if (j<10)j=10
          if (j>n-10) j=n-10
          w(j+1)=(-1.d0+w(j+1)/2)
          do i=j+2,n
             w(i)=-2.d0
          enddo
           write(std_out,*) 'Reset tmpPot ', j
           write(std_out,*) '   Last points '
              write(std_out,'(1p,20e15.7)') Gridwk%r(n), dum(n),w(n)

           CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,w,success)
           write(std_out,*) 'after updatwfn from reset ',success;
           call flush_unit(std_out)
       Endif

    !if FC core calc , restore core info backinto tmpOrbit

    IF(frozencorecalculation) THEN
       DO io = 1 , Orbitwk%norbit
          IF(Orbitwk%iscore(io)) THEN
             tmpOrbit%eig(io)=Orbitwk%eig(io)
             tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
             if(diracrelativistic) tmpOrbit%lwfn(:,io)=Orbitwk%lwfn(:,io)
             tmpOrbit%otau(:,io)=Orbitwk%otau(:,io)
          ENDIF
       ENDDO
    ENDIF

    IF (.NOT.success) THEN
       WRITE(STD_OUT,*) 'Bad luck in Sub'
    ENDIF

    !write(std_out,*) 'in LDAGGAsub before Get'; call flush_unit(std_out)
    CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)

    CALL Fixdensity(Gridwk,tmpOrbit%den)
    
    !write(std_out,*) 'in LDAGGAsub before EXC'; call flush_unit(std_out)
    CALL Get_EXC(Gridwk,tmpPot,tmpOrbit,SCFwk)
    !write(std_out,*) 'after Get_EXC'; call flush_unit(std_out)
    dum(1:n)=tmpPot%rvh(1:n)+tmpPot%rvx(1:n)-w(1:n)
    !write(std_out,*) 'after Get_EXC'; call flush_unit(std_out)
    residue=dum
    err=Dot_Product(residue,residue)
    write(STD_OUT,*) 'in LDAGGASub   err ', err;call flush_unit(STD_OUT)
    
   if(frozencorecalculation) then
     Call Get_FCKinCoul(Gridwk,tmpPot,tmpOrbit,FCwk,SCFwk)
     CALL Get_FCEXC(SCFwk)
     energy=SCFwk%evale
     CALL Total_FCEnergy_Report(SCFwk,STD_OUT)
   else
     energy=SCFwk%etot
     CALL Total_Energy_Report(SCFwk,STD_OUT)
   endif


    IF (update) THEN
       Potwk%rv=w+tmpPot%rvn
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx
       if (needvtau) Potwk%vtau=tmpPot%vtau
    !!!!   if (needvtau) Potwk%vtau=0.5d0*tmpPot%vtau+0.5d0*Potwk%vtau
    !!!!   if (needvtau) Potwk%vtau=0.1d0*tmpPot%vtau+0.9d0*Potwk%vtau
    !!!!     !   needed to stabilize calculation -- may need adjustments
       Orbitwk%wfn=tmpOrbit%wfn
       If(diracrelativistic)Orbitwk%lwfn=tmpOrbit%lwfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den
       Orbitwk%otau=tmpOrbit%otau
       Orbitwk%tau=tmpOrbit%tau
       Orbitwk%deltatau=tmpOrbit%deltatau
       Call One_electron_energy_Report(Orbitwk,std_out)
       !!! testing
       if(needvtau.and.fcount==0) then
      OPEN (unit=1001,file='potwithkedR2SCAN',form='formatted')
      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau   vtau'      
      DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i),Potwk%vtau(i)
      ENDDO
      CLOSE(1001)
      fcount=fcount+1
      endif
       if(needvtau.and.fcount==1) then
      OPEN (unit=1001,file='potwithkedR2SCAN1',form='formatted')
      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau   vtau'      
      DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i),Potwk%vtau(i)
      ENDDO
      CLOSE(1001)
      fcount=fcount+1
      endif
       if(needvtau.and.fcount==2) then
      OPEN (unit=1001,file='potwithkedR2SCAN2',form='formatted')
      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau   vtau'      
      DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i),Potwk%vtau(i)
      ENDDO
      CLOSE(1001)
      fcount=fcount+1
      endif
    !!! end testing
    ENDIF


    !write(std_out,*) 'in LDAGGAsub before end'; call flush_unit(std_out)
    CALL DestroyOrbit(tmpOrbit)
    CALL DestroyPot(tmpPot)
    DEALLOCATE (dum)

  END SUBROUTINE  LDAGGASub

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   For reasons that are not completely clear, the electron
  !    density near the origin occasionally becomes erratic
  !    this routine checks and fixes this behavior for r<0.01 bohr
  !  SUBROUTINE fixdensity(Grid,den)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fixdensity(Grid,den)
    TYPE (GridInfo), INTENT(IN) :: Grid
    REAL(8), INTENT(INOUT) :: den(:)

    REAL(8), ALLOCATABLE :: tmpd(:),dum(:),grad(:)
    INTEGER :: n,i,j,k
    REAL(8) :: fpi,con,a,b
    REAL(8), parameter :: range=0.01d0, tol=0.01d0
    LOGICAL :: problem
    INTEGER :: counter=0

    fpi=4*pi
    n=Grid%n
    allocate(tmpd(n),dum(n),grad(n))

       tmpd=0
       DO i=2,n
          tmpd(i)=den(i)/(fpi*(Grid%r(i)**2))
       ENDDO
       CALL extrapolate(Grid,tmpd)
       do i=1,n
         dum(i)=ddlog(abs(tmpd(i)))
       enddo  
       call derivative(Grid,dum,grad,1,n)

       problem=.false.; j=0
       do i=2,n
          if (Grid%r(i)>range) exit
          if (abs(grad(i)-grad(i-1))>tol) then
            problem=.true.
            j=i      
            exit
          endif
       enddo  

       if (problem) then
         write(std_out,*) 'density problem   #', counter      
         k=0
         do i=j,n
          if (Grid%r(i)>range) exit
          if (abs(grad(i)-grad(i-1))<tol) then
            k=i      
            exit
          endif
         enddo  
         write(std_out,*) 'Resetting density for < ', k,Grid%r(k)
         con=(Grid%r(k+2)*dum(k+1)-Grid%r(k+1)*dum(k+2))/(Grid%r(k+2)-Grid%r(k+1))
         a=(dum(k+1)-dum(k+2))/(Grid%r(k+2)-Grid%r(k+1))
         write(std_out,*) 'con, a ', con,a
         con=fpi*exp(con)
         do i=1,k+2
           den(i)=con*exp(-a*Grid%r(i))*(Grid%r(i)**2)
         enddo
        ! do i=1,n
        !   write(5000+counter,'(1p,10e15.7)') Grid%r(i),den(i),tmpd(i),dum(i),grad(i)      
        ! enddo  
       endif  
         counter=counter+1
    deallocate(tmpd,dum,grad)
  END SUBROUTINE fixdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get_EXC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_EXC(Grid,Pot,Orbit,SCF)
    !  program to calculate exc energy and potential
    !     assumes Orbit%den already known
    !     also assume kinetic and coulomb energies
    !       calculated and stored in SCF
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(INOUT) :: Pot
    TYPE(Orbitinfo), INTENT(INOUT) :: Orbit
    TYPE(SCFInfo), INTENT(OUT) :: SCF

    REAL(8) :: eex,etot,etxc
    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: k,n

    n=Grid%n
    !write(std_out,*) 'In Get_EXC', n; call flush_unit(std_out)
    CALL exch(Grid,Orbit%den,Pot%rvx,etxc,eex,&
&       tau=Orbit%tau,vtau=Pot%vtau)
    !write(std_out,*) 'After exch', etxc,eex; call flush_unit(std_out)

    SCF%eexc=eex
    etot = SCF%ekin+SCF%estatic+SCF%eexc
    SCF%etot=etot
    WRITE(STD_OUT,*) '    Total                    :  ',etot

    !write(std_out,*) 'before allocate'; call flush_unit(std_out)
    ALLOCATE(dum(n))
    dum=0
    !write(std_out,*) 'before dum'; call flush_unit(std_out)
    dum(2:n)=Pot%rvx(2:n)*Orbit%den(2:n)/Grid%r(2:n)
    !write(std_out,*) 'after dum'; call flush_unit(std_out)
    !write(std_out,*) SCF%eone;call flush_unit(std_out)
    !write(std_out,*) SCF%ecoul;call flush_unit(std_out)
    !write(std_out,*) eex;call flush_unit(std_out)
    !write(std_out,*) integrator(Grid,dum);call flush_unit(std_out)
    WRITE(STD_OUT,*) '    Total   (DC form)        :  ',&
&        SCF%eone-SCF%ecoul+eex-integrator(Grid,dum)
    call flush_unit(std_out)
    DEALLOCATE(dum)
  END SUBROUTINE Get_EXC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Get_FCEXC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Get_FCEXC(SCF)
    !  program to get exc from SCF%eexc and to
    !    valence kinetic and couloub terms
    !    to accumulate SCF%evale
    !     assumes Orbit%den already known
    TYPE(SCFInfo), INTENT(INOUT) :: SCF

    SCF%valeexc=SCF%eexc   !due to nonlinear dependence
                           ! of core and valence contributions
                           ! in lda and gga
    SCF%evale=SCF%valekin+SCF%valecoul+SCF%valeexc

  END SUBROUTINE Get_FCEXC



  SUBROUTINE Report_LDAGGA_functions(sub)
    CHARACTER(2) :: sub
    INTEGER :: i,j,k,n,nmap

    INTEGER, SAVE :: counter=0
    CHARACTER(4) :: stuff

    CALL mkname(counter,stuff)

     write(std_out,*) 'in Report ', stuff; call flush_unit(std_out)

    OPEN (unit=1001,file='pot'//sub//TRIM(stuff),form='formatted')
    WRITE(1001,*)'#    r         rv               rvh           rvx       den                  tau          vtau          deltatau '
   
    n=Gridwk%n
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i),&
&            Potwk%vtau(i),Orbitwk%deltatau(i)
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
    if (.not.diracrelativistic) then
    WRITE(1001,*) '#         r          wfn in order of s, p, d ... '
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&           (Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
    ENDDO
    endif
    if (diracrelativistic) then
    WRITE(1001,*) '#         r          wfn, lwfn in order of s, p-1/2, p+1/2 ... '
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&           (Orbitwk%wfn(i,j),Orbitwk%lwfn(i,j),j=1,Orbitwk%norbit)
    ENDDO
    endif

    CLOSE(1001)


    counter=counter+1

  END SUBROUTINE Report_LDAGGA_functions

END MODULE ldagga_mod
