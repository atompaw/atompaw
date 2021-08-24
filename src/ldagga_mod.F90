!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!     LDAGGA_SCF, LDAGGASub, Get_EXC, Get_FCEXC, Report_LDAGGA_functions,
!          DENITERSub
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
!!!      OPEN (unit=1001,file='potinitLDA',form='formatted')
!!!      !OPEN (unit=1001,file='potinitGGA',form='formatted')
!!!      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau  '      
!!!      DO i = 1,n
!!!       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
!!!&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i)
!!!      ENDDO
!!!      CLOSE(1001)
      CALL FreeAnderson(AC)
      write(STD_OUT,*) 'Completed initial iteration ' ; call flush_unit(std_out)
      needvtau=.true.
      exctype=exctypesave
      call initexch
      counter=1

!!!    CALL exch(Gridwk,Orbitwk%den,Potwk%rvx,etxc,eex,&
!!!&       tau=Orbitwk%tau,vtau=Potwk%vtau)
!!!     !!! testing
!!!     if (needvtau) then
!!!      OPEN (unit=1001,file='potinitR2SCAN0',form='formatted')
!!!      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau   vtau '      
!!!      DO i = 1,n
!!!       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i),Potwk%rvh(i),&
!!!&         Potwk%rvx(i),Orbitwk%den(i),Orbitwk%tau(i),Potwk%vtau(i)
!!!      ENDDO
!!!      CLOSE(1001)
!!!      endif 
!!!      counter=1
    Endif  
    
    CALL exch(Gridwk,Orbitwk%den,Potwk%rvx,etxc,eex,&
&       tau=Orbitwk%tau,vtau=Potwk%vtau)

    Potwk%rv=Potwk%rvh+Potwk%rvx-Potwk%rvx(1)
    CALL zeropot(Gridwk,Potwk%rv,Potwk%v0,Potwk%v0p)
    Potwk%rv=Potwk%rv+Potwk%rvn+Potwk%rvx(1)

      SCFwk%iter=0
      SCFwk%delta=0

      CALL InitAnderson_dr(AC,6,5,n,0.5d0,1.d3,2000,1.d-11,1.d-16,.true.)
      If (needvtau) then
         arg=Orbitwk%den   ! iterating on density
         CALL DoAndersonMix(AC,arg,en1,DENITERsub,success)
         !!!!   evaluate vxc and vtau on universal grid
         !!!!    because it may be bumpy at intermediate range from spline
         !!!!     evaluation
         CALL exch(Gridwk,Orbitwk%den,Potwk%rvx,etxc,eex,&
&           tau=Orbitwk%tau,vtau=Potwk%vtau)
         Potwk%rv=Potwk%rvn+Potwk%rvh+Potwk%rvx   
      else   
         arg=Potwk%rvh+Potwk%rvx  ! iterating only on electronic part of pot
         CALL DoAndersonMix(AC,arg,en1,LDAGGAsub,success)
      endif   
      SCFwk%iter=SCFwk%iter+AC%CurIter
      SCFwk%delta=AC%res

      WRITE(STD_OUT,*) 'Anderson Mix ',success,AC%res ,' iter = ',AC%CurIter
      write(std_out,*) 'line 136 ldagga ', scftype; call flush_unit(std_out)
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
    INTEGER:: fcount=0,dcount=0
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
    write(std_out,*) 'completed updatewfn with success ', success
    !tmpPot%rv=w         ! no longer used
      If (.not.success) then   !  attempt to stablize solution
              !   probably needs rethinking....
          write(std_out,*) 'Current eigs', (Orbitwk%eig(io),io=1,Orbitwk%norbit)
          
!!!!!!!!!!!!!!!!!!!    The following code was taken out --
!!!!   but is apparently needed
!!!!
!!!!          write(std_out,*) 'Eigensolver failed and recovery algorithm is not worthy -- program stopping' 
!!!!          stop
!!!!!!!!!!!!!!!!!!    end of  modification

          j=n
          x=Orbitwk%eig(1)
          if (Orbitwk%norbit>1) then
             do io = 2, Orbitwk%norbit
                 if (Orbitwk%eig(io)<0.d0.and.x<Orbitwk%eig(io)) &
&                      x=Orbitwk%eig(io)
             enddo
          endif
          !write(std_out,*) x,1.d0/sqrt(abs(x))
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

    !CALL Fixdensity(Gridwk,tmpOrbit%den)
    
    !write(std_out,*) 'in LDAGGAsub before EXC'; call flush_unit(std_out)
    !write(std_out,*)  'Check tau ', integrator(Gridwk,tmpOrbit%tau)
    CALL Get_EXC(Gridwk,tmpPot,tmpOrbit,SCFwk)
    !write(std_out,*) 'after Get_EXC'; call flush_unit(std_out)
    dum(1:n)=tmpPot%rvh(1:n)+tmpPot%rvx(1:n)-w(1:n)
    !write(std_out,*) 'after Get_EXC'; call flush_unit(std_out)
    residue=dum
    err=Dot_Product(residue,residue)
    write(STD_OUT,*) 'in LDAGGASub   err ', err;call flush_unit(STD_OUT)
!!    write(std_out,*) 'write out residue ', 5000+dcount
!!    do i=1,n
!!       write(5000+dcount,'(1P,10E15.7)') Gridwk%r(i),tmpPot%rvh(i),tmpPot%rvx(i),w(i),residue(i)
!!    enddo
 !!   close(5000+dcount);dcount=dcount+1   
 !!   if(needvtau.and.err>1.d-5) then
 !!       do i=1,n
 !!         write(500+fcount,'(1p,5e15.7)') Gridwk%r(i),w(i),residue(i)
 !!       enddo
 !!       close(500+fcount)
 !!       fcount=fcount+1
 !!       if (AC%CurIter==1000) then
 !!         write(std_out,*) 'Iteration 1000 -- try to stabilize'
 !!         w=tmpPot%rvh+tmpPot%rvx
 !!       endif       
 !!    endif

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
    !!! end testing
    ENDIF


    !write(std_out,*) 'in LDAGGAsub before end'; call flush_unit(std_out)
    CALL DestroyOrbit(tmpOrbit)
    CALL DestroyPot(tmpPot)
    DEALLOCATE (dum)

  END SUBROUTINE  LDAGGASub


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  DENITERSub	-- w is the electron density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DENITERSub(w,energy,residue,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: residue(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io,nw
    REAL(8),ALLOCATABLE :: dum(:)
    REAL(8) :: x,y
    INTEGER:: fcount=0,dcount=0
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot

    n=Gridwk%n
    nw=SIZE(w)
    if(n/=nw) then
      write(std_out,*) 'problem in DENITERsub n,nw', n,nw
      stop
    endif
    ALLOCATE(dum(nw),stat=i)
    IF (i/=0) THEN
        WRITE(6,*) 'Error in DENITERsub allocation' ,nw
        STOP
    ENDIF

    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    !  w enters subroutine as next iteration density
    !  It needs to be renormalized 
    !  tmpPot%rvh, tmpPot%rvx,  tmpPot%vtau need to be generated 

    x=integrator(Gridwk,w)
    write(std_out,*) 'In DENITERsub norm(m) adjust ', x,Potwk%q
    w=w*tmpPot%q/x
    call  poisson(Gridwk,x,w,tmpPot%rvh,y)
    write(std_out,*) 'after poisson  q,ecoul ',x,y
    call exch(Gridwk,w,tmpPot%rvx,x,y,tau=tmpOrbit%tau,vtau=tmpPot%vtau)
    write(std_out,*) 'after exch   exvct, eexc ',x,y
    tmpPot%rv=tmpPot%rvn+tmpPot%rvh+tmpPot%rvx
    tmpOrbit%den=w

    CALL Updatewfnwden(Gridwk,tmpPot,tmpOrbit,w,success)
    write(std_out,*) 'completed updatewfnwden with success ', success

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

    CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)

    write(std_out,*)  'Check tau ', integrator(Gridwk,tmpOrbit%tau)
    CALL Get_EXC(Gridwk,tmpPot,tmpOrbit,SCFwk)
    dum(1:n)=tmpOrbit%den(1:n)-w(1:n)
    residue=dum
    err=Dot_Product(residue,residue)
    write(STD_OUT,*) 'in DENITERSub   err ', err;call flush_unit(STD_OUT)
!!!    write(std_out,*) 'write out residue ', 5000+dcount
!!!    do i=1,n
!!!       write(5000+dcount,'(1P,10E15.7)') Gridwk%r(i),tmpPot%rvh(i),tmpPot%rvx(i),w(i),residue(i)
!!!    enddo
!!!    close(5000+dcount);dcount=dcount+1   

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
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx
       if (needvtau) Potwk%vtau=tmpPot%vtau
       Orbitwk%wfn=tmpOrbit%wfn
       If(diracrelativistic)Orbitwk%lwfn=tmpOrbit%lwfn
       Orbitwk%eig=tmpOrbit%eig
       !Orbitwk%den=tmpOrbit%den
       Orbitwk%den=w
       Orbitwk%otau=tmpOrbit%otau
       Orbitwk%tau=tmpOrbit%tau
       Orbitwk%deltatau=tmpOrbit%deltatau
       Call One_electron_energy_Report(Orbitwk,std_out)
!!!       !!! testing
!!!       if(needvtau.and.fcount==0) then
!!!      OPEN (unit=1001,file='potwithkedR2SCAN',form='formatted')
!!!      WRITE(1001,*) '#    r         rv               rvh           rvx       den    tau   vtau'      
!!!      DO i = 1,n
!!!       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
!!!&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i),Potwk%vtau(i)
!!!      ENDDO
!!!      CLOSE(1001)
!!!      fcount=fcount+1
!!!      endif
!!!    !!! end testing
    ENDIF


    CALL DestroyOrbit(tmpOrbit)
    CALL DestroyPot(tmpPot)
    DEALLOCATE (dum)

  END SUBROUTINE  DENITERSub



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
    REAL(8), ALLOCATABLE :: dum(:),dum1(:)
    INTEGER :: k,n

    n=Grid%n
    !write(std_out,*) 'In Get_EXC', n; call flush_unit(std_out)
    !write(std_out,*) 'In Get_EXC check ekin', integrator(Grid,Orbit%tau)
    CALL exch(Grid,Orbit%den,Pot%rvx,etxc,eex,&
&       tau=Orbit%tau,vtau=Pot%vtau)
    !write(std_out,*) 'After exch', etxc,eex; call flush_unit(std_out)

    SCF%eexc=eex
    etot = SCF%ekin+SCF%estatic+SCF%eexc
    SCF%etot=etot
    WRITE(STD_OUT,*) '    Total                    :  ',etot

    !write(std_out,*) 'before allocate'; call flush_unit(std_out)
    ALLOCATE(dum(n),dum1(n))
    dum=0;dum1=0
    !write(std_out,*) 'before dum'; call flush_unit(std_out)
    dum(2:n)=Pot%rvx(2:n)*Orbit%den(2:n)/Grid%r(2:n)
    if (needvtau) then
       dum=dum+Orbit%tau*Pot%vtau       !Kinetic energy correction    
       dum1=Orbit%tau*Pot%vtau
    endif   
    !write(std_out,*) 'after dum'; call flush_unit(std_out)
    !write(std_out,*) 'one',SCF%eone;call flush_unit(std_out)
    !write(std_out,*) 'coul',SCF%ecoul;call flush_unit(std_out)
    !write(std_out,*) 'exc',eex;call flush_unit(std_out)
    !write(std_out,*) 'DC',integrator(Grid,dum);call flush_unit(std_out)
    !write(std_out,*) 'vt',integrator(Grid,dum1);call flush_unit(std_out)
    WRITE(STD_OUT,*) '    Total   (DC form)        :  ',&
&        SCF%eone-SCF%ecoul+eex-integrator(Grid,dum)
    call flush_unit(std_out)
    DEALLOCATE(dum,dum1)
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
    n=Gridwk%n

     write(std_out,*) 'in Report ', stuff; call flush_unit(std_out)

    IF (frozencorecalculation) THEN
      OPEN (unit=1001,file='pot'//sub//TRIM(stuff),form='formatted')
      WRITE(1001,'(2a)') &
&         '#    r         rv               rvh           rvx       den',&
&         '                  tau          vtau          deltatau '   
      n=Gridwk%n
      DO i = 1,n
        WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&            Potwk%rvh(i),Potwk%rvx(i),Orbitwk%den(i), Orbitwk%tau(i),&
&            Potwk%vtau(i),Orbitwk%deltatau(i)
      ENDDO
      CLOSE(1001)
    ENDIF

    IF (.not.diracrelativistic) THEN
      OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
      WRITE(1001,'(a)') '#         r          wfn in order of s, p, d ... '
!      write(std_out,*) 'line 539 in ldagga ',Orbitwk%norbit;call flush_unit(std_out)
      DO i = 1,n
        WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&           (Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
!      write(std_out,*) 'line 544 in ldagga ',i;call flush_unit(std_out)
      ENDDO
      CLOSE(1001)
!      write(std_out,*) 'line 547 in ldagga ',Orbitwk%norbit;call flush_unit(std_out)
    ELSEIF (diracrelativistic) THEN
      OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
      WRITE(1001,'(a)') '#         r          wfn, lwfn in order of s, p-1/2, p+1/2 ... '
      DO i = 1,n
        WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&            (Orbitwk%wfn(i,j),Orbitwk%lwfn(i,j),j=1,Orbitwk%norbit)
      ENDDO
      CLOSE(1001)
    ENDIF
!      write(std_out,*) 'line 556 in ldagga ',Orbitwk%norbit;call flush_unit(std_out)

    counter=counter+1

  END SUBROUTINE Report_LDAGGA_functions

END MODULE ldagga_mod

