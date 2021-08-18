!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!     LDAGGA_SCF, LDAGGASub, Get_EXC, Get_FCEXC, Report_LDAGGA_functions
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

    INTEGER :: n,i
    REAL(8) :: en1,etxc,eex
    REAL(8), ALLOCATABLE :: arg(:)
    LOGICAL :: success

    Gridwk=>Gridin
    Orbitwk=>Orbitin
    Potwk=>Potin
    FCwk=>FCin
    SCFwk=>SCFin

    !write(std_out,*) 'in ldagga '; call flush_unit(std_out)
    n=Gridwk%n     ;!write(std_out,*) 'n ', n; call flush_unit(std_out)
    ALLOCATE(arg(n))

    !write(std_out,*) 'Before exch'; call flush_unit(std_out)
    CALL exch(Gridwk,Orbitwk%den,Potwk%rvx,etxc,eex,&
&       needvtau=Potwk%needvtau,tau=Orbitwk%tau,vtau=Potwk%vtau)
    
    !open(1001,file='testpot',form='formatted')
    !do i=1,n
    !   write(1001,'(1p,50e15.6)') Gridwk%r(i),Orbitwk%den(i),Potwk%rvx(i), &
&   !    Orbitwk%tau(i),Potwk%vtau(i)
    !enddo
    !close(1001)
    !stop
               

    Potwk%rv=Potwk%rvh+Potwk%rvx
    CALL zeropot(Gridwk,Potwk%rv,Potwk%v0,Potwk%v0p)
    Potwk%rv=Potwk%rv+Potwk%rvn

    !write(std_out,*) 'in ldagga before arg ',Potwk%v0,Potwk%v0p; 
    !call flush_unit(std_out)
    arg=Potwk%rv

    CALL InitAnderson_dr(AC,6,5,n,0.5d0,1.d3,1000,1.d-11,1.d-16,.true.)
    !write(std_out,*) 'in ldagga before Doand '; call flush_unit(std_out)
    !write(std_out,*) arg(1),arg(2)
    CALL DoAndersonMix(AC,arg,en1,LDAGGAsub,success)

    SCFwk%iter=AC%CurIter
    SCFwk%delta=AC%res

    CALL Report_LDAGGA_functions(scftype)

    CALL FreeAnderson(AC)
    WRITE(STD_OUT,*) 'Finished Anderson Mix', en1 ,' success = ', success
    DEALLOCATE(arg)
  END SUBROUTINE LDAGGA_SCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  LDAGGASub									!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE LDAGGASub(w,energy,residue,err,success,update)
    REAL(8), INTENT(INOUT) :: w(:)
    REAL(8), INTENT(OUT) :: energy
    REAL(8), INTENT(OUT) :: residue(:)
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    LOGICAL, INTENT(IN) :: update

    INTEGER :: i,j,k,n,io
    REAL(8),ALLOCATABLE :: dum(:)
    REAL(8) :: x
    INTEGER:: fcount=0
    TYPE (OrbitInfo) :: tmpOrbit
    TYPE (PotentialInfo) :: tmpPot

    ALLOCATE(dum(Gridwk%n),stat=i)
    IF (i/=0) THEN
       WRITE(STD_OUT,*) 'Error in LDAGGAsub allocation' ,Gridwk%n
       STOP
    ENDIF

    n=Gridwk%n

    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,w,success)
    tmpPot%rv=w
      If (.not.success) then   !  attempt to stablize solution
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
          tmpPot%rv(j+1)=(-1.d0+w(j+1)/2)
          do i=j+2,n
             tmpPot%rv(i)=-2.d0
          enddo
           write(std_out,*) 'Reset tmpPot ', j
           write(std_out,*) '   Last points '
              write(std_out,'(1p,20e15.7)') Gridwk%r(n), tmpPot%rv(n),w(n)

           CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,tmpPot%rv,success)
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
          ENDIF
       ENDDO
    ENDIF

    IF (.NOT.success) THEN
       WRITE(STD_OUT,*) 'Bad luck in Sub'
    ENDIF

    !write(std_out,*) 'in LDAGGAsub before Get'; call flush_unit(std_out)
    CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)
    !write(std_out,*) 'in LDAGGAsub before EXC'; call flush_unit(std_out)
    CALL Get_EXC(Gridwk,tmpPot,tmpOrbit,SCFwk)
    !write(std_out,*) 'after Get_EXC'; call flush_unit(std_out)
    dum=tmpPot%rvh+tmpPot%rvx+tmpPot%rvn-tmpPot%rv
    !write(std_out,*) 'after Get_EXC'; call flush_unit(std_out)
    residue=dum
    err=Dot_Product(residue,residue)
    w=tmpPot%rvh+tmpPot%rvx+tmpPot%rvn


    IF (update) THEN
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%otau=tmpOrbit%otau
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den
       Orbitwk%tau=tmpOrbit%tau
       If(diracrelativistic)Orbitwk%lwfn=tmpOrbit%lwfn

       Call One_electron_energy_Report(Orbitwk,std_out)
    ENDIF


    if(frozencorecalculation) then
         Call Get_FCKinCoul(Gridwk,tmpPot,tmpOrbit,FCwk,SCFwk)
         CALL Get_FCEXC(SCFwk)
         energy=SCFwk%evale
         CALL Total_FCEnergy_Report(SCFwk,std_out)
    else
         energy=SCFwk%etot
         CALL Total_Energy_Report(SCFwk,std_out)
    endif
    !write(std_out,*) 'in LDAGGAsub before end'; call flush_unit(std_out)
    CALL DestroyOrbit(tmpOrbit)
    DEALLOCATE (dum,tmpPot%rvn,tmpPot%rv,tmpPot%rvh,tmpPot%rvx)

  END SUBROUTINE  LDAGGASub


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
&       needvtau=Pot%needvtau,tau=Orbit%tau,vtau=Pot%vtau)
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
    INTEGER, ALLOCATABLE :: mmap(:)

    INTEGER, SAVE :: counter=0
    CHARACTER(4) :: stuff

    CALL mkname(counter,stuff)

     write(std_out,*) 'in Report ', stuff; call flush_unit(std_out)

    OPEN (unit=1001,file='pot'//sub//TRIM(stuff),form='formatted')
    n=Gridwk%n
    DO i = 1,n
      if (frozencorecalculation) then
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i)
      else
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i),Orbitwk%tau(i),Potwk%vtau(i)
      endif
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
    if (.not.diracrelativistic) then
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&           (Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
    ENDDO
    endif
    if (diracrelativistic) then
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&           (Orbitwk%wfn(i,j),Orbitwk%lwfn(i,j),j=1,Orbitwk%norbit)
    ENDDO
    endif

    CLOSE(1001)

    IF (ALLOCATED(mmap)) DEALLOCATE(mmap)

    counter=counter+1

  END SUBROUTINE Report_LDAGGA_functions

END MODULE ldagga_mod
