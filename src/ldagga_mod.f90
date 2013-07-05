MODULE ldagga_mod
  USE atomdata
  USE anderson_driver
  USE excor
  USE general_mod
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

    INTEGER :: n
    REAL(8) :: en1,etxc,eex
    REAL(8), ALLOCATABLE :: arg(:)
    LOGICAL :: success

    Gridwk=>Gridin
    Orbitwk=>Orbitin
    Potwk=>Potin
    FCwk=>FCin
    SCFwk=>SCFin

    !write(6,*) 'in ldagga '; call flush(6)
    n=Gridwk%n     ;write(6,*) 'n ', n; call flush(6)
    ALLOCATE(arg(n))

    CALL exch(Gridwk,Orbitwk%den,Potwk%rvx,etxc,eex)

    Potwk%rv=Potwk%rvh+Potwk%rvx
    CALL zeropot(Gridwk,Potwk%rv,Potwk%v0,Potwk%v0p)
    Potwk%rv=Potwk%rv+Potwk%rvn

    !write(6,*) 'in ldagga before arg '; call flush(6)
    arg=Potwk%rv

    CALL InitAnderson_dr(AC,6,5,n,0.5d0,1.d3,1000,1.d-11,1.d-16,.true.)
    !write(6,*) 'in ldagga before Doand '; call flush(6)
    CALL DoAndersonMix(AC,arg,en1,LDAGGAsub,success)

    SCFwk%iter=AC%CurIter
    SCFwk%delta=AC%res

    CALL Report_LDAGGA_functions(scftype)

    CALL FreeAnderson(AC)
    WRITE(6,*) 'Finished Anderson Mix', en1 ,' success = ', success
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

    write(6,*) 'in LDAGGAsub '; call flush(6)
    ALLOCATE(dum(Gridwk%n),stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in LDAGGAsub allocation' ,Gridwk%n
       STOP
    ENDIF

    n=Gridwk%n

    CALL CopyOrbit(Orbitwk,tmpOrbit)
    CALL CopyPot(Potwk,tmpPot)

    write(6,*) 'in LDAGGAsub before Updatewfn'; call flush(6)
    CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,w,success)
    tmpPot%rv=w
      If (.not.success) then   !  attempt to stablize solution
          write(6,*) 'Current eigs', (Orbitwk%eig(io),io=1,Orbitwk%norbit)
          j=n
          x=Orbitwk%eig(1)
          if (Orbitwk%norbit>1) then
             do io = 2, Orbitwk%norbit
                 if (Orbitwk%eig(io)<0.d0.and.x<Orbitwk%eig(io)) &
&                      x=Orbitwk%eig(io)
             enddo
          endif
          write(6,*) x,1.d0/sqrt(abs(x))
          x=1.d0/sqrt(abs(x))
          j=FindGridIndex(Gridwk,x)
          write(6,*) 'index', x,j,Gridwk%r(j)
          if (j<10)j=10
          if (j>n-10) j=n-10
          tmpPot%rv(j+1)=(-1.d0+w(j+1)/2)
          do i=j+2,n
             tmpPot%rv(i)=-2.d0
          enddo
           write(6,*) 'Reset tmpPot ', j
           !write(6,*) '   Last points '
           !   write(6,'(1p,20e15.7)') Gridwk%r(n), tmpPot%rv(n),w(n)

           CALL Updatewfn(Gridwk,tmpPot,tmpOrbit,tmpPot%rv,success)
       Endif

    !if FC core calc , restore core info backinto tmpOrbit

    IF(frozencorecalculation) THEN
       DO io = 1 , Orbitwk%norbit
          IF(Orbitwk%iscore(io)) THEN
             tmpOrbit%eig(io)=Orbitwk%eig(io)
             tmpOrbit%wfn(:,io)=Orbitwk%wfn(:,io)
          ENDIF
       ENDDO
    ENDIF

    IF (.NOT.success) THEN
       WRITE(6,*) 'Bad luck in Sub'
    ENDIF

    !write(6,*) 'in LDAGGAsub before Get'; call flush(6)
    CALL Get_KinCoul(Gridwk,tmpPot,tmpOrbit,SCFwk)
    CALL Get_EXC(Gridwk,tmpPot,tmpOrbit,SCFwk)
    dum=tmpPot%rvh+tmpPot%rvx+tmpPot%rvn-tmpPot%rv


    if(frozencorecalculation) then
         Call Get_FCKinCoul(Gridwk,tmpPot,tmpOrbit,FCwk,SCFwk)
         CALL Get_FCEXC(SCFwk)
         energy=SCFwk%evale
         CALL Total_FCEnergy_Report(SCFwk,6)
    else
         energy=SCFwk%etot
         CALL Total_Energy_Report(SCFwk,6)
    endif
    residue=dum
    err=Dot_Product(residue,residue)
    w=tmpPot%rv

    IF (update) THEN
       Potwk%rv=tmpPot%rv
       Potwk%rvh=tmpPot%rvh
       Potwk%rvx=tmpPot%rvx

       Orbitwk%wfn=tmpOrbit%wfn
       Orbitwk%eig=tmpOrbit%eig
       Orbitwk%den=tmpOrbit%den

       Call One_electron_energy_Report(Orbitwk,6)
    ENDIF

    !write(6,*) 'in LDAGGAsub before end'; call flush(6)
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
    TYPE(Gridinfo), INTENT(INOUT) :: Grid
    TYPE(Potentialinfo), INTENT(INOUT) :: Pot
    TYPE(Orbitinfo), INTENT(INOUT) :: Orbit
    TYPE(SCFInfo), INTENT(OUT) :: SCF

    REAL(8) :: eex,etot,etxc
    REAL(8), ALLOCATABLE :: dum(:)
    INTEGER :: k,n

    n=Grid%n
    CALL exch(Grid,Orbit%den,Pot%rvx,etxc,eex)

    SCF%eexc=eex
    etot = SCF%ekin+SCF%estatic+SCF%eexc
    SCF%etot=etot
    WRITE(6,*) '    Total                    :  ',etot

    ALLOCATE(dum(n))
    dum=0
    dum(2:n)=Pot%rvx(2:n)*Orbit%den(2:n)/Grid%r(2:n)
    WRITE(6,*) '    Total   (DC form)        :  ',&
&        SCF%eone-SCF%ecoul+eex-integrator(Grid,dum)
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

     write(6,*) 'in Report ', stuff; call flush(6)

    OPEN (unit=1001,file='pot'//sub//TRIM(stuff),form='formatted')
    n=Gridwk%n
    DO i = 1,n
      if (frozencorecalculation) then
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i)
      else
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i), &
&           Potwk%rvh(i),Potwk%rvx(i)
      endif
    ENDDO
    CLOSE(1001)
    OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
    DO i = 1,n
       WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&           (Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
    ENDDO
    CLOSE(1001)

    IF (ALLOCATED(mmap)) DEALLOCATE(mmap)

    counter=counter+1

  END SUBROUTINE Report_LDAGGA_functions

END MODULE ldagga_mod
