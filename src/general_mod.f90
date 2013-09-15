MODULE general_mod
  USE atomdata
  USE numerov_mod
  USE globalmath
  USE gridmod
  USE radialsr

  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  SUBROUTINE Updatewfn(Grid,Pot,Orbit,rvin,success)
  !      Given new potential rvin, generate new Orbit%wfn,Orbit%eig,Pot%den
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Updatewfn(Grid,Pot,Orbit,rvin,success)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    TYPE (PotentialInfo), INTENT(INOUT) :: Pot
    TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
    REAL(8), INTENT(IN) ::  rvin(:)
    LOGICAL :: success

    !  program to calculate wavefunctions given potential rvin

    INTEGER :: icount,i,j,k,n,it,start,np,ierr,nroot,s1,s2
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir,nzeff,jierr,nz
    REAL(8) :: h,emin,zz
    REAL(8), ALLOCATABLE :: dum(:)
    LOGICAL :: OK

    n=Grid%n; h=Grid%h;    nz=Pot%nz;   zz=Pot%zz
    success=.TRUE.

    allocate(dum(n))

    Pot%rv=rvin
    dum=rvin-Pot%rvn
    CALL zeropot(Grid,dum,Pot%v0,Pot%v0p)
    IF (ABS(Pot%v0)> 100.d0) Pot%v0=0
    IF (ABS(Pot%v0p)> 100.d0) Pot%v0p=0

    !  solve for bound states of Schroedinger equation
    !
    icount=0
    jierr=0
    it=0
    !  s states :
    IF (Orbit%nps.GT.0) THEN
       it=it+1
       emin=-nz*nz-0.1d0
       l=0
       nroot=Orbit%nps
       start=1;s1=start;s2=start+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  p states :
    IF (Orbit%npp.GT.1) THEN
       it=it+1
       emin=-nz*nz/4.d0-0.5d0
       l=1
       nroot=Orbit%npp-1
       start=start+Orbit%nps
       s1=start;s2=start+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  d states :
    IF (Orbit%npd.GT.2) THEN
       it=it+1
       emin=-nz*nz/9.d0-0.5d0
       l=2
       nroot=Orbit%npd-2
       start=start+Orbit%npp-1
       s1=start;s2=start+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  f states :
    IF (Orbit%npf.GT.3) THEN
       it=it+1
       emin=-nz*nz/16.d0-0.5d0
       l=3
       nroot=Orbit%npf-3
       start=start+Orbit%npd-2
       s1=start;s2=start+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  g states :
    IF (Orbit%npg.GT.4) THEN
       it=it+1
       emin=-nz*nz/25.d0-0.5d0
       l=4
       nroot=Orbit%npg-4
       start=start+Orbit%npf-3
       s1=start;s2=start+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
       IF (.NOT.OK) THEN
          success=.FALSE.
       ENDIF
    ENDIF

    deallocate(dum)

  END SUBROUTINE Updatewfn

  SUBROUTINE Get_KinCoul(Grid,Pot,Orbit,SCF,noalt)
    !  program to calculate Kinetic energy and Coulomb Energies from Orbit%wfn
    !   also update Pot%rvh
    TYPE(GridInfo), INTENT(INOUT) :: Grid
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot
    TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
    TYPE(SCFInfo), INTENT(INOUT) :: SCF
    LOGICAL, OPTIONAL :: noalt

    REAL(8) :: ecoul,eex,etot,ekin,eone,h,x,v0,qcal,etxc,small,rescale,electrons
    INTEGER :: icount,i,j,k,l,m,n,io
    REAL(8), ALLOCATABLE :: dum(:)
    REAL(8) :: small0=1.d-6
    INTEGER :: counter=1


    n=Grid%n; h=Grid%h
    small=small0

    ALLOCATE(dum(n),STAT=k)
    IF (k /= 0) THEN
       WRITE(6,*) 'Error in Get_KinCoul allocation ', n,k
       STOP
    ENDIF

    !update density
    Orbit%den(1:n)=0.d0
    DO io=1,Orbit%norbit
       IF (Orbit%occ(io).GT.small) THEN
          DO i=1,n
             IF (ABS(Orbit%wfn(i,io))<machine_zero)Orbit%wfn(i,io)=0
             Orbit%den(i)=Orbit%den(i)+ Orbit%occ(io)*(Orbit%wfn(i,io)**2)
          ENDDO
       ENDIF

    ENDDO
    qcal=integrator(Grid,Orbit%den)
    !WRITE(6,*) 'qcal = ', qcal

    ! for FC, Big change here , electrons are total , not valence anymore!
    electrons=Pot%q
    rescale=electrons/qcal
    Orbit%den(1:n)=Orbit%den(1:n)*rescale
    !WRITE(6,*) 'rescaled qcal = ', integrator(Grid,Orbit%den), Pot%q

    CALL poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul,v0)

    dum=0
    dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/Grid%r(2:n)
    SCF%estatic=integrator(Grid,dum)+ecoul

    !WRITE(6,*) ' n  l     occupancy       energy'
    ekin=0.0d0 ; if (frozencorecalculation) ekin=SCF%corekin
    eone=0.0d0
    DO io=1,Orbit%norbit
       if(.not.frozencorecalculation &
&          .or.frozencorecalculation.and.(.not.Orbit%iscore(io))) then
         !WRITE(6,'(i2,1x,i2,4x,1p,2e15.7)') &
         !&  Orbit%np(io),Orbit%l(io),&
         !&  Orbit%occ(io),Orbit%eig(io)
         eone=eone+Orbit%occ(io)*Orbit%eig(io)
         IF (counter>1.and..not.present(noalt)) THEN
            CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io),Pot%rv,x)
         ELSE
            CALL kinetic(Grid,Orbit%wfn(:,io),Orbit%l(io),x)
         ENDIF
         ekin=ekin+Orbit%occ(io)*x
       endif
    ENDDO

    SCF%eone=eone
    SCF%ekin=ekin
    SCF%ecoul=ecoul
    counter=counter+1

    DEALLOCATE(dum)
  END SUBROUTINE Get_KinCoul

  SUBROUTINE Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF,noalt)
    !  program to calculate Frozencore part of
    !    Kinetic energy and Coulomb Energies from Orbit%wfn
    !  Valence energy defined in paper --
    !    U. von Barth and C. D. Gelatt, Phys. Rev. B 21, 2222(1980)
    !  Exchange terms treated later
    TYPE(GridInfo), INTENT(INOUT) :: Grid
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot
    TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
    TYPE(FCInfo), INTENT(INOUT) :: FC
    TYPE(SCFInfo), INTENT(INOUT) :: SCF
    LOGICAL, OPTIONAL :: noalt

    REAL(8) :: tv,tc,x,y,electrons,vcoul,ccoul,vnucl,rescale
    INTEGER :: icount,i,j,k,l,m,n,io
    REAL(8), ALLOCATABLE :: dum(:)
    REAL(8) :: qcal,small,small0=1.d-6
    INTEGER :: firsttime=0

    n=Grid%n
    small=small0

    ALLOCATE(dum(n),STAT=k)
    IF (k /= 0) THEN
       WRITE(6,*) 'Error in Get_FCKinCoul allocation ', n,k
       STOP
    ENDIF


     write(6,*) 'In Get_FCKinCoul ', firsttime
    !update density  and calculated  kinetic energy
    Orbit%den(1:n)=0.d0
    FC%valeden(1:n)=0.d0     ;   tv=0 ; tc=0; SCF%eone=0
    DO io=1,Orbit%norbit
       IF (Orbit%occ(io).GT.small) THEN
          If(Orbit%iscore(io)) THEN
            If (firsttime==0) then
              If(present(noalt)) then
                  CALL kinetic(Grid,Orbit%wfn(:,io),Orbit%l(io),x)
              else
                  CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io),Pot%rv,x)
              endif
              tc=tc+Orbit%occ(io)*x
            Endif
          Else
              FC%valeden(:)=FC%valeden(:)+ Orbit%occ(io)*(Orbit%wfn(:,io)**2)
              If(present(noalt)) then
                  CALL kinetic(Grid,Orbit%wfn(:,io),Orbit%l(io),x)
              else
                  CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io),Pot%rv,x)
              endif
              tv=tv+Orbit%occ(io)*x
          ENDIF
          Orbit%den(:)=Orbit%den(:)+ Orbit%occ(io)*(Orbit%wfn(:,io)**2)
          SCF%eone=SCF%eone+Orbit%occ(io)*Orbit%eig(io)
       ENDIF
    ENDDO

    write(6,*) 'eone = ', SCF%eone

    if (firsttime==0) SCF%corekin=tc
    firsttime=1
    SCF%valekin=tv
    SCF%ekin=SCF%corekin+tv
    qcal=integrator(Grid,FC%valeden)
    !WRITE(6,*) 'qcalval = ', qcal
    Write(6,*) 'Core kin ', SCF%corekin,'  Vale kin  ', SCF%valekin

    electrons=FC%zvale
    rescale=electrons/qcal
    FC%valeden(1:n)=FC%valeden(1:n)*rescale
    !WRITE(6,*) 'rescaled qcalval = ', integrator(Grid,FC%valeden), FC%zvale

    x=FC%zvale
    CALL poisson(Grid,x,FC%valeden,dum,vcoul,y)  !valence-valence

    dum(2:n)=(dum(2:n)*FC%coreden(2:n)+Pot%rvn(2:n)*FC%valeden(2:n))/Grid%r(2:n)
    dum(1)=0.d0
    SCF%valecoul=vcoul+integrator(Grid,dum)

    !Update complete Hartree potential
    Call poisson(Grid,Pot%q,Orbit%den,Pot%rvh,x,Pot%v0)
        dum=0
    dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/Grid%r(2:n)
    SCF%estatic=integrator(Grid,dum)+x
    SCF%ecoul=x


    DEALLOCATE(dum)
  END SUBROUTINE Get_FCKinCoul

  SUBROUTINE Get_Nuclearpotential(Grid,Pot)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot

    INTEGER :: n
    REAL(8) :: h,q,v0,v0p
    REAL(8) :: r,RR
    INTEGER :: i,j,k

    IF (.NOT.finitenucleus) THEN
       DO i=1,Grid%n
          Pot%rvn(i)=-2*Pot%nz
       ENDDO
    ELSE
       RR=Pot%nz
       RR=2.9d-5*(RR**0.3333333333333333333333333d0)
       DO i=1,n
          Pot%rvn(i)=-2*Pot%nz*derf(Grid%r(i)/RR)
       ENDDO
    ENDIF
  END SUBROUTINE Get_Nuclearpotential

  SUBROUTINE ORTHONORMALIZE(Grid,Orbit)
    TYPE(GridInfo) ,INTENT(IN):: Grid
    TYPE(OrbitInfo) ,INTENT(INOUT):: Orbit

    Integer :: l,lmax,many,io,n
    REAL(8), ALLOCATABLE :: wfn(:,:)

    lmax=MAXVAL(Orbit%l)
    n=Grid%n
    allocate (wfn(n,Orbit%norbit))

    lmax=MAXVAL(Orbit%l)
    DO l=0,lmax
       wfn=0;many=0
       DO io=1,Orbit%norbit
          IF(Orbit%l(io)==l) THEN
            if (frozencorecalculation.AND.Orbit%iscore(io)) then
               many=many+1
               wfn(:,many)=Orbit%wfn(:,io)
            else
               CALL gramschmidt(Grid,many,wfn,Orbit%wfn(:,io))
               many=many+1
               Orbit%wfn(:,io)=wfn(:,many)
               CALL ADJUSTSIGN(Orbit%wfn(:,io),3)
            endif
          ENDIF
       ENDDO
    ENDDO

    deallocate (wfn)

  END SUBROUTINE ORTHONORMALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Subroutine to input f(:) and return f(:) if f(index)>=0 or
!                                return -f(:) if f(index)<0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ADJUSTSIGN(f,index)
     REAL(8),INTENT(INOUT) :: f(:)
     INTEGER,INTENT(IN) :: index

     IF (f(index)<0.d0) then
        f=-f
     ENDIF
  END SUBROUTINE ADJUSTSIGN

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Adjustnodes(Grid,Orbit,whichone)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE Adjustnodes(Grid,Orbit,whichone)
!    TYPE(GridInfo),INTENT(IN)::Grid
!    TYPE(OrbitInfo),INTENT(INOUT)::Orbit
!    INTEGER, INTENT(IN), OPTIONAL :: whichone
!
!    INTEGER :: io, i,j,n,num,nodes
!
!    n=Grid%n
!
!    if (PRESENT(whichone)) then
!       io=whichone
!       num=0;j=n+1
!       nodes=Orbit%np(io)-Orbit%l(io)-1
!       do i=3,n
!          if(Orbit%wfn(i-1,io)*Orbit%wfn(i,io)<0.d0) then
!            num=num+1
!            if (num>nodes) then
!               j=i
!               exit
!            endif
!          endif
!       enddo
!       If (j<n+1) then
!            Orbit%wfn(j:n,io)=0.d0
!            write(6,*) 'Adjusted orbit ', io,'  at r = ', &
!&                  Grid%r(j), Orbit%wfn(j-1,io)
!       endif
!     else
!       do io=1,Orbit%norbit
!          num=0;j=n+1
!          nodes=Orbit%np(io)-Orbit%l(io)-1
!          do i=3,n
!             if(Orbit%wfn(i-1,io)*Orbit%wfn(i,io)<0.d0) then
!               num=num+1
!               if (num>nodes) then
!                  j=i
!                  exit
!               endif
!             endif
!          enddo
!          If (j<n+1) Orbit%wfn(j:n,io)=0.d0
!       enddo
!     endif
!  END SUBROUTINE Adjustnodes
!

END  MODULE general_mod

