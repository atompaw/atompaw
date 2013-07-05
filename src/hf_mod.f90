MODULE hf_mod
  USE anderson_driver
  USE atomdata
  USE fock
  USE general_mod
  USE globalmath
  USE gridmod
  USE report_mod

  IMPLICIT NONE

  TYPE HFInfo
     REAL(8), POINTER :: lmbd(:,:),SumY(:,:),emin(:),emax(:)
     REAL(8), POINTER :: SumYCV(:,:),SumYVV(:,:),rVxcore(:)
     INTEGER, POINTER :: lmany(:),lmap(:,:)
     INTEGER, POINTER :: CSlmany(:),CSlmap(:,:)     !! closed shell blocks
     INTEGER  :: lmax
  ENDTYPE HFInfo

  TYPE(HFInfo), TARGET :: HF

  TYPE (Gridinfo) ,PRIVATE, POINTER :: Gridwk
  TYPE (PotentialInfo) ,PRIVATE, POINTER :: Potwk
  TYPE (OrbitInfo) ,PRIVATE, POINTER :: Orbitwk
  TYPE (OrbitInfo) ,PRIVATE, POINTER  :: EigOrbitwk
  TYPE (FCInfo) ,PRIVATE, POINTER :: FCwk
  TYPE (SCFInfo) ,PRIVATE, POINTER :: SCFwk

  LOGICAL, PRIVATE :: verboseoutput
  REAL(8), PRIVATE :: previouserror, currenterror
  REAL(8), PARAMETER, PRIVATE :: thrsh=1.d-7,tol=1.d-10

CONTAINS

  SUBROUTINE HF_tmp(Gridin,Orbitin,Potin,FCin,SCFin)
    TYPE(GridInfo) ,TARGET:: Gridin
    TYPE(OrbitInfo) ,TARGET:: Orbitin
    TYPE(PotentialInfo) ,TARGET:: Potin
    TYPE(FCInfo) ,TARGET:: FCin
    TYPE(SCFInfo) ,TARGET:: SCFin

    Gridwk=>Gridin
    Orbitwk=>Orbitin
    Potwk=>Potin
    FCwk=>FCin
    SCFwk=>SCFin

  END SUBROUTINE HF_tmp

  SUBROUTINE HF_SCF(scftype,lotsofoutput,Gridin,Orbitin,Potin,FCin,SCFin)
    CHARACTER(2), INTENT(IN) :: scftype
    LOGICAL, INTENT(IN) :: lotsofoutput
    TYPE(GridInfo) ,TARGET:: Gridin
    TYPE(OrbitInfo) ,TARGET:: Orbitin
    TYPE(PotentialInfo) ,TARGET:: Potin
    TYPE(FCInfo) ,TARGET:: FCin
    TYPE(SCFInfo) ,TARGET:: SCFin

    INTEGER :: i,n,io,icount,ok,loop,lng,l,many,np
    REAL(8) :: en1,x,ecoul,v0,err
    REAL(8) :: errv0=1.d-7,errv
    REAL(8), PARAMETER :: mix=0.4d0
    LOGICAL :: success,done
    INTEGER :: nrestart=5,mxloop1=20,mxloop2=100,mxloop=400
    INTEGER :: firsttime=0

    Gridwk=>Gridin
    Orbitwk=>Orbitin
    Potwk=>Potin
    FCwk=>FCin
    SCFwk=>SCFin

    verboseoutput=lotsofoutput

    n=Gridwk%n
    IF (firsttime<1) THEN
            write(6,*) 'in HF_SCF '; call flush(6);
       io=Orbitwk%norbit
       ALLOCATE(HF%lmbd(io,io),HF%SumY(n,io),HF%lmany(io), &
&           HF%lmap(io,io), HF%emin(io),HF%emax(io),   &
&           HF%CSlmany(io), HF%CSlmap(io,io),stat=ok)
       IF (ok /=0 ) THEN
          WRITE(6,*) 'Error in allocate HFdata ', io,ok
          STOP
       ENDIF
       HF%lmbd=0

       HF%lmax=0
       DO io=1,Orbitwk%norbit
          IF (Orbitwk%l(io)>HF%lmax) &
&              HF%lmax=Orbitwk%l(io)
          HF%emin(io)=-REAL(Potwk%nz**2)/Orbitwk%np(io)**2
          HF%emax(io)=-1.d-5
       ENDDO

       HF%lmap=0
       DO l=0,HF%lmax
          HF%lmany(l+1)=0
          DO io=1,Orbitwk%norbit
             IF (Orbitwk%l(io)==l) THEN
                HF%lmany(l+1)=HF%lmany(l+1)+1
                HF%lmap(l+1,HF%lmany(l+1))=io
             ENDIF
          ENDDO
       ENDDO

       HF%CSlmap=0
       DO l=0,HF%lmax
          HF%CSlmany(l+1)=0
          x=2*(2*l+1)
          DO io=1,Orbitwk%norbit
             IF (Orbitwk%l(io)==l.AND.ABS(x-Orbitwk%occ(io))<1.d-5) THEN
                HF%CSlmany(l+1)=HF%CSlmany(l+1)+1
                HF%CSlmap(l+1,HF%CSlmany(l+1))=io
             ENDIF
          ENDDO
       ENDDO

       firsttime=1
    ENDIF

    IF(TRIM(scftype)=='AE'.OR.TRIM(scftype)=='NC') THEN
       err=1.d10; currenterror=err
       DO loop=1,mxloop
          WRITE(6,*) '----loop---- ',loop
          !IF (err<1.d0) THEN
          !   DO io=1,Orbitwk%norbit
          !      CALL Adjustnodes(Gridwk,Orbitwk,io)
          !   ENDDO
          !ENDIF
          previouserror=currenterror
          CALL  HFIter(mix,tol,err,success,'ALL')
          IF (success) THEN
             WRITE(6,*) ' wfn iteration converged ', loop
             EXIT
          ENDIF
       ENDDO

       !DO io=1,Orbitwk%norbit
       !   CALL Adjustnodes(Gridwk,Orbitwk,io)
       !ENDDO

       currenterror=1.d10
       WRITE(6,*) 'Final  iterations '
       DO loop=1,mxloop
          previouserror=currenterror
          WRITE(6,*) '----Last loops---- ',loop
          CALL  HFIter(mix,tol,err,success,'CLO')
          IF (success) THEN
             WRITE(6,*) ' wfn iteration converged ', loop
             EXIT
          ENDIF
       ENDDO

       CALL GetLambda(Gridwk,Orbitwk,Potwk)
       CALL Total_Energy_Report(SCFwk,6)
       CALL One_electron_energy_Report(Orbitwk,6)
       CALL  Report_HF_functions(scftype)

       RETURN
    END IF

    IF (TRIM(scftype)=='SC') THEN
       IF(Orbitwk%exctype=='HFV') THEN
          io=Orbitwk%norbit
          ALLOCATE(HF%SumYVV(n,io),HF%SumYCV(n,io),HF%rVxcore(n))
          CALL HFVterms(Gridwk,Orbitwk,FCwk)
          CALL Get_FCEnergy_HFV(Gridwk,Potwk,Orbitwk,FCwk,SCFwk)
          CALL Total_FCEnergy_Report(SCFwk,6)
       ELSE
          CALL GetLambda(Gridwk,Orbitwk,Potwk,FCwk,SCFwk)
          CALL Total_Energy_Report(SCFwk,6)
          CALL Total_FCEnergy_Report(SCFwk,6)
          CALL One_electron_energy_Report(Orbitwk,6)
       ENDIF
       CALL  Report_HF_functions(scftype)

       RETURN
    ENDIF

    IF (TRIM(scftype)=='FC') THEN
       currenterror=1.d10;
       DO loop=1,mxloop
          WRITE(6,*) '----FC loops---- ',loop
          !DO io=1,Orbitwk%norbit
          !   IF (.NOT.Orbitwk%iscore(io)) CALL Adjustnodes(Gridwk,Orbitwk,io)
          !ENDDO
          previouserror=currenterror
          if (Orbitwk%exctype=='HF') then
             CALL  HFFCIter(mix,tol,err,success)
          elseif (Orbitwk%exctype=='HFV') then
             CALL  HFVFCIter(mix,tol,err,success)
          endif
          IF (success) THEN
             WRITE(6,*) ' wfn iteration converged ', loop
             EXIT
          ENDIF
       ENDDO
       DO loop=1,mxloop
          WRITE(6,*) '----FC loops again ---- ',loop
          previouserror=currenterror
          if (Orbitwk%exctype=='HF') then
             CALL  HFFCIter(mix,tol,err,success)
          elseif (Orbitwk%exctype=='HFV') then
             CALL  HFVFCIter(mix,tol,err,success)
          endif
          IF (success) THEN
             WRITE(6,*) ' wfn iteration converged ', loop
             EXIT
          ENDIF
       ENDDO
       !CALL GetLambda(Gridwk,Orbitwk,Potwk,FCwk,SCFwk)
       !CALL Total_Energy_Report(SCFwk,6)
       !CALL Total_FCEnergy_Report(SCFwk,6)
       !CALL One_electron_energy_Report(Orbitwk,6)
       CALL  Report_HF_functions(scftype)

       RETURN
    END IF


  END SUBROUTINE HF_SCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  HFIter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE HFIter(mix,tol,err,success,whichtype)
    REAL(8), INTENT(IN) :: mix,tol
    REAL(8), INTENT(OUT) :: err
    LOGICAL, INTENT(OUT) :: success
    CHARACTER(3), INTENT(IN) :: whichtype

    TYPE(GridInfo) ,POINTER:: Grid
    TYPE(OrbitInfo) ,POINTER:: Orbit
    TYPE(PotentialInfo) ,POINTER:: Pot


    INTEGER :: i,j,k,n,io,jo,lng,icount,many,lmax,l,norbit,iter,nodes,irmax
    REAL(8) :: ecoul,v0,accum,val,en,x
    REAL(8), POINTER :: r(:),rv(:),lmbd(:,:)
    REAL(8),ALLOCATABLE :: res(:),dum(:)
    INTEGER:: fcount=0
    INTEGER:: mxiter=20
    TYPE (OrbitInfo) :: tmpOrbit
    CHARACTER(4) :: stuff
    REAL(8), PARAMETER :: threshold=1.d-8, fixenergy=0.01d0,Rmax=10.d0

    Grid=>Gridwk
    Orbit=>Orbitwk
    Pot=>Potwk

    success=.FALSE.
    CALL ORTHONORMALIZE(Grid,Orbit)
    CALL CopyOrbit(Orbit,tmpOrbit)
    r=>Grid%r
    rv=>Pot%rv
    lmbd=>HF%lmbd
    norbit=Orbit%norbit
    n=Grid%n
    ALLOCATE(res(n),dum(n))
    irmax=FindGridIndex(Grid,Rmax)
    lmax=MAXVAL(Orbit%l)

    write(6,*) 'Before HFDiag '; call flush(6)
    CALL HFDiag(Grid,Orbit,Pot,whichtype)
    write(6,*) 'After HFDiag '; call flush(6)
    CALL Total_Energy_Report(SCFwk,6)
    CALL One_electron_energy_Report(Orbit,6)

    ! Solve inhomogeneous diffeq. and store result in tmpOrbit
    rv=Pot%rvn+Pot%rvh
    err=0
    DO io=1,norbit
       l=Orbit%l(io)
       res=0
       res(2:n)=-HF%SumY(2:n,io)/r(2:n)
       DO jo=1,norbit
          IF (jo==io.OR.(Orbit%l(jo)==l.AND.Orbit%occ(jo)>threshold)) THEN
             res=res+lmbd(jo,io)*Orbit%wfn(:,jo)
          ENDIF
       ENDDO
       CALL extrapolate(Grid,res)
       en=Orbit%eig(io)
       !IF((Orbit%occ(io)>threshold).AND.(en>HF%emax(io).OR.en<HF%emin(io)))THEN
       IF(en>HF%emax(io).OR.en<HF%emin(io))THEN
          en=(HF%emin(io)+HF%emax(io))/2
          if (io>1) then
             if (Orbit%l(io)==Orbit%l(io-1)) then
                 en=Orbit%eig(io-1)+0.1d0
             endif
          endif
          WRITE(6,*) 'wARNING -- en reset ', Orbit%eig(io),en
       ENDIF
       x=en
       do iter=1,mxiter
          dum=-(res-x*Orbit%wfn(:,io))
          tmpOrbit%wfn(:,io)=0.d0
          CALL inhomo_bound_numerov(Grid,l,n,x,rv,dum,tmpOrbit%wfn(:,io))
          j=countnodes(2,irmax,tmpOrbit%wfn(:,io))
          nodes=Orbit%np(io)-Orbit%l(io)-1
          if (x>0.or.nodes==j)exit
          if (j>nodes)x=x-fixenergy
          if (j<nodes)x=x+fixenergy
          write(6,*) 'Rerunning inhomo with x = ',&
&                 x,j,nodes
       enddo
       dum=(Orbit%wfn(:,io)-tmpOrbit%wfn(:,io))**2
       err=err+Orbit%occ(io)*Integrator(Grid,dum)
    ENDDO

    WRITE(6,*) 'COmpleted iter ', fcount,err
    SCFwk%delta=err; SCFwk%iter=fcount; currenterror=err

   ! CALL mkname(fcount,stuff)
   ! OPEN (unit=1001,file='hfwfn.'//TRIM(stuff),form='formatted')
   ! DO i=1,n
   !    WRITE(1001,'(1p,60E15.7)') Grid%r(i),(Orbit%wfn(i,io),&
   !&        tmpOrbit%wfn(i,io),io=1,norbit)
   ! ENDDO
   ! CLOSE(1001)

   ! OPEN (unit=1001,file='pot.'//TRIM(stuff),form='formatted')
   ! DO i=1,n
   !    WRITE(1001,'(1p,60E15.7)') Grid%r(i),Pot%rv(i),Pot%rvh(i),Orbit%den(i)
   ! ENDDO
   ! CLOSE(1001)

    ! update wfn if tolerance not satisfied
    IF (err>tol.AND.ABS(err-previouserror)>tol) THEN
       val=1.d0-mix
       WRITE(6,*) 'mixing wfns ', val
       DO io=1,norbit
          Orbit%wfn(:,io)=val*Orbit%wfn(:,io)+mix*tmpOrbit%wfn(:,io)
       ENDDO
       CALL ORTHONORMALIZE(Grid,Orbit)
    ELSE
       success=.TRUE.
    ENDIF

    fcount=fcount+1

    DEALLOCATE(res,dum)
    CALL  DestroyOrbit(tmpOrbit)

  END  SUBROUTINE HFIter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  On input, Orbit contains orthonormal wfns
  !    Orthonormality is not checked, but it is essential
  !     On output, wfns are "rotated" so that lambda is diagonal within closed
  !              shell blocks diagonal
  !      Pot%rvh, Pot$rv, and HF%SumY are updated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE HFDiag(Grid,Orbit,Pot,whichtype)
    TYPE(GridInfo) ,INTENT(IN):: Grid
    TYPE(OrbitInfo) ,INTENT(INOUT):: Orbit
    TYPE(PotentialInfo) ,INTENT(INOUT):: Pot
    CHARACTER(3), INTENT(IN) :: whichtype

    INTEGER :: i,j,k,n,io,jo,lng,icount,many,lmax,l,norbit,iter,lwork
    REAL(8) :: ecoul,v0,accum,val,en,x
    REAL(8), POINTER :: r(:),rv(:),lmbd(:,:)
    INTEGER, POINTER :: lmany(:),lmap(:,:)
    REAL(8),ALLOCATABLE :: A(:,:),w(:),work(:)
    INTEGER:: fcount=0
    TYPE (OrbitInfo) :: tmpOrbit
    CHARACTER(4) :: stuff

    r=>Grid%r
    rv=>Pot%rv
    lmbd=>HF%lmbd
    norbit=Orbit%norbit;lwork=norbit**2+10
    n=Grid%n
    CALL CopyOrbit(Orbit,tmpOrbit)


    ALLOCATE(A(norbit,norbit),w(norbit),work(lwork))

    IF (whichtype=="ALL") THEN
       lmany=>HF%lmany
       lmap=>HF%lmap
    ELSE IF (whichtype=="CLO") THEN
       lmany=>HF%CSlmany
       lmap=>HF%CSlmap
    ELSE
       WRITE(6,*) 'Error in HFDiag no matching type ', whichtype
       STOP
    ENDIF

    WRITE(6,*) ' Beginning HFDiag with fcount ', fcount
    CALL GetLambda(Grid,Orbit,Pot)

!!! now done in GetLambda
    !DO io=1,norbit
    !   Orbit%eig(io)=HF%lmbd(io,io)
    !ENDDO
    !  Change  blocks of lmbd:
    lmax=MAXVAL(Orbit%l)
    DO l=0,lmax
       A=0;w=0;work=0
       k=lmany(l+1)
       WRITE(6,*) 'l k ', l,k
       IF(k>1) THEN
          DO i=1,k
             DO j=1,k
                A(i,j)=0.5d0*lmbd(lmap(l+1,i),lmap(l+1,j)) &
&                    +0.5d0*lmbd(lmap(l+1,j),lmap(l+1,i))
             ENDDO
          ENDDO
          CALL dsyev('V','U',k,A,norbit,w,work,lwork,i)
          WRITE(6,*) 'returning form dsyev with info = ',i
          DO i=1,k
             WRITE(6,'(i5,2x,1p,15e15.7)') i,(A(i,j),j=1,k)
          ENDDO
          DO i=1,k
             io=lmap(l+1,i)
             WRITE(6,'("l,i,io, eig = ",3i5,2x,1p,5E15.7)')  l,i,io,w(i)
             Orbit%eig(io)=w(i)
             Orbit%wfn(:,io)=0.d0
             DO j=1,k
                jo=lmap(l+1,j)
                Orbit%wfn(:,io)=Orbit%wfn(:,io)+ &
&                    A(j,i)*tmpOrbit%wfn(:,jo)
             ENDDO
          ENDDO
       ENDIF
    ENDDO   ! l
    CALL ORTHONORMALIZE(Grid,Orbit)
    DO io=1,norbit
       CALL ADJUSTSIGN(Orbit%wfn(:,io),3)
    ENDDO

    CALL GetLambda(Grid,Orbit,Pot)

    WRITE(6,*) 'Checking orthonormality'
    DO io=1,norbit
       DO jo=1,norbit
          IF (Orbit%l(io)==Orbit%l(jo)) THEN
             WRITE(6,*) io,jo,&
&                 overlap(Grid,Orbit%wfn(:,io),Orbit%wfn(:,jo))
          ENDIF
       ENDDO
    ENDDO
    WRITE(6,*) ' Ending HFDiag with fcount ', fcount; call flush(6)

    fcount=fcount+1
    CALL DestroyOrbit(tmpOrbit)
    DEALLOCATE(A,w,work)

  END SUBROUTINE HFDiag

  SUBROUTINE GetLambda(Grid,Orbit,Pot,FC,SCF)
    TYPE(GridInfo) ,INTENT(IN):: Grid
    TYPE(OrbitInfo) ,INTENT(INOUT):: Orbit
    TYPE(PotentialInfo) ,INTENT(INOUT):: Pot
    !REAL(8), INTENT(OUT) :: SumY(:,:),lmbd(:,:)
    TYPE(FCInfo), INTENT(INOUT), OPTIONAL,TARGET :: FC
    TYPE(SCFInfo), INTENT(INOUT), OPTIONAL,TARGET :: SCF

    INTEGER :: i,j,k,n,io,jo,lng,icount,many,l,norbit
    REAL(8) :: ecoul,v0,accum,val,en,x,y,t1,t2,t3
    REAL(8), POINTER :: r(:),rv(:)
    REAL(8),ALLOCATABLE :: res(:),dum(:)
    INTEGER :: counter=0
    REAL(8), POINTER :: SumY(:,:),lmbd(:,:)
    TYPE(FCInfo), POINTER :: FCloc
    TYPE(SCFInfo), POINTER :: SCFloc

    SumY=>HF%SumY
    lmbd=>HF%lmbd

    IF(PRESENT(FC)) THEN
       FCloc=>FC
    ELSE
       FCloc=>FCwk
    ENDIF
    IF(PRESENT(SCF)) THEN
       SCFloc=>SCF
    ELSE
       SCFloc=>SCFwk
    ENDIF
    SCFloc%eone=0;   SCFloc%ekin=0; SCFloc%estatic=0;    SCFloc%ecoul=0
    SCFloc%eexc=0;   SCFloc%etot=0
    SCFloc%valekin=0; SCFloc%valecoul=0; SCFloc%valeexc=0;   SCFloc%evale=0
    SumY=0.d0; lmbd=0.d0
    r=>Grid%r
    norbit=Orbit%norbit
    n=Grid%n
    ALLOCATE(res(n),dum(n))

    !  recalculate electron density
    Orbit%den=0
    DO io=1,Orbit%norbit
       Orbit%den=Orbit%den+Orbit%occ(io)*Orbit%wfn(:,io)**2
    ENDDO

    IF (frozencorecalculation) THEN
       FCloc%valeden=0
       DO io=1,Orbit%norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             FCloc%valeden=FCloc%valeden+Orbit%occ(io)*Orbit%wfn(:,io)**2
          ENDIF
       ENDDO
    ENDIF


    dum=0; dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/r(2:n)
    SCFloc%estatic=integrator(Grid,dum)
    CALL poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul,v0) ;Pot%v0=v0
    SCFloc%ecoul=ecoul ; SCFloc%estatic=SCFloc%estatic+SCFloc%ecoul
    Pot%rv=Pot%rvn+Pot%rvh
    IF (frozencorecalculation) THEN
       CALL poisson(Grid,x,FCloc%valeden,dum,val,y)
       dum(1)=0;dum(2:n)=(Pot%rvn(2:n)*FCloc%valeden(2:n)+&
&           dum(2:n)*FCloc%coreden(2:n))/r(2:n)
       SCFloc%valecoul=integrator(Grid,dum)+val
    ENDIF

    lmbd=0
    IF(frozencorecalculation) THEN
       CALL Get_Energy_EXX_VC(Grid,Orbit,x)
       CALL Get_Energy_EXX_VV(Grid,Orbit,y)
       WRITE(6,*) '--EXX_VC = ',x,'--EXX_VV = ',y
       SCFloc%valeexc=x+y
       CALL Get_Energy_EXX(Grid,Orbit,x)
       SCFloc%eexc=x
    ELSE
       CALL Get_Energy_EXX(Grid,Orbit,x)
       SCFloc%eexc=x
    ENDIF

    WRITE(6,*) 'Checking Lambda '
    DO io=1,norbit
       CALL Calc_dexdphi_io(Grid,Orbit,io,sumY(:,io));l=Orbit%l(io)
       res(2:n)=sumY(2:n,io)/r(2:n); res(1)=0.d0
       DO jo=1,norbit
          IF (Orbit%l(io)==Orbit%l(jo)) THEN
             x=0
             CALL kinetic_ij(Grid,Orbit%wfn(:,io),Orbit%wfn(:,jo),&
&                 Orbit%l(io),y)
             x=x+y ;t1=y
             IF(io==jo) THEN
                SCFloc%ekin=SCFloc%ekin+Orbit%occ(io)*y
                IF(frozencorecalculation.AND..NOT.Orbit%iscore(io)) THEN
                   SCFloc%valekin=SCFloc%valekin+Orbit%occ(io)*y
                ENDIF
             ENDIF
             dum=Pot%rv*Orbit%wfn(:,io)*Orbit%wfn(:,jo)
             dum(2:n)=dum(2:n)/r(2:n); dum(1)=0
             x=x+integrator(Grid,dum); t2=integrator(Grid,dum)
             x=x+overlap(Grid,Orbit%wfn(:,jo),res)
             t3=overlap(Grid,Orbit%wfn(:,jo),res)
             WRITE(6,'("CHK lmbdD  ",3i5,2x,1p,5e15.7)') l,jo,io,x,t1,t2,t3
             lmbd(jo,io)=x
          ENDIF
       ENDDO
    ENDDO

    DO io=1,norbit
       DO jo=1,norbit
          WRITE(6,'("CHK lmbdagain ", 2i5,1Pe15.7)') io,jo,lmbd(io,jo)
       ENDDO
    ENDDO
    DO io=1,norbit
       Orbit%eig(io)=lmbd(io,io)
       SCFloc%eone=SCFloc%eone+Orbit%occ(io)*Orbit%eig(io)
    ENDDO

    SCFloc%etot=SCFloc%ekin+SCFloc%estatic+SCFloc%eexc
    SCFloc%evale=SCFloc%valekin+SCFloc%valecoul+SCFloc%valeexc
    DEALLOCATE(res,dum)
  END SUBROUTINE GetLambda

  SUBROUTINE hf_energy_only(Grid,Orbit,Pot,SCF)
    TYPE(GridInfo) ,INTENT(IN):: Grid
    TYPE(OrbitInfo) ,INTENT(IN):: Orbit
    TYPE(PotentialInfo) ,INTENT(INOUT):: Pot
    TYPE(SCFInfo), INTENT(INOUT) :: SCF


    INTEGER :: i,j,k,n,io,jo,lng,icount,many,l,norbit
    REAL(8) :: ecoul,v0,accum,val,en,x,y,t1,t2,t3
    REAL(8), POINTER :: r(:)
    REAL(8),ALLOCATABLE :: res(:),dum(:),den(:),coreden(:)

    SCF%eone=0;   SCF%ekin=0; SCF%estatic=0;    SCF%ecoul=0
    SCF%eexc=0;   SCF%etot=0
    SCF%valekin=0; SCF%valecoul=0; SCF%valeexc=0;   SCF%evale=0

    n=Grid%n
    r=>Grid%r
    allocate(res(n),dum(n),den(n),coreden(n))

    !  recalculate electron density
    den=0; coreden=0
    DO io=1,Orbit%norbit
       den=den+Orbit%occ(io)*Orbit%wfn(:,io)**2
       If (frozencorecalculation.AND.Orbit%iscore(io)) then
           coreden=coreden+Orbit%occ(io)*Orbit%wfn(:,io)**2
       endif
    ENDDO

    dum=0; dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/r(2:n)
    SCF%estatic=integrator(Grid,dum)
    CALL poisson(Grid,x,Orbit%den,dum,ecoul,v0)  ; Pot%v0=v0
    SCF%ecoul=ecoul ; SCF%estatic=SCF%estatic+SCF%ecoul
    IF (frozencorecalculation) THEN
       CALL poisson(Grid,x,coreden,dum,val,y)
       dum(1)=0;dum(2:n)=Pot%rvn(2:n)*coreden(2:n)/r(2:n)
       SCF%valecoul=SCF%estatic-(integrator(Grid,dum)+y)
    ENDIF

    IF(frozencorecalculation) THEN
       CALL Get_Energy_EXX_VC(Grid,Orbit,x)
       CALL Get_Energy_EXX_VV(Grid,Orbit,y)
       WRITE(6,*) '--EXX_VC = ',x,'--EXX_VV = ',y
       SCF%valeexc=x+y
       CALL Get_Energy_EXX(Grid,Orbit,x)
       SCF%eexc=x
    ELSE
       CALL Get_Energy_EXX(Grid,Orbit,x)
       SCF%eexc=x
    ENDIF

    Do io=1,Orbit%norbit
       l=Orbit%l(io)
       CALL kinetic_ij(Grid,Orbit%wfn(:,io),Orbit%wfn(:,io),l,y)
       SCF%ekin=SCF%ekin+Orbit%occ(io)*y
       if(frozencorecalculation.and..not.Orbit%iscore(io))  &
&            SCF%valekin=SCF%valekin+Orbit%occ(io)*y
    enddo

    CALL Total_Energy_Report(SCF,6)
    if (frozencorecalculation) CALL Total_FCEnergy_Report(SCF,6)

    deallocate(res,dum,den)
  END SUBROUTINE hf_energy_only

  SUBROUTINE GetLambdaHFV(Grid,Orbit,Pot,FC,SCF)
    TYPE(GridInfo) ,INTENT(IN):: Grid
    TYPE(OrbitInfo) ,INTENT(INOUT):: Orbit
    TYPE(PotentialInfo) ,INTENT(INOUT):: Pot
    !REAL(8), INTENT(OUT) :: SumY(:,:),lmbd(:,:)
    TYPE(FCInfo), INTENT(INOUT), OPTIONAL,TARGET :: FC
    TYPE(SCFInfo), INTENT(INOUT), OPTIONAL,TARGET :: SCF

    INTEGER :: i,j,k,n,io,jo,lng,icount,many,l,norbit
    REAL(8) :: ecoul,v0,accum,val,en,x,y,z,t1,t2,t3
    REAL(8), POINTER :: r(:),rv(:)
    REAL(8),ALLOCATABLE :: res(:),dum(:)
    INTEGER :: counter=0
    REAL(8), POINTER :: SumY(:,:),lmbd(:,:)
    TYPE(FCInfo), POINTER :: FCloc
    TYPE(SCFInfo), POINTER :: SCFloc

    lmbd=>HF%lmbd

    IF(PRESENT(FC)) THEN
       FCloc=>FC
    ELSE
       FCloc=>FCwk
    ENDIF
    IF(PRESENT(SCF)) THEN
       SCFloc=>SCF
    ELSE
       SCFloc=>SCFwk
    ENDIF
    SCFloc%eone=0;   SCFloc%ekin=0; SCFloc%estatic=0;    SCFloc%ecoul=0
    SCFloc%eexc=0;   SCFloc%etot=0
    SCFloc%valekin=0; SCFloc%valecoul=0; SCFloc%valeexc=0;   SCFloc%evale=0
    HF%SumYVV=0.d0; lmbd=0.d0
    r=>Grid%r
    norbit=Orbit%norbit
    n=Grid%n
    ALLOCATE(res(n),dum(n))

    !  recalculate electron density
    Orbit%den=0
    DO io=1,Orbit%norbit
       Orbit%den=Orbit%den+Orbit%occ(io)*Orbit%wfn(:,io)**2
    ENDDO

    !IF (frozencorecalculation) THEN    (assumed)
    FCloc%valeden=0
    DO io=1,Orbit%norbit
       IF (.NOT.Orbit%iscore(io)) THEN
          FCloc%valeden=FCloc%valeden+Orbit%occ(io)*Orbit%wfn(:,io)**2
       ENDIF
    ENDDO
    !ENDIF


    dum=0; dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/r(2:n)
    SCFloc%estatic=integrator(Grid,dum)
    CALL poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul,v0)  ; Pot%v0=v0
    SCFloc%ecoul=ecoul ; SCFloc%estatic=SCFloc%estatic+SCFloc%ecoul
    Pot%rv=Pot%rvn+Pot%rvh
    !IF (frozencorecalculation) THEN    (assumed)
    CALL poisson(Grid,x,FCloc%valeden,dum,val,y)
    dum(1)=0;dum(2:n)=(Pot%rvn(2:n)*FCloc%valeden(2:n)+&
&        dum(2:n)*FCloc%coreden(2:n))/r(2:n)
    SCFloc%valecoul=integrator(Grid,dum)+val
    !ENDIF

    lmbd=0
    !If(frozencorecalculation) then   (assumed)
    CALL Get_Energy_EXX_VC(Grid,Orbit,x)
    CALL Get_Energy_EXX_VV(Grid,Orbit,y)
    WRITE(6,*) '--EXX_VC = ',x,'--EXX_VV = ',y
    !SCFloc%valeexc=x+y
    dum=0; dum(2:n)=HF%rvxcore(2:n)/Grid%r(2:n)
    z=overlap(Grid,dum,FCloc%valeden)
    WRITE(6,*) '--EXX_VCV = ',z ; CALL flush(6)
    SCFloc%valeexc=y+z
    !call Get_Energy_EXX(Grid,Orbit,x)
    !SCFloc%eexc=x
    !else
    !  call Get_Energy_EXX(Grid,Orbit,x)
    !  SCFloc%eexc=x
    !endif

    WRITE(6,*) 'Checking Lambda '; CALL flush(6)
    DO io=1,norbit
       IF (.NOT.Orbit%iscore(io)) THEN
          CALL Calc_dexdphi_io_v(Grid,Orbit,io,HF%sumYVV(:,io));l=Orbit%l(io)
          res(2:n)=HF%sumYVV(2:n,io)/r(2:n); res(1)=0.d0
          DO jo=1,norbit
             IF (Orbit%l(io)==Orbit%l(jo)) THEN
                x=0
                CALL kinetic_ij(Grid,Orbit%wfn(:,io),Orbit%wfn(:,jo),l,y)
                x=x+y ;t1=y
                IF(io==jo) THEN
                   SCFloc%valekin=SCFloc%valekin+Orbit%occ(io)*y
                ENDIF
                dum=(Pot%rv+HF%rvxcore)*Orbit%wfn(:,io)*Orbit%wfn(:,jo)
                dum(2:n)=dum(2:n)/r(2:n); dum(1)=0
                x=x+integrator(Grid,dum); t2=integrator(Grid,dum)
                x=x+overlap(Grid,Orbit%wfn(:,jo),res)
                t3=overlap(Grid,Orbit%wfn(:,jo),res)
                WRITE(6,'("CHK lmbdV  ",3i5,2x,1p,5e15.7)') l,io,jo,x,t1,t2,t3
                call flush(6)
                lmbd(jo,io)=x
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    DO io=1,norbit
       DO jo=1,norbit
          WRITE(6,'("CHK lmbdagain ", 2i5,1Pe15.7)') io,jo,lmbd(io,jo)
          CALL flush(6)
       ENDDO
    ENDDO
    DO io=1,norbit
       IF (.NOT.Orbit%iscore(io)) THEN
          Orbit%eig(io)=lmbd(io,io)
          SCFloc%eone=SCFloc%eone+Orbit%occ(io)*Orbit%eig(io)
       ENDIF
    ENDDO
       !SCFloc%etot=SCFloc%ekin+SCFloc%estatic+SCFloc%eexc
    SCFloc%evale=SCFloc%valekin+SCFloc%valecoul+SCFloc%valeexc
    DEALLOCATE(res,dum)
    END SUBROUTINE GetLambdaHFV

     SUBROUTINE HFVterms(Grid,Orbit,FC)
       TYPE(GridInfo), INTENT(IN) :: Grid
       TYPE(OrbitInfo), INTENT(IN) :: Orbit
       TYPE(FCInfo), INTENT(IN) :: FC

       INTEGER :: i,j,k,n,io,norbit

       norbit=Orbit%norbit
       n=Grid%n

       HF%SumYVV=0;HF%SumYCV=0;HF%rVxcore=0
       DO io=1,norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             CALL Calc_dexdphi_io_v(Grid,Orbit,io,HF%sumYVV(:,io))
             CALL Calc_dexdphi_io_c(Grid,Orbit,io,HF%sumYCV(:,io))
             HF%rvxcore=HF%rvxcore+Orbit%occ(io)*HF%sumYCV(:,io)*Orbit%wfn(:,io)
          ENDIF
       ENDDO

       HF%rvxcore(2:n)= HF%rvxcore(2:n)/FC%valeden(2:n)
       HF%rvxcore(1)=0
     END SUBROUTINE HFVterms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  HFFCIter   -- Frozen core version;
!!!!   assume not more than 1 valence shell as actually
!!!!           close shell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE HFFCIter(mix,tol,err,success)
       REAL(8), INTENT(IN) :: mix,tol
       REAL(8), INTENT(OUT) :: err
       LOGICAL, INTENT(OUT) :: success

       TYPE(GridInfo) ,POINTER:: Grid
       TYPE(OrbitInfo) ,POINTER:: Orbit
       TYPE(PotentialInfo) ,POINTER:: Pot


       INTEGER :: i,j,k,n,io,jo,lng,icount,many,lmax,l,norbit,iter,nodes
       REAL(8) :: ecoul,v0,accum,val,en,x
       REAL(8), POINTER :: r(:),rv(:),lmbd(:,:)
       REAL(8),ALLOCATABLE :: res(:),dum(:)
       INTEGER:: fcount=0
       INTEGER:: mxiter=20
       TYPE (OrbitInfo) :: tmpOrbit
       CHARACTER(4) :: stuff
       REAL(8), PARAMETER :: threshold=1.d-8

       IF (.NOT.frozencorecalculation) THEN
          WRITE(6,*) 'Error -- subroutine HFFCIter is only for Frozen core '
          STOP
       ENDIF

       Grid=>Gridwk
       Orbit=>Orbitwk
       Pot=>Potwk

       success=.FALSE.
       CALL ORTHONORMALIZE(Grid,Orbit)
       CALL CopyOrbit(Orbit,tmpOrbit)
       r=>Grid%r
       rv=>Pot%rv
       lmbd=>HF%lmbd
       norbit=Orbit%norbit
       n=Grid%n
       ALLOCATE(res(n),dum(n))

       lmax=MAXVAL(Orbit%l)

       CALL GetLambda(Grid,Orbit,Pot)
       CALL Total_FCEnergy_Report(SCFwk,6)
       CALL Total_Energy_Report(SCFwk,6)
       CALL One_electron_energy_Report(Orbit,6)

       ! Solve inhomogeneous diffeq. and store result in tmpOrbit
       rv=Pot%rvn+Pot%rvh
       err=0
       DO io=1,norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             Orbit%eig(io)=lmbd(io,io)
             l=Orbit%l(io)
             res=0
             res(2:n)=-HF%SumY(2:n,io)/r(2:n)
             DO jo=1,norbit
                IF (jo==io.OR.(Orbit%l(jo)==l.AND.Orbit%occ(jo)>threshold)) THEN
                   res=res+lmbd(jo,io)*Orbit%wfn(:,jo)
                ENDIF
             ENDDO
             CALL extrapolate(Grid,res)
             en=Orbit%eig(io)
             IF ((Orbit%occ(io)>threshold).OR.(en>HF%emax(io).OR.en<HF%emin(io))) THEN
                en=(HF%emin(io)+HF%emax(io))/2
             ENDIF
             dum=-(res-en*Orbit%wfn(:,io))
             tmpOrbit%wfn(:,io)=0.d0
             CALL inhomo_bound_numerov(Grid,l,n,en,rv,dum,tmpOrbit%wfn(:,io))
             dum=(Orbit%wfn(:,io)-tmpOrbit%wfn(:,io))**2
             err=err+Orbit%occ(io)*Integrator(Grid,dum)
          ENDIF
       ENDDO

       WRITE(6,*) 'COmpleted iter ', fcount,err
       SCFwk%delta=err; SCFwk%iter=fcount;currenterror=err

      ! CALL mkname(fcount,stuff)
      ! OPEN (unit=1001,file='FChfwfn.'//TRIM(stuff),form='formatted')
      ! DO i=1,n
      !    WRITE(1001,'(1p,60E15.7)') Grid%r(i),(Orbit%wfn(i,io),&
      !&        tmpOrbit%wfn(i,io),io=1,norbit)
      ! ENDDO
      ! CLOSE(1001)

      ! OPEN (unit=1001,file='FCpot.'//TRIM(stuff),form='formatted')
      ! DO i=1,n
      !    WRITE(1001,'(1p,60E15.7)') Grid%r(i),Pot%rv(i),Pot%rvh(i),Orbit%den(i)
      ! ENDDO
      ! CLOSE(1001)

       ! update wfn if tolerance not satisfied
       IF (err>tol.AND.ABS(err-previouserror)>tol) THEN
          val=(1.d0-mix)
          WRITE(6,*) 'mixing wfns ', val
          DO io=1,norbit
             IF (.NOT.Orbit%iscore(io)) THEN
                Orbit%wfn(:,io)=val*Orbit%wfn(:,io)+mix*tmpOrbit%wfn(:,io)
             ENDIF
          ENDDO
          CALL ORTHONORMALIZE(Grid,Orbit)
       ELSE
          success=.TRUE.
       ENDIF

       fcount=fcount+1

       DEALLOCATE(res,dum)
       CALL  DestroyOrbit(tmpOrbit)

     END  SUBROUTINE HFFCIter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  HFVFCIter   -- Frozen core version with potential;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE HFVFCIter(mix,tol,err,success)
       REAL(8), INTENT(IN) :: mix,tol
       REAL(8), INTENT(OUT) :: err
       LOGICAL, INTENT(OUT) :: success

       TYPE(GridInfo) ,POINTER:: Grid
       TYPE(OrbitInfo) ,POINTER:: Orbit
       TYPE(PotentialInfo) ,POINTER:: Pot


       INTEGER :: i,j,k,n,io,jo,lng,icount,many,lmax,l,norbit,iter,nodes
       REAL(8) :: ecoul,v0,accum,val,en,x
       REAL(8), POINTER :: r(:),rv(:),lmbd(:,:)
       REAL(8),ALLOCATABLE :: res(:),dum(:)
       INTEGER:: fcount=0
       INTEGER:: mxiter=20
       TYPE (OrbitInfo) :: tmpOrbit
       CHARACTER(4) :: stuff
       REAL(8), PARAMETER :: threshold=1.d-8

       Grid=>Gridwk
       Orbit=>Orbitwk
       Pot=>Potwk

       IF (.NOT.frozencorecalculation.OR.Orbit%exctype/='HFV') THEN
          WRITE(6,*) 'Error -- subroutine HFVFCIter is only for Frozen core '
          STOP
       ENDIF


       success=.FALSE.
       CALL ORTHONORMALIZE(Grid,Orbit)
       CALL CopyOrbit(Orbit,tmpOrbit)
       r=>Grid%r
       rv=>Pot%rv
       lmbd=>HF%lmbd
       norbit=Orbit%norbit
       n=Grid%n
       ALLOCATE(res(n),dum(n))

       lmax=MAXVAL(Orbit%l)


       CALL GetLambdaHFV(Grid,Orbit,Pot)
       CALL Total_FCEnergy_Report(SCFwk,6)
       !CALL Total_Energy_Report(SCFwk,6)
       CALL One_electron_energy_Report(Orbit,6)

       ! Solve inhomogeneous diffeq. and store result in tmpOrbit
       rv=Pot%rvn+Pot%rvh + HF%rvxcore
       err=0
       DO io=1,norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             Orbit%eig(io)=lmbd(io,io)
             l=Orbit%l(io)
             res=0
             res(2:n)=-HF%SumYVV(2:n,io)/r(2:n)
             DO jo=1,norbit
                IF (jo==io.OR.(Orbit%l(jo)==l.AND.Orbit%occ(jo)>threshold)) THEN
                   res=res+lmbd(jo,io)*Orbit%wfn(:,jo)
                ENDIF
             ENDDO
             CALL extrapolate(Grid,res)
             en=Orbit%eig(io)
             IF ((Orbit%occ(io)>threshold).OR.(en>HF%emax(io).OR.en<HF%emin(io))) THEN
                en=(HF%emin(io)+HF%emax(io))/2
             ENDIF
             dum=-(res-en*Orbit%wfn(:,io))
             tmpOrbit%wfn(:,io)=0.d0
             CALL inhomo_bound_numerov(Grid,l,n,en,rv,dum,tmpOrbit%wfn(:,io))
             dum=(Orbit%wfn(:,io)-tmpOrbit%wfn(:,io))**2
             err=err+Orbit%occ(io)*Integrator(Grid,dum)
          ENDIF
       ENDDO

       WRITE(6,*) 'COmpleted iter ', fcount,err
       SCFwk%delta=err; SCFwk%iter=fcount;currenterror=err

      ! CALL mkname(fcount,stuff)
      ! OPEN (unit=1001,file='FChfwfn.'//TRIM(stuff),form='formatted')
      ! DO i=1,n
      !    WRITE(1001,'(1p,60E15.7)') Grid%r(i),(Orbit%wfn(i,io),&
      !&        tmpOrbit%wfn(i,io),io=1,norbit)
      ! ENDDO
      ! CLOSE(1001)

      ! OPEN (unit=1001,file='FCpot.'//TRIM(stuff),form='formatted')
      ! DO i=1,n
      !    WRITE(1001,'(1p,60E15.7)') Grid%r(i),Pot%rv(i),Pot%rvh(i),&
      !&        HF%rvxcore(i),Orbit%den(i)
      ! ENDDO
      ! CLOSE(1001)

       ! update wfn if tolerance not satisfied
       IF (err>tol.AND.ABS(err-previouserror)>tol) THEN
          val=(1.d0-mix)
          WRITE(6,*) 'mixing wfns ', val
          DO io=1,norbit
             IF (.NOT.Orbit%iscore(io)) THEN
                Orbit%wfn(:,io)=val*Orbit%wfn(:,io)+mix*tmpOrbit%wfn(:,io)
             ENDIF
          ENDDO
          CALL ORTHONORMALIZE(Grid,Orbit)
       ELSE
          success=.TRUE.
       ENDIF

       fcount=fcount+1

       DEALLOCATE(res,dum)
       CALL  DestroyOrbit(tmpOrbit)

     END  SUBROUTINE HFVFCIter

     SUBROUTINE Get_FCEnergy_HF(Grid,Pot,Orbit,FC,SCF)
       TYPE(GridInfo), INTENT(INOUT) :: Grid
       TYPE(PotentialInfo), INTENT(INOUT) :: Pot
       TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
       TYPE(FCInfo), INTENT(INOUT) :: FC
       TYPE(SCFInfo), INTENT(INOUT) :: SCF

       CALL GetLambda(Grid,Orbit,Pot,FC,SCF)

     END SUBROUTINE Get_FCEnergy_HF

     SUBROUTINE Get_FCEnergy_HFV(Grid,Pot,Orbit,FC,SCF)
       TYPE(GridInfo), INTENT(INOUT) :: Grid
       TYPE(PotentialInfo), INTENT(INOUT) :: Pot
       TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
       TYPE(FCInfo), INTENT(INOUT) :: FC
       TYPE(SCFInfo), INTENT(INOUT) :: SCF

       REAL(8), ALLOCATABLE :: dum(:)
       REAL(8) :: x,y,z,val
       INTEGER :: i,j,k,n,norbit,io

       n=Grid%n
       ALLOCATE(dum(n))
       CALL poisson(Grid,x,FC%valeden,dum,val,y)
       dum(1)=0;dum(2:n)=(Pot%rvn(2:n)*FC%valeden(2:n)+&
&           dum(2:n)*FC%coreden(2:n))/Grid%r(2:n)
       SCF%valecoul=integrator(Grid,dum)+val
       CALL Get_Energy_EXX_VC(Grid,Orbit,x)
       CALL Get_Energy_EXX_VV(Grid,Orbit,y)
       WRITE(6,*) '--EXX_VC = ',x,'--EXX_VV = ',y    ; CALL flush(6)
       dum=0; dum(2:n)=HF%rvxcore(2:n)/Grid%r(2:n)
       z=overlap(Grid,dum,FC%valeden)
       WRITE(6,*) '--EXX_VCV = ',z ; CALL flush(6)
       SCF%valeexc=y+z

       norbit=Orbit%norbit; SCF%valekin=0
       DO io=1,norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             CALL kinetic_ij(Grid,Orbit%wfn(:,io),Orbit%wfn(:,io),&
&                 Orbit%l(io),y)
             SCF%valekin=SCF%valekin+Orbit%occ(io)*y
          ENDIF
       ENDDO

       DEALLOCATE(dum)
     END SUBROUTINE Get_FCEnergy_HFV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    HFunocc
!      solve HF equation for unoccupied continuum state
!            No Lagrange multipliers used
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE HFunocc(Grid,Orbit,l,energy,rv,v0,v0p,wfn,lng,Xout)
       TYPE(GridInfo) , INTENT(IN) :: Grid
       TYPE(OrbitInfo) , INTENT(IN) :: Orbit
       INTEGER, INTENT(IN) :: l
       REAL(8), INTENT(IN) :: energy
       REAL(8), INTENT(IN) :: rv(:),v0,v0p
       REAL(8), INTENT(INOUT) :: wfn(:)        ! unbound wfn
       INTEGER, INTENT(OUT) :: lng
       REAL(8), OPTIONAL, INTENT(OUT) :: Xout(:)

       INTEGER :: i,j,k,n,io,jo,many,lmax,norbit,iter,nodes,flag,ZZ
       INTEGER :: istart,imtch
       REAL(8) :: ecoul,accum,val,en,xx,A,B,C,D,Coeff,kappa,jl,nl,jlp,nlp,norm
       REAL(8) :: jjl,nnl,jjlp,nnlp,previous,residue,range,r0,c0,c1,c2,scale
       REAL(8), POINTER :: r(:)
       REAL(8),ALLOCATABLE :: wfn0(:),ve(:),X(:),dum(:)
       INTEGER:: fcount=0
       INTEGER:: mxiter=2000,ni=999
       REAL(8), PARAMETER :: threshold=1.d-8, tol=1.d-8,mix=0.5d0,rch=1.d0
       LOGICAL :: success
       INTEGER :: icount=0

       IF (energy<0.d0) THEN
          WRITE(6,*) 'Error in HFunocc -- assumes energy>0 ',energy
          STOP
       ENDIF

       write(6,*) 'Entering HFunocc with energy l = ', energy,l
       n=Grid%n
       r=>Grid%r
       r0=Grid%drdu(1)

       range=0;lng=1
       norbit=Orbit%norbit
       do io=1,norbit
            !Call Filter(n,Orbit%wfn(:,io),machine_zero)
               !!! done in calling program
            k=Maxloc(abs(Orbit%wfn(:,io)),1)
            do i=n,k,-1
               j=i
               if (abs(Orbit%wfn(i,io))>1.d-5) exit
            enddo
            if (j>lng) lng=j
            range=Grid%r(lng)
       end do
       write(6,*) 'HFunocc: range ', lng, range
       write(6,*) 'Calculating HFunocc for range ', range, lng

       imtch=1
       do i=1,lng
          imtch=i
          if(Grid%r(i)>rch) exit
       enddo
       write(6,*) 'Calculating HFunocc for imtch ', rch , &
&            Grid%r(imtch), imtch

       ALLOCATE(wfn0(n),ve(n),X(n),dum(n))

       istart=6

       ZZ=(-rv(1)+0.001d0)/2
       write(6,*) 'ZZ = ', ZZ
       ve(2:n)=(rv(2:n)+2*ZZ)/Grid%r(2:n)
       CALL extrapolate(Grid,ve)
       write(6,*) 've ', ve(1),ve(2),ve(3), ve(4),ve(5)
       wfn0=0
       wfn0(1:lng)=((Grid%r(1:lng)/r0)**(l+1))*(1-&
&            (float(ZZ)/(l+1))*Grid%r(1:lng))
       wfn0(1:lng)=wfn0(1:lng)/wfn0(lng)

       Call Calc_Xp(Grid,Orbit,ni,l,wfn0,X,lng)
       X(2:lng)=X(2:lng)/Grid%r(2:lng)
       Call extrapolate(Grid,X)
       !write(6,*) 'X ', X(1),X(2),X(3), X(4),X(5)

       dum(2:lng)=X(2:lng)/(Grid%r(2:lng)**(l+1))
       call extrapolate(Grid,dum)
       write(6,*) 'dum ', dum(1),dum(2),dum(3), dum(4),dum(5)

       wfn=0
       do i=2,istart
          xx=Grid%r(i);c0=1.d0;c1=-float(ZZ)/(l+1)
           c2=(dum(1)+ve(1)-energy+2*float(ZZ*ZZ)/(l+1))/(4*l+6)
          wfn(i)=((xx)**(l+1))*(c0+xx*(c1+xx*c2))
       enddo
       scale=wfn0(istart)/wfn(istart)
       wfn(1:istart)=wfn(1:istart)*scale

       !!!!Call Milne(Grid,istart,energy,l,ZZ,ve,X,wfn)
       CALL midrange_numerov(Grid,l,istart,lng,energy,rv,X,wfn)
       wfn0(1:lng)=wfn(1:lng)/wfn(lng)

       success=.false.
       Do iter=1,mxiter
          Call Calc_Xp(Grid,Orbit,ni,l,wfn0,X,lng)
          X(2:lng)=X(2:lng)/Grid%r(2:lng)
          Call extrapolate(Grid,X)
          ! write(6,*) 'X ', X(1),X(2),X(3), X(4),X(5)
          !Call Milne(Grid,istart,energy,l,ZZ,ve,X,wfn)
          dum(2:lng)=X(2:lng)/(Grid%r(2:lng)**(l+1))
          call extrapolate(Grid,dum)
           write(6,*) 'dum ', dum(1),dum(2),dum(3), dum(4),dum(5)
          wfn=0
          do i=2,istart
             xx=Grid%r(i);c0=1.d0;c1=-float(ZZ)/(l+1)
              c2=(dum(1)+ve(1)-energy+2*float(ZZ*ZZ)/(l+1))/(4*l+6)
             wfn(i)=((xx/r0)**(l+1))*(c0+xx*(c1+xx*c2))
          enddo
          scale=wfn0(istart)/wfn(istart)
          wfn(1:istart)=wfn(1:istart)*scale
          CALL midrange_numerov(Grid,l,istart,lng,energy,rv,X,wfn)
          !wfn(1:lng)=wfn(1:lng)/wfn(lng)
          !do i=1,lng
          !   write(600+iter,'(1p,8e16.7)') Grid%r(i),wfn(i),wfn0(i),X(i)
          !enddo
          residue=sum(abs(wfn(1:lng)-wfn0(1:lng)))
          write(6,*) 'iter ', iter,residue
          if (residue<tol) then
                  success=.true.
                  exit
          endif
          wfn0(1:lng)=wfn0(1:lng)*(1.d0-mix)+wfn(1:lng)*(mix)
          wfn0=wfn0/wfn0(lng)
       enddo


       IF (.NOT.success) THEN
          WRITE(6,*) 'Warning in HFunocc not converged ',l,energy,range,residue
          WRITE(6,*)  'stopping after dump ' , 800+icount
          do i=1,lng
           write(800+icount,'(1P3e15.7)')Grid%r(i),wfn0(i),wfn(i)
          enddo
          stop

       ENDIF

       !do i=2,lng
       !    write(800+icount,'(1P9e15.7)')Grid%r(i),wfn0(i),wfn(i),X(i), &
       !&        -Gsecondderiv(Grid,i,wfn0)+&
       !&         (l*(l+1)/r(i)**2+rv(i)/r(i)-energy)*wfn0(i)+X(i)
       !enddo

       wfn=0
       wfn(1:lng)=wfn0(1:lng)
       icount=icount+1

       If (present(Xout)) then
               Xout=0
               Xout(1:lng)=X(1:lng)
       ENdif

       DEALLOCATE(wfn0,ve,X)

     END  SUBROUTINE HFunocc

     SUBROUTINE Report_HF_functions(sub)
       CHARACTER(2) :: sub
       INTEGER :: i,j,k,n,nmap,io,nc,nv
       INTEGER, ALLOCATABLE :: vmap(:),cmap(:)

       INTEGER, SAVE :: counter=0
       CHARACTER(4) :: stuff

       !!  store variables in Orbit datastructure
       Orbitwk%X=HF%SumY
       Orbitwk%lqp=HF%lmbd
       WRITE(6,*)
       WRITE(6,*) ' Summary of HF orbitals '
       WRITE(6,"(' n  l     occupancy            energy      ')")
       DO io=1,Orbitwk%norbit
          WRITE(6,'(i2,1x,i2,4x,1p,3e15.7)') &
&              Orbitwk%np(io),Orbitwk%l(io),Orbitwk%occ(io),Orbitwk%eig(io)
       ENDDO

       CALL mkname(counter,stuff)

       n=Gridwk%n
       IF (frozencorecalculation) THEN
          ALLOCATE(vmap(Orbitwk%norbit), cmap(Orbitwk%norbit))
          j=0; k=0
          DO i=1,Orbitwk%norbit
             IF (.NOT.Orbitwk%iscore(i)) THEN
                j=j+1
                vmap(j)=i
             ELSE
                k=k+1
                cmap(k)=i
             ENDIF
          ENDDO
          nc=k; nv=j
          IF (nc>0) THEN
             OPEN (unit=1001,file='cwfn'//sub//TRIM(stuff),form='formatted')
             DO i = 1,n
                WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&                    (Orbitwk%wfn(i,cmap(k)),k=1,nc)
             ENDDO
             CLOSE(1001)
          ENDIF
          IF (nv>0) THEN
             OPEN (unit=1001,file='vwfn'//sub//TRIM(stuff),form='formatted')
             DO i = 1,n
                WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&                    (Orbitwk%wfn(i,vmap(k)),k=1,nv)
             ENDDO
             CLOSE(1001)
          ENDIF
       ELSE
          OPEN (unit=1001,file='wfn'//sub//TRIM(stuff),form='formatted')
          DO i = 1,n
             WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&                 (Orbitwk%wfn(i,j),j=1,Orbitwk%norbit)
          ENDDO
          CLOSE(1001)
       ENDIF
       OPEN (unit=1001,file='pot'//sub//TRIM(stuff),form='formatted')
       DO i = 1,n
          WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),Potwk%rv(i),Potwk%rvh(i)
       ENDDO
       CLOSE(1001)
       OPEN (unit=1001,file='den'//sub//TRIM(stuff),form='formatted')
       IF (frozencorecalculation) then
          DO i = 1,n
             WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), Orbitwk%den(i), &
&               FCwk%valeden(i), FCwk%coreden(i)
          ENDDO
       ELSE
          DO i = 1,n
             WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), Orbitwk%den(i)
          ENDDO
       ENDIF
       CLOSE(1001)

       IF (frozencorecalculation.AND.Orbitwk%exctype/='HFV') THEN
          IF (nc>0) THEN
             OPEN (unit=1001,file='cfx'//sub//TRIM(stuff),form='formatted')
             DO i = 1,n
                WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&                    (HF%sumY(i,cmap(k)),k=1,nc)
             ENDDO
             CLOSE(1001)
          ENDIF
          IF (nv>0) THEN
             OPEN (unit=1001,file='vfx'//sub//TRIM(stuff),form='formatted')
             DO i = 1,n
                WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),&
&                    (HF%sumY(i,vmap(k)),k=1,nv)
             ENDDO
             CLOSE(1001)
          ENDIF
       ELSEIF (frozencorecalculation.AND.Orbitwk%exctype=='HFV') THEN
          IF (nv>0) THEN
             OPEN (unit=1001,file='cvfx'//sub//TRIM(stuff),form='formatted')
             OPEN (unit=1002,file='vvfx'//sub//TRIM(stuff),form='formatted')
             DO i = 1,n
                WRITE(1001,'(1p,50e15.7)') Gridwk%r(i),HF%rvxcore(i),&
&                    (HF%sumYCV(i,vmap(k)),k=1,nv)
                WRITE(1002,'(1p,50e15.7)') Gridwk%r(i),&
&                    (HF%sumYVV(i,vmap(k)),k=1,nv)
             ENDDO
             CLOSE(1001)
             CLOSE(1002)
          ENDIF
       ELSE
          OPEN (unit=1001,file='fx'//sub//TRIM(stuff),form='formatted')
          DO i = 1,n
             WRITE(1001,'(1p,50e15.7)') Gridwk%r(i), &
&                 (HF%sumY(i,j),j=1,Orbitwk%norbit)
          ENDDO
          CLOSE(1001)
       ENDIF

       IF (ALLOCATED(vmap)) DEALLOCATE(vmap)
       IF (ALLOCATED(cmap)) DEALLOCATE(cmap)

       counter=counter+1

     END SUBROUTINE Report_HF_functions

     SUBROUTINE HFdump(Grid,Orbit,ifo)
       TYPE(GridInfo), INTENT(IN) :: Grid
       TYPE(OrbitInfo), INTENT(IN) :: Orbit
       INTEGER, INTENT(IN) :: ifo

       INTEGER :: i,io,jo,n,norbit

       n=Grid%n;norbit=Orbit%norbit

       WRITE(ifo,*) HF%lmax
       WRITE(ifo,*) (HF%emin(io),HF%emax(io),io=1,norbit)
       WRITE(ifo,*) ((HF%lmbd(io,jo),io=1,norbit),jo=1,norbit)
       WRITE(ifo,*) ((HF%SumY(i,io),i=1,n),io=1,norbit)
       WRITE(ifo,*) (HF%lmany(io),io=1,norbit),(HF%CSlmany(io),io=1,norbit)
       WRITE(ifo,*) ((HF%lmap(io,jo),io=1,norbit),jo=1,norbit)
       WRITE(ifo,*) ((HF%CSlmap(io,jo),io=1,norbit),jo=1,norbit)

     END SUBROUTINE HFdump

     SUBROUTINE HFload(Grid,Orbit,ifo)
       TYPE(GridInfo), INTENT(IN) :: Grid
       TYPE(OrbitInfo), INTENT(IN) :: Orbit
       INTEGER, INTENT(IN) :: ifo

       INTEGER :: i,io,jo,n,norbit

       n=Grid%n;norbit=Orbit%norbit

       READ(ifo,*) HF%lmax

       ALLOCATE(HF%emin(norbit),HF%emax(norbit),HF%lmbd(norbit,norbit),&
&           HF%SumY(n,norbit),HF%lmany(norbit),HF%CSlmany(norbit), &
&           HF%lmap(norbit,norbit),HF%CSlmap(norbit,norbit))

       READ(ifo,*) (HF%emin(io),HF%emax(io),io=1,norbit)
       READ(ifo,*) ((HF%lmbd(io,jo),io=1,norbit),jo=1,norbit)
       READ(ifo,*) ((HF%SumY(i,io),i=1,n),io=1,norbit)
       READ(ifo,*) (HF%lmany(io),io=1,norbit),(HF%CSlmany(io),io=1,norbit)
       READ(ifo,*) ((HF%lmap(io,jo),io=1,norbit),jo=1,norbit)
       READ(ifo,*) ((HF%CSlmap(io,jo),io=1,norbit),jo=1,norbit)

     END SUBROUTINE HFload


     Subroutine DestroyHF(HF)
       Type(HFInfo), INTENT(INOUT) :: HF

       If (ASSOCIATED(HF%lmbd)) DEALLOCATE(HF%lmbd)
       If (ASSOCIATED(HF%SumY)) DEALLOCATE(HF%SumY)
       If (ASSOCIATED(HF%emin)) DEALLOCATE(HF%emin)
       If (ASSOCIATED(HF%emax)) DEALLOCATE(HF%emax)
       If (ASSOCIATED(HF%SumYCV)) DEALLOCATE(HF%SumYCV)
       If (ASSOCIATED(HF%SumYVV)) DEALLOCATE(HF%SumYVV)
       If (ASSOCIATED(HF%rVxcore)) DEALLOCATE(HF%rVxcore)
       If (ASSOCIATED(HF%lmany)) DEALLOCATE(HF%lmany)
       If (ASSOCIATED(HF%lmap)) DEALLOCATE(HF%lmap)
       If (ASSOCIATED(HF%CSlmany)) DEALLOCATE(HF%CSlmany)
       If (ASSOCIATED(HF%CSlmap)) DEALLOCATE(HF%CSlmap)

     End Subroutine DestroyHF
   END MODULE hf_mod
