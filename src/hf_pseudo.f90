MODULE hf_pseudo
  USE atomdata
  USE hf_mod
  USE fock
  USE globalmath
  USE gridmod
  USE paw_sub
  USE pseudodata
  USE pseudo_sub

  IMPLICIT NONE

CONTAINS


  !modified version of SUBROUTINE makebasis_custom which was
  !           written by Marc Torrent from earlier version by NAWH
  !           only one option programmed here
  ! This version does not calculate tp because local pseudopotential is not
  !     yet known
  SUBROUTINE make_hf_basis_only(Grid,Pot,PAW,ifinput)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    INTEGER, INTENT(IN) :: ifinput

    INTEGER ::  i,j,k,l,n,ib,jb,io,jo,norbit,nbase,irc,np,thisrc,lmax,ni
    INTEGER :: icount,ioo,joo,lng
    INTEGER, ALLOCATABLE :: omap(:)
    REAL(8) :: rc,gp,gpp,xx,occ,stuff
    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: Ci(:),dum(:),tdum(:),hat(:),arg(:),aa(:,:),ai(:,:)
    CHARACTER(132) :: inputline
    REAL(8), PARAMETER :: threshold=1.d-8
    LOGICAL :: last

    write(6,*) 'Entering make_hf_basis '; call flush(6)
    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    lmax=PAW%lmax
    nbase=PAW%nbase
    np=5
    ALLOCATE(Ci(np),arg(n),dum(n),tdum(n),hat(n))

    WRITE(6,*) 'For each of the following basis functions enter rc'

    ! Loop on basis elements
    DO ib=1,nbase
       l=PAW%l(ib)
       lng=min(n,PAW%rng(ib))
       ! Read matching radius
       WRITE(6,'(3i5,1p,1e15.7)') ib,PAW%np(ib),PAW%l(ib),PAW%eig(ib)
       READ(5,'(a)') inputline
       WRITE(ifinput,'(a)') TRIM(inputline)
       READ(inputline,*) rc
       thisrc=FindGridIndex(Grid,rc)
       thisrc=MIN(thisrc,irc)       ! make sure rc<total rc
       rc=r(thisrc)
       WRITE(6,*) 'rc for this wfn', rc
       IF (thisrc<3.OR.thisrc>irc.OR.thisrc>lng-3) THEN
          WRITE(6,*) 'rc out of range', thisrc,lng,irc
          STOP
       ENDIF
       CALL pspolyn(PAW%phi(:,ib),Ci,r,l,np,thisrc,lng)
       ! Compute pseudized partial wave and part of  projector
       PAW%tphi(:,ib)=0;PAW%tphi(1:lng,ib)=PAW%phi(1:lng,ib)
       PAW%tp(:,ib)=0.d0 ; PAW%Kop(:,ib)=0
       DO i=1,thisrc-1
          xx=r(i)*r(i)
          PAW%tphi(i,ib)=Ci(1)+Ci(2)*xx
          gp=2.d0*Ci(2)
          gpp=2.d0*Ci(2)
          DO j=3,np
             PAW%tphi(i,ib)=PAW%tphi(i,ib)+Ci(j)*xx**(j-1)
             gp=gp+DBLE(2*j-2)*Ci( j)*xx**(j-2)
             gpp=gpp+DBLE((2*j-2)*(2*j-3))*Ci(j)*xx**(j-2)
          ENDDO
          PAW%tphi(i,ib)=PAW%tphi(i,ib)*r(i)**(l+1)
          PAW%Kop(i,ib)=-(DBLE(2*(l+1))*gp+gpp)*r(i)**(l+1)  !kinetic pt.
       ENDDO
       DO i=thisrc,lng
          gpp=Gsecondderiv(Grid,i,PAW%tphi(:,ib))
          PAW%Kop(i,ib)=-gpp+DBLE(l*(l+1))/(r(i)**2)*PAW%tphi(i,ib)
       ENDDO
    ! test
      !do i=1,n
      !   write(200+ib,'(1p,20e15.7)') r(i),PAW%Kop(i,ib),PAW%tp(i,ib),&
      !&      -Gsecondderiv(Grid,i,PAW%phi(:,ib))+&
      !&       DBLE(l*(l+1))/(r(i)**2)*PAW%phi(i,ib)
      !enddo
    !end test
    ENDDO     !nbase

    ! loop on all occupied states  -- assume all core states localized
    PAW%tcore=0
    norbit=PAW%OCCWFN%norbit
    thisrc=irc; rc=r(irc)
    DO io=1,norbit
       IF(PAW%valencemap(io)>0) THEN
          ib=PAW%valencemap(io)
          PAW%TOCCWFN%wfn(:,io)=PAW%tphi(:,ib)
          !PAW%OCCWFN%wfn(:,io)=PAW%phi(:,ib)    !already set in setbasis
          write(6,*) 'io = ', io,ib
       ELSE    ! presumably a core state
          PAW%TOCCWFN%wfn(:,io)=0
          write(6,*) 'io =  core ', io
          !!! comment out code for extended core states
          !If (PAW%TOCCWFN%iscore(io)) then
          !   last=.true. ; l=PAW%TOCCWFN%l(io)
          ! ! check to find if this is the outermost core state for this l
          !   if (io<norbit) then
          !      do jo=io+1,norbit
          !         if (PAW%TOCCWFN%iscore(jo).and.PAW%TOCCWFN%l(jo)==l) &
          !&              last=.false.
          !      enddo
          !   endif
          !   if (last) then
          !      IF (MAXVAL(ABS(PAW%OCCWFN%wfn(irc-np/2:irc+np/2,io))) &
          !&            >PAW%coretol) THEN
          !         CALL Smoothfunc(r,l,np,thisrc,n,PAW%OCCWFN%wfn(:,io), &
          !&                PAW%TOCCWFN%wfn(:,io))
          !      ENDIF
          !   endif
          !else
          !    write(6,*) 'Something wrong -- should be core state ',&
          !&      io,PAW%TOCCWFN%l(io),PAW%TOCCWFN%eig(io),PAW%TOCCWFN%occ(io)
          !    stop
          !endif
       ENDIF
    ENDDO
    DEALLOCATE(Ci,arg,dum,tdum,hat)
  END SUBROUTINE make_hf_basis_only

  !Calculate tp for hf case;  assume tphi, Kop and vreff known
  SUBROUTINE make_hf_tp_only(Grid,Pot,PAW,ifinput)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    INTEGER, INTENT(IN) :: ifinput

    INTEGER ::  i,j,k,l,n,ib,jb,io,jo,norbit,nbase,irc,np,thisrc,lmax,ni
    INTEGER :: icount,ioo,joo,lng
    INTEGER, ALLOCATABLE :: omap(:)
    REAL(8) :: rc,gp,gpp,xx,occ,stuff
    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: Ci(:),dum(:),tdum(:),hat(:),arg(:),aa(:,:),ai(:,:)
    CHARACTER(132) :: inputline
    REAL(8), PARAMETER :: threshold=1.d-8
    LOGICAL :: last

    write(6,*) 'Entering make_hf_tp '; call flush(6)
    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    lmax=PAW%lmax
    nbase=PAW%nbase
    np=5
    ALLOCATE(Ci(np),arg(n),dum(n),tdum(n),hat(n))

    norbit=PAW%OCCWFN%norbit
    ! Find smooth Coulombic potential for V_N + V_H
    dum=0; tdum=0; PAW%core=0;   PAW%tcore=0 ; PAW%den=0; PAW%tden=0
    CALL hatL(Grid,PAW,0,hat)
    DO io=1,norbit
       occ=PAW%OCCwfn%occ(io)
       if (occ>1.d-5) then
          dum=dum+occ*PAW%OCCwfn%wfn(:,io)**2
          tdum=tdum+occ*PAW%TOCCwfn%wfn(:,io)**2
       endif
       IF (PAW%OCCwfn%iscore(io)) THEN
          PAW%core=PAW%core+occ*PAW%OCCwfn%wfn(:,io)**2
          !PAW%tcore=PAW%tcore+occ*PAW%TOCCwfn%wfn(:,io)**2
       ELSE
          PAW%den=PAW%den+occ*PAW%OCCwfn%wfn(:,io)**2
          PAW%tden=PAW%tden+occ*PAW%TOCCwfn%wfn(:,io)**2
       ENDIF
    ENDDO

    ! tdum loads PAW%rveff/r
    tdum=0; tdum(2:n)=PAW%rveff(2:n)/r(2:n)
    call extrapolate(Grid,tdum)

    ! complete projector functions  ;
    DO ib=1,nbase
       PAW%tp(:,ib)=PAW%Kop(:,ib)
       l=PAW%l(ib)
       lng=PAW%rng(ib)
       CALL Calc_tXp_basis(Grid,PAW,ib,dum,lng);
       ! testing
       !arg=0
       !ni=PAW%np(ib)
       !DO jo=1,norbit
       !   IF(PAW%valencemap(jo)==ib) THEN
       !      ni=jo
       !   ENDIF
       !ENDDO
       !CALL Calc_XP(Grid,PAW%OCCwfn,ni,PAW%l(ib),PAW%phi(:,ib),hat)

       !DO i=1,n
       !   WRITE(300+ib,'(1p,20e15.7)') r(i),dum(i),hat(i)
       !   arg(i)=   -Gsecondderiv(Grid,i,PAW%phi(:,ib))+&
       !&                DBLE(l*(l+1))/(r(i)**2)*PAW%phi(i,ib)
       !ENDDO
       ! end testing

       dum(1)=0; dum(2:lng)=dum(2:lng)/r(2:lng)
       PAW%tp(1:lng,ib)=PAW%tp(1:lng,ib)+tdum(1:lng)*PAW%tphi(1:lng,ib)+&
&                dum(1:lng)
       arg(2:n)=arg(2:n)+Pot%rv(2:n)*PAW%phi(2:n,ib)/r(2:n)+hat(2:n)/r(2:n)
       if (PAW%eig(ib)>0.d0) then
          if (PAW%occ(ib)> 0.001d0) then
            write(6,*) 'Problem with basis ib', &
&           ib,PAW%l(ib),PAW%eig(ib),PAW%occ(ib)
            stop
          endif
          PAW%tp(1:lng,ib)=PAW%tp(1:lng,ib)-PAW%eig(ib)*PAW%tphi(1:lng,ib)
          arg(:)=arg(:)-PAW%eig(ib)*PAW%phi(:,ib)
       else
          ioo=0;joo=0
          DO io=1,norbit
             IF(PAW%valencemap(io)==ib) ioo=io
          ENDDO
          Do jo=1,norbit
             IF(PAW%OCCWFN%l(jo)==l.and..not.PAW%OCCWFN%iscore(jo)) then
              joo=jo
              PAW%tp(1:lng,ib)=PAW%tp(1:lng,ib)-&
&                    PAW%OCCWFN%lqp(joo,ioo)*PAW%TOCCWFN%wfn(1:lng,joo)
              arg(:)=arg(:)-&
&                    PAW%OCCWFN%lqp(joo,ioo)*PAW%OCCWFN%wfn(:,joo)
              write(6,*) 'ioo,joo,lam', ioo,joo,PAW%OCCWFN%lqp(joo,ioo)
             ENDIF
          Enddo
       endif
       !do i=1,n
       !    write(900+ib,'(1p,20e15.7)') r(i),PAW%tp(i,ib),arg(i)
       !enddo
    ENDDO

    open(123,file="tp",form="formatted")
    do i=1,n
       write(123,'(1p,50e15.7)') Grid%r(i),(PAW%tp(i,ib),ib=1,nbase)
    enddo
    close(123)
    ! Form orthogonalized projector functions
    DO ib=1,nbase
       PAW%ophi(:,ib)=PAW%phi(:,ib)
       PAW%otphi(:,ib)=PAW%tphi(:,ib)
    ENDDO

    DO l=0,lmax
       icount=0
       DO ib=1,nbase
          IF (PAW%l(ib)==l) icount=icount+1
       ENDDO
       WRITE(6,*) 'For l = ', l,icount,' basis functions'
       ALLOCATE(aa(icount,icount),ai(icount,icount),omap(icount))
       aa=0;icount=0;lng=n
       DO ib=1,nbase
          IF (PAW%l(ib)==l) THEN
             icount=icount+1
             omap(icount)=ib
             lng=min(lng,PAW%rng(ib))
          ENDIF
       ENDDO
       DO i=1,icount
          ib=omap(i)
          DO j=1,icount
             jb=omap(j)
             aa(i,j)=overlap(Grid,PAW%otphi(:,ib),PAW%tp(:,jb),1,lng)
             write(6,*) 'Overlap', ib,jb,aa(i,j)
          ENDDO
       ENDDO
       ai=aa;CALL minverse(ai,icount,icount,icount)

       DO i=1,icount
          ib=omap(i)
          PAW%ck(ib)=ai(i,i)
          PAW%otp(:,ib)=0
          DO j=1,icount
             jb=omap(j)
             PAW%otp(1:lng,ib)=PAW%otp(1:lng,ib)+PAW%tp(1:lng,jb)*ai(j,i)
          ENDDO
       ENDDO
       DEALLOCATE(aa,ai,omap)
    ENDDO

    DEALLOCATE(Ci,arg,dum,tdum,hat)
  END SUBROUTINE make_hf_tp_only


  SUBROUTINE Smoothfunc(r,l,np,thisrc,n,wfnin,wfnout)
  REAL(8), INTENT(IN):: r(:),wfnin(:)
  INTEGER, INTENT(IN) :: l,np,thisrc,n
  REAL(8), INTENT(OUT) :: wfnout(:)

  REAL(8), ALLOCATABLE :: Ci(:)
  REAL(8) :: xx
  INTEGER :: i,j

  ALLOCATE(Ci(np))
  wfnout(:)=0
  CALL pspolyn(wfnin(:),Ci,r,l,np,thisrc,n)
  wfnout(:)=wfnin(:)
  DO i=1,thisrc-1
     xx=r(i)*r(i)
     wfnout(i)=Ci(1)+Ci(2)*xx
     DO j=3,np
        wfnout(i)=wfnout(i)+Ci(j)*xx**(j-1)
     ENDDO
     wfnout(i)=wfnout(i)*r(i)**(l+1)
  ENDDO

  DEALLOCATE(Ci)
END SUBROUTINE Smoothfunc

  ! Version of Calc_tXp for PAW basis functions only
  SUBROUTINE Calc_tXp_basis(Grid,PAW,ib,tres,lng)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(PseudoInfo), INTENT(IN):: PAW
    INTEGER, INTENT(IN) :: ib,lng
    REAL(8), INTENT(OUT) :: tres(:)

    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:),arg(:),f(:),hat(:)
    INTEGER :: i,j,k,n,m,ll,li,lj,lmin,lmax,io,jo,ok,norbit,ko
    INTEGER :: nu,nup,nj,ni
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo
    REAL(8), PARAMETER :: threshold=1.d-8

    tres=0.d0
    norbit=PAW%OCCwfn%norbit
    n=Grid%n
    !lng=min(n,PAW%rng(ib))

    ALLOCATE(wfp(n),vl(n),dum(n),arg(n),f(n),hat(n))
    r=>Grid%r
    li=PAW%l(ib)
    ni=PAW%np(ib)
    DO jo=1,norbit
       IF(PAW%valencemap(jo)==ib) THEN
          ni=jo
       ENDIF
    ENDDO

    occ=1.d0

    DO jo=1,norbit
       occj=PAW%OCCwfn%occ(jo); vl=0.d0
       IF (occj>threshold) THEN
          lj=PAW%OCCwfn%l(jo);nj=PAW%OCCwfn%np(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:lng)=PAW%tphi(1:lng,ib)*PAW%TOCCwfn%wfn(1:lng,jo)
          arg(1:lng)=PAW%phi(1:lng,ib)*PAW%OCCwfn%wfn(1:lng,jo)-wfp(1:lng)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,ni,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL hatL(Grid,PAW,ll,dum)
                f=arg*(r**ll)
                f=wfp+integrator(Grid,f,1,lng)*dum
                CALL apoisson(Grid,ll,lng,f,dum)
                vl(1:lng)=vl(1:lng)-wgt*dum(1:lng)*PAW%TOCCwfn%wfn(1:lng,jo)/occ
             ENDIF
          ENDDO
       ENDIF
       tres=tres+vl
    ENDDO   !jo

    DEALLOCATE(wfp,vl,dum,arg,f,hat)
  END SUBROUTINE Calc_tXp_basis

  SUBROUTINE HFOrthotocore(Grid,l,wfn,PAW)
    TYPE(GridInfo) , INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(INOUT) :: wfn(:)
    TYPE(PseudoInfo), INTENT(IN) :: PAW

    INTEGER :: io,ib
    REAL(8) :: x,y
    REAL(8) :: small=1.d-9

    DO io=1,PAW%TOCCWFN%norbit
       IF (PAW%TOCCWFN%iscore(io).AND.PAW%TOCCWFN%l(io)==l) THEN
          x=CoreOverlap(Grid,PAW,l,io,wfn)
          y=CoreOverlap(Grid,PAW,l,io,PAW%TOCCWFN%wfn(:,io))
          if (ABS(x)>small) wfn=wfn-PAW%TOCCWFN%wfn(:,io)*x/y
       ENDIF
    ENDDO
  END SUBROUTINE HFOrthotocore

  FUNCTION CoreOverlap(Grid,PAW,l,ic,wfn)
    REAL(8) :: CoreOverlap
    TYPE(GridInfo) , INTENT(IN) :: Grid
    INTEGER, INTENT(IN) :: l,ic
    REAL(8), INTENT(IN) :: wfn(:)
    TYPE(PseudoInfo), INTENT(IN) :: PAW

    INTEGER :: io,ib
    REAL(8) :: x
    REAL(8) :: small=1.d-9

    CoreOverlap=0.d0 ; x=0.d0
    write(6,*) 'In CoreOverlap ',l,ic
    if(PAW%TOCCWFN%iscore(ic).and.PAW%TOCCWFN%l(ic)==l) then
          x=overlap(Grid,wfn,PAW%TOCCWFN%wfn(:,ic))
          write(6,*) 'In CoreOverlap ', x
          DO ib=1,PAW%nbase
             IF (PAW%l(ib)==l.AND.ABS(PAW%mLic(ib,ic,1))>small) THEN
                x=x+overlap(Grid,wfn,PAW%otp(:,ib))*PAW%mLic(ib,ic,1)
                write(6,*) 'In CoreOverlap ', ib,x
             ENDIF
          ENDDO
          WRITE(6,*) 'overlap with core ', l,ic,x
     ENDIF
     CoreOverlap=x

  END FUNCTION CoreOverlap

  ! HF version of generalized Graham-Schmidt orthogonalization of wfn1 to wfn2
  ! on output <wfn1|O|wfn2>=0
  SUBROUTINE HFgenOrthog(Grid,PAW,l,wfn1,wfn2)
    !REAL(8) :: genoverlap
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(INOUT) :: wfn1(:)
    REAL(8), INTENT(IN) :: wfn2(:)

    REAL(8) :: x,y

    !   Orthogonalize
    x=HFgenoverlap(Grid,PAW,l,wfn1,wfn2)
    y=HFgenoverlap(Grid,PAW,l,wfn2,wfn2)

    WRITE(6,*) 'overlap ', x,y
    wfn1(:)=wfn1(:)-(x/y)*wfn2(:)

  END SUBROUTINE HFgenOrthog

  !*******************************************************************
  !  HF version of genoverlap
  !*******************************************************************
  FUNCTION HFgenoverlap(Grid,PAW,l,wfn1,wfn2)
    REAL(8) :: HFgenoverlap
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: l
    REAL(8), INTENT(IN) :: wfn1(:),wfn2(:)

    INTEGER :: n,ib,ic,nbase
    REAL(8) :: h
    REAL(8), ALLOCATABLE :: a(:),b(:)

    n=Grid%n; h=Grid%h; nbase=PAW%nbase;
    ALLOCATE(a(nbase),b(nbase),stat=ib)
    IF (ib/=0) THEN
       WRITE(6,*) 'Error in HFgenoverlap allocation', nbase,ib
       STOP
    ENDIF

    HFgenoverlap=overlap(Grid,wfn1,wfn2)
    a=0; b=0
    DO ib=1,nbase
       IF (l==PAW%l(ib)) THEN
          a(ib)=overlap(Grid,wfn1,PAW%otp(:,ib))
          b(ib)=overlap(Grid,wfn2,PAW%otp(:,ib))
       ENDIF
    ENDDO
    DO ib=1,nbase
       DO ic=1,nbase
          HFgenoverlap=HFgenoverlap+a(ib)*PAW%oij(ib,ic)*b(ic)
       ENDDO
    ENDDO

    DEALLOCATE(a,b)
  END FUNCTION HFgenoverlap

  SUBROUTINE GetMoment(Grid,PAW,io,jo,ll,o,ml,zero)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(IN) :: PAW
    INTEGER, INTENT(IN) :: io,jo,ll
    REAL(8), INTENT(IN) :: o(:,:)   ! o(io,ib) for valence states
    REAL(8), INTENT(OUT) :: ml(:)
    LOGICAL, INTENT(OUT) :: zero

    INTEGER :: ib, jb,l,li,lj
    REAL(8) :: x
    REAL(8) :: thres=1.d-9
    REAL(8) , ALLOCATABLE :: arg(:)

    ml=0; zero=.TRUE.

    IF (PAW%OCCWFN%iscore(io)) THEN
       WRITE(6,*) 'In GetMoment -- io must not be core state !'
       STOP
    ENDIF

    li=PAW%OCCWFN%l(io)
    lj=PAW%OCCWFN%l(jo)
    !write(6,*) 'In GetMoment ', li,lj,ll
    IF (ll<ABS(li-lj).OR.ll>(li+lj)) THEN
       RETURN
    ENDIF

    zero=.FALSE.
    IF (PAW%OCCWFN%iscore(jo)) THEN
       DO ib=1,PAW%nbase
          IF (PAW%l(ib)==li) THEN
             ml=ml+o(io,ib)*PAW%mLic(ib,jo,ll+1)*PAW%g(:,ll+1)
          ENDIF
       ENDDO
    ELSE
       DO ib=1,PAW%nbase
          IF (PAW%l(ib)==li) THEN
             DO jb=1,PAW%nbase
                IF (PAW%l(jb)==lj) THEN
                   ml=ml+&
&                       o(io,ib)*o(jo,jb)*PAW%mLij(ib,jb,ll+1)*PAW%g(:,ll+1)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDIF

    !allocate(arg(Grid%n))
    !arg=ml*(Grid%r**ll)
    !write(6,*) 'GetMoment',io,jo,ll,integrator(Grid,arg)
    !deallocate(arg)
  END SUBROUTINE GetMoment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Calc_tXv(Grid,PAW,io,o,tres,Xiv,tEx,Eax)
  !      Calculate \tilde{X}_v for valence orbital io as well as
  !        one-center matrix elements Xiv(ib,io), smooth energy tEx
  !           and one-center contributions Eax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Calc_tXv(Grid,PAW,io,o,tres,Xiv,tEx,Eax)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(PseudoInfo), INTENT(IN):: PAW
    INTEGER, INTENT(IN) :: io
    REAL(8), INTENT(IN) :: o(:,:)      ! o(io,ib) for valence states
    REAL(8), INTENT(OUT) :: tres(:) ,Xiv(:) ,tEx, Eax

    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:),arg(:),f(:),hat(:)
    INTEGER :: i,j,k,n,m,ll,li,lj,lmin,lmax,jo,ok,norbit,ko,ib,jb,kb,lb
    INTEGER :: nu,nup,nj,ni
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo,fac,x
    REAL(8), PARAMETER :: threshold=1.d-8
    LOGICAL :: iszero


    tres=0.d0; Xiv=0.d0; tEx=0.d0; Eax=0.d0
    norbit=PAW%OCCwfn%norbit
    n=Grid%n

    ALLOCATE(wfp(n),vl(n),dum(n),arg(n),f(n),hat(n))
    r=>Grid%r
    li=PAW%OCCwfn%l(io)

    occ=1.d0

    DO jo=1,norbit
       occj=PAW%OCCwfn%occ(jo); vl=0.d0
       IF (occj>threshold) THEN
          lj=PAW%OCCwfn%l(jo)
          IF (PAW%OCCwfn%iscore(jo)) THEN
             fac=1.d0
          ELSE
             fac=0.5d0
          ENDIF
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=PAW%TOCCwfn%wfn(1:n,io)*PAW%TOCCwfn%wfn(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL GetMoment(Grid,PAW,io,jo,ll,o,arg,iszero)
                IF (iszero) THEN
                   WRITE(6,*) 'Zero moment ', io,jo,li,lj,ll
                   STOP
                ENDIF
                arg=arg+wfp
                CALL apoisson(Grid,ll,n,arg,dum)
                vl(1:n)=vl(1:n)-wgt*dum(1:n)*PAW%TOCCwfn%wfn(1:n,jo)/occ
                dum(1)=0; dum(2:n)=dum(2:n)/Grid%r(2:n)
                tEx=tEx-fac*wgt*overlap(Grid,arg,dum)/occ
                IF (PAW%OCCwfn%iscore(jo)) THEN
                   DO ib=1,PAW%nbase      ! add Z terms to Xiv
                      x=PAW%mLic(ib,jo,ll+1)
                      IF (ABS(x)>threshold) THEN
                         f=x*PAW%g(:,ll+1)
                         Xiv(ib)=Xiv(ib)-wgt*overlap(Grid,dum,f)/occ
                      ENDIF
                   ENDDO
                ELSE
                   DO ib=1,PAW%nbase      ! add Z terms to Xiv
                      x=0
                      DO jb=1,PAW%nbase
                         IF (PAW%l(jb)==lj) THEN
                            x=x+o(jo,jb)*PAW%mLij(ib,jb,ll+1)
                         ENDIF
                      ENDDO
                      f=x*PAW%g(:,ll+1)
                      Xiv(ib)=Xiv(ib)-wgt*overlap(Grid,dum,f)/occ
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF
       tres=tres+vl
    ENDDO   !jo

    !  complete Xiv terms
    DO jo=1,norbit
       occj=PAW%OCCwfn%occ(jo); vl=0.d0
       IF (occj>threshold) THEN
          lj=PAW%OCCwfn%l(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          IF (PAW%OCCwfn%iscore(jo)) THEN
             DO ll=lmin,lmax,2
                CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                IF (wgt>threshold) THEN
                   wgt=wgt*(2*ll+1)   ! because of apoisson convention
                   DO ib=1,PAW%nbase
                      x=0
                      DO jb=1,PAW%nbase
                         IF (PAW%l(jb)==li) THEN
                            x=x+o(io,jb)*PAW%DRC(ib,jb,jo,ll+1)
                         ENDIF
                      ENDDO
                      Xiv(ib)=Xiv(ib)-wgt*x/occ
                      Eax=Eax-wgt*x*o(io,ib)/occ
                   ENDDO
                ENDIF
             ENDDO
          ELSE
             DO ll=lmin,lmax,2
                CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                IF (wgt>threshold) THEN
                   wgt=wgt*(2*ll+1)   ! because of apoisson convention
                   DO ib=1,PAW%nbase
                      x=0
                      DO jb=1,PAW%nbase
                         IF (PAW%l(jb)==lj) THEN
                            DO kb=1,PAW%nbase
                               IF (PAW%l(kb)==lj) THEN
                                  DO lb=1,PAW%nbase
                                     IF (PAW%l(lb)==li) THEN
                                        x=x+o(jo,jb)*o(jo,kb)*o(io,lb)*&
&                                            PAW%DR(ib,jb,kb,lb,ll+1)
                                     ENDIF
                                  ENDDO
                               ENDIF
                            ENDDO
                         ENDIF
                      ENDDO
                      Xiv(ib)=Xiv(ib)-wgt*x/occ
                      Eax=Eax-0.5d0*wgt*x*o(io,ib)/occ
                   ENDDO
                ENDIF
             ENDDO
          ENDIF
       ENDIF
    ENDDO

    DEALLOCATE(wfp,vl,dum,arg,f,hat)
  END SUBROUTINE Calc_tXv

  SUBROUTINE Calc_Xc(Grid,PAW,io,ic,o,Xcv)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(PseudoInfo), INTENT(IN):: PAW
    INTEGER, INTENT(IN) :: io ,ic   !valence, core
    REAL(8), INTENT(IN) :: o(:,:)      ! o(io,ib) for valence states
    REAL(8), INTENT(OUT) :: Xcv

    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:),arg(:),f(:),hat(:)
    INTEGER :: i,j,k,n,m,ll,li,lj,lmin,lmax,jo,ok,norbit,ko,ib,jb,kb,lb
    INTEGER :: nu,nup,nj,ni
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo,fac,x
    REAL(8), PARAMETER :: threshold=1.d-8
    LOGICAL :: iszero

    Xcv=0.d0
    IF (PAW%OCCWFN%iscore(io).OR..NOT.PAW%OCCWFN%iscore(ic)) THEN
       WRITE(6,*) 'Calc_Xc assumes io is valence and ic core ', ic
       STOP
    ENDIF

    norbit=PAW%OCCwfn%norbit
    n=Grid%n

    ALLOCATE(wfp(n),vl(n),dum(n),arg(n),f(n),hat(n))
    r=>Grid%r
    li=PAW%OCCwfn%l(io)

    occ=1.d0

    DO jo=1,PAW%TOCCwfn%norbit
       occj=PAW%OCCwfn%occ(jo)
       IF (occj>threshold) THEN
          lj=PAW%OCCwfn%l(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:n)=PAW%TOCCwfn%wfn(1:n,io)*PAW%TOCCwfn%wfn(1:n,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                CALL GetMoment(Grid,PAW,io,jo,ll,o,arg,iszero)
                IF (iszero) THEN
                   WRITE(6,*) 'zero moment', io,jo,li,lj,ll
                   STOP
                ENDIF
                arg=arg+wfp
                CALL apoisson(Grid,ll,n,arg,dum)
                dum(1)=0; dum(2:n)=dum(2:n)/Grid%r(2:n)
                IF (PAW%OCCwfn%iscore(jo)) THEN
                   x=PAW%mLcc(jo,ic,ll+1)  ! add Z terms to Xcv
                   IF (ABS(x)>threshold) THEN
                      f=x*PAW%g(:,ll+1)
                      Xcv=Xcv-wgt*overlap(Grid,dum,f)/occ
                   ENDIF
                ELSE
                   x=0              ! add Z terms to Xcv
                   DO jb=1,PAW%nbase
                      IF (PAW%l(jb)==lj) THEN
                         x=x+o(jo,jb)*PAW%mLic(jb,ic,ll+1)
                      ENDIF
                   ENDDO
                   f=x*PAW%g(:,ll+1)
                   Xcv=Xcv-wgt*overlap(Grid,dum,f)/occ
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO   !jo

    !  complete Xcv terms
    DO jo=1,norbit
       occj=PAW%OCCwfn%occ(jo)
       IF (occj>threshold) THEN
          lj=PAW%OCCwfn%l(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          IF (PAW%OCCwfn%iscore(jo)) THEN
             DO ll=lmin,lmax,2
                CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                IF (wgt>threshold) THEN
                   wgt=wgt*(2*ll+1)   ! because of apoisson convention
                   x=0
                   DO jb=1,PAW%nbase
                      IF (PAW%l(jb)==li) THEN
                         x=x+o(io,jb)*PAW%DRCC(jb,ic,jo,ll+1) !jo last index
                      ENDIF
                   ENDDO
                   Xcv=Xcv-wgt*x/occ
                ENDIF
             ENDDO
          ELSE
             DO ll=lmin,lmax,2
                CALL EXXwgt(occ,occj,io,li,jo,lj,ll,wgt)
                IF (wgt>threshold) THEN
                   wgt=wgt*(2*ll+1)   ! because of apoisson convention
                   x=0
                   DO jb=1,PAW%nbase
                      IF (PAW%l(jb)==lj) THEN
                         DO kb=1,PAW%nbase
                            IF (PAW%l(kb)==lj) THEN
                               DO lb=1,PAW%nbase
                                  IF (PAW%l(lb)==li) THEN
                                     x=x+o(jo,jb)*o(jo,kb)*o(io,lb)*&
&                                         PAW%DRCjkl(ic,jb,kb,lb,ll+1)
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                   Xcv=Xcv-wgt*x/occ
                ENDIF
             ENDDO
          ENDIF
       ENDIF
    ENDDO

    DEALLOCATE(wfp,vl,dum,arg,f,hat)
  END SUBROUTINE Calc_Xc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  PAWIter_HF(Grid,PAW,err0,err,success)
  !     On input PAW%TOCCWFN  contains initial guess of smooth wfn's
  !        On output the guess is updated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PAWIter_HF(Grid,PAW,mix,err0,err,success)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    REAL(8) , INTENT(IN) :: mix,err0
    REAL(8) , INTENT(OUT) :: err
    LOGICAL , INTENT(OUT) :: success

    INTEGER :: i,j,k,l,li,lj,n,io,jo,ib,jb,kb,lb,norbit,irc,nbase,nocc
    INTEGER, ALLOCATABLE :: tmap(:)
    REAL(8) , ALLOCATABLE :: arg(:),rhs(:),rv(:),aden(:),v1(:),v2(:),o(:,:)
    REAL(8) , ALLOCATABLE :: Xv(:,:),Xiv(:,:)
    TYPE(OrbitInfo), POINTER :: PSO
    TYPE(OrbitInfo) :: tmpOrbit
    REAL(8) :: occ,x,y,q,v0term,en,val
    INTEGER :: fcount=0
    REAL(8), PARAMETER :: threshold=1.d-8
    CHARACTER(4) :: stuff

    success=.FALSE.
    n=Grid%n; nbase=PAW%nbase; norbit=PAW%OCCWFN%norbit
    ALLOCATE(arg(n),rhs(n),rv(n),aden(n),v1(n),v2(n),tmap(norbit),&
&        o(norbit,nbase),Xv(n,norbit),Xiv(nbase,norbit))

    PAW%tkin=0; PAW%tion=0; PAW%tvale=0; PAW%txc=0; PAW%Ea=0; PAW%Etotal=0
    PSO=>PAW%TOCCWFN
    ! orthonormalize
    nocc=0
    DO io=1,PSO%norbit
       IF (.NOT.PSO%iscore(io)) THEN
          CALL  HFOrthotocore(Grid,PSO%l(io),PSO%wfn(:,io),PAW)
          IF (io>1) THEN
             DO jo=1,io-1
                IF (PSO%l(jo)==PSO%l(io).AND..NOT.PSO%iscore(jo)) THEN
                   CALL HFgenOrthog(Grid,PAW,PSO%l(io),&
&                       PSO%wfn(:,io),PSO%wfn(:,jo))
                   WRITE(6,*) 'orthog', io,jo
                ENDIF
             ENDDO
          ENDIF
          x=HFgenoverlap(Grid,PAW,PSO%l(io),PSO%wfn(:,io),PSO%wfn(:,io))
          PSO%wfn(:,io)=PSO%wfn(:,io)/SQRT(x)
          WRITE(6,*) 'normalize ', io, x
          nocc=nocc+1
          tmap(nocc)=io
       ENDIF
    ENDDO

    !   checking overlap
    do i=1,nocc
       io=tmap(i)
       do jo=1,PSO%norbit
          if (PSO%l(jo)==PSO%l(io))  then
             if (PSO%iscore(jo)) then
                x=CoreOverlap(Grid,PAW,PSO%l(io),jo,PSO%wfn(:,io))
             else
                x=HFgenoverlap(Grid,PAW,PSO%l(io),PSO%wfn(:,jo),PSO%wfn(:,io))
             endif
            write(6,*) 'overlap ',io,jo,x
         endif
       enddo
    enddo

    !  prepare overlap terms
    o=0
    DO k=1,nocc
       io=tmap(k); l=PSO%l(io)
       DO ib=1,PAW%nbase
          IF (l==PAW%l(ib)) THEN
             o(io,ib)=overlap(Grid,PSO%wfn(:,io),PAW%otp(:,ib))
             WRITE(6,'("<p|psi> ", 2i5,1p,1e15.7)') io,ib,o(io,ib)
          ENDIF
       ENDDO
    ENDDO

    CALL CopyOrbit(PSO,tmpOrbit)

    rv=PAW%rtVf; PAW%tkin=0;   PAW%wij=0;    PSO%lqp=0
    PAW%tden=0;rhs=0;aden=0;x=0
    DO k=1,nocc
       io=tmap(k) ; l=PSO%l(io)
       occ=PSO%occ(io)
       PAW%tden=PAW%tden+occ*(PSO%wfn(:,io))**2
       DO jo=1,PSO%norbit
          IF(PSO%l(jo)==l) THEN
             CALL kinetic_ij(Grid,PSO%wfn(:,io),PSO%wfn(:,jo),l,y)
             PSO%lqp(jo,io)=y
             !WRITE(6,*) 'Kinetic ',jo,io,y
             IF (io==jo) THEN
                PAW%tkin= PAW%tkin + occ*y
                !PSO%eig(io)=y
             ENDIF
          ENDIF
       ENDDO
       DO ib=1,PAW%nbase
          DO jb=1,PAW%nbase
             IF (PAW%l(ib)==l.AND.PAW%l(jb)==l) THEN
                x=x+occ*o(io,ib)*o(io,jb)*PAW%mLij(ib,jb,1)
                PAW%wij(ib,jb)=PAW%wij(ib,jb)+occ*o(io,ib)*o(io,jb)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    aden=PAW%tden+x*PAW%g(:,1)

    ! now have kinetic part of lqv and valence pseudodensity
    arg=0; arg(2:n)=(aden(2:n)*PAW%hatpot(2:n))/Grid%r(2:n)
    v0term=integrator(Grid,arg)
    WRITE(6,*) 'v0term, aden', v0term,integrator(Grid,aden)

    CALL poisson(Grid,q,aden,rhs,x)
    WRITE(6,*) 'PAW poisson ', q,x
    PAW%tvale=x
    arg=0; arg(2:n)=PAW%rtVf(2:n)/Grid%r(2:n)
    PAW%tion=overlap(Grid,arg,PAW%tden)
    rv=rv+rhs    ! vion + valence-Hartree

    DO k=1,nocc
       io=tmap(k);l=PSO%l(io)
       arg=rv*(PSO%wfn(:,io));
       arg(1)=0; arg(2:n)=arg(2:n)/Grid%r(2:n)
       DO jo=1,PSO%norbit
          IF (PSO%l(jo)==l) THEN
             x=overlap(Grid,PSO%wfn(:,jo),arg)
             !WRITE(6,*) ' potential term ',jo,io, x
             PSO%lqp(jo,io)= PSO%lqp(jo,io)+x
          ENDIF
       ENDDO
    ENDDO

    PAW%dij=0; PAW%Ea=0
    DO ib=1,PAW%nbase
       DO jb=1,PAW%nbase
          IF (PAW%l(ib)==PAW%l(jb)) THEN
             PAW%dij(ib,jb)=PAW%dij(ib,jb) + PAW%Kij(ib,jb) + &
&                 PAW%VFij(ib,jb)+PAW%mLij(ib,jb,1)*v0term
             PAW%Ea=PAW%Ea + PAW%wij(ib,jb)*(PAW%Kij(ib,jb) + &
&                 PAW%VFij(ib,jb))
             x=0;    !  accumulate Hartree term
             DO kb=1,PAW%nbase
                DO lb=1,PAW%nbase
                   IF (PAW%l(kb)==PAW%l(lb)) THEN
                      x=x+PAW%wij(kb,lb)*PAW%DR(ib,jb,kb,lb,1)
                   ENDIF
                ENDDO
             ENDDO
             PAW%dij(ib,jb)=PAW%dij(ib,jb)+x
             PAW%Ea=PAW%Ea + 0.5d0*PAW%wij(ib,jb)*x
             write(6,'(" Dij chk again ", 3I5,1p,1e20.13)')ib,jb,PAW%l(ib),&
&               PAW%dij(ib,jb)
          ENDIF
       ENDDO
    ENDDO

    ! valence Dij terms are completed; update valence lqp
    DO i=1,nocc
       io=tmap(i); li=PSO%l(io)
       DO j=1,nocc
          jo=tmap(j); lj=PSO%l(jo)
          !WRITE(6,*) 'lqv before dij ', jo,io,PSO%lqp(jo,io)
          DO ib=1,PAW%nbase
             IF (PAW%l(ib)==li) THEN
                DO jb=1,PAW%nbase
                   IF (PAW%l(jb)==lj) THEN
                      PSO%lqp(jo,io)=PSO%lqp(jo,io)+&
&                          PAW%dij(ib,jb)*o(io,ib)*o(jo,jb)
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
          !WRITE(6,*) 'lqv after dij ', jo,io,PSO%lqp(jo,io)
       ENDDO
    ENDDO

    !  update core lqp
    DO io=1,PSO%norbit
       IF (PSO%iscore(io)) THEN
          DO jb=1,PAW%nbase
             IF (PAW%l(jb)==PSO%l(io)) THEN
                x=PAW%Dcj(io,jb)+PAW%mLic(jb,io,1)*v0term
                DO kb=1,PAW%nbase
                   DO lb=1,PAW%nbase
                      if (PAW%l(kb)==PAW%l(lb))  &
&                        x=x+PAW%wij(kb,lb)*PAW%DRCjkl(io,jb,kb,lb,1)
                   ENDDO
                ENDDO
                DO j=1,nocc
                   jo=tmap(j); lj=PSO%l(jo)
                   IF (lj==PAW%l(jb)) THEN
                      PSO%lqp(io,jo)=PSO%lqp(io,jo)+x*o(jo,jb)
                      !WRITE(6,*) 'Update lqv ', io,jo,jb,PSO%lqp(io,jo)
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO

    ! setup exchange terms
    PAW%txc=0
    DO i=1,nocc
       io=tmap(i);l=PSO%l(io)
       CALL Calc_tXv(Grid,PAW,io,o,Xv(:,io),Xiv(:,io),x,y)
       WRITE(6,*) 'Xiv ',io, (Xiv(ib,io),ib=1,PAW%nbase)
       !WRITE(6,*) 'x y', x,y
       PAW%txc=PAW%txc+PSO%occ(io)*x
       PAW%Ea=PAW%Ea+PSO%occ(io)*y
       arg=0;   arg(2:n)=Xv(2:n,io)/Grid%r(2:n)
       DO jo=1,PSO%norbit
          IF (PSO%l(jo)==l) THEN
             x=overlap(Grid,PSO%wfn(:,jo),arg)
             PSO%lqp(jo,io)=PSO%lqp(jo,io)+x
             !WRITE(6,*) 'Smooth exchange ' , jo,io,x
             IF (PSO%iscore(jo)) THEN
                CALL Calc_Xc(Grid,PAW,io,jo,o,x)
                PSO%lqp(jo,io)=PSO%lqp(jo,io)+x
                !WRITE(6,*) 'One center exchange',x
             ELSE
                x=0
                DO ib=1,PAW%nbase
                   IF(PAW%l(ib)==l) x=x+o(jo,ib)*Xiv(ib,io)
                ENDDO
                !WRITE(6,*) 'One center exchange',x
                PSO%lqp(jo,io)=PSO%lqp(jo,io)+x
             ENDIF
             WRITE(6,'("lqp jo,io ", 2i5,1p,1e15.7)') jo,io, PSO%lqp(jo,io)
          ENDIF
       ENDDO
       PSO%eig(io)=PSO%lqp(io,io)
    ENDDO
     !do i=1,Grid%n
     !    write(601,'(1P15e15.7)') Grid%r(i),(Xv(i,tmap(k)),k=1,nocc)
     !enddo

    PAW%Etotal=PAW%tkin+PAW%tion+PAW%tvale+PAW%txc+PAW%Ea
    WRITE(6,*) '*******Total energy*********', PAW%Etotal
    WRITE(6,*) 'PAW%tkin                    ',PAW%tkin
    WRITE(6,*) 'PAW%tion                    ',PAW%tion
    WRITE(6,*) 'PAW%tvale                   ',PAW%tvale
    WRITE(6,*) 'PAW%txc                     ',PAW%txc
    WRITE(6,*) 'PAW%Ea                      ',PAW%Ea

    ! Solve inhomogeneous diffeq. and store result in tmpOrbit
    err=0;
    DO k=1,nocc
       io=tmap(k); l=PSO%l(io)
       en=PSO%eig(io)
       arg=0; arg(2:n)=Xv(2:n,io)/Grid%r(2:n)
       rhs=-arg-en*PSO%wfn(:,io)
       DO jo=1,PSO%norbit
          if (io==jo.or.(PSO%l(jo)==l.and.PSO%occ(jo)>threshold)) &
&                  rhs=rhs+PSO%lqp(jo,io)*PSO%wfn(:,jo)
       ENDDO
       DO ib=1,PAW%nbase
          IF(PAW%l(ib)==l) THEN
             x=-Xiv(ib,io)
             DO jb=1,PAW%nbase
                IF (PAW%l(jb)==l) x=x-PAW%dij(ib,jb)*o(io,jb)
             ENDDO
             DO jo=1,PSO%norbit
                if (io==jo.or.(PSO%l(jo)==l.and.PSO%occ(jo)>threshold)) then
                   IF (PSO%iscore(jo)) THEN
                       x=x+PSO%lqp(jo,io)*PAW%mLic(ib,jo,1)
                   ELSE
                       DO jb=1,PAW%nbase
                          IF(PAW%l(jb)==l) THEN
                            x=x+PSO%lqp(jo,io)*PAW%oij(ib,jb)*o(jo,jb)
                          ENDIF
                       ENDDO
                   ENDIF
                endif
             ENDDO
          rhs=rhs+PAW%otp(:,ib)*x
          ENDIF
       ENDDO

       tmpOrbit%wfn(:,io)=0   ; rhs=-rhs
       !  note rhs is -(what we expect)
       CALL inhomo_bound_numerov(Grid,l,n,en,rv,rhs,tmpOrbit%wfn(:,io))
       arg=(PSO%wfn(:,io)-tmpOrbit%wfn(:,io))**2
       err=err+PSO%occ(io)*Integrator(Grid,arg)
       !do i=1,Grid%n
       !   write(800+k,'(1p,20e15.7)') Grid%r(i),PSO%wfn(i,io),&
       !&     tmpOrbit%wfn(i,io),rhs(i),rv(i)
       !enddo
    ENDDO

    !stop
    WRITE(6,*) 'PAWIter ', fcount,err
    call mkname(fcount,stuff)
     open(1001, file='hfpswfn.'//stuff,form='formatted')
       do i=1,Grid%n
          write(1001,'(1p,20e15.7)') Grid%r(i),(tmpOrbit%wfn(i,io),&
&                   PSO%wfn(i,io),io=1,PSO%norbit)
       enddo
     close(1001)

    ! update wfn if tolerance not satisfied
    IF (err>err0) THEN
       val=(1.d0-mix)
       WRITE(6,*) 'mixing wfns ', val
       DO k=1,nocc
          io=tmap(k)
          PSO%wfn(:,io)=val*PSO%wfn(:,io)+mix*tmpOrbit%wfn(:,io)
       ENDDO
    ELSE
       success=.TRUE.
    ENDIF
    fcount=fcount+1
    DEALLOCATE(arg,rhs,rv,aden,v1,v2,tmap,o,Xv,Xiv)
    CALL DestroyOrbit(tmpOrbit)
  END SUBROUTINE PAWIter_HF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Calc_tXp(Grid,PAW,li,wfn,twfn,tres)
  !      Calculate \tilde{X}_p for continuum pseudostate wfn
  !        Used for finding local pseudopotential in HF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Calc_tXp(Grid,PAW,ni,li,wfn,twfn,tres,lng)
    TYPE(GridInfo), INTENT(IN):: Grid
    TYPE(PseudoInfo), INTENT(IN):: PAW
    INTEGER, INTENT(IN) :: ni,li
    REAL(8), INTENT(IN) :: wfn(:),twfn(:)
    REAL(8), INTENT(INOUT) :: tres(:)
    INTEGER, OPTIONAL, INTENT(IN) :: lng

    REAL(8), POINTER :: r(:)
    REAL(8), ALLOCATABLE :: wfp(:),vl(:),dum(:),arg(:),f(:),hat(:,:)
    INTEGER :: i,j,k,l,n,m,ll,lj,lmin,lmax,jo,ok,norbit,ko,ib,jb,kb,lb,last
    INTEGER :: nu,nup,nj
    REAL(8) :: occ,term,term1,wgt,occi,occj,rc,uvjo,fac,x
    REAL(8), PARAMETER :: threshold=1.d-8
    LOGICAL :: iszero


    norbit=PAW%TOCCwfn%norbit
    n=Grid%n

    write(6,*) 'In Calc_tXp ', n,norbit; call flush(6)
    ll=MAXVAL(PAW%TOCCWFN%l(:)); ll=li+ll

    ALLOCATE(wfp(n),vl(n),dum(n),arg(n),f(n),hat(n,ll+1))
    r=>Grid%r

    write(6,*) 'After alloc ',n;call flush(6)
    tres(:)=0.d0;
    occ=1.d0

    DO l=0,ll
       CALL hatL(Grid,PAW,l,hat(:,l+1))
       f=hat(:,l+1)*(Grid%r(:)**l)
       write(6,*) 'test hat l', l,integrator(Grid,f)
    ENDDO

    Open(33,file='hattest',form='formatted')
    do i=1,n
      write(33,'(1p,50e15.7)') Grid%r(i),(hat(i,l+1),l=0,ll)
    enddo
    close(33)

    last=n
    if (present(lng)) last=min(last,lng)
    r=>Grid%r

    DO jo=1,norbit
       write(6,*) 'jo : ', jo; call flush(6)
       occj=PAW%OCCwfn%occ(jo); vl=0.d0
       IF (occj>threshold) THEN
          lj=PAW%OCCwfn%l(jo)
          lmax=li+lj
          lmin=ABS(li-lj)
          wfp(1:last)=twfn(1:last)*PAW%TOCCwfn%wfn(1:last,jo)
          arg(1:last)=wfn(1:last)*PAW%OCCwfn%wfn(1:last,jo) - &
&              twfn(1:last)*PAW%TOCCwfn%wfn(1:last,jo)
          DO ll=lmin,lmax,2
             CALL EXXwgt(occ,occj,ni,li,jo,lj,ll,wgt)
             IF (wgt>threshold) THEN
                wgt=wgt*(2*ll+1)   ! because of apoisson convention
                f(1:last)=arg(1:last)*(r(1:last)**ll)
                x=integrator(Grid,f,1,last); write(6,*) 'moment',ll,x;
                call flush(6)
                f(1:last)=wfp(1:last)+x*hat(1:last,ll+1)
                CALL apoisson(Grid,ll,last,f,dum)
                vl(1:last)=vl(1:last)- &
&                   wgt*dum(1:last)*PAW%TOCCwfn%wfn(1:last,jo)/occ
             ENDIF
          ENDDO
       tres=tres+vl
       ENDIF
    ENDDO   !jo

    write(6,*) 'Finished tXP '; call flush(6)

    DEALLOCATE(wfp,vl,dum,arg,f,hat)
  END SUBROUTINE Calc_tXp

  !***************************************************************
  ! SUBROUTINE troullier_HF(lmax,Grid,Pot)
  !  Creates  screened norm-conserving pseudopotential following
  !    approach of N. Troullier and J. L. Martins, PRB 43, 1993 (1991)
  !    Uses p(r)=a0+f(r); f(r)=SUMm(Coef(m)*r^(2*m), where
  !          m=1,2..6
  !    Psi(r) = r^(l+1)*exp(p(r))
  !      HF version
  !***************************************************************
  SUBROUTINE Troullier_HF(Grid,Pot,PAW,l,e)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
    INTEGER,INTENT(IN) :: l
    REAL(8),INTENT(IN) :: e

    REAL(8) :: A0,A,B,B0,C,C0,D,F,S
    REAL(8) :: Coef(6),Coef0,Coef0old
    REAL(8) :: h,rc,delta,x,pp,dpp,ddpp,dddpp,ddddpp
    REAL(8) :: gam,bet
    INTEGER :: i,j,k,n,iter,nr,nodes,irc,ok,m,wavetype,lng
    INTEGER, PARAMETER :: niter=5000,ni=999
    REAL(8), PARAMETER :: small=1.0d-9
    REAL(8), ALLOCATABLE ::  VNC(:),XX(:),twfn(:),wfn(:),p(:),dum(:),tX(:)
    REAL(8), ALLOCATABLE :: rvx(:),trvx(:)
    REAL(8), POINTER :: r(:),rv(:)
    CHARACTER(132) :: line

    n=Grid%n
    h=Grid%h
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    write(6,*) 'nr = ', nr , n, PAW%irc_vloc
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc
    write(6,*) ' Check rc', irc,rc,Grid%r(irc)

    ALLOCATE(VNC(n),wfn(n),twfn(n),p(n),dum(n),XX(n),tX(n),&
&                 rvx(n),trvx(n),stat=ok)
    IF (ok /=0) THEN
       WRITE(6,*) 'Error in troullier  -- in allocating wfn,p', nr,ok
       STOP
    ENDIF

    !write(6,*) ' Troullier ', n,nr,irc
    !call flush(6)
    !if (scalarrelativistic) then
    !   CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    !else
    !   CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
    !endif

    If (PAW%exctype=='HF') then
        CALL HFunocc(Grid,PAW%OCCWFN,l,e,Pot%rv,Pot%v0,Pot%v0p,wfn,lng,XX)
    Else
       Write(6,*) 'Error -- this version of Troulliers is for HF only'
       stop
    endif

    rvx(1)=0.d0; rvx(2:lng)=XX(2:lng)*Grid%r(2:lng)/wfn(2:lng)
    VNC(1:lng)=Pot%rv(1:lng)+rvx(1:lng)

    !do i=1,lng
    !    write(701,'(1P8E16.7)') Grid%r(i),wfn(i),XX(i),Pot%rv(i),VNC(i)
    !enddo
    IF (wfn(irc)<0) then
             wfn=-wfn;
     EndIF
    dum(1:irc)=(wfn(1:irc)**2)
    S=integrator(Grid,dum(1:irc),1,irc)
    A0=LOG(wfn(irc)/(rc**(l+1)))
    B0=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
    C0=rc*(VNC(irc)-rc*e)-B0*(B0+2*l+2)
    D=-rc*(VNC(irc)-rc*Gfirstderiv(Grid,irc,VNC))-2*B0*C0-2*(l+1)*(C0-B0)
    F=rc*(2*VNC(irc)-rc*(2*Gfirstderiv(Grid,irc,VNC) &
&        -rc*Gsecondderiv(Grid,irc,VNC)))+&
&        4*(l+1)*(C0-B0)-2*(l+1)*D-2*C0**2-2*B0*D

    WRITE(6,*) 'In troullier -- matching parameters',S,A0,B0,C0,D,F

    delta=1.d10
    iter=0
    Coef0=0

    DO WHILE(delta>small.AND.iter<=niter)
       iter=iter+1
       A=A0-Coef0
       B=B0
       C=C0
       CALL EvaluateTp(l,A,B,C,D,F,coef)

       dum=0
       DO  i=1,irc
          x=(r(i)/rc)**2
          p(i)=x*(Coef(1)+x*(Coef(2)+x*(Coef(3)+&
&              x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
          dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
       ENDDO
       Coef0old=Coef0

       x=integrator(Grid,dum(1:irc),1,irc)
       Coef0=(LOG(S/x))/2

       delta=ABS(Coef0-Coef0old)
       !WRITE(6,'(" VNC: iter Coef0 delta",i5,1p,2e15.7)') iter,Coef0,delta
    ENDDO

    WRITE(6,*) '  VNC converged in ', iter,'  iterations'
    WRITE(6,*) '  Coefficients  -- ', Coef0,Coef(1:6)
    !
    ! Now  calculate VNC
    OPEN(88,file='NC',form='formatted')
    !
    dum=VNC; VNC=0
    DO  i=2,nr
       x=(r(i)/rc)**2
       p(i)=Coef0+x*(Coef(1)+x*(Coef(2)+&
&           x*(Coef(3)+x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
       dpp=2*r(i)/(rc**2)*(Coef(1)+x*(2*Coef(2)+x*(3*Coef(3)+&
&           x*(4*Coef(4)+x*(5*Coef(5)+x*6*Coef(6))))))
       ddpp=(1/(rc**2))*(2*Coef(1)+x*(12*Coef(2)+x*(30*Coef(3)+&
&           x*(56*Coef(4)+x*(90*Coef(5)+x*132*Coef(6))))))
       dddpp=(r(i)/rc**4)*(24*Coef(2)+x*(120*Coef(3)+x*(336*Coef(4)+&
&           x*(720*Coef(5)+x*1320*Coef(6)))))
       ddddpp=(1/(rc**4)*(24*Coef(2)+x*(360*Coef(3)+x*(1680*Coef(4)+&
&           x*(5040*Coef(5)+x*11880*Coef(6))))))
       IF (i==irc) THEN
          WRITE(6,*) 'check  dp ', dpp,  B0/rc
          WRITE(6,*) 'check ddp ', ddpp, C0/rc**2
          WRITE(6,*) 'check dddp', dddpp, D/rc**3
          WRITE(6,*) 'check ddddp', ddddpp, F/rc**4
       ENDIF
       VNC(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
       twfn(i)=(r(i)**(l+1))*EXP(p(i))
       WRITE(88,'(1p,6e15.7)') r(i),wfn(i),twfn(i),VNC(i)*r(i),rv(i),dum(i)
    ENDDO
    CLOSE(88)
    x=overlap(Grid,twfn(1:irc),twfn(1:irc),1,irc)
    WRITE(6,*) 'check norm ',x,S; call flush(6)

    twfn(irc:lng)=wfn(irc:lng)
    call Calc_tXp(Grid,PAW,ni,l,wfn,twfn,tX,lng)
    trvx(2:lng)=tX(2:lng)/twfn(2:lng)
    call extrapolate(Grid,trvx)

    !   need to remove exchange part
    PAW%rveff=Pot%rv
    PAW%rveff(1:irc)=VNC(1:irc)*r(1:irc)-trvx(1:irc)

    open(88,file='checkvxc',form='formatted')
    do i=1,lng
       write(88,'(1p,50e15.7)') r(i),XX(i)*r(i),tX(i),rvx(i),trvx(i),VNC(i),&
&              rv(i),PAW%rveff(i)
    enddo
    close(88)

    DEALLOCATE(VNC,wfn,p,dum,XX,twfn,tX,rvx,trvx)
  END SUBROUTINE troullier_HF

END MODULE hf_pseudo

