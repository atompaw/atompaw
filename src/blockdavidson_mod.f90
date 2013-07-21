MODULE BlockDavidson_mod
  USE globalmath
  USE search_sort

  IMPLICIT NONE

  REAL(8), ALLOCATABLE, PRIVATE :: A(:,:),O(:,:),v(:,:),f(:,:),w(:),hv(:),ov(:)
  REAL(8), PARAMETER, PRIVATE :: base_eps=1.d-6
  REAL(8), PRIVATE :: eps
  REAL(8), PARAMETER, PRIVATE :: conv1=4.d13,conv2=3.d13,conv3=2.d13,conv4=1.d13
  INTEGER, PARAMETER, PRIVATE :: mxiter=300

CONTAINS

  SUBROUTINE InitBlockDavidson(nvec,vec,dup)
    INTEGER, INTENT(IN) :: nvec,dup
    REAL(8), INTENT(IN) :: vec(:,:)

    INTEGER :: ndim,i,ns

    ns=nvec*dup
    ndim=SIZE(vec,1)
    eps=MIN(SQRT(base_eps),base_eps*REAL(ndim))

    ALLOCATE(A(ns,ns),O(ns,ns),w(ns),v(ns,ns), &
&        f(ndim,ns),hv(ndim),ov(ndim), stat=i)
    IF (i/=0) THEN
       WRITE(6,*) 'Error in InitBlock...', ns,ndim,i
       STOP
    ENDIF
  END SUBROUTINE InitBlockDavidson

  SUBROUTINE EndBlockDavidson
    DEALLOCATE(A,O,v,f,w,hv,ov)
  END SUBROUTINE EndBlockDavidson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  On input vec(ndim,nvec) contains initial guesses for eigenfunctions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE BlockDavidson(ntest,nvec,vec,eig,success,dup,ressub,multsub)
    INTEGER, INTENT(IN) :: ntest,nvec,dup
    REAL(8), INTENT(INOUT) :: vec(:,:),eig(:)
    LOGICAL :: success

    INTEGER :: ndim,i,j,k,n,ns,start,finish,last,iter
    REAL(8) :: delta,ee,v1,v2,v3,v4

    INTERFACE
       SUBROUTINE ressub(vin,hv,ov,ee)
         REAL(8), INTENT(IN) :: vin(:)
         REAL(8), INTENT(OUT) :: hv(:),ov(:),ee
       END SUBROUTINE ressub
    END INTERFACE

    INTERFACE
       FUNCTION multsub(v1,v2)
         REAL(8) :: multsub
         REAL(8), INTENT(IN) :: v1(:),v2(:)
       END FUNCTION multsub
    END INTERFACE

    success=.false.
    !WRITE(6,*) 'In BlockDavidson ', nvec,dup
    !CALL flush(6)
    CALL InitBlockDavidson(nvec,vec,dup)
    delta=1.d10; f=0; eig=0;   ns=SIZE(A,1)
    v1=conv1;v2=conv2;v3=conv3;v4=conv4
    DO iter=1,mxiter
       CALL shift4(v1,v2,v3,v4,delta)
       IF (iter>=4.AND.(v4<eps.AND.v4>v3))THEN
          CALL EndBlockDavidson
          success=.true.
          RETURN

       ELSE
          start=1; finish=nvec; last=finish

          DO i=1,nvec
             f(:,i)=vec(:,i)
             CALL ressub(f(:,i),hv,ov,eig(i))
          ENDDO

          DO k=1,dup-1
             last=finish
             DO i=start,finish
                CALL ressub(f(:,i),hv,ov,ee)
                last=last+1
                f(:,last)=hv-ee*ov
                f(:,last)=f(:,last)/SQRT(multsub(f(:,last),f(:,last)))
             ENDDO
             start=finish+1; finish=last
          ENDDO

          !open (8,file='starting', form='formatted')
          !do i=1,Size(hv)
          !   write(8,'(i5, 1p,50e15.7)') i,(f(i,j),j=1,finish)
          !enddo
          !close(8)
          !stop

          A=0; O=0
          DO i=1,finish
             CALL ressub(f(:,i),hv,ov,ee)
             DO j=1,finish
                A(j,i)=multsub(f(:,j),hv)
                O(j,i)=multsub(f(:,j),ov)
             ENDDO
          ENDDO

          CALL Diagonalizer(finish,ns,i,A,O,w,v)

          IF (i < nvec) THEN
             WRITE(6,*) 'Too few eigenvalues',i,nvec
             STOP
          ENDIF

          delta=0
          DO i=1,nvec
             IF(i<=ntest) delta=delta+ABS(w(i)-eig(i))
             Eig(i)=w(i)
          ENDDO

          Vec=0
          DO i=1,nvec
             Eig(i)=w(i)
             DO j=1,finish
                Vec(:,i)=Vec(:,i)+f(:,j)*v(j,i)
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL EndBlockDavidson
    WRITE(6,*) ' BlockDavidson did not converge '

  END SUBROUTINE BlockDavidson


  SUBROUTINE Diagonalizer(VecSize,ArraySize,NewSize,Hbase,Obase,Eigen,Vec)
    INTEGER,          INTENT(IN)  :: VecSize
    INTEGER,          INTENT(IN)  :: ArraySize
    INTEGER,          INTENT(OUT) :: NewSize
    REAL(8),          INTENT(INOUT) :: Hbase(:,:)
    REAL(8),          INTENT(INOUT) :: Obase(:,:)
    REAL(8),          INTENT(OUT) :: Eigen(:)
    REAL(8),          INTENT(OUT) :: Vec(:,:)

    INTEGER ::  i, ii, j, k, LWork,  LSize
    INTEGER ::  Info
    REAL(8), ALLOCATABLE  :: C(:,:),U(:,:),VT(:,:),WORK(:),S(:)
    INTEGER, ALLOCATABLE  :: LUT(:)
    REAL(8)     :: tol,val

    LWork = MAX(200,VecSize**2)
    k=VecSize
    ALLOCATE(C(k,k),U(k,k),VT(k,k),WORK(LWORK),S(k),LUT(k),STAT=i)
    IF (i /= 0) THEN
       WRITE(6,*) 'Diagonalizer: allocation ', k,Lwork,i
       STOP
    ENDIF

    tol=1.d-8

    NewSize=VecSize

    C=0

    C(1:VecSize,1:VecSize)=Obase(1:VecSize,1:VecSize)

    CALL DGESVD('A','A',k,k,C,k,S,U,k,VT,k,WORK,LWORK,i)

    ii=0
    DO i=1,VecSize
       IF (S(i)>tol) THEN
          ii=ii+1
       ELSE
          EXIT
       ENDIF
    ENDDO
    IF (ii>0) THEN
       NewSize=ii
    ELSE
       WRITE(6,*) 'Error in Diagonalizer -- Obase is singular '
       STOP
    ENDIF

    C=0 ;
    DO i=1,NewSize
       DO j=1,Newsize
          C(i,j)=C(i,j)+&
&              DOT_PRODUCT(U(:,i),MATMUL(Hbase,VT(j,:)))/S(i)
       ENDDO
    ENDDO

    CALL DGEEV('N', 'V', NewSize, C(1,1), ArraySize, Eigen, &
&        S, U, ArraySize, U, ArraySize, Work, LWork, Info)

    !WRITE(6,*) ' completed Hmat diagonalization with Info=',Info
    !WRITE(6,*) 'h ', Eigen(1:NewSize)

    !WRITE(6,*) ' completed Hmat diagonalization with Info=',Info
    IF (info /= 0) THEN
       WRITE(6,*) 'Stopping due to diagonalizer error'
       STOP
    ENDIF

    Work=1.e10
    Work(1:NewSize)=Eigen(1:NewSize)
    CALL Insertion_Sort(work(1:Newsize),LUT(1:Newsize),.TRUE.)
    Eigen(1:Newsize)=work(LUT(1:NewSize))
    Vec=0
    DO i=1,NewSize
       Vec(1:VecSize,i)=&
&           MATMUL(TRANSPOSE(VT(1:NewSize,1:VecSize)),U(1:NewSize,LUT(i)))
    ENDDO

    DEALLOCATE(C,U,VT,WORK,S,LUT)

  END SUBROUTINE  Diagonalizer


END MODULE BlockDavidson_mod
