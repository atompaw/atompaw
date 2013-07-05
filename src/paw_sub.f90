Module paw_sub
   USE Globalmath
   USE gridmod

   IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Support routines for PAW calculations;
!     Several of these routines were written by Marc Torrent of CEA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CONTAINS
!******************************************************************************
! Pseudization routine: PSPOLYN
!  Pseudize a function with a polynom
!                   tfunc(r)=r^(l+1).Sum[Ci.r^2i]  0<=i<=np-1  if r<=rc
!                   tfunc(r)=func(r)                           if r>rc
!  Ci coefficients are returned
!******************************************************************************

 SUBROUTINE pspolyn(func,Ci,r,l,np,irc,n)

   INTEGER,INTENT(IN) :: irc,l,n,np
   REAL(8),INTENT(IN) :: func(n),r(n)
   REAL(8),INTENT(OUT) :: Ci(np)

   INTEGER :: i,j,np2
   REAL(8) :: rc,xx,y,scale
   REAL(8),ALLOCATABLE :: A(:,:),X(:)

   if (irc<3.or.irc>irc.or.irc>n-3) stop 'pspolyn: rc out of range'
   if (np<1) stop 'pspolyn: p out of range'

   allocate(A(np,np),X(np),stat=i)
   if (i/=0) stop 'allocation error in pspolyn'

   rc=r(irc);np2=np/2
   scale=(rc/(rc-r(irc-1)))**2 ! Scale to limit rounding error in linsol
   do i=1,np
    xx=r(i+irc-np2-1)/rc
    y=xx*xx
    A(i,1)=scale
    do j=2,np
     A(i,j)=A(i,j-1)*y
    enddo
    X(i)=scale*func(i+irc-np2-1)/(xx**(l+1))
   enddo

   call linsol(A,X,np,np,np,np)
   !write(6,*) 'Completed linsol with coefficients'
   !write(6,'(1p,10e15.7)') (X(i),i=1,np)

   do i=1,np
    Ci(i)=X(i)/rc**(l+2*i-1)
   enddo

   deallocate(A,X)

 END SUBROUTINE pspolyn

!******************************************************************************
! Pseudization routine: PSUSPOLYN
!  Pseudize a function with a polynom
!                   tfunc(r)=r^(l+1).Sum[Ci.r^2i]  0<=i<=np-1  if r<=rc
!                   tfunc(r)=func(r)                           if r>rc
!      For i>3, Ci coefficients are computed so that to minimize
!      Fourier coefficients of pseudized function for q>qcut
!  Ci coefficients are returned
!******************************************************************************

 SUBROUTINE psuspolyn(func,Ci,r,l,np,irc,n,qcut)

   INTEGER,INTENT(IN) :: irc,l,n,np
   REAL(8),INTENT(IN) :: qcut
   REAL(8),INTENT(IN) :: func(n),r(n)
   REAL(8),INTENT(OUT) :: Ci(np)

   INTEGER,PARAMETER :: nq=2001
   REAL(8),PARAMETER :: lfact(0:5)=(/1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/)
   INTEGER :: i,j,k,ip,iq,jq,il,ix,jx
   REAL(8) :: qh,qrc2,rc,xx,scale,yy(6,6),zz(6)
   REAL(8),ALLOCATABLE :: A(:,:),X(:),Q(:,:),qq(:),ff(:),gg(:)

   if (irc<3.or.irc>irc.or.irc>n-6) stop 'psuspolyn: rc out of range'
   if (np<4) stop 'psuspolyn: polynomial degree for pseudization was too small'

   allocate(A(np+4,np+4),X(np+4),Q(nq,np+4),qq(nq),ff(nq),gg(nq),stat=i)
   if (i/=0) stop 'allocation error in psuspolyn'

   il=l+1;rc=r(irc)

   qh=qcut/dfloat(nq-1)
   do iq=1,nq
    qq(iq)=0.5d0*qcut+qh*dble(iq-1)
   enddo

   Q(:,:)=0.d0
   do ip=1,np
    do iq=1,nq
     qrc2=0.5d0*(qq(iq)*rc)**2
     xx=1.d0;ix=0
     do while (abs(xx)>1.d-20)
      ix=ix+1
      xx=xx*qrc2/dble(ix)/dble(2*(ix+l)+1)
     enddo
     xx=0.d0
     do jx=ix,1,-1
      xx=qrc2/dble(jx)/dble(2*(jx+l)+1)*(1.d0/dble(2*(jx+ip+l)+1)-xx)
     enddo
     Q(iq,ip)=4.d0*pi*rc**dble(2*ip+il)*(qq(iq)*rc)**dble(l)/lfact(l) &
&                 *(1.0/dble(2*(l+ip)+1)-xx)
    enddo
   enddo

   A(:,:)=0.d0
   do iq=1,np
    do ip=1,np
     do jq=1,nq
      ff(jq)=qq(jq)**4*Q(jq,iq)*Q(jq,ip)
     enddo
     A(iq,ip)=overint(nq,qh,ff)
    enddo
   enddo
   do iq=1,np
    ix=2*(iq-1)+il
    A(np+1,iq)=rc**dble(ix)
    A(np+2,iq)=dble(ix)*rc**(ix-1)
    A(np+3,iq)=dble(ix*(ix-1))*rc**(ix-2)
    A(np+4,iq)=dble(ix*(ix-1)*(ix-2))*rc**(ix-3)
    A(iq,np+1)=A(np+1,iq)
    A(iq,np+2)=A(np+2,iq)
    A(iq,np+3)=A(np+3,iq)
    A(iq,np+4)=A(np+4,iq)
   enddo

   yy(:,:)=0.d0;zz(:)=0.d0
   do i=1,6
    do j=1,6
     if (i==1.and.j==1) then
      yy(i,j)=12.d0
     else
      do k=1,12
       yy(i,j)=yy(i,j)+(r(irc+k-6)-rc)**(i+j-2)
      enddo
     endif
    enddo
   enddo
   do k=1,12
    zz(1)=zz(1)+func(irc+k-6)
   end do
   do i=2,6
    do k=1,12
     zz(i)=zz(i)+func(irc+k-6)*(r(irc+k-6)-rc)**(i-1)
    end do
   end do
   scale=1/(rc-r(irc-1))**3;yy=yy*scale;zz=zz*scale ! Scale to limit rounding error in linsol
   call linsol(yy,zz,6,6,6,6)
   zz(3)=2.d0*zz(3);zz(4)=6.d0*zz(4)
   X(np+1:np+4)=zz(1:4)

   do iq=1,nq
    gg(iq)=4.d0*pi*intjl(rc,qq(iq),zz,l)
   enddo
   do ip=1,np
    do iq=1,nq
     ff(iq)=qq(iq)**4*gg(iq)*Q(iq,ip)
    enddo
    X(ip)=-overint(nq,qh,ff)
   enddo

   scale=(rc/(rc-r(irc-1)))**2;A=A*scale;X=X*scale ! Scale to limit rounding error in linsol
   call linsol(A,X,np+4,np+4,np+4,np+4)
   !write(6,*) 'Completed linsol with coefficients'
   !write(6,'(1p,10e15.7)') (X(i),i=1,np+4)

   Ci(1:np)=X(1:np)

   deallocate(A,X,Q,qq,ff,gg)

 END SUBROUTINE psuspolyn

!******************************************************************************
! Pseudization routine: PSBES
!  Pseudize a function with a sum of 2 Bessel functions
!                (following PHYS REV B 41,1227 (1990))
!                   tfunc(r)=[al(1)*jl(ql(1)*r)+al(2)*jl(ql(2)*r)]*r if r<=rc
!                   tfunc(r)=func(r)                                 if r>rc
!  al and ql coefficients are returned
!******************************************************************************

 SUBROUTINE psbes(func,al,ql,Grid,l,irc,n)

   INTEGER,INTENT(IN) :: irc,l,n
   REAL(8),INTENT(IN) :: func(n)
   REAL(8),INTENT(OUT) :: al(2),ql(2)
   TYPE(GridInfo),INTENT(IN) :: Grid

   INTEGER :: i
   REAL(8) :: alpha,beta,det,qr,jbes,jbesp,jbespp,rc
   REAL(8) :: amat(2,2),bb(2)

   rc=Grid%r(irc)

   beta=1.D0
   alpha=1.D0-Gfirstderiv(Grid,irc,func)*rc/func(irc)
   call solvbes(ql,alpha,beta,l,2)
   ql(1:2)=ql(1:2)/rc

   do i=1,2
    qr=ql(i)*rc
    call jbessel(jbes,jbesp,jbespp,l,2,qr)
    jbespp=2.d0*ql(i)*jbesp+jbespp*ql(i)*ql(i)*rc
    jbesp=jbes+jbesp*ql(i)*rc
    jbes=jbes*rc
    amat(1,i)=jbes
    amat(2,i)=jbespp
   enddo

   bb(1)=func(irc)
   bb(2)=Gsecondderiv(Grid,irc,func)

   det=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
   al(1)=(amat(2,2)*bb(1)-amat(1,2)*bb(2))/det
   al(2)=(amat(1,1)*bb(2)-amat(2,1)*bb(1))/det

 END SUBROUTINE psbes

 SUBROUTINE trunk(Grid,f,rstart,rend)
   TYPE (GridInfo), INTENT(IN) :: Grid
   REAL(8),  INTENT(INOUT) :: f(:)
   REAL(8),  INTENT(IN) :: rstart,rend

   INTEGER :: i,n
   REAL(8) :: r,delta,arg

   delta=rend-rstart

   DO i=1,Grid%n
      r=Grid%r(i)
      IF (r>rstart .AND. r<=rend) THEN
         arg=pi*(r-rstart)/delta
         f(i)=f(i)*(SIN(arg)/arg)**2
      ENDIF
      IF (r>rend) f(i)=0
   ENDDO
 END SUBROUTINE trunk

END MODULE paw_sub
