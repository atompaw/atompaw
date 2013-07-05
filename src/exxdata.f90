Module exxdata

  Implicit none

  TYPE HSZInfo
     REAL(8), POINTER :: psi(:,:),psiref(:,:)
     REAL(8), POINTER :: shift(:),grad(:),U(:),Uref(:),Ucore(:),Uvale(:)
     REAL(8), POINTER :: LMBD(:,:),LMBDref(:,:),LMBDcore(:,:),LMBDvale(:,:)
     LOGICAL :: Fixed_Zero  ! If true, read in constant zero_index
     ! If false (default) adjust zero_index each iteration
     INTEGER :: zero_index      ! index of most extended wave function
     INTEGER :: lmax    ! maximum l value of bound states
     REAL(8) :: betaL   ! expontial decay of zero_index wfn
     REAL(8) :: grad2   ! Dot_Product(grad,grad)
     INTEGER :: matchpoint  ! beyond this grid point use asymptotic form
     REAL(8), POINTER :: rVxref(:),coreshift(:)
     REAL(8), POINTER :: rVxKLI(:),rDVxKLI(:)
  END TYPE HSZInfo

  TYPE(HSZInfo), TARGET :: HSZ

End module
