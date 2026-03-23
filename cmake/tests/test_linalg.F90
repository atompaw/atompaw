program test_linalg
  implicit none
  
  integer, parameter :: n = 2
  real(8) :: A(2,2), b(2), x(2)
  integer :: ipiv(2), info
  
  A = reshape([1.0d0, 2.0d0, 3.0d0, 4.0d0], [2,2])
  b = [5.0d0, 6.0d0]
  
  call dcopy(n, b, 1, x, 1)
    call dgesv(n, 1, A, n, ipiv, x, n, info)
  
  stop (info == 0)  ! 0=OK, 1=échec

end program test_linalg
