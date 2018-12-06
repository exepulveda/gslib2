module solvers
use kinds_kit
implicit none
contains
subroutine system_solver(a,b,symmetric,info)
use lapack95
implicit none
  INTEGER, PARAMETER :: WP = KIND(1.0E0)
  real(kind=WP), intent(inout) :: a(:,:),b(:)
  logical, intent(in) :: symmetric
  integer, intent(out), optional :: info
  !locals
  info = 0
  call ssytrf_f95(a)
  !call ssytrs_f95(a,b)
end subroutine

end module
