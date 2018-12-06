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
  integer ipiv(size(b)),i,j,neq,info_tmp
  real(kind=WP) y(size(b),1)
  
  neq = size(b)
  info = 0
  y(:,1) = b
  call ssytrf_f95(a,ipiv=ipiv,info=info_tmp)
  if (info_tmp>0) then
    if (present(info)) info = info_tmp
    return
  end if
  call ssytrs_f95(a,y,ipiv=ipiv,info=info_tmp)
  if (info_tmp>0) then
    if (present(info)) info = info_tmp
    return
  end if
  b = y(:,1)

end subroutine

end module
