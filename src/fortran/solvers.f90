module solvers
use kinds_kit
implicit none
contains
subroutine system_solver(a,b,x,symmetric,info)
use lapack95
implicit none
  real(kind=fk), intent(in) :: a(:,:),b(:)
  real(kind=fk), intent(out) :: x(:)
  logical, intent(in) :: symmetric
  integer, intent(out), optional :: info
  !locals
  real(kind=8) int_a(size(b),size(b))
  integer ipiv(size(b)),i,j,neq,info_tmp
  real(kind=8) y(size(b),1)
  logical failed
  
  int_a = real(a,kind=8)
  y(:,1) = real(b,kind=8)
  
  failed = .false.
  
  neq = size(b)
  
  
  if (neq<= 1) stop "neq<= 1"
  if (neq /= 11) stop "neq /= 11"
  
  !check symmetry and negative
  do j=1,neq
      do i=1,neq
        if (int_a(i,j) /= int_a(j,i) .or. int_a(i,j) < 0) then
            print *,i,j,int_a(i,j),int_a(j,i)
            failed = .true.
            go to 100
        end if
      end do
  end do


 100 if (failed) then
    !debug krigin system
    do i=1,neq
      write(*,888) i,b(i),(int_a(i,j),j=1,neq)
    end do

888        format('SOLVER ERROR b(',i3,') =',f7.4,'  a= ',99f7.4)  

    stop "symmetry failed"
  end if

  info = 0
  
  call dposv_f95(int_a,y,info=info_tmp)
  print *,'dposv_f95::info=',info_tmp
  
  !call ssytrf_f95(int_a,ipiv=ipiv,info=info_tmp,uplo='U')
  if (info_tmp>0) then
    do i=1,neq
      write(*,888) i,b(i),(int_a(i,j),j=1,neq)
    end do

    if (present(info)) info = info_tmp
    return
  end if
  !call ssytrs_f95(int_a,y,ipiv=ipiv,info=info_tmp,uplo='U')
  !if (info_tmp>0) then
  !  if (present(info)) info = info_tmp
  !  return
  !end if
  x = real(y(:,1),kind=fk)

end subroutine

end module
