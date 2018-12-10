module solvers
use kinds_kit
implicit none
contains
#ifdef USE_MKL
subroutine system_solver(a,b,x,option,info)
use lapack95
implicit none
  real(kind=fk), intent(in) :: a(:,:),b(:)
  real(kind=fk), intent(out) :: x(:)
  integer, intent(in), optional :: option
  integer, intent(out), optional :: info
  !locals
  real(kind=8) int_a(size(b),size(b))
  integer ipiv(size(b)),i,j,neq,info_tmp
  real(kind=8) y(size(b),1)
  integer local_option !0=LU, 1=LDL', 2=Cholesky
  logical failed
  
  int_a = real(a,kind=8)
  y(:,1) = real(b,kind=8)
  
  local_option = 0
  if (present(option)) local_option = option
  if (local_option < 0 .or. local_option > 1 ) local_option = 0
  
  neq = size(b)

#ifdef DEBUG_SOLVER  
  failed = .false.
  !debug krigin system
  do i=1,neq
    write(*,999) i,b(i),(a(i,j),j=1,neq)
  end do

999        format('SOLVER INICIAL b(',i3,') =',f7.4,'  a= ',99f7.4)  
#endif

  
#ifdef DEBUG_SOLVER  
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
#endif
  info = 0

  if (option == 1) then
    call dsysv_f95(int_a,y,ipiv=ipiv,info=info_tmp)
  else if (option == 2) then
    call dposv_f95(int_a,y,info=info_tmp)
  else
    call dgesv_f95(int_a,y,ipiv=ipiv,info=info_tmp)
  end if  
#ifdef DEBUG_SOLVER  
  print *,'dgesv_f95::info=',info_tmp
  
  !call ssytrf_f95(int_a,ipiv=ipiv,info=info_tmp,uplo='U')
  if (info_tmp>0) then
    do i=1,neq
      write(*,888) i,b(i),(int_a(i,j),j=1,neq)
    end do

    if (present(info)) info = info_tmp
    return
  end if
#endif
  x = real(y(:,1),kind=fk)

end subroutine
#else
subroutine system_solver(a,b,x,option,info)
implicit none
  real(kind=fk), intent(in) :: a(:,:),b(:)
  real(kind=fk), intent(out) :: x(:)
  integer, intent(in), optional :: option
  integer, intent(out), optional :: info
  !locals
  real(kind=8) int_a(size(b),size(b)),int_b(size(b))
  integer INDX(size(b))
  real(kind=8) y(size(b))
  integer n,D,rc
  logical failed
  
  int_a = real(a,kind=8)
  int_b = real(b,kind=8)
  n = size(b)

  call gauss_2(int_a,int_b,y,n)
  
  x = real(y,kind=fk)
  
  info = 0
  
end subroutine

subroutine gauss_2(a,b,x,n)
!===========================================================
! Solutions to a system of linear equations A*x=b
! Method: Gauss elimination (with scaling and pivoting)
! Alex G. (November 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - array of the right hand coefficients b
! n      - number of equations (size of matrix A)
! output ...
! x(n)   - solutions
! coments ...
! the original arrays a(n,n) and b(n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), b(n), x(n)
double precision s(n)
double precision c, pivot, store
integer i, j, k, l

! step 1: begin forward elimination
do k=1, n-1

! step 2: "scaling"
! s(i) will have the largest element from row i 
  do i=k,n                       ! loop over rows
    s(i) = 0.0
    do j=k,n                    ! loop over elements of row i
      s(i) = max(s(i),abs(a(i,j)))
    end do
  end do

! step 3: "pivoting 1" 
! find a row with the largest pivoting element
  pivot = abs(a(k,k)/s(k))
  l = k
  do j=k+1,n
    if(abs(a(j,k)/s(j)) > pivot) then
      pivot = abs(a(j,k)/s(j))
      l = j
    end if
  end do

! Check if the system has a sigular matrix
  if(pivot == 0.0) then
    write(*,*) ' The matrix is sigular '
    return
  end if

! step 4: "pivoting 2" interchange rows k and l (if needed)
if (l /= k) then
  do j=k,n
     store = a(k,j)
     a(k,j) = a(l,j)
     a(l,j) = store
  end do
  store = b(k)
  b(k) = b(l)
  b(l) = store
end if

! step 5: the elimination (after scaling and pivoting)
   do i=k+1,n
      c=a(i,k)/a(k,k)
      a(i,k) = 0.0
      b(i)=b(i)- c*b(k)
      do j=k+1,n
         a(i,j) = a(i,j)-c*a(k,j)
      end do
   end do
end do

! step 6: back substiturion 
x(n) = b(n)/a(n,n)
do i=n-1,1,-1
   c=0.0
   do j=i+1,n
     c= c + a(i,j)*x(j)
   end do 
   x(i) = (b(i)- c)/a(i,i)
end do

end subroutine gauss_2
#endif

end module
