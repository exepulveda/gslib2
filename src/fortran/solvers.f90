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
  real(kind=8) int_a(size(b),size(b))
  integer INDX(size(b))
  real(kind=8) y(size(b))
  integer n,D,rc
  logical failed
  
  int_a = real(a,kind=8)
  n = size(b)

  
  call LUDCMP(int_a,n,INDX,D,rc)
  
  if (rc.eq.0) then
    call LUBKSB(int_a,n,INDX,y)
    x = real(y,kind=fk)
  endif
  info=int(rc,kind=gik)

end subroutine


    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the LU *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. INDX is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
     Subroutine LUDCMP(A,N,INDX,D,CODE)
     IMPLICIT NONE
     integer, parameter :: nmax = 100
     real, parameter :: tiny = 1.5D-16

     real*8, intent(inout), dimension(N,N) :: A
     integer, intent(in) :: N
     integer, intent(out) :: D, CODE
     integer, intent(out), dimension(N) :: INDX
     !f2py depend(N) A, indx

     REAL*8  :: AMAX, DUM, SUMM, VV(NMAX)
     INTEGER :: i, j, k, imax

     D=1; CODE=0

     DO I=1,N
       AMAX=0.d0
       DO J=1,N
         IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
       END DO ! j loop
       IF(AMAX.LT.TINY) THEN
         CODE = 1
         RETURN
       END IF
       VV(I) = 1.d0 / AMAX
     END DO ! i loop

     DO J=1,N
       DO I=1,J-1
         SUMM = A(I,J)
         DO K=1,I-1
           SUMM = SUMM - A(I,K)*A(K,J) 
         END DO ! k loop
         A(I,J) = SUMM
       END DO ! i loop
       AMAX = 0.d0
       DO I=J,N
         SUMM = A(I,J)
         DO K=1,J-1
           SUMM = SUMM - A(I,K)*A(K,J) 
         END DO ! k loop
         A(I,J) = SUMM
         DUM = VV(I)*DABS(SUMM)
         IF(DUM.GE.AMAX) THEN
           IMAX = I
           AMAX = DUM
         END IF
       END DO ! i loop  
       
       IF(J.NE.IMAX) THEN
         DO K=1,N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
         END DO ! k loop
         D = -D
         VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

       IF(J.NE.N) THEN
         DUM = 1.d0 / A(J,J)
         DO I=J+1,N
           A(I,J) = A(I,J)*DUM
         END DO ! i loop
       END IF 
     END DO ! j loop

     RETURN
 END subroutine LUDCMP
 
!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************

!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 subroutine LUBKSB(A, N, INDX, B)
 implicit none
 integer, intent(in) :: N 
 real*8, intent(in), dimension(N,N) :: A
 integer, intent(in), dimension(N) :: INDX
 real*8, intent(inout), dimension(N) :: B
 !locals
 integer ii,i,j,ll

 REAL*8  SUMM

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUMM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUMM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUMM
 END DO ! i loop

 DO I=N,1,-1
   SUMM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUMM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine LUBKSB
#endif

end module
