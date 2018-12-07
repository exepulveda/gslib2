module interpolation
use kinds_kit
implicit none

real(kind=fk), parameter :: EPSLON=1.0e-20


contains

!-----------------------------------------------------------------------
!
! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
!                 for a value of x and a power pow.
!
!-----------------------------------------------------------------------
elemental function powint(xlow,xhigh,ylow,yhigh,xval,pow)
  implicit none
  real(kind=fk), intent(in) :: xlow
  real(kind=fk), intent(in) :: xhigh
  real(kind=fk), intent(in) :: ylow
  real(kind=fk), intent(in) :: yhigh
  real(kind=fk), intent(in) :: xval
  real(kind=fk), intent(in) :: pow
  !return
  real(kind=fk) :: powint

      if((xhigh-xlow).lt.EPSLON) then
            powint = (yhigh+ylow)/2.0
      else
            powint = ylow + (yhigh-ylow)* & 
                    (((xval-xlow)/(xhigh-xlow))**pow)
      end if

      return
end function

!-----------------------------------------------------------------------
!
! Given an array "xx" of length "n", and given a value "x", this routine
! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
! returned to indicate that x is out of range.
!
! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
!-----------------------------------------------------------------------
pure function locate(xx,is,ie,x) result(j)
  implicit none
  real(kind=fk), intent(in) :: xx(:)
  integer(kind=gik), intent(in) :: is
  integer(kind=gik), intent(in) :: ie
  real(kind=fk), intent(in) :: x
  integer(kind=gik) :: j
  !locals
  integer(kind=gik) :: n,jl,ju,jm
  
  n = size(xx)
!
! Initialize lower and upper methods:
!
      if(is.le.0) then
        jl = 0
      else
       jl = is-1
      end if
      
      ju = ie
      if(xx(n).le.x) then
            j = ie
            return
      end if
!
! If we are not done then compute a midpoint:
!
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
!
! Replace the lower or upper limit with the midpoint:
!
            if((xx(ie).gt.xx(is)).eqv.(x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
!
! Return with the array index:
!
      j = jl
      return
end function

      
end module
