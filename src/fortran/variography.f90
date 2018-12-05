module variography
use kinds_kit

real(8), parameter :: PI=3.14159265
integer, parameter :: PMX=999.
real(4), parameter :: EPSLON=1.e-5

type vmodel_structure_type
  integer(kind=gik) type
  real(kind=fk) sill
  real(kind=fk) range(3)
  real(kind=fk) angle(3)
  real(kind=fk) rotmat(3,3)
end type

type vmodel_type
  integer(kind=gik) :: nst
  real(kind=fk) :: nugget
  type(vmodel_structure_type),allocatable,dimension(:) :: structure
end type

integer(4), parameter :: SPHERICAL = 1
integer(4), parameter :: EXPONENTIAL = 2
integer(4), parameter :: GAUSSIAN = 3
integer(4), parameter :: POWER = 4
integer(4), parameter :: HOLE_EFFECT = 5

public setup_vmodel_structure,cova3

contains

subroutine setup_vmodel_structure(vmodel,is,type,sill,r1,r2,r3,a1,a2,a3)
use geometry, only: setrot
!arguments
type(vmodel_type), intent(inout) :: vmodel
integer(kind=gik), intent(in) :: is
integer(kind=gik), intent(in) :: type
real(kind=fk), intent(in) :: sill
real(kind=fk), intent(in) :: r1,r2,r3
real(kind=fk), intent(in) :: a1,a2,a3
!locals
  vmodel%structure(is)%type = type
  vmodel%structure(is)%sill = sill
  vmodel%structure(is)%angle(1) = a1
  vmodel%structure(is)%angle(2) = a2
  vmodel%structure(is)%angle(3) = a3

  vmodel%structure(is)%range(1) = r1
  vmodel%structure(is)%range(2) = r2
  vmodel%structure(is)%range(3) = r3

  call setrot(vmodel%structure(is)%angle(1), &
              vmodel%structure(is)%angle(2), &
              vmodel%structure(is)%angle(3), &
              vmodel%structure(is)%range(1)/vmodel%structure(is)%range(2), &
              vmodel%structure(is)%range(1)/vmodel%structure(is)%range(3), &
              vmodel%structure(is)%rotmat)
end subroutine

pure function max_cova(p1,p2,vmodel) result(cmax)
!dummies
real(kind=fk), intent(in),dimension(:) :: p1,p2
type(vmodel_type), intent(in) :: vmodel
real(kind=fk) :: cmax
!locals
integer(4) i
real(kind=dk) hsqd,h,aa,cc
!
! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):
!
  cmax   = vmodel%nugget
  do i=1,vmodel%nst
    if (vmodel%structure(i)%type  == POWER_MODEL) then
      cmax = cmax + PMX
    else
      cmax = cmax + vmodel%structure(i)%sill
    endif
  end do
end function

pure function cova3_new(p1,p2,vmodel) result(cova)
use geometry, only : sqdist
!dummies
real(kind=fk), intent(in),dimension(:) :: p1,p2
type(vmodel_type), intent(in) :: vmodel
real(kind=fk) :: cova
!locals
real(kind=fk) :: cmax
integer(4) i
real(kind=dk) hsqd,h,aa,cc
logical :: cmax_calculated

!
! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):
!
  cmax_calculated = .false.
!
! Check for "zero" distance, return with cmax if so:
!
  hsqd = sqdist(p1(1),p1(2),p1(3),p2(1),p2(2),p2(3),vmodel%structure(1)%rotmat)
  if (hsqd <= EPSILON) then
    cova = max_cova(p1,p2,vmodel)
    return
  endif
!
! Loop over all the structures:
!
  cova = 0.0
  do i=1,vmodel%nst
!
! Compute the appropriate distance:
!
    if (i /= 1) then
      hsqd=sqdist(p1(1),p1(2),p1(3),p2(1),p2(2),p2(3),vmodel%structure(i)%rotmat)
    end if
    h = sqrt(hsqd)

    aa = vmodel%structure(i)%range(1)
    cc = vmodel%structure(i)%sill

    select case(vmodel%structure(i)%type)
    case (SPHERICAL)
      hr = h/aa
      if(hr.lt.1.) cova=cova+cc*(1.-hr*(1.5-.5*hr*hr))
    case (EXPONENTIAL)
      cova = cova + cc*exp(-3.0*h/aa)
    case (GAUSSIAN)
      cova = cova + cc*exp(-3.*(h/aa)*(h/aa))
    case (POWER)
      if (.not. cmax_calculated) then
        cmax = max_cova(p1,p2,vmodel)
        cmax_calculated = .true.
      endif

      cova = cova + cmax - cc*(h**aa)
    case (HOLE_EFFECT)
      cova = cova + cc*cos(h/aa*PI)
    end select
  end do
end function

subroutine cova3(x1,y1,z1,x2,y2,z2,nst,c0,it,cc,aa,vrotmat,cmax,cova)
use geometry, only : sqdist
!dummies
real(kind=fk), intent(in) :: x1,y1,z1,x2,y2,z2
real(kind=fk), intent(in) :: c0,cc(:),aa(:),vrotmat(:,:,:)
integer(kind=gik), intent(in) :: it(:),nst
real(kind=fk), intent(out) :: cmax,cova
!locals
integer(4) i,is,ist
real(kind=dk) hsqd,h

!
! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):
!
cmax   = c0
do is=1,nst
      ist = is
      if(it(ist).eq.4) then
            cmax = cmax + PMX
      else
            cmax = cmax + cc(ist)
      endif
end do
!
! Check for "zero" distance, return with cmax if so:
!
hsqd = sqdist(x1,y1,z1,x2,y2,z2,vrotmat(1,:,:))
if(real(hsqd).lt.EPSLON) then
      cova = cmax
      return
endif
!
! Loop over all the structures:
!
cova = 0.0
do is=1,nst
      ist = is
!
! Compute the appropriate distance:
!
      if(ist.ne.1) then
            ir = min((irot+is-1),MAXROT)
            hsqd=sqdist(x1,y1,z1,x2,y2,z2,vrotmat(:,:,ist))
      end if
      h = real(dsqrt(hsqd))
!
! Spherical Variogram Model?
!
      if(it(ist).eq.SPHERICAL) then
            hr = h/aa(ist)
            if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
!
! Exponential Variogram Model?
!
else if(it(ist).eq.EXPONENTIAL) then
            cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
!
! Gaussian Variogram Model?
!
else if(it(ist).eq.GAUSSIAN) then
            cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
!
! Power Variogram Model?
!
else if(it(ist).eq.POWER_MODEL) then
            cova = cova + cmax - cc(ist)*(h**aa(ist))
!
! Hole Effect Model?
!
else if(it(ist).eq.HOLE_EFFECT) then
!                 d = 10.0 * aa(ist)
!                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
            cova = cova + cc(ist)*cos(h/aa(ist)*PI)
      endif
end do
!
! Finished:
!
return
end subroutine

end module variography
