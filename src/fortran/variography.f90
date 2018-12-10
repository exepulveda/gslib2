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

integer(4), parameter :: SPHERICAL_MODEL = 1
integer(4), parameter :: EXPONENTIAL_MODEL = 2
integer(4), parameter :: GAUSSIAN_MODEL = 3
integer(4), parameter :: POWER_MODEL = 4
integer(4), parameter :: HOLE_EFFECT_MODEL = 5

public setup_vmodel_structure,cova3

contains

subroutine setup_vmodel(vmodel,nst,c0,it,cc,aa,ang1,ang2,ang3,anis1,anis2)
  use geometry, only: setrot
  implicit none
  !arguments
  type(vmodel_type), intent(inout) :: vmodel
  integer(kind=gik), intent(in) :: nst
  integer(kind=gik), intent(in) :: it(:)
  real(kind=fk), intent(in) :: c0,cc(:),aa(:)
  real(kind=fk), intent(in) :: ang1(:),ang2(:),ang3(:)
  real(kind=fk), intent(in) :: anis1(:),anis2(:)
  !locals 
  integer is

  vmodel%nugget = c0
  vmodel%nst = nst

  allocate(vmodel%structure(nst))

  do is=1,nst
    vmodel%structure(is)%type = it(is)
    vmodel%structure(is)%sill = cc(is)
    vmodel%structure(is)%angle(1) = ang1(is)
    vmodel%structure(is)%angle(2) = ang2(is)
    vmodel%structure(is)%angle(3) = ang3(is)

    vmodel%structure(is)%range(1) = aa(is)
    !vmodel%structure(is)%range(2) = range2(is)
    !vmodel%structure(is)%range(3) = range3(is)

    call setrot(vmodel%structure(is)%angle(1), &
                vmodel%structure(is)%angle(2), &
                vmodel%structure(is)%angle(3), &
                anis1(is), &
                anis2(is), &
                vmodel%structure(is)%rotmat)
                
  end do
end subroutine

subroutine print_vmodel(vmodel)
implicit none
  !arguments
  type(vmodel_type), intent(in) :: vmodel
  !locals 
  integer is

  print *, "nugget=",vmodel%nugget
  print *, "nst=",vmodel%nst

  do is=1,vmodel%nst
    print *, "st=",1,"type=",vmodel%structure(is)%type
    print *, "st=",1,"sill=",vmodel%structure(is)%sill
    print *, "st=",1,"range=",vmodel%structure(is)%range(1)
    print *, "st=",1,"rotmat=",vmodel%structure(is)%rotmat
  end do
end subroutine


pure function max_cova(vmodel) result(cmax)
implicit none
!dummies
type(vmodel_type), intent(in) :: vmodel
!return
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

subroutine cova3_vmodel(x1,y1,z1,x2,y2,z2,vmodel,cmax,cova)
use geometry, only : sqdist
implicit none
!dummies
real(kind=fk), intent(in) :: x1,y1,z1,x2,y2,z2
type(vmodel_type), intent(in) :: vmodel
real(kind=fk),intent(out) :: cova,cmax
!locals
integer(4) i
real(kind=dk) hsqd,h,aa,cc,hr
logical :: cmax_calculated

!
! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):
!
  cmax_calculated = .false.
!
! Check for "zero" distance, return with cmax if so:
!
  hsqd = sqdist(x1,y1,z1,x2,y2,z2,vmodel%structure(1)%rotmat)
  if (hsqd <= EPSLON) then
    cova = max_cova(vmodel)
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
      hsqd=sqdist(x1,y1,z1,x2,y2,z2,vmodel%structure(i)%rotmat)
    end if
    h = sqrt(hsqd)

    aa = vmodel%structure(i)%range(1)
    cc = vmodel%structure(i)%sill

    select case(vmodel%structure(i)%type)
    case (SPHERICAL_MODEL)
      hr = h/aa
      if(hr.lt.1.) cova=cova+cc*(1.-hr*(1.5-.5*hr*hr))
    case (EXPONENTIAL_MODEL)
      cova = cova + cc*exp(-3.0*h/aa)
    case (GAUSSIAN_MODEL)
      cova = cova + cc*exp(-3.*(h/aa)*(h/aa))
    case (POWER_MODEL)
      if (.not. cmax_calculated) then
        cmax = max_cova(vmodel)
        cmax_calculated = .true.
      endif

      cova = cova + cmax - cc*(h**aa)
    case (HOLE_EFFECT_MODEL)
      cova = cova + cc*cos(h/aa*PI)
    end select
  end do
end subroutine

subroutine cova3(x1,y1,z1,x2,y2,z2,nst,c0,it,cc,aa,vrotmat,cmax,cova)
  use geometry, only : sqdist
  implicit none
!dummies
real(kind=fk), intent(in) :: x1,y1,z1,x2,y2,z2
real(kind=fk), intent(in) :: c0,cc(:),aa(:),vrotmat(:,:,:)
integer(kind=gik), intent(in) :: it(:),nst
real(kind=fk), intent(out) :: cmax,cova
!locals
integer(4) i,is,ist
real(kind=dk) hsqd,h,hr

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
      if (cova <0) stop "COVA < 0"
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
            hsqd=sqdist(x1,y1,z1,x2,y2,z2,vrotmat(ist,:,:))
      end if
      h = real(dsqrt(hsqd))
!
! Spherical Variogram Model?
!
      if(it(ist).eq.SPHERICAL_MODEL) then
            hr = h/aa(ist)
            if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
!
! Exponential Variogram Model?
!
else if(it(ist).eq.EXPONENTIAL_MODEL) then
            cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
!
! Gaussian Variogram Model?
!
else if(it(ist).eq.GAUSSIAN_MODEL) then
            cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
!
! Power Variogram Model?
!
else if(it(ist).eq.POWER_MODEL) then
            cova = cova + cmax - cc(ist)*(h**aa(ist))
!
! Hole Effect Model?
!
else if(it(ist).eq.HOLE_EFFECT_MODEL) then
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
