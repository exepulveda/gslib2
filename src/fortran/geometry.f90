module geometry
use kinds_kit
implicit none

real(8), parameter :: DEG2RAD=3.141592654/180.0
real(8), parameter :: EPSLON=1.e-20

real(fk), parameter :: zero_3d(3) = (/ 0.0,0.0,0.0 /)
real(fk), parameter :: zero_2d(2) = (/ 0.0,0.0 /)

type grid_type
  integer(kind=ik) :: nodes(3)
  real(kind=fk) :: sizes(3)
  real(kind=fk) :: starts(3)
end type

contains

pure function sqdist(x1,y1,z1,x2,y2,z2,rotmat)
implicit none
real(kind=fk), intent(in) :: x1,y1,z1
real(kind=fk), intent(in) :: x2,y2,z2
real(kind=fk), intent(in) :: rotmat(3,3)

real(kind=dk) :: sqdist
!-----------------------------------------------------------------------
!
!    Squared Anisotropic Distance Calculation Given Matrix Indicator
!    ***************************************************************
!
! This routine calculates the anisotropic distance between two points
!  given the coordinates of each point and a definition of the
!  anisotropy.
!
!
! INPUT VARIABLES:
!
!   x1,y1,z1         Coordinates of first point
!   x2,y2,z2         Coordinates of second point
!   ind              The rotation matrix to use
!   MAXROT           The maximum number of rotation matrices dimensioned
!   rotmat           The rotation matrices
!
!
!
! OUTPUT VARIABLES:
!
!   sqdist           The squared distance accounting for the anisotropy
!                      and the rotation of coordinates (if any).
!
!
! NO EXTERNAL REFERENCES
!
!
!-----------------------------------------------------------------------

real(8) cont,dx,dy,dz
integer i
!
! Compute component distance vectors and the squared distance:
!
  dx = dble(x1 - x2)
  dy = dble(y1 - y2)
  dz = dble(z1 - z2)
  sqdist = 0.0
  do i=1,3
        cont   = rotmat(i,1) * dx &
               + rotmat(i,2) * dy &
               + rotmat(i,3) * dz
        sqdist = sqdist + cont * cont
  end do
  return
end function

subroutine setrot(ang1,ang2,ang3,anis1,anis2,rotmat)
implicit none

real(kind=fk), intent(in) :: ang1,ang2,ang3
real(kind=fk), intent(in) :: anis1,anis2
real(kind=fk), intent(inout) :: rotmat(3,3)
!-----------------------------------------------------------------------
!
!              Sets up an Anisotropic Rotation Matrix
!              **************************************
!
! Sets up the matrix to transform cartesian coordinates to coordinates
! accounting for angles and anisotropy (see manual for a detailed
! definition):
!
!
! INPUT PARAMETERS:
!
!   ang1             Azimuth angle for principal direction
!   ang2             Dip angle for principal direction
!   ang3             Third rotation angle
!   anis1            First anisotropy ratio
!   anis2            Second anisotropy ratio
!   ind              matrix indicator to initialize
!   MAXROT           maximum number of rotation matrices dimensioned
!   rotmat           rotation matrices
!
!
! NO EXTERNAL REFERENCES
!
!
!-----------------------------------------------------------------------

real(8) afac1,afac2,sina,sinb,sint,cosa,cosb,cost
real(8) alpha,beta,theta
!
! Converts the input angles to three angles which make more
!  mathematical sense:
!
!         alpha   angle between the major axis of anisotropy and the
!                 E-W axis. Note: Counter clockwise is positive.
!         beta    angle between major axis and the horizontal plane.
!                 (The dip of the ellipsoid measured positive down)
!         theta   Angle of rotation of minor axis about the major axis
!                 of the ellipsoid.
!
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
!
! Get the required sines and cosines:
!
      sina  = dble(sin(alpha))
      sinb  = dble(sin(beta))
      sint  = dble(sin(theta))
      cosa  = dble(cos(alpha))
      cosb  = dble(cos(beta))
      cost  = dble(cos(theta))
!
! Construct the rotation matrix in the required memory:
!
      afac1 = 1.0 / dble(max(anis1,EPSLON))
      afac2 = 1.0 / dble(max(anis2,EPSLON))
      rotmat(1,1) =       (cosb * cosa)
      rotmat(1,2) =       (cosb * sina)
      rotmat(1,3) =       (-sinb)
      rotmat(2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(2,3) = afac1*( sint * cosb)
      rotmat(3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(3,3) = afac2*(cost * cosb)
!
! Return to calling program:
!
return
end subroutine

subroutine getindx(n,min,siz,loc,index,inflag)
  integer(kind=ik),intent(in) ::   n
  integer(kind=ik),intent(out) ::   index
  real(kind=fk),intent(in) ::      min,siz,loc
  logical,intent(out) ::    inflag
  !
  ! Compute the index of "loc":
  !
  index = int( (loc-min)/siz + 1.5 )
  !
  ! Check to see if in or out:
  !
  if(index.lt.1) then
        index  = 1
        inflag = .false.
  else if(index.gt.n) then
        index  = n
        inflag = .false.
  else
        inflag = .true.
  end if
end subroutine

end module
