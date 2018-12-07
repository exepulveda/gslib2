module gaussian_utils
use kinds_kit
implicit none

real(kind=fk), parameter :: EPSLON=1.0e-20


contains

subroutine nscore(vr,tmin,tmax,iwt,wt,tmp,vrg,ierror)
!-----------------------------------------------------------------------
!
!              Transform Univariate Data to Normal Scores
!              ******************************************
!
! This subroutibe takes "nd" data "vr(i),i=1,...,nd" possibly weighted
! by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
! as "vrg(i),i=1,...,nd".  The extra storage array "tmp" is required
! so that the data can be returned in the same order (just in case
! there are associated arrays like the coordinate location).
!
!
!
! INPUT VARIABLES:
!
!   nd               Number of data (no missing values)
!   vr(nd)           Data values to be transformed
!   tmin,tmax        data trimming limits
!   iwt              =0, equal weighted; =1, then apply weight
!   wt(nd)           Weight for each data (don't have to sum to 1.0)
!   tmp(nd)          Temporary storage space for sorting
!   lout             if > 0 then transformation table will be written
!
!
!
! OUTPUT VARIABLES:
!
!   vrg(nd)          normal scores
!   ierror           error flag (0=error free,1=problem)
!
!
!
! EXTERNAL REFERENCES:
!
!   gauinv           Calculates the inverse of a Gaussian cdf
!   sortem           sorts a number of arrays according to a key array
!
!
!
!-----------------------------------------------------------------------
use sorting, only: sortem_original
implicit none
!arguments
real(kind=fk), intent(inout)    :: vr(:)
real(kind=fk), intent(in)    :: tmin,tmax
integer(kind=gik), intent(in) :: iwt
real(kind=fk), intent(inout)    :: wt(:)
real(kind=fk), intent(out)   :: vrg(:)
integer(kind=gik), intent(out) :: ierror
!locals
real(kind=fk) tmp(size(vr)),oldcp,cp
real*8    inv_val,twt
integer i,nd
integer(kind=gik) ierr

nd = size(vr)
!
! Sort the data in ascending order and calculate total weight:
!
      ierror = 0
      twt    = 0.0
      do i=1,nd
            tmp(i) = real(i)
            if(vr(i).ge.tmin.and.vr(i).lt.tmax) then
                  if(iwt.eq.0) then
                        twt = twt + 1.
                  else
                        twt = twt + wt(i)
                  end if
            end if
      end do
      if(nd.lt.1.or.twt.lt.EPSLON) then
            ierror = 1
            return
      end if
      call sortem_original(1,nd,vr,wt,tmp)
!
! Compute the cumulative probabilities:
!
      oldcp = 0.0
      cp    = 0.0
      do i=1,nd
            cp     =  cp + wt(i) / twt
            wt(i)  = (cp + oldcp)/ 2.0
            oldcp  =  cp
            call gauinv(dble(wt(i)),inv_val,ierr)
            vrg(i) = real(inv_val,kind=fk)
!            if(lout.gt.0) write(lout,'(f12.5,1x,f12.5)') vr(i),vrg(i)
      end do
!
! Get the arrays back in original order:
!
      call sortem_original(1,nd,tmp,wt,vr,vrg)
!
! Finished:
!
      return
end subroutine

subroutine gauinv(p,xp,ierr)
real(kind=dk),intent(in) :: p
real(kind=dk),intent(out) :: xp
integer(kind=gik),intent(out) ::  ierr
!locals
real*8 p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim
save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
!
! Coefficients of approximation:
!
data lim/1.0e-10/
data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/, &
     p3/-0.0204231210245/,p4/-0.0000453642210148/
data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/, &
     q3/0.103537752850/,q4/0.0038560700634/
!
! Check for an error situation:
!
ierr = 1
if(p.lt.lim) then
      xp = -1.0e10
      return
end if
if(p.gt.(1.0-lim)) then
      xp =  1.0e10
      return
end if
ierr = 0
!
! Get k for an error situation:
!
pp   = p
if(p.gt.0.5) pp = 1 - pp
xp   = 0.0
if(p.eq.0.5) return
!
! Approximate the function:
!
y  = dsqrt(dlog(1.0/(pp*pp)))
xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / &
               ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
if(real(p).eq.real(pp)) xp = -xp
!
! Return with G^-1(p):
!
return
end subroutine

!-----------------------------------------------------------------------
!
!           Back Transform Univariate Data from Normal Scores
!           *************************************************
!
! This subroutine backtransforms a standard normal deviate from a
! specified back transform table and option for the tails of the
! distribution.  Call once with "first" set to true then set to false
! unless one of the options for the tail changes.
!
!
!
! INPUT VARIABLES:
!
!   vrgs             normal score value to be back transformed
!   vr(:)           original data values that were transformed
!   vrg(:)          the corresponding transformed values
!   zmin,zmax        limits possibly used for linear or power model
!   ltail            option to handle values less than vrg(1):
!   ltpar            parameter required for option ltail
!   utail            option to handle values greater than vrg(nt):
!   utpar            parameter required for option utail
!
!
!
!-----------------------------------------------------------------------
real(kind=fk) function backtr(vrgs,vr,vrg,zmin,zmax,ltail,ltpar,utail,utpar)
  use interpolation, only: powint,locate
  implicit none
  !arguments
  real(kind=fk), intent(in)    :: vrgs
  real(kind=fk), intent(in)    :: vr(:)
  real(kind=fk), intent(in)    :: vrg(:)
  real(kind=fk), intent(in)    :: zmin,zmax
  integer(kind=gik), intent(in) :: ltail,utail
  real(kind=fk), intent(in)    :: ltpar,utpar
  !locals
  integer(kind=gik) nt,j
  real(kind=fk) lambda,cdflo,cdfbt,cdfhi,cpow
!
! Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
!
      if(vrgs.le.vrg(1)) then
            backtr = vr(1)
            cdflo  = gcum(vrg(1))
            cdfbt  = gcum(vrgs)
            if(ltail.eq.1) then
                  backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,1.0)
            else if(ltail.eq.2) then
                  cpow   = 1.0 / ltpar
                  backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,cpow)
            endif
!
! Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
!
      else if(vrgs.ge.vrg(nt)) then
            backtr = vr(nt)
            cdfhi  = gcum(vrg(nt))
            cdfbt  = gcum(vrgs)
            if(utail.eq.1) then
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,1.0)
            else if(utail.eq.2) then
                  cpow   = 1.0 / utpar
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,cpow)
            else if(utail.eq.4) then
                  lambda = (vr(nt)**utpar)*(1.0-gcum(vrg(nt)))
                  backtr = (lambda/(1.0-gcum(vrgs)))**(1.0/utpar)
            endif
      else
!
! Value within the transformation table:
!
            j = locate(vrg,1,nt,vrgs)
            j = max(min((nt-1),j),1)
            backtr = powint(vrg(j),vrg(j+1),vr(j),vr(j+1),vrgs,1.0)
      endif
      return
end function

!-----------------------------------------------------------------------
!
! Evaluate the standard normal cdf given a normal deviate x.  gcum is
! the area under a unit normal curve to the left of x.  The results are
! accurate only to about 5 decimal places.
!
!
!-----------------------------------------------------------------------
real(kind=fk) function gcum(x)
real(kind=fk), intent(in) :: x
real(kind=fk) z,t,e2
      z = x
      if(z.lt.0.) z = -z
      t    = 1./(1.+ 0.2316419*z)
      gcum = t*(0.31938153   + t*(-0.356563782 + t*(1.781477937 + &
             t*(-1.821255978 + t*1.330274429))))
      e2   = 0.
!
!  6 standard deviations out gets treated as infinity:
!
      if(z.le.6.) e2 = exp(-z*z/2.)*0.3989422803
      gcum = 1.0- e2 * gcum
      if(x.ge.0.) return
      gcum = 1.0 - gcum
      return
end function
end module
