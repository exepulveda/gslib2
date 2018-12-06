module sorting
use kinds_kit
implicit none

contains

function swap(a,i,j)
  implicit none
  !arguments
  integer(kind=gik), intent(in) :: i,j
  real(kind=fk), intent(inout) :: a(:)
  real(kind=fk) :: swap
  !locals
  real(kind=fk) :: tmp

  tmp = a(j)
  a(j) = a(i)
  a(i) = tmp


  swap = tmp
end function

subroutine swap_put(a,m,n,ta)
  implicit none
  !arguments
  integer(kind=gik), intent(in) :: m,n
  real(kind=fk), intent(inout) :: a(:)
  real(kind=fk), intent(in) :: ta

  !locals
  real(kind=fk) :: tmp

  a(n) = a(m)
  a(m) = ta

end subroutine


function swap_set(a,m,n)
  implicit none
  !arguments
  integer(kind=gik), intent(in) :: m,n
  real(kind=fk), intent(inout) :: a(:)
  real(kind=fk) :: swap_set
  !locals
  real(kind=fk) :: tmp

  swap_set = a(m)
  a(m) = a(n)
end function

!original sortem from gslib
subroutine sortem_original(ib,ie,a,b,c,d,e,f,g,h)
implicit none
!arguments
integer(kind=gik), intent(in) :: ib,ie
real(kind=fk), intent(inout) :: a(:)
real(kind=fk), intent(inout), dimension(:), optional :: b,c,d,e,f,g,h
!
! The dimensions for lt and ut have to be at least log (base 2) n
!
integer   lt(64),ut(64),i,j,k,m,p,q
real(kind=fk) ta,tb,tc,td,te,tf,tg,th,xa,xb,xc,xd,xe,xf,xg,xh
!
! Initialize:
!
j     = ie
m     = 1
i     = ib
!
! If this segment has more than two elements  we split it
!
10   if (j-i-1) 100,90,15
!
! p is the position of an arbitrary element in the segment we choose the
! middle element. Under certain circumstances it may be advantageous
! to choose p at random.
!
15   p    = (j+i)/2

ta = swap_set(a,p,i)

if (present(b)) tb = swap_set(b,p,i)
if (present(c)) tc = swap_set(c,p,i)
if (present(d)) td = swap_set(d,p,i)
if (present(e)) te = swap_set(e,p,i)
if (present(f)) tf = swap_set(f,p,i)
if (present(g)) tg = swap_set(g,p,i)
if (present(h)) th = swap_set(h,p,i)

!
! Start at the beginning of the segment, search for k such that a(k)>t
!
q = j
k = i
20   k = k+1
if(k.gt.q)     go to 60
if(a(k).le.ta) go to 20
!
! Such an element has now been found now search for a q such that a(q)<t
! starting at the end of the segment.
!
30   continue
if(a(q).lt.ta) go to 40
q = q-1
if(q.gt.k)     go to 30
go to 50
!
! a(q) has now been found. we interchange a(q) and a(k)
!
!40   xa   = a(k)


40   xa = swap(a,q,k)

if (present(b)) xb = swap(b,q,k)
if (present(c)) xc = swap(c,q,k)
if (present(d)) xd = swap(d,q,k)
if (present(e)) xe = swap(e,q,k)
if (present(f)) xf = swap(f,q,k)
if (present(g)) xg = swap(g,q,k)
if (present(h)) xh = swap(h,q,k)

!
! Update q and search for another pair to interchange:
!
q = q-1
go to 20
50   q = k-1
60   continue
!
! The upwards search has now met the downwards search:
!
call swap_put(a,i,q,ta)

if (present(b)) call swap_put(b,q,i,tb)
if (present(c)) call swap_put(c,i,q,tc)
if (present(d)) call swap_put(d,i,q,td)
if (present(e)) call swap_put(e,i,q,te)
if (present(f)) call swap_put(f,i,q,tf)
if (present(g)) call swap_put(g,i,q,tg)
if (present(h)) call swap_put(h,i,q,th)

!
! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
! store the position of the largest segment in lt and ut
!
if (2*q.le.i+j) go to 70
lt(m) = i
ut(m) = q-1
i = q+1
go to 80
70   lt(m) = q+1
ut(m) = j
j = q-1
!
! Update m and split the new smaller segment
!
80   m = m+1
go to 10
!
! We arrive here if the segment has  two elements we test to see if
! the segment is properly ordered if not, we perform an interchange
!
90   continue
if (a(i).le.a(j)) go to 100

xa = swap(a,i,j)

if (present(b)) xa = swap(b,i,j)
if (present(c)) xa = swap(c,i,j)
if (present(d)) xa = swap(d,i,j)
if (present(e)) xa = swap(e,i,j)
if (present(f)) xa = swap(f,i,j)
if (present(g)) xa = swap(g,i,j)
if (present(h)) xa = swap(h,i,j)


!
! If lt and ut contain more segments to be sorted repeat process:
!
100  m = m-1
if (m.le.0) go to 110
i = lt(m)
j = ut(m)
go to 10
110  continue
return
end subroutine

subroutine sortem(ib,ie,a,sec)
!                      Quickersort Subroutine
!                      **********************
implicit none
!arguments
integer(kind=gik), intent(in) :: ib,ie
real(kind=fk), intent(inout) :: a(:)
real(kind=fk), intent(inout), optional :: sec(:,:)
!locals
integer   lt(64),ut(64),i,j,k,m,p,q
real(kind=fk),allocatable :: tmp_sec(:)
real(kind=fk) :: ta,xa

if (present(sec)) then
  allocate(tmp_sec(size(sec,2)))
end if
!
! The dimensions for lt and ut have to be at least log (base 2) n
!
!
! Initialize:
!
      j     = ie
      m     = 1
      i     = ib
!
! If this segment has more than two elements  we split it
!
 10   if (j-i-1) 100,90,15
!
! p is the position of an arbitrary element in the segment we choose the
! middle element. Under certain circumstances it may be advantageous
! to choose p at random.
!
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      !secondary variables swap
      if (present(sec)) then
        tmp_sec   = sec(k,:)
        sec(k,:) = sec(q,:)
        sec(q,:) = tmp_sec
      end if
!
! Start at the beginning of the segment, search for k such that a(k)>t
!
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
!
! Such an element has now been found now search for a q such that a(q)<t
! starting at the end of the segment.
!
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
!
! a(q) has now been found. we interchange a(q) and a(k)
!
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      if (present(sec)) then
        tmp_sec   = sec(k,:)
        sec(k,:) = sec(q,:)
        sec(q,:) = tmp_sec
      end if

!
! Update q and search for another pair to interchange:
!
      q = q-1
      go to 20
 50   q = k-1
 60   continue
!
! The upwards search has now met the downwards search:
!
      a(i)=a(q)
      a(q)=ta

      if (present(sec)) then
        tmp_sec   = sec(k,:)
        sec(k,:) = sec(q,:)
        sec(q,:) = tmp_sec
      end if
!
! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
! store the position of the largest segment in lt and ut
!
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
!
! Update m and split the new smaller segment
!
 80   m = m+1
      go to 10
!
! We arrive here if the segment has  two elements we test to see if
! the segment is properly ordered if not, we perform an interchange
!
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      if (present(sec)) then
        tmp_sec   = sec(k,:)
        sec(k,:) = sec(q,:)
        sec(q,:) = tmp_sec
      end if

!
! If lt and ut contain more segments to be sorted repeat process:
!
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
end subroutine

end module
