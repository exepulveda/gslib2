module sorting
use kinds_kit
implicit none

contains

subroutine sortem_real_int(ib,ie,a,b)
implicit none
!arguments
integer(kind=gik), intent(in) :: ib,ie
real(kind=fk), intent(inout) :: a(:)
integer(kind=gik), intent(inout), dimension(:) :: b
!locals 
integer(kind=gik) iring
integer(kind=gik)   lt(64),ut(64),i,j,k,m,p,q
real(kind=fk) ta,xa
integer(kind=gik) tb,xb


!
! Initialize:
!
      j     = ie
      m     = 1
      i     = ib
      iring = 2
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
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
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
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
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
 64      b(i) = b(q)
         b(q) = tb
 65   continue
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
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
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


!original sortem from gslib
subroutine sortem_original(ib,ie,a,b,c,d,e,f,g,h)
implicit none
!arguments
integer(kind=gik), intent(in) :: ib,ie
real(kind=fk), intent(inout) :: a(:)
real(kind=fk), intent(inout), dimension(:), optional :: b,c,d,e,f,g,h
!locals 
integer(kind=gik) iperm,iring
integer(kind=gik)   lt(64),ut(64),i,j,k,m,p,q
real(kind=fk) ta,tb,tc,td,te,tf,tg,th,xa,xb,xc,xd,xe,xf,xg,xh

iperm = 0
if (present(b)) iperm = 1
if (present(c)) iperm = 2
if (present(d)) iperm = 3
if (present(e)) iperm = 4
if (present(f)) iperm = 5
if (present(g)) iperm = 6
if (present(h)) iperm = 7

!
! Initialize:
!
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
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
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
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
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
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
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
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
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
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
