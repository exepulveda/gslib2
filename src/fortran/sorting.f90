module sorting
use kinds_kit
implicit none

contains

!original sortem from gslib
subroutine sortem(ib,ie,a,sec)
!                      Quickersort Subroutine
!                      **********************
implicit none
!arguments
integer(kind=ik), intent(in) :: ib,ie
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
