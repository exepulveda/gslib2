module searching
use kinds_kit
implicit none

type sb_structure_type
  integer(kind=gik) :: nxsup,nysup,nzsup
  real(kind=fk)    :: xsizsup,ysizsup,zsizsup
  real(kind=fk)    :: xmnsup,ymnsup,zmnsup
  integer(kind=gik), allocatable :: nisb(:)
  integer(kind=gik) , allocatable :: ixsbtosr(:),iysbtosr(:),izsbtosr(:)
  integer(kind=gik)  :: nsbtosr
end type

contains

!high level functions
function setup_superblock(x,y,z,vr,sec,grid,sb_grid,rotmat,radsqd,MAXSBX,MAXSBY,MAXSBZ)
use geometry, only: grid_type
implicit none
!arguments
  real(kind=fk), intent(in) :: x(:),y(:),z(:)
  real(kind=fk), intent(in) :: vr(:)
  real(kind=fk), intent(in),optional :: sec(:,:)
  type(grid_type), intent(in)  ::grid
  type(grid_type), intent(out) ::sb_grid
  real(kind=fk), intent(in) :: rotmat(:,:),radsqd
  integer(kind=gik), intent(in) :: MAXSBX,MAXSBy,MAXSBZ

!return
  type(sb_structure_type) :: setup_superblock
!locals
  type(sb_structure_type) :: sb_strcture
  integer MAXSB
  integer test

  MAXSB = min(grid%nodes(1),MAXSBX)*min(grid%nodes(2),MAXSBY)*min(grid%nodes(3),MAXSBZ)
  !determine
  allocate(sb_strcture%nisb(MAXSB),stat = test)

  call setsupr( &
      grid%nodes(1),grid%starts(1),grid%sizes(1), &
      grid%nodes(2),grid%starts(2),grid%sizes(2), &
      grid%nodes(3),grid%starts(3),grid%sizes(3), &
      x,y,z, &
      vr,sec,MAXSBX,MAXSBY,MAXSBZ,sb_strcture%nisb, &
      sb_grid%nodes(1),sb_grid%starts(1),sb_grid%sizes(1), &
      sb_grid%nodes(2),sb_grid%starts(2),sb_grid%sizes(2), &
      sb_grid%nodes(3),sb_grid%starts(3),sb_grid%sizes(3))

  allocate(sb_strcture%ixsbtosr(8 * MAXSB),stat = test)
  allocate(sb_strcture%ixsbtosr(8 * MAXSB),stat = test)
  allocate(sb_strcture%ixsbtosr(8 * MAXSB),stat = test)

  call picksup(sb_grid%nodes(1),sb_grid%sizes(1), &
               sb_grid%nodes(2),sb_grid%sizes(2), &
               sb_grid%nodes(3),sb_grid%sizes(3), &
               rotmat,radsqd,sb_strcture%nsbtosr,sb_strcture%ixsbtosr, &
               sb_strcture%iysbtosr,sb_strcture%izsbtosr)

  setup_superblock = sb_strcture
end function

subroutine search_super_block( &
    xloc,yloc,zloc,radsqd,rotmat, &
    sb_strcture, &
    ndmax,noct, &
    x,y,z, &
    nclose,close,infoct)
use geometry, only: grid_type
implicit none
!arguments
  real(kind=fk), intent(in) :: xloc,yloc,zloc
  real(kind=fk), intent(in) :: rotmat(:,:),radsqd
  type(sb_structure_type) :: sb_strcture
  integer(kind=gik), intent(in) :: ndmax,noct
  real(kind=fk), intent(in) :: x(:),y(:),z(:)
  integer(kind=gik), intent(out) :: nclose
  integer(kind=gik), intent(inout) :: close(:)
  integer(kind=gik), intent(inout) :: infoct(:)

!locals

  nclose = ndmax

  call srchsupr( &
      xloc,yloc,zloc,radsqd,rotmat, &
      sb_strcture%nsbtosr, &
      sb_strcture%ixsbtosr,sb_strcture%iysbtosr,sb_strcture%izsbtosr, &
      noct, &
      x,y,z,&
      sb_strcture%nisb, &
      sb_strcture%nxsup,sb_strcture%xmnsup,sb_strcture%xsizsup, &
      sb_strcture%nysup,sb_strcture%ymnsup,sb_strcture%ysizsup, &
      sb_strcture%nzsup,sb_strcture%zmnsup,sb_strcture%zsizsup, &
      nclose,close,infoct)

end subroutine
!low level
subroutine picksup( &
    nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup, &
    rotmat,radsqd,nsbtosr, &
    ixsbtosr,iysbtosr,izsbtosr)
!-----------------------------------------------------------------------
!
!             Establish Which Super Blocks to Search
!             **************************************
!
! This subroutine establishes which super blocks must be searched given
! that a point being estimated/simulated falls within a super block
! centered at 0,0,0.
!
!
!
! INPUT VARIABLES:
!
!   nxsup,xsizsup    Definition of the X super block grid
!   nysup,ysizsup    Definition of the Y super block grid
!   nzsup,zsizsup    Definition of the Z super block grid
!   rotmat           rotation matrices
!   radsqd           squared search radius
!
!
!
! OUTPUT VARIABLES:
!
!   nsbtosr          Number of super blocks to search
!   ixsbtosr         X offsets for super blocks to search
!   iysbtosr         Y offsets for super blocks to search
!   izsbtosr         Z offsets for super blocks to search
!
!
!
! EXTERNAL REFERENCES:
!
!   sqdist           Computes anisotropic squared distance
!
!
!
!-----------------------------------------------------------------------
use geometry, only: sqdist
implicit none
!arguments
integer(kind=gik), intent(in) :: nxsup,nysup,nzsup
real(kind=fk), intent(in) :: xsizsup,ysizsup,zsizsup
real(kind=fk), intent(in) :: rotmat(:,:),radsqd
integer(kind=gik), intent(inout) :: ixsbtosr(:),iysbtosr(:),izsbtosr(:)
integer(kind=gik), intent(inout) :: nsbtosr

!locals
integer i,j,k,i1,j1,k1,i2,j2,k2
real xo,yo,zo,xdis,ydis,zdis,hsqd,shortest
!
! MAIN Loop over all possible super blocks:
!
      nsbtosr = 0
      do i=-(nxsup-1),(nxsup-1)
      do j=-(nysup-1),(nysup-1)
      do k=-(nzsup-1),(nzsup-1)
            xo = real(i)*xsizsup
            yo = real(j)*ysizsup
            zo = real(k)*zsizsup
!
! Find the closest distance between the corners of the super blocks:
!
            shortest = 1.0e21
            do i1=-1,1
            do j1=-1,1
            do k1=-1,1
                  do i2=-1,1
                  do j2=-1,1
                  do k2=-1,1
                        if(i1.ne.0.and.j1.ne.0.and.k1.ne.0.and. &
                           i2.ne.0.and.j2.ne.0.and.k2.ne.0) then
                              xdis = real(i1-i2)*0.5*xsizsup + xo
                              ydis = real(j1-j2)*0.5*ysizsup + yo
                              zdis = real(k1-k2)*0.5*zsizsup + zo
                              hsqd = sqdist(0.0,0.0,0.0,xdis,ydis,zdis,rotmat)
                              if(hsqd.lt.shortest) shortest = hsqd
                        end if
                  end do
                  end do
                  end do
            end do
            end do
            end do
!
! Keep this super block if it is close enoutgh:
!
            if(real(shortest).le.radsqd) then
                  nsbtosr = nsbtosr + 1
                  ixsbtosr(nsbtosr) = i
                  iysbtosr(nsbtosr) = j
                  izsbtosr(nsbtosr) = k
            end if
      end do
      end do
      end do
!
! Finished:
!
end subroutine

subroutine setsupr( &
    nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,x,y,z, &
    vr,sec,MAXSBX,MAXSBY, &
    MAXSBZ,nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup, &
    ysizsup,nzsup,zmnsup,zsizsup)
!-----------------------------------------------------------------------
!
!           Establish Super Block Search Limits and Sort Data
!           *************************************************
!
! This subroutine sets up a 3-D "super block" model and orders the data
! by super block number.  The limits of the super block is set to the
! minimum and maximum limits of the grid; data outside are assigned to
! the nearest edge block.
!
! The idea is to establish a 3-D block network that contains all the
! relevant data.  The data are then sorted by their index location in
! the search network, i.e., the index location is given after knowing
! the block index in each coordinate direction (ix,iy,iz):
!          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
! An array, the same size as the number of super blocks, is constructed
! that contains the cumulative number of data in the model.  With this
! array it is easy to quickly check what data are located near any given
! location.
!
!
!
! INPUT VARIABLES:
!
!   nx,xmn,xsiz      Definition of the X grid being considered
!   ny,ymn,ysiz      Definition of the Y grid being considered
!   nz,zmn,zsiz      Definition of the Z grid being considered
!   nd               Number of data
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   vr(nd)           Variable at each location.
!   tmp(nd)          Temporary storage to keep track of the super block
!                      index associated to each data (uses the same
!                      storage already allocated for the simulation)
!   nsec             Number of secondary variables to carry with vr
!   sec1(nd)         First secondary variable (if nsec >= 1)
!   sec2(nd)         Second secondary variable (if nsec >= 2)
!   sec3(nd)         Third secondary variable (if nsec = 3)
!   MAXSB[X,Y,Z]     Maximum size of super block network
!
!
!
! OUTPUT VARIABLES:
!
!   nisb()                Array with cumulative number of data in each
!                           super block.
!   nxsup,xmnsup,xsizsup  Definition of the X super block grid
!   nysup,ymnsup,ysizsup  Definition of the Y super block grid
!   nzsup,zmnsup,zsizsup  Definition of the Z super block grid
!
!
!
! EXTERNAL REFERENCES:
!
!   sortem           Sorting routine to sort the data
!
!
!
!-----------------------------------------------------------------------
use geometry, only: getindx
implicit none
!arguments
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
real(kind=fk), intent(in) :: x(:),y(:),z(:)
real(kind=fk), intent(in) :: vr(:)
real(kind=fk), intent(in),optional :: sec(:,:)
integer(kind=gik), intent(in) :: MAXSBX,MAXSBy,MAXSBZ

integer(kind=gik), allocatable, intent(out) :: nisb(:)
integer(kind=gik), intent(out) :: nxsup,nysup,nzsup
real(kind=fk), intent(out) :: xmnsup,ymnsup,zmnsup
real(kind=fk), intent(out) :: xsizsup,ysizsup,zsizsup

!locals
real(kind=fk) tmp(size(vr))
integer i,nd,nsec,ii,ix,iy,iz,nsort
logical inflag

nd = size(x)
!
! Establish the number and size of the super blocks:
!
      nxsup   = min(nx,MAXSBX)
      nysup   = min(ny,MAXSBY)
      nzsup   = min(nz,MAXSBZ)
      xsizsup = real(nx)*xsiz/real(nxsup)
      ysizsup = real(ny)*ysiz/real(nysup)
      zsizsup = real(nz)*zsiz/real(nzsup)
      xmnsup  = (xmn-0.5*xsiz)+0.5*xsizsup
      ymnsup  = (ymn-0.5*ysiz)+0.5*ysizsup
      zmnsup  = (zmn-0.5*zsiz)+0.5*zsizsup
!
! Initialize the extra super block array to zeros:
!
  allocate(nisb(nxsup*nysup*nzsup))

  nisb = 0

!
! Loop over all the data assigning the data to a super block and
! accumulating how many data are in each super block:
!
      do i=1,nd
            call getindx(nxsup,xmnsup,xsizsup,x(i),ix,inflag)
            call getindx(nysup,ymnsup,ysizsup,y(i),iy,inflag)
            call getindx(nzsup,zmnsup,zsizsup,z(i),iz,inflag)
            ii = ix + (iy-1)*nxsup + (iz-1)*nxsup*nysup
            tmp(i)   = ii
            nisb(ii) = nisb(ii) + 1
      end do
!
! Sort the data by ascending super block number:
!
      nsort = 4 + nsec
      call sortem(1,nd,tmp,nsort,x,y,z,vr,sec)
!
! Set up array nisb with the starting address of the block data:
!
      do i=1,(nxsup*nysup*nzsup-1)
            nisb(i+1) = nisb(i) + nisb(i+1)
      end do
!
! Finished:
!
end subroutine

subroutine srchsupr( &
    xloc,yloc,zloc,radsqd,rotmat, &
    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct, &
    x,y,z,nisb,nxsup,xmnsup,xsizsup, &
    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup, &
    nclose,close,infoct)
!-----------------------------------------------------------------------
!
!              Search Within Super Block Search Limits
!              ***************************************
!
!
! This subroutine searches through all the data that have been tagged in
! the super block subroutine.  The close data are passed back in the
! index array "close".  An octant search is allowed.
!
!
!
! INPUT VARIABLES:
!
!   xloc,yloc,zloc   location of point being estimated/simulated
!   radsqd           squared search radius
!   irot             index of the rotation matrix for searching
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!   nsbtosr          Number of super blocks to search
!   ixsbtosr         X offsets for super blocks to search
!   iysbtosr         Y offsets for super blocks to search
!   izsbtosr         Z offsets for super blocks to search
!   noct             If >0 then data will be partitioned into octants
!   nd               Number of data
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   tmp(nd)          Temporary storage to keep track of the squared
!                      distance associated with each data
!   nisb()                Array with cumulative number of data in each
!                           super block.
!   nxsup,xmnsup,xsizsup  Definition of the X super block grid
!   nysup,ymnsup,ysizsup  Definition of the X super block grid
!   nzsup,zmnsup,zsizsup  Definition of the X super block grid
!
!
!
! OUTPUT VARIABLES:
!
!   nclose           Number of close data
!   close()          Index of close data
!   infoct           Number of informed octants (only computes if
!                      performing an octant search)
!
!
!
! EXTERNAL REFERENCES:
!
!   sqdist           Computes anisotropic squared distance
!   sortem           Sorts multiple arrays in ascending order
!
!
!
!-----------------------------------------------------------------------
use geometry, only: getindx, sqdist
implicit none
!arguments
integer(kind=gik), intent(in) :: nxsup,nysup,nzsup
real(kind=fk), intent(in) :: xmnsup,ymnsup,zmnsup
real(kind=fk), intent(in) :: xsizsup,ysizsup,zsizsup
integer(kind=gik), intent(in) :: nisb(:)
integer(kind=gik), intent(in) :: noct


real(kind=fk), intent(in) :: x(:),y(:),z(:)
real(kind=fk), intent(in) :: xloc,yloc,zloc
real(kind=fk), intent(in) :: rotmat(:,:),radsqd
integer(kind=gik),intent(in) :: nsbtosr, ixsbtosr(:),iysbtosr(:),izsbtosr(:)

integer(kind=gik),intent(inout) :: nclose
integer(kind=gik),intent(inout) :: infoct(:)
integer(kind=gik),intent(inout) :: close(:)
!locals
real(kind=dk) hsqd,tmp(nclose),dx,dy,dz,h
integer inoct(8),ix,iy,iz,ixsup,iysup,izsup,i,ii,isup,j,iq,na,nt,nums,nclose_int
logical inflag
!
! Determine the super block location of point being estimated:
!
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
!
! Loop over all the possible Super Blocks:
!
      nclose_int = 0
      do 1 isup=1,nsbtosr
!
! Is this super block within the grid system:
!
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
            if(ixsup.le.0.or.ixsup.gt.nxsup.or. &
               iysup.le.0.or.iysup.gt.nysup.or. &
               izsup.le.0.or.izsup.gt.nzsup) go to 1
!
! Figure out how many samples in this super block:
!
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
!
! Loop over all the data in this super block:
!
            do 2 ii=1,nums
                  i = i + 1
!
! Check squared distance:
!
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
!
! Accept this sample:
!
                  nclose_int = nclose_int + 1
                  close(nclose_int) = real(i)
                  tmp(nclose_int)  = real(hsqd)
 2          continue
 1    continue
!
! Sort the nearby samples by distance to point being estimated:
!
      call sortem(nclose_int,tmp,close)
!
! If we aren't doing an octant search then just return:
!
      if(noct.le.0) return
!
! PARTITION THE DATA INTO OCTANTS:
!
      do i=1,8
            inoct(i) = 0
      end do
!
! Now pick up the closest samples in each octant:
!
      nt = 8*noct
      na = 0
      do j=1,nclose_int
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
!
! Keep this sample if the maximum has not been exceeded:
!
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
!
! End of data selection. Compute number of informed octants and return:
!
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
!
! Finished:
!
end subroutine

end module
