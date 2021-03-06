module searching
use kinds_kit
use iso_c_binding
implicit none

type sb_structure_type
  integer(kind=gik) :: nxsup,nysup,nzsup
  real(kind=fk)    :: xsizsup,ysizsup,zsizsup
  real(kind=fk)    :: xmnsup,ymnsup,zmnsup
  real(kind=fk) rotmat(3,3)  
  integer(kind=gik), allocatable :: nisb(:)
  integer(kind=gik) , allocatable :: ixsbtosr(:),iysbtosr(:),izsbtosr(:)
  integer(kind=gik)  :: nsbtosr, ndata
  
  !for checking sizes
  integer(kind=gik)  :: nisb_len,ixsbtosr_len,iysbtosr_len,izsbtosr_len
  
end type

contains

!C_API
subroutine array_example(x,n) bind(c, name='array_example')
  implicit none
  real(kind=fk), intent(inout) :: x(n)
  integer(kind=gik), intent(in), value :: n
  print *,"x",size(x)
  print *,"n",n
  
  x(1) = 999

end subroutine


!C_API
subroutine super_block_print(sb_p) bind(c, name='super_block_print')
  implicit none
  type(c_ptr), intent(in), value :: sb_p
  !
  type(sb_structure_type),pointer :: sb_structure
  call c_f_pointer(sb_p, sb_structure)  
  
  print *, sb_structure%nxsup,sb_structure%nysup,sb_structure%nzsup

end subroutine


function setup_superblock_c(x,y,z,vr,sec,n,nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,sang1,sang2,sang3,sanis1,sanis2,radsqd,MAXSBX,MAXSBY,MAXSBZ) bind(c, name='setup_superblock')
  implicit none
  real(kind=fk), intent(inout) :: x(n),y(n),z(n)
  real(kind=fk), intent(inout) :: vr(n)
  real(kind=fk), intent(inout),optional :: sec(n)
  integer(kind=gik), intent(in),value :: n
  integer(kind=gik), intent(in),value :: nx,ny,nz
  real(kind=fk), intent(in),value    :: xsiz,ysiz,zsiz
  real(kind=fk), intent(in),value    :: xmn,ymn,zmn
  real(kind=fk), intent(in) :: sang1,sang2,sang3,sanis1,sanis2
  real(kind=fk), intent(in),value :: radsqd
  integer(kind=gik), intent(in), value :: MAXSBX,MAXSBy,MAXSBZ
!return
  type(c_ptr) :: setup_superblock_c
!locals
  type(sb_structure_type) :: sb_structure

  print *,"grid",nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz
  print *,"x",size(x)
  print *,"y",size(y)
  print *,"z",size(z)
  
  sb_structure = setup_superblock(x,y,z,vr,sec,nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,sang1,sang2,sang3,sanis1,sanis2,radsqd,MAXSBX,MAXSBY,MAXSBZ)
  
  setup_superblock_c = c_loc(sb_structure)
end function

integer(kind=gik) function search_super_block_c(sb_p, &
        xloc,yloc,zloc,radsqd,&
        ndmax,noct,x,y,z,nd,nclose,close,infoct) bind(c, name='search_super_block')
implicit none
!arguments
  type(c_ptr), intent(in), value :: sb_p
  real(kind=fk), intent(in), value :: xloc,yloc,zloc
  real(kind=fk), intent(in), value :: radsqd
  integer(kind=gik), intent(in), value :: ndmax,noct
  real(kind=fk), intent(in) :: x(nd),y(nd),z(nd)
  integer(kind=gik), intent(in),value :: nd
  integer(kind=gik), intent(in) :: nclose
  integer(kind=gik), intent(inout) :: close(nclose),infoct(8)
!locals
  type(sb_structure_type),pointer :: sb_structure
  integer(kind=gik) :: buffer(nd)
  integer(kind=gik) n_found
!begin
  call c_f_pointer(sb_p, sb_structure)  
  
  call search_super_block(xloc,yloc,zloc,radsqd,sb_structure,ndmax,noct,x,y,z,n_found,buffer,infoct)
  close(1:n_found) = buffer(1:n_found)
  
  search_super_block_c = n_found
  
end function

!=======================================================================
!Fortran side
function setup_superblock(x,y,z,vr,sec,nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,sang1,sang2,sang3,sanis1,sanis2,radsqd,MAXSBX,MAXSBY,MAXSBZ)
use geometry, only: grid_type,setrot
implicit none
!arguments
  real(kind=fk), intent(inout) :: x(:),y(:),z(:)
  real(kind=fk), intent(inout) :: vr(:)
  real(kind=fk), intent(inout),optional :: sec(:)
  integer(kind=gik) :: nx,ny,nz
  real(kind=fk)    :: xsiz,ysiz,zsiz
  real(kind=fk)    :: xmn,ymn,zmn
  real(kind=fk), intent(in) :: sang1,sang2,sang3,sanis1,sanis2,radsqd
  integer(kind=gik), intent(in) :: MAXSBX,MAXSBy,MAXSBZ

!return
  type(sb_structure_type) :: setup_superblock
!locals
  type(sb_structure_type) :: sb_structure
  integer(kind=gik) MAXSB
  integer(kind=gik) test
  
  call setrot(sang1,sang2,sang3,sanis1,sanis2,sb_structure%rotmat)  

  MAXSB = min(nx,MAXSBX)*min(ny,MAXSBY)*min(nz,MAXSBZ)
  !determine
  allocate(sb_structure%nisb(MAXSB),stat = test)
  
  sb_structure%nisb_len = MAXSB
  sb_structure%ndata = size(vr)
  

#ifdef TRACE                  
  print *,'setup_superblock::','before setsupr'
#endif
  call setsupr( &
      nx,xmn,xsiz, &
      ny,ymn,ysiz, &
      nz,zmn,zsiz, &
      x,y,z, &
      vr,sec,MAXSBX,MAXSBY,MAXSBZ, &
      sb_structure%nisb, &
      sb_structure%nxsup,sb_structure%xmnsup,sb_structure%xsizsup, &
      sb_structure%nysup,sb_structure%ymnsup,sb_structure%ysizsup, &
      sb_structure%nzsup,sb_structure%zmnsup,sb_structure%zsizsup)

#ifdef TRACE                  
  print *,'setup_superblock::','after setsupr'
  print *,'setup_superblock xdim::',sb_structure%nxsup,sb_structure%xmnsup,sb_structure%xsizsup
  print *,'setup_superblock ydim::',sb_structure%nysup,sb_structure%ymnsup,sb_structure%ysizsup
  print *,'setup_superblock zdim::',sb_structure%nzsup,sb_structure%zmnsup,sb_structure%zsizsup
  print *,'setup_superblock::','MAXSB=',MAXSB
#endif
  allocate(sb_structure%ixsbtosr(8 * MAXSB),stat = test)
  allocate(sb_structure%iysbtosr(8 * MAXSB),stat = test)
  allocate(sb_structure%izsbtosr(8 * MAXSB),stat = test)

  sb_structure%ixsbtosr_len = size(sb_structure%ixsbtosr)
  sb_structure%iysbtosr_len = size(sb_structure%iysbtosr)
  sb_structure%izsbtosr_len = size(sb_structure%izsbtosr)

#ifdef TRACE                  
  print *,'setup_superblock::','before picksup'
#endif
  call picksup(sb_structure%nxsup,sb_structure%xsizsup, &
               sb_structure%nysup,sb_structure%ysizsup, &
               sb_structure%nzsup,sb_structure%zsizsup, &
               sb_structure%rotmat,radsqd, &
               sb_structure%nsbtosr,sb_structure%ixsbtosr, &
               sb_structure%iysbtosr,sb_structure%izsbtosr)
#ifdef TRACE                  
  print *,'setup_superblock::','after picksup'
  print *,'setup_superblock::','nsbtosr',sb_structure%nsbtosr
#endif
  setup_superblock = sb_structure
end function

subroutine search_super_block( &
    xloc,yloc,zloc,radsqd, &
    sb_structure, &
    ndmax,noct, &
    x,y,z, &
    nclose,close,infoct)
use geometry, only: grid_type
implicit none
!arguments
  real(kind=fk), intent(in) :: xloc,yloc,zloc
  real(kind=fk), intent(in) :: radsqd
  type(sb_structure_type), intent(in) :: sb_structure
  integer(kind=gik), intent(in) :: ndmax,noct
  real(kind=fk), intent(in) :: x(:),y(:),z(:)
  integer(kind=gik), intent(out) :: nclose
  integer(kind=gik), intent(inout) :: close(:),infoct(:)

!locals

  nclose = ndmax

  call srchsupr( &
      xloc,yloc,zloc,radsqd,sb_structure%rotmat, &
      sb_structure%nsbtosr, &
      sb_structure%ixsbtosr,sb_structure%iysbtosr,sb_structure%izsbtosr, &
      noct, &
      x,y,z,&
      sb_structure%nisb, &
      sb_structure%nxsup,sb_structure%xmnsup,sb_structure%xsizsup, &
      sb_structure%nysup,sb_structure%ymnsup,sb_structure%ysizsup, &
      sb_structure%nzsup,sb_structure%zmnsup,sb_structure%zsizsup, &
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
integer(kind=gik) i,j,k,i1,j1,k1,i2,j2,k2
real(kind=fk) xo,yo,zo,xdis,ydis,zdis,hsqd,shortest
!
! MAIN Loop over all possible super blocks:
!
#ifdef TRACE
    print *, "picksup::n_sup",nxsup,nysup,nzsup
    print *, "picksup::_sizsup",xsizsup,ysizsup,zsizsup
    print *, "picksup::radsqd",radsqd,rotmat
#endif

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

#ifdef TRACE
            print *, "picksup::shortest",shortest
#endif
!
! Keep this super block if it is close enoutgh:
!
            if(real(shortest).le.radsqd) then
                  nsbtosr = nsbtosr + 1
                  if (nsbtosr > size(ixsbtosr)) STOP "nsbtosr > size(ixsbtosr)"
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
    nx,xmn,xsiz,&
    ny,ymn,ysiz,&
    nz,zmn,zsiz,&
    x,y,z, &
    vr,sec,MAXSBX,MAXSBY, &
    MAXSBZ,nisb,&
    nxsup,xmnsup,xsizsup,&
    nysup,ymnsup,ysizsup,&
    nzsup,zmnsup,zsizsup)
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
use sorting, only: sortem_original
implicit none
!arguments
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
real(kind=fk), intent(inout) :: x(:),y(:),z(:)
real(kind=fk), intent(inout) :: vr(:)
real(kind=fk), intent(inout),optional :: sec(:)
integer(kind=gik), intent(in) :: MAXSBX,MAXSBy,MAXSBZ

integer(kind=gik), intent(inout) :: nisb(:)
integer(kind=gik), intent(out) :: nxsup,nysup,nzsup
real(kind=fk), intent(out) :: xmnsup,ymnsup,zmnsup
real(kind=fk), intent(out) :: xsizsup,ysizsup,zsizsup

!locals
real(kind=fk) tmp(size(vr))
integer(kind=gik) i,nd,nsec,ii,ix,iy,iz,nsort
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
      do i=1,nxsup*nysup*nzsup
        nisb(i) = 0
      end do

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
            sec(i) = i
      end do
!
! Sort the data by ascending super block number:
!
      nsort = 4 + nsec
      if (size(sec) > 0) then
        call sortem_original(1,nd,tmp,x,y,z,vr,sec)
      else
        call sortem_original(1,nd,tmp,x,y,z,vr)
      end if
#ifdef TRACE_SORTED
      print *,"SORTED XYZ"
      do i=1,nd
        print *,x(i),y(i),z(i),vr(i)
      end do
#endif
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
    x,y,z,nisb,&
    nxsup,xmnsup,xsizsup, &
    nysup,ymnsup,ysizsup, &
    nzsup,zmnsup,zsizsup, &
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
use sorting, only: sortem_real_int
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

integer(kind=gik),intent(inout) :: nclose,infoct(:),close(:)
!locals
real(kind=dk) hsqd,dx,dy,dz,h
real(kind=fk) tmp(size(x)),tmp_close(size(x)) !size 
integer(kind=gik) inoct(8),ix,iy,iz,ixsup,iysup,izsup,i,ii,isup,j,iq,na,nt,nums,nclose_int
logical inflag
!
! Determine the super block location of point being estimated:
!
  ix=0
  iy=0
  iz=0
  
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
!
! Loop over all the possible Super Blocks:
!
#ifdef TRACE                  
    print *,"srchsupr::nsbtosr",nsbtosr,xloc,yloc,zloc,ix,iy,iz
    print *,"srchsupr::n_sup",nxsup,nysup,nzsup
#endif


      nclose_int = 0
      do 1 isup=1,nsbtosr
!
! Is this super block within the grid system:
!
#ifdef CHECK_BOUNDS
            if (isup > size(ixsbtosr)) stop "isup > size(ixsbtosr)"
            if (isup > size(iysbtosr)) stop "isup > size(ixsbtosr)"
            if (isup > size(izsbtosr)) stop "isup > size(ixsbtosr)"
#endif

            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)


#ifdef TRACE                  
            print *,"srchsupr::isup",isup
            print *,"srchsupr::i_sup",ixsup,iysup,izsup
            print *,"srchsupr::i_sbtosr",ixsbtosr(isup),iysbtosr(isup),izsbtosr(isup)
#endif

            if(ixsup.le.0.or.ixsup.gt.nxsup.or. &
               iysup.le.0.or.iysup.gt.nysup.or. &
               izsup.le.0.or.izsup.gt.nzsup) go to 1
!
! Figure out how many samples in this super block:
!
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup

#ifdef CHECK_BOUNDS
            if (ii > size(nisb)) stop "ii > size(nisb)"
#endif

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
                  
#ifdef TRACE                  
                  
                  print *,"loc",xloc,yloc,zloc
                  print *,"x(i)",i,x(i),y(i),z(i)
                  print *,"rotmat",rotmat
                  print *,"hsqd",hsqd
                  print *,"radsqd",radsqd
                  print *,"nclose_int",nclose_int
#endif
                  
                  if(real(hsqd).gt.radsqd) go to 2
!
! Accept this sample:
!
                  nclose_int = nclose_int + 1

#ifdef CHECK_BOUNDS
                  if (nclose_int > size(close)) stop "nclose_int > size(close)"
                  if (nclose_int > size(tmp)) stop "nclose_int > size(tmp)"
#endif

                  close(nclose_int) = i
                  tmp(nclose_int)  = real(hsqd)
 2          continue
 1    continue
!
! Sort the nearby samples by distance to point being estimated:
!

      if (size(tmp) /= size(close)) stop "size(tmp) /= size(close)"
      
      !tmp_close = real(close,kind=fk)

      call sortem_real_int(1,nclose_int,tmp,close)

      !close(1:nclose_int) = int(tmp_close(1:nclose_int),kind=gik)
      
      nclose = nclose_int


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
            i  = close(j)
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
