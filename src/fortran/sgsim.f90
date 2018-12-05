module sgsimulation
use kinds_kit

real,parameter :: TINY = 0.0001
real,parameter :: UNEST = -999

type ctable_structure_type
  integer(kind=gik) :: ncnode,nlooku
  integer(kind=gik),allocatable :: icnode(:),ixnode(:),iynode(:),iznode(:)
  real(kind=fk),allocatable :: covtab(:,:,:)
  integer(kind=gik) :: MAXCTX,MAXCTY,MAXCTZ
end type

contains

subroutine setup_ctable(ctable_st, &
  nx,ny,nz, &
  xmn,ymn,zmn, &
  xsiz,ysiz,zsiz, &
  nst,c0,cc,aa,it,vrotmat, &
  radsqd,nodmax,&
  MAXCTX,MAXCTY,MAXCTZ)
!arguments
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
real(kind=fk), intent(in) :: radsqd
integer(kind=gik), intent(in) ::nst,it(:),nodmax
real(kind=fk), intent(in) :: c0,cc(:),aa(:),vrotmat(:,:,:)
integer(kind=gik), intent(in) ::MAXCTX,MAXCTY,MAXCTZ
type(ctable_structure_type),intent(inout) :: ctable_st
!locals
  integer(kind=gik) MAXXYZ,MAXCXY
  integer test

  MAXCXY = MAXCTX * MAXCTY
  MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
  allocate(ctable_st%covtab(MAXCTX,MAXCTY,MAXCTZ))
  allocate(ctable_st%icnode(MAXNOD),stat = test)
  allocate(ctable_st%ixnode(MAXXYZ),stat = test)
  allocate(ctable_st%iynode(MAXXYZ),stat = test)
  allocate(ctable_st%iznode(MAXXYZ),stat = test)

  ctable_st%MAXCTX = MAXCTX
  ctable_st%MAXCTY = MAXCTY
  ctable_st%MAXCTZ = MAXCTZ

  call ctable( &
      nx,ny,nz, &
      xmn,ymn,zmn, &
      xsiz,ysiz,zsiz, &
      nst,c0,cc,aa,it,vrotmat, &
      radsqd,nodmax,&
      MAXCTX,MAXCTY,MAXCTZ,&
      ctable_st%covtab,ctable_st%nlooku,ctable_st%ixnode,ctable_st%iynode,ctable_st%iznode)
end subroutine

subroutine search_ctable(ix,iy,iz,sim, &
    nx,ny,nz, xmn,ymn,zmn, xsiz,ysiz,zsiz, &
    noct,nodmax, &
    ctable_st, &
    ncnode,icnode,cnodex,cnodey,cnodez,cnodev)
!arguments
  integer(kind=gik), intent(in) :: ix,iy,iz
  real(kind=fk), intent(in) ::sim(:)
  integer(kind=gik), intent(in) :: nx,ny,nz
  real(kind=fk), intent(in) :: xmn,ymn,zmn
  real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
  integer(kind=gik), intent(in) ::nodmax,noct
  type(ctable_structure_type),intent(in) :: ctable_st
  integer(kind=gik), intent(inout) ::ncnode,icnode(:)
  real(kind=fk), intent(inout) :: cnodex(:),cnodey(:),cnodez(:),cnodev(:)
!locals
  integer(kind=gik) :: nctx,ncty,nctz


  nctx = min(((MAXCTX-1)/2),(nx-1))
  ncty = min(((MAXCTY-1)/2),(ny-1))
  nctz = min(((MAXCTZ-1)/2),(nz-1))


  call srchnd(ix,iy,iz, &
    sim, &
    nx,ny,nz, &
    xmn,ymn,zmn, &
    xsiz,ysiz,zsiz, &
    noct,nodmax, &
    nctx,ncty,nctz,&
    ctable_st%nlooku, &
    ctable_st%ixnode,ctable_st%iynode,ctable_st%iznode, &
    ncnode,icnode,cnodex,cnodey,cnodez,cnodev)
end subroutine

subroutine sgsim_grid( &
    x,y,z, & ! coordinates
    vr, &    ! variable
    nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
    ktype,radius,srotmat, & !search configuration
    ndmin,ndmax,noct,nodmax, & !simulated data search
    nst,c0,cc,aa,it,vrotmat, & !variogram
    ctable_st, & !cache for variogram
    seed, &
    sim & !output
  )
use geometry
use searching, only: sb_structure_type,setup_superblock,search_super_block
implicit none
!arguments
real(kind=fk), intent(in) :: x(:),y(:),z(:)
real(kind=fk), intent(in) :: vr(:)
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=sk), intent(in) :: ktype
real(kind=fk), intent(in) :: radius,srotmat(:,:)
integer(kind=sk), intent(in) :: ndmin,ndmax,noct,nodmax
integer(kind=gik), intent(in) ::nst,it(:)
real(kind=fk), intent(in) :: c0,cc(:),aa(:),vrotmat(:,:,:)
type(ctable_structure_type),intent(in) :: ctable_st
integer(kind=sk), intent(in) :: seed
real(kind=fk), intent(in) :: sim(:)

!locals
type(sb_structure_type) :: sb_strcture


!sb_strcture = setup_superblock(x,y,z,vr,sec,grid,sb_grid,rotmat,radsqd,MAXSBX,MAXSBY,MAXSBZ)

end subroutine

subroutine sgsim_realizarion_grid( &
    x,y,z, & ! coordinates
    vr, &    ! variable
    nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
    ktype,radius,srotmat, & !search configuration
    ndmin,ndmax,noct, & !
    nodmax,nmult, sstrat, & !simulated data search
    sb_strcture, &
    nst,c0,cc,aa,it,vrotmat, & !variogram
    ctable_st,&
    seed,&
    sim) !output
use geometry
use searching, only: sb_structure_type,setup_superblock,search_super_block
implicit none
!arguments
real(kind=fk), intent(in) :: x(:),y(:),z(:)
real(kind=fk), intent(in) :: vr(:)
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=gik), intent(in) :: ktype
real(kind=fk), intent(in) :: radius,srotmat(:,:)
integer(kind=gik), intent(in) :: ndmin,ndmax,noct,nodmax,nmult,sstrat
type(sb_structure_type),optional :: sb_strcture
integer(kind=gik), intent(in) ::nst,it(:)
real(kind=fk), intent(in) :: c0,cc(:),aa(:),vrotmat(:,:,:)
type(ctable_structure_type),intent(in) :: ctable_st
integer(kind=gik), intent(in) :: seed
real(kind=fk), intent(inout) :: sim(:)
!locals
real      var(10)
real*8    p,acorni,cp,oldcp,w
logical   testind, trace
integer nxyz,imult,ind,id,in
integer ix,iy,iz,nnx,nny,nnz,nd,nxy
real    xx,yy,zz,test,test2
integer(kind=gik) order(nx*ny*nz)
real(kind=dk),allocatable :: A(:,:),b(:),s(:)
real(kind=fk) :: radsqd
real(kind=dk) :: xp
real    gmean,cmean,cstdev
integer id2,idbg,idum,index,irepo,jx,jy,jz,ldbg,ierr
integer(kind=gik) nclose,close(ndmax+nodmax),infoct(8),lktype
integer(kind=gik) ncnode,icnode(nodmax) !results of searching simulated nodes
real(kind=fk) :: cnodex(nodmax),cnodey(nodmax),cnodez(nodmax),cnodev(nodmax)

radsqd = radius  * radius
idbg = 0

nd = size(vr)
nxy = nx*ny
nxyz = nxy*nz

!
! Work out a random path for this realization:
!
      do ind=1,nxyz
            sim(ind)   = real(acorni(idum))
            order(ind) = ind
      end do
!
! The multiple grid search works with multiples of 4 (yes, that is
! somewhat arbitrary):
!
      if(nmult > 0) then
            do imult=1,nmult
                  nnz = max(1,nz/(imult*4))
                  nny = max(1,ny/(imult*4))
                  nnx = max(1,nx/(imult*4))
                  jz  = 1
                  jy  = 1
                  jx  = 1
                  do iz=1,nnz
                     if(nnz.gt.1) jz = iz*imult*4
                     do iy=1,nny
                        if(nny.gt.1) jy = iy*imult*4
                        do ix=1,nnx
                           if(nnx.gt.1) jx = ix*imult*4
                           index = jx + (jy-1)*nx + (jz-1)*nxy
                           sim(index) = sim(index) - imult
                        end do
                     end do
                  end do
            end do
      end if
      call sortem(1,nxyz,sim,order)
!
! Initialize the simulation:
!
      do ind=1,nxyz
            sim(ind) = UNEST
      end do
!
! Assign the data to the closest grid node:
!

      do id=1,nd
            call getindx(nx,xmn,xsiz,x(id),ix,testind)
            call getindx(ny,ymn,ysiz,y(id),iy,testind)
            call getindx(nz,zmn,zsiz,z(id),iz,testind)
            ind = ix + (iy-1)*nx + (iz-1)*nxy
            xx  = xmn + real(ix-1)*xsiz
            yy  = ymn + real(iy-1)*ysiz
            zz  = zmn + real(iz-1)*zsiz
            test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
!
! Assign this data to the node (unless there is a closer data):
!
            if(sstrat.eq.1) then
                  if(sim(ind).ge.0.0) then
                        id2 = int(sim(ind)+0.5)
                        test2 = abs(xx-x(id2)) + abs(yy-y(id2)) + abs(zz-z(id2))
                        if(test.le.test2) sim(ind) = real(id)
                        write(ldbg,102) id,id2
                  else
                        sim(ind) = real(id)
                  end if
            end if
!
! Assign a flag so that this node does not get simulated:
!
            if(sstrat.eq.0.and.test.le.TINY) sim(ind)=10.0*UNEST
      end do
102        format(' WARNING data values ',2i5,' are both assigned to ', &
            /,'         the same node - taking the closest')
!
! Now, enter data values into the simulated grid:
!
      do ind=1,nxyz
            id = int(sim(ind)+0.5)
            if(id.gt.0) sim(ind) = vr(id)
      end do
      irepo = max(1,min((nxyz/10),10000))
!
! MAIN LOOP OVER ALL THE NODES:
!
  do in=1,nxyz
    if((int(in/irepo)*irepo).eq.in) write(*,103) in
103              format('   currently on node ',i9)
!
! Figure out the location of this point and make sure it has
! not been assigned a value already:
!
    index = order(in)
    if(sim(index).gt.(UNEST+EPSLON) .or. sim(index).lt.(UNEST*2.0)) go to 5
    iz = int((index-1)/nxy) + 1
    iy = int((index-(iz-1)*nxy-1)/nx) + 1
    ix = index - (iz-1)*nxy - (iy-1)*nx
    xx = xmn + real(ix-1)*xsiz
    yy = ymn + real(iy-1)*ysiz
    zz = zmn + real(iz-1)*zsiz
!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously
! simulated grid nodes:
!
    if(sstrat.eq.0) then
      call search_super_block( &
          xx,yy,zz,radsqd,srotmat, &
          sb_strcture, &
          ndmax,noct, &
          x,y,z, &
          nclose,close,infoct)

      if(nclose.lt.ndmin) go to 5
      if(nclose.gt.ndmax) nclose = ndmax
    endif

    call search_ctable(ix,iy,iz,sim, &
        nx,ny,nz, xmn,ymn,zmn, xsiz,ysiz,zsiz, &
        noct,nodmax, &
        ctable_st, &
        ncnode,icnode,cnodex,cnodey,cnodez,cnodev)
!
! Calculate the conditional mean and standard deviation.  This will be
! done with kriging if there are data, otherwise, the global mean and
! standard deviation will be used:
!
    if(ktype.eq.2) then
      !TODO gmean = lvm(index)
    else
      gmean = 0.0
    end if
    if((nclose+ncnode).lt.1) then
      cmean  = gmean
      cstdev = 1.0
    else
!
! Perform the kriging.  Note that if there are fewer than four data
! then simple kriging is prefered so that the variance of the
! realization does not become artificially inflated:
!
      lktype = ktype
      if(ktype.eq.1.and.(nclose+ncnode).lt.4)lktype=0
      call krige(ix,iy,iz,xx,yy,zz,lktype,gmean,nclose,close,ncnode,ctable_st,cmean,cstdev)
    endif
!
! Draw a random number and assign a value to this node:
!
    p = acorni(idum)
    call gauinv(p,xp,ierr)
    sim(index) = xp * cstdev + cmean
    if(idbg.ge.3) write(ldbg,141) p,sim(index)
141              format(' random number ',f6.4,' realization ',f7.4)
!
! Quick check for far out results:
!
    if(abs(cmean).gt.5.0.or.abs(cstdev).gt.5.0.or. abs(sim(index)).gt.6.0) then
            write(ldbg,104) ix,iy,iz,cmean,cstdev,sim(index)
104             format('WARNING: grid node location: ',3i5,/, &
                    '         conditional mean:   ',f12.5,/,  &
                    '         conditional stdev:  ',f12.5,/,  &
                    '         simulated value:    ',f12.5)
    endif
!
! END MAIN LOOP OVER NODES:
!
5                continue
  end do
!
! Do we need to reassign the data to the grid nodes?
!
  if(sstrat.eq.0) then
        do id=1,nd
              call getindx(nx,xmn,xsiz,x(id),ix,testind)
              call getindx(ny,ymn,ysiz,y(id),iy,testind)
              call getindx(nz,zmn,zsiz,z(id),iz,testind)
              xx  = xmn + real(ix-1)*xsiz
              yy  = ymn + real(iy-1)*ysiz
              zz  = zmn + real(iz-1)*zsiz
              ind = ix + (iy-1)*nx + (iz-1)*nxy
              test=abs(xx-x(id))+abs(yy-y(id))+abs(zz-z(id))
              if(test.le.TINY) sim(ind) = vr(id)
        end do
  end if
end subroutine

subroutine krige( &
      ix,iy,iz,xx,yy,zz, &
      ktype, &
      gmean, &
      nclose,close,ncnode, & !conditioning data
      ctable_st, &
      cmean, &
      cstdev)
!arguments
integer(kind=gik), intent(in) :: ix,iy,iz
real(kind=fk), intent(in) :: xx,yy,zz
real(kind=fk), intent(out) :: gmean
integer(kind=gik), intent(in) :: ktype
integer(kind=gik), intent(in) :: nclose,close(:),ncnode
type(ctable_structure_type),intent(in) :: ctable_st
real(kind=fk), intent(out) :: cmean,cstdev
!locals
logical first
integer na,neq,i,j,k,ierr
real(kind=fk), allocatable :: a(:,:),b(:),s(:),r(:),rr(:)
integer(kind=gik) :: lktype,ind
!
! Size of the kriging system:
!
lktype = ktype
first = .false.
na    = nclose + ncnode
33   continue
if(lktype.eq.0) neq = na
if(lktype.eq.1) neq = na + 1
if(lktype.eq.2) neq = na
if(lktype.eq.3) neq = na + 2
if(lktype.eq.4) neq = na + 1
if(lktype.ge.3) then
      ind = ix + (iy-1)*nx + (iz-1)*nxy
      if(lvm(ind).le.-6.0.or.lvm(ind).ge.6.0) then
            lktype = 0
            go to 33
      end if
end if

allocate(a(na,na),b(na),s(na),r(na),rr(na))
!
! Set up kriging matrices:
!
!samples part
do j=1,nclose
  do i=j,nclose
    call cova3(x(close(i)),y(close(i)),z(close(i)),&
               x(close(j)),y(close(j)),z(close(j)),&
               nst,c0,it,cc,aa,rotmat,cmax,cov)
    a(i,j) = dble(cov)
    a(j,i) = dble(cov)

  end do

  call cova3(xx,yy,zz,x(close(i)),y(close(i)),z(close(i)),&
             nst,c0,it,cc,aa,rotmat,cmax,cov)
  b(j) = dble(cov)

end do

!simulated part
do j=nclose+1,na
  do i=1,nclose+1,na
    call cova3(x(close(i)),y(close(i)),z(close(i)),&
               x(close(j)),y(close(j)),z(close(j)),&
               nst,c0,it,cc,aa,rotmat,cmax,cov)
    a(i,j) = dble(cov)
    a(j,i) = dble(cov)

  end do
  call cova3(xx,yy,zz,xx,yy,zz,&
             nst,c0,it,cc,aa,rotmat,cmax,cov)
  b(j) = dble(cov)
end do

!samples against simulated part
do j=1,nclose
  do i=1,nclose+1,na
    call cova3(x(close(i)),y(close(i)),z(close(i)),&
               x(close(j)),y(close(j)),z(close(j)),&
               nst,c0,it,cc,aa,rotmat,cmax,cov)
    a(i,j) = dble(cov)
    a(j,i) = dble(cov)
  end do
  call cova3(xx,yy,zz,xx,yy,zz,&
             nst,c0,it,cc,aa,rotmat,cmax,cov)
  b(j) = dble(cov)
end do

!
! Addition of OK constraint:
!
if(lktype.eq.1.or.lktype.eq.3) then
      do i=1,na
        a(i,na) = 1.0
        a(na,i) = 1.0
      end do
      a(na,na)    = 0.0
      r(na+1)  = 1.0
      rr(na+1) = 1.0
endif
!
! Addition of the External Drift Constraint:
!
if(lktype.eq.3) then
      ! edmin =  999999.
      ! edmax = -999999.
      ! do i=1,na
      !       in    = in + 1
      !       a(in) = vrea(i)
      !       if(a(in).lt.edmin) edmin = a(in)
      !       if(a(in).gt.edmax) edmax = a(in)
      ! end do
      ! in       = in + 1
      ! a(in)    = 0.0
      ! in       = in + 1
      ! a(in)    = 0.0
      ! ind      = ix + (iy-1)*nx + (iz-1)*nxy
      ! r(na+2)  = dble(lvm(ind))
      ! rr(na+2) = r(na+2)
      ! if((edmax-edmin).lt.EPSLON) neq = neq - 1
endif
!
! Addition of Collocated Cosimulation Constraint:
!
if(lktype.eq.4) then
      ! sfmin =  1.0e21
      ! sfmax = -1.0e21
      ! do i=1,na
      !       in    = in + 1
      !       a(in) = dble(colocorr)*r(i)
      !       if(a(in).lt.sfmin) sfmin = a(in)
      !       if(a(in).gt.sfmax) sfmax = a(in)
      ! end do
      ! in    = in + 1
      ! a(in) = 1.0
      ! ii    = na + 1
      ! r(ii) = dble(colocorr)
      ! rr(ii)= r(ii)
!           if((sfmax-sfmin).lt.EPSLON) neq = neq - 1
end if
!
! Write out the kriging Matrix if Seriously Debugging:
!
if(idbg.ge.3) then
      write(ldbg,100) ix,iy,iz
      is = 1
      do i=1,neq
            ie = is + i - 1
            !!!write(ldbg,101) i,r(i),(a(j),j=is,ie)
            is = is + i
      end do
100        format(/,'Kriging Matrices for Node: ',3i4,' RHS first')
101        format('    r(',i2,') =',f7.4,'  a= ',99f7.4)
endif
!
! Solve the Kriging System:
!
if(neq.eq.1.and.lktype.ne.3) then
      s(1)  = r(1) / a(1,1)
      ising = 0
else
      call ksol(a,r,s,ising)
endif
!
! Write a warning if the matrix is singular:
!
if(ising.ne.0) then
      if(idbg.ge.1) then
            write(ldbg,*) 'WARNING SGSIM: singular matrix'
            write(ldbg,*) '               for node',ix,iy,iz
      endif
      cmean  = gmean
      cstdev = 1.0
      return
endif
!
! Compute the estimate and kriging variance.  Recall that kriging type
!     0 = Simple Kriging:
!     1 = Ordinary Kriging:
!     2 = Locally Varying Mean:
!     3 = External Drift:
!     4 = Collocated Cosimulation:
!
cmean  = 0.0
cstdev = cbb
sumwts = 0.0
do i=1,na
      cmean  = cmean  + real(s(i))*vra(i)
      cstdev = cstdev - real(s(i)*rr(i))
      sumwts = sumwts + real(s(i))
end do

if(lktype.eq.1) cstdev = cstdev - real(s(na+1))

if(lktype.eq.2) cmean  = cmean + gmean

if(lktype.eq.4) then
      ind    = ix + (iy-1)*nx + (iz-1)*nxy
      cmean  = cmean  + real(s(na+1))*lvm(ind)
      cstdev = cstdev - real(s(na+1) *rr(na+1))
end if
!
! Error message if negative variance:
!
if(cstdev.lt.0.0) then
      write(ldbg,*) 'ERROR: Negative Variance: ',cstdev
      cstdev = 0.0
endif
cstdev = sqrt(max(cstdev,0.0))
!
! Write out the kriging Weights if Seriously Debugging:
!
if(idbg.ge.3) then
      do i=1,na
            write(ldbg,140) i,vra(i),s(i)
      end do
140        format(' Data ',i4,' value ',f8.4,' weight ',f8.4)
      if(lktype.eq.4) write(ldbg,141) lvm(ind),s(na+1)
141        format(' Sec Data  value ',f8.4,' weight ',f8.4)
      write(ldbg,142) gmean,cmean,cstdev
142        format(' Global mean ',f8.4,' conditional ',f8.4,' std dev ',f8.4)
end if
!
! Finished Here:
!
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


subroutine ctable( &
    nx,ny,nz, &
    xmn,ymn,zmn, &
    xsiz,ysiz,zsiz, &
    nst,c0,cc,aa,it,vrotmat, &
    radsqd,nodmax,&
    MAXCTX,MAXCTY,MAXCTZ,&
    covtab,nlooku,ixnode,iynode,iznode)
use variography, only: cova3
implicit none
!arguments
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=gik), intent(in) ::nst,it(:)
real(kind=fk), intent(in) :: c0,cc(:),aa(:),vrotmat(:,:,:)
real(kind=fk), intent(in) :: radsqd
real(kind=fk), intent(inout) :: covtab(:,:,:)
integer(kind=gik), intent(out) ::nlooku
integer(kind=gik), intent(inout) ::ixnode(:),iynode(:),iznode(:)
integer(kind=gik), intent(in) ::nodmax
integer(kind=gik), intent(in) ::MAXCTX,MAXCTY,MAXCTZ

!locals
real,parameter :: TINY=1.0e-10
real(kind=fk)::tmp(MAXCTX * MAXCTY * MAXCTZ),order(MAXCTX * MAXCTY * MAXCTZ)
real(kind=fk)::cmax,cbb
integer nctx,ncty,nctz,i,j,k,MAXCXY,ic,jc,kc,il,ix,iy,iz,loc
real xx,yy,zz

real*8    hsqd,sqdist

MAXCXY = MAXCTX * MAXCTY


!
! Size of the look-up table:
!
nctx = min(((MAXCTX-1)/2),(nx-1))
ncty = min(((MAXCTY-1)/2),(ny-1))
nctz = min(((MAXCTZ-1)/2),(nz-1))

!
! NOTE: If dynamically allocating memory, and if there is no shortage
!       it would a good idea to go at least as far as the radius and
!       twice that far if you wanted to be sure that all covariances
!       in the left hand covariance matrix are within the table look-up.
!
! Initialize the covariance subroutine and cbb at the same time:
!
call cova3(0.0,0.0,0.0,0.0,0.0,0.0,nst,c0,it,cc,aa,vrotmat,cmax,cbb)
!
! Now, set up the table and keep track of the node offsets that are
! within the search radius:
!
nlooku = 0
do i=-nctx,nctx
xx = i * xsiz
ic = nctx + 1 + i
do j=-ncty,ncty
yy = j * ysiz
jc = ncty + 1 + j
do k=-nctz,nctz
zz = k * zsiz
kc = nctz + 1 + k
      call cova3(0.0,0.0,0.0,xx,yy,zz,nst,c0,it,cc,aa,vrotmat,cmax,covtab(ic,jc,kc))
      hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,vrotmat(:,:,1))
      if(real(hsqd).le.radsqd) then
            nlooku         = nlooku + 1
!
! We want to search by closest variogram distance (and use the
! anisotropic Euclidean distance to break ties:
!
            tmp(nlooku)   = - (covtab(ic,jc,kc)-TINY*real(hsqd))
            order(nlooku) = real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
      endif
end do
end do
end do
!
! Finished setting up the look-up table, now order the nodes such
! that the closest ones, according to variogram distance, are searched
! first. Note: the "loc" array is used because I didn't want to make
! special allowance for 2 byte integers in the sorting subroutine:
!
call sortem(1,nlooku,tmp,order)
do il=1,nlooku
      loc = int(order(il))
      iz  = int((loc-1)/MAXCXY) + 1
      iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
      ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
      iznode(il) = int(iz)
      iynode(il) = int(iy)
      ixnode(il) = int(ix)
end do

return
end subroutine



subroutine srchnd(ix,iy,iz, &
  sim, &
  nx,ny,nz, &
  xmn,ymn,zmn, &
  xsiz,ysiz,zsiz, &
  noct,nodmax, &
  nctx,ncty,nctz,&
  nlooku,ixnode,iynode,iznode, &
  ncnode,icnode,cnodex,cnodey,cnodez,cnodev)
implicit none
!arguments
integer(kind=gik), intent(in) :: ix,iy,iz
real(kind=fk), intent(in) ::sim(:)
integer(kind=gik), intent(in) ::ixnode(:),iynode(:),iznode(:)
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=gik), intent(in) :: nctx,ncty,nctz,noct,nodmax
integer(kind=gik), intent(in) ::nlooku
integer(kind=gik), intent(inout) ::ncnode,icnode(:)
real(kind=fk), intent(inout) ::cnodex(:),cnodey(:),cnodez(:),cnodev(:)
!locals
integer   ninoct(8),i,j,k,il,ind,idx,idy,idz,iq,nxy

nxy = nx*ny
!
! Consider all the nearby nodes until enough have been found:
!
ncnode = 0
if(noct.gt.0) then
      do i=1,8
            ninoct(i) = 0
      end do
end if
do 2 il=2,nlooku
      if(ncnode.eq.nodmax) return
      i = ix + (int(ixnode(il))-nctx-1)
      j = iy + (int(iynode(il))-ncty-1)
      k = iz + (int(iznode(il))-nctz-1)
      if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
      if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
      ind = i + (j-1)*nx + (k-1)*nxy
      if(sim(ind).gt.UNEST) then
!
! Check the number of data already taken from this octant:
!
            if(noct.gt.0) then
                  idx = ix - i
                  idy = iy - j
                  idz = iz - k
                  if(idz.gt.0) then
                        iq = 4
                        if(idx.le.0 .and. idy.gt.0) iq = 1
                        if(idx.gt.0 .and. idy.ge.0) iq = 2
                        if(idx.lt.0 .and. idy.le.0) iq = 3
                  else
                        iq = 8
                        if(idx.le.0 .and. idy.gt.0) iq = 5
                        if(idx.gt.0 .and. idy.ge.0) iq = 6
                        if(idx.lt.0 .and. idy.le.0) iq = 7
                  end if
                  ninoct(iq) = ninoct(iq) + 1
                  if(ninoct(iq).gt.noct) go to 2
            end if
            ncnode = ncnode + 1
            icnode(ncnode) = il
            cnodex(ncnode) = xmn + real(i-1)*xsiz
            cnodey(ncnode) = ymn + real(j-1)*ysiz
            cnodez(ncnode) = zmn + real(k-1)*zsiz
            cnodev(ncnode) = sim(ind)
      endif
2    continue
!
! Return to calling program:
!
return
end subroutine
end module