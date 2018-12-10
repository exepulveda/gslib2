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
  vmodel, & !variogram
  radsqd, & !search
  nodmax,&
  MAXCTX,MAXCTY,MAXCTZ,cbb)
use variography, only: vmodel_type
implicit none
!arguments
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
real(kind=fk), intent(in) :: radsqd
type(vmodel_type),intent(in) :: vmodel
integer(kind=gik), intent(in) ::nodmax
integer(kind=gik), intent(in) ::MAXCTX,MAXCTY,MAXCTZ
type(ctable_structure_type),intent(inout) :: ctable_st
real(kind=fk),intent(out) :: cbb
!locals
  integer(kind=gik) MAXXYZ,MAXCXY
  integer test

  MAXCXY = MAXCTX * MAXCTY
  MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
  allocate(ctable_st%covtab(MAXCTX,MAXCTY,MAXCTZ))
  allocate(ctable_st%icnode(nodmax),stat = test)
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
      vmodel, &
      radsqd,nodmax,&
      MAXCTX,MAXCTY,MAXCTZ,&
      ctable_st%covtab,ctable_st%nlooku,ctable_st%ixnode,ctable_st%iynode,ctable_st%iznode,cbb)
end subroutine

subroutine search_ctable(ix,iy,iz,sim, &
    nx,ny,nz, xmn,ymn,zmn, xsiz,ysiz,zsiz, &
    noct,nodmax, &
    ctable_st, &
    ncnode,icnode,cnodex,cnodey,cnodez,cnodev)
implicit none
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


  nctx = min(((ctable_st%MAXCTX-1)/2),(nx-1))
  ncty = min(((ctable_st%MAXCTY-1)/2),(ny-1))
  nctz = min(((ctable_st%MAXCTZ-1)/2),(nz-1))


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
    vr, sec,&    ! variable
    nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
    ktype,radius,sang1,sang2,sang3,sanis1,sanis2, & !search configuration
    ndmin,ndmax,noct,nodmax, & !simulated data search
    nst,c0,cc,aa,it,ang1,ang2,ang3,anis1,anis2, & !variogram
    seed,MAXSBX,MAXSBY,MAXSBZ,mxctx,mxcty,mxctz,nsim, &
    sim & !output
  )
use geometry, only: grid_type, setrot
use searching, only: sb_structure_type,setup_superblock,search_super_block
use variography, only: max_cova,setup_vmodel,vmodel_type,cova3_vmodel,print_vmodel
implicit none
!arguments
real(kind=fk), intent(inout) :: x(:),y(:),z(:)
real(kind=fk), intent(inout) :: vr(:),sec(:)
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=gik), intent(in) :: ktype
real(kind=fk), intent(in) :: radius,sang1,sang2,sang3,sanis1,sanis2
integer(kind=gik), intent(in) :: ndmin,ndmax,noct,nodmax
integer(kind=gik), intent(in) ::nst,it(:)
real(kind=fk), intent(in) :: c0,cc(:),aa(:),ang1(:),ang2(:),ang3(:),anis1(:),anis2(:)
integer(kind=gik), intent(in) :: seed,MAXSBX,MAXSBY,MAXSBZ,mxctx,mxcty,mxctz,nsim
real(kind=fk), intent(inout) :: sim(:)

!locals
type(vmodel_type) :: vmodel
type(sb_structure_type) :: sb_structure
type(ctable_structure_type) :: ctable_st
type(grid_type) :: sb_grid
integer(kind=gik) :: nmult = 0
integer(kind=gik) :: sstrat = 0
integer(kind=gik) :: nr,seed_size,max_int,rseed,i,j
real(kind=fk) :: realisation_seed(nsim),cbb

 !trivial checkings
 !1.- samples duplicated
 do i=1,size(x)
   do j=1,size(x)
     if (i/=j .and. x(i) == x(j) .and. y(i)==y(j) .and. z(i) == z(j)) then
      stop "DATA DUPLICATED"
      end if
   end do
  end do

  max_int = huge(seed)
  !with the given seed, generate seeds for each realisation
  call random_number(realisation_seed)


  sb_structure = setup_superblock(x,y,z,vr,sec,nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,sang1,sang2,sang3,sanis1,sanis2,radius*radius,MAXSBX,MAXSBY,MAXSBZ)

  call setup_vmodel(vmodel,nst,c0,it,cc,aa,ang1,ang2,ang3,anis1,anis2)

  !ctable
  call setup_ctable(ctable_st, &
    nx,ny,nz, &
    xmn,ymn,zmn, &
    xsiz,ysiz,zsiz, &
    vmodel, & !variogram
    radius*radius, & !search
    nodmax,&
    mxctx,mxcty,mxctz,cbb)

  do nr=1,nsim
    rseed = int(realisation_seed(nr)*max_int)
    print *,'realisation=',nr,'seed=',rseed
    call sgsim_realizarion_grid(nr, &
      x,y,z, & ! coordinates
      vr, &    ! variable
      nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
      ktype,radius, & !search configuration
      ndmin,ndmax,noct, & !
      nodmax,nmult, sstrat, & !simulated data search
      sb_structure, &
      vmodel,cbb, & !variogram
      ctable_st,&
      rseed,&
      sim) !output
  end do
end subroutine

subroutine sgsim_realizarion_grid( isim, &
    x,y,z, & ! coordinates
    vr, &    ! variable
    nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
    ktype,radius, & !search configuration
    ndmin,ndmax,noct, & !
    nodmax,nmult, sstrat, & !simulated data search
    sb_structure, &
    vmodel,cbb, & !variogram
    ctable_st,&
    seed,&
    sim) !output
use geometry
use searching, only: sb_structure_type,setup_superblock,search_super_block
use sorting, only: sortem_original
use gaussian_utils, only: gauinv
use variography, only: vmodel_type
implicit none
!arguments
integer(kind=gik), intent(in) :: isim
real(kind=fk), intent(in) :: x(:),y(:),z(:)
real(kind=fk), intent(in) :: vr(:)
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=gik), intent(in) :: ktype
real(kind=fk), intent(in) :: radius
integer(kind=gik), intent(in) :: ndmin,ndmax,noct,nodmax,nmult,sstrat
type(sb_structure_type),intent(in) :: sb_structure
type(vmodel_type),intent(in) :: vmodel
real(kind=fk), intent(in) :: cbb
type(ctable_structure_type),intent(in) :: ctable_st
integer(kind=gik), intent(in) :: seed
real(kind=fk), intent(inout) :: sim(:)
!locals
real(kind=fk)     var(10)
real(kind=dk)    p,cp,oldcp,w
logical   testind, trace
integer(kind=gik) nxyz,imult,ind,id,in
integer(kind=gik) ix,iy,iz,nnx,nny,nnz,nd,nxy
real(kind=fk)    xx,yy,zz,test,test2
real(kind=fk) order(nx*ny*nz)
real(kind=dk),allocatable :: A(:,:),b(:),s(:)
real(kind=fk) :: radsqd
real(kind=dk) :: xp
real(kind=fk)    cmean,cstdev,gmean
integer(kind=gik) id2,idbg,idum,index,irepo,jx,jy,jz,ldbg,ierr
integer(kind=gik) nclose,infoct(8),lktype
integer(kind=gik) close(size(x))
integer(kind=gik) ncnode,icnode(nodmax) !results of searching simulated nodes
real(kind=fk) :: cnodex(nodmax),cnodey(nodmax),cnodez(nodmax),cnodev(nodmax)
integer(kind=gik) :: ne,seed_size,i,j
real(kind=fk) :: av,ss,simval

radsqd = radius  * radius
idbg = 0
ldbg = 0
gmean = 0

nd = size(vr)
nxy = nx*ny
nxyz = nxy*nz

print *,"seed for realization:",seed,size(sim)
!set seed
seed_size=1
call random_seed(size=seed_size)
call random_seed(put=(/ seed /))

!
! Work out a random path for this realization:
!
      call random_number(sim)
      do ind=1,nxyz
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
      call sortem_original(1,nxyz,sim,order)
      print *,'FIRST PATH=',order(1)

#ifdef TRACE                  
      do ind=1,nxyz
        print *,'RANDOM PATH=',ind,sim(ind),order(ind)
      end do
#endif

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
    
    
#ifdef TRACE                  
    print *,in,index,"xx,yy,zz=",xx,yy,zz
    print *,"sstrat=",sstrat
    print *,"radsqd=",radsqd
    print *,"ndmax=",ndmax
    print *,"noct=",noct
    print *,"nclose=",nclose
#endif

!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously
! simulated grid nodes:
!
    if(sstrat.eq.0) then
      call search_super_block( &
          xx,yy,zz,radsqd, &
          sb_structure, &
          ndmax,noct, &
          x,y,z, &
          nclose,close,infoct)
          
      !check close are not duplicated
      do i=1,nclose
       do j=1,nclose
         if (i/=j .and. close(i) == close(j)) then
          stop "search_super_block DUPLICATED"
          end if
       end do
      end do


#ifdef TRACE                  
      print *,"ndmax=",ndmax
      print *,"noct=",noct
      print *,"nclose=",nclose
      print *,"close=",close(1:nclose)
      print *,"infoct=",infoct
#endif

      if(nclose.lt.ndmin) go to 5
      if(nclose.gt.ndmax) nclose = ndmax
    endif

    call search_ctable(ix,iy,iz,sim, &
        nx,ny,nz, xmn,ymn,zmn, xsiz,ysiz,zsiz, &
        noct,nodmax, &
        ctable_st, &
        ncnode,icnode,cnodex,cnodey,cnodez,cnodev)

#ifdef TRACE                  
    print *,"ncnode=",ncnode
    print *,"icnode=",icnode(1:ncnode)
    print *,"cnodex=",cnodex(1:ncnode)
    print *,"cnodey=",cnodey(1:ncnode)
    print *,"cnodez=",cnodez(1:ncnode)
    print *,"cnodev=",cnodev(1:ncnode)
    print *,"nclose,ncnode=",nclose,ncnode
#endif
        
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
      
      !hack
      !ncnode = 0
      call krige(x,y,z,vr,ix,iy,iz,xx,yy,zz,lktype,gmean,&
                 nclose,close, &
                 ncnode,icnode,cnodex,cnodey,cnodez,cnodev, &
                 vmodel,cbb,ctable_st,cmean,cstdev)
                 
      cmean = 0.0
      cstdev = 1.0
    endif
!
! Draw a random number and assign a value to this node:
!
    !p = acorni(idum)
    call random_number(p)
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

  ne = 0
  av = 0.0
  ss = 0.0
  do ind=1,nxyz
        simval = sim(ind)
        if(simval.gt.-9.0.and.simval.lt.9.0) then
              ne = ne + 1
              av = av + simval
              ss = ss + simval*simval
        end if
        !write(lout,'(g14.8)') simval
  end do

  av = av / max(real(ne),1.0)
  ss =(ss / max(real(ne),1.0)) - av * av
  write(ldbg,111) isim,ne,av,ss
  write(*,   111) isim,ne,av,ss
111        format(/,' Realization ',i3,': number   = ',i8,/, &
              '                  mean     = ',f12.4,       &
              ' (close to 0.0?)',/,                        &
              '                  variance = ',f12.4,       &
              ' (close to gammabar(V,V)? approx. 1.0)',/)

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

subroutine krige(x,y,z,vra, &
      ix,iy,iz,xx,yy,zz, &
      ktype, &
      gmean, &
      nclose,close,& !conditioning data
      ncnode,icnode,cnodex,cnodey,cnodez,cnodev, & !results of searching simulated nodes
      vmodel,cbb, & !variogram
      ctable_st, &
      cmean, &
      cstdev)
use variography, only: cova3_vmodel,vmodel_type
use geometry, only: sqdist
use solvers, only: system_solver
implicit none
!arguments
real(kind=fk), intent(in) :: x(:),y(:),z(:),vra(:)
integer(kind=gik), intent(in) :: ix,iy,iz
real(kind=fk), intent(in) :: xx,yy,zz
real(kind=fk), intent(in) :: gmean
integer(kind=gik), intent(in) :: ktype
type(vmodel_type), intent(in) :: vmodel
integer(kind=gik), intent(in) :: nclose,close(:)
integer(kind=gik) ncnode,icnode(:) !results of searching simulated nodes
real(kind=fk) :: cnodex(:),cnodey(:),cnodez(:),cnodev(:)
type(ctable_structure_type),intent(in) :: ctable_st
real(kind=fk), intent(in) :: cbb
real(kind=fk), intent(out) :: cmean,cstdev
!locals
logical first
integer na,neq,i,j,k,ierr,idbg,ie,is,ising,ldbg,indexi,indexj
real(kind=fk), allocatable :: a(:,:),b(:),s(:)
integer(kind=gik) :: lktype,ind
integer :: info,ss_option
real(kind=fk) :: cmax,cov,sumwts,x1,y1,z1,x2,y2,z2
!
! Size of the kriging system:
!
idbg = 0

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
      ! ind = ix + (iy-1)*nx + (iz-1)*nxy
      ! if(lvm(ind).le.-6.0.or.lvm(ind).ge.6.0) then
      !       lktype = 0
      !       go to 33
      ! end if
end if

allocate(a(neq,neq),b(neq),s(neq))

a = 0
b = 0
!
! Set up kriging matrices:
!
!samples part
do j=1,nclose
  indexj  = int(close(j))
  x2 = x(indexj)
  y2 = y(indexj)
  z2 = z(indexj)
  do i=j,nclose
    indexi  = int(close(i))
    x1 = x(indexi)
    y1 = y(indexi)
    z1 = z(indexi)
    if (indexi /=indexj .and. (x1 == x2 .and. y1 == y2 .and. z1 == z2)) then
    print *, indexi,indexj,x1,y1,z1,x2,y2,z2,close(i),close(j)
    stop "SAMPLES REPETEAD"
    end if
    call cova3_vmodel(x1,y1,z1,x2,y2,z2,vmodel,cmax,cov)
    a(i,j) = cov
    a(j,i) = cov

  end do

  call cova3_vmodel(xx,yy,zz,x2,y2,z2,vmodel,cmax,cov)
  b(j) = cov

end do

!simulated part
do j=nclose+1,na
  indexj  = j-nclose
  x2 = cnodex(indexj)
  y2 = cnodey(indexj)
  z2 = cnodez(indexj)
  do i=j,na
    indexi  = i-nclose
    x1 = cnodex(indexi)
    y1 = cnodey(indexi)
    z1 = cnodez(indexi)
    
    call cova3_vmodel(x1,y1,z1,x2,y2,z2,vmodel,cmax,cov)
#ifdef TRACE                  
    print *,x1,y1,z1,x2,y2,z2,cov
#endif
    a(i,j) = dble(cov)
    a(j,i) = dble(cov)

  end do
  call cova3_vmodel(xx,yy,zz,x2,y2,z2,vmodel,cmax,cov)
  
  !print *,'sim par',j,xx,yy,zz,x2,y2,z2,cov,sqdist(xx,yy,zz,x2,y2,z2,vrotmat(1,:,:))
  b(j) = dble(cov)
end do

!samples against simulated part
do j=nclose+1,na
  indexj  = j-nclose
  x2 = cnodex(indexj)
  y2 = cnodey(indexj)
  z2 = cnodez(indexj)
  do i=1,nclose
    indexi  = int(close(i))
    x1 = x(indexi)
    y1 = y(indexi)
    z1 = z(indexi)
    call cova3_vmodel(x1,y1,z1,x2,y2,z2,vmodel,cmax,cov)
    a(i,j) = dble(cov)
    a(j,i) = dble(cov)
  end do
end do


!
! Addition of OK constraint:
!
if(lktype.eq.1.or.lktype.eq.3) then
      do i=1,na
        a(i,na+1) = 1.0
        a(na+1,i) = 1.0
      end do
      a(na+1,na+1)    = 0.0
      b(na+1)  = 1.0
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


!debug krigin system
#ifdef TRACE_KS         
do i=1,neq
  write(*,*) 'TRACE_KS:',i,b(i),(a(i,j),j=1,neq)
end do
#endif

101        format('BEFORE SE: r(',i3,') =',f7.4,'  a= ',99f7.4)
!
! Solve the Kriging System:
!
if(neq.eq.1.and.lktype.ne.3) then
      b(1)  = b(1) / a(1,1)
      ising = 0
else
  if (lktype == 0) then
    ss_option = 1
  else
    ss_option = 1
  end if
  call system_solver(a,b,s,ss_option,info)
  if (info /= 0) then
    print *,"info=",info
    stop "BAD SYSTEM"
  end if
  !r must sum up 1
  sumwts = 0.0
  do i=1,na
    sumwts = sumwts + s(i)
  end do
  
  if (sumwts > 1.1 ) then
    print *, 'sumwts',sumwts
    stop "Weights must sum up 1.0"
  end if
endif

!
! Write a warning if the matrix is singular:
!
if(info.ne.0) then
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
#ifdef TRACE                  
write(ldbg,*) 'cbb: ',cbb
#endif

cmean  = 0.0
cstdev = cbb
sumwts = 0.0
do i=1,na
#ifdef TRACE                  
      write(ldbg,*) 'EST/VAR: ',i,s(i),vra(i),b(i)
#endif
      cmean  = cmean  + real(s(i))*vra(i)
      cstdev = cstdev - real(s(i)*b(i))
      sumwts = sumwts + real(s(i))
end do

if(lktype.eq.1) cstdev = cstdev - real(s(na+1))

if(lktype.eq.2) cmean  = cmean + gmean

if(lktype.eq.4) then
      ! ind    = ix + (iy-1)*nx + (iz-1)*nxy
      ! cmean  = cmean  + real(s(na+1))*lvm(ind)
      ! cstdev = cstdev - real(s(na+1) *rr(na+1))
end if
!
! Error message if negative variance:
!
if(cstdev.lt.0.0) then
      write(ldbg,*) 'ERROR: Negative Variance: ',cstdev,cbb
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
!      if(lktype.eq.4) write(ldbg,141) lvm(ind),s(na+1)
!141        format(' Sec Data  value ',f8.4,' weight ',f8.4)
      write(ldbg,142) gmean,cmean,cstdev
142        format(' Global mean ',f8.4,' conditional ',f8.4,' std dev ',f8.4)
end if
!
! Finished Here:
!
deallocate(a,s,b)

end subroutine

subroutine ctable( &
    nx,ny,nz, &
    xmn,ymn,zmn, &
    xsiz,ysiz,zsiz, &
    vmodel, &
    radsqd,nodmax,&
    MAXCTX,MAXCTY,MAXCTZ,&
    covtab,nlooku,ixnode,iynode,iznode,cbb)
use variography, only: cova3_vmodel,vmodel_type
use sorting, only: sortem_original
use geometry, only: sqdist
implicit none
!arguments
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
type(vmodel_type), intent(in) :: vmodel
real(kind=fk), intent(in) :: radsqd
real(kind=fk), intent(inout) :: covtab(:,:,:)
integer(kind=gik), intent(out) ::nlooku
integer(kind=gik), intent(inout) ::ixnode(:),iynode(:),iznode(:)
integer(kind=gik), intent(in) ::nodmax
integer(kind=gik), intent(in) ::MAXCTX,MAXCTY,MAXCTZ
real(kind=fk),intent(out) :: cbb
!locals
real,parameter :: TINY=1.0e-10
real(kind=fk)::tmp(MAXCTX * MAXCTY * MAXCTZ),order(MAXCTX * MAXCTY * MAXCTZ)
real(kind=fk)::cmax
integer nctx,ncty,nctz,i,j,k,MAXCXY,ic,jc,kc,il,ix,iy,iz,loc
real xx,yy,zz

real(kind=dk)    hsqd

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
call cova3_vmodel(0.0,0.0,0.0,0.0,0.0,0.0,vmodel,cmax,cbb)
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
      call cova3_vmodel(0.0,0.0,0.0,xx,yy,zz,vmodel,cmax,covtab(ic,jc,kc))
      hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,vmodel%structure(1)%rotmat)
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
call sortem_original(1,nlooku,tmp,order)
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
integer(kind=gik) :: ninoct(8),i,j,k,il,ind,idx,idy,idz,iq,nxy

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

#ifdef TRACE                  
print *, "srchnd::nlooku",nlooku
print *, "srchnd::(ix,iy,iz)",ix,iy,iz
print *, "srchnd::(nctx,ncty,nctz)",nctx,ncty,nctz
#endif

do 2 il=2,nlooku
      if(ncnode.eq.nodmax) return
      i = ix + (int(ixnode(il))-nctx-1)
      j = iy + (int(iynode(il))-ncty-1)
      k = iz + (int(iznode(il))-nctz-1)

#ifdef TRACE                  
      print *, "srchnd::do2.(il,ijk)",il,i,j,k
#endif


      if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
      if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
      ind = i + (j-1)*nx + (k-1)*nxy

#ifdef TRACE                  
      print *, "srchnd::do2.ind",ind,sim(ind)
      print *, "srchnd::ncnode",ncnode
#endif

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
