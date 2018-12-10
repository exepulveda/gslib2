module kriging
use kinds_kit

real,parameter :: TINY = 0.0001
real,parameter :: UNEST = -999
real(4), parameter :: EPSLON=1.e-5

contains

subroutine kriging_estimation( &
    koption, & !kriging option 0=grid, 1=cross, 2=jackknife
    xest,yest,zest, & !locations to estimate if koption /= 0
    x,y,z, & ! coordinates samples
    vr, ve, &    ! variable samples
    dh,extest, & !ddh, if koption==1, extest for drift
    nxdis,nydis,nzdis, & !block support
    nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
    ktype, &
    radius,sang1,sang2,sang3,sanis1,sanis2, & !search configuration
    ndmin,ndmax,noct, & !data search
    nst,c0,cc,aa,it,ang1,ang2,ang3,anis1,anis2, & !variogram
    itrend,idrif, &!drifts
    MAXSBX,MAXSBY,MAXSBZ, &
    idbg,ldbg, & !trace & debug
    est_val,est_var & !output
  )
use geometry, only: grid_type, setrot
use searching, only: sb_structure_type,setup_superblock,search_super_block
use variography, only: max_cova,setup_vmodel,vmodel_type,cova3_vmodel,print_vmodel
implicit none
!arguments
real(kind=fk), intent(in) :: xest(:),yest(:),zest(:)
real(kind=fk), intent(inout) :: x(:),y(:),z(:)
real(kind=fk), intent(inout) :: vr(:),ve(:),extest(:)
integer(kind=gik), intent(inout) :: dh(:)
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=gik), intent(in) :: ktype,koption,nxdis,nydis,nzdis
real(kind=fk), intent(in) :: radius,sang1,sang2,sang3,sanis1,sanis2
integer(kind=gik), intent(in) :: ndmin,ndmax,noct
integer(kind=gik), intent(in) ::nst,it(:)
real(kind=fk), intent(in) :: c0,cc(:),aa(:),ang1(:),ang2(:),ang3(:),anis1(:),anis2(:)
integer(kind=gik), intent(in) :: MAXSBX,MAXSBY,MAXSBZ
real(kind=fk), intent(inout) :: est_val(:),est_var(:)
integer(kind=gik), intent(in) :: idrif(:), idbg,ldbg,itrend

!locals
type(vmodel_type) :: vmodel
type(sb_structure_type) :: sb_structure
integer(kind=gik) :: nr,i,j,ix,iy,iz,index
real(kind=fk) :: cbb,radsqd,covmax,resc,cov,cmax
real(kind=fk) :: xdb(nxdis*nydis*nzdis),ydb(nxdis*nydis*nzdis),zdb(nxdis*nydis*nzdis)
real(kind=fk) :: xdis,ydis,zdis,xloc,yloc,zloc,resce
integer(kind=gik) :: nk,xk,vk,xkmae,xkmse,nsec,mdt,ndb,nxy,nxyz,nloop,irepo,nd,ddh
real(kind=fk) :: est,estv,bv(9),extest_iter!drifts

  if (koption /= 0) nd = size(xest)

  nxy = nx*ny
  nxyz = nxy*nz
  

!
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search.  Also compute the maximum covariance for
! the rescaling factor:
!
      write(*,*) 'Setting up rotation matrices for variogram and search'
      call setup_vmodel(vmodel,nst,c0,it,cc,aa,ang1,ang2,ang3,anis1,anis2)
      
      radsqd = radius * radius
      covmax = max_cova(vmodel)
!
! Finish computing the rescaling factor and stop if unacceptable:
!
      if(radsqd.lt.1.0) then
            resc = 2.0 * radius / max(covmax,0.0001)
      else
            resc =(4.0 * radsqd)/ max(covmax,0.0001)
      endif
      if(resc.le.0.0) then
            write(*,*) 'ERROR KT3D: The rescaling value is wrong ',resc
            write(*,*) '            Maximum covariance: ',covmax
            write(*,*) '            search radius:      ',radius
            stop
      endif
      resc = 1.0 / resc
!
! Set up for super block searching:
!
      write(*,*) 'Setting up super block search strategy'
      nsec = 2

      sb_structure = setup_superblock(x,y,z,vr,ve,nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,sang1,sang2,sang3,sanis1,sanis2,radsqd,MAXSBX,MAXSBY,MAXSBZ)


!
! Compute the number of drift terms, if an external drift is being
! considered then it is one more drift term, if SK is being considered
! then we will set all the drift terms off and mdt to 0):
!
      mdt = 1
      do i=1,9
!            if(ktype.eq.0.or.ktype.eq.2) idrif(i) = 0
!            if(idrif(i).lt.0.or.idrif(i).gt.1) then
!                  write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
!                  stop
!            endif
            mdt = mdt + idrif(i)
      end do
      if(ktype.eq.3) mdt = mdt + 1
      if(ktype.eq.0) mdt = 0
      if(ktype.eq.2) mdt = 0
!
! Set up the discretization points per block.  Figure out how many
! are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
! the offsets relative to the block center (this only gets done once):
!
! In all cases the offsets are relative to the lower left corner.
! This is done for rescaling the drift terms in the kriging matrix.
!
!      if(nxdis.lt.1) nxdis = 1
!      if(nydis.lt.1) nydis = 1
!      if(nzdis.lt.1) nzdis = 1
      ndb = max(nxdis,1) * max(nydis,1) * max(nzdis,1)
! MAXDIS is not necessary
!      if(ndb.gt.MAXDIS) then
!            write(*,*) 'ERROR KT3D: Too many discretization points',ndb
!            write(*,*) '            Increase MAXDIS or lower n[xyz]dis'
!            stop
!      endif
      xdis = xsiz  / max(real(nxdis),1.0)
      ydis = ysiz  / max(real(nydis),1.0)
      zdis = zsiz  / max(real(nzdis),1.0)
      i    = 0
      xloc = -0.5*(xsiz+xdis)
      do ix =1,max(nxdis,1)
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do iy=1,max(nydis,1)
                  yloc = yloc + ydis
                  zloc = -0.5*(zsiz+zdis)
                  do iz=1,max(nzdis,1)
                        zloc = zloc + zdis
                        i = i+1
                        xdb(i) = xloc + 0.5*xsiz
                        ydb(i) = yloc + 0.5*ysiz
                        zdb(i) = zloc + 0.5*zsiz
                  end do
            end do
      end do
!
! Initialize accumulators:
!
      nk    = 0
      xk    = 0.0
      vk    = 0.0
      xkmae = 0.0
      xkmse = 0.0
!
! Calculate Block Covariance. Check for point kriging.
!
      call cova3_vmodel(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),vmodel,cmax,cov)
!
! Set the ``unbias'' variable so that the matrix solution is more stable
!
!      unbias = cov
      cbb    = dble(cov)
      if(ndb.gt.1) then
            cbb = 0.0
            do i=1,ndb
               do j=1,ndb
                  call cova3_vmodel(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j),vmodel,cmax,cov)
                  if(i.eq.j) cov = cov - vmodel%nugget
                  cbb = cbb + dble(cov)
               end do
            end do
            cbb = cbb/dble(real(ndb*ndb))
      end if
      if(idbg.gt.1) then
            write(ldbg,*) ' '
            write(ldbg,*) 'Block Covariance: ',cbb
            write(ldbg,*) ' '
      end if
!
! Mean values of the drift functions:
!
      do i=1,9
            bv(i) = 0.0
      end do
      do i=1,ndb
            bv(1) = bv(1) + xdb(i)
            bv(2) = bv(2) + ydb(i)
            bv(3) = bv(3) + zdb(i)
            bv(4) = bv(4) + xdb(i)*xdb(i)
            bv(5) = bv(5) + ydb(i)*ydb(i)
            bv(6) = bv(6) + zdb(i)*zdb(i)
            bv(7) = bv(7) + xdb(i)*ydb(i)
            bv(8) = bv(8) + xdb(i)*zdb(i)
            bv(9) = bv(9) + ydb(i)*zdb(i)
      end do  
      do i=1,9
            bv(i) = (bv(i) / real(ndb)) * resc
      end do  
!
! Report on progress from time to time:
!
      if(koption.eq.0) then
            nxy   = nx*ny
            nxyz  = nxy*nz
            nloop = nxyz
            irepo = max(1,min((nxyz/10),10000))
      else
            nloop = 10000000
            irepo = max(1,min((nd/10),10000))
      end if
!      ddh = 0.0
!      write(*,*)
!      write(*,*) 'Working on the kriging '
!
! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
!

!$omp parallel do private(index,iz,iy,ix,xloc,yloc,zloc,ddh,extest_iter,est,estv)
  do index=1,nloop
      if((int(index/irepo)*irepo).eq.index) write(*,103) index
 103  format('   currently on estimate ',i9)
!
! Where are we making an estimate?
!
      if(koption.eq.0) then
            iz   = int((index-1)/nxy) + 1
            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
            ix   = index - (iz-1)*nxy - (iy-1)*nx
            xloc = xmn + real(ix-1)*xsiz
            yloc = ymn + real(iy-1)*ysiz
            zloc = zmn + real(iz-1)*zsiz
            extest_iter = 0
      else
            xloc = xest(index)
            yloc = yest(index)
            zloc = zest(index)
            if (size(dh) >= index) then
              ddh = dh(index)
            else
              ddh = 0
            end if
            
            if (size(extest) >= index) then
              extest_iter = extest(index)
            else
              extest_iter = 0
            end if
      end if

!
! Read in the external drift variable for this grid node if needed:
!
      if(ktype.eq.2.or.ktype.eq.3) then
            if(koption.eq.0) then
!                  read(lext,*) (var(i),i=1,iextve)
!                  extest = var(iextve)
            end if
!            if(extest(index).lt.tmin.or.extest(index).ge.tmax) then
!                  est  = UNEST
!                  estv = UNEST
!                  go to 1
!            end if
            resce  = covmax / max(extest_iter,0.0001)
      endif


    call krige(koption,xloc,yloc,zloc, ddh, & !location to krige
               x,y,z, & ! coordinates
               vr, ve, dh,&    ! variable
               nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
               ktype,radsqd,sb_structure, & !search configuration
               ndmin,ndmax,noct, & !data search
               vmodel, &
               xdb,ydb,zdb,cbb,& !discretization
               mdt, extest_iter, itrend, & !drift
               idbg,ldbg, &
               est,estv) !output
               
    est_val(index) = est
    est_var(index) = estv
!
! END OF MAIN KRIGING LOOP:
!
 1          continue
!            if(iktype.eq.0) then
!                  if(koption.eq.0) then
!                        write(lout,'(g14.8,1x,g14.8)') est,estv
!                  else
!                        err = UNEST
!                        if(koption.eq.1) then
!                        if(true(index).ne.UNEST.and.est.ne.UNEST) then
!                              err=est-true(index)
!                              xkmae = xkmae + abs(err)
!                              xkmse = xkmse + err*err
!                        end if
!                        write(lout,'(7(g14.8,1x))') xloc,yloc,zloc,true(index),est,estv,err
!                  end if
!            else
!
! Work out the IK-type distribution implicit to this data configuration
! and kriging weights:
!
!                  do icut=1,ncut
!                        cdf(icut) = -1.0
!                  end do
!                  wtmin = 1.0
!                  do i=1,na
!                        if(s(i).lt.wtmin) wtmin = s(i)
!                  end do
!                  sumwt = 0.0
!                  do i=1,na
!                        s(i)  = s(i) - wtmin
!                        sumwt = sumwt + s(i)
!                  end do
!                  do i=1,na
!                        s(i) = s(i) / max(0.00001,sumwt)
!                  end do
!                  if(na.gt.1.and.sumwt.gt.0.00001) then
!                        do icut=1,ncut
!                              cdf(icut) = 0.0
!                              do i=1,na
!                                    if(vra(i).le.cut(icut)) cdf(icut)=cdf(icut)+s(i)
!                              end do
!                        end do
!                  end if
!                  if(koption.eq.0) then
!                        write(lout,'(30(f8.4))') (cdf(i),i=1,ncut)
!                  else
!                        write(lout,'(30(f8.4))') (cdf(i),i=1,ncut),true
!                  end if
!            end if
  end do
!$omp end parallel do  
 2    continue
!     if(koption.gt.0) close(ljack)
!
! Write statistics of kriged values:
!
 
!      if(nk.gt.0.and.idbg.gt.0) then
!            xk    = xk/real(nk)
!            vk    = vk/real(nk) - xk*xk
!            xkmae = xkmae/real(nk)
!            xkmse = xkmse/real(nk)
!            write(ldbg,105) nk,xk,vk
!            write(*,   105) nk,xk,vk
! 105        format(/,'Estimated   ',i8,' blocks ',/,'  average   ',g14.8,/,'  variance  ',g14.8,/)
!            if(koption.ne.0) then
!                  write(*,106) xkmae,xkmse
! 106              format(/,'  mean error',g14.8,/,'  mean sqd e',g14.8)
!            end if
!      endif
!
! All finished the kriging:
!
      return
 96   stop 'ERROR in jackknife file!'
end subroutine


subroutine krige(koption,xloc,yloc,zloc,ddh, & !location to krige
    x,y,z, & ! coordinates
    vr, ve, dh, &    ! variable
    nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, & !grid configuration
    ktype,radsqd,sb_structure, & !search configuration
    ndmin,ndmax,noct, & !data search
    vmodel, &
    xdb,ydb,zdb,cbb, & !discretization
    mdt,extest,itrend, & !drift
    idbg,ldbg, & ! trace & debug
    est,estv & !output
  )
use geometry, only: grid_type, setrot
use searching, only: sb_structure_type,setup_superblock,search_super_block
use variography, only: vmodel_type,cova3_vmodel
use solvers, only: system_solver
implicit none
!arguments
integer(kind=gik), intent(in) :: koption,ddh,dh(:)
real(kind=fk), intent(in) :: xloc,yloc,zloc
real(kind=fk), intent(in) :: x(:),y(:),z(:)
real(kind=fk), intent(in) :: vr(:),ve(:)
integer(kind=gik), intent(in) :: nx,ny,nz
real(kind=fk), intent(in) :: xmn,ymn,zmn
real(kind=fk), intent(in) :: xsiz,ysiz,zsiz
integer(kind=gik), intent(in) :: ktype
real(kind=fk), intent(in) :: radsqd
type(sb_structure_type), intent(in) :: sb_structure
integer(kind=gik), intent(in) :: ndmin,ndmax,noct,mdt,ldbg,itrend,idbg
type(vmodel_type), intent(in) :: vmodel
real(kind=fk), intent(in) :: xdb(:),ydb(:),zdb(:),extest,cbb
!real(kind=fk), intent(in) :: xa(:),ya(:),za(:),vra(:),vea(:)
real(kind=fk), intent(inout) :: est,estv

!locals
integer(kind=gik) :: ix,iy,iz,i,j,ndb,nk,neq
real(kind=fk) :: xa(ndmax),ya(ndmax),za(ndmax),vra(ndmax),vea(ndmax)
real(kind=fk) :: dx,dy,dz,cmax,cb1,cb,cov,wt,skmean
integer(kind=gik) :: nclose,close(size(x)),infoct(8),na,ind
integer ss_info
logical accept,fircon,first
!equations
real(kind=fk),allocatable :: a(:,:),b(:),s(:)

  ndb = size(xdb)

  fircon = .true.
  est  = UNEST
  estv = UNEST

!
! Where are we making an estimate?
!
      if(koption.eq.0) then
!            iz   = int((index-1)/nxy) + 1
!            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
!            ix   = index - (iz-1)*nxy - (iy-1)*nx
!            xloc = xmn + real(ix-1)*xsiz
!            yloc = ymn + real(iy-1)*ysiz
!            zloc = zmn + real(iz-1)*zsiz
      else
!            read(ljack,*,err=96,end=2) (var(i),i=1,nvarij)
!            ddh  = 0.0
!            xloc = xmn
!            yloc = ymn
!            zloc = zmn
!            true = UNEST
!            secj = UNEST
!            if(idhlj.gt.0)  ddh    = var(idhlj)
!            if(ixlj.gt.0)   xloc   = var(ixlj)
!            if(iylj.gt.0)   yloc   = var(iylj)
!            if(izlj.gt.0)   zloc   = var(izlj)
!            if(ivrlj.gt.0)  true   = var(ivrlj)
!            if(iextvj.gt.0) extest = var(iextvj)
!            if(true.lt.tmin.or.true.ge.tmax) true = UNEST
      end if

!
! Read in the external drift variable for this grid node if needed:
!
      if(ktype.eq.2.or.ktype.eq.3) then
!            if(koption.eq.0) then
!                  read(lext,*) (var(i),i=1,iextve)
!                  extest = var(iextve)
!           end if
!            if(extest.lt.tmin.or.extest.ge.tmax) then
!                  est  = UNEST
!                  estv = UNEST
!                  go to 1
!            end if
!            resce  = covmax / max(extest,0.0001)
      endif
!
! Find the nearest samples:
!
      call search_super_block( &
          xloc,yloc,zloc,radsqd, &
          sb_structure, &
          ndmax,noct, &
          x,y,z, &
          nclose,close,infoct)

!
! Load the nearest data in xa,ya,za,vra,vea:
!
      na = 0
      do i=1,nclose
            ind    = close(i)
            accept = .true.
            if(koption.ne.0.and.(abs(x(ind)-xloc)+abs(y(ind)-yloc)+ abs(z(ind)-zloc)).lt.EPSLON) accept = .false.
            if(koption.ne.0.and.(abs(dh(ind)-ddh)).lt.EPSLON) accept = .false.
            if(accept) then
                  if(na.lt.ndmax) then
                        na = na + 1
                        xa(na)  = x(ind) - xloc + 0.5*xsiz
                        ya(na)  = y(ind) - yloc + 0.5*ysiz
                        za(na)  = z(ind) - zloc + 0.5*zsiz
                        vra(na) = vr(ind)
                        if (size(ve) > 0) vea(na) = ve(ind)
                  end if
            end if
      end do
!
! Test number of samples found:
!
      if(na.lt.ndmin) then
            est  = UNEST
            estv = UNEST
            return
      end if
!
! Test if there are enough samples to estimate all drift terms:
!
      if(na.ge.1.and.na.le.mdt) then
            if(fircon) then
                  write(ldbg,999)
                  fircon = .false.
            end if
            est  = UNEST
            estv = UNEST
            return
      end if
 999  format(' Encountered a location where there were too few data ',/, &
             ' to estimate all of the drift terms but there would be',/, &
             ' enough data for OK or SK.   KT3D currently leaves ',/, &
             ' these locations unestimated.',/, &
             ' This message is only written once - the first time.',/)
!
! There are enough samples - proceed with estimation.
!
      if(na.le.1) then
!
! Handle the situation of only one sample:
!
            call cova3_vmodel(xa(1),ya(1),za(1),xa(1),ya(1),za(1),vmodel,cmax,cb1)
!
! Establish Right Hand Side Covariance:
!
            if(ndb.le.1) then
                  call cova3_vmodel(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),vmodel,cmax,cb)
            else
                  cb  = 0.0
                  do i=1,ndb
                        call cova3_vmodel(xa(1),ya(1),za(1),xdb(i),ydb(i),zdb(i),vmodel,cmax,cov)
                        cb = cb + cov
                        dx = xa(1) - xdb(i)
                        dy = ya(1) - ydb(i)
                        dz = za(1) - zdb(i)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-vmodel%nugget
                  end do
                  cb = cb / real(ndb)
            end if
!vd
!
! Early bug - always did OK in presence of one data.
!
!vd
            if(ktype.eq.2) skmean = extest
            if(ktype.eq.0.or.ktype.eq.2) then
                  wt   = cb / cb1
                  est  = wt * vra(1) + (1.0-wt) * skmean
                  estv = real(cbb) - wt*cb
            else
                  est  = vra(1)
                  estv = real(cbb) - 2.0*cb + cb1
            end if
            return
      end if
!
! Go ahead and set up the OK portion of the kriging matrix:
!
      neq = mdt+na
      
      allocate(a(neq,neq),b(neq),s(neq))      
      
!
! Initialize the main kriging matrix:
!

      first = .false.
      a = 0.0
      b = 0.0
!
! Fill in the kriging matrix:
!
      do i=1,na
      do j=i,na
            call cova3_vmodel(xa(i),ya(i),za(i),xa(j),ya(j),za(j),vmodel,cmax,cov)
            a(i,j) = cov
            a(j,i) = cov
      end do
      end do
!
! Fill in the OK unbiasedness portion of the matrix (if not doing SK):
!
      if(neq.gt.na) then
            do i=1,na
                  a(na+1,i) = 1.0
                  a(i,na+1) = 1.0
            end do
      endif
!
! Set up the right hand side:
!
      do i=1,na
            if(ndb.le.1) then
                  call cova3_vmodel(xa(i),ya(i),za(i),xdb(1),ydb(1),zdb(1),vmodel,cmax,cb)
            else
                  cb  = 0.0
                  do j=1,ndb
                        call cova3_vmodel(xa(i),ya(i),za(i),xdb(j),ydb(j),zdb(j),vmodel,cmax,cov)
                        cb = cb + cov
                        dx = xa(i) - xdb(j)
                        dy = ya(i) - ydb(j)
                        dz = za(i) - zdb(j)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-vmodel%nugget
                  end do
                  cb = cb / real(ndb)
            end if
            b(i) = cb
      end do
      if(neq.gt.na) b(na+1) = 1.0
!COMMENT EXE
!!
!! Add the additional unbiasedness constraints:
!!
!      im = na + 1
!!
!! First drift term (linear in "x"):
!!
!      if(idrif(1).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(xa(k)*resc)
!                  a(neq*(k-1)+im) = dble(xa(k)*resc)
!            end do
!            r(im) = dble(bv(1))
!      endif
!!
!! Second drift term (linear in "y"):
!!
!      if(idrif(2).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(ya(k)*resc)
!                  a(neq*(k-1)+im) = dble(ya(k)*resc)
!            end do
!            r(im) = dble(bv(2))
!      endif
!!
!! Third drift term (linear in "z"):
!!
!      if(idrif(3).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(za(k)*resc)
!                  a(neq*(k-1)+im) = dble(za(k)*resc)
!            end do
!            r(im) = dble(bv(3))
!      endif
!!
!! Fourth drift term (quadratic in "x"):
!!
!      if(idrif(4).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(xa(k)*xa(k)*resc)
!                  a(neq*(k-1)+im) = dble(xa(k)*xa(k)*resc)
!            end do
!            r(im) = dble(bv(4))
!      endif
!!
!! Fifth drift term (quadratic in "y"):
!!
!      if(idrif(5).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(ya(k)*ya(k)*resc)
!                  a(neq*(k-1)+im) = dble(ya(k)*ya(k)*resc)
!            end do
!            r(im) = dble(bv(5))
!      endif
!!
!! Sixth drift term (quadratic in "z"):
!!
!      if(idrif(6).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(za(k)*za(k)*resc)
!                  a(neq*(k-1)+im) = dble(za(k)*za(k)*resc)
!            end do
!            r(im) = dble(bv(6))
!      endif
!!
!! Seventh drift term (quadratic in "xy"):
!!
!      if(idrif(7).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(xa(k)*ya(k)*resc)
!                  a(neq*(k-1)+im) = dble(xa(k)*ya(k)*resc)
!            end do
!            r(im) = dble(bv(7))
!      endif
!!
!! Eighth drift term (quadratic in "xz"):
!!
!      if(idrif(8).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(xa(k)*za(k)*resc)
!                  a(neq*(k-1)+im) = dble(xa(k)*za(k)*resc)
!            end do
!            r(im) = dble(bv(8))
!      endif
!!
!! Ninth drift term (quadratic in "yz"):
!!
!      if(idrif(9).eq.1) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(ya(k)*za(k)*resc)
!                  a(neq*(k-1)+im) = dble(ya(k)*za(k)*resc)
!            end do
!            r(im) = dble(bv(9))
!      endif
!!
!! External drift term (specified by external variable):
!!
!      if(ktype.eq.3) then
!            im=im+1
!            do k=1,na
!                  a(neq*(im-1)+k) = dble(vea(k)*resce)
!                  a(neq*(k-1)+im) = dble(vea(k)*resce)
!            end do
!            r(im) = dble(extest*resce)
!      endif
!
! Copy the right hand side to compute the kriging variance later:
!
!      do k=1,neq
!            rr(k) = r(k)
!      end do
!      kadim = neq * neq
!      ksdim = neq
!      nrhs  = 1
!      nv    = 1
!
! If estimating the trend then reset all the right hand side terms=0.0:
!
      if(itrend.ge.1) then
            do i=1,na
                  b(i)  = 0.0
                  !rr(i) = 0.0
            end do
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
      if(idbg.eq.3) then
            write(ldbg,*) 'Estimating node index : ',xloc,yloc,zloc
            do i=1,neq
                  write(ldbg,100) i,b(i),(a(i,j),j=1,neq)
 100              format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
            end do
      endif
!
! Solve the kriging system:
!
      !call ktsol(neq,nrhs,nv,a,r,s,ising,maxeq)
      call system_solver(a,b,s,1,ss_info)
!
! Compute the solution:
!
      if(ss_info.ne.0) then
            if(idbg.ge.3) write(ldbg,*) ' Singular Matrix ',xloc,yloc,zloc
            est  = UNEST
            estv = UNEST
      else
            est  = 0.0
            estv = real(cbb)
            if(ktype.eq.2) skmean = extest
            do j=1,neq
                  estv = estv - real(s(j))*b(j)
                  if(j.le.na) then
                        if(ktype.eq.0) then
                              est = est + real(s(j))*(vra(j)-skmean)
                        else if(ktype.eq.2) then
                              est = est + real(s(j))*(vra(j)-vea(j))
                        else
                              est = est + real(s(j))*vra(j)
                        endif
                  endif
            end do
            if(ktype.eq.0.or.ktype.eq.2) est = est + skmean
!
! Write the kriging weights and data if debugging level is above 2:
!
            if(idbg.ge.2) then
                  write(ldbg,*) '       '
                  write(ldbg,*) 'BLOCK: ',ix,iy,iz,' at ',xloc,yloc,zloc
                  write(ldbg,*) '       '
                  if(ktype.ne.0) write(ldbg,*) '  Lagrange : ',s(na+1)
                  write(ldbg,*) '  BLOCK EST: x,y,z,vr,wt '
                  do i=1,na
                        xa(i) = xa(i) + xloc - 0.5*xsiz
                        ya(i) = ya(i) + yloc - 0.5*ysiz
                        za(i) = za(i) + zloc - 0.5*zsiz
                        write(ldbg,'(5f12.3)') xa(i),ya(i),za(i),vra(i),s(i)
                  end do
                  write(ldbg,*) '  estimate, variance  ',est,estv
            endif
      endif
end subroutine


end module
