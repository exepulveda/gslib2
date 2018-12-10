C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 2003, Statios Software and Services Incorporated.  All %
C rights reserved.                                                     %
C                                                                      %
C This program has been modified from the one distributed in 1996 (see %
C below).  This version is also distributed in the hope that it will   %
C be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
C code may be redistributed without restriction; however, this code is %
C for one developer only. Each developer or user of this source code   %
C must purchase a separate copy from Statios.                          %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
C Junior University.  All rights reserved.                             %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Module to declare dynamic arrays in multiple subroutines:
c
      module geostat
      
      integer,allocatable :: nisb(:),ixsbtosr(:),iysbtosr(:),izsbtosr(:)
      real,allocatable    :: x(:),y(:),z(:),vr(:),ve(:),dh(:),tmp(:),
     +          close(:),xa(:),ya(:),za(:),vra(:),vea(:),xdb(:),ydb(:),
     +          zdb(:),cut(:),cdf(:)
      real*8,allocatable  :: r(:),rr(:),s(:),a(:)
      
      end module
c
c
c
      program main
c-----------------------------------------------------------------------
c
c             Kriging (SK,OK,KT) of a 3-D Rectangular Grid
c             ********************************************
c
c The program is executed with no command line arguments.  The user
c will be prompted for the name of a parameter file.  The parameter
c file is described in the documentation (see the example kt3d.par)
c and should contain the following information:
c
c
c
c AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
c-----------------------------------------------------------------------
      use       geostat
      use       omp_lib
      include  'kt3d.inc'
      
      IF (OMP_GET_THREAD_NUM() .EQ. 0) THEN
        PRINT *, 'Number of threads = ', omp_get_max_threads()
        PRINT *, 'Number of processors = ', omp_get_num_procs()
      END IF      
      
c
c Read the parameters, the data, and open the output files:
c
      call readparm(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
c
c Call kt3d to krige the grid:
c
      call kt3d(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
c
c Finished:
c
      close(ldbg)
      close(lout)
      write(*,9998) VERSION
 9998 format(/' KT3D Version: ',f5.3, ' Finished'/)
      stop
      end
 
 
 
      subroutine readparm(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
c-----------------------------------------------------------------------
c
c                  Initialization and Read Parameters
c                  **********************************
c
c The input parameters and data are read in from their files. Some quick
c error checking is performed and the statistics of all the variables
c being considered are written to standard output.
c
c
c
c-----------------------------------------------------------------------
      use       geostat
      include  'kt3d.inc'
      parameter(MV=100)
      real      var(MV)
      character datafl*512,jackfl*512,extfl*512,outfl*512,dbgfl*512,
     +          str*512,title*80
      logical   testfl
c
c FORTRAN Units:
c
      lin   = 1
      ldbg  = 3
      lout  = 4
      lext  = 7
      ljack = 8
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' KT3D Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'kt3d.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'kt3d.par            ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) idhl,ixl,iyl,izl,ivrl,iextv
      write(*,*) ' columns = ',idhl,ixl,iyl,izl,ivrl,iextv

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,*,err=98) koption
      write(*,*) ' kriging option = ',koption

c
c This is an undocumented feature to have kt3d construct an IK-type
c distribution:
c
      iktype = 0
      if(koption.lt.0) then
            iktype  = 1
            koption = -koption
      end if
      if(iktype.eq.1) then

            read(lin,*,err=98) ncut
            write(*,*) ' number of cutoffs = ',ncut
c
c Find the needed parameter:
c
            MAXCUT = ncut
c
c Allocate the needed memory:
c21
            allocate(cut(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to',
     +                        ' insufficient memory.'
                        stop
                  end if
c22
            allocate(cdf(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to',
     +                        ' insufficient memory.'
                        stop
                  end if
c
            read(lin,*,err=98) (cut(i),i=1,ncut)
            write(*,*) ' cutoffs = ',(cut(i),i=1,ncut)

      end if

      read(lin,'(a512)',err=98) jackfl
      call chknam(jackfl,512)

      write(*,*) ' jackknife data file = ',jackfl(1:40)

      read(lin,*,err=98) ixlj,iylj,izlj,ivrlj,iextvj
      write(*,*) ' columns = ',ixlj,iylj,izlj,ivrlj,iextvj

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      
      write(*,*) ' debugging file = ',dbgfl(1:40)

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' nx, xmn, xsiz = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' ny, ymn, ysiz = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' nz, zmn, zsiz = ',nz,zmn,zsiz

      read(lin,*,err=98) nxdis,nydis,nzdis
      write(*,*) ' block discretization:',nxdis,nydis,nzdis

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' ndmin,ndmax = ',ndmin,ndmax

      read(lin,*,err=98) noct
      write(*,*) ' max per octant = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) ktype,skmean
      write(*,*) ' ktype, skmean =',ktype,skmean

      read(lin,*,err=98) (idrif(i),i=1,9)
      write(*,*) ' drift terms = ',(idrif(i),i=1,9)

      read(lin,*,err=98) itrend
      write(*,*) ' itrend = ',itrend

      read(lin,'(a512)',err=98) extfl
      call chknam(extfl,40)

      write(*,*) ' external drift file = ',extfl(1:40)

      read(lin,*,err=98) iextve
      write(*,*) ' variable in external drift file = ',iextve

      read(lin,*,err=98) nst(1),c0(1)
      write(*,*) ' nst, c0 = ',nst(1),c0(1)

      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/,
     +             ' The c or a values can be set to zero')
            stop
      endif

      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)
            read(lin,*,err=98) aa(i),aa1,aa2
            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i),
     +                   ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
      end do

      close(lin)
c
c Find the needed parameters:
c
      MAXDIS = nxdis*nydis*nzdis
      MAXSAM = ndmax + 1
      MAXEQ = MAXSAM + MAXDT + 2
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2.00)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
c
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2.00)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
c
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2.00)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
c
      MAXSB = MAXSBX*MAXSBY*MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX  = 2 * MAXSBX
c
c Allocate the needed memory:
c1
      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c2
      allocate(ixsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c3
      allocate(iysbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c4
      allocate(izsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c13
      allocate(xa(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c14
      allocate(ya(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c15
      allocate(za(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c16
      allocate(vra(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c17
      allocate(vea(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c18
      allocate(xdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c19
      allocate(ydb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c20
      allocate(zdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c23
      allocate(r(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c24
      allocate(rr(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c25
      allocate(s(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c26
      allocate(a(MAXEQ * MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c
c Perform some quick error checking:
c
      if(ndmax.gt.MAXSAM) stop 'ndmax is too big - modify .inc file'
      if(ktype.eq.3.and.iextv.le.0) stop 'must have external variable'
      if(ixl.le.0.and.nx.gt.1) write(*,*) ' WARNING: ixl=0 and nx>1 ! '
      if(iyl.le.0.and.ny.gt.1) write(*,*) ' WARNING: iyl=0 and ny>1 ! '
      if(izl.le.0.and.nz.gt.1) write(*,*) ' WARNING: izl=0 and nz>1 ! '
c
c Check to make sure the data file exists, then either read in the
c data or write an error message and stop:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR data file ',datafl,' does not exist!'
            stop
      endif
c
c The data file exists so open the file and read in the header
c information. Initialize the storage that will be used to summarize
c the data found in the file:
c
      title(1:22) = 'KT3D ESTIMATES WITH: '
      open(lin,file=datafl,status='OLD')
      read(lin,*)
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,*)
      end do
      MAXDAT = 0
 22   read(lin,*,end=33,err=99) (var(j),j=1,nvari)
      if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) go to 22
      MAXDAT = MAXDAT + 1
      go to 22
 33   continue
c
c Allocate the needed memory:
c5
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c6
      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c7
      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c8
      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c9
      allocate(ve(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c10
      allocate(dh(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c11
      allocate(tmp(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c12
      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to',
     +                  ' insufficient memory.'
                  stop
            end if
c
      rewind(lin)
      read(lin,'(a58)') title(23:80)
      read(lin,*,err=99)       nvari
      nd = 0
      av = 0.0
      ss = 0.0
      do i=1,nvari
            read(lin,'(a40)',err=99) str
      end do
c
c Some tests on column numbers:
c
      if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.ivrl.gt.nvari)
     +      then
            write(*,*) 'There are only ',nvari,' columns in input data'
            write(*,*) '  your specification is out of range'
            stop
      end if
c
c Read all the data until the end of the file:
c
 2    read(lin,*,end=3,err=99) (var(j),j=1,nvari)
      vrt = var(ivrl)
      if(vrt.lt.tmin.or.vrt.ge.tmax) go to 2
      nd = nd + 1
      if(nd.gt.MAXDAT) then
            write(*,*) ' ERROR: Exceeded available memory for data'
            stop
      end if
c
c Establish the location of this datum:
c
      if(idhl.le.0) then
            dh(nd) = -99
      else
            dh(nd) = var(idhl)
      endif
      if(ixl.le.0) then
            x(nd) = xmn
      else
            x(nd) = var(ixl)
      endif
      if(iyl.le.0) then
            y(nd) = ymn
      else
            y(nd) = var(iyl)
      endif
      if(izl.le.0) then
            z(nd) = zmn
      else
            z(nd) = var(izl)
      endif
c
c Establish the external drift variable (if needed):
c
      ve(nd) = 1.0
      if(ktype.eq.3.or.ktype.eq.2) then
            ve(nd) = var(iextv)
            if(ve(nd).lt.tmin.or.ve(nd).ge.tmax) then
                  write(*,*) ' External drift variable must be present',
     +                       ' at all data locations!'
                  write(*,*) ' Encountered at data number ',nd
                  stop
            end if
      end if
      vr(nd) = vrt
      av     = av + vrt
      ss     = ss + vrt*vrt
      go to 2
 3    close(lin)
c
c Compute the averages and variances as an error check for the user:
c
      av = av / max(real(nd),1.0)
      ss =(ss / max(real(nd),1.0)) - av * av
      write(*,*) 'Data for KT3D: Variable number ',ivrl
      write(*,*) '  Number   = ',nd
      write(*,*) '  Average  = ',av
      write(*,*) '  Variance = ',ss
      if(nd.lt.1) then
            write(*,*) ' ERROR: there are no data'
            stop
      end if
c
c Open the debugging and output files:
c
      open(ldbg,file=dbgfl,status='UNKNOWN')
      open(lout,file=outfl,status='UNKNOWN')
      write(lout,'(a80)') title

      if(iktype.eq.0.and.koption.eq.0) then
           write(lout,201) 2,nx,ny,nz
           write(lout,102)
 102       format('Estimate',/,'EstimationVariance')
      end if
      if(iktype.eq.0.and.koption.ge.1) then
           write(lout,201) 7
           write(lout,103)
 103       format('X',/,'Y',/,'Z',/,'True',/,'Estimate',/,
     +            'EstimationVariance',/,'Error: est-true')
      end if
 201  format(4(1x,i4))

      if(iktype.eq.1) then
            if(koption.eq.0) then
                  write(lout,201) ncut,nx,ny,nz
            else
                  write(lout,201) ncut+1
            end if
            do i=1,ncut
                  write(lout,104) i,cut(i)
 104              format('Threshold: ',i2,' = ',f12.5)
            end do
            if(koption.eq.1) write(lout,105)
 105        format('true value')
      end if
c
c Open the external drift file if needed and position it at the
c first grid node in the file:
c
      if((ktype.eq.2.or.ktype.eq.3).and.koption.eq.0) then
            inquire(file=extfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR file ',extfl,' does not exist!'
                  stop
            endif
            open(lext,file=extfl,status='UNKNOWN')
            read(lext,'(a40)',err=97) str
            read(lext,*,err=97)       nvari
            do i=1,nvari
                  read(lext,'(a40)',err=97) str
            end do
            if(idbg.ge.3) write(ldbg,100) iextve
 100        format('A secondary variable is being used.  The gridded '
     +             'file',/,'must have the same grid specifications '
     +             'as the grid you are kriging.',/,'The external '
     +             'drift variable was taken from column ',i2)
      endif
c
c Set up for cross validation:
c
      if(koption.eq.1) then
            jackfl = datafl
            idhlj  = idhl
            ixlj   = ixl
            iylj   = iyl
            izlj   = izl
            ivrlj  = ivrl
            iextvj = iextv
      end if
c
c Open the file with the jackknife data?
c
      if(koption.gt.0) then
            inquire(file=jackfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR file ',jackfl,' does not exist!'
                  stop
            endif
            open(ljack,file=jackfl,status='OLD')
            read(ljack,*,err=96)
            read(ljack,*,err=96) nvarij
            do i=1,nvarij
                  read(ljack,*,err=96)
            end do
      end if
c
c Finished here:
c
      return
c
c Error in an Input File Somewhere:
c
 96   stop 'ERROR in jackknife file!'
 97   stop 'ERROR in external drift file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      subroutine kt3d(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
c-----------------------------------------------------------------------
c
c                Krige a 3-D Grid of Rectangular Blocks
c                **************************************
c
c This subroutine estimates point or block values of one variable by
c simple, ordinary, or kriging with a trend model.  It is also possible
c to estimate the trend directly.
c
c
c
c PROGRAM NOTES:
c
c   1. The data and parameters are passed in common blocks defined
c      in kt3d.inc.  Local storage is allocated in the subroutine
c      for kriging matrices, i.e.,
c         - xa,ya,za,vra   arrays for data within search neighborhood
c         - a,r,rr,s       kriging arrays
c         - xdb,ydb,zdb    relative position of discretization points
c         - cbb            block covariance
c   2. The kriged value and the kriging variance is written to Fortran
c      unit number "lout".
c
c
c
c
c Original:  A.G. Journel and C. Lemmer                             1981
c Revisions: A.G. Journel and C. Kostov                             1984
c-----------------------------------------------------------------------
      use        geostat
      use        kriging, only: kriging_estimation
      use        kinds_kit
      include   'kt3d.inc'
      real*8     cbb
      real       var(0)
      logical    first,fircon,accept
      data       fircon/.true./
      real(fk)   xest(0),yest(0),zest(0),extest(0)
      real(fk),allocatable :: est_var(:),est_val(:)
      integer(gik)  dh_tmp(0),true(0),nest,i
c
c Set up the rotation/anisotropy matrices that are needed for the
c variogram and search.  Also compute the maximum covariance for
c the rescaling factor:
c
      if (koption == 0) then
        nest = nx*ny*nz
      else
        nest = nd
      end if
      
      allocate(est_var(nest),est_val(nest))

      call kriging_estimation(
     +   koption, 
     +   xest,yest,zest,
     +   x,y,z, 
     +   vr, ve,
     +   dh_tmp,extest, 
     +   nxdis,nydis,nzdis, 
     +   nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz, 
     +   ktype, 
     +   radius,sang1,sang2,sang3,sanis1,sanis2, 
     +   ndmin,ndmax,noct, 
     +   nst(1),c0(1),cc,aa,it,ang1,ang2,ang3,anis1,anis2,
     +   itrend,idrif, 
     +   MAXSBX,MAXSBY,MAXSBZ,
     +   idbg,ldbg, 
     +   est_val,est_var)

      do i=1,nest
         write(lout,'(g14.8,1x,g14.8)') est_val(i),est_var(i)
      end do
      
 2    continue
c      if(koption.gt.0) close(ljack)
c
c Write statistics of kriged values:
c
 
      if(nk.gt.0.and.idbg.gt.0) then
            xk    = xk/real(nk)
            vk    = vk/real(nk) - xk*xk
            xkmae = xkmae/real(nk)
            xkmse = xkmse/real(nk)
            write(ldbg,105) nk,xk,vk
            write(*,   105) nk,xk,vk
 105        format(/,'Estimated   ',i8,' blocks ',/,
     +               '  average   ',g14.8,/,'  variance  ',g14.8,/)
            if(koption.ne.0) then
                  write(*,106) xkmae,xkmse
 106              format(/,'  mean error',g14.8,/,'  mean sqd e',g14.8)
            end if
      endif
c
c All finished the kriging:
c
      return
 96   stop 'ERROR in jackknife file!'
      end



      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='kt3d.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for KT3D',/,
     +       '                  *******************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat              ',
     +       '-file with data')
      write(lun,12)
 12   format('0  1  2  0  3  0                 ',
     +       '-   columns for DH,X,Y,Z,var,sec var')
      write(lun,13)
 13   format('-1.0e21   1.0e21                 ',
     +       '-   trimming limits')
      write(lun,14)
 14   format('0                                ',
     +       '-option: 0=grid, 1=cross, 2=jackknife')
      write(lun,15)
 15   format('xvk.dat                          ',
     +       '-file with jackknife data')
      write(lun,16)
 16   format('1   2   0    3    0              ',
     +       '-   columns for X,Y,Z,vr and sec var')
      write(lun,17)
 17   format('3                                ',
     +       '-debugging level: 0,1,2,3')
      write(lun,18)
 18   format('kt3d.dbg                         ',
     +       '-file for debugging output')
      write(lun,19)
 19   format('kt3d.out                         ',
     +       '-file for kriged output')
      write(lun,20)
 20   format('50   0.5    1.0                  ',
     +       '-nx,xmn,xsiz')
      write(lun,21)
 21   format('50   0.5    1.0                  ',
     +       '-ny,ymn,ysiz')
      write(lun,22)
 22   format('1    0.5    1.0                  ',
     +       '-nz,zmn,zsiz')
      write(lun,23)
 23   format('1    1      1                    ',
     +       '-x,y and z block discretization')
      write(lun,24)
 24   format('4    8                           ',
     +       '-min, max data for kriging')
      write(lun,25)
 25   format('0                                ',
     +       '-max per octant (0-> not used)')
      write(lun,26)
 26   format('20.0  20.0  20.0                 ',
     +       '-maximum search radii')
      write(lun,27)
 27   format(' 0.0   0.0   0.0                 ',
     +       '-angles for search ellipsoid')
      write(lun,28)
 28   format('0     2.302                      ',
     +       '-0=SK,1=OK,2=non-st SK,3=exdrift')
      write(lun,29)
 29   format('0 0 0 0 0 0 0 0 0                ',
     +       '-drift: x,y,z,xx,yy,zz,xy,xz,zy')
      write(lun,30)
 30   format('0                                ',
     +       '-0, variable; 1, estimate trend')
      write(lun,31)
 31   format('extdrift.dat                     ',
     +       '-gridded file with drift/mean')
      write(lun,32)
 32   format('4                                ',
     +       '-  column number in gridded file')
      write(lun,33)
 33   format('1    0.2                         ',
     +       '-nst, nugget effect')
      write(lun,34)
 34   format('1    0.8  0.0   0.0   0.0        ',
     +       '-it,cc,ang1,ang2,ang3')
      write(lun,35)
 35   format('         10.0  10.0  10.0        ',
     +       '-a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end

      subroutine chknam(str,len)
c-----------------------------------------------------------------------
c
c                   Check for a Valid File Name
c                   ***************************
c
c This subroutine takes the character string "str" of length "len" and
c removes all leading blanks and blanks out all characters after the
c first blank found in the string (leading blanks are removed first).
c
c
c
c-----------------------------------------------------------------------
      parameter (MAXLEN=512)
      character str(MAXLEN)*1
c
c Find first two blanks and blank out remaining characters:
c
      do i=1,len-1
            if(str(i)  .eq.' '.and.
     +         str(i+1).eq.' ') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 2
            end if
      end do
 2    continue
c
c Look for "-fi" for file
c
      do i=1,len-2
            if(str(i)  .eq.'-'.and.
     +         str(i+1).eq.'f'.and.
     +         str(i+2).eq.'i') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 3
            end if
      end do
 3    continue
c
c Look for "\fi" for file
c
      do i=1,len-2
            if(str(i)  .eq.'\'.and.
     +         str(i+1).eq.'f'.and.
     +         str(i+2).eq.'i') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 4
            end if
      end do
 4    continue
c
c Return with modified file name:
c
      return
      end
