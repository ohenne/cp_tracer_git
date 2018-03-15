! identifies objects and establishes links
! needs input data as SRV file
! Compile: ifort -no-wrap-margin -o irt_objects_release.x irt_objects_release.f90
! to dos:
! timestep handling> at timestep x new tracers are set and at the same timestep
! tracers are updated alread. DOes the new position corresponds to timestep x or
! x+1?
PROGRAM irt_cp_tracking

USE netcdf

USE irt_parameters, ONLY: domainsize_x, domainsize_y, lperiodic, &
    time_steps, nt_bins, nx_bins, ny_bins, &
    max_no_of_cells, miss, edge_fraction, dt, &
    tracer_steps, resolution, dt, max_tfields, max_tracers
! resolution in m
! dt in sec

IMPLICIT NONE

INTEGER              :: domsize_x, domsize_y, domsize_z
INTEGER              :: ii, ij,ix,iy,idx,idy 
REAL, ALLOCATABLE    :: input_field(:,:)       ! eg precip field to find COGs
REAL, ALLOCATABLE    :: vel(:,:,:,:)           ! velocity field
REAL, ALLOCATABLE    :: QC(:,:,:), QG(:,:,:)   ! passive tracer 
REAL, ALLOCATABLE    :: track_numbers(:,:)     ! ID for precip and coldpools objects
REAL, ALLOCATABLE    :: nneighb(:,:)           ! identifier for cell boundaries
REAL, ALLOCATABLE    :: traced(:,:)            ! traced information
REAL, ALLOCATABLE    :: tracerfield(:,:,:)     ! field of tracer to check results
INTEGER              :: i,j,t,it,fileid,zi
INTEGER              :: counter                ! cell counter
REAL                 :: COMx(max_no_of_cells), COMy(max_no_of_cells)
REAL                 :: zm(75), zt(74)
INTEGER              :: already_tracked(max_no_of_cells) ! memory of cell counter
INTEGER              :: srv_header_input(8)

! time handling
INTEGER              :: date, time, timestep, time_step_event, Skip

! auxiliary variables
REAL                 :: vx,vy
INTEGER              :: nevent
INTEGER              :: ntracer, count_tracer
INTEGER              :: ntr
INTEGER              :: onset

! file names
CHARACTER (len=90)   :: input_fn_tracks
CHARACTER (len=90)   :: mask_nneighb
CHARACTER (len=90)   :: mask_u
CHARACTER (len=90)   :: mask_v
CHARACTER (len=90)   :: mask_rv
CHARACTER (len=90)   :: mask_tracer

INTEGER              :: ierr

    real :: start, finish
    call cpu_time(start)
! open the info file to get information on timestep starting
999 FORMAT(12X, I3)
!open(unit=100,file='info.txt',status='old',action=&
!     'read', iostat=ierr)
!if ( ierr == 0) then
!    read(100,999) onset 
!else
!     write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', &
!                 ierr,' aufgetreten'
!end if
onset = 136
WRITE(input_fn_tracks,"(A19)"),    "input/irt_tracks_mask.srv"


! open the input filenames
OPEN(1,FILE='input/irt_objects_input_00.srv',FORM='unformatted', ACTION='read')
OPEN(2,FILE='input/irt_tracks_mask.srv',    FORM='unformatted', ACTION='read')
OPEN(3,FILE='input/irt_tracks_sorted.txt',FORM='formatted', ACTION='read')
OPEN(4,FILE='input/level.txt',FORM='formatted', ACTION='read')

! open the output filenames
OPEN(36,FILE='irt_objects_tracer.srv',FORM='unformatted', ACTION='write')
OPEN(40,FILE='cp_3Dhistory.txt',FORM='formatted', ACTION='write')

! text output file, write header
151 FORMAT  (2X,A4,   1X,A6,   3X,A3, 2X,A7,   1X,A4,  2(8X, A4),    7X,A6,   3(2X,A4),2X, &
  A1, 5(X,A10))
WRITE(40,151) 'time','tstart','age','traceID','cpID', 'Xpos','Ypos','Height','Xpos','Ypos','Zpos', &
  'F', ' LWC/gkg-1',' GWC/gkg-1', '    w/ms-1', 'rad' , 'dist'
! get the levels
DO 
  READ(4,*,END=300) zi, zm(zi)
END DO
300 CONTINUE
zt = (zm(1:74)+zm(2:75))/2.

ntracer      = 1 ! counts tracer patches
count_tracer = 1 ! counts individual pixels !OCH was ist mit pixeln gemeint? der
!tracer selber?

counter = 0
timestep=-1 !OCH why -1 !-1
already_tracked(:) = 0

! OCH: DO WE NEED THIS??
! If periodic boundary conditions are switched off, a 1-gridbix thick 
! frame of missing values will be laid around the field, and the domainsize
! has to be increased by 2 in both dimensions
domsize_z = 74
IF (lperiodic) THEN
  domsize_x = domainsize_x
  domsize_y = domainsize_y
ELSE
  domsize_x = domainsize_x+2
  domsize_y = domainsize_y+2
ENDIF

! allocating fields for the different quantities
ALLOCATE(input_field(domsize_x,domsize_y))
ALLOCATE(nneighb(domsize_x,domsize_y)) ! 1=nn, 2=u wind of boundary, 3=v wind of boundary, 4=event number.
ALLOCATE(vel(domsize_z,domsize_x,domsize_y,3))
ALLOCATE(track_numbers(domsize_x,domsize_y))
! allocate array to hold all traced particles
ALLOCATE(traced(INT(max_tracers),11))
ALLOCATE(tracerfield(domsize_x,domsize_y,max_tfields))
ALLOCATE(QC(domsize_z,domsize_x,domsize_y))
ALLOCATE(QG(domsize_z,domsize_x,domsize_y))

count_tracer=1
READ(3,*) Skip,time_step_event

! beginning of main loop
WRITE(*,*) "beginning main loop "
DO
 date=srv_header_input(3)
 time=srv_header_input(4)
 timestep=timestep+1 ! OCH: starts with 0 ?
 WRITE(*,*) "at",timestep
 
 ! reading the standard input fields
 READ(1,END=200) srv_header_input
 WRITE(*,*) "RAED 1 input/irt_objects_input_00.srv header"
 IF (lperiodic) THEN
   READ(1) input_field(:,:)
 ELSE
   input_field(:,:) = miss-1.
   READ(1) input_field(2:domsize_x-1,2:domsize_y-1)
 ENDIF
 WRITE(*,*) "RAED 1 input/irt_objects_input_00.srv file"
 
 track_numbers(:,:)=0
 IF (timestep .GE. max(onset,time_step_event)) THEN
    ! reading velocity field and passive tracer
    CALL read_in_3d(vel(:,:,:,1), 'u', 'input/test1plus4K.out.vol.u.nc',timestep,2,domsize_z)
    CALL read_in_3d(vel(:,:,:,2), 'v', 'input/test1plus4K.out.vol.v.nc',timestep,2,domsize_z)
    CALL read_in_3d(vel(:,:,:,3), 'w', 'input/test1plus4K.out.vol.w.nc',timestep,1,domsize_z)
    CALL read_in_3d(QC, 'l', '/nbi/ac/conv1/henneb/modeloutput/test1plus4K/level1/test1plus4K.out.vol.l.nc'&
                   ,timestep,2,domsize_z)
    CALL read_in_3d(QG, 'rgrp','/nbi/ac/conv1/henneb/modeloutput//test1plus4K/level1/test1plus4K.out.vol.rgrp.nc'&
                   ,timestep,2,domsize_z)
    ! reading the track input files
    READ(2,END=200) srv_header_input   

    IF (lperiodic) THEN
       READ(2) track_numbers(:,:)
    ELSE
       track_numbers(:,:) = miss-1.
       READ(2) track_numbers(2:domsize_x-1,2:domsize_y-1)
    ENDIF
 ENDIF
 ! initializing/clearing several fields
 nneighb(:,:)=0
 tracerfield(:,:,:)=0
 counter=1
 
 ! identification of edges
 CALL neigh(domsize_x,domsize_y,track_numbers,nneighb,COMx,COMy,input_field,max_no_of_cells,already_tracked)
 CALL set_tracer(counter,domsize_x,domsize_y,nneighb, track_numbers, vel(1:2,:,:,:), & 
      COMx, COMy, timestep,traced,edge_fraction,max_no_of_cells, &
      ntracer,count_tracer,max_tracers,already_tracked,zt(1))
 CALL update_tracer(vel(:,:,:,1),vel(:,:,:,2),vel(:,:,:,3),domsize_x,domsize_y, domsize_z,&
      timestep,traced,resolution,dt,count_tracer,max_tracers,track_numbers,zt,zm,QC,QG)
 ! find closed lines around COG by connecting tracer
! CALL find_outlines(timestep,traced,COMx,COMy,already_tracked,max_no_of_cells,max_tracers,count_tracer) 
! OCH: might be better to separate some routines: 
 CALL outlines(traced,max_tracers,COMx,COMy,already_tracked,max_no_of_cells)
 CALL write_output(traced,max_tracers,count_tracer,QC,QG,zt,domsize_z,domsize_x,domsize_y,timestep) 
 CALL fill_tracer(timestep,traced,COMx,COMy,already_tracked,max_no_of_cells,max_tracers,count_tracer)
 ! writing tracers to a grid (for time step output)
 tracerfield(:,:,:) = 0
 ntr = 1
 DO WHILE (ntr .LT. count_tracer)
    ix = INT(traced(ntr,1))
    iy = INT(traced(ntr,2))
 
    IF ((ix .GT. 0) .AND. (iy .GT. 0)) THEN
       ! "speed" field
       ! tracerfield(ix,iy,1) = sqrt(traced(ntr,3)**2+traced(ntr,4)**2)
       tracerfield(ix,iy,1) = traced(ntr,9)
       ! "index" field
       tracerfield(ix,iy,2) = traced(ntr, 8)
    ENDIF
    ntr=ntr+1
 ENDDO
 
 IF (lperiodic) THEN
!   CALL write_srv(domainsize_x,domainsize_y,REAL(nneighb(:,:)),date,time,31)
   CALL write_srv(domainsize_x,domainsize_y,REAL(tracerfield(:,:,1)),date,time,36)
 !  CALL write_srv(domainsize_x,domainsize_y,REAL(nneighb(:,:,3)),date,time,36)
 ENDIF
 
! end main loop (time)
!WRITE(*,*) timestep
ENDDO
WRITE(*,*) 'finished main loop'
200 CONTINUE

! close input/output files
CLOSE(11)
CLOSE(36)
CLOSE(1)
call cpu_time(finish)
print '("Time = ",f10.3," seconds.")',finish-start

CONTAINS
! -----------------------------------------------------------------------
!subroutines for reading and writing netcdf
! -----------------------------------------------------------------------

  SUBROUTINE read_in_3d (poutput,varname, filename,ctime,zstart,zend)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  real, dimension(:,:,:), intent(inout) :: poutput
  integer :: ncId, rhVarId, nx, ny, nz, nt, ctime
  integer, intent(in) :: zstart,zend
  integer, dimension(nf90_max_var_dims) :: dimIDs
!  real, allocatable, dimension(:,:,:) ::  zvar

    CALL check(nf90_open(filename, nf90_NoWrite, ncid))
    CALL check(nf90_inq_varid(ncid,varname, rhVarId))
    CALL check(nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(4), len = nt))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(3), len = nx))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(2), len = ny))
    CALL check(nf90_inquire_dimension(ncid, dimIDs(1), len = nz))
!    allocate(zvar(nx,ny,nz)) !,nx, ny, nz, nt))
  !  status = nf90_get_var(ncid, rhVarId, zvar,start = (/ 1, 2, 1,1 /))
    CALL check(nf90_get_var(ncid, rhVarId, poutput, start =(/zstart,1,1,ctime/), &
                                                    count= (/ zend,nx, ny, 1 /)))
!    poutput = zvar(:,:,:)
!WRITE(*,*) varname," nz",nz,"ny",ny,"nx",nx,"nt",nt, zstart, zend, ctime  
    CALL check(nf90_close(ncid))

!    deallocate(zvar)
!    write(*,*) "read data: ",filename
  END SUBROUTINE read_in_3d

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check


END PROGRAM irt_cp_tracking

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                           SUBROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE set_tracer(counter,domsize_x,domsize_y,nneighb,track_numbers,vel,center_of_mass_x, &
     center_of_mass_y,timestep,traced,edge_fraction,& 
     max_no_of_cells,ntracer,count_tracer,max_tracers,already_tracked,zt)

  INTEGER, INTENT(IN)	    :: domsize_x,domsize_y, counter
  REAL, INTENT(INOUT)       :: nneighb(domsize_x,domsize_y)
  REAL, INTENT(IN)          :: track_numbers(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: timestep
  INTEGER, INTENT(IN)       :: max_no_of_cells
  INTEGER, INTENT(IN)       :: max_tracers
  REAL, INTENT(IN)          :: zt
  INTEGER, INTENT(INOUT)    :: ntracer
  REAL, INTENT(IN)          :: edge_fraction
  REAL                      :: center_of_mass_x(max_no_of_cells),center_of_mass_y(max_no_of_cells)
  REAL, INTENT(INOUT)       :: traced(max_tracers,11)
  REAL, INTENT(IN)          :: vel(2,domsize_x,domsize_y,3)
  INTEGER, INTENT(INOUT)    :: count_tracer
  INTEGER                   :: ix, iy
  REAL                      :: vx, vy
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  DO iy=1, domsize_y
     DO ix=1, domsize_x
        IF (nneighb(ix,iy) .EQ. 1 ) THEN 
           vx             = .5*(vel(1,ix,iy,1)+vel(1,MOD(ix+1-1,domsize_x)+1,iy,1)) ! average vx
           vy             = .5*(vel(1,ix,iy,2)+vel(1,ix,MOD(iy+1-1,domsize_y)+1,2)) ! average vy

           IF (already_tracked(INT(track_numbers(ix,iy))) .LT. 300) THEN
              IF (track_numbers(ix,iy) .GT. 0 .AND. track_numbers(ix,iy) .GE. 1) THEN 
                 traced(count_tracer,1) = ix 
                 traced(count_tracer,2) = iy
                 traced(count_tracer,3) = zt    ! in meter not gridpoints, start
                                                ! at full level 
!now used for rad and dist
!nicht genutzt    traced(count_tracer,4) = vx 
!                 traced(count_tracer,5) = vy 
                 traced(count_tracer,6) = timestep     ! timestep when tracking begins
                 traced(count_tracer,7) = 0   ! tracer age ???timestep     ! 
!nichtgenutzt                 traced(count_tracer,8) = counter              ! ! OCH ist das die
!track ID Tracer ID
!nein? was fuer ein counter ist das?
                 traced(count_tracer,9) = track_numbers(ix,iy) ! TRACK ID 
                 traced(count_tracer,10) = 1.                   ! has tracer been involved in other precip events
                 traced(count_tracer,11) = 1.
                 count_tracer           = count_tracer + 1     ! Das ist die Tracer ID ??? 
!>>OCH 
!150 format (2X,I10,2X,I10,2X,I4,2X, I4,2X,I3,2X,I3,2X,I1)  ! firstout put
!timetsep macht mehr sinn
!                 WRITE(40,150) ix,iy,timestep,0,counter,INT(track_numbers(ix,iy)), 1
!<<OCH
!ix, iy, timestep, counter, track_numbers(ix,iy)

              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO
! updating track numbers which have already been used for tracking -> set already_tracked=TRUE
  DO iy=1, domsize_y
     DO ix=1, domsize_x
        IF (track_numbers(ix,iy) .GE. 1 .AND. nneighb(ix,iy) .EQ. 1 ) THEN 
           already_tracked(INT(track_numbers(ix,iy))) = already_tracked(INT(track_numbers(ix,iy))) +1
        ENDIF
     ENDDO
  ENDDO
  
  RETURN
  
END SUBROUTINE set_tracer

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE update_tracer(velx,vely,velz,domsize_x,domsize_y,domsize_z, &
     timestep,traced,resolution, &
     dt, count_tracer, max_tracers,track_numbers,zt,zm, QC,QG)
  INTEGER, INTENT(IN)       :: domsize_x, domsize_y, domsize_z, timestep
  REAL, INTENT(IN)          :: velx(domsize_z,domsize_x,domsize_y), &
                               vely(domsize_z,domsize_x,domsize_y), &
                               velz(domsize_z,domsize_x,domsize_y),&
                               QC(domsize_z,domsize_x,domsize_y), &
                               QG(domsize_z,domsize_x,domsize_y)
  REAL, INTENT(IN)          :: resolution
  REAL, INTENT(IN)          :: track_numbers(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: count_tracer,max_tracers
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_tracers,11)
  REAL, INTENT(IN)          :: dt
  REAL, INTENT(IN)          :: zt(domsize_z), zm(domsize_z+1)
  REAL                      :: ix_new, iy_new, hh_new, vx_intp, vy_intp, vz_intp
  REAL                      :: ix, iy, hh
  INTEGER                   :: iz_dummy(1)
  INTEGER                   :: ix_round, iy_round,iz_round, tottracer, start_time, ix_org, iy_org
  INTEGER                   :: ix_m, iy_m 
  INTEGER                   :: ix_round_new, iy_round_new, iz_new
  REAL                      :: vx_lft, vx_rgt, vy_bot, vy_top,  wgt_x, wgt_y, wgt_zh, wgt_zf
  REAL                      :: LWC, GWC
  it=1 
  ! updating previous tracers
  DO WHILE (it .LT. count_tracer)        
    IF (traced(it,11)  .eq. 1.) THEN  !trace only if tracer is active
      ! determining the first time step of the event
      start_time=traced(it,6)
      ! determining how many timesteps have passed since then
      tracer_ts =timestep-start_time+1  
     
      IF (start_time .GT. 0) THEN
        ! getting the previous positions
        ix=traced(it,1)
        iy=traced(it,2) ! position in gridpoints
        hh=traced(it,3) ! height in meter
        IF (hh .lt. 10000.) THEN
          ix_round=MOD(INT(ix)-1+domsize_x,domsize_x)+1
          iy_round=MOD(INT(iy)-1+domsize_y,domsize_y)+1
          ix_m=MOD(INT(ix)-2+domsize_x,domsize_x)+1
          iy_m=MOD(INT(iy)-2+domsize_y,domsize_y)+1

          ! get weights for vertical interpolation 
          wgt_x  = MOD(ix,1.)
          wgt_y  = MOD(iy,1.)

          ! in z for half level
          izm = 1
          DO WHILE (zm(izm)-hh .lt. 0. )  ! Loop untill levelheight(iz) is below hh
            izm= izm+1                   ! First level above hh
          END DO
          ! and full level
          izt = 1
          DO WHILE (zt(izt)-hh .lt. 0. )  
            izt= izt+1
          END DO
          ! distance to the next upper level =weight to the upper level(izm)          
          wgt_zh = ((zm(izm)-hh))/(zm(izm)-zm(izm-1))   
          wgt_zf = ((zt(izt)-hh))/(zt(izt)-zt(izt-1))   

          ! 3D inter
          CALL trilininterpol(vx_intp,velx(izt-1:izt,(/ix_m,ix_round/),(/iy_m,iy_round/)),wgt_x,wgt_y,wgt_z)
          CALL trilininterpol(vy_intp,vely(izt-1:izt,(/ix_m,ix_round/),(/iy_m,iy_round/)),wgt_x,wgt_y,wgt_z)
          CALL trilininterpol(vz_intp,velz(izm-1:izm,(/ix_m,ix_round/),(/iy_m,iy_round/)),wgt_x,wgt_y,wgt_z)

! kann vlt weg    IF (ix .GT. 0 .AND. iy .GT. 0) THEN ! bogus now?

          ! get new location as decimal
          ix_new = MOD((ix + dt*vx_intp/resolution)-1.+FLOAT(domsize_x),FLOAT(domsize_x))+1.
          iy_new = MOD((iy + dt*vy_intp/resolution)-1.+FLOAT(domsize_y),FLOAT(domsize_y))+1.
          hh_new = MAX(hh + dt*vz_intp,50.) !keep tracer minimum on first full level zt(1)  

          ! and as gridded values
          ix_round_new= INT(ix_new) !MOD(INT(ix_new)-1+domsize_x,domsize_x)+1
          iy_round_new=INT(iy_new)  !MOD(INT(iy_new)-1+domsize_y,domsize_y)+1
          iz_dummy=MINLOC(ABS(zt(:)-hh_new)) 
          iz_new = iz_dummy(1)

          ! save new values for next loop
          traced(it,1) = ix_new
          traced(it,2) = iy_new
          traced(it,3) = hh_new
!          traced(it,4) = vx_new do we need to save this values? 
!          traced(it,5) = vy_new
          traced(it,7) = tracer_ts
          traced(it,11) = 1
          !IF (traced(it,10) .eq. 0) THEN    ! if tracer is dead cant get back 
          !  traced(it,9) = 0  ! stays dead
          !  ! else (if it was alive) it keeps living as long no other precip
          !  ! event is above
          !ELSEIF (track_numbers(ix_round,iy_round) .eq. traced(it,9) .or. track_numbers(ix_round,iy_round) .le. 0. ) THEN
          !  traced(it,10) = 1
          !ELSE
          !  traced(it,10) = 0
          !ENDIF 
!          CALL outlines(traced,max_tracers,it,COMx,COMy)
! call outline routine 
! make separate out put routine: 
!          LWC = QC(iz_new,ix_round_new,iy_round_new) *1000.
!          GWC = QG(iz_new,ix_round_new,iy_round_new) *1000.
!          150 FORMAT   (2X,I4,   3X,I4,            2X,I4    ,3X,I5, 2X,I4, &
!                        2(2X,F10.5),  2X,F11.5,3(2X,I4),                        2X,I1,&
!                        3(2X,F10.6))
!          WRITE(40,150) timestep,INT(traced(it,6)),tracer_ts,it    ,INT(traced(it,9)),&
!                        ix_new,iy_new,hh_new ,ix_round_new,iy_round_new,iz_new, INT(traced(it,10)), &
!                         LWC, GWC, vz_intp 
        ELSE ! set tracer inactive if it is to high 
          traced(it,11) = 0.
        ENDIF ! end if: check if tracer is too high
      ENDIF  ! end if start time .gt. 0
    END IF ! check if tracer is active   
    it = it + 1
  ENDDO
  RETURN
END SUBROUTINE update_tracer

SUBROUTINE write_output(traced,max_tracers,count_tracer,QC,QG,zt,domsize_z,domsize_x,domsize_y,timestep)
  INTEGER, INTENT(IN)       :: count_tracer,max_tracers, timestep, &
                               domsize_z, &
                               domsize_x,domsize_y
  INTEGER                   :: it, ix_round_new, iy_round_new, iz_dummy(1), iz_new 
  REAL, INTENT(INOUT)       :: traced(max_tracers,11)
  REAL, INTENT(IN)          :: QC(domsize_z,domsize_x,domsize_y), &
                               QG(domsize_z,domsize_x,domsize_y)
  REAL                      :: LWC, GWC 
  REAL, INTENT(IN)          :: zt(domsize_z)
  it=1
  ! updating previous tracers
  DO WHILE (it .LT. count_tracer)
    IF (traced(it,11)  .eq. 1.) THEN  !trace only if tracer is active
        IF (traced(it,3) .lt. 10000.) THEN
          ! and as gridded values
          ix_round_new= INT(traced(it,1)) !MOD(INT(ix_new)-1+domsize_x,domsize_x)+1
          iy_round_new=INT(traced(it,2))  !MOD(INT(iy_new)-1+domsize_y,domsize_y)+1
          iz_dummy=MINLOC(ABS(zt(:)-traced(it,3)))
          iz_new = iz_dummy(1)

          LWC = QC(iz_new,ix_round_new,iy_round_new) *1000.
          GWC = QG(iz_new,ix_round_new,iy_round_new) *1000.
          150 FORMAT   (2X,I4,   3X,I4,              3X,I5, 2X,I4, &
                        2(2X,F10.5),  2X,F11.5,3(2X,I4),2X,I1,&
                        5(2X,F10.6))
          WRITE(40,150) timestep,INT(traced(it,6)),it,INT(traced(it,9)),&
                        traced(it,1),traced(it,2),traced(it,3),ix_round_new,iy_round_new,iz_new, INT(traced(it,10)), &
                         LWC, GWC, vz_intp, traced(it,4),traced(it,5)


        END IF
    END IF
    it = it+1
  END DO

END SUBROUTINE write_output


!SUBROUTINE outlines(traced,max_tracers,it)
!  INTEGER, INTENT(IN)       :: max_tracers,it
!  REAL, INTENT(INOUT)       :: traced(max_tracers,11)
!  REAL                      :: DELTAx, DELTAy, pi
!
!   pi = 2.*asin(1.)
!
!
!
!  DELTAx = traced(it,1) - COMx(INT(traced(it,9)))
!  DELTAy = traced(it,2) - COMy(INT(traced(it,9)))
!
!  ! calculate the angle of every tracer to COG
!  IF (DELTAx .gt. 0 .and. DELTAy .gt. 0) THEN ! 1st quadrant
!     traced(it,4) = tan(DELTAy/DELTAx)
!  ELSE IF (DELTAx .gt. 0 .and. DELTAy .lt. 0) THEN ! 2nd
!     traced(it,4) = tan(DELTAy/DELTAx) + pi
!  ELSE IF (DELTAx .lt. 0 .and. DELTAy .lt. 0) THEN ! 3rd
!     traced(it,4) = tan(DELTAy/DELTAx) + pi
!  ELSE IF (DELTAx .gt. 0 .and. DELTAy .gt. 0) THEN ! 4th
!     traced(it,4) = tan(DELTAy/DELTAx) +2.* pi
!  ELSE IF (DELTAx .gt. 0 .and. DELTAy .eq. 0) THEN
!     traced(it,4) = pi/2.
!  ELSE IF (DELTAx .lt. 0 .and. DELTAy .eq. 0) THEN
!     traced(it,4) = (3./2.) * pi
!  ELSE IF (DELTAx .eq. 0 .and. DELTAy .gt. 0) THEN
!     traced(it,4) = 0.
!  ELSE IF (DELTAx .eq. 0 .and. DELTAy .lt. 0) THEN
!     traced(it,4) =  pi
!  ENDIF
!  ! and the distance do the COG
!  traced(it,5) = sqrt(DELTAx**2+DELTAy**2)
!
!END SUBROUTINE outlines

SUBROUTINE outlines(traced,max_tracers,COMx,COMy,already_tracked,max_no_of_cells)
  INTEGER, INTENT(IN)       :: max_tracers,max_no_of_cells
  REAL, INTENT(INOUT)       :: traced(max_tracers,11)
  REAL                      :: DELTAx(max_tracers), DELTAy(max_tracers), pi
  INTEGER                   :: i
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  REAL, INTENT(INOUT)    :: COMx(max_no_of_cells),COMy(max_no_of_cells)

   pi = 2.*asin(1.)

  i = 1
  DO WHILE (already_tracked(i) .gt. 0) 
    WHERE (i .eq. INT(traced(:,9)))
      DELTAx = traced(:,1) - COMx(i)
      DELTAy = traced(:,2) - COMy(i)
    END WHERE
     ! calculate the angle of every tracer to COG
    WHERE (DELTAx .gt. 0 .and. DELTAy .gt. 0)  ! 1st quadrant
       traced(:,4) = tan(DELTAy/DELTAx)
    ELSE WHERE (DELTAx .gt. 0 .and. DELTAy .lt. 0)  ! 2nd
       traced(:,4) = tan(DELTAy/DELTAx) + pi
    ELSE WHERE (DELTAx .lt. 0 .and. DELTAy .lt. 0)  ! 3rd
       traced(:,4) = tan(DELTAy/DELTAx) + pi
    ELSE WHERE (DELTAx .gt. 0 .and. DELTAy .gt. 0)  ! 4th
       traced(:,4) = tan(DELTAy/DELTAx) +2.* pi
    ELSE WHERE (DELTAx .gt. 0 .and. DELTAy .eq. 0) 
       traced(:,4) = pi/2.
    ELSE WHERE (DELTAx .lt. 0 .and. DELTAy .eq. 0) 
       traced(:,4) = (3./2.) * pi
    ELSE WHERE (DELTAx .eq. 0 .and. DELTAy .gt. 0) 
       traced(:,4) = 0.
    ELSE WHERE (DELTAx .eq. 0 .and. DELTAy .lt. 0) 
       traced(:,4) =  pi
    ENDWHERE
  END DO 
  ! and the distance do the COG
  traced(:,5) = sqrt(DELTAx**2+DELTAy**2)

END SUBROUTINE outlines

! ++++++++++++++++++++++++++++++

SUBROUTINE fill_tracer(timestep,traced,COMx,COMy,already_tracked,max_no_of_cells,max_tracers, count_tracer)
!SUBROUTINE add_tracer(timestep,traced,COMx,COMy,already_tracked,max_no_of_cells,max_tracers,count_tracer)


  INTEGER, INTENT(IN)       :: already_tracked(max_no_of_cells)
  INTEGER, INTENT(IN)       :: max_tracers
  INTEGER, INTENT(IN)       :: timestep
  INTEGER, INTENT(IN)       :: max_no_of_cells
  INTEGER, INTENT(INOUT)    :: count_tracer
  REAL, INTENT(INOUT)       :: COMx(max_no_of_cells),COMy(max_no_of_cells)
  REAL                      :: DELTAx(max_tracers), DELTAy(max_tracers)
  REAL, INTENT(INOUT)       :: traced(max_tracers,11) 
  REAL                      :: traced_sort(max_tracers,11)
  INTEGER                   :: i, j, k, l, counter, start
  REAL                      :: rad(max_tracers), dist(max_tracers), &
                               xpos(max_tracers), ypos(max_tracers), cpID(max_tracers) 
  REAL                      :: pi, a, dst

   rad = traced(1:count_tracer,4)
   dist = traced(1:count_tracer,5)
 ! sort by angle for every cpID
  ! initialize
  i = 1 ! rain cell resp cp ID
  start =2
  DO WHILE (already_tracked(i) .gt. 0) ! loop trough all cps with tracer
    do j=start,start+already_tracked(i) !   count_tracer!  already_tracked(i)
      a=rad(j)  !save value at j (eg 2)
      do k=j-1,1,-1 ! go backwarts from value below current
        ! if  val before current value is smaller than actual value
        if (rad(k)<=a) goto 10
        ! set this value at the current position 
        rad(k+1)=rad(k)
        dist(k+1)=dist(k)
        cpID(k+1)=traced(k,9)
        xpos(k+1)=traced(k,1)
        ypos(k+1)=traced(k,2)
        traced_sort(k+1,1) = traced(k,1)
        traced_sort(k+1,2) = traced(k,2)
        traced_sort(k+1,4) = dist(k)
        traced_sort(k+1,5) = rad(k)
        traced_sort(k+1,9) = traced(k,9)

        ! DELTA X abspeichern
      end do
      k=0
      10 rad(k+1)=a
         dist(k+1)=b
    end do
    start = already_tracked(i)+2
    i = i +1
  end do 
  
 ! OCH: first smooth the distance to cog before checking the distance 
!   loop durcj cps
!   WHERE(cp = cp)
!     ratio = rad/mean(rad)  !wird dann mean nur ueber rad where genommen?
!   ENDWHERE
!   xpos = xpos + deltax  mal ratio
!   y pos = y pos +deltay mal ratio
  counter = 1
  do l =2,count_tracer
    if (count_tracer+counter .lt. max_tracers) THEN
     if (cpID(l)  .eq. cpID(l-1)) THEN! only if tracer belong to same CP
       dst =  sqrt((xpos(l)- xpos(l-1))**2 + (ypos(l)- ypos(l-1))**2)
       if (dst .gt. 4.)  then! neighbouring tracers have a distance larger than 4gp 
         ! set new tracer
         traced(count_tracer+counter,1) = MOD(xpos(l-1) + (xpos(l)- xpos(l-1))+320.-1.,320.)+1.
         traced(count_tracer+counter,2) = MOD(ypos(l-1) + (ypos(l)- ypos(l-1))+320.-1.,320.)+1.
         traced(count_tracer+counter,3) = 50.
         traced(count_tracer+counter,6) = timestep
         traced(count_tracer+counter,7) = 0
!         traced(count_tracer+counter,8) = count_tracer+counter           
         traced(count_tracer+counter,9) = cpID(l)
         traced(count_tracer+counter,10) = 1. 
         traced(count_tracer+counter,11) = 1. 
         counter = counter+1
       end if
     end if
    end if
  end do
  count_tracer = count_tracer+counter
END SUBROUTINE fill_tracer !find_outlines

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE write_srv(domainsize_x,domainsize_y,field,date,time,file_id)

  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: domainsize_x,domainsize_y
  REAL, INTENT(IN)      :: field(domainsize_x,domainsize_y)
  INTEGER, INTENT(IN)   :: date,time
  INTEGER, INTENT(IN)   :: file_id
  INTEGER               :: srv_header(8)

  srv_header(1) = 1	      ! Code
  srv_header(2) = 1	      ! Level
  srv_header(3) = date        ! Datum
  srv_header(4) = time        ! Zeitinkrement
  srv_header(5) = domainsize_x
  srv_header(6) = domainsize_y
  srv_header(7) = 0
  srv_header(8) = 0

  WRITE (file_id) srv_header
  WRITE (file_id) field
  
  RETURN

END SUBROUTINE write_srv


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! identify boundaries of the precipitation cells
SUBROUTINE neigh(domsize_x,domsize_y,track_numbers,nneighb,COMx,COMy,input_field,max_no_of_cells,already_tracked)

  IMPLICIT NONE
  INTEGER                :: i,j, imodm, imodp, jmodm,jmodp
  INTEGER, INTENT(IN)    :: domsize_x,domsize_y
  REAL, INTENT(IN)       :: track_numbers(domsize_x,domsize_y)
  REAL, INTENT(INOUT)    :: nneighb(domsize_x,domsize_y) 
  REAL, INTENT(INOUT)    :: input_field(domsize_x,domsize_y)
  INTEGER, INTENT(IN)    :: max_no_of_cells
  REAL, INTENT(INOUT)    :: COMx(max_no_of_cells),COMy(max_no_of_cells) 
  REAL                   :: COUNTER(max_no_of_cells)
  INTEGER, INTENT(IN)    :: already_tracked(max_no_of_cells)

  COMx = 0
  COMy = 0
  COUNTER = 0
  nneighb = 0
  DO i =2,domsize_x-1
    DO j =2,domsize_y-1
     IF (already_tracked(INT(track_numbers(i,j))) .lt. 30 ) THEN 
     ! already_tracked counts how many tracer exist for a single CP, 
     ! CPID is gien in track_numbers(i,j)
     ! if already 30 tracers are set, no more new tracers will be started
     ! associated with the precipitation object, oly by interpolating
      imodp = mod(i+domsize_x,domsize_x)+1
      imodm = mod(i-2+domsize_x,domsize_x)+1
      jmodp = mod(j+domsize_y,domsize_y)+1
      jmodm = mod(j-2+domsize_y,domsize_y)+1

      IF (track_numbers(i,j) .ne. track_numbers(imodp,j)) nneighb((/i,imodp/),j) =1  
      IF (track_numbers(i,j) .ne. track_numbers(imodm,j)) nneighb((/i,imodm/),j) =1 
      IF (track_numbers(i,j) .ne. track_numbers(i,jmodp)) nneighb(i,(/j,jmodp/)) =1
      IF (track_numbers(i,j) .ne. track_numbers(i,jmodm)) nneighb(i,(/j,jmodm/)) =1 

      ! diagonal neigh
      IF (track_numbers(i,j) .ne. track_numbers(imodp,jmodp)) nneighb(imodp,jmodp) =1
      IF (track_numbers(i,j) .ne. track_numbers(imodm,jmodp)) nneighb(imodm,jmodp) =1
      IF (track_numbers(i,j) .ne. track_numbers(imodp,jmodm)) nneighb(imodp,jmodm) =1
      IF (track_numbers(i,j) .ne. track_numbers(imodm,jmodm)) nneighb(imodm,jmodm) =1  

      IF (track_numbers(i,j) .gt. 0 .and. track_numbers(i,j) .lt. max_no_of_cells ) THEN
        COMx(INT(track_numbers(i,j))) = COMx(INT(track_numbers(i,j)))  + i*input_field(i,j) 
        COMy(INT(track_numbers(i,j))) = COMy(INT(track_numbers(i,j)))  + j*input_field(i,j)
        COUNTER(INT(track_numbers(i,j))) = COUNTER(INT(track_numbers(i,j)))+ input_field(i,j)
      END IF
     END IF
    END DO
  END DO
  ! dont forget boundaries
  COMx=COMx/COUNTER
  COMy=COMy/COUNTER

END SUBROUTINE neigh


SUBROUTINE trilininterpol(v_intp,vel,wgt_x,wgt_y,wgt_z)

  IMPLICIT NONE
  REAL, INTENT(INOUT) :: v_intp
  REAL, INTENT(IN) :: vel(2,2,2)
  REAL, INTENT(IN) :: wgt_x,wgt_y,wgt_z


  v_intp = vel(1,1,1) *(1.-wgt_x)*(1.-wgt_y)*wgt_z &
         + vel(1,1,2) *(1.-wgt_x)*(   wgt_y)*wgt_z &
         + vel(1,2,2) *(   wgt_x)*(   wgt_y)*wgt_z &
         + vel(1,2,1) *(   wgt_x)*(1.-wgt_y)*wgt_z &
         + vel(2,1,1) *(1.-wgt_x)*(1.-wgt_y)*(1.-wgt_z) &
         + vel(2,1,2) *(1.-wgt_x)*(   wgt_y)*(1.-wgt_z) &
         + vel(2,2,2) *(   wgt_x)*(   wgt_y)*(1.-wgt_z) &
         + vel(2,2,1) *(   wgt_x)*(1.-wgt_y)*(1.-wgt_z) 

END SUBROUTINE trilininterpol
