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
    miss, edge_fraction, dt, &
    tracer_steps, resolution, dt, max_tfields
!) resolution in m
! dt in sec

IMPLICIT NONE

INTEGER              :: domsize_x, domsize_y, domsize_z
INTEGER              :: ii, ij,ix,iy,idx,idy 
REAL, ALLOCATABLE    :: input_field(:,:)       ! eg precip field to find COGs
REAL, ALLOCATABLE    :: vel(:,:,:,:)           ! velocity field
REAL, ALLOCATABLE    :: QC(:,:,:), QG(:,:,:)   ! passive tracer 
REAL, ALLOCATABLE    :: track_numbers(:,:)     ! ID for precip and coldpools objects
REAL, ALLOCATABLE    :: nneighb(:,:)           ! identifier for cell boundaries
REAL, ALLOCATABLE    :: traced(:,:,:)            ! traced information
REAL, ALLOCATABLE    :: tracerfield(:,:,:)     ! field of tracer to check results
INTEGER              :: i,j,t,it,fileid,zi
INTEGER              :: counter                ! cell counter
INTEGER              :: max_no_of_cells
REAL, ALLOCATABLE    :: COMx(:), COMy(:)
REAL                 :: zm(75), zt(74)
INTEGER, ALLOCATABLE :: already_tracked(:) ! memory of cell counter
INTEGER              :: srv_header_input(8)
INTEGER              :: max_tracers
! time handling
INTEGER              :: date, time, timestep, time_step_event, Skip

! auxiliary variables
REAL                 :: vx,vy
INTEGER              :: nevent
INTEGER              :: ntracer, count_tracer
INTEGER              :: ntr
INTEGER              :: onset
INTEGER              :: max_tracer_CP
INTEGER,ALLOCATABLE  :: tracpo(:,:)
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
WRITE(input_fn_tracks,"(A19)"),    "input/irt_tracks_mask.srv"


! open the input filenames
OPEN(1,FILE='input/irt_objects_input_00.srv',FORM='unformatted', ACTION='read')
OPEN(2,FILE='input/irt_tracks_mask.srv',    FORM='unformatted', ACTION='read')
OPEN(3,FILE='input/irt_tracks_sorted.txt',FORM='formatted', ACTION='read')
OPEN(4,FILE='input/level.txt',FORM='formatted', ACTION='read')

! open the output filenames
OPEN(36,FILE='irt_objects_tracer.srv',FORM='unformatted', ACTION='write')
OPEN(40,FILE='cp_3Dhistory.txt',FORM='formatted', ACTION='write')
OPEN(50,FILE='cp_3Dsorted.txt',FORM='formatted', ACTION='write')

! text output file, write header
151 FORMAT  (2X,A4,   1X,A6,  2X,A7,   1X,A4,  2(8X, A4),    7X,A6,   3(2X,A4),2X, &
  A1, 5(X,A10))
152 FORMAT (  2X,A4 , X,A5   ,2(8X,A5),         11X,A2, 2X,A3, 8X,A5, 9X,A4)
WRITE(40,151) 'time','tstart','traceID','cpID', 'Xpos','Ypos','Height','Xpos','Ypos','Zpos', &
  'F', ' LWC/gkg-1',' GWC/gkg-1', '    w/ms-1', 'rad' , 'dist'
WRITE(50,152) 'time', 'CP ID', 'x pos', 'y pos', 'hh' ,'phi', 'D2cog','ravg' 
! get the levels
DO 
  READ(4,*,END=300) zi, zm(zi)
END DO
300 CONTINUE
zt = (zm(1:74)+zm(2:75))/2.

max_tracer_CP = 600 ! automate
! OCH Automate onset again, automate max no of cells according to rain tracking
! output
onset = 136
max_no_of_cells = 290
max_tracers = max_no_of_cells*max_tracer_CP
ntracer      = 1 ! counts tracer patches
count_tracer = 1 ! counts individual pixels !OCH was ist mit pixeln gemeint? der
!tracer selber?

counter = 0
timestep=-1 !OCH why -1 !-1

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
ALLOCATE(traced(max_no_of_cells,max_tracer_CP,11))
ALLOCATE(tracerfield(domsize_x,domsize_y,max_tfields))
ALLOCATE(QC(domsize_z,domsize_x,domsize_y))
ALLOCATE(QG(domsize_z,domsize_x,domsize_y))
ALLOCATE(tracpo(2,max_tracers))
ALLOCATE(COMx(max_no_of_cells))
ALLOCATE(COMy(max_no_of_cells))
ALLOCATE(already_tracked(max_no_of_cells))
already_tracked(:) = 0
READ(3,*) Skip,time_step_event

! beginning of main loop
WRITE(*,*) "beginning main loop "
DO
IF (timestep .gt. 240) GOTO 5000
 date=srv_header_input(3)
 time=srv_header_input(4)
 timestep=timestep+1 ! OCH: starts with 0 ?
 WRITE(*,*) "at",timestep
 
 ! reading the standard input fields
 READ(1,END=200) srv_header_input
 IF (lperiodic) THEN
   READ(1) input_field(:,:)
 ELSE
   input_field(:,:) = miss-1.
   READ(1) input_field(2:domsize_x-1,2:domsize_y-1)
 ENDIF
 
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
      timestep,traced,edge_fraction,max_no_of_cells, &
      ntracer,count_tracer,max_tracer_CP,already_tracked,zt(1),tracpo, max_tracers)
 CALL update_tracer(vel(:,:,:,1),vel(:,:,:,2),vel(:,:,:,3),domsize_x,domsize_y, domsize_z,&
      timestep,traced,resolution,dt,count_tracer,max_tracer_CP,track_numbers,zt,zm,max_no_of_cells,&
      tracpo,max_tracers)
 ! find closed lines around COG by connecting tracer
! CALL find_outlines(timestep,traced,COMx,COMy,already_tracked,max_no_of_cells,max_tracers,count_tracer) 
! OCH: might be better to separate some routines: 
!IF (timestep .gt. 136) then
!write(*,*) "before outl ", already_tracked
!end if
 CALL outlines(traced,max_tracer_CP,COMx,COMy,already_tracked,max_no_of_cells,REAL(domsize_x),REAL(domsize_y))
 CALL write_output(traced,max_tracers,count_tracer,QC,QG,zt,domsize_z,domsize_x,domsize_y,timestep, &
                   tracpo,max_no_of_cells, max_tracer_CP) 
write(*,*) "vorher",sum(already_tracked)
 CALL fill_tracer(timestep,traced,already_tracked,max_no_of_cells,max_tracer_CP,count_tracer)
write(*,*) "nacher",sum(already_tracked)
 ! writing tracers to a grid (for time step output)
 tracerfield(:,:,:) = 0
 ntr = 1
! DO WHILE (ntr .LT. count_tracer)
!    ix = INT(traced(ntr,1))
!    iy = INT(traced(ntr,2))
! 
!    IF ((ix .GT. 0) .AND. (iy .GT. 0)) THEN
!       ! "speed" field
!       ! tracerfield(ix,iy,1) = sqrt(traced(ntr,3)**2+traced(ntr,4)**2)
!       tracerfield(ix,iy,1) = traced(ntr,9)
!       ! "index" field
!       tracerfield(ix,iy,2) = traced(ntr, 8)
!    ENDIF
!    ntr=ntr+1
! ENDDO
 
 IF (lperiodic) THEN
!   CALL write_srv(domainsize_x,domainsize_y,REAL(nneighb(:,:)),date,time,31)
   CALL write_srv(domainsize_x,domainsize_y,REAL(tracerfield(:,:,1)),date,time,36)
 !  CALL write_srv(domainsize_x,domainsize_y,REAL(nneighb(:,:,3)),date,time,36)
 ENDIF
 
! end main loop (time)
!WRITE(*,*) timestep
ENDDO
5000 CONTINUE
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

SUBROUTINE set_tracer(counter,domsize_x,domsize_y,nneighb,track_numbers,vel, &
     timestep,traced,edge_fraction,& 
     max_no_of_cells,ntracer,count_tracer,max_tracer_CP,already_tracked,zt,tracpo,max_tracers)

  INTEGER, INTENT(IN)	    :: domsize_x,domsize_y, counter
  REAL, INTENT(INOUT)       :: nneighb(domsize_x,domsize_y)
  REAL, INTENT(IN)          :: track_numbers(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: timestep
  INTEGER, INTENT(IN)       :: max_no_of_cells , max_tracers
  INTEGER, INTENT(IN)       :: max_tracer_CP
  REAL, INTENT(IN)          :: zt
  INTEGER, INTENT(INOUT)    :: tracpo(2,max_tracers)
  INTEGER, INTENT(INOUT)    :: ntracer
  REAL, INTENT(IN)          :: edge_fraction
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,11)
  REAL, INTENT(IN)          :: vel(2,domsize_x,domsize_y,3)
  INTEGER, INTENT(INOUT)    :: count_tracer
  INTEGER                   :: ix, iy
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  
  DO iy=1, domsize_y
     DO ix=1, domsize_x
        IF (nneighb(ix,iy) .EQ. 1 ) THEN 
! min size of rain cell 256 (16^2) -> quadratic cell has 4*16 boundary gp
!bisher werden nicht genug tracer durch auffuellen gesetzt ....
            IF (already_tracked(INT(track_numbers(ix,iy))) .LT. 150) THEN !max_tracer_CP ) THEN !64) THEN
              IF (track_numbers(ix,iy) .GT. 0 .AND. track_numbers(ix,iy) .GE. 1) THEN 
                 already_tracked(INT(track_numbers(ix,iy)))=already_tracked(INT(track_numbers(ix,iy)))+1

                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),1) = ix 
                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),2) = iy
                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),3) = zt    ! in meter not gridpoints, start
!Write(*,*) "set height ", zt
                                                ! at full level 
                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),6) = timestep     ! timestep when tracking begins
                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),7) = 0   ! tracer age ???timestep     ! 
!nichtgenutzt                 traced(count_tracer,8) = counter              ! ! OCH ist das die
!track ID Tracer ID
!nein? was fuer ein counter ist das?
                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),9) = track_numbers(ix,iy) ! TRACK ID 
                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),10) = 1.                   ! has tracer been involved in other precip events
                 traced(INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy))),11) = 1.
!                 tracpo(:,count_tracer)=(/INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))/)
                 tracpo(:,count_tracer)=(/INT(track_numbers(ix,iy)),already_tracked(INT(track_numbers(ix,iy)))/)
                 count_tracer           = count_tracer + 1     ! Das ist die Tracer ID ??? 
!write(*,*) 'endif'
              ENDIF
           ENDIF
     c = c+1
        ENDIF
     ENDDO
  ENDDO
! updating track numbers which have already been used for tracking -> set already_tracked=TRUE
!  DO iy=1, domsize_y
!     DO ix=1, domsize_x
!        IF (track_numbers(ix,iy) .GE. 1 .AND. nneighb(ix,iy) .EQ. 1 ) THEN 
!           already_tracked(INT(track_numbers(ix,iy))) = already_tracked(INT(track_numbers(ix,iy))) +1
!        ENDIF
!     ENDDO
!  ENDDO
  
  RETURN
  
END SUBROUTINE set_tracer

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE update_tracer(velx,vely,velz,domsize_x,domsize_y,domsize_z, &
     timestep,traced,resolution, &
     dt, count_tracer, max_tracer_CP,track_numbers,zt,zm,max_no_of_cells,tracpo,max_tracers)
  INTEGER, INTENT(IN)       :: domsize_x, domsize_y, domsize_z, &
                               timestep,max_no_of_cells, max_tracers
  REAL, INTENT(IN)          :: velx(domsize_z,domsize_x,domsize_y), &
                               vely(domsize_z,domsize_x,domsize_y), &
                               velz(domsize_z,domsize_x,domsize_y)
  REAL, INTENT(IN)          :: resolution
  REAL, INTENT(IN)          :: track_numbers(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: count_tracer,max_tracer_CP
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,11)
  REAL, INTENT(IN)          :: dt
  REAL, INTENT(IN)          :: zt(domsize_z), zm(domsize_z+1)
  REAL                      :: ix_new, iy_new, hh_new, vx_intp, vy_intp, vz_intp
  REAL                      :: ix, iy, hh
  INTEGER                   :: iz_dummy(1)
  INTEGER                   :: ix_round, iy_round,iz_round, tottracer, start_time, ix_org, iy_org
  INTEGER                   :: ix_m, iy_m 
  INTEGER                   :: ix_round_new, iy_round_new, iz_new
  REAL                      :: vx_lft, vx_rgt, vy_bot, vy_top,  wgt_x, wgt_y, wgt_zh, wgt_zf
  INTEGER, INTENT(IN)       :: tracpo(2,max_tracers)
  it=1 
  ! updating previous tracers
  DO WHILE (it .LT. count_tracer)       

    IF (traced(tracpo(1,it),tracpo(2,it),11)  .eq. 1.) THEN  !trace only if tracer is active
      ! determining the first time step of the event
      start_time=traced(tracpo(1,it),tracpo(2,it),6)
      ! determining how many timesteps have passed since then
      tracer_ts =timestep-start_time+1  
     
      IF (start_time .GT. 0) THEN
        ! getting the previous positions
        ix=traced(tracpo(1,it),tracpo(2,it),1)
        iy=traced(tracpo(1,it),tracpo(2,it),2) ! position in gridpoints
        hh=traced(tracpo(1,it),tracpo(2,it),3) ! height in meter
        IF (hh .lt. 2500.) THEN
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
!write(*,*) izm-1, izm, velz(izm-1:izm,(/ix_m/),(/iy_round/))
!write(*,*) izm-1, izm, velz(izm-1:izm,(/ix_m/),(/iy_m/))
!write(*,*) izm-1, izm, velz(izm-1:izm,(/ix_round/),(/iy_m/))
!write(*,*) izm-1, izm, velz(izm-1:izm,(/ix_round/),(/iy_round/))

!         make heavy tracer 
          traced(tracpo(1,it),tracpo(2,it),8) = vz_intp
!IF (vz_intp .gt. 5.) THEN
! write(*,*) timestep, vz_intp ,tracpo(1,it)
!END IF
          vz_intp = max(0.,vz_intp - 5.)   ! add fall velocity of 3m/s to reduce vertical advection
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
          traced(tracpo(1,it),tracpo(2,it),1) = ix_new
          traced(tracpo(1,it),tracpo(2,it),2) = iy_new
          traced(tracpo(1,it),tracpo(2,it),3) = hh_new
!          traced(tracpo(1,it),tracpo(2,it),4) = vx_new do we need to save this values? 
!          traced(tracpo(1,it),tracpo(2,it),5) = vy_new
          traced(tracpo(1,it),tracpo(2,it),6) = timestep
          traced(tracpo(1,it),tracpo(2,it),7) = tracer_ts
          traced(tracpo(1,it),tracpo(2,it),11) = 1
          IF (hh .gt. 1000) THEN
            traced(tracpo(1,it),tracpo(2,it),10) = 0
          ELSE
            traced(tracpo(1,it),tracpo(2,it),10) = 1
          END IF

          !IF (traced(tracpo(1,it),tracpo(2,it),10) .eq. 0) THEN    ! if tracer is dead cant get back 
          !  traced(tracpo(1,it),tracpo(2,it),9) = 0  ! stays dead
          !  ! else (if it was alive) it keeps living as long no other precip
          !  ! event is above
          !ELSEIF (track_numbers(ix_round,iy_round) .eq. traced(tracpo(1,it),tracpo(2,it),9) .or. track_numbers(ix_round,iy_round) .le. 0. ) THEN
          !  traced(tracpo(1,it),tracpo(2,it),10) = 1
          !ELSE
          !  traced(tracpo(1,it),tracpo(2,it),10) = 0
          !ENDIF 
!          CALL outlines(traced,max_tracers,it,COMx,COMy)
! call outline routine 
! make separate out put routine: 
!          LWC = QC(iz_new,ix_round_new,iy_round_new) *1000.
!          GWC = QG(iz_new,ix_round_new,iy_round_new) *1000.
!          150 FORMAT   (2X,I4,   3X,I4,            2X,I4    ,3X,I5, 2X,I4, &
!                        2(2X,F10.5),  2X,F11.5,3(2X,I4),                        2X,I1,&
!                        3(2X,F10.6))
!          WRITE(40,150) timestep,INT(traced(tracpo(1,it),tracpo(2,it),6)),tracer_ts,it    ,INT(traced(tracpo(1,it),tracpo(2,it),9)),&
!                        ix_new,iy_new,hh_new ,ix_round_new,iy_round_new,iz_new, INT(traced(tracpo(1,it),tracpo(2,it),10)), &
!                         LWC, GWC, vz_intp 
        ELSE ! set tracer inactive if it is to high 
          traced(tracpo(1,it),tracpo(2,it),11) = 0.
        ENDIF ! end if: check if tracer is too high
      ENDIF  ! end if start time .gt. 0
    END IF ! check if tracer is active   
    it = it + 1
  ENDDO
  RETURN
END SUBROUTINE update_tracer


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE write_output(traced,max_tracers,count_tracer,QC,QG,zt,domsize_z,domsize_x,domsize_y,timestep,tracpo,&
                       max_no_of_cells, max_tracer_CP)
  INTEGER, INTENT(IN)       :: count_tracer,max_tracers, timestep, &
                               domsize_z, &
                               domsize_x,domsize_y
  INTEGER                   :: it, ix_round_new, iy_round_new, iz_dummy(1), iz_new
  INTEGER,INTENT(IN)        :: tracpo(2,max_tracers) 
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells, max_tracer_CP,11)
  REAL, INTENT(IN)          :: QC(domsize_z,domsize_x,domsize_y), &
                               QG(domsize_z,domsize_x,domsize_y)
  REAL                      :: LWC, GWC 
  REAL, INTENT(IN)          :: zt(domsize_z)
  it=1
  ! updating previous tracers
  DO WHILE (it .LT. count_tracer)
    IF (traced(tracpo(1,it),tracpo(2,it),11)  .eq. 1.) THEN  !trace only if tracer is active
        IF (traced(tracpo(1,it),tracpo(2,it),3) .lt. 10000.) THEN
          ! and as gridded values
          ix_round_new= INT(traced(tracpo(1,it),tracpo(2,it),1)) !MOD(INT(ix_new)-1+domsize_x,domsize_x)+1
          iy_round_new=INT(traced(tracpo(1,it),tracpo(2,it),2))  !MOD(INT(iy_new)-1+domsize_y,domsize_y)+1
          iz_dummy=MINLOC(ABS(zt(:)-traced(tracpo(1,it),tracpo(2,it),3)))
          iz_new = iz_dummy(1)

          LWC = QC(iz_new,ix_round_new,iy_round_new) *1000.
          GWC = QG(iz_new,ix_round_new,iy_round_new) *1000.
          150 FORMAT   (2X,I4,  3X,I4,   &
                        2X,I6, 2X,I4, &
                        2(2X,F10.5),  2X,F11.5,3(2X,I4),2X,I1,&
                        4(2X,F10.6), 2X,F10.3)
          WRITE(40,150) timestep,INT(traced(tracpo(1,it),tracpo(2,it),6)),& !timestep, timestart
                        it,INT(traced(tracpo(1,it),tracpo(2,it),9)),&       !tarcer and CP ID
                        traced(tracpo(1,it),tracpo(2,it),1),traced(tracpo(1,it),tracpo(2,it),2),& !pos 1-2
                        traced(tracpo(1,it),tracpo(2,it),3),ix_round_new,iy_round_new,iz_new, & !pos 3 and rounded
                        INT(traced(tracpo(1,it),tracpo(2,it),10)), &
                        LWC, GWC, traced(tracpo(1,it),tracpo(2,it),8), traced(tracpo(1,it),tracpo(2,it),4),&
                        traced(tracpo(1,it),tracpo(2,it),5)


        END IF
    END IF
    it = it+1
  END DO

END SUBROUTINE write_output

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE outlines(traced,max_tracer_CP,COMx,COMy,already_tracked,max_no_of_cells,domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: max_tracer_CP,max_no_of_cells
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,11)
  REAL                      :: DELTAx(max_no_of_cells,max_tracer_CP), &
                               DELTAy(max_no_of_cells,max_tracer_CP), pi
  INTEGER                   :: i,j
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  REAL, INTENT(INOUT)    :: COMx(max_no_of_cells),COMy(max_no_of_cells)
  REAL, INTENT(IN)          :: domsize_x,domsize_y
   pi = 2.*asin(1.)

  DO i = 1,max_no_of_cells ! loop trough every cp
!write(*,*) "traced ", already_tracked(i)
   IF (already_tracked(i) .gt. 0) THEN ! do only sth if there are already tracerfor the CP

    DO j=1,already_tracked(i)

      DELTAx(i,j) = mod(traced(i,j,1) -COMx(i)-1.+domsize_x,domsize_x)+1.
      DELTAy(i,j) = mod(traced(i,j,2) - COMy(i)-1.+domsize_y,domsize_y)+1.
!write(*,*) traced(i,j,2), COMy(i), domsize_y, mod(traced(i,j,2) - COMy(i)-1.+domsize_y,domsize_y)+1.
!write(*,*) "delta y",DELTAy(i,j)
!write(*,*) "at CP", i,"traced", already_tracked(i), traced(i,j,1),traced(i,j,2), "COGS", COMx(i), COMy(i),"deltas",DELTAx(i,j),DELTAy(i,j)
!write(*,*) domsize_y, domsize_x
      IF      (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .gt. 0) THEN
       traced(i,j,4) = atan(DELTAx(i,j)/DELTAy(i,j))
      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 2nd
       traced(i,j,4) = atan(DELTAx(i,j)/DELTAy(i,j)) + pi
      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .lt. 0) THEN  ! 3nd
       traced(i,j,4) = atan(DELTAx(i,j)/DELTAy(i,j)) + pi
      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .gt. 0) THEN  ! 4nd
       traced(i,j,4) = atan(DELTAx(i,j)/DELTAy(i,j)) +2.* pi
      ELSE IF (DELTAx(i,j) .gt. 0 .and. DELTAy(i,j) .eq. 0)THEN
       traced(i,j,4) = pi/2.
      ELSE IF (DELTAx(i,j) .lt. 0 .and. DELTAy(i,j) .eq. 0) THEN
       traced(i,j,4) = (3./2.) * pi
      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .gt. 0) THEN
       traced(i,j,4) = 0.
      ELSE IF (DELTAx(i,j) .eq. 0 .and. DELTAy(i,j) .lt. 0) THEN
       traced(i,j,4) =  pi
      END IF
      traced(i,j,5) = sqrt(DELTAx(i,j)**2.+DELTAy(i,j)**2.)
     END DO
    END IF
   END DO



!    END DO
!   END IF
!  END DO
!  DO i = 1,max_no_of_cells ! loop trough every cp
!   IF (already_tracked(i) .gt. 0) THEN ! do only sth if there are already tracerfor the CP
!     ! calculate the angle of every tracer to COG
!    WHERE (DELTAx .gt. 0 .and. DELTAy .gt. 0)  ! 1st quadrant
!       traced(:,:,4) = atan(DELTAy/DELTAx)
!    ELSE WHERE (DELTAx .gt. 0 .and. DELTAy .lt. 0)  ! 2nd
!       traced(:,:,4) = atan(DELTAy/DELTAx) + pi
!    ELSE WHERE (DELTAx .lt. 0 .and. DELTAy .lt. 0)  ! 3rd
!       traced(:,:,4) = atan(DELTAy/DELTAx) + pi
!    ELSE WHERE (DELTAx .gt. 0 .and. DELTAy .gt. 0)  ! 4th
!       traced(:,:,4) = atan(DELTAy/DELTAx) +2.* pi
!    ELSE WHERE (DELTAx .gt. 0 .and. DELTAy .eq. 0) 
!       traced(:,:,4) = pi/2.
!    ELSE WHERE (DELTAx .lt. 0 .and. DELTAy .eq. 0) 
!       traced(:,:,4) = (3./2.) * pi
!    ELSE WHERE (DELTAx .eq. 0 .and. DELTAy .gt. 0) 
!       traced(:,:,4) = 0.
!    ELSE WHERE (DELTAx .eq. 0 .and. DELTAy .lt. 0) 
!       traced(:,:,4) =  pi
!    ENDWHERE
!   END IF
!  END DO 
  ! and the distance do the COG
!  traced(:,:,5) = sqrt(DELTAx**2.+DELTAy**2.)

END SUBROUTINE outlines

! ++++++++++++++++++++++++++++++

SUBROUTINE fill_tracer(timestep,traced,already_tracked,max_no_of_cells,max_tracer_CP, count_tracer)
!SUBROUTINE add_tracer(timestep,traced,COMx,COMy,already_tracked,max_no_of_cells,max_tracers,count_tracer)


  INTEGER, INTENT(INOUT)       :: already_tracked(max_no_of_cells)
  INTEGER, INTENT(IN)       :: max_tracer_CP
  INTEGER, INTENT(IN)       :: timestep
  INTEGER, INTENT(IN)       :: max_no_of_cells
  INTEGER, INTENT(INOUT)    :: count_tracer
  INTEGER                   :: count_tracer_new
  REAL, INTENT(INOUT)       :: traced(max_no_of_cells,max_tracer_CP,11) 
  REAL                      :: traced_dummy(max_tracer_CP,11)
  INTEGER                   :: i,jc, j, k, l, counter, start

   count_tracer_new = count_tracer
   DO i =1,max_no_of_cells ! loop trough all cps with tracer
!    jc = 1
    IF (already_tracked(i) .gt. 0 ) THEN
!    IF (already_tracked(i) .gt. 0 .and. already_tracked(i) .lt. max_tracer_CP) THEN ! loop trough all cps with tracer
!!     DO j = 1,count_tracer ! loop trough current number of tracer
!!       IF (traced(j,9) .eq. i .and. traced(j,11) .eq. 1) THEN !active tracer belonging to CPi  
!!       traced_dummy(jc,:) = traced(jc,(/1,2,4,5/)) 
!!       jc = jc+1 
!!      END IF
!!     END DO
     traced_dummy(:,:)=traced(i,:,:)
     CALL sort(traced_dummy(1:already_tracked(i),:),already_tracked(i))
write(*,*) "vorher",already_tracked(i)
     CALL add_tracer(traced_dummy(1:already_tracked(i),:),traced(i,:,:),max_no_of_cells,&
                     already_tracked(i),i,timestep,count_tracer_new,&
                     counter,max_tracer_CP)
write(*,*) "nachher",already_tracked(i)
     count_tracer = count_tracer+counter
    END IF
    
   END DO
!   count_tracer = count_tracer_new
END SUBROUTINE fill_tracer 


SUBROUTINE sort(traced_dummy,jc)

  INTEGER, INTENT(IN)      :: jc
  REAL, INTENT(INOUT)      :: traced_dummy(jc,11)
  REAL                     :: a(11)

 ! sort by angle for every cpID:
  ! initialize
    do j=1,jc !   count_tracer!  
      a=traced_dummy(j,:)  !save value at j 
      do k=j-1,1,-1 ! go backwarts from value below current
        ! if  val before current value is smaller than actual value
        if (traced_dummy(k,4)<= a(4)) goto 10
        ! set this value at the current position 
        traced_dummy(k+1,:)=traced_dummy(k,:)

        ! DELTA X abspeichern
      end do
      k=0
      10 traced_dummy(k+1,:)=a
    end do
END SUBROUTINE sort


SUBROUTINE add_tracer(traced_dummy,traced,max_no_of_cells,already_tracked,CPI,timestep,count_tracer_new,counter,max_tracer_CP)

  INTEGER, INTENT(IN) :: CPI, timestep,count_tracer_new,max_no_of_cells,&
                         max_tracer_CP
  INTEGER, INTENT(INOUT) :: already_tracked
  INTEGER, INTENT(OUT) :: counter
  REAL, INTENT(IN)    :: traced_dummy(already_tracked,11)
  REAL, INTENT(INOUT) :: traced(max_tracer_CP,11)
  REAL                :: xpos(max_tracer_CP), ypos(max_tracer_CP), &
                         phi(max_tracer_CP), r(max_tracer_CP), rthr(max_tracer_CP), &
                         maxdist, dst, pi2, abstand, weighted, sumweight, &
                         xposnew, yposnew, phip, cruc_dst
  INTEGER             :: l,lm,lp, ind, k ,c, il(already_tracked)
   pi2 = 4.*asin(1.)
! take only tracer below 150m 
  c = 0
  DO l =1,already_tracked 
   IF (traced_dummy(l,3) .lt. 150.) then
     c = c+1
     phi(c) = traced_dummy(l,4)
     r(c) = traced_dummy(l,5)
     il(c) = l 
   END IF 
  END DO
! CALCULATE OUTER radius thr
  do l =1,c !already_tracked !2,count_tracer
    lm = mod(l-5+already_tracked-1,already_tracked)+1
    lp = mod(l+5+already_tracked-1,already_tracked)+1

    maxdist = mod(phi(l)- phi(lm)+pi2,pi2)+mod(phi(lp)-phi(l)+pi2,pi2)
    weighted = 0.
    sumweight = 0.

    do k = -5,5,1
      ind = mod(l+k+already_tracked-1,already_tracked)+1
      if (k .lt. 0 ) then
        abstand = mod(phi(l) - phi(ind) +pi2,pi2)
      else if (k .gt. 0) then
        abstand = mod(phi(ind) - phi(l) +pi2,pi2)
      else if (k .eq. 0) then
        abstand =0.
      end if
      weighted = weighted + r(ind)* (maxdist-abstand  )*r(ind) 
      sumweight = sumweight + ((maxdist-abstand)*r(ind))
    end do
    rthr(l) = weighted/sumweight
    ! sorted tracer as output
    if (traced_dummy(il(l),6) .eq. timestep) then
    WRITE(50,200) INT(traced_dummy(il(l),6)), INT(traced_dummy(il(l),9)), &             !timestep, CP ID
                  traced_dummy(il(l),1), traced_dummy(il(l),2), traced_dummy(il(l),3), & !position
                  traced_dummy(il(l),4), traced_dummy(il(l),5), &             !angle and dist to cog
                  rthr(l) ! averaged radius  
    end if
  end do
  
  counter = 0
  xpos = traced_dummy(:,1)
  ypos = traced_dummy(:,2)


! add tracer in between
write(*,*) "c, already tracked ",c ,already_tracked 
  do l =1,c !already_tracked !2,count_tracer
    if (already_tracked .lt. max_tracer_CP) THEN ! maximal 300 tracer allowed per CP
      ! take only outer tracer 
      if (r(il(l)) .gt. rthr(l)) then
        ii = 1
        do while (r(il(l)+ii) .lt. rthr(l+ii) ) 
         ii = ii+1
        end do
        ! compare to the next outer tracer
        if (r(il(l)+ii) .gt. rthr(l+ii) ) then
           dst = mod(phi(il(l)+ii)- phi(il(l)) + pi2,pi2) ! radial differences
           !!tracer distance should be less than one gp: 
           cruc_dst = (pi2/8.)/max(rthr(l+ii),1.) !tracer distance should be less
           xposnew = xpos(il(l))
           yposnew = ypos(il(l))
           phip    = phi(il(l))
           do while (dst .lt. cruc_dst)
            counter = counter +1
            already_tracked = already_tracked + 1 !counter
            xposnew= mod(xposnew-sin(phip)*rthr(l)+sin(phip+cruc_dst)*rthr(l)+320.-1.,320.)
            yposnew= mod(yposnew-cos(phip)*rthr(l)+cos(phip+cruc_dst)*rthr(l)+320.-1.,320.)
            ! set new tracer
            traced(already_tracked,1) = xposnew
            traced(already_tracked,2) = yposnew 
            traced(already_tracked,3) = 50.
            traced(already_tracked,6) = timestep
            traced(already_tracked,7) = 0
!            traced(count_tracer+counter,8) = count_tracer+counter           
            traced(already_tracked,9) = CPI
            traced(already_tracked,10) = 1.
            traced(already_tracked,11) = 1.
!write(*,*) "problem", l, phip, dst, cruc_dst, rthr(l+ii)
            phip = mod(phip +cruc_dst + pi2,pi2)
            dst = mod(phi(il(l)+ii)- phip + pi2,pi2)            
           end do
             
!          dst =  sqrt((xpos(il(l))- xpos(il(l)+ii))**2 + (ypos(il(l))- ypos(il(l)+ii))**2)
!          ! if next outer tracer is to far, set new in btw
!          if (dst .gt. 0.5)  then! neighbouring tracers have a distance larger than 4gp 
!write(*,*) "set new tracer for", CPI, "between", xpos(il(l)+ii) , "and", xpos(il(l)) 
!            counter = counter +1
!            ! set new tracer
!            already_tracked = already_tracked + 1 !counter
!            traced(already_tracked,1) = MOD(xpos(il(l)+ii) + (xpos(il(l))- xpos(il(l)+ii))+320.-1.,320.)+1.
!            traced(already_tracked,2) = MOD(ypos(il(l)+ii) + (ypos(il(l))- ypos(il(l)+ii))+320.-1.,320.)+1.
!            traced(already_tracked,3) = 50.
!            traced(already_tracked,6) = timestep
!            traced(already_tracked,7) = 0
!            traced(count_tracer+counter,8) = count_tracer+counter           
!            traced(already_tracked,9) = CPI
!            traced(already_tracked,10) = 1. 
!            traced(already_tracked,11) = 1.
!          end if
        end if
      end if
    end if
  end do
 200 FORMAT   (2(2X,I4) ,   3(2X,F11.5), 2X,F5.3, 2(2X,F11.3))

END SUBROUTINE add_tracer 

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
  REAL                   :: COGx(max_no_of_cells),COGy(max_no_of_cells)
  COGx = COMx ! save COMs already calculated
  COGy = COMy
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
  ! if preicpitation is disipated, keep the old COG
  WHERE (COUNTER .eq. 0) COMx = COGx
  WHERE (COUNTER .eq. 0) COMy = COGy

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
