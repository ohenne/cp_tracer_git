! identifies objects and establishes links
! needs input data as SRV file
! Compile: ifort -no-wrap-margin -o irt_objects_release.x irt_objects_release.f90
! to dos:
! timestep handling> at timestep x new tracers are set and at the same timestep
! tracers are updated alread. DOes the new position corresponds to timestep x or
! x+1?
PROGRAM irt_cp_tracking

USE irt_parameters, ONLY: domainsize_x, domainsize_y, lperiodic, &
    time_steps, nt_bins, nx_bins, ny_bins, &
    max_no_of_cells, miss, edge_fraction, dt, &
    tracer_steps, resolution, dt, max_tfields, max_tracers
! resolution in m
! dt in sec

IMPLICIT NONE

INTEGER              :: domsize_x, domsize_y
INTEGER              :: ii, ij,ix,iy,idx,idy 
REAL, ALLOCATABLE    :: input_field(:,:)     ! eg precip field
REAL, ALLOCATABLE    :: vel(:,:,:)             ! velocity field
REAL, ALLOCATABLE    :: track_numbers(:,:)     ! ID for precip and coldpools objects
REAL, ALLOCATABLE    :: nneighb(:,:)           ! identifier for cell boundaries
REAL, ALLOCATABLE    :: traced(:,:)            ! traced information
REAL, ALLOCATABLE    :: tracerfield(:,:,:)     ! field of tracer to check results
INTEGER              :: i,j,t,it,fileid
INTEGER              :: counter                ! cell counter
REAL                 :: COMx(max_no_of_cells), COMy(max_no_of_cells)
INTEGER              :: already_tracked(max_no_of_cells) ! memory of cell counter

! WINDFIELD INPUT
INTEGER              :: srv_header_input(8)
INTEGER              :: srv_header_coarse(8)
INTEGER              :: headlen

! time handling
INTEGER              :: date, time, timestep, time_step_event, Skip

! auxiliary variables
REAL                 :: vx,vy
REAL                 :: cm_x,cm_y
REAL                 :: vel_x,vel_y
REAL                 :: rel_x,rel_y
INTEGER              :: nevent
INTEGER              :: ntracer, count_tracer
INTEGER              :: ntr
INTEGER              :: onset

! file names
CHARACTER (len=90)   :: input_filename
CHARACTER (len=90)   :: input_filename_v(2)
CHARACTER (len=90)   :: input_fn_tracks
CHARACTER (len=90)   :: output_filename
CHARACTER (len=90)   :: mask_filename
CHARACTER (len=90)   :: mask_nneighb
CHARACTER (len=90)   :: mask_u
CHARACTER (len=90)   :: mask_v
CHARACTER (len=90)   :: mask_rv
CHARACTER (len=90)   :: mask_tracer
CHARACTER (len=90)   :: coarsevel_filename

INTEGER              :: ierr

    real :: start, finish
    call cpu_time(start)
!
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
WRITE(*,*) "OCH onset", onset
WRITE(input_filename,"(A24)") "input/irt_objects_input_00.srv"
WRITE(input_filename_v(1),"(A29)"),"input/irt_objects_input_u.srv"
WRITE(input_filename_v(2),"(A29)"),"input/irt_objects_input_v.srv"
WRITE(input_fn_tracks,"(A25)"),    "input/irt_tracks_mask.srv"

! defining the output filenames
!output_filename    = "irt_objects_output.txt"
!mask_filename      = "irt_objects_mask.srv"
mask_tracer        = "irt_objects_tracer.srv"
!coarsevel_filename = "irt_advection_field.srv"

! open the input filenames
OPEN(1,FILE='input/irt_objects_input_00.srv',FORM='unformatted', ACTION='read')
OPEN(10,FILE=trim(input_filename_v(1)),FORM='unformatted', ACTION='read')
OPEN(11,FILE=trim(input_filename_v(2)),FORM='unformatted', ACTION='read')
OPEN(12,FILE=trim(input_fn_tracks),    FORM='unformatted', ACTION='read')
OPEN(13,FILE='input/irt_tracks_sorted.txt',FORM='formatted', ACTION='read')

! open the output filenames
OPEN(36,FILE=trim(mask_tracer),FORM='unformatted', ACTION='write')
OPEN(40,FILE='cp_history.txt',FORM='formatted', ACTION='write')

! text output file, write header
151 format  (2X,A4,     1X,A6,3X,A3,    2X,A7, 1X,A4,     2(8X, A4),     2(2X,A4),2X, A1)
WRITE(40,151) 'time','tstart','age','traceID','cpID', 'Xpos','Ypos','Xpos','Ypos', 'F' 

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
ALLOCATE(vel(domsize_x,domsize_y,2))
ALLOCATE(track_numbers(domsize_x,domsize_y))
! allocate array to hold all traced particles
ALLOCATE(traced(INT(max_tracers),9))
ALLOCATE(tracerfield(domsize_x,domsize_y,max_tfields))

count_tracer=1
READ (13,*) Skip,time_step_event
! beginning of main loop
WRITE(*,*) "beginning main loop of irt_objects_release"
DO

 date=srv_header_input(3)
 time=srv_header_input(4)
 timestep=timestep+1 ! OCH: starts with 0 ?
 
 ! reading the velocity input files
 DO fileid=10,11
    READ (fileid,END=200) srv_header_input   
    IF (lperiodic) THEN
       READ (fileid) vel(:,:,fileid-9)
    ELSE
       vel(:,:,fileid-9) = miss-1.
       READ (fileid) vel(2:domsize_x-1,2:domsize_y-1,fileid-9)
    ENDIF
 ENDDO
 
 ! reading the standard input fields
 READ (1,END=200) srv_header_input
 IF (lperiodic) THEN
   READ (1) input_field(:,:)
 ELSE
   input_field(:,:) = miss-1.
   READ (1) input_field(2:domsize_x-1,2:domsize_y-1)
 ENDIF
 
 track_numbers(:,:)=0
 IF (timestep .GE. max(onset,time_step_event)) THEN
    ! reading the track input files
    READ (12,END=200) srv_header_input   
    IF (lperiodic) THEN
       READ (12) track_numbers(:,:)
    ELSE
       track_numbers(:,:) = miss-1.
       READ (12) track_numbers(2:domsize_x-1,2:domsize_y-1)
    ENDIF
 ENDIF
 
 ! initializing/clearing several fields
 nneighb(:,:)=0
 tracerfield(:,:,:)=0
 counter=1
 
 ! identification of edges
 CALL neigh(domsize_x,domsize_y,track_numbers,nneighb,COMx,COMy,input_field,max_no_of_cells)
       CALL set_tracer(counter,domsize_x,domsize_y,nneighb, track_numbers, vel, & 
            COMx, COMy, timestep,traced,edge_fraction,max_no_of_cells, &
            ntracer,count_tracer,max_tracers,already_tracked)
 CALL update_tracer(vel(:,:,1),vel(:,:,2),domsize_x,domsize_y, &
      timestep,traced,resolution,dt,count_tracer,max_tracers,track_numbers)
 ! writing tracers to a grid (for time step output)
 tracerfield(:,:,:) = 0
 ntr = 1
 DO WHILE (ntr .LT. count_tracer)
    ix = INT(traced(ntr,1))
    iy = INT(traced(ntr,2))
 
    IF ((ix .GT. 0) .AND. (iy .GT. 0)) THEN
       ! "speed" field
       ! tracerfield(ix,iy,1) = sqrt(traced(ntr,3)**2+traced(ntr,4)**2)
       tracerfield(ix,iy,1) = traced(ntr,8)
       ! "index" field
       tracerfield(ix,iy,2) = traced(ntr, 7)
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
CLOSE(10)
CLOSE(11)
CLOSE(36)
CLOSE(1)
call cpu_time(finish)
print '("Time = ",f10.3," seconds.")',finish-start


END PROGRAM irt_cp_tracking

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                           SUBROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE set_tracer(counter,domsize_x,domsize_y,nneighb,track_numbers,vel,center_of_mass_x, &
     center_of_mass_y,timestep,traced,edge_fraction,& 
     max_no_of_cells,ntracer,count_tracer,max_tracers,already_tracked)

  INTEGER, INTENT(IN)	    :: domsize_x,domsize_y, counter
  REAL, INTENT(INOUT)       :: nneighb(domsize_x,domsize_y)
  REAL, INTENT(IN)       :: track_numbers(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: timestep
  INTEGER, INTENT(IN)       :: max_no_of_cells
  INTEGER, INTENT(IN)       :: max_tracers
  INTEGER, INTENT(INOUT)    :: ntracer
  REAL, INTENT(IN)          :: edge_fraction
  REAL                      :: center_of_mass_x(max_no_of_cells),center_of_mass_y(max_no_of_cells)
  REAL, INTENT(INOUT)       :: traced(max_tracers,9)
  REAL, INTENT(IN)          :: vel(domsize_x,domsize_y,2)
  INTEGER, INTENT(INOUT)    :: count_tracer
  INTEGER                   :: ix, iy
  REAL                      :: vx, vy
  INTEGER, INTENT(INOUT)    :: already_tracked(max_no_of_cells)
  DO iy=1, domsize_y
     DO ix=1, domsize_x
        IF (nneighb(ix,iy) .EQ. 1 ) THEN 
           vx             = .5*(vel(ix,iy,1)+vel(MOD(ix+1-1,domsize_x)+1,iy,1)) ! average vx
           vy             = .5*(vel(ix,iy,2)+vel(ix,MOD(iy+1-1,domsize_y)+1,2)) ! average vy

           IF (already_tracked(INT(track_numbers(ix,iy))) .LT. 300) THEN
              IF (track_numbers(ix,iy) .GT. 0 .AND. track_numbers(ix,iy) .GE. 1) THEN 
                 traced(count_tracer,1) = ix 
                 traced(count_tracer,2) = iy 
                 traced(count_tracer,3) = vx 
                 traced(count_tracer,4) = vy 
                 traced(count_tracer,5) = timestep     ! timestep when tracking begins
                 traced(count_tracer,6) = 0   ! tracer age ???timestep     ! 
                 traced(count_tracer,7) = counter              ! ! OCH ist das die
!track ID Tracer ID
!nein? was fuer ein counter ist das?
                 traced(count_tracer,8) = track_numbers(ix,iy) ! TRACK ID 
                 traced(count_tracer,9) = 1.                   ! has tracer been involved in other precip events
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

SUBROUTINE update_tracer(velx,vely,domsize_x,domsize_y, &
     timestep,traced,resolution, &
     dt, count_tracer, max_tracers,track_numbers)
  INTEGER, INTENT(IN)       :: domsize_x, domsize_y, timestep
  REAL, INTENT(IN)          :: velx(domsize_x,domsize_y),vely(domsize_x,domsize_y)
  REAL, INTENT(IN)          :: resolution
  REAL, INTENT(IN)          :: track_numbers(domsize_x,domsize_y)
  INTEGER, INTENT(IN)       :: count_tracer,max_tracers
  INTEGER                   :: tracer_ts, it
  REAL, INTENT(INOUT)       :: traced(max_tracers,9)
  REAL, INTENT(IN)          :: dt
  REAL                      :: ix_new, iy_new, vx_new, vy_new
  REAL                      :: ix, iy, vx_intp, vy_intp
  INTEGER                   :: ix_round, iy_round, tottracer, start_time, ix_org, iy_org
  INTEGER                   :: ix_round_new, iy_round_new
  REAL                      :: vx_lft, vx_rgt, vx_mn, vy_bot, vy_top, vy_mn, wgt_x, wgt_y
     it=1 
     ! updating previous tracers
     DO WHILE (it .LT. count_tracer)        
       ! determining the first time step of the event
       start_time=traced(it,5)
       ! determining how many timesteps have passed since then
       tracer_ts =timestep-start_time+1  
        
       IF (start_time .GT. 0) THEN
         ! getting the previous positions
         ix=traced(it,1)
         iy=traced(it,2)
         ix_round=MOD(INT(ix)-1+domsize_x,domsize_x)+1
         iy_round=MOD(INT(iy)-1+domsize_y,domsize_y)+1
         ix_m=MOD(INT(ix)-2+domsize_x,domsize_x)+1
         iy_m=MOD(INT(iy)-2+domsize_y,domsize_y)+1

         ! get weights for vertical interpolation 
         wgt_x  = MOD(ix,1.)
         wgt_y  = MOD(iy,1.)

         ! bilinear interpolation
         vx_intp = velx(ix_m,    iy_m    )*(1-wgt_x)*(1-wgt_y) &
                 + velx(ix_m,    iy_round)*(1-wgt_x)*(  wgt_y) &
                 + velx(ix_round,    iy_m)*(  wgt_x)*(1-wgt_y) &
                 + velx(ix_round,iy_round)*(  wgt_x)*(  wgt_y) 
         vy_intp = vely(ix_m,        iy_m)*(1-wgt_x)*(1-wgt_y) &
                 + vely(ix_m,   iy_round )*(1-wgt_x)*(  wgt_y) &
                 + vely(ix_round,    iy_m)*(  wgt_x)*(1-wgt_y) &
                 + vely(ix_round,iy_round)*(  wgt_x)*(  wgt_y) 


! kann vlt weg       IF (ix .GT. 0 .AND. iy .GT. 0) THEN ! bogus now?

         ! get new location as decimal
         ix_new = ix + dt*vx_intpi/resolution
         iy_new = iy + dt*vely(ix_round,iy_round)/resolution

         ! and as gridded values
         ix_round_new= INT(ix_new) !MOD(INT(ix_new)-1+domsize_x,domsize_x)+1
         iy_round_new=INT(iy_new)  !MOD(INT(iy_new)-1+domsize_y,domsize_y)+1

         traced(it,1) = ix_new   ! new position at current timestep
         traced(it,2) = iy_new
         traced(it,3) = vx_intp  ! old interpolated velocity at prev timestep
         traced(it,4) = vy_intp
         traced(it,6) = tracer_ts ! current timestep

         IF (traced(it,9) .eq. 0) THEN    ! if tracer is dead cant get back 
           traced(it,9) = 0  ! stays dead
         ! else (if it was alive) it keeps living as long no other precip
         ! event is above
         ELSEIF (track_numbers(ix_round,iy_round) .eq. traced(it,8) .or. track_numbers(ix_round,iy_round) .le. 0. ) THEN
           traced(it,9) = 1
         ELSE
           traced(it,9) = 0
         ENDIF 

150 format (           2X,I4,   3X,I4,            2X,I4    ,2X,I6, 2X,I4, &
                       2(2X,F10.6),               2(2X,I4)    ,2X,I1)
         WRITE(40,150) timestep,INT(traced(it,5)),tracer_ts,it,INT(traced(it,8)),&
                       ix_new,iy_new,ix_round_new,iy_round_new,INT(traced(it,9))
!<<OCH             
        ENDIF
     it = it + 1
  ENDDO
! WRITE(40,*) ix_new, iy+new, timestep, counter, track_numbers(ix,iy)

  RETURN
 WRITE(*,*) "end f routine" 
END SUBROUTINE update_tracer


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
SUBROUTINE neigh(domsize_x,domsize_y,track_numbers,nneighb,COMx,COMy,input_field,max_no_of_cells)

  IMPLICIT NONE
  INTEGER                :: i,j, imodm, imodp, jmodm,jmodp
  INTEGER, INTENT(IN)    :: domsize_x,domsize_y
  REAL, INTENT(IN)    :: track_numbers(domsize_x,domsize_y)
  REAL, INTENT(INOUT)    :: nneighb(domsize_x,domsize_y) 
  REAL, INTENT(INOUT)    :: input_field(domsize_x,domsize_y)
  INTEGER, INTENT(IN)    :: max_no_of_cells
  REAL, INTENT(INOUT)    :: COMx(max_no_of_cells),COMy(max_no_of_cells) 
  REAL                   :: COUNTER(max_no_of_cells)
  COMx = 0
  COMy = 0
  COUNTER = 0
  nneighb = 0
  DO i =2,domsize_x-1
    DO j =2,domsize_y-1
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
      IF (track_numbers(i,j) .ne. track_numbers(imodp,jmodm)) nneighb(imodp,jmodp) =1
      IF (track_numbers(i,j) .ne. track_numbers(imodm,jmodm)) nneighb(imodm,jmodm) =1  

!       IF ((track_numbers(i,j) .ne. track_numbers(imodp,j) .or. &
!         (track_numbers(i,j) .ne. track_numbers(imodm,j)) .or. &
!         (track_numbers(i,j) .ne. track_numbers(i,jmodp)) .or. &
!         (track_numbers(i,j) .ne. track_numbers(i,jmodm)) ) ) THEN
!          nneighb(i,j) =1
!      END IF     
      IF (track_numbers(i,j) .gt. 0 .and. track_numbers(i,j) .lt. max_no_of_cells ) THEN
        COMx(INT(track_numbers(i,j))) = COMx(INT(track_numbers(i,j)))  + i*input_field(i,j) 
        COMy(INT(track_numbers(i,j))) = COMy(INT(track_numbers(i,j)))  + j*input_field(i,j)
        COUNTER(INT(track_numbers(i,j))) = COUNTER(INT(track_numbers(i,j)))+ input_field(i,j)
      END IF
    END DO
  END DO
  ! dont forget boundaries
  COMx=COMx/COUNTER
  COMy=COMy/COUNTER

END SUBROUTINE neigh

