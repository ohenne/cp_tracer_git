! OCH March 2018
!
!
!
PROGRAM celloverlay
USE irt_parameters, ONLY: domainsize_x, domainsize_y, lperiodic, &
    resolution, dt

IMPLICIT NONE
INTEGER                 cell_id_var1, cell_ts_var1, &
                        cell_id_var2, cell_ts_var2, &
                        readdummy1, readdummy2, readdummy3, &
                        prev_time, current_time, counter, time_main, time_2nd, &
                        ii, ij, cmerker
INTEGER              :: srv_header_input_main(8), srv_header_input_2nd(8)
REAL, ALLOCATABLE    :: track_numbers1(:,:), track_numbers2(:,:)     ! ID for first and second tracks
REAL                 :: merker(100) 
CHARACTER (len=90)   :: input_filename0  = "irt_tracks_header_sorted.txt"
CHARACTER (len=90)   :: input_filename2  = "irt_tracks_mask_main.srv"
CHARACTER (len=90)   :: input_filename3  = "irt_tracks_mask_2nd.srv"
CHARACTER (len=90)   :: output_filename0 = "cloud_rain_overlap.txt"

! OPEN INPUT FILES
OPEN(10,FILE=trim(input_filename0),FORM='formatted', ACTION='read')
!OPEN(11,FILE=trim(input_filename1),FORM='formatted', ACTION='read')
OPEN(12,FILE=trim(input_filename2),    FORM='unformatted', ACTION='read')
OPEN(13,FILE=trim(input_filename3),    FORM='unformatted', ACTION='read')

! OPEN OUTPUT FILES
OPEN(20,FILE=trim(output_filename0),FORM='formatted', ACTION='write')
501 FORMAT (5X,A3,2(X,A7))
WRITE(20,501) 'ts ','IDmain','IDprev'
500 FORMAT (3(3X,I4))

! ALLOCATE
ALLOCATE(track_numbers1(domainsize_x,domainsize_y))
ALLOCATE(track_numbers2(domainsize_x,domainsize_y))

! INITIALIZE
prev_time = 0
counter   = 1

! read main tracks for the first time
READ (12,END=200) srv_header_input_main
time_main = srv_header_input_main(4)
READ (12) track_numbers1(:,:)

! read secondary tracks for the first time
READ (13,END=200) srv_header_input_2nd
time_2nd = srv_header_input_2nd(4)
READ (13) track_numbers2(:,:)

! read the first timesteps from cloud tracks
! while no rian trakc sexist yet
DO WHILE (time_2nd .lt. time_main)
  READ (13,END=200) srv_header_input_2nd
  time_2nd = srv_header_input_2nd(4)
  READ (13) track_numbers2(:,:)
END DO

! START MAIN LOOP trough every cell track
DO  
  READ(10,*,END=200) cell_id_var1, cell_ts_var1, readdummy1, readdummy2, readdummy3
  !find overlaping cloud
  cmerker = 1
  merker(:) = 0
  DO ii =1,domainsize_x 
   DO ij =1,domainsize_y
      ! find current cellID1 in field and check which cloud track can be found
      IF (track_numbers1(ii,ij) .eq. cell_id_var1 &
         .and. track_numbers2(ii,ij) .gt. 0 ) THEN 
        IF (.not. ANY(merker .eq. track_numbers2(ii,ij)))  THEN
          WRITE(20,500) cell_ts_var1, cell_id_var1, INT(track_numbers2(ii,ij))
          merker(cmerker) = track_numbers2(ii,ij)
        END IF
      END IF
   END DO
  END DO
!  300 CONTINUE
  IF (cell_ts_var1 .gt. prev_time) THEN
    ! read the next time step
    READ (12,END=200) srv_header_input_main
    time_main = srv_header_input_main(4)
    READ (12) track_numbers1(:,:)
    prev_time = cell_ts_var1 

    READ (13,END=200) srv_header_input_2nd
    time_2nd = srv_header_input_2nd(4)
    READ (13) track_numbers2(:,:) 
  END IF
! 
!  counter = counter+1 
!
END DO
200 CONTINUE
END PROGRAM celloverlay
