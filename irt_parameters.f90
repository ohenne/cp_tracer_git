MODULE irt_parameters

INTEGER, PARAMETER    :: domainsize_x = 320
INTEGER, PARAMETER    :: domainsize_y = 320
REAL,    PARAMETER    :: resolution = 200 ! resolution in m
REAL,    PARAMETER    :: dt = 300 ! time step in sec

LOGICAL, PARAMETER    :: lperiodic = .TRUE.

INTEGER, PARAMETER    :: n_fields = 0   ! number of additional averaging fields


! bins of coarse velocity field
INTEGER, PARAMETER    :: time_steps = 267    ! total number of timesteps
INTEGER, PARAMETER    :: nt_bins = 10
INTEGER, PARAMETER    :: nx_bins = 2
INTEGER, PARAMETER    :: ny_bins = 2

!REAL, PARAMETER       :: cutoff = 0.5 !1.0            ! for rain intensity
!REAL, PARAMETER       :: minimum_size = 100.       ! events smaller than that will be sorted out
! for cp_tracking: We are interested in deep convective cps. Thus precipi event
! causing cp should be larger than RB scale (2km) ~ 10x10 gps 
REAL, PARAMETER       :: threshold_ratio=1      ! Choose value between 0.0 and 1.0

REAL, PARAMETER       :: max_velocity = 10.   ! adjust acordingly
                                              ! velocities>max_velocity will be ignored to remove outliers
! define a minimal number of cells required for a coarse grained coordinate to be evaluated 
! if there are less, missing value will be assigned to that coarse cell
REAL, PARAMETER       :: edge_fraction = .01  ! the estimated fraction of the domain located on the edge of cold pools
INTEGER, PARAMETER    :: min_cells = 10
INTEGER, PARAMETER    :: tracer_steps = 20
INTEGER, PARAMETER    :: max_tfields = 20     ! maximum number of distinct tracer fields to keep track of

INTEGER, PARAMETER    :: max_no_of_cells=10000  ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_no_of_tracks=10000    ! buffer size, increase if necessary
INTEGER, PARAMETER    :: max_length_of_track=267  ! buffer size, increase if necessary

REAL, PARAMETER       :: miss=-9999.           ! value<miss ==> missing_value
INTEGER, PARAMETER    :: max_tracers=500000     ! max number of tracers to track


END MODULE irt_parameters
