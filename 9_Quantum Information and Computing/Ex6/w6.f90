!Week 5 -  Faorlin Tommaso (2021857)
!compile with gfortran w6.f90  -I/usr/local/Cellar/fftw/3.3.10/include -L/usr/local/Cellar/fftw/3.3.10/lib error_handling.f90 my_mat_mul.f90 qm_utils.f90 -o w6.out -lfftw3 -lfftw3f -llapack -lm  -ffree-line-length-0

PROGRAM week_6

    USE my_mat_mul
    USE error_handling
    USE qm_utils
    
    IMPLICIT NONE

    !discretization of space variables
    INTEGER(4)                                          :: local_dim, sys_dim
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE           :: wf
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE         :: dmat, dred
    CHARACTER(10),  DIMENSION(:),    ALLOCATABLE        :: commandl_args
    LOGICAL                                             :: custom, debug
    LOGICAL                                             :: trace
    INTEGER(4)                                          :: aa, num_commandl_args, iarg, red_sub
    !time tracking
    REAL(4)                                             :: t_i, t_f    

    debug = .FALSE.
    
    !header prints
    IF (debug.EQV..TRUE.) THEN
        PRINT *, "================================================================="
        PRINT *, "-> Initialization of generic and separable quantum states <-"
        PRINT *, ''
        PRINT *, "Author:  Tommaso Faorlin (2021857)"
        PRINT *, ''
        PRINT *, "Run with '-custom system_dimension local_dimension -trace reduced_subsystem' for non-default values"
    ENDIF

        !define number of command line arguments
    num_commandl_args = COMMAND_ARGUMENT_COUNT()
    !allocate in memory the vector
    ALLOCATE(commandl_args(num_commandl_args))  

    !loop on the number of cl arguments and fill the vector
    DO aa = 1, num_commandl_args
        CALL GET_COMMAND_ARGUMENT(aa, commandl_args(aa))
    END DO

    !arguments are passed
    IF(num_commandl_args.EQ.3.OR.num_commandl_args.EQ.5) THEN
        !if argument is -custom insert manually discretization parameters
        IF(commandl_args(1).EQ.'-custom') THEN
            custom = .TRUE.
            READ(commandl_args(2),"(I8)") iarg
            local_dim = iarg
            READ(commandl_args(3),"(I8)") iarg
            sys_dim = iarg
            IF (commandl_args(4).EQ.'-trace') THEN
                trace = .TRUE.
                READ(commandl_args(5),"(I8)") iarg
                red_sub = iarg
            ELSE 
                trace = .FALSE.
            ENDIF
        ELSE
            !other arguments for the moment are not accepted
            CALL checkpoint(.TRUE., 'Invalid argument inserted', 'e')
            STOP
        ENDIF
    !no arguments are passed, so no debug mode
    ELSE IF (num_commandl_args.EQ.0) THEN
        custom    = .FALSE.
        local_dim = 2
        sys_dim   = 2
        red_sub   = 1
        trace = .FALSE.
    ELSE
        !having more than one argument is not accepted at the moment
        CALL checkpoint(.TRUE., 'Invalid number of arguments in the command line', 'e')
        STOP 
    ENDIF
    
    CALL CPU_TIME(t_i)
    wf = sep_wf_init(local_dim, sys_dim)
    CALL CPU_TIME(t_f)
    !WRITE (*, '(F10.7)') (t_f-t_i)

    ! CALL CPU_TIME(t_i)
    ! wf = gen_wf_init(local_dim, sys_dim)
    ! CALL CPU_TIME(t_f)
    ! WRITE (*, '(F10.7)') (t_f-t_i)

    IF(trace.EQV..TRUE.) THEN
        dmat = dmatrix(wf)
        dred = p_trace(local_dim,dmat, red_sub)
    ENDIF

    PRINT *, dmat

END PROGRAM week_6