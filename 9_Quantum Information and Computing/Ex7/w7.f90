!Week 7 -  Faorlin Tommaso (2021857)
!compile with gfortran w7.f90  -I/usr/local/Cellar/fftw/3.3.10/include -L/usr/local/Cellar/fftw/3.3.10/lib error_handling.f90 my_mat_mul.f90 qm_utils.f90 -o w7.out -lfftw3 -lfftw3f -llapack -lm  -ffree-line-length-0

PROGRAM week_7

    USE my_mat_mul
    USE error_handling
    USE qm_utils
    USE dc_matrix

    
    IMPLICIT NONE

    !discretization of space variables
    INTEGER(4)                                          :: N
    REAL(8)                                             :: lambda, farg, t_i, t_f
    CHARACTER(10),  DIMENSION(:),    ALLOCATABLE        :: commandl_args
    LOGICAL                                             :: custom, debug
    INTEGER(4)                                          :: aa, num_commandl_args, iarg
    DOUBLE COMPLEX, DIMENSION(:,:),  ALLOCATABLE        :: TFI_ham
    REAL(8),  DIMENSION(:),    ALLOCATABLE              :: eigs
    TYPE(mymatrix)                                      :: TFI_ham_dcmatrix

    !time tracking   

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
    ALLOCATE(eigs(N))

    !loop on the number of cl arguments and fill the vector
    DO aa = 1, num_commandl_args
        CALL GET_COMMAND_ARGUMENT(aa, commandl_args(aa))
    END DO

    !arguments are passed
    IF(num_commandl_args.EQ.3) THEN
        !if argument is -custom insert manually discretization parameters
        IF(commandl_args(1).EQ.'-custom') THEN
            custom = .TRUE.
            READ(commandl_args(2),"(I8)") iarg
            N = iarg
            READ(commandl_args(3),"(F10.2)") farg
            lambda = farg
        ELSE
            !other arguments for the moment are not accepted
            CALL checkpoint(.TRUE., 'Invalid argument inserted', 'e')
            STOP
        ENDIF
    !no arguments are passed, so no debug mode
    ELSE IF (num_commandl_args.EQ.0) THEN
        N      = 12
        lambda = 2.0
    ELSE
        !having more than one argument is not accepted at the moment
        CALL checkpoint(.TRUE., 'Invalid number of arguments in the command line', 'e')
        STOP 
    ENDIF

    ALLOCATE(TFI_ham(2**N, 2**N))

    CALL CPU_TIME(t_i)
    TFI_ham = TFI_ham_init(N, lambda)
    CALL CPU_TIME(t_f)
    WRITE (*, '(F14.7)') (t_f-t_i)
    
    TFI_ham_dcmatrix%m_dims(1)= 2**N
    TFI_ham_dcmatrix%m_dims(2)= 2**N
    TFI_ham_dcmatrix%m_el = TFI_ham
    PRINT *, size(TFI_ham_dcmatrix%m_el)
    CALL CPU_TIME(t_i)
    eigs = eigenvalues(TFI_ham_dcmatrix)
    CALL CPU_TIME(t_f)
    WRITE (*, '(F14.7)') (t_f-t_i)

    !PRINT *, eigs(1), eigs(2), eigs(3), eigs(4)

    !PRINT *, eigs(2)-eigs(1)

END PROGRAM week_7