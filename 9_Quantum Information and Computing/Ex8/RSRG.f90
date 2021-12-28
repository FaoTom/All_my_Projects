!Week 8 -  Faorlin Tommaso (2021857)
!compile with gfortran RSRG.f90  -I/usr/local/Cellar/fftw/3.3.10/include -L/usr/local/Cellar/fftw/3.3.10/lib error_handling.f90 my_mat_mul.f90 qm_utils.f90 -o RSRG.out -lfftw3 -lfftw3f -llapack -lm  -ffree-line-length-0

PROGRAM RSRG

    USE my_mat_mul
    USE error_handling
    USE qm_utils
    USE dc_matrix

    
    IMPLICIT NONE

    !discretization of space variables
    INTEGER(4)                                          :: N, k, tot_it
    REAL(8)                                             :: lambda, farg, delta, gs, est_gs, sys_size
    CHARACTER(10),  DIMENSION(:),    ALLOCATABLE        :: commandl_args
    LOGICAL                                             :: custom, debug
    INTEGER(4)                                          :: aa, num_commandl_args, iarg
    DOUBLE COMPLEX, DIMENSION(:,:),  ALLOCATABLE        :: H_N_T
    REAL(8),  DIMENSION(:),    ALLOCATABLE              :: eigvls
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE        :: H_2N, P, eigvct, HabA_T, HabB_T, HabA_2N, HabB_2N
    INTEGER                                            :: info, lwork
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE          :: work(:)
    REAL(8),    DIMENSION(:), ALLOCATABLE              :: rwork(:)
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

    !____________________________________________________________________________________
    !____________________________________________________________________________________
    !____________________________________________________________________________________

    !define number of command line arguments
    num_commandl_args = COMMAND_ARGUMENT_COUNT()
    !allocate in memory the vector
    ALLOCATE(commandl_args(num_commandl_args))  
    

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
        N      = 2
        lambda = 2.0
    ELSE
        !having more than one argument is not accepted at the moment
        CALL checkpoint(.TRUE., 'Invalid number of arguments in the command line', 'e')
        STOP 
    ENDIF

    !____________________________________________________________________________________
    !___________________________RSRG_____________________________________________________
    !____________________________________________________________________________________

    ALLOCATE(H_2N(2**(2*N), 2**(2*N)))
    ALLOCATE(H_N_T(2**N, 2**N))
    ALLOCATE(eigvls(2**(2*N)))
    
    !Hamiltonian operator on the doubled system H=2**2Nx2**2N
    ! ALLOCATE(HabA_T(2**(N), 2**(N)))          !interaction term on subsystem AA
    ! ALLOCATE(HabB_T(2**(N), 2**(N)))          !interaction term on subsystem B
    ! ALLOCATE(HabA_2N(2**(2*N), 2**(2*N)))          !interaction term on subsystem AA
    ! ALLOCATE(HabB_2N(2**(2*N), 2**(2*N)))          !interaction term on subsystem B
    
    k = 2**N
    delta = 1e-12
    tot_it = 0
    sys_size = DBLE(N)

    ALLOCATE(P(2**(2*N),k))

    H_N_T = TFI_ham_init(N, lambda)
    HabA_T = tp(id_mat_init(N-1), sigma_x())
    HabB_T = tp(sigma_x(), id_mat_init(N-1))

    gs = -1.0d0

    DO WHILE(ABS(gs-est_gs)>delta)

        gs = est_gs

        H_2N = tp(H_N_T, id_mat_init(N)) + tp(id_mat_init(N), H_N_T) + tp(HabA_T, HabB_T)
        HabA_2N = tp(HabA_T, id_mat_init(N))
        HabB_2N = tp(id_mat_init(N),HabB_T)

        eigvct = H_2N

        lwork = MAX(1,2*(2**(2*N))-1)
        ALLOCATE(work(MAX(1, lwork)))
        ALLOCATE(rwork(MAX(1, 3*(2**(2*N))-2)))
        
        CALL ZHEEV('V','U', 2**(2*N), eigvct, 2**(2*N), eigvls, work, lwork, rwork, info)

        DEALLOCATE(work)
        DEALLOCATE(rwork)
        sys_size = sys_size * 2
        
        est_gs = eigvls(1)/DBLE(sys_size)
        
        DO ii=1,SIZE(P,2)
            P(:,ii) = eigvct(:,ii)
        ENDDO

        H_N_T  =  MATMUL(MATMUL(TRANSPOSE(P),H_2N),P)
        HabA_T =  MATMUL(MATMUL(TRANSPOSE(P),HabA_2N),P)
        HabB_T =  MATMUL(MATMUL(TRANSPOSE(P),HabB_2N),P)

        tot_it = tot_it + 1
    
    ENDDO

    PRINT * , tot_it, est_gs

END PROGRAM RSRG