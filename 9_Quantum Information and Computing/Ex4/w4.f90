!Week 4 -  Faorlin Tommaso (2021857)
!compile with gfortran -Wall dc_matrix.f90 error_handling.f90 my_mat_mul.f90 qm_utils.f90 w4.f90 -o w4.out -llapack

PROGRAM harmonic_osc
    !This program solves the time independent Schrödinger equation
    !for a Harmonic oscillator

    USE my_mat_mul
    USE error_handling
    USE qm_utils
    
    IMPLICIT NONE

    !discretization variables
    REAL(8)                                      :: Lmin, Lmax, m, omega
    REAL(8),    DIMENSION(:),   ALLOCATABLE      :: eigs
    REAL(8),    DIMENSION(:,:), ALLOCATABLE      :: ham_op
    INTEGER(4)                                   :: N, aa, num_commandl_args, iarg
    CHARACTER(10),  DIMENSION(:),    ALLOCATABLE :: commandl_args
    !boolean debug variable
    LOGICAL                                      :: scripting
    !variables for string formatting
    CHARACTER(len=1024) :: filename

    !header prints
    PRINT *, "================================================================="
    PRINT *, "-> Numerical solution of the quantum harmonic oscillator <-"
    PRINT *, ''
    PRINT *, "Author:  Tommaso Faorlin (2021857)"
    PRINT *, ''
    PRINT *, "Run with '-custom N Lmin Lmax omega' for non-default values"
    
    !define number of command line arguments
    num_commandl_args = COMMAND_ARGUMENT_COUNT()
    !allocate in memory the vector
    ALLOCATE(commandl_args(num_commandl_args))  

    !loop on the number of cl arguments and fill the vector
    DO aa = 1, num_commandl_args
        CALL GET_COMMAND_ARGUMENT(aa, commandl_args(aa))
    END DO

    !arguments are passed
    IF(num_commandl_args.EQ.5) THEN
        !if argument is -custom insert manually discretization parameters
        IF(commandl_args(1).EQ.'-custom') THEN
            scripting = .TRUE.
            READ(commandl_args(2),"(I8)") iarg
            N = iarg
            READ(commandl_args(3),"(I8)") iarg
            Lmin = iarg
            READ(commandl_args(4),"(I8)") iarg
            Lmax = iarg
            READ(commandl_args(5),"(I8)") iarg
            omega = iarg
        ELSE
            !other arguments for the moment are not accepted
            CALL checkpoint(.TRUE., 'Invalid argument inserted', 'e')
            STOP
        ENDIF
    !no arguments are passed: set values to default
    ELSE IF (num_commandl_args.EQ.0) THEN
        scripting = .FALSE.
        N     = 2500
        Lmin  = -10
        Lmax  = +10
        omega = 1
    ELSE
        !having more than one argument is not accepted at the moment
        CALL checkpoint(.TRUE., 'Invalid number of arguments in the command line', 'e')
        STOP 
    ENDIF

    !set mass and angular frequency equal to 1
    m = 1

    !simulation parameters print
    PRINT *, "================================================================="
    PRINT *, "SIMULATION PARAMETERS: "
    PRINT *, "-N:", N
    PRINT *, "-Lmax:  ", Lmin
    PRINT *, "-Lmax:  ", Lmax
    PRINT *, "-dx:    ",(Lmax-Lmin)/N
    PRINT *, "-omega: ", omega
    PRINT *, "-m:     ", m
    PRINT *, "================================================================="

    !N+1 because after dividing the interval in N steps we obtain N+1 points
    ALLOCATE(ham_op(N+1,N+1))
    ALLOCATE(eigs(N+1))
    
    ham_op = harmosc_H_op_init(N, Lmin, Lmax, m, omega)
    print *, ham_op
    CALL diagonalize_H_op(ham_op, eigs)

    WRITE(filename, "(ai4af4.1af4.1a)") "data/eigenfunc_N", N, "_Lmax", Lmax, "_omega", omega, ".dat"
    OPEN(UNIT = 2, FILE = filename, ACTION = 'write', STATUS = 'replace')
    DO aa = 1, SIZE(ham_op,1)
        WRITE(2, *)  ham_op(aa,10)
    ENDDO
    CLOSE(2)

    WRITE(filename, "(ai4af4.1af4.1a)") "data/eigenvalues_N", N, "_Lmax", Lmax, "_omega", omega, ".dat"
    OPEN(UNIT = 2, FILE = filename, ACTION = 'write', STATUS = 'replace')
    DO aa = 1, SIZE(eigs)
        WRITE(2, *)  eigs(aa)
    ENDDO
    CLOSE(2)

END PROGRAM harmonic_osc
