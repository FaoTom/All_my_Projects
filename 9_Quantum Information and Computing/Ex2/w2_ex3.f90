!Week 2 - Exercise 3, Faorlin Tommaso (2021857)
!compile with gfortran -W w2_ex3.f90 dc_matrix.f90 error_handling.f90 -o w2_ex3.out

PROGRAM test_program

    USE dc_matrix
    USE error_handling

    IMPLICIT NONE

    !boolean debug variable
    LOGICAL                                      :: debug
    !variables for command line arguments (iterator and quantity)
    INTEGER(4)                                   :: aa, num_commandl_args
    !dynamic declaration of a vector of command line inputs
    CHARACTER(10),  DIMENSION(:),    ALLOCATABLE :: commandl_args
    !initialize test and test adjoint matrices
    TYPE(mymatrix)                               :: test_matrix, test_matrix_adj
    !dimensions passed via command line
    INTEGER(4)                                   :: rows_mat_A, cols_mat_A
    !dimensions array for mymatrix type
    INTEGER(4), DIMENSION(2) :: m_dims

    !header prints
    PRINT *, "Operations on mymatrix user-defined type"
    PRINT *, "Author:  Tommaso Faorlin (2021857)"
    PRINT *, "Add -debug to the executable to enable debug mode"

    !define number of command line arguments
    num_commandl_args = COMMAND_ARGUMENT_COUNT()
    !allocate in memory the vector
    ALLOCATE(commandl_args(num_commandl_args))  

    !loop on the number of cl arguments and fill the vector
    DO aa = 1, num_commandl_args
        CALL GET_COMMAND_ARGUMENT(aa, commandl_args(aa))
    END DO

    !arguments are passed
    IF(num_commandl_args.EQ.1) THEN
        !if argument is -debug then set debug mode
        IF(commandl_args(1).EQ.'-debug') THEN
            
            debug = .TRUE. 
        ELSE
            !other arguments for the moment are not accepted
            CALL checkpoint(.TRUE., 'Invalid argument inserted', 'e')
            STOP
        ENDIF
    !no arguments are passed, so no debug mode
    ELSE IF (num_commandl_args.EQ.0) THEN
        debug = .FALSE.
    ELSE
        !having more than one argument is not accepted at the moment
        CALL checkpoint(.TRUE., 'Invalid number of arguments in the command line', 'e')
        STOP 
    ENDIF

    PRINT *, "================================================================="
    PRINT *, "-> Insert the dimensions of the two matrices <-"
    WRITE (*, '(A)', ADVANCE='no') "Insert rows of matrix A: "
    READ(*,'(I8)') rows_mat_A
    WRITE (*, '(A)', ADVANCE='no') "Insert cols of matrix A: "
    READ(*,'(I8)') cols_mat_A

    !PRE-CONDITION: check if rows and columns of matrix A are .GT. 0
    IF(rows_mat_A.LE.0.OR.cols_mat_A.LE.0) THEN
        CALL checkpoint(.TRUE., 'Invalid dimension of matrix A: insert a&
                                & number greater than 0.', 'e')
        STOP
    ENDIF

    !PRE-CONDITION: check if matrix A size is valid (rows and cols less than 10000)
    IF(rows_mat_A.GT.10000.OR.cols_mat_A.GT.10000) THEN
        CALL checkpoint(.TRUE., 'Matrix A is too big: Only sizes below 10000&
                                & are valid.', 'e')
        STOP
    ENDIF

    !set fixed dimension of the matrices
    m_dims(1) = rows_mat_A
    m_dims(2) = cols_mat_A

    test_matrix = random_init(m_dims)    
    test_matrix_adj = .ADJ.test_matrix

    CALL print_to_file(test_matrix,     20, 'test.txt'   )
    CALL print_to_file(test_matrix_adj, 20, 'test_adj.txt')

END PROGRAM test_program