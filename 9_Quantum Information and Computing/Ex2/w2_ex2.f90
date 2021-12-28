!Week 2 - Exercise 2, Faorlin Tommaso (2021857)
!compile with gfortran -W w2_ex2.f90 my_mat_mul.f90 error_handling.f90 -o w2_ex2.out

PROGRAM test_performance_squared

    USE my_mat_mul
    USE error_handling

    IMPLICIT NONE

    !variables for command line arguments (iterator and quantity)
    INTEGER(4)                                   :: aa, num_commandl_args
    !matrices dimensions
    INTEGER(8)                                   :: rows_mat_A, rows_mat_B, cols_mat_A, cols_mat_B
    !boolean debug variable
    LOGICAL                                      :: debug
    !time tracking
    REAL(4)                                      :: t_i, t_f
    !dynamic declaration for matrices
    REAL(4),        DIMENSION(:,:),  ALLOCATABLE :: mat_A, mat_B, mat_C
    !dynamic declaration of a vector of command line inputs
    CHARACTER(10),  DIMENSION(:),    ALLOCATABLE :: commandl_args

    !header prints
    PRINT *, "Matrix by matrix multiplication benchmarking"
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

    WRITE (*, '(A)', ADVANCE='no') "Insert rows of matrix B: "
    READ(*,'(I8)') rows_mat_B
    WRITE (*, '(A)', ADVANCE='no') "Insert cols of matrix B: "
    READ(*,'(I8)') cols_mat_B

    !PRE-CONDITION: check if rows and columns of matrix B are .GT. 0
    IF(rows_mat_B.LE.0.OR.cols_mat_B.LE.0) THEN
        CALL checkpoint(.TRUE., 'Invalid dimension of matrix B: insert a&
                                & number greater than 0.', 'e')
        STOP
    ENDIF

    !PRE-CONDITION: check if matrix B size is valid (rows and cols less than 10000)
    IF(rows_mat_B.GT.10000.OR.cols_mat_B.GT.10000) THEN
        CALL checkpoint(.TRUE., 'Matrix B is too big: Only sizes below 10000&
                                & are valid.', 'e')
        STOP
    ENDIF

    !PRE-CONDITION: check if matrix multiplication is possible
    IF(cols_mat_A.NE.rows_mat_B) THEN
        CALL checkpoint(.TRUE., 'Dimensions mismatch for matrix multiplication:&
                                & cols_mat_A and rows_mat_B are differing.', 'e')
        STOP
    ENDIF

    PRINT *, "================================================================="

    !if debug, print a checkpoint on matrices dimensions
    CALL checkpoint(debug, 'Rows of matrix A:',    'd', input = rows_mat_A)
    CALL checkpoint(debug, 'Columns of matrix A:', 'd', input = cols_mat_A)
    CALL checkpoint(debug, 'Rows of matrix B:',    'd', input = rows_mat_B)
    CALL checkpoint(debug, 'Columns of matrix B:', 'd', input = cols_mat_B)

    !finally, occupy space in memory for the matrices
    ALLOCATE( mat_A(rows_mat_A, cols_mat_A))
    ALLOCATE( mat_B(rows_mat_B, cols_mat_B))
    ALLOCATE( mat_C(rows_mat_A, cols_mat_B))

    !initialize mat_A and mat_B at random
    CALL RANDOM_NUMBER(mat_A)
    CALL RANDOM_NUMBER(mat_B)

    !results..
    PRINT *, "================================================================="
    PRINT *, 'RESULTS:'
    
    !..for matmul_bycol
    CALL CPU_TIME(t_i)
    mat_C = matmul_bycol(mat_A, mat_B)
    CALL CPU_TIME(t_f)
    PRINT *, "Time for moltiplication by columns", t_f-t_i, 's'
    OPEN (10, FILE = "matmul_bycols.txt", POSITION="append", ACTION="WRITE")
    WRITE (10,'(F10.7)') (t_f-t_i)
    CLOSE (10)

    !POST-CONDITION: check if dimensions are correct in the output matrix (matmul_bycol)
    IF(rows_mat_A.NE.SIZE(mat_C,1).OR.cols_mat_B.NE.SIZE(mat_C,2)) THEN
        CALL checkpoint(.TRUE., 'Dimension mismatch in the output matrix', 'e')
        STOP
    ENDIF

    !..for matmul_byrow
    CALL CPU_TIME(t_i)
    mat_C = matmul_byrow(mat_A, mat_B)
    CALL CPU_TIME(t_f)
    PRINT *, "Time for moltiplication by rows: ", t_f-t_i, 's'
    OPEN (10, FILE="matmul_byrows.txt", POSITION="append", ACTION="WRITE")
    WRITE (10,'(F10.7)') (t_f-t_i)
    CLOSE (10)

    !POST-CONDITION: check if dimensions are correct in the output matrix (matmul_byrow)
    IF(rows_mat_A.NE.SIZE(mat_C,1).OR.cols_mat_B.NE.SIZE(mat_C,2)) THEN
        CALL checkpoint(.TRUE., 'Dimension mismatch in the output matrix', 'e')
        STOP
    ENDIF

    !..for MATMUL
    CALL CPU_TIME(t_i)
    mat_C = MATMUL(mat_A, mat_B)
    CALL CPU_TIME(t_f)
    PRINT *, "Time for moltiplication using 'MATMUL' subroutine: ", t_f-t_i, 's'
    OPEN (10, FILE="matmul_fortran.txt", POSITION="append", ACTION="WRITE")
    WRITE (10,'(F10.7)') (t_f-t_i)

    !POST-CONDITION: check if dimensions are correct in the output matrix (MATMUL)
    IF(rows_mat_A.NE.SIZE(mat_C,1).OR.cols_mat_B.NE.SIZE(mat_C,2)) THEN
        CALL checkpoint(.TRUE., 'Dimension mismatch in the output matrix: ', 'e')
        STOP
    ENDIF

    PRINT *, "================================================================="
    
    CLOSE (10)

    !empty space in memory
    DEALLOCATE(mat_A)
    DEALLOCATE(mat_B)
    DEALLOCATE(mat_C)

END PROGRAM test_performance_squared