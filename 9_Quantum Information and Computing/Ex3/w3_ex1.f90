!Week 3 - Exercise 1, Faorlin Tommaso (2021857)

PROGRAM test_performance_squared

USE my_mat_mul
USE error_handling

IMPLICIT NONE

!matrices dimensions
INTEGER(8)                                   :: rows_mat_A, rows_mat_B, cols_mat_A, cols_mat_B
!boolean debug variable
LOGICAL                                      :: debug
!time tracking
REAL(4)                                      :: t_i, t_f
!dynamic declaration for matrices
REAL(4),        DIMENSION(:,:),  ALLOCATABLE :: mat_A, mat_B, mat_C
!dynamic declaration of a vector of command line inputs
CHARACTER(10)                                :: command_line_arg

!set .TRUE. to enter debug mode
debug = .FALSE.

CALL GET_COMMAND_ARGUMENT(1, command_line_arg)
READ(command_line_arg, FMT="(I8)") rows_mat_A

cols_mat_A = rows_mat_A
rows_mat_B = rows_mat_A
cols_mat_B = rows_mat_A

!if debug, print a checkpoint on matrices dimensions
CALL checkpoint(debug, 'Rows of matrix A:',    'd', input = rows_mat_A)
CALL checkpoint(debug, 'Columns of matrix A:', 'd', input = cols_mat_A)
CALL checkpoint(debug, 'Rows of matrix B:',    'd', input = rows_mat_B)
CALL checkpoint(debug, 'Columns of matrix B:', 'd', input = cols_mat_B)

!PRE-CONDITION: check if rows and columns of matrix A are .GT. 0
IF(rows_mat_A.LE.0.OR.cols_mat_A.LE.0) THEN
    CALL checkpoint(debug, 'Invalid dimension of matrix A: insert a&
                            & number greater than 0.', 'e')
    STOP
ENDIF

!PRE-CONDITION: check if rows and columns of matrix B are .GT. 0
IF(rows_mat_B.LE.0.OR.cols_mat_B.LE.0) THEN
    CALL checkpoint(debug, 'Invalid dimension of matrix B: insert a&
                            & number greater than 0.', 'e')
    STOP
ENDIF

!PRE-CONDITION: check if matrix A size is valid (rows and cols less than 10000)
IF(rows_mat_A.GT.10000.OR.cols_mat_A.GT.10000) THEN
    CALL checkpoint(debug, 'Matrix A is too big: Only sizes below 10000&
                            & are valid.', 'e')
    STOP
ENDIF

!PRE-CONDITION: check if matrix B size is valid (rows and cols less than 10000)
IF(rows_mat_B.GT.10000.OR.cols_mat_B.GT.10000) THEN
    CALL checkpoint(debug, 'Matrix B is too big: Only sizes below 10000&
                            & are valid.', 'e')
    STOP
ENDIF

!PRE-CONDITION: check if matrix multiplication is possible
IF(cols_mat_A.NE.rows_mat_B) THEN
    CALL checkpoint(debug, 'Dimensions mismatch for matrix multiplication:&
                            & cols_mat_A and rows_mat_B are differing.', 'e')
    STOP
ENDIF

!finally, occupy space in memory for the matrices
ALLOCATE( mat_A(rows_mat_A, cols_mat_A))
ALLOCATE( mat_B(rows_mat_B, cols_mat_B))
ALLOCATE( mat_C(rows_mat_A, cols_mat_B))

!initialize mat_A and mat_B at random
CALL RANDOM_NUMBER(mat_A)
CALL RANDOM_NUMBER(mat_B)

!..for matmul_bycol
CALL CPU_TIME(t_i)
mat_C = matmul_bycol(mat_A, mat_B)
CALL CPU_TIME(t_f)
!WRITE (*, '(F10.7)', ADVANCE='no') (t_f-t_i)
WRITE (*, '(F10.7)') (t_f-t_i)

!POST-CONDITION: check if dimensions are correct in the output matrix (matmul_bycol)
IF(rows_mat_A.NE.SIZE(mat_C,1).OR.cols_mat_B.NE.SIZE(mat_C,2)) THEN
    CALL checkpoint(debug, 'Dimension mismatch in the output matrix', 'e')
    STOP
ENDIF

!..for matmul_byrow
CALL CPU_TIME(t_i)
mat_C = matmul_byrow(mat_A, mat_B)
CALL CPU_TIME(t_f)
!WRITE (*, '(F10.7)', ADVANCE='no') (t_f-t_i)
WRITE (*, '(F10.7)') (t_f-t_i)

!POST-CONDITION: check if dimensions are correct in the output matrix (matmul_byrow)
IF(rows_mat_A.NE.SIZE(mat_C,1).OR.cols_mat_B.NE.SIZE(mat_C,2)) THEN
    CALL checkpoint(debug, 'Dimension mismatch in the output matrix', 'e')
    STOP
ENDIF

!..for MATMUL
CALL CPU_TIME(t_i)
mat_C = MATMUL(mat_A, mat_B)
CALL CPU_TIME(t_f)
!WRITE (*, '(F10.7)', ADVANCE='no') (t_f-t_i)
WRITE (*, '(F10.7)') (t_f-t_i)

!POST-CONDITION: check if dimensions are correct in the output matrix (MATMUL)
IF(rows_mat_A.NE.SIZE(mat_C,1).OR.cols_mat_B.NE.SIZE(mat_C,2)) THEN
    CALL checkpoint(debug, 'Dimension mismatch in the output matrix', 'e')
    STOP
ENDIF

!empty space in memory
DEALLOCATE(mat_A)
DEALLOCATE(mat_B)
DEALLOCATE(mat_C)

END PROGRAM test_performance_squared