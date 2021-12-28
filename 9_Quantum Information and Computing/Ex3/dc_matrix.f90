!some other functions have been added to the original module used in Week 2 - Exercise 3

MODULE dc_matrix

        ! MODULE dc_matrix
        !  -> This module dc_matrix (doublecomplex_matrix) can be used
        !     to deal comfortably with matrices with double precision complex
        !     numbers as entries.
        !
        ! =================================================================
        !
        ! FUNCTIONS AND SUBROUTINES:
        ! - FUNCTION   random_init(input_dims)       RESULT(matrix_out)
        ! - FUNCTION   adjoint(matrix_in)            RESULT(matrix_out)
        ! - FUNCTION   trace(matrix_in)              RESULT(trace_value)
        ! - FUNCTION   determinant(matrix_in)        RESULT(determinant_value)
        ! - SUBROUTINE print_to_file(matrix_in, unit, filename)
        !
        ! ** NEW WRT PREVIOUS VERSION OF dcmatrix.f90 **
        ! - FUNCTION   hermitian_init(input_dims)     RESULT(matrix_out)
        ! - FUNCTION   diagonal_init(input_dims)      RESULT(matrix_out)
        ! - FUNCTION   eigenvalues(matrix_in)         RESULT(vector_output)
        ! - FUNCTION   norm_eigenvalues_spacing(eigs) RESULT(vector_output)

    use error_handling

    TYPE mymatrix
        !matrix elements
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE  :: m_el
        !matrix dimensions
        INTEGER(4),  DIMENSION(2)                :: m_dims
        !trace and determinant variables
        COMPLEX*16                               :: m_trace, m_det
    END TYPE mymatrix 

    !define the interface operator for trace of a mymatrix TYPE
    INTERFACE OPERATOR(.TR.)
        MODULE PROCEDURE trace
    END INTERFACE
    !define the interface operator for adjoint of a mymatrix TYPE
    INTERFACE OPERATOR(.ADJ.)
        MODULE PROCEDURE adjoint
    END INTERFACE
    !define the interface operator for determinant of a mymatrix TYPE
    INTERFACE OPERATOR(.DET.)
        MODULE PROCEDURE determinant
    END INTERFACE

    CONTAINS

    FUNCTION random_init(input_dims) RESULT(matrix_out)

        ! FUNCTION random_init
        !  -> This function randomly initializes a mymatrix TYPE 
        !     of given dimensions with COMPLEX*16 numbers and 'if 
        !     possible' computes also its determinant and trace.
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - input_dims: vector containing the number of rows and of columns
        ! - ii, jj:     integer iterators 
        ! - matrix_out: mymatrix TYPE in output 
        !

        INTEGER(4), DIMENSION(2) :: input_dims
        INTEGER(4)               :: ii, jj
        TYPE(mymatrix)           :: matrix_out

        IF ((input_dims(1) > 0).AND.(input_dims(2) > 0)) THEN

            !initialize the dimensions of matrix_out
            matrix_out%m_dims(1) = input_dims(1)
            matrix_out%m_dims(2) = input_dims(2)

            !allocate space in memory for matrix elements
            ALLOCATE(matrix_out%m_el(matrix_out%m_dims(1),matrix_out%m_dims(2)))

            !loop over the entries of matrix_out and assign random values
            DO ii = 1,input_dims(1)
                DO jj = 1,input_dims(2)
                    matrix_out%m_el(ii,jj) = COMPLEX(RAND(),RAND())
                ENDDO 
            ENDDO  
        ENDIF 

        !compute the trace
        matrix_out%m_trace = trace(matrix_out)
        !compute the determinant
        matrix_out%m_det = determinant(matrix_out)
        
    END FUNCTION random_init

    FUNCTION hermitian_init(input_dims) RESULT(matrix_out)

        ! FUNCTION hermitian_init
        !  -> 
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - input_dims: vector containing the number of rows and of columns
        !

        INTEGER(4), DIMENSION(2) :: input_dims
        INTEGER(4)               :: ii, jj
        TYPE(mymatrix)           :: matrix_out

        !initialize the dimensions of matrix_out
        matrix_out%m_dims(1) = input_dims(1)
        matrix_out%m_dims(2) = input_dims(2)

        !allocate space in memory for matrix elements
        ALLOCATE(matrix_out%m_el(matrix_out%m_dims(1),matrix_out%m_dims(2)))

        DO jj = 1, input_dims(1)
            DO ii = 1, input_dims(2)
                IF (ii.EQ.jj) THEN
                    matrix_out%m_el(ii, jj) = COMPLEX(-1+2*RAND(),0d0)
                ELSE IF (jj.GT.ii) THEN
                    matrix_out%m_el(ii, jj) = COMPLEX(-1+2*RAND(),-1+2*RAND())
                    matrix_out%m_el(jj, ii) = CONJG(matrix_out%m_el(ii, jj))
                ELSE 
                    CONTINUE
                ENDIF           
            ENDDO
        ENDDO


        !compute the trace
        matrix_out%m_trace = trace(matrix_out)
        !compute the determinant
        matrix_out%m_det = determinant(matrix_out)
        
    END FUNCTION hermitian_init

    FUNCTION diagonal_init(input_dims) RESULT(matrix_out)

        ! FUNCTION hermitian_init
        !  -> 
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - input_dims: vector containing the number of rows and of columns
        !

        INTEGER(4), DIMENSION(2) :: input_dims
        INTEGER(4)               :: ii, jj
        TYPE(mymatrix)           :: matrix_out

        !initialize the dimensions of matrix_out
        matrix_out%m_dims(1) = input_dims(1)
        matrix_out%m_dims(2) = input_dims(2)

        !allocate space in memory for matrix elements
        ALLOCATE(matrix_out%m_el(matrix_out%m_dims(1),matrix_out%m_dims(2)))

        matrix_out%m_el = 0

        DO ii = 1, input_dims(1)
            DO jj = ii, input_dims(2)
                IF(ii.EQ.jj) THEN
                    matrix_out%m_el(ii, jj) = COMPLEX(-1+2*RAND(),0d0)
                ENDIF
            ENDDO
        ENDDO

        !compute the trace
        matrix_out%m_trace = trace(matrix_out)
        !compute the determinant
        matrix_out%m_det = determinant(matrix_out)
        
    END FUNCTION diagonal_init


    FUNCTION adjoint(matrix_in) RESULT(matrix_out)

        ! FUNCTION adjoint
        !  -> This function computes the adjoint of a mymatrix TYPE.
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - matrix_in:  input mymatrix
        ! - ii, jj:     integer iterators 
        ! - matrix_out: adjoint of matrix_in
        !

        TYPE(mymatrix), INTENT(IN) :: matrix_in
        INTEGER(4)                 :: ii, jj
        TYPE(mymatrix)             :: matrix_out

        !initialize the dimensions of matrix_out
        matrix_out%m_dims(1) = matrix_in%m_dims(2)
        matrix_out%m_dims(2) = matrix_in%m_dims(1)

        !allocate space in memory for matrix elements
        ALLOCATE(matrix_out%m_el(matrix_out%m_dims(1), matrix_out%m_dims(2)))

        !loop over matrix_out(ii,jj) and assign matrix_out(jj,ii) to each
        DO ii = 1, matrix_out%m_dims(1)
            DO jj = 1, matrix_out%m_dims(2)
                matrix_out%m_el(ii,jj) = CONJG(matrix_in%m_el(jj,ii))
            ENDDO
        ENDDO
        
        ! or, in just one line:
        ! matrix_out%m_el = CONJG(TRANSPOSE(matrix_in%m_el)) 
        
        !compute the trace
        matrix_out%m_trace = .TR.matrix_out
        !compute the determinant
        matrix_out%m_det = .DET.matrix_out

    END FUNCTION adjoint

    FUNCTION trace(matrix_in) RESULT(trace_value)

        ! FUNCTION trace
        !  -> This function computes the trace of a mymatrix TYPE.
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - matrix_in:   input mymatrix
        ! - trace_value: trace of matrix_in
        ! - ii:          integer iterator
        !
        ! ERROR HANDLING:
        ! This function does not compute the trace if matrix_in is rectangular 

        TYPE(mymatrix), INTENT(IN) :: matrix_in
        COMPLEX*16                :: trace_value
        INTEGER(4)                 :: ii

        IF (matrix_in%m_dims(1).EQ.matrix_in%m_dims(2)) THEN
            !initialize to 0+0i the output value
            trace_value = COMPLEX(0d0,0d0) 
            !iterate and sum over diagonal elements
            DO ii = 1, matrix_in%m_dims(1)
                trace_value = trace_value+matrix_in%m_el(ii,ii)      
            ENDDO
        ELSE 
            CALL checkpoint(.TRUE., 'Rectangular matrix: cannot calculate the trace (set to -1)', 'e')
            trace_value = COMPLEX(-1,-1) 
        ENDIF

    END FUNCTION

    FUNCTION determinant(matrix_in) RESULT(determinant_value)

        ! FUNCTION determinant [TODO]
        !  -> This function computes the determinant of a mymatrix TYPE.
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - matrix_in:          input mymatrix
        ! - determinant_value:  determinant of matrix_in
        !
        ! ERROR HANDLING:
        ! This function does not compute the determinant if matrix_in is rectangular 

        TYPE(mymatrix), INTENT(IN) :: matrix_in
        COMPLEX*16                 :: determinant_value

        IF (matrix_in%m_dims(1).EQ.matrix_in%m_dims(2)) THEN
            !TODO
            determinant_value = (0,0)
        ELSE
            CALL checkpoint(.TRUE., 'Rectangular matrix: cannot calculate the determinant (set to -1)', 'e')
            determinant_value = COMPLEX(-1,-1) 
        ENDIF

    END FUNCTION determinant
    
    FUNCTION eigenvalues(matrix_in) RESULT(vector_output)

        ! FUNCTION eigenvalues 
        !  -> This function computes the eigenvalues of a mymatrix TYPE
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - matrix_in:     input mymatrix
        ! - vector_output: vector of eigenvalues in output
        ! - work(:):
        ! - rwork(:):
        ! - info:
        ! - lwork(:):
        ! - N:             dimension of the rectangular matrix
        !
        ! ERROR HANDLING:
        ! This function does not compute the eigenvalues if matrix_in is rectangular 

        TYPE(mymatrix), INTENT(IN)                   :: matrix_in
        !eigenvalues vector in double precision
        REAL(8),    DIMENSION(:), ALLOCATABLE        :: vector_output
        !auxiliary variables for ZHEEV subroutine
        COMPLEX*16, DIMENSION(:), ALLOCATABLE        :: work(:)
        REAL(8),    DIMENSION(:), ALLOCATABLE        :: rwork(:)
        INTEGER(4)                                   :: info, lwork, N
        
        IF (matrix_in%m_dims(1).EQ.matrix_in%m_dims(2)) THEN
            N = matrix_in%m_dims(1)
            ALLOCATE(vector_output(N))
    
            !allocate auxiliary variables for ZHEEV subroutine
            lwork = MAX(1,2*N-1)
            ALLOCATE(work(MAX(1, lwork)))
            ALLOCATE(rwork(MAX(1, 3*N-2)))
    
            CALL ZHEEV('N','U', N, matrix_in%m_el, N, vector_output, work, lwork, rwork, info)
            
            CALL dlasrt('I', N, vector_output, INFO)

            DEALLOCATE(work)
            DEALLOCATE(rwork)
        ELSE
            CALL checkpoint(.TRUE., 'Rectangular matrix: cannot calculate the eigenvalues', 'e')
            STOP
        ENDIF

    END FUNCTION eigenvalues

    FUNCTION norm_eigenvalues_spacing(eigs) RESULT(vector_output)

        !eigenvalues vector in double precision
        REAL(8),    DIMENSION(:), INTENT(IN)         :: eigs
        REAL(8),    DIMENSION(:), ALLOCATABLE        :: vector_output
        REAL(8)                                      :: mean
        INTEGER(4)                                   :: ee
    
        ALLOCATE(vector_output(SIZE(eigs)-1))

        vector_output = eigs(2:SIZE(eigs)) - eigs(1:(SIZE(eigs)-1))
        mean = SUM(vector_output)/SIZE(vector_output)

        do ee = 1, SIZE(vector_output)
            vector_output(ee) = vector_output(ee)/mean
        END DO
        
    END FUNCTION norm_eigenvalues_spacing

    SUBROUTINE print_to_file(matrix_in, unit, filename)

        ! SUBROUTINE print_to_file
        !  -> This function prints to  a file, in a formatted fashion,
        !      a mymatrix TYPE
        !
        ! =================================================================
        !
        ! PARAMETERS:
        ! - matrix_in:           input mymatrix
        ! - unit:                specify the unit number to use in writing to a file
        ! - ii, jj:              integer iterators
        ! - filename (optional): name of file in output + extension
        !                        if not specified is set to 'matrix.txt'
        !

        TYPE(mymatrix), INTENT(IN)         :: matrix_in
        INTEGER(4)                         :: unit, ii, jj
        CHARACTER(*), INTENT(IN), OPTIONAL :: filename

        IF(PRESENT(filename)) THEN
            OPEN(UNIT = unit, FILE = filename, ACTION = 'write', STATUS = 'replace')
        ELSE
            OPEN(UNIT = unit, FILE = 'matrix.txt', ACTION = 'write', STATUS = 'replace')
        ENDIF
        
        WRITE(unit, *) "Matrix rows:    ", matrix_in%m_dims(1)
        WRITE(unit, *) "Matrix columns: ", matrix_in%m_dims(2)

        WRITE(unit, *) "Matrix elements row-by-row:"
        DO ii = 1,matrix_in%m_dims(1)
            write (unit, "(*('('sf10.7xspf10.7x'i)':x))") (matrix_in%m_el(ii,jj),jj=1,matrix_in%m_dims(2))
        END DO

        WRITE(unit, *) "Matrix trace:       ", matrix_in%m_trace
        WRITE(unit, *) "Matrix determinant: ", matrix_in%m_det

        CLOSE(unit)

    END SUBROUTINE print_to_file
    
END MODULE dc_matrix