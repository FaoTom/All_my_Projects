!Week 3 - Exercise 2, Faorlin Tommaso (2021857)
!compile with gfortran -W w3_ex2.f90 dc_matrix.f90 error_handling.f90 my_mat_mul.f90 -o w3_ex2.o -llapack 

PROGRAM eigenproblem

    USE my_mat_mul
    USE error_handling
    USE dc_matrix
    
    IMPLICIT NONE

    INTEGER(4)                                   :: N, ee, statistic, ss
    !boolean debug variable
    LOGICAL                                      :: debug
    !dimensions array for mymatrix type
    INTEGER(4), DIMENSION(2)                     :: m_dims
    !test mymatrix
    TYPE(mymatrix)                               :: test_hermitian, test_diagonal
    !eigenvalues vector in double precision
    REAL(8),    DIMENSION(:),   ALLOCATABLE      :: herm_eigenvalues_vec, diag_eigenvalues_vec&
                                                    &,herm_norm_spacing, diag_norm_spacing

    !set .TRUE. to enter debug mode
    debug = .TRUE.
    statistic = 80

    !set fixed dimension of the matrices
    N = 2500
    m_dims(1) = N
    m_dims(2) = N
    
    DO ss = 1, statistic
    !allocate space for the hermitian random matrix
        ALLOCATE(herm_eigenvalues_vec(N))
        !allocate space for the diagonal random matrix
        ALLOCATE(diag_eigenvalues_vec(N))
        !allocate space for eigenvalues spacing for hermitian matrix
        ALLOCATE(herm_norm_spacing(N-1))
        !allocate space for eigenvalues spacing for diagonal matrix
        ALLOCATE(diag_norm_spacing(N-1))
    
        !initialize the matrices
        test_hermitian = hermitian_init(m_dims)
        test_diagonal  = diagonal_init(m_dims)
    
        !compute eigenvalues
        herm_eigenvalues_vec = eigenvalues(test_hermitian)
        diag_eigenvalues_vec = eigenvalues(test_diagonal)
    
        !compute eigenvalues spacing
        herm_norm_spacing = norm_eigenvalues_spacing(herm_eigenvalues_vec)
        diag_norm_spacing = norm_eigenvalues_spacing(diag_eigenvalues_vec)
    
        !write to file hermitian matrix eigenvalues spacing
        OPEN (10, FILE="data/herm_eigenvalues_spacing_norm.dat", POSITION = "append", ACTION = "WRITE")
        DO ee = 1, SIZE(herm_norm_spacing)
            WRITE (10,'(F15.7)') herm_norm_spacing(ee)
        ENDDO
        CLOSE (10)
        
        !write to file diagonal matrix eigenvalues spacing
        OPEN (10, FILE="data/diag_eigenvalues_spacing_norm.dat", POSITION = "append", ACTION = "WRITE")
        DO ee = 1, SIZE(diag_norm_spacing)
            WRITE (10,'(F15.7)') diag_norm_spacing(ee)
        ENDDO
        CLOSE (10)
    
        !deallocate space for the hermitian random matrix
        DEALLOCATE(herm_eigenvalues_vec)
        !deallocate space for the diagonal random matrix
        DEALLOCATE(diag_eigenvalues_vec)
        !deallocate space for eigenvalues spacing for hermitian matrix
        DEALLOCATE(herm_norm_spacing)
        !deallocate space for eigenvalues spacing for diagonal matrix
        DEALLOCATE(diag_norm_spacing)

        CALL checkpoint(debug, 'Iteration ended', 'd')

    ENDDO

    END PROGRAM eigenproblem