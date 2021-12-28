MODULE my_mat_mul

    ! MODULE my_mat_mul
    !  -> This module my_mat_mul can be used to perform matrix 
    !     matrix matrix multiplication in two different ways:
    !     matmul_byrow and matmul_bycol. The last two methods will
    !     show different results in terms of CPU time due to
    !     memory caching.
    !
    ! =================================================================
    !
    ! FUNCTIONS AND SUBROUTINES:
    ! - FUNCTION  matmul_byrow(mat_A, mat_B)  RESULT(mat_C1)
    ! - FUNCTION  matmul_bycol(mat_A, mat_B)  RESULT(mat_C1)
    !

IMPLICIT NONE

!define iterators
INTEGER(4) :: ii, jj, kk

CONTAINS

FUNCTION matmul_byrow(mat_A, mat_B) RESULT(mat_C)

    ! FUNCTION matmul_byrow
    !  -> This function multiplies two matrices of single precision
    !     real variables as entries by rows.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - mat_A: first input matrix
    ! - mat_B: second input matrix
    ! - mat_C: output matrix
    !

    REAL(4), DIMENSION(:,:)                         :: mat_A, mat_B
    REAL(4), DIMENSION(SIZE(mat_A,1),SIZE(mat_B,2)) :: mat_C

    DO ii = 1, SIZE(mat_A,1)
        DO kk = 1, SIZE(mat_B,2)
            DO jj = 1, SIZE(mat_A,2)
                mat_C(ii,jj) = mat_C(ii,jj) + mat_A(ii,kk) * mat_B(kk,jj)
            ENDDO
        ENDDO
    ENDDO 
    
END FUNCTION matmul_byrow

FUNCTION matmul_bycol(mat_A, mat_B) RESULT(mat_C)

    ! FUNCTION matmul_bycol
    !  -> This function multiplies two matrices of single precision
    !     real variables as entries by columns.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - mat_A: first input matrix
    ! - mat_B: second input matrix
    ! - mat_C: output matrix
    !

    REAL(4), DIMENSION(:,:)                         :: mat_A, mat_B
    REAL(4), DIMENSION(SIZE(mat_A,1),SIZE(mat_B,2)) :: mat_C

    DO jj = 1, SIZE(mat_A,2)
        DO kk = 1, SIZE(mat_B,2)
            DO ii = 1, SIZE(mat_A,1)
                mat_C(ii,jj) = mat_C(ii,jj) + mat_A(ii,kk) * mat_B(kk,jj)
            ENDDO
        ENDDO
    ENDDO 
    
END FUNCTION matmul_bycol

END MODULE my_mat_mul
