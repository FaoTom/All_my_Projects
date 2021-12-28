MODULE QM_utils

    ! MODULE QM_utils
    !  -> This module contains some useful functions and subroutine for solving computational
    !     problems involving quantum mechanics.
    !
    ! =================================================================
    !
    ! FUNCTIONS AND SUBROUTINES:
    ! - FUNCTION   harmosc_H_op_init(N, Lmin, Lmax, m, omega) result(ham_op)
    ! - SUBROUTINE diagonalize_H_op(ham_op, eigs)

USE error_handling

IMPLICIT NONE

CONTAINS

FUNCTION harmosc_H_op_init(N, Lmin, Lmax, m, omega) result(ham_op)

    ! FUNCTION harmosc_H_op_init
    !  -> This function returns a discretized Hamiltonian operator
    !     for the quantum harmonic oscillator in space representation
    !     ψ(x)=<ψ|x>, up to O(dx^4) precision, where dx is the discretization step.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - N:      number of discretization step
    ! - Lmin:   left bound of discretized interval 
    ! - Lmax:   right bound of discretized interval
    ! - m:      mass of harmonic oscillator
    ! - omega:  harmonic oscillator angular frequency
    ! - ham_op: discretized Hamiltonian operator in output
    !

    !variables for the one dimensional lattice
    REAL(8)                                      :: Lmin, Lmax, dx
    !coefficients that enters the harm osc hamiltonian
    REAL(8)                                      :: m, omega, h_bar
    INTEGER(4)                                   :: N, aa, bb
    !harmonic oscillator hamiltonian operator
    REAL(8),    DIMENSION(:,:), ALLOCATABLE      :: ham_op

    !set reduced Planck constant to 1
    h_bar = 1

    !allocate and initialize the operator
    ALLOCATE(ham_op(N+1,N+1))

    !initialize all to zero and then fill it with for loops
    ham_op=0

    !compute the spacing of the discretized one dimensional lattice
    dx = (Lmax-Lmin)/N

    !filling the matrix
    DO aa = 1 ,N+1
        DO bb = 1,N+1
            IF(aa.EQ.bb) THEN
                ham_op(aa,bb) = (h_bar**2)/(dx**2)
            ELSE IF(aa.EQ.bb+1) THEN
                ham_op(aa,bb) = -(h_bar**2)/(2*dx**2)
                ham_op(bb,aa) = ham_op(aa,bb)
            ENDIF
        ENDDO
    ENDDO

    !adding the lattice-position dependent value to the diagonal elements
    DO aa = 1, N+1
        ham_op(aa,aa) = ham_op(aa,aa) + 0.5*m*omega**2*((Lmin+(aa-1)*dx)**2)
    ENDDO

END FUNCTION harmosc_H_op_init

SUBROUTINE diagonalize_H_op(ham_op, eigs)

    ! SUBROUTINE diagonalize_H_op
    !  -> This subroutine diagonalizes a tri-diagonal symmetric
    !      matrix with real entries and returns eigenvectors, and
    !      eigenvalues. Here we use DSTEV subroutine.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - ham_op: discretized Hamiltonian operator in input
    ! - eigs:   vector with eigenvalues in output
    !
    
    !hamiltonian operator in input
    REAL(8),    DIMENSION(:,:), ALLOCATABLE      :: ham_op
    REAL(8),    DIMENSION(:),   ALLOCATABLE      :: eigs, diag_A, subd_A
    !auxiliary variables for DSTEV subroutine
    REAL(8),    DIMENSION(:),   ALLOCATABLE      :: work(:)
    INTEGER(4)                                   :: info, lwork, N, aa, bb
    REAL(8)                                      :: norm

    !PRE-CONDITION: check if diagonalization is possible
    IF (size(ham_op, 1).NE.size(ham_op, 2)) THEN
        CALL checkpoint(.TRUE.,'Rectangular matrix: cannot diagonalize', 'e')
        STOP
    ENDIF

    N = SIZE(ham_op, 1)
    lwork = 2*(N+1)-2
    
    info = 0 !in this way, on exit of DSTEV, diag_A will contain the eigenvalues

    ALLOCATE(diag_A(N))
    ALLOCATE(subd_A(N-1))
    ALLOCATE(work(lwork))

    !filling the array of elements on the diagonal
    DO aa = 1, N
        diag_A(aa) = ham_op(aa,aa)
    ENDDO

    !filling the array of elements on the subdiagonal
    DO aa = 1, N-1
        subd_A(aa) = ham_op(aa,aa+1)
    ENDDO

    !compute eigenvalues and eigenvectors
    CALL DSTEV('V', N, diag_A, subd_A, ham_op, N, work, info)

    IF(info.NE.0) THEN
        CALL checkpoint(.TRUE., 'DSTEV exited with some error', 'e')
        STOP
    ENDIF

    !if info = 0, on exit, diag_A contains the eigenvalues in ascending order
    eigs = diag_A

    !normalize k-th eigenfunciton ψk dividing by SQRT(|<ψk|ψk>|^2)
    DO bb = 1,SIZE(ham_op,2)
        norm = DOT_PRODUCT(ham_op(:,bb), ham_op(:,bb))
        ham_op(:,bb) = ham_op(:,bb) / SQRT(norm)
    ENDDO

    DEALLOCATE(work)
    DEALLOCATE(diag_A)
    DEALLOCATE(subd_A)

END SUBROUTINE diagonalize_H_op

END MODULE QM_utils