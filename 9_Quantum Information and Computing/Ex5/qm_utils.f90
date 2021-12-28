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
    ! - FUNCTION H_kinetic_part(N, Lmin, Lmax, m) RESULT(K)
    ! - FUNCTION H_potential_part(N, Lmin, Lmax, m, omega, time, T) RESULT(V)
    ! - FUNCTION myFFTW(input) RESULT(output)
    ! - FUNCTION myAFFTW(input) RESULT(output)
    ! - SUBROUTINE evolve_dt(wf, dt, Lmin, Lmax, m, omega, time, T)

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
                ham_op(aa,bb) = (h_bar**2)/(m*dx**2)
            ELSE IF(aa.EQ.bb+1) THEN
                ham_op(aa,bb) = -(h_bar**2)/(2*m*dx**2)
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
    
    REAL(8),    DIMENSION(:,:), ALLOCATABLE      :: ham_op
    REAL(8),    DIMENSION(:),   ALLOCATABLE      :: eigs, diag_A, subd_A
    !auxiliary variables for DSTEV subroutine
    REAL(8),    DIMENSION(:),   ALLOCATABLE      :: work(:)
    INTEGER(4)                                   :: info, lwork, N, aa, bb
    REAL(8)                                      :: norm

    !PRE-CONDITION: check if diagonalization is possible
    IF (SIZE(ham_op, 1).NE.SIZE(ham_op, 2)) THEN
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

    !compute eigenvalues and eigenvectors with the suitable routine
    CALL DSTEV('V', N, diag_A, subd_A, ham_op, N, work, info)

    !if info = 0, on exit, diag_A contains the eigenvalues in ascending order
    eigs = diag_A

    DO bb=1,SIZE(ham_op,2)
        norm = dot_product(ham_op(:,bb),ham_op(:,bb))
        ham_op(:,bb) = ham_op(:,bb)/norm
    ENDDO

    DEALLOCATE(work)
    DEALLOCATE(diag_A)
    DEALLOCATE(subd_A)

END SUBROUTINE diagonalize_H_op

FUNCTION H_kinetic_part(N, Lmin, Lmax, m) RESULT(K)

    ! FUNCTION H_kinetic_part
    !  -> This function returns the discretized kinetic part for the
    !     Hamiltonian operator of the quantum harmonic oscillator in momentum
    !     representation up to O(dx^4) precision, where dx is the discretization step. 
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - Lmin:   left bound of discretized interval 
    ! - Lmax:   right bound of discretized interval
    ! - N:      number of discretization step
    ! - aa:     iteratore
    ! - m:      mass of harmonic oscillator
    ! - K:      discretized kinetic term of harmonic oscillator
    !
    
    REAL(8)                                      :: Lmin, Lmax
    INTEGER(4)                                   :: N, aa
    REAL(8)                                      :: m
    REAL(8),   DIMENSION(N+1)                    :: K

    DO aa=1,INT(N/2)
        K(aa) = (0.5d0/m) * (4.0d0*ASIN(1.0d0)*(aa)/(Lmax-Lmin))**2
    END DO

    DO aa=INT(N/2)+1,N+1
        K(aa) = (0.5d0/m) * (4.0d0*ASIN(1.0d0)*(aa-N-1)/(Lmax-Lmin))**2
    ENDDO

END FUNCTION H_kinetic_part

FUNCTION H_potential_part(N, Lmin, Lmax, m, omega, time, T) RESULT(V)

    ! FUNCTION H_potential_part
    !  -> This function returns the discretized potential part for the
    !     Hamiltonian operator of the quantum harmonic oscillator in 
    !     space representation ψ(x)=<ψ|x>, up to O(dx^4) precision, 
    !     where dx is the discretization step and for a simulation time T. 
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - Lmin:   left bound of discretized interval 
    ! - Lmax:   right bound of discretized interval
    ! - dx:     discretization step
    ! - T:      characteristic time of the potential
    ! - time:   current time instant of the simulation
    ! - N:      number of discretization step
    ! - aa:     iterator
    ! - m:      mass of harmonic oscillator
    ! - omega:  harmonic oscillator angular frequency
    ! - V:      discretized potential term of harmonic oscillator
    !

    REAL(8)                                      :: Lmin, Lmax, dx, T, time
    INTEGER(4)                                   :: N, aa
    REAL(8)                                      :: m, omega
    REAL(8),    DIMENSION(N+1)                   :: V

    dx = (Lmax-Lmin)/N
    DO aa = 1, N+1
        V(aa) = 0.5*m*omega**2*(Lmin+(aa-1)*dx-time/T)**2
    ENDDO

END FUNCTION H_potential_part

FUNCTION myFFTW(input) RESULT(output)

    ! FUNCTION myFFTW
    !  -> This function myFFTW calls some wrappers routines, written in C,
    !     to the original FFTW,  in order to compute the one dimensional 
    !     discrete Fourier transform.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - input:           wavefunction in the coordinate space
    ! - output:          wavefunction in the momentum space
    ! - temp:            copy of the input
    ! - plan:            best algorithm to perform the FFTW
    ! - N:               size of the input wavefunction
    ! - FFTW_FORWARD:    specify to go forward (sign -1 in the exponential of Fourier transform)
    ! - FFTW_MEASURE:    specify the dimension of the set of algorithm in which to look for the best 
    !

    DOUBLE COMPLEX, DIMENSION(:)                 :: input
    DOUBLE COMPLEX, DIMENSION(SIZE(input, 1))    :: output, temp
    !at least as big as a pointer (address) on your machine (integer*8 recommended)
    INTEGER(8)                                   :: plan
    INTEGER(4)                                   :: N
    INTEGER                                      :: FFTW_FORWARD
    PARAMETER(FFTW_FORWARD=-1)
    INTEGER                                      :: FFTW_MEASURE 
    PARAMETER(FFTW_MEASURE=64) 

    N = SIZE(input, 1)
    
    !initialize the array once the plan is created
    temp = input

    CALL dfftw_plan_dft_1d( plan, N, input, output, FFTW_FORWARD, FFTW_MEASURE)
    CALL dfftw_execute_dft( plan, temp, output)
    output = output/SQRT(1.0d0*N)
    CALL dfftw_destroy_plan(plan)

END FUNCTION myFFTW

FUNCTION myAFFTW(input) RESULT(output)

    ! FUNCTION myAFFTW
    !  -> This function myAFFTW calls some wrappers routines written in C
    !     in order to compute the one dimensional discrete Fourier anti-transform.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - input:           wavefunction in the momentum space
    ! - output:          wavefunction in the coordinate space
    ! - temp:            copy of the input
    ! - plan:            best algorithm to perform the FFTW
    ! - N:               size of the input wavefunction
    ! - FFTW_BACKWARD:   specify to go forward (sign -1 in the exponential of Fourier transform)
    ! - FFTW_MEASURE:    specify the dimension of the set of algorithm in which to look for the best 
    !
    
    COMPLEX(8), DIMENSION(:)                     :: input
    COMPLEX(8), DIMENSION(SIZE(input, 1))        :: output, temp
    INTEGER(8)                                   :: plan
    INTEGER(4)                                   :: N
    INTEGER                                      :: FFTW_BACKWARD
    PARAMETER(FFTW_BACKWARD=1)
    INTEGER                                      :: FFTW_MEASURE       
    PARAMETER(FFTW_MEASURE=64)

    N = SIZE(input, 1)

    !initialize the array once the plan is created
    temp = input

    CALL dfftw_plan_dft_1d( plan, N, input, output, FFTW_BACKWARD, FFTW_MEASURE)
    CALL dfftw_execute_dft( plan, temp,  output)
    output = output/SQRT(1.0d0*N)
    CALL dfftw_destroy_plan(plan)

END FUNCTION myAFFTW

SUBROUTINE evolve_dt(wf, dt, Lmin, Lmax, m, omega, time, T)

    ! SUBROUTINE evolve_dt
    !  -> This subroutine evolve_dt applies the splitting operator
    !     method, used to evolve the function for a time step dt.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - wf:              wavefunction in the momentum space to be evolved
    ! - Lmin:            left bound of discretized interval 
    ! - Lmax:            right bound of discretized interval
    ! - dt:              time discretization step
    ! - T:               characteristic time of the potential
    ! - N:               number of discretization step
    ! - aa:              iterator
    ! - m:               mass of harmonic oscillator
    ! - omega:           harmonic oscillator angular frequency
    ! - V:               discretized potential term of harmonic oscillator
    ! - K:               discretized kinetic term of harmonic oscillator
    !

    DOUBLE COMPLEX, DIMENSION(:)                 :: wf
    REAL(8)                                      :: Lmin, Lmax, dt, T, time
    INTEGER(4)                                   :: N, aa
    REAL(8)                                      :: m, omega
    REAL(8),    DIMENSION(:),  ALLOCATABLE       :: V, K

    !recall that the wavefunction is discretize over N+1 points
    N = size(wf, 1)-1

    ALLOCATE(V(N+1))
    ALLOCATE(K(N+1))

    V = H_potential_part(N, Lmin, Lmax, m, omega, time, T)
    K = H_kinetic_part(  N, Lmin, Lmax, m)
    
    DO aa = 1, N+1
        wf(aa) = EXP(COMPLEX(0.0d0, - 1.0d0 / 2.0d0 * dt * V(aa))) * wf(aa)
    END DO
    wf = myFFTW(wf)
    DO aa = 1, N+1
        wf(aa) = EXP(COMPLEX(0.0d0, -1.0d0 * dt * K(aa))) * wf(aa)
    END DO
    wf = myAFFTW(wf)
    DO aa = 1, N+1
        wf(aa) = EXP(COMPLEX(0.0d0, -0.5d0 * dt * V(aa))) * wf(aa)
    END DO
    
END SUBROUTINE evolve_dt

END MODULE QM_utils