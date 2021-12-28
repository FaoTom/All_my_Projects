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
    ! - FUNCTION sep_wf_init(local_dimension, system_dimension) RESULT(output)
    ! - FUNCTION gen_wf_init(local_dimension, system_dimension) RESULT(output)
    ! - FUNCTION p_trace(local_dimension, dmatrix, i_subsystem) RESULT(red_dmatrix)

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
    ! - N:                number of discretization step
    ! - Lmin:             left bound of discretized interval 
    ! - Lmax:             right bound of discretized interval
    ! - m:                mass of harmonic oscillator
    ! - omega:            harmonic oscillator angular frequency
    ! - ham_op:           discretized Hamiltonian operator in output
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
    !  -> This subroutine 
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - N:                number of discretization step
    ! - ham_op:           discretized Hamiltonian operator in input
    ! - eigs:             vector with eigenvalues in output
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
    call DSTEV('V', N, diag_A, subd_A, ham_op, N, work, info)

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
    ! - N:                number of discretization step
    ! - Lmin:             left bound of discretized interval 
    ! - Lmax:             right bound of discretized interval
    ! - m:                mass of harmonic oscillator
    ! - K:                discretized kinetic term of harmonic oscillator
    !
    
    REAL(8)                                      :: Lmin, Lmax
    INTEGER(4)                                   :: N, aa
    REAL(8)                                      :: m
    REAL(8),   DIMENSION(N+1)                    :: K

    DO aa=1,INT(N/2)
        K(aa) = (0.5d0/m) * (4.0d0*ASIN(1.0d0)*(aa)/(Lmax-Lmin))**2
    END DO

    do aa=INT(N/2)+1,N+1
        K(aa) = (0.5d0/m) * (4.0d0*ASIN(1.0d0)*(aa-N-1)/(Lmax-Lmin))**2
    end do

END FUNCTION H_kinetic_part

FUNCTION H_potential_part(N, Lmin, Lmax, m, omega, time, T) RESULT(V)

    ! FUNCTION H_potential_part
    !  -> This function returns the discretized potential part for the
    !     space representation ψ(x)=<ψ|x>, up to O(dx^4) precision, 
    !     Hamiltonian operator of the quantum harmonic oscillator in 
    !     where dx is the discretization step and for a simulation time T. 
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - N:                number of discretization step
    ! - Lmin:             left bound of discretized interval 
    ! - Lmax:             right bound of discretized interval
    ! - m:                mass of harmonic oscillator
    ! - omega:            harmonic oscillator angular frequency
    ! - time:             initial instant of the simulation
    ! - T:                characteristic time of the simulation
    ! - V:                discretized kinetic term of harmonic oscillator
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
    !  -> This function myFFTW calls some wrappers routines written in C
    !     in order to compute the one dimensional discrete Fourier transform.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - wf:               wavefunction in the coordinate space
    ! - fftw_wf:          wavefunction in the momentum space
    ! - plan:             container for data useful to perform the FFT
    ! - FFTW_FORWARD:     specify to go forward (sign -1 in the exponential of Fourier transform)
    ! - FFTW_ESTIMATE:   
    !

    DOUBLE COMPLEX, DIMENSION(:)                 :: input
    DOUBLE COMPLEX, DIMENSION(SIZE(input, 1))    :: output, temp
    !at least as big as a pointer (address) on your machine (integer*8 recommended)
    INTEGER(8)                                   :: plan
    INTEGER(4)                                   :: N
    INTEGER                                  :: FFTW_FORWARD
    PARAMETER(FFTW_FORWARD=-1)
    INTEGER                                   :: FFTW_MEASURE 
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
    ! - fftw_wf:          wavefunction in the momentum space
    ! - wf:               wavefunction in the coordinate space
    ! - plan:             container for data useful to perform the FFT
    ! - FFTW_BACKWARD:    specify to go backward (sign +1 in the exponential of Fourier transform)
    ! - FFTW_ESTIMATE:   
    !
    
    DOUBLE COMPLEX, DIMENSION(:)                 :: input
    DOUBLE COMPLEX, DIMENSION(SIZE(input, 1))    :: output, temp
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
    !  -> s
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - fftw_wf:          wavefunction in the momentum space
    ! - wf:               wavefunction in the coordinate space
    ! - plan:             container for data useful to perform the FFT
    ! - FFTW_BACKWARD:    specify to go backward (sign +1 in the exponential of Fourier transform)
    ! - FFTW_ESTIMATE:   
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

! sostituire qui sopra non appena ho commentato per bene la settimana precedente.

FUNCTION sep_wf_init(local_dimension, system_dimension) RESULT(output)

    ! FUNCTION sep_wf_init
    !  -> This function initializes a separable wave function of a 
    !     system_dimension-body quantum system where each of them is 
    !     embedded in a space of dimension local_dimension 
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - output:           output wavefunction
    ! - re:               real part of wavefunction
    ! - im:               imaginary part of wavefunction
    ! - norm:             norm of the wavefunction
    ! - local_dimension:  local dimension
    ! - system_dimension: system dimension
    !

    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE   :: output
    REAL(8),        DIMENSION(:), ALLOCATABLE   :: re, im
    DOUBLE COMPLEX                              :: norm
    INTEGER(4)                                  :: local_dimension, system_dimension

    ALLOCATE(re(local_dimension*system_dimension))
    ALLOCATE(im(local_dimension*system_dimension))
    ALLOCATE(output(local_dimension*system_dimension))

    call RANDOM_NUMBER(re)
    call RANDOM_NUMBER(im)

    output = ((1.0d0-2.0d0)*re*COMPLEX(1.0d0,0.0d0) + (1.0d0-2.0d0)*im*COMPLEX(0.0d0,1.0d0))
    norm = dot_product(output,output)
    output=output/SQRT(norm)

END FUNCTION sep_wf_init

FUNCTION gen_wf_init(local_dimension, system_dimension) RESULT(output)

    ! FUNCTION gen_wf_init
    !  -> This function initializes a pure wave function of a 
    !     system_dimension-body quantum system where each of them is 
    !     embedded in a space of dimension local_dimension 
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - output:           output wavefunction
    ! - re:               real part of wavefunction
    ! - im:               imaginary part of wavefunction
    ! - norm:             norm of the wavefunction
    ! - local_dimension:  local dimension
    ! - system_dimension: system dimension
    !

    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE   :: output
    REAL(8),        DIMENSION(:), ALLOCATABLE   :: re, im
    DOUBLE COMPLEX                              :: norm
    INTEGER(4)                                  :: local_dimension, system_dimension

    ALLOCATE(re(local_dimension**system_dimension))
    ALLOCATE(im(local_dimension**system_dimension))
    ALLOCATE(output(local_dimension**system_dimension))

    call RANDOM_NUMBER(re)
    call RANDOM_NUMBER(im)

    output = ((1.0d0-2.0d0)*re*COMPLEX(1.0d0,0.0d0) + (1.0d0-2.0d0)*im*COMPLEX(0.0d0,1.0d0))
    norm = dot_product(output,output)
    output=output/SQRT(norm)

END FUNCTION gen_wf_init

FUNCTION dmatrix(state) RESULT(output)

    ! FUNCTION dmatrix
    !  -> This function computes the density matrix of a generic quantum state
    !     The dimension is 2^Nx2^N given a 2^N dimensional state vector.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - output:           output matrix
    ! - state:            real part of wavefunction
    ! - b:                BRA
    ! - k:                KET
    ! - N:                size of the vector state
    !
    
    DOUBLE COMPLEX, DIMENSION(:),   ALLOCATABLE   :: state
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: output
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: b, k
    INTEGER(4)                                    :: N

    N=size(state,1)

    ALLOCATE(output(N,N), k(N,1), b(1,N))

    k(:,1) = state
    b(1,:) = CONJG(state)

    output = MATMUL(k, b)

END FUNCTION dmatrix

FUNCTION p_trace(local_dimension, dmatrix, i_subsystem) RESULT(red_dmatrix)

    ! FUNCTION p_trace
    !  -> This function computes the partial trace of a given density
    !     matrix with respect to the i_subsystem-th subsystem.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - dmatrix:          input density matrix
    ! - red_dmatrix:      output reduced density matrix
    ! - local_dimension:  local dimension of each subsystem
    ! - dN:               size of input density matrix
    ! - d:                local dimension of each subsystem
    ! - aa, bb, cc:       iterators
    ! - i_subsystem:      subsystem to trace out
    ! - row, col          useful vector used for the computation
    !

    DOUBLE COMPLEX, DIMENSION(:,:),   ALLOCATABLE   :: dmatrix, red_dmatrix
    INTEGER(4)                                      :: local_dimension, dN, d
    INTEGER(4)                                      :: aa, bb, cc, i_subsystem
    INTEGER(4), DIMENSION(:), ALLOCATABLE           :: row, col

    !in this case dN = d^N
    dN = SIZE(dmatrix,1)
    d = local_dimension

    !given a local dimension d and the size of the complete system N, 
    !the reduce density matrix has dimension (d^(N-1),d^(N-1))
    ALLOCATE(red_dmatrix(dN/d, dN/d))
    ALLOCATE(row(d), col(d))

    !intialize the reduced dansity matrix
    DO aa = 1, dN/d
        DO bb = 1, dN/d
            red_dmatrix(aa, bb) = COMPLEX(0.0d0, 0.0d0)
        ENDDO
    ENDDO

    IF (i_subsystem == 2) THEN
        DO aa=1, d
            row(aa) = 2*aa-1
            col(aa) = 2*aa-1
        ENDDO

        DO aa = 1, SIZE(row)
            DO bb = 1, SIZE(col)
                DO cc = 1, d
                    red_dmatrix(aa, bb) = red_dmatrix(aa, bb) + dmatrix(row(aa)+(cc-1), row(bb)+(cc-1))
                ENDDO 
            ENDDO
        ENDDO

    ELSE IF(i_subsystem == 1) THEN
        DO aa = 1, dN/d
            DO bb = 1, dN/d
                DO cc = 1, d
                    red_dmatrix(aa, bb) = red_dmatrix(aa, bb) + dmatrix(aa+(cc-1)*d, bb+(cc-1)*d)
                ENDDO
            ENDDO 
        ENDDO

    ELSE
        call checkpoint(.TRUE., 'The function is not built for more than two subsystems', 'e')
    ENDIF

END FUNCTION p_trace

FUNCTION tp(m1 ,m2) RESULT(output)

    ! FUNCTION dmatrix
    !  -> This function computes the density matrix of a generic quantum state
    !     The dimension is 2^Nx2^N given a 2^N dimensional state vector.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - output:           output matrix
    ! - state:            real part of wavefunction
    ! - b:                BRA
    ! - k:                KET
    ! - N:                size of the vector state
    !
    
    DOUBLE COMPLEX, DIMENSION(:,:)                :: m1, m2
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: output
    INTEGER(4)                                    :: ii, jj
    INTEGER(4)                                    :: n_col1, n_col2
    INTEGER(4)                                    :: n_row1, n_row2
    INTEGER(4)                                    :: rf, cf, rl, cl

    n_row1 = SIZE(m1, 1)
    n_col1 = SIZE(m1, 2)
    n_row2 = SIZE(m2, 1)
    n_col2 = SIZE(m2, 2)

    ALLOCATE(output(n_row1*n_row2,n_col1*n_col2))

    !cycle on the elements of the first matrix
    DO ii = 1, n_row1
        DO jj = 1, n_col1
            !building the indices of the output matrix
            rf = (ii - 1) * n_row2 + 1
            cf = (jj - 1) * n_col2 + 1
            rl = ii * n_row2
            cl = jj * n_col2
            output(rf:rl,cf:cl) = m1(ii,jj)*m2
        ENDDO
    ENDDO

END FUNCTION tp

FUNCTION sigma_z() RESULT(output)

    ! FUNCTION sigma_z
    !  -> This function returns the sigma_z Pauli matrix
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - output:           output matrix
    !

    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: output

    ALLOCATE(output(2,2))

    output      = COMPLEX( 0.0d0, 0.0d0)
    output(1,1) = COMPLEX( 1.0d0, 0.0d0)
    output(2,2) = COMPLEX(-1.0d0, 0.0d0)

END FUNCTION sigma_z

FUNCTION sigma_x() RESULT(output)

    ! FUNCTION sigma_x
    !  -> This function returns the sigma_x Pauli matrix
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - output:           output matrix
    !

    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: output

    ALLOCATE(output(2,2))

    output      = COMPLEX( 0.0d0, 0.0d0)
    output(1,2) = COMPLEX( 1.0d0, 0.0d0)
    output(2,1) = COMPLEX( 1.0d0, 0.0d0)

END FUNCTION sigma_x

FUNCTION id_mat_init(N) RESULT(output)

    ! FUNCTION id_mat_init
    !  -> This function initializes an identity matrix for N particles
    !     so of dimension (2**N, 2**N).
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - N:                number of particle
    ! - ii:               iterator
    ! - output:           identity matrix in output
    !

    INTEGER(4)                            :: N, ii
    DOUBLE COMPLEX, dimension(2**N, 2**N) :: output

    output = COMPLEX(0.0d0, 0.0d0)

    DO ii = 1, 2**N
        output(ii,ii) = COMPLEX(1.0d0, 0.0d0)
    ENDDO

END FUNCTION id_mat_init

FUNCTION TFI_ham_init(N, lambda) RESULT(output)

    ! FUNCTION TFI_ham_init
    !  -> This function initializes the hamiltonian of the TFI.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - output:           TFI hamiltonian in output
    ! - N:                system size
    ! - ii:               iterator
    ! - lambda:           external field
    !
    
    DOUBLE COMPLEX, DIMENSION(2**N,2**N) :: output
    INTEGER(4)                           :: N, ii
    REAL(8)                              :: lambda

    output = COMPLEX(0.0d0, 0.0d0)

    DO ii = 1, N
        !single site component
        output = output + tp(tp(id_mat_init(ii-1), sigma_z()),id_mat_init(N-ii))
    ENDDO

    output = output * lambda* COMPLEX(1.0d0, 0.0d0)

    DO ii = 1, N-1
        output = output - tp(tp(tp(id_mat_init(ii-1), sigma_x()),sigma_x()),id_mat_init(N-ii-1))
    ENDDO

END FUNCTION TFI_ham_init

END MODULE QM_utils