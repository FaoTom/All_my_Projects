!Week 5 -  Faorlin Tommaso (2021857)
!compile with 'gfortran w5.f90  -I/usr/local/Cellar/fftw/3.3.10/include -L/usr/local/Cellar/fftw/3.3.10/lib dc_matrix.f90 error_handling.f90 my_mat_mul.f90 qm_utils.f90 -o w5.out -lfftw3 -lfftw3f -llapack -lm  -ffree-line-length-0'

!run with './w5.out -custom 2000 -10 10 3000 0 40 1 0'
!N_x Lmin Lmax N_t Tmin Tmax mass omega

PROGRAM time_dep_schrod

    USE my_mat_mul
    USE error_handling
    USE qm_utils
    
    IMPLICIT NONE

    !discretization of space variables
    REAL(8)                                             :: Lmin, Lmax, dx, expval_x, expval1sigma_x
    INTEGER(4)                                          :: N_x
    REAL(8),        DIMENSION(:), ALLOCATABLE           :: x_latt
    !discretization of time variables       
    REAL(8)                                             :: Tmin, Tmax, dt, T
    INTEGER(4)                                          :: N_t
    REAL(8),        DIMENSION(:), ALLOCATABLE           :: t_latt
    !discretization of momentum variables
    REAL(8),        DIMENSION(:), ALLOCATABLE           :: p_latt
    REAL(8)                                             :: expval_p, expval1sigma_p, PI
    REAL(8),        DIMENSION(:),   ALLOCATABLE         :: expval_p_vector, expval1sigma_p_vector
    !variables for the harmonic oscillator       
    REAL(8)                                             :: m, omega
    !eigenvalues
    REAL(8),        DIMENSION(:),   ALLOCATABLE         :: eigs, expval_x_vector, expval1sigma_x_vector
    DOUBLE COMPLEX, DIMENSION(:),   ALLOCATABLE         :: wf
    !hamiltonian operator       
    REAL(8),        DIMENSION(:,:), ALLOCATABLE         :: ham_op
    INTEGER(4)                                          :: aa, num_commandl_args, iarg, bb
    CHARACTER(10),  DIMENSION(:),    ALLOCATABLE        :: commandl_args
    LOGICAL                                             :: custom
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE         :: total_evolution_x, total_evolution_p

    !header prints
    PRINT *, "================================================================="
    PRINT *, "-> Numerical solution of the time depedent Schrödinger equation <-"
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
    IF(num_commandl_args.EQ.9) THEN
        !if argument is -custom insert manually discretization parameters
        IF(commandl_args(1).EQ.'-custom') THEN
            custom = .TRUE.
            READ(commandl_args(2),"(I8)") iarg
            N_x = iarg
            READ(commandl_args(3),"(I8)") iarg
            Lmin = iarg
            READ(commandl_args(4),"(I8)") iarg
            Lmax = iarg
            READ(commandl_args(5),"(I8)") iarg
            N_t = iarg
            READ(commandl_args(6),"(I8)") iarg
            Tmin = iarg
            READ(commandl_args(7),"(I8)") iarg
            Tmax = iarg
            READ(commandl_args(8),"(I8)") iarg
            m = iarg
            READ(commandl_args(9),"(I8)") iarg
            omega = iarg
        ELSE
            !other arguments for the moment are not accepted
            CALL checkpoint(.TRUE., 'Invalid argument inserted', 'e')
            STOP
        ENDIF
    !no arguments are passed, so no debug mode
    ELSE IF (num_commandl_args.EQ.0) THEN
        custom = .FALSE.
        N_x = 1000
        Lmin = -10
        Lmax = +10
        N_t = 200
        Tmin = 0
        Tmax = 16
        m = 1
        omega = 1
    ELSE
        !having more than one argument is not accepted at the moment
        CALL checkpoint(.TRUE., 'Invalid number of arguments in the command line', 'e')
        STOP 
    ENDIF

    !potential intrinsic time
    T = 1

    IF(Tmax.LE.0) THEN
        !Tmax should be greater than zero
        CALL checkpoint(.TRUE., 'Invalid value for Tmax inserted (Insert Tmax>0)', 'e')
        STOP 
    ENDIF

    !simulation parameters print
    PRINT *, "================================================================="
    PRINT *, "SIMULATION PARAMETERS: "
    PRINT *, "___________ SPACE ___________"
    PRINT *, "-N_x:", N_x
    PRINT *, "-Lmin:  ", Lmin
    PRINT *, "-Lmax:  ", Lmax
    PRINT *, "-dx:    ",(Lmax-Lmin)/N_x
    PRINT *, "___________ TIME ___________"
    PRINT *, "-N_t:", N_t
    PRINT *, "-Tmin:  ", Tmin
    PRINT *, "-Tmax:  ", Tmax
    PRINT *, "-dx:    ", (Tmax-Tmin)/N_t
    PRINT *, "___________ SIMULATION ___________"
    PRINT *, "-m:     ", m
    PRINT *, "-omega: ", omega
    PRINT *, "================================================================="

    !================================================
    !MEMORY ALLOCATION
    !================================================
    
    ALLOCATE(t_latt(N_t+1))
    ALLOCATE(ham_op(N_x+1,N_x+1))
    ALLOCATE(eigs(N_x+1))
    ALLOCATE(wf(N_x+1))
    ALLOCATE(total_evolution_x(N_x+1, N_t+1))
    ALLOCATE(total_evolution_p(N_x+1, N_t+1))
    !expectation value of position
    ALLOCATE(x_latt(N_x+1))
    ALLOCATE(expval_x_vector(N_t+1))
    ALLOCATE(expval1sigma_x_vector(N_t+1))
    !expectation value of momentum
    ALLOCATE(p_latt(N_t+1))
    ALLOCATE(expval_p_vector(N_t+1))
    ALLOCATE(expval1sigma_p_vector(N_t+1))

    !================================================
    !LATTICES
    !================================================
    
    !discretization variables
    dx=(Lmax-Lmin)/N_x
    dt=(Tmax-Tmin)/N_t

    !defining the position lattice
    DO aa = 1, N_x+1
        x_latt(aa) = Lmin + (aa-1)*dx
    ENDDO

    !defining the time lattice
    DO aa = 1, N_t+1
        t_latt(aa) = Tmin + (aa-1)*dt
    ENDDO

    !defining the momentum lattice
    PI = 2.0d0*ASIN(1.0d0)

    DO aa=1,INT(N_x/2)
        p_latt(aa) = 2.0d0*PI*(aa)/(Lmax-Lmin)*dx
    ENDDO

    DO aa=INT(N_x/2)+1,N_x+1
        p_latt(aa) = 2.0d0*PI*(aa-N_x-1)/(Lmax-Lmin)*dx
    ENDDO

    !================================================

    !initialize hamiltonian of harmonic oscillator
    ham_op = harmosc_H_op_init(N_x, Lmin, Lmax, m, omega)

    !find eigenstates and save the ground state
    CALL diagonalize_H_op(ham_op, eigs)

    !print *, eigs
    DO aa = 1, SIZE(ham_op, 1)
        wf(aa) = COMPLEX(ham_op(aa,1), 0.0d0)
    ENDDO

    DO aa = 1, SIZE(wf)
        wf(aa) = wf(aa)*SQRT(N_x/(Lmax-Lmin))
    ENDDO
    
    DO aa = 1, N_t + 1

        expval_x = 0.0d0
        expval1sigma_x = 0.0d0
        expval_p = 0.0d0
        expval1sigma_p = 0.0d0

        CALL evolve_dt(wf, dt, Lmin, Lmax, m, omega, t_latt(aa), T)

        total_evolution_x(:,aa) = wf
        total_evolution_p(:,aa) = myFFTW(wf)

        DO bb = 1, N_x + 1
            expval_x = expval_x + REAL(total_evolution_x(bb,aa)*CONJG(total_evolution_x(bb,aa)))*x_latt(bb)*dx
            expval_p = expval_p + REAL(total_evolution_p(bb,aa)*CONJG(total_evolution_p(bb,aa)))*p_latt(bb)*dx
        ENDDO

        DO bb = 1, N_x + 1
            expval1sigma_x = expval1sigma_x + REAL(total_evolution_x(bb,aa)*CONJG(total_evolution_x(bb,aa)))*(x_latt(bb)-expval_x)**2*dx
            expval1sigma_p = expval1sigma_p + REAL(total_evolution_p(bb,aa)*CONJG(total_evolution_p(bb,aa)))*(p_latt(bb)-expval_p)**2*dx
        ENDDO

        expval_x_vector(aa) = expval_x
        expval_p_vector(aa) = expval_p
        expval1sigma_x_vector(aa) = SQRT(expval1sigma_x)
        expval1sigma_p_vector(aa) = SQRT(expval1sigma_p)

    ENDDO
 
    !PRINT TO FILE: color map for the probability distribution of the wavefunction evolution
    OPEN(2, file='x.dat', status='replace')
    DO aa=1, SIZE(x_latt,1)
        WRITE(2,*) x_latt(aa)
    END DO
    CLOSE(2)

    OPEN(2, file='p.dat', status='replace')
    DO aa=1, SIZE(p_latt,1)
        WRITE(2,*) p_latt(aa)
    END DO
    CLOSE(2)

    OPEN(2, file='y.dat', status='replace')
    DO aa=1,size(t_latt,1)
        WRITE(2,*) t_latt(aa)
    END DO
    CLOSE(2)

    OPEN(2, file='zx.dat', status='replace')
    DO bb = 1, size(total_evolution_x,2)
        WRITE(2,*) REAL(total_evolution_x(:,bb)*CONJG(total_evolution_x(:,bb)))
    END DO
    CLOSE(2)

    OPEN(2, file='zp.dat', status='replace')
    DO bb = 1, size(total_evolution_p,2)
        WRITE(2,*) REAL(total_evolution_p(:,bb)*CONJG(total_evolution_p(:,bb)))
    END DO
    CLOSE(2)

    OPEN(2, file='expvalue1sigma.dat', status='replace')
    DO aa = 1, size(expval1sigma_x_vector,1)
        WRITE(2,*) expval1sigma_x_vector(aa)
    END DO
    CLOSE(2)

    OPEN(2, file='expvalue.dat', status='replace')
    DO aa = 1, size(expval_x_vector,1)
        WRITE(2,*) expval_x_vector(aa)
    END DO
    CLOSE(2)

    OPEN(2, file='expvaluep1sigma.dat', status='replace')
    DO aa = 1, size(expval1sigma_p_vector,1)
        WRITE(2,*) expval1sigma_p_vector(aa)
    END DO
    CLOSE(2)

    OPEN(2, file='expvaluep.dat', status='replace')
    DO aa = 1, size(expval_p_vector,1)
        WRITE(2,*) expval_p_vector(aa)
    END DO
    CLOSE(2)

END PROGRAM time_dep_schrod