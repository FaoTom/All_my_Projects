import subprocess as sp
import numpy as np

N_vector = [2, 3, 4, 5, 6, 7, 8, 9 , 10, 11, 12]

llambda =  1

src = 'w7.f90'
exe = 'w7.out'

sp.run(['gfortran',
            src,
            '-I/usr/local/Cellar/fftw/3.3.10/include', 
            '-L/usr/local/Cellar/fftw/3.3.10/lib',
            'qm_utils.f90', 
            'my_mat_mul.f90',
            'error_handling.f90',
            'dc_matrix.f90',
            '-lfftw3',
            '-lfftw3f',
            '-llapack',
            '-lm',
            '-ffree-line-length-0',
            '-w',
            '-o',
            exe]
            )

ham = np.empty(len(N_vector), dtype=float)
eig   = np.empty(len(N_vector), dtype=float)

for nn, N in enumerate(N_vector):
    output = sp.Popen( [('./' + exe), '-custom', str(int(N)), str(llambda)], stdout=sp.PIPE ).communicate()[0]
    times = output.decode('utf-8').split()
    
    ham[nn] = float(times[0])
    eig[nn] = float(times[1])

np.savetxt('data/ham_times'+str(N)+'.txt',ham)
np.savetxt('data/eig_times'+str(N)+'.txt',eig)
