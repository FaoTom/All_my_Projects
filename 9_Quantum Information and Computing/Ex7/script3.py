import subprocess as sp
import numpy as np

N = 10

lambda_vector = np.linspace(0,3,200)

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

diff = np.empty(len(lambda_vector), dtype=float)

for ll, llambda in enumerate(lambda_vector):
    output = sp.Popen( [('./' + exe), '-custom', str(int(N)), str(llambda)], stdout=sp.PIPE ).communicate()[0]
    times = output.decode('utf-8').split()
    
    diff[ll] = float(times[0])

np.savetxt('data/diff'+str(N)+'.txt',diff)
