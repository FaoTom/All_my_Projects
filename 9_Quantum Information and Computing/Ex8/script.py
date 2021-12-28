import subprocess as sp
import numpy as np

N = 4

lambda_vector = np.linspace(0,3,200)

src = 'RSRG.f90'
exe = 'RSRG.out'

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

gs = np.empty(len(lambda_vector), dtype=float)
it = np.empty(len(lambda_vector), dtype=float)


for ll, llambda in enumerate(lambda_vector):
    output = sp.Popen( [('./' + exe), '-custom', str(int(N)), str(llambda)], stdout=sp.PIPE ).communicate()[0]
    times = output.decode('utf-8').split()
    
    it[ll] = float(times[0])
    gs[ll] = float(times[1])

np.savetxt('data/gs'+str(N)+'.txt',gs)
np.savetxt('data/it'+str(N)+'.txt',it)

