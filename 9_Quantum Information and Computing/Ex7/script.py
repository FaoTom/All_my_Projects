import subprocess as sp
import numpy as np

N = 2

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

k1 = np.empty(len(lambda_vector), dtype=float)
k2   = np.empty(len(lambda_vector), dtype=float)
k3   = np.empty(len(lambda_vector), dtype=float)
k4   = np.empty(len(lambda_vector), dtype=float)

for ll, llambda in enumerate(lambda_vector):
    output = sp.Popen( [('./' + exe), '-custom', str(int(N)), str(llambda)], stdout=sp.PIPE ).communicate()[0]
    times = output.decode('utf-8').split()
    
    k1[ll] = float(times[0])
    k2[ll] = float(times[1])
    k3[ll] = float(times[2])
    k4[ll] = float(times[3])

np.savetxt('data/k1_N'+str(N)+'.txt',k1)
np.savetxt('data/k2_N'+str(N)+'.txt',k2)
np.savetxt('data/k3_N'+str(N)+'.txt',k3)
np.savetxt('data/k4_N'+str(N)+'.txt',k4)
