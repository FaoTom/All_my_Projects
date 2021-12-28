import subprocess as sp
import numpy as np

N_body_min = 10
N_body_max = 1000000
points  = 200
local_dimension = 5


#N_body_vector = np.linspace(N_body_min, N_body_max, points)
N_body_vector = np.arange(10,1000000,10000)
#N_body_vector = [2, 5, 10, 100, 200, 500, 1000, 2000, 5000,  10000, 20000, 50000, 100000, 200000, 5000000, 1000000]
print(len(N_body_vector))

src = 'w6.f90'
exe = 'w6.out'

sp.run(['gfortran',
            src,
            '-I/usr/local/Cellar/fftw/3.3.10/include', 
            '-L/usr/local/Cellar/fftw/3.3.10/lib',
            'qm_utils.f90', 
            'my_mat_mul.f90',
            'error_handling.f90',
            '-lfftw3',
            '-lfftw3f',
            '-llapack',
            '-lm',
            '-ffree-line-length-0',
            '-w',
            '-o',
            exe]
            )

separable = np.empty(len(N_body_vector), dtype=float)
generic   = np.empty(len(N_body_vector), dtype=float)


for nn, N in enumerate(N_body_vector):
    output = sp.Popen( [('./' + exe), '-custom', str(int(N)), str(local_dimension)], stdout=sp.PIPE ).communicate()[0]
    times = output.decode('utf-8').split()
    
    separable[nn] = float(times[0])
    #generic[nn]   = float(times[0])

np.savetxt('data/sep_'+str(local_dimension)+'.txt',separable)
#np.savetxt('data/gen_'+str(local_dimension)+'txt',generic)
