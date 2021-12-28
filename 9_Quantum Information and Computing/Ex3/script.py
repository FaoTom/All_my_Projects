import subprocess as sp
import numpy as np

N_min = 20
N_max = 2020
points  = 20

N_vector = np.arange(N_min, N_max, points)

src = 'w3_ex1.f90'
exe = 'w3_ex1.o'
sp.run(['gfortran', src , 'my_mat_mul.f90','error_handling.f90','-o', exe])

matmul_col = np.empty(len(N_vector), dtype=float)
matmul_row = np.empty(len(N_vector), dtype=float)
matmul_for = np.empty(len(N_vector), dtype=float)

for nn, N in enumerate(N_vector):
    output = sp.Popen( [('./' + exe), str(N)], stdout=sp.PIPE ).communicate()[0]
    times = output.decode('utf-8').split()
    
    matmul_col[nn] = float(times[0])
    matmul_row[nn] = float(times[1])
    matmul_for[nn] = float(times[2])

np.savetxt('data/matmul_bycol.txt',matmul_col)
np.savetxt('data/matmul_byrow.txt',matmul_row)
np.savetxt('data/matmul_fortran.txt',matmul_for)