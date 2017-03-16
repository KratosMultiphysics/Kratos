import numpy as np
import h5py
import matplotlib.pyplot as plt


def MakeCode(method, dt, Nq, quad_order, tw = None, m = None):
    code = '/' + method + '_' + 'dt=' + str(dt) + '_' + 'Nq=' + str(Nq) + '_' + 'quadrature_order=' + str(quad_order)
    if method == 'Hinsberg':
        code += '_tw=' + str(tw) + '_m=' + str(m)
    return code


methods = {'Daitche'}
time_steps = {0.01}
N_qs = {1, 2, 4, 8, 16, 32, 64, 128}
quadrature_orders = {2}
numbers_of_exp = {0}
tws = {1}
codes = set()

for method in methods:
    for dt in time_steps:
        for Nq in N_qs:
            for order in quadrature_orders:
                for tw in tws:
                    for m in numbers_of_exp:
                        codes.add(MakeCode(method, dt, Nq, order, tw, m))

with h5py.File('candelier_results.h5py', 'r') as f:
    # result_code = '/Hinsberg_dt=0.01_Nq=1_quadrature_order=2_tw=1_m=0'
    for code in codes:
        times = np.array(f[code + '/time'])
        errors = np.array(f[code + '/E(t)'])
        plt.plot(times, errors)

plt.semilogy()
plt.show()
