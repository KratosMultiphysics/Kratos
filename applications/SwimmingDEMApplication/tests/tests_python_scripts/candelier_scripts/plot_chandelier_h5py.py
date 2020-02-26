import numpy as np
import h5py
import matplotlib.pyplot as plt


methods = {'Daitche'}
time_steps = {0.01, 0.02, 0.04, 0.08, 0.16}
N_qs = {1, 2, 4, 8, 16, 32, 64, 128}
quadrature_orders = {2}
numbers_of_exp = {1}
tws = {1}
cases = []

class Case:
    def __init__(self, method, dt, Nq, quad_order, tw = None, m = None):
        self.method = method
        self.dt = dt
        self.Nq = Nq
        self.quad_order = quad_order
        self.tw = tw
        self.m = m
        self.MakeCode()

    def MakeCode(self):
        self.code = '/' + self.method + '_' + 'dt=' + str(self.dt) + '_' + 'Nq=' + str(self.Nq) + '_' + 'quadrature_order=' + str(self.quad_order)
        if self.method == 'Hinsberg':
            self.code += '_tw=' + str(self.tw) + '_m=' + str(self.m)


for method in methods:
    for dt in time_steps:
        for Nq in N_qs:
            for order in quadrature_orders:
                for tw in tws:
                    for m in numbers_of_exp:
                        if Nq == 1 or dt == 0.01:
                            cases.append(Case(method, dt, Nq, order, tw, m))

cases.sort(key=lambda x: (x.dt, x.Nq), reverse=False)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
font_size = 40
ax.tick_params(axis='x', labelsize=0.5 * font_size)
ax.tick_params(axis='y', labelsize=0.5 * font_size)
ax.tick_params(axis='both', which='major', pad=10)
fig.tight_layout()

with h5py.File('candelier_results.h5py', 'r') as f:
    # result_code = '/Hinsberg_dt=0.01_Nq=1_quadrature_order=2_tw=1_m=0'
    for case in cases:
        times = np.array(f[case.code + '/time'])
        errors = np.array(f[case.code + '/E(t)'])

        if case.Nq == 1 and case.dt == 0.01:
            line_width = 4
            plt.plot(times, errors, color = 'k', label = 'Nq = ' + str(case.Nq), linewidth = line_width)
        elif case.dt == 0.01:
            line_width = 2
            plt.plot(times, errors, linestyle = '--', dashes = (0.3 * case.Nq, case.Nq ** 0.5), color = 'k', label = 'Nq = ' + str(case.Nq), linewidth = line_width)
        else:
            line_width = 1
            times_dt = int(case.dt / 0.01)
            plt.plot(times, errors, linestyle = '--', dashes = (0.3 * times_dt, times_dt ** 0.5), color = 'r', label = '$\Delta t \\times $' + str(times_dt), linewidth = line_width)

plt.xlabel('$t$', fontsize = 0.75 * font_size)
plt.ylabel('$E(t)$', fontsize = 0.75 * font_size)
plt.semilogy()
plt.legend(loc = 'lower right', prop={'size':20}, ncol = 2)
fig.set_size_inches(14, 11)
plt.savefig('EvsSubstepping.pdf', dpi = 1000, format='pdf', bbox_inches = 'tight')
plt.savefig('EvsSubstepping.eps', dpi = 1000, format='eps', bbox_inches = 'tight')
plt.show()
