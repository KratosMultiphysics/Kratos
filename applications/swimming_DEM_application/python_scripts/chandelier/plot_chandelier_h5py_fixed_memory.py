import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

n_points_in_memory = 10
methods = {'Hinsberg'}
time_steps = [0.01, 0.005, 0.001, 0.0005]
Nqs = {1, 2, 4, 8, 16, 32}
quadrature_orders = {2}
numbers_of_exp = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
tws = set()
for Nq in Nqs:
    tws.add(Nq * n_points_in_memory)
cases_list = []

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

for i_plot in range(len(list(time_steps))):
    cases_list_i = []
    for method in methods:
        for order in quadrature_orders:
            for m in numbers_of_exp:
                cases_m = []
                dt = list(time_steps)[i_plot]
                for Nq in Nqs:
                    tw = 10 * Nq * dt
                    cases_m.append(Case(method, dt, Nq, order, tw, m))
                cases_m.sort(key=lambda x: x.Nq, reverse=False)
                cases_list_i.append(cases_m)
    cases_list_i.sort(key=lambda x: x[0].Nq, reverse=False)
    cases_list.append(cases_list_i)

fig = plt.figure()
font_size = 4
line_width = 1
gs1 = gridspec.GridSpec(4, 4)
gs1.update(wspace=0.025, hspace=0.05) # set the spacing between axes.
with h5py.File('candelier_results.h5py', 'r') as f:
    # result_code = '/Hinsberg_dt=0.01_Nq=1_quadrature_order=2_tw=1_m=0'
    for i_plot, ref_time_step in enumerate(list(time_steps)):


        for case_m in cases_list[i_plot]:
            m = case_m[0].m
            Nqs_m = []
            E20m = []
            for case in case_m:
                print(case.code)
                times = np.array(f[case.code + '/time'])
                errors = np.array(f[case.code + '/E(t)'])
                Nqs_m.append(case.Nq)
                E20m.append(errors[-1])

            ax = plt.subplot(2, 2, i_plot + 1)
            plt.semilogy()
            plt.xlabel('$N_q$', fontsize = 4 * font_size)
            plt.ylabel('$E(20)$', fontsize = 4 * font_size)
            ax.tick_params(axis='x', labelsize= 3 * font_size)
            ax.tick_params(axis='y', labelsize= 3 * font_size)
            ax.tick_params(axis='both', which='major', pad=4)
            ax.plot(Nqs_m, E20m, color = 'k', label = '$m =$' + str(m), linewidth = 0.5 * m)
            ax.set_title('$\\Delta t = $' + str(time_steps[i_plot]))
            handles, labels = ax.get_legend_handles_labels()

fig.subplots_adjust(wspace=0, hspace=0)
# plt.savefig('EvsNqFixedMemory.pdf', dpi = 1000, format='pdf', bbox_inches = 'tight')
# plt.savefig('EvsNqFixedMemory.eps', dpi = 1000, format='eps', bbox_inches = 'tight')
# plt.figlegend(handles, labels, loc = (0.2, - 0.01), ncol=4)
plt.legend( handles, labels, loc = 'lower center', bbox_to_anchor = (0,-0.12,1,1),
            bbox_transform = plt.gcf().transFigure, ncol = 5 )
plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
plt.savefig('EvsNqFixedMemory.pdf', dpi = 1000, format='pdf', bbox_inches = 'tight')
plt.savefig('EvsNqFixedMemory.eps', dpi = 1000, format='eps', bbox_inches = 'tight')
plt.show()
