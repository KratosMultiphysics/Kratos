import matplotlib.pyplot as plt
import math
import chandelier_parameters as ch_pp
import chandelier as ch
import DEM_explicit_solver_var as DEM_parameters
import os

def Radius(a):
    a_2 = [x ** 2 for x in a]
    return math.sqrt(sum(a_2))

def ComputeL1NormOfError(method, j_var, k_var):
    E = 0.
    E_r = 0.
    main_path = os.getcwd()
    with open(main_path + '/coors_' + str(method) + str(j_var) + str(k_var) + '.txt','r') as f_coors:
        ch_pp.include_history_force = 1
        sim = ch.AnalyticSimulator(ch_pp)
        i = 0
        for line in f_coors:
            i += 1
            coor = line.split()
            coor = [float(comp) for comp in coor]
            exact_coor = [0.] * 3
            sim.CalculatePosition(exact_coor, coor[0])
            x = exact_coor[0] * ch_pp.R
            y = exact_coor[1] * ch_pp.R
            r = Radius(coor[1:3])
            exact_r = Radius([x, y])
            delta_coor = [coor[1] - x, coor[2] - y]
            #E = max(E, Radius(delta_coor) / exact_r)
            #E_r = max(E_r, (r - exact_r) / exact_r)
            E += Radius(delta_coor) / exact_r
            E_r += (r - exact_r) / exact_r
    return E/i, E_r/i

#sim = ch.AnalyticSimulator(ch_pp)
#coors = [None] * 3
#sim.CalculatePosition(coors, 100.0)
#x_final = coors[0] * ch_pp.R
#y_final = coors[1] * ch_pp.R
#r_final = Radius([x_final, y_final])

final_positions = []
with open("final_positions_decreasing_DT", 'r') as f:
    for line in f:
        position = line.split()
        position = [float(position[0]), float(position[1])]
        final_positions.append(position)

positions_table = []

n_exponentials = 12
n_doubling = 8

#for m in range(n_exponentials):
    #if m == 0:
        #print('Reading Daitche...')
    #else:
        #print('Reading Hinsberg, m = ' + str(m - 1) + '...')
    #positions_line = []
    #for k in range(n_doubling):
        #positions_line.append(final_positions[m * n_doubling + k])
    #positions_table.append(positions_line)

r_table = [[0.0 for k in range(len(positions_table[0]))] for m in range(len(positions_table))]
E = [[0. for k in range(n_doubling)] for m in range(n_exponentials)]
E_r = [[0. for k in range(n_doubling)] for m in range(n_exponentials)]

for m in range(0, n_exponentials):
    if m == 0:
        method = 2
        print('Calculating error for Daitche...')
    else:
        method = 4
        print('Calculating error for Hinsberg, m = ' + str(m - 1) + '...')
    for k in range(0, n_doubling):
        error, radius_error = ComputeL1NormOfError(method, max(m - 1, 0), k)
        E[m][k] = error
        E_r[m][k] = radius_error

dts = [DEM_parameters["MaxTimeStep"].GetDouble() / 2 ** k for k in range(n_doubling)]

#E_r = [[abs(r_final - r_table[m][k]) / r_final for k in range(n_doubling)] for m in range(n_exponentials)]
#E   = [[Radius([positions_table[m][k][0] - x_final, positions_table[m][k][1] - y_final]) / r_final for k in range(n_doubling)] for m in range(len(r_table))]
index_short = len([value for value in dts if value < 2])


#from matplotlib import rc
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
font_size = 40
ax.tick_params(axis='x', labelsize=0.5 * font_size)
ax.tick_params(axis='y', labelsize=0.5 * font_size)
ax.tick_params(axis='both', which='major', pad=10)
fig.tight_layout()

m = 0

for line in reversed(E):
    if m == n_exponentials - 1:
        plt.plot(dts, line, '--', marker='o', markersize=10, label='Daitche', color='k')
    else:
        plt.plot(dts, line, '-o', markersize=10, label='m = ' + str(10 - m), linewidth=0.5*(11-m), ms=10)
    m += 1


plt.xlabel('$h$', fontsize = 0.75 * font_size)
plt.ylabel('$E(100)$', fontsize = 0.75 * font_size)
plt.legend(loc='upper left', prop={'size':25})
plt.semilogx()
plt.semilogy()
figure = plt.gcf() # get current figure
figure.set_size_inches(14, 11)
plt.savefig('EvsDt.pdf', dpi = 1000, bbox_inches = 'tight')
plt.show()
