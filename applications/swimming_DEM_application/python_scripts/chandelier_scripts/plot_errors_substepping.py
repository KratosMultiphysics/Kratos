import matplotlib.pyplot as plt
import math
import chandelier_parameters as ch_pp
import chandelier as ch
import DEM_explicit_solver_var as DEM_parameters

def Radius(a):
    a_2 = [x ** 2 for x in a]
    return math.sqrt(sum(a_2))


sim = ch.AnalyticSimulator(ch_pp)
coors = [None] * 3
sim.CalculatePosition(coors, 100.0)
x_final = coors[0] * ch_pp.R
y_final = coors[1] * ch_pp.R
r_final = Radius([x_final, y_final])

final_positions = []
with open("final_positions_substepping", 'r') as f:
    for line in f:
        position = line.split()
        position = [float(position[0]), float(position[1])]
        final_positions.append(position)

positions_table = []

n_exponentials = 12
n_doubling = 6
print('len(final_positions)',len(final_positions))
for m in range(n_exponentials):
    positions_line = []
    for k in range(n_doubling):
        positions_line.append(final_positions[m * n_doubling + k])
    positions_table.append(positions_line)

r_table = [[0.0 for k in range(len(positions_table[0]))] for m in range(len(positions_table))]

for m in range(0, n_exponentials):
    for k in range(0, n_doubling):
        r_table[m][k] = Radius(positions_table[m][k])

dts_quad = [DEM_parameters["MaxTimeStep"].GetDouble() * 2 ** k for k in range(n_doubling)]

E_r = [[abs(r_final - r_table[m][k]) / r_final for k in range(n_doubling)] for m in range(n_exponentials)]
E   = [[Radius([positions_table[m][k][0] - x_final, positions_table[m][k][1] - y_final]) / r_final for k in range(n_doubling)] for m in range(len(r_table))]
index_short = len([value for value in dts_quad if value < 2])

m = 0
for line in reversed(E):
    if m == 11:
        plt.plot(dts_quad, line, '--', marker='o', markersize=10, label='Daitche', color='k')
    else:
        plt.plot(dts_quad, line, '-o', markersize=10, label='Hinsberg, m = ' + str(10 - m), linewidth=0.5*(11-m), ms=10)
    m += 1
plt.xlabel('$\mathrm{d}t_q$')
plt.ylabel('$E$')
plt.legend(loc='lower right')
plt.semilogx()
plt.semilogy()
plt.savefig('EvsNq.pdf')
plt.show()
