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
ch_pp.include_history_force = 0
sim = ch.AnalyticSimulator(ch_pp)
sim.CalculatePosition(coors, 100.0)
x_final_NH = coors[0] * ch_pp.R
y_final_NH = coors[1] * ch_pp.R
r_final_NH = Radius([x_final_NH, y_final_NH])
E_NH = Radius([abs(x_final_NH - x_final), abs(y_final - y_final_NH)]) / r_final
E_RNH = abs(r_final_NH - r_final) / r_final
final_positions = []
with open("final_positions_Dt=0.0025", 'r') as f:
    for line in f:
        position = line.split()
        position = [float(position[0]), float(position[1])]
        final_positions.append(position)

positions_table = []

n_exponentials = 11
n_doubling = 9

for m in range(n_exponentials):
    positions_line = []
    for k in range(n_doubling):
        positions_line.append(final_positions[m * n_doubling + k])
    positions_table.append(positions_line)

r_table = [[0.0 for k in range(len(positions_table[0]))] for m in range(len(positions_table))]

for m in range(0, n_exponentials):
    for k in range(0, n_doubling):
        r_table[m][k] = Radius(positions_table[m][k])

E_r = [[abs(r_final - r_table[m][k]) / r_final for k in range(n_doubling)] for m in range(n_exponentials)]
E   = [[Radius([positions_table[m][k][0] - x_final, positions_table[m][k][1] - y_final]) / r_final for k in range(n_doubling)] for m in range(n_exponentials)]
bytes_per_vector = 24
memory = [[(m + 5 * 2 ** k) * 24 for k in range(n_doubling)] for m in range(n_exponentials)]

m = 0
for line in reversed(E):
    plt.plot(memory[10-m], line, '-o', label='m = ' + str(10 - m), linewidth=0.5*(11-m), ms=10)
    m += 1
plt.xlabel('bytes stored per particle')
plt.ylabel('$E$')
plt.legend()
plt.semilogy()
plt.savefig('EvsMemory.pdf')
plt.show()

