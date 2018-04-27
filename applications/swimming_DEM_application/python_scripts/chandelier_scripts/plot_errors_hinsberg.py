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

t_wins = [5 * DEM_parameters["MaxTimeStep"].GetDouble() * 2 ** k for k in range(n_doubling)]

E_r = [[abs(r_final - r_table[m][k]) / r_final for k in range(n_doubling)] for m in range(n_exponentials)]
E   = [[Radius([positions_table[m][k][0] - x_final, positions_table[m][k][1] - y_final]) / r_final for k in range(n_doubling)] for m in range(len(r_table))]
index_short = len([value for value in t_wins if value < 2])


m = 0
for line in reversed(E_r):
    plt.plot(t_wins, line, '-o', label='m = ' + str(10 - m), linewidth=0.5*(11-m), ms=10)
    m += 1
plt.xlabel('r$t_w$')
plt.ylabel('r$E_R$')
plt.legend()
plt.semilogy()
plt.axhline(y = E_RNH, color='k',linestyle='--', linewidth = 2)
plt.savefig('E_r.pdf')
plt.clf()

m = 0
plt.clf()
for line in reversed(E_r):
    plt.plot(t_wins[:index_short], line[:index_short], '-o', label='m = ' + str(10 - m), linewidth=0.5*(11-m), ms=10)
    m += 1
plt.xlabel('r$t_w$')
plt.ylabel('r$E_R$')
plt.legend()
plt.semilogy()
plt.axhline(y = E_RNH, color='k',linestyle='--', linewidth = 2)
plt.savefig('E_r_short.pdf')

m = 0
for line in reversed(E):
    plt.plot(t_wins, line, '-o', label='m = ' + str(10 - m), linewidth=0.5*(11-m), ms=10)
    m += 1
plt.xlabel('r$t_w$')
plt.ylabel('r$E$')
plt.legend()
plt.semilogy()
plt.axhline(y = E_NH, color='k',linestyle='--', linewidth = 2)
plt.savefig('E.pdf')

m = 0
plt.clf()
for line in reversed(E):
    plt.plot(t_wins[:index_short], line[:index_short], '-o', label='m = ' + str(10 - m), linewidth=0.5*(11-m), ms=10)
    m += 1
plt.xlabel('r$t_w$')
plt.ylabel('r$E$')
plt.legend()
plt.semilogy()
plt.axhline(y = E_NH, color='k',linestyle='--', linewidth = 2)
plt.savefig('E_short.pdf')
plt.show()

