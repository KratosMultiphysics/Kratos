from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
kratos_path = '../../../..'
import sys
sys.path.append(kratos_path)

f1 = open("../../../../applications/DEM_application/custom_problemtype/DEM_explicit_solver.gid/spheric_particle_script.py", 'r')
f2 = open("benchmark_lines.py", 'r')
fileOne = f1.readlines()
fileTwo = f2.readlines()
f1.close()
f2.close()
outFile = open("two_balls_no_damp.py", 'w')
x = 0
for line in fileOne:
    if line.find("BENCHMARK") == -1:
        outFile.write(line)
    else:
        outFile.write(fileTwo[x])
        x += 1
outFile.close()
