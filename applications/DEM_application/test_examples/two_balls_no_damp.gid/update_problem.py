from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
kratos_path = '../../../..'
import sys
import os
sys.path.append(kratos_path)

os.system("gid -n -b two_balls_no_damp.gid/batch.bch")

# os.system("mv ../../../../applications/DEM_application/test_examples/two_balls_no_damp.gid/two_balls_no_damp-3.dat ../../../../applications/DEM_application/test_examples/two_balls_no_damp.gid/DEM_explicit_solver_var.py")

# os.system("cat python_variables.py >> DEM_explicit_solver_var.py")

# os.system("rm ../../../../applications/DEM_application/test_examples/two_balls_no_damp.gid/two_balls_no_damp-1.dat ../../../../applications/DEM_application/test_examples/two_balls_no_damp.gid/two_balls_no_damp-2.dat ../../../../applications/DEM_application/test_examples/two_balls_no_damp.gid/two_balls_no_damp_aux.unix.bat")

# os.system("cp ../../../../applications/DEM_application/custom_problemtype/DEM_explicit_solver.gid/DEM_procedures.py ../../../../applications/DEM_application/test_examples/two_balls_no_damp.gid/DEM_procedures.py")
