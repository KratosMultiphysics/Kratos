# import modules for TAU
import PyPara
import PyPrep
import PySolv
import shutil

# Definition of the parameter file
para_path='airfoil_Structured.cntl'
para_path_mod = para_path + ".mod"
shutil.copy(para_path, para_path_mod)

# Init Tau python classes
Para = PyPara.Parafile(para_path_mod)
Prep = PyPrep.Preprocessing(para_path_mod)
Solver = PySolv.Solver(para_path_mod)

Prep.run(write_dualgrid=0,free_primgrid=False)
Solver.init(verbose = 1, reset_steps = True, n_time_steps = 1)
Solver.outer_loop()
Solver.output()
Solver.finalize()
Para.free_parameters()
tau("exit")