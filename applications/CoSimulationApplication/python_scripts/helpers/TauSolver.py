# import modules for TAU
import PyPara, PyPrep, PySolv
import shutil, sys

control_signal = ""
while True:
    try:
        control_signal = sys.stdin.readline()
        if control_signal == "Initialize\n":
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
        elif control_signal == "SolveSolutionStep\n":
            Solver.outer_loop()
            Solver.output()
        elif control_signal == "Finalize\n":
            Solver.finalize()
            Para.free_parameters()
            tau("exit")
            break
        elif control_signal == "":
            pass
        else:
            raise Exception("Unknown control signal received: {}\n".format(control_signal))
    except IOError:
        raise IOError("TAU Solver cannot recieve signal")
