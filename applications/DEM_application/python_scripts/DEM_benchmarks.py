from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
print("import main_script")
import main_script

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import plot_variables                # Related to benchmarks in Chung, Ooi
import DEM_benchmarks_class as DBC   # Related to benchmarks in Chung, Ooi
#from glob import glob

sys.path.insert(0,'')
# DEM Application

benchmark_number, nodeplotter = int(sys.argv[1]), 0
benchmark = getattr(DBC, 'Benchmark' + str(benchmark_number))()

listDISCONT   = list(range(1,12))
listROLLFR    = list(range(12,13))
listDEMFEM    = list(range(13,18))
listCONT      = list(range(20,27))
listDISclZHAO = [30,32]
listDISclRK   = [31,33]

'''
listDISCONT   = [] #list(range(1,12))
listROLLFR    = [] #list(range(12,13))
listDEMFEM    = [] #list(range(13,18))
listCONT      = [] #list(range(20,27))
listDISclZHAO = [] #[30,32]
listDISclRK   = [] #[31,33]
'''
if benchmark_number in listDISCONT:
    nodeplotter = 1
    import DEM_explicit_solver_var as DEM_parameters
elif benchmark_number in listROLLFR:
    import DEM_explicit_solver_var_ROLLFR as DEM_parameters
elif benchmark_number in listDEMFEM:
    import DEM_explicit_solver_var_DEMFEM as DEM_parameters
elif benchmark_number in listCONT:
    import DEM_explicit_solver_var_CONT as DEM_parameters
elif benchmark_number == 27:
    import DEM_explicit_solver_var_UCS as DEM_parameters
elif benchmark_number == 28:
    import DEM_explicit_solver_var_PENDULO3D as DEM_parameters
elif benchmark_number in listDISclZHAO:
    import DEM_explicit_solver_var_DISclZHAO as DEM_parameters
elif benchmark_number in listDISclRK:
    import DEM_explicit_solver_var_DISclRK as DEM_parameters
else:
    print('Benchmark number does not exist')
    sys.exit()

class Solution(main_script.Solution):

    def __init__(self):
        super().__init__()
        os.chdir('..')
        self.main_path = os.getcwd()
        print("-----execute init")
        self.dt = self.output_time_step = self.number_of_coeffs_of_restitution = self.number_of_points_in_the_graphic = 0
        self.PreviousOperations()
         
    def Initialize(self):
        print("execute initialize") 
        super().Initialize()
        
    def SetInitialData(self, iteration, coeff_of_restitution_iteration):
        benchmark.set_initial_data(self.spheres_model_part, self.rigid_face_model_part, iteration, self.number_of_points_in_the_graphic, coeff_of_restitution_iteration)

    def GetMpFilename(self):
        print("getmpfilename-benchmark")
        return 'benchmark' + str(benchmark_number) + "DEM"
    
    def GetInletFilename(self):
        return 'benchmarkDEM_Inlet'
      
    def GetProblemTypeFilename(self):
        return 'benchmark' + str(benchmark_number)

    def BeforeSolveOperations(self, time):        
        benchmark.ApplyNodalRotation(self.spheres_model_part, time, self.dt)

    def BeforePrintingOperations(self, time):
        benchmark.generate_graph_points(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part, time, self.output_time_step, self.dt)

    def BeforeFinalizeOperations(self):
        benchmark.get_final_data(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part)

    def InitializeTimeStep(self):
        return self.dt

    def PreviousOperations(self):
        DEM_parameters.FinalTime, self.dt, self.output_time_step, self.number_of_points_in_the_graphic, self.number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
        DEM_parameters.problem_name = 'benchmark' + str(benchmark_number)

    def BeforeRunMainTemporalLoopOperations(self, iteration):
        print("Computing points in the curve...", 1 + self.number_of_points_in_the_graphic - iteration, "point(s) left to finish....",'\n')
        list_of_nodes_ids = [1]
        if nodeplotter:
            os.chdir(self.main_path)
            self.plotter = plot_variables.variable_plotter(self.spheres_model_part, list_of_nodes_ids)
            self.tang_plotter = plot_variables.tangential_force_plotter(self.spheres_model_part, list_of_nodes_ids, iteration)

    def FinalizeTimeStep(self, time):
            if nodeplotter:
                os.chdir(self.main_path)
                self.plotter.plot_variables(time) #Related to the benchmark in Chung, Ooi
                self.tang_plotter.plot_tangential_force(time)

    def Run(self, iteration):  
        print("execute run")      
        self.Initialize()
        self.BeforeRunMainTemporalLoopOperations(iteration)
        self.RunMainTemporalLoop()
        self.BeforeFinalizeOperations()
        self.Finalize()
        self.CleanUpOperations()

    def AdditionalFinalizeOperations(self):
        if nodeplotter:
            os.chdir(self.main_path)
            self.plotter.close_files()
            self.tang_plotter.close_files()

    def RunBenchmarks(self):
        print("RunBenchmarks")
        for coeff_of_restitution_iteration in range(1, self.number_of_coeffs_of_restitution + 1):
            print("coeff_of_restitution_iteration")
            for iteration in range(1, self.number_of_points_in_the_graphic + 1):
                print("iteration")
                Solution().Run(iteration)
            benchmark.print_results(self.number_of_points_in_the_graphic, self.dt)
        
    def CleanUpOperations(self):
        print("execute CleanUpOperations") 
        super().CleanUpOperations()
    
Solution().RunBenchmarks()