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
         
    def Initialize(self):
        print("execute initialize") 
        DEM_parameters.problem_name = 'benchmark' + str(benchmark_number)
        super().Initialize()
        

        print("Computing points in the curve...", 1 + self.number_of_points_in_the_graphic - self.iteration, "point(s) left to finish....",'\n')
        list_of_nodes_ids = [1]
        if nodeplotter:
            os.chdir(self.main_path)
            self.plotter = plot_variables.variable_plotter(self.spheres_model_part, list_of_nodes_ids)
            self.tang_plotter = plot_variables.tangential_force_plotter(self.spheres_model_part, list_of_nodes_ids, self.iteration)

    def ReadModelParts(self):
        super().ReadModelParts()
        benchmark.set_initial_data(self.spheres_model_part, self.rigid_face_model_part, self.iteration, self.number_of_points_in_the_graphic, coeff_of_restitution_iteration)

    def GetMpFilename(self):
        print("getmpfilename-benchmark")
        return 'benchmark' + str(benchmark_number) + "DEM"
    
    def GetInletFilename(self):
        return 'benchmarkDEM_Inlet'
      
    def GetProblemTypeFilename(self):
        return 'benchmark' + str(benchmark_number)

    def BeforeSolveOperations(self, time): 
        super().BeforeSolveOperations(time)     
        benchmark.ApplyNodalRotation(self.spheres_model_part, time, self.dt)

    def BeforePrintingOperations(self, time):
        super().BeforePrintingOperations(time)
        benchmark.generate_graph_points(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part, time, self.output_time_step, self.dt)

    def Finalize(self):
        benchmark.get_final_data(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part)
        super().Finalize()
        if nodeplotter:
            os.chdir(self.main_path)
            self.plotter.close_files()
            self.tang_plotter.close_files()

    def FinalizeTimeStep(self, time):
        super().FinalizeTimeStep(time)
        if nodeplotter:
            os.chdir(self.main_path)
            self.plotter.plot_variables(time) #Related to the benchmark in Chung, Ooi
            self.tang_plotter.plot_tangential_force(time)   
        
    def CleanUpOperations(self):
        print("execute CleanUpOperations") 
        DBC.delete_archives() #.......Removing some unuseful files 
        super().CleanUpOperations()
    

final_time, dt, output_time_step, number_of_points_in_the_graphic, number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
for coeff_of_restitution_iteration in range(1, number_of_coeffs_of_restitution + 1):
    for iteration in range(1, number_of_points_in_the_graphic + 1):
        sol = Solution()
        sol.iteration = iteration
        sol.dt = dt
        DEM_parameters.FinalTime = final_time
        sol.output_time_step = output_time_step
        print(DEM_parameters.FinalTime)        
        sol.number_of_points_in_the_graphic = number_of_points_in_the_graphic
        sol.number_of_coeffs_of_restitution = number_of_coeffs_of_restitution
        sol.Run()
    benchmark.print_results(number_of_points_in_the_graphic, dt)