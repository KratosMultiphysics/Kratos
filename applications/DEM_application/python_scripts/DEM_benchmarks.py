from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
import main_script

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import plot_variables                # Related to benchmarks in Chung, Ooi
import DEM_benchmarks_class as DBC   # Related to benchmarks in Chung, Ooi

sys.path.insert(0,'')
benchmark_number, nodeplotter = int(sys.argv[1]), 0
benchmark = getattr(DBC, 'Benchmark' + str(benchmark_number))()

listDISCONT   = list(range(1,12))
listROLLFR    = list(range(12,13))
listDEMFEM    = list(range(13,18))
listCONT      = list(range(20,27))
listDISclZHAO = [30,32]
listDISclRK   = [31,33]

if benchmark_number in listDISCONT:
    nodeplotter = 1
    import DEM_explicit_solver_var_DISCONT as DEM_parameters
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
        super().__init__(DEM_parameters)
        os.chdir('..')
        self.main_path = os.getcwd()

    def GetProblemTypeFilename(self):
        return benchmark

    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):        
        return ModelPartIO(modelpart)

    def SetSolverStrategy(self):
        # Strategy object
        if (DEM_parameters.ElementType == "SphericPartDEMElement3D" or DEM_parameters.ElementType == "CylinderPartDEMElement2D"):
            import sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "SphericContPartDEMElement3D" or DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
            import continuum_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "ThermalSphericContPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "ThermalSphericPartDEMElement3D"):
            import thermal_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "SinteringSphericConPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "IceContPartDEMElement3D"):
            import ice_continuum_sphere_strategy as SolverStrategy
        else:
            self.KRATOSprint('Error: Strategy unavailable. Select a different scheme-element')

        return SolverStrategy    
         
    def SetSolver(self):
        return self.solver_strategy.ExplicitStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, DEM_parameters, self.procedures)

    def SetFinalTime(self):        
        self.final_time = final_time  
        #return self.final_time

    def Setdt(self):       
        self.dt = dt
        #print("setdt", self.dt )
        #return self.dt

    def Initialize(self):
        DEM_parameters.problem_name = 'benchmark' + str(benchmark_number)
        #self.final_time = slt.final_time  
        #self.dt = slt.dt
        #self.graph_print_interval = slt.graph_print_interval
        super().Initialize()           

        print("Computing points in the curve...", 1 + self.number_of_points_in_the_graphic - self.iteration, "point(s) left to finish....",'\n')
        list_of_nodes_ids = [benchmark_number]
        if nodeplotter:
            os.chdir(self.main_path)
            self.plotter = plot_variables.variable_plotter(self.spheres_model_part, list_of_nodes_ids)
            self.tang_plotter = plot_variables.tangential_force_plotter(self.spheres_model_part, list_of_nodes_ids, self.iteration)

    def ReadModelParts(self):
        super().ReadModelParts()
        benchmark.set_initial_data(self.spheres_model_part, self.rigid_face_model_part, self.iteration, self.number_of_points_in_the_graphic, coeff_of_restitution_iteration)       

    def GetMpFilename(self):
        return 'benchmark' + str(benchmark_number) + "DEM"
    
    def GetInletFilename(self):
        return 'benchmarkDEM_Inlet'
        #return 'benchmark' + str(benchmark_number) + "DEM_Inlet"

    def GetFemFilename(self):
        return 'benchmark' + str(benchmark_number) + "DEM_FEM_boundary"   

    def GetClusterFilename(self):
        return 'benchmark' + str(benchmark_number) + "DEM_Clusters"   
      
    def GetProblemTypeFilename(self):
        return 'benchmark' + str(benchmark_number)

    def BeforeSolveOperations(self, time): 
        super().BeforeSolveOperations(time)     
        benchmark.ApplyNodalRotation(time, self.dt, self.spheres_model_part)

    def BeforePrintingOperations(self, time):
        super().BeforePrintingOperations(time)
        self.Setdt()
        benchmark.generate_graph_points(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part, time, self.graph_print_interval, self.dt)        
    
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
        print("running CleanUpOperations") 
        #DBC.delete_archives() #.......Removing some unuseful files 
        super().CleanUpOperations() 
                  

final_time, dt, graph_print_interval, number_of_points_in_the_graphic, number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
for coeff_of_restitution_iteration in range(1, number_of_coeffs_of_restitution + 1):
    for iteration in range(1, number_of_points_in_the_graphic + 1):
        slt = Solution()
        slt.iteration = iteration
        slt.dt = dt
        slt.final_time = final_time    
        slt.graph_print_interval = graph_print_interval     
        slt.number_of_points_in_the_graphic = number_of_points_in_the_graphic
        slt.number_of_coeffs_of_restitution = number_of_coeffs_of_restitution        
        slt.Run()
    benchmark.print_results(number_of_points_in_the_graphic, dt)
#DBC.delete_archives() #.......Removing some unuseful files     