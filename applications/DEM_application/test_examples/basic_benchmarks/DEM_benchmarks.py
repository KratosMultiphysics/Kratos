from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import main_script
import plot_variables                # Related to benchmarks in Chung, Ooi
import DEM_benchmarks_class as DBC   # Related to benchmarks in Chung, Ooi

sys.path.insert(0,'')
start = timer.time()
benchmark_number = int(sys.argv[1])
benchmark = getattr(DBC, 'Benchmark' + str(benchmark_number))()

listDISCONT   = list(range(1,12))
listROLLFR    = list(range(12,13))
listDEMFEM    = list(range(13,18))
listCONT      = list(range(20,27))
listDISclZHAO = [30,32]
listDISclRK   = [31,33]


class Solution(main_script.Solution):

    def LoadParametersFile(self):
        file_name = None
        if benchmark_number in listDISCONT:
            self.nodeplotter = True
            file_name = "ProjectParametersDISCONT.json"
        elif benchmark_number in listROLLFR:
            file_name = "ProjectParametersROLLFR.json"
        elif benchmark_number in listDEMFEM:
            file_name = "ProjectParametersDEMFEM.json"
        elif benchmark_number in listCONT:
            file_name = "ProjectParametersDEMCONT.json"
        elif benchmark_number == 27:
            file_name = "ProjectParametersUCS.json"
        #elif benchmark_number == 28:
        #    import DEM_explicit_solver_var_PENDULO3D as DEM_parameters  #disappeared?
        #    parameters_file = open("ProjectParametersDEM.json",'r')
        elif benchmark_number in listDISclZHAO:
            file_name = "ProjectParametersDISclZHAO.json"
        elif benchmark_number in listDISclRK:
            file_name = "ProjectParametersDISclRK.json"
        else:
            Logger.PrintInfo("DEM",'Benchmark number does not exist')
            sys.exit()

        with open(file_name, 'r') as parameters_file:
            self.DEM_parameters = Parameters(parameters_file.read())

    def __init__(self):
        super(Solution, self).__init__()
        self.nodeplotter = False
        self.LoadParametersFile()
        self.main_path = os.getcwd()

    def GetProblemTypeFilename(self):
        return benchmark

    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):
        return ModelPartIO(modelpart)

    def SetSolverStrategy(self):
        # Strategy object
        element_type = self.DEM_parameters["ElementType"].GetString()
        if (element_type == "SphericPartDEMElement3D" or element_type == "CylinderPartDEMElement2D"):
            import sphere_strategy as SolverStrategy
        elif (element_type == "SphericContPartDEMElement3D" or element_type == "CylinderContPartDEMElement2D"):
            import continuum_sphere_strategy as SolverStrategy
        elif (element_type == "ThermalSphericContPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (element_type == "ThermalSphericPartDEMElement3D"):
            import thermal_sphere_strategy as SolverStrategy
        elif (element_type == "SinteringSphericConPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (element_type == "IceContPartDEMElement3D"):
            import ice_continuum_sphere_strategy as SolverStrategy
        else:
            self.KRATOSprint('Error: Strategy unavailable. Select a different scheme-element')

        return SolverStrategy

    def SetSolver(self):
        return self.solver_strategy.ExplicitStrategy(self.all_model_parts,
                                                     self.creator_destructor,
                                                     self.dem_fem_search,
                                                     self.DEM_parameters,
                                                     self.procedures)

    def SetFinalTime(self):
        self.final_time = final_time

    def Setdt(self):
        self.dt = dt

    def Initialize(self):
        self.DEM_parameters["problem_name"].SetString('benchmark' + str(benchmark_number))
        #self.final_time = slt.final_time
        #self.dt = slt.dt
        #self.graph_print_interval = slt.graph_print_interval
        super(Solution, self).Initialize()

        Logger.PrintInfo("DEM","Computing points in the curve...", 1 + self.number_of_points_in_the_graphic - self.iteration, "point(s) left to finish....",'\n')
        list_of_nodes_ids = [1]
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter = plot_variables.variable_plotter(self.spheres_model_part, list_of_nodes_ids, benchmark_number)
            self.tang_plotter = plot_variables.tangential_force_plotter(self.spheres_model_part, list_of_nodes_ids, self.iteration)

    def ReadModelParts(self):
        super(Solution, self).ReadModelParts()
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
        super(Solution, self).BeforeSolveOperations(time)
        benchmark.ApplyNodalRotation(time, self.dt, self.spheres_model_part)

    def BeforePrintingOperations(self, time):
        super(Solution, self).BeforePrintingOperations(time)
        self.Setdt()
        benchmark.generate_graph_points(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part, time, self.graph_print_interval, self.dt)

    def Finalize(self):
        benchmark.get_final_data(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part)
        super(Solution, self).Finalize()
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter.close_files()
            self.tang_plotter.close_files()

        self.procedures.RemoveFoldersWithResults(self.main_path, self.problem_name)

    def FinalizeTimeStep(self, time):
        super(Solution, self).FinalizeTimeStep(time)
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter.plot_variables(time) #Related to the benchmark in Chung, Ooi
            self.tang_plotter.plot_tangential_force(time)

    def CleanUpOperations(self):
        Logger.PrintInfo("DEM","running CleanUpOperations")
        #DBC.delete_archives() #.......Removing some unuseful files
        super(Solution, self).CleanUpOperations()


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
        del slt
    end = timer.time()
    benchmark.print_results(number_of_points_in_the_graphic, dt, elapsed_time = end - start)
#DBC.delete_archives() #.......Removing some unuseful files
