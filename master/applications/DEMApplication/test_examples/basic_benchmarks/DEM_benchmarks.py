import time as timer
import os
import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import KratosMultiphysics.DEMApplication.main_script as main_script
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables # Related to benchmarks in Chung, Ooi
import KratosMultiphysics.DEMApplication.DEM_benchmarks_class as DBC   # Related to benchmarks in Chung, Ooi

Logger.PrintInfo("DEM", "WARNING: DEM_benchmarks.py is is deprecated since 20/03/2019")
Logger.PrintInfo("DEM", "WARNING: Please use DEM_benchmarks_analysis.py")

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
listGeneric   = [40]



class Solution(main_script.Solution):

    def GetInputParameters(self):
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
        elif benchmark_number in listGeneric:
            file_name = "ProjectParametersDEMGeneric.json"
        else:
            raise Exception("Benchmark number does not exist")

        with open(file_name, 'r') as parameters_file:
            parameters = Parameters(parameters_file.read())

        return parameters


    def __init__(self, model):
        super().__init__(model)
        self.nodeplotter = False
        self.LoadParametersFile()
        self.main_path = os.getcwd()

    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):
        return ModelPartIO(modelpart)

    def SetSolverStrategy(self):
        # Strategy object
        element_type = self.DEM_parameters["ElementType"].GetString()
        if (element_type == "SphericPartDEMElement3D" or element_type == "CylinderPartDEMElement2D"):
            import KratosMultiphysics.DEMApplication.sphere_strategy as SolverStrategy
        elif (element_type == "SphericContPartDEMElement3D" or element_type == "CylinderContPartDEMElement2D"):
            import KratosMultiphysics.DEMApplication.continuum_sphere_strategy as SolverStrategy
        elif (element_type == "IceContPartDEMElement3D"):
            import KratosMultiphysics.DEMApplication.ice_continuum_sphere_strategy as SolverStrategy
        else:
            self.KratosPrintWarning('Error: Strategy unavailable. Select a different scheme-element')

        return SolverStrategy

    def SetSolver(self):
        return self.solver_strategy.ExplicitStrategy(self.all_model_parts,
                                                     self.creator_destructor,
                                                     self.dem_fem_search,
                                                     self.DEM_parameters,
                                                     self.procedures)

    def SetFinalTime(self):
        self.end_time = end_time

    def SetDt(self):
        self.solver.dt = dt

    def Initialize(self):
        self.DEM_parameters["problem_name"].SetString('benchmark' + str(benchmark_number))
        #self.end_time = slt.end_time
        #self.dt = slt.dt
        #self.graph_print_interval = slt.graph_print_interval
        super().Initialize()

        Logger.PrintInfo("DEM","Computing points in the curve...", 1 + self.number_of_points_in_the_graphic - self.iteration, "point(s) left to finish....",'\n')
        list_of_nodes_ids = [1]
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter = plot_variables.variable_plotter(self.spheres_model_part, list_of_nodes_ids, benchmark_number)
            self.tang_plotter = plot_variables.tangential_force_plotter(self.spheres_model_part, list_of_nodes_ids, self.iteration)

    def ReadModelParts(self):
        super().ReadModelParts()
        benchmark.set_initial_data(self.spheres_model_part, self.rigid_face_model_part, self.iteration, self.number_of_points_in_the_graphic, coeff_of_restitution_iteration)

    def GetMpFilename(self):
        return 'benchmark' + str(benchmark_number) + "DEM"

    def GetInletFilename(self):
        if benchmark_number == 40:
            return 'benchmark' + str(benchmark_number) + "DEM_Inlet"
        else:
            return 'benchmarkDEM_Inlet'

    def GetFemFilename(self):
        return 'benchmark' + str(benchmark_number) + "DEM_FEM_boundary"

    def GetClusterFilename(self):
        return 'benchmark' + str(benchmark_number) + "DEM_Clusters"

    def GetProblemTypeFilename(self):
        return 'benchmark' + str(benchmark_number)

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        benchmark.ApplyNodalRotation(self.time, self.dt, self.spheres_model_part)

    def BeforePrintingOperations(self, time):
        super().BeforePrintingOperations(time)
        self.SetDt()
        benchmark.generate_graph_points(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part, time, self.graph_print_interval, self.dt)

    def Finalize(self):
        benchmark.get_final_data(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part)
        super().Finalize()
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter.close_files()
            self.tang_plotter.close_files()

        self.procedures.RemoveFoldersWithResults(self.main_path, self.problem_name)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter.plot_variables(time) #Related to the benchmark in Chung, Ooi
            self.tang_plotter.plot_tangential_force(time)

    def CleanUpOperations(self):
        Logger.PrintInfo("DEM","running CleanUpOperations")
        super().CleanUpOperations()


end_time, dt, graph_print_interval, number_of_points_in_the_graphic, number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
for coeff_of_restitution_iteration in range(1, number_of_coeffs_of_restitution + 1):
    for iteration in range(1, number_of_points_in_the_graphic + 1):
        model = Model()
        slt = Solution(model)
        slt.iteration = iteration
        slt.dt = dt
        slt.end_time = end_time
        slt.graph_print_interval = graph_print_interval
        slt.number_of_points_in_the_graphic = number_of_points_in_the_graphic
        slt.number_of_coeffs_of_restitution = number_of_coeffs_of_restitution
        slt.Run()
        del slt
    end = timer.time()
    benchmark.print_results(number_of_points_in_the_graphic, dt, elapsed_time = end - start)

