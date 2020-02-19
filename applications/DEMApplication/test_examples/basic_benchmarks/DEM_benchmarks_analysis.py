from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
import KratosMultiphysics.DEMApplication.plot_variables as plot_variables        # Related to benchmarks in Chung, Ooi
import KratosMultiphysics.DEMApplication.DEM_benchmarks_class as DBC

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

def GetInputParameters():
    file_name = None
    if benchmark_number in listDISCONT:
        file_name = "ProjectParametersDISCONT.json"
    elif benchmark_number in listROLLFR:
        file_name = "ProjectParametersROLLFR.json"
    elif benchmark_number in listDEMFEM:
        file_name = "ProjectParametersDEMFEM.json"
    elif benchmark_number in listCONT:
        file_name = "ProjectParametersDEMCONT.json"
    elif benchmark_number == 27:
        file_name = "ProjectParametersUCS.json"
    elif benchmark_number in listDISclZHAO:
        file_name = "ProjectParametersDISclZHAO.json"
    elif benchmark_number in listDISclRK:
        file_name = "ProjectParametersDISclRK.json"
    elif benchmark_number in listGeneric:
            file_name = "ProjectParametersDEMGeneric.json"
    else:
        Logger.PrintInfo("DEM",'Benchmark number does not exist')
        sys.exit()

    with open(file_name, 'r') as parameters_file:
        parameters = Parameters(parameters_file.read())

    return parameters


def partial_delete_archives():
    from glob import glob

    files_to_delete_list = glob('*.time')
    files_to_delete_list.extend(glob('*.gp'))
    files_to_delete_list.extend(glob('*.txt'))
    files_to_delete_list.extend(glob('*.lst'))
    files_to_delete_list.extend(glob('*.info'))
    files_to_delete_list.extend(glob('*.hdf5'))

    for to_erase_file in files_to_delete_list:
        try:
            os.remove(to_erase_file)
        except OSError:
            pass

    # folders_to_delete_list      = glob('*Data')
    # folders_to_delete_list.extend(glob('*ists'))
    # folders_to_delete_list.extend(glob('*ults'))
    # folders_to_delete_list.extend(glob('*he__'))
    # folders_to_delete_list.extend(glob('*aphs'))
    # folders_to_delete_list.extend(glob('*iles'))

    # for to_erase_folder in folders_to_delete_list:
    #     try:
    #         shutil.rmtree(to_erase_folder)
    #     except OSError:
    #         pass


class DEMBenchmarksAnalysisStage(DEMAnalysisStage):

    def __init__(self, model, DEM_parameters):
        self.DEM_parameters = DEM_parameters
        self.project_parameters = DEM_parameters
        default_input_parameters = self.GetDefaultInputParameters()
        self.DEM_parameters.ValidateAndAssignDefaults(default_input_parameters)
        self.FixParametersInconsistencies()
        self.main_path = os.getcwd()
        self.nodeplotter = False
        if benchmark_number in listDISCONT:
            self.nodeplotter = True

        self.DEM_parameters["problem_name"].SetString("benchmark_" + str(benchmark_number))

        super(DEMBenchmarksAnalysisStage, self).__init__(model, DEM_parameters)

    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):
        return ModelPartIO(modelpart)

    def PrintResultsForGid(self, time):
        pass

    def SetDt(self):
        self._GetSolver().dt = dt

    def Initialize(self):
        self.DEM_parameters["problem_name"].SetString('benchmark' + str(benchmark_number))
        #self.end_time = slt.end_time
        #self.dt = slt.dt
        #self.graph_print_interval = slt.graph_print_interval
        super(DEMBenchmarksAnalysisStage, self).Initialize()

        Logger.PrintInfo("DEM","Computing points in the curve...", 1 + self.number_of_points_in_the_graphic - self.iteration, "point(s) left to finish....",'\n')
        list_of_nodes_ids = [1]
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter = plot_variables.variable_plotter(self.spheres_model_part, list_of_nodes_ids, benchmark_number)
            self.tang_plotter = plot_variables.tangential_force_plotter(self.spheres_model_part, list_of_nodes_ids, self.iteration)

    def ReadModelParts(self):
        super(DEMBenchmarksAnalysisStage, self).ReadModelParts()
        benchmark.set_initial_data(self.spheres_model_part, self.rigid_face_model_part, self.iteration, self.number_of_points_in_the_graphic, coeff_of_restitution_iteration)

    def GetInputFilePath(self, file_name):
        if file_name == "DEM_Inlet" and benchmark_number != 40: # TODO: get rid of this exception
            return 'benchmarkDEM_Inlet'
        else:
            return 'benchmark' + str(benchmark_number) + file_name

    def InitializeSolutionStep(self):
        super(DEMBenchmarksAnalysisStage, self).InitializeSolutionStep()
        benchmark.ApplyNodalRotation(self.time, self._GetSolver().dt, self.spheres_model_part)

    def BeforePrintingOperations(self, time):
        super(DEMBenchmarksAnalysisStage, self).BeforePrintingOperations(time)
        self.SetDt()
        benchmark.generate_graph_points(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part, time, self.graph_print_interval, self._GetSolver().dt)

    def Finalize(self):
        benchmark.get_final_data(self.spheres_model_part, self.rigid_face_model_part, self.cluster_model_part)
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter.close_files()
            self.tang_plotter.close_files()

        self.procedures.RemoveFoldersWithResults(self.main_path, self.problem_name)
        super(DEMBenchmarksAnalysisStage, self).Finalize()

    def FinalizeTimeStep(self, time):
        super(DEMBenchmarksAnalysisStage, self).FinalizeTimeStep(time)
        if self.nodeplotter:
            os.chdir(self.main_path)
            self.plotter.plot_variables(time) #Related to the benchmark in Chung, Ooi
            self.tang_plotter.plot_tangential_force(time)

    def CleanUpOperations(self):
        Logger.PrintInfo("DEM","running CleanUpOperations")
        super(DEMBenchmarksAnalysisStage, self).CleanUpOperations()

end_time, dt, graph_print_interval, number_of_points_in_the_graphic, number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
for coeff_of_restitution_iteration in range(1, number_of_coeffs_of_restitution + 1):
    for iteration in range(1, number_of_points_in_the_graphic + 1):
        model = Model()
        parameters = GetInputParameters()
        slt = DEMBenchmarksAnalysisStage(model, parameters)
        slt.iteration = iteration
        slt.dt = dt
        slt.graph_print_interval = graph_print_interval
        slt.number_of_points_in_the_graphic = number_of_points_in_the_graphic
        slt.number_of_coeffs_of_restitution = number_of_coeffs_of_restitution
        slt.DEM_parameters["FinalTime"].SetDouble(end_time)
        slt.project_parameters["problem_data"]["end_time"].SetDouble(end_time)
        slt.DEM_parameters["MaxTimeStep"].SetDouble(dt)
        slt.Run()
        del slt
    end = timer.time()
    benchmark.print_results(number_of_points_in_the_graphic, dt, elapsed_time = end - start)
    #partial_delete_archives()


