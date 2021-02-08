import os
import KratosMultiphysics
from KratosMultiphysics import Logger
import multiprocessing

def CreateAndRunStageInSelectedNumberOfOpenMPThreads(my_obj, model, parameters_file_name, number_of_threads):
    omp_utils = KratosMultiphysics.OpenMPUtils()
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        omp_utils.SetNumThreads(number_of_threads)
    else:
        Logger.PrintInfo("cores detected:", multiprocessing.cpu_count())
        omp_utils.SetNumThreads(number_of_threads)

    with open(parameters_file_name,'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    my_obj(model, project_parameters).Run()

    if "OMP_NUM_THREADS" in os.environ:
        omp_utils.SetNumThreads(int(initial_number_of_threads))
    else:
        omp_utils.SetNumThreads(multiprocessing.cpu_count())
