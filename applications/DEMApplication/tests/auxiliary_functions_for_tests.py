import os
import KratosMultiphysics
from KratosMultiphysics import Logger
import multiprocessing


def RunStageInSelectedNumberOfOpenMPThreads(my_obj, number_of_threads):
    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        KratosMultiphysics.ParallelUtilities.SetNumThreads(number_of_threads)
    else:
        Logger.PrintInfo("cores detected:", multiprocessing.cpu_count())
        KratosMultiphysics.ParallelUtilities.SetNumThreads(number_of_threads)

    my_obj.Run()

    if "OMP_NUM_THREADS" in os.environ:
        KratosMultiphysics.ParallelUtilities.SetNumThreads(int(initial_number_of_threads))
    else:
        KratosMultiphysics.ParallelUtilities.SetNumThreads(multiprocessing.cpu_count())

def CreateAndRunStageInSelectedNumberOfOpenMPThreads(my_class, model, parameters, number_of_threads):

    if "OMP_NUM_THREADS" in os.environ:
        initial_number_of_threads = os.environ['OMP_NUM_THREADS']
        KratosMultiphysics.ParallelUtilities.SetNumThreads(number_of_threads)
    else:
        Logger.PrintInfo("cores detected:", multiprocessing.cpu_count())
        KratosMultiphysics.ParallelUtilities.SetNumThreads(number_of_threads)

    if isinstance(parameters, str):
        with open(parameters,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        my_class(model, project_parameters).Run()
    else:
        my_class(model, parameters).Run()

    if "OMP_NUM_THREADS" in os.environ:
        KratosMultiphysics.ParallelUtilities.SetNumThreads(int(initial_number_of_threads))
    else:
        KratosMultiphysics.ParallelUtilities.SetNumThreads(multiprocessing.cpu_count())


class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

def  GetHardcodedNumberOfThreads():
    return 1
