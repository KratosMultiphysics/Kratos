from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time monitoring
import time as timer
print(timer.ctime())
initial_time = timer.perf_counter()

import os
import sys

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam

import dam_main
import Dam_benchmarks_class as DamBenchmarksClass   # Related to benchmarks in Chung, Ooi

sys.path.insert(0,'')
start = timer.time()
benchmark_number = int(sys.argv[1])
benchmark = getattr(DamBenchmarksClass, 'Benchmark' + str(benchmark_number))()

list2D   = list(range(201,213))
list3D   = list(range(301,307))

class Solution(dam_main.Solution):

    def __init__(self):
        super(Solution, self).__init__()
        self.main_path = os.getcwd()
        benchmark.set_initial_data()

    def LoadParametersFile(self):
        file_name = None
        if benchmark_number == 201:
            file_name = "ProjectParameters2D01.json"
        elif benchmark_number == 202:
            file_name = "ProjectParameters2D02.json"
        elif benchmark_number == 203:
            file_name = "ProjectParameters2D03.json"
        elif benchmark_number == 204:
            file_name = "ProjectParameters2D04.json"
        elif benchmark_number == 205:
            file_name = "ProjectParameters2D05.json"
        elif benchmark_number == 206:
            file_name = "ProjectParameters2D06.json"
        elif benchmark_number == 207:
            file_name = "ProjectParameters2D07.json"
        elif benchmark_number == 208:
            file_name = "ProjectParameters2D08.json"
        elif benchmark_number == 209:
            file_name = "ProjectParameters2D09.json"
        elif benchmark_number == 210:
            file_name = "ProjectParameters2D10.json"
        elif benchmark_number == 211:
            file_name = "ProjectParameters2D11.json"
        elif benchmark_number == 212:
            file_name = "ProjectParameters2D12.json"
        elif benchmark_number == 301:
            file_name = "ProjectParameters3D01.json"
        elif benchmark_number == 302:
            file_name = "ProjectParameters3D02.json"
        elif benchmark_number == 303:
            file_name = "ProjectParameters3D03.json"
        elif benchmark_number == 304:
            file_name = "ProjectParameters3D04.json"
        elif benchmark_number == 305:
            file_name = "ProjectParameters3D03.json"
        elif benchmark_number == 306:
            file_name = "ProjectParameters3D04.json"
        else:
            print('Benchmark number does not exist')
            sys.exit()

        with open(file_name, 'r') as parameters_file:
            self.ProjectParameters = KratosMultiphysics.Parameters(parameters_file.read())

    def Initialize(self):
        super(Solution, self).Initialize()
        benchmark.set_initial_data()

    def PrintOutput(self):
        super(Solution, self).PrintOutput()
        benchmark.generate_graph_points(self.main_model_part, self.time)

solution = Solution()
solution.Run()
benchmark.print_results()
benchmark.finalize_data()
DamBenchmarksClass.delete_archives() #.......Removing some unuseful files
