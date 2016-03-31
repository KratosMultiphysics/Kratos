from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# benchmarking...
import sys
kratos_benchmarking_path = '../../../../../../benchmarking'  # kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking


def FindMinDisp(nodes):
    dx_min = 0.0
    dy_min = 0.0
    dz_min = 0.0
    for node in nodes:
        dx = node.GetSolutionStepValue(DISPLACEMENT_X)
        dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
        dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
        if(dx_min > dx):
            dx_min = dx
        if(dy_min > dy):
            dy_min = dy
        if(dz_min > dz):
            dz_min = dz

    return [dx_min, dy_min, dz_min]


def FindMaxDisp(nodes):
    dx_max = 0.0
    dy_max = 0.0
    dz_max = 0.0
    for node in nodes:
        dx = node.GetSolutionStepValue(DISPLACEMENT_X)
        dy = node.GetSolutionStepValue(DISPLACEMENT_Y)
        dz = node.GetSolutionStepValue(DISPLACEMENT_Z)
        if(dx_max < dx):
            dx_max = dx
        if(dy_max < dy):
            dy_max = dy
        if(dz_max < dz):
            dz_max = dz

    return [dx_max, dy_max, dz_max]


def WriteBenchmarkResults(model_part):

    # write
    abs_tol = 1e-9
    rel_tol = 1e-5

    [dx_min, dy_min, dz_min] = FindMinDisp(model_part.Nodes)
    print([dx_min, dy_min, dz_min])

    if (benchmarking.InBenchmarkingMode()):
        # write
        abs_tol = 1e-9
        rel_tol = 1e-5
        benchmarking.Output(dx_min, "min dx", abs_tol, rel_tol)
        benchmarking.Output(dy_min, "min dy", abs_tol, rel_tol)
        benchmarking.Output(dz_min, "min dz", abs_tol, rel_tol)
 
    [dx_max, dy_max, dz_max] = FindMaxDisp(model_part.Nodes)
    print([dx_max, dy_max, dz_max])

    if (benchmarking.InBenchmarkingMode()):
        benchmarking.Output(dx_max, "max dx", abs_tol, rel_tol)
        benchmarking.Output(dy_max, "max dy", abs_tol, rel_tol)
        benchmarking.Output(dz_max, "max dz", abs_tol, rel_tol)
