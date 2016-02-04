from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import sys
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking


def FindMinDispY(nodes):
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


def WriteBenchmarkResults(model_part):

    [dx_min, dy_min, dz_min] = FindMinDispY(model_part.Nodes)
    print([dx_min, dy_min, dz_min])
    
    if (benchmarking.InBenchmarkingMode()):
        abs_tol = 1e-9
        rel_tol = 1e-5
        benchmarking.Output(dx_min, "min dx", abs_tol, rel_tol)
        benchmarking.Output(dy_min, "min dy", abs_tol, rel_tol)
        benchmarking.Output(dz_min, "min dz", abs_tol, rel_tol)
        
