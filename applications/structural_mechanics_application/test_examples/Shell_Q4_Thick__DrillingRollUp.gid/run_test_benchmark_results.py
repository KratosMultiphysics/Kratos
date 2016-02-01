from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# benchmarking...
import sys
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking

def WriteBenchmarkResults(model_part):
    
    if (benchmarking.InBenchmarkingMode()):
        # find max Y rotation
        r_max = 0.0
        for node in model_part.Nodes:
            ir = node.GetSolutionStepValue(ROTATION_Z)
            if(ir > r_max):
                r_max = ir
        # write
        abs_tol = 1e-9
        rel_tol = 1e-5
        benchmarking.Output(r_max, "Z Rotation", abs_tol, rel_tol)
