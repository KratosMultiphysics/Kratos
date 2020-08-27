from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_adjoint_dynamic_analysis import StructuralMechanicsAdjointDynamicAnalysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

def _get_test_working_dir():
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(this_file_dir)

with controlledExecutionScope(_get_test_working_dir()):

    # primal analysis
    with open("ProjectParameters.json",'r') as parameter_file:
        primal_parameters = KratosMultiphysics.Parameters(parameter_file.read())
    primal_model = KratosMultiphysics.Model()
    primal_simulation = StructuralMechanicsAnalysis(primal_model,primal_parameters)
    primal_simulation.Run()
    print("PRIMAL PROBLEM SOLVED")

    with open("AdjointParameters.json",'r') as parameter_file:
       adjoint_parameters = KratosMultiphysics.Parameters(parameter_file.read())
    adjoint_model = KratosMultiphysics.Model()
    adjoint_simulation = StructuralMechanicsAdjointDynamicAnalysis(adjoint_model,adjoint_parameters)
    adjoint_simulation.Run()
 
    model_part_name = primal_parameters["solver_settings"]["model_part_name"].GetString()
    adjoint_model_part = adjoint_simulation.model.GetModelPart(model_part_name)
    sensitivity = adjoint_model_part.Elements[1].GetValue(KratosMultiphysics.StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS_SENSITIVITY)
    print("\n sensitivity:", sensitivity[0])     
    print("\n ** finished computations ** \n")
 