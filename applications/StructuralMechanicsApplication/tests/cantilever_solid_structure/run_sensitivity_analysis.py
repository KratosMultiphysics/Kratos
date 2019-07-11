# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_response_function_factory
import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

if kratos_utilities.CheckIfApplicationsAvailable("EigenSolversApplication"):
    has_eigensolvers_application = True
else:
    has_eigensolvers_application = False


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

    # -----------------------------------------------
    ## computation using python structural response
    # -----------------------------------------------
    #with open("adjoint_strain_energy_response_parameters.json",'r') as parameter_file:
    #    parameters = KratosMultiphysics.Parameters( parameter_file.read())
    #
    #model = KratosMultiphysics.Model()
    #response_function = structural_response_function_factory.CreateResponseFunction("strain_energy", parameters["kratos_response_settings"], model)
    #response_function.RunCalculation(True)

    # -----------------------------------------------
    ## computation structural mechanics analysis
    # -----------------------------------------------

    # primal analysis
    with open("ProjectParameters.json",'r') as parameter_file:
        ProjectParametersPrimal = KratosMultiphysics.Parameters( parameter_file.read())
    model_primal = KratosMultiphysics.Model()
    primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
    primal_analysis.Run()

    # adjoint analysis
    with open("AdjointParameters.json",'r') as parameter_file:
        ProjectParametersAdjoint = KratosMultiphysics.Parameters( parameter_file.read())
    model_adjoint = KratosMultiphysics.Model()
    adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
    adjoint_analysis.Run()


    # -----------------------------------------------
    ## delete files
    # -----------------------------------------------
    kratos_utilities.DeleteFileIfExisting("cantilever_structure.time")
    for file_name in os.listdir():
        if file_name.endswith(".h5"):
            kratos_utilities.DeleteFileIfExisting(file_name)

    print("\n *** Finished adjoint sensitivity analysis! *** \n")