# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication import *
import structural_response_function_factory
import structural_mechanics_analysis
import structural_response
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

def ReadYoungsModulusDict(model_part):
    youngs_modulus_dict = {}

    test_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 29, 30, 31, 32, 39, 40, 41, 42, 49, 50, 51, 52, 59, 60, 61, 62, 69, 70, 71, 72, 79, 80, 81, 82, 89, 90, 91, 92, 99, 100, 101, 102, 109, 110, 111, 112, 119, 120, 121, 122, 129, 130, 131, 132, 139, 140, 141, 142, 149, 150, 151, 152, 159, 160, 161, 162, 169, 170, 171, 172, 179, 180, 181, 182, 189, 190, 191, 192, 199, 200]

    for element in model_part.Elements:

        if element.Id in test_list:
            youngs_modulus = 10000.0
        else:
            youngs_modulus = 0.001

        # youngs_modulus = 4600 # constant youngs modulus for all elements

        youngs_modulus_dict[element.Id] = youngs_modulus
    return youngs_modulus_dict

def ModifyYoungsModulus(model_part, youngs_modulus_dict):

    original_properties = model_part.GetProperties()

    for element in model_part.Elements:
        youngs_modulus = youngs_modulus_dict[element.Id]

        old_element_property = element.Properties

        original_property = original_properties[old_element_property.Id]

        new_property = KratosMultiphysics.Properties(old_element_property.Id)

        # TODO use a clone method and only overwrite the youngs modulus
        new_property.SetValue(KratosMultiphysics.DENSITY, original_property.GetValue(KratosMultiphysics.DENSITY) )
        new_property.SetValue(KratosMultiphysics.POISSON_RATIO, original_property.GetValue(KratosMultiphysics.POISSON_RATIO) )
        new_property.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, original_property.GetValue(KratosMultiphysics.CONSTITUTIVE_LAW) )
        new_property.SetValue(KratosMultiphysics.YOUNG_MODULUS, youngs_modulus)

        element.Properties = new_property

        # for output
        element.SetValue(KratosMultiphysics.YOUNG_MODULUS, youngs_modulus)

class SiemensTopologyResponseFunction(structural_response.AdjointResponseFunction):

    def Initialize(self):
        super().Initialize()

        youngs_modulus_dict = ReadYoungsModulusDict(self.primal_model_part)

        ModifyYoungsModulus(self.primal_model_part, youngs_modulus_dict)
        ModifyYoungsModulus(self.adjoint_model_part, youngs_modulus_dict)

        print(self.primal_model_part.GetElement(1).Properties)


with controlledExecutionScope(_get_test_working_dir()):

    # -----------------------------------------------
    # # computation using python structural response
    # -----------------------------------------------
    with open("adjoint_strain_energy_response_parameters.json",'r') as parameter_file:
       parameters = KratosMultiphysics.Parameters( parameter_file.read())

    model = KratosMultiphysics.Model()
    response_function = SiemensTopologyResponseFunction("adjoint_linear_strain_energy", parameters["kratos_response_settings"], model)
    response_function.RunCalculation(True)

    # -----------------------------------------------
    ## computation structural mechanics analysis
    # -----------------------------------------------

    # # primal analysis
    # with open("ProjectParameters.json",'r') as parameter_file:
    #     ProjectParametersPrimal = KratosMultiphysics.Parameters( parameter_file.read())
    # model_primal = KratosMultiphysics.Model()
    # primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
    # primal_analysis.Run()

    # # adjoint analysis
    # with open("AdjointParameters.json",'r') as parameter_file:
    #     ProjectParametersAdjoint = KratosMultiphysics.Parameters( parameter_file.read())
    # model_adjoint = KratosMultiphysics.Model()
    # adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
    # adjoint_analysis.Run()


    # -----------------------------------------------
    ## delete files
    # -----------------------------------------------
    kratos_utilities.DeleteFileIfExisting("cantilever_structure.time")
    for file_name in os.listdir():
        if file_name.endswith(".h5"):
            kratos_utilities.DeleteFileIfExisting(file_name)

    print("\n *** Finished adjoint sensitivity analysis! *** \n")