# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# --- STD Imports ---
from typing import Union


def CreateSolver(
    model: KratosMultiphysics.Model,
    parameters: KratosMultiphysics.Parameters) -> "StructuralMechanicsAdjointSolver":
        return StructuralMechanicsAdjointSolver(
            model,
            parameters)


class StructuralMechanicsAdjointSolver(MechanicalSolver):

    def __init__(
        self,
        model: KratosMultiphysics.Model,
        parameters: KratosMultiphysics.Parameters):
            super().__init__(model, parameters)


    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        additional_parameters: KratosMultiphysics.Parameters = KratosMultiphysics.Parameters("""{
            "response_function_settings" : {},
            "sensitivity_settings" : {},
            "sensitivity_variables" : ["SHAPE_SENSITIVITY", "TEMPERATURE_SENSITIVITY"]
        }""")
        additional_parameters.AddMissingParameters(super().GetDefaultParameters())
        return additional_parameters


    def AddVariables(self) -> None:
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_DISPLACEMENT)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_ROTATION)

        kernel: KratosMultiphysics.Kernel = KratosMultiphysics.Kernel()
        for sensitivity_variable_name in self.settings["sensitivity_variables"].GetStringArray():
            variable: Union[None,KratosMultiphysics.DoubleVariable,KratosMultiphysics.Array1DVariable3] = None
            if kernel.HasDoubleVariable(sensitivity_variable_name):
                variable = kernel.GetDoubleVariable(sensitivity_variable_name)
            elif kernel.HasArrayVariable(sensitivity_variable_name):
                variable = kernel.GetArrayVariable(sensitivity_variable_name)
            if variable is None:
                raise RuntimeError(f"\"{sensitivity_variable_name}\" is not a registered variable")


    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_DISPLACEMENT_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_DISPLACEMENT_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_DISPLACEMENT_Z, self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_ROTATION_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_ROTATION_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_ROTATION_Z, self.main_model_part)


    def Initialize(self):
        response_type = self.settings["response_function_settings"]["response_type"].GetString()
        if response_type == "adjoint_local_stress":
            self.response_function = StructuralMechanicsApplication.AdjointLocalStressResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        elif response_type == "adjoint_max_stress":
            self.response_function = StructuralMechanicsApplication.AdjointMaxStressResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        elif response_type == "adjoint_nodal_displacement":
            self.response_function = StructuralMechanicsApplication.AdjointNodalDisplacementResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        elif response_type == "adjoint_linear_strain_energy":
            self.response_function = StructuralMechanicsApplication.AdjointLinearStrainEnergyResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        elif response_type == "adjoint_nodal_reaction":
            self.response_function = StructuralMechanicsApplication.AdjointNodalReactionResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        else:
            raise Exception("invalid response_type: " + response_type)

        self.sensitivity_builder: KratosMultiphysics.SensitivityBuilder = KratosMultiphysics.SensitivityBuilder(
            self.settings["sensitivity_settings"],
            self.main_model_part,
            self.response_function)
        self.sensitivity_builder.Initialize()

        super().Initialize()
        self.response_function.Initialize()


    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()


    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.sensitivity_builder.UpdateSensitivities()


    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            if self.settings["compute_reactions"].GetBool():
                raise Exception("\"compute_reactions\" is not allowed for adjoint models parts")
            if self.settings["move_mesh_flag"].GetBool():
                raise Exception("\"move_mesh_flag\" is not allowed for adjoint models parts")
            mechanical_solution_strategy = self._create_linear_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available for adjoints!\n"
            err_msg += "Available options are: \"linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy


    def _CreateScheme(self):
        return KratosMultiphysics.StaticAdjointScheme(self.response_function)
