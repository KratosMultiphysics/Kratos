from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return StructuralMechanicsAdjointStaticSolver(model, custom_settings)

class StructuralMechanicsAdjointStaticSolver(MechanicalSolver):

    def __init__(self, model, custom_settings):

        self.response_function_settings = custom_settings["response_function_settings"].Clone()
        self.sensitivity_settings = custom_settings["sensitivity_settings"].Clone()
        custom_settings.RemoveValue("response_function_settings")
        custom_settings.RemoveValue("sensitivity_settings")
        # Construct the base solver.
        super(StructuralMechanicsAdjointStaticSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Construction finished")

    def AddVariables(self):
        super(StructuralMechanicsAdjointStaticSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
        # TODO evaluate if these variables should be historical
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Variables ADDED")

    def PrepareModelPart(self):
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]!= 3):
            raise Exception("there are currently only 3D adjoint elements available")
        super(StructuralMechanicsAdjointStaticSolver, self).PrepareModelPart()
        # TODO Why does replacement need to happen after reading materials?

        process_info = self.main_model_part.ProcessInfo
        if (process_info.Has(StructuralMechanicsApplication.IS_ADJOINT) and
            process_info.GetValue(StructuralMechanicsApplication.IS_ADJOINT)):
            raise RuntimeError("Modelpart '{}' is already adjoint modelpart!".format(self.main_model_part.Name))

        # defines how the primal elements should be replaced with their adjoint counterparts
        replacement_settings = KratosMultiphysics.Parameters("""
            {
                "element_name_table" :
                {
                    "ShellThinElement3D3N"           : "AdjointFiniteDifferencingShellThinElement3D3N",
                    "CrLinearBeamElement3D2N"        : "AdjointFiniteDifferenceCrBeamElementLinear3D2N",
                    "TrussLinearElement3D2N"         : "AdjointFiniteDifferenceTrussLinearElement3D2N",
                    "TrussElement3D2N"               : "AdjointFiniteDifferenceTrussElement3D2N"
                },
                "condition_name_table" :
                {
                    "PointLoadCondition2D1N"         : "AdjointSemiAnalyticPointLoadCondition2D1N",
                    "PointLoadCondition3D1N"         : "AdjointSemiAnalyticPointLoadCondition3D1N"
                }
            }
        """)

        StructuralMechanicsApplication.ReplaceMultipleElementsAndConditionsProcess(self.main_model_part, replacement_settings).Execute()
        process_info.SetValue(StructuralMechanicsApplication.IS_ADJOINT, True)

        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "ModelPart prepared for Solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z, self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "DOF's ADDED.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        if self.response_function_settings["response_type"].GetString() == "adjoint_local_stress":
            self.response_function = StructuralMechanicsApplication.AdjointLocalStressResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_nodal_displacement":
            self.response_function = StructuralMechanicsApplication.AdjointNodalDisplacementResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_linear_strain_energy":
            self.response_function = StructuralMechanicsApplication.AdjointLinearStrainEnergyResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_nodal_reaction":
            self.response_function = StructuralMechanicsApplication.AdjointNodalReactionResponseFunction(self.main_model_part, self.response_function_settings)
        else:
            raise Exception("invalid response_type: " + self.response_function_settings["response_type"].GetString())

        self.sensitivity_builder = KratosMultiphysics.SensitivityBuilder(self.sensitivity_settings, self.main_model_part, self.response_function)
        self.sensitivity_builder.Initialize()

        super(StructuralMechanicsAdjointStaticSolver, self).Initialize()
        self.response_function.Initialize()

        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Finished initialization.")

    def InitializeSolutionStep(self):
        super(StructuralMechanicsAdjointStaticSolver, self).InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(StructuralMechanicsAdjointStaticSolver, self).FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.sensitivity_builder.UpdateSensitivities()

    def SolveSolutionStep(self):
        if self.response_function_settings["response_type"].GetString() == "adjoint_linear_strain_energy":
            self._SolveSolutionStepSpecialLinearStrainEnergy()
        else:
            super(StructuralMechanicsAdjointStaticSolver, self).SolveSolutionStep()

    def _SolveSolutionStepSpecialLinearStrainEnergy(self):
        for node in self.main_model_part.Nodes:
            adjoint_displacement = 0.5 * node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            node.SetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT, adjoint_displacement )
            if self.settings["rotation_dofs"].GetBool():
                adjoint_rotation = 0.5 * node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
                node.SetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_ROTATION, adjoint_rotation )

    def _create_mechanical_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            if self.settings["compute_reactions"].GetBool():
                raise Exception("\"compute_reactions\" is not possible for adjoint models parts")
            if self.settings["move_mesh_flag"].GetBool():
                raise Exception("\"move_mesh_flag\" is not allowed for adjoint models parts")
            mechanical_solution_strategy = self._create_linear_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available for adjoints!\n"
            err_msg += "Available options are: \"linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedAdjointStaticScheme(self.response_function)
