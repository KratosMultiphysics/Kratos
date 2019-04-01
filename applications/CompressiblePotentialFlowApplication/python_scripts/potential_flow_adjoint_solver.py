from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_solver import LaplacianSolver
import KratosMultiphysics.StructuralMechanicsApplication

def CreateSolver(model, custom_settings):
    return PotentialAdjointSolver(model, custom_settings)

class PotentialAdjointSolver(LaplacianSolver):
    def __init__(self, model, custom_settings):
        adjoint_settings = KratosMultiphysics.Parameters("""
        {
            "element_replace_settings" : {
                "element_name":"IncompressibleAdjointPotentialFlowElement",
                "condition_name": "IncompressibleAdjointPotentialWallCondition"
            }
        }
        """)

        self.validate_and_transfer_matching_settings(custom_settings, adjoint_settings)
        self.element_replace_settings=adjoint_settings["element_replace_settings"]

        self.response_function_settings = custom_settings["response_function_settings"].Clone()
        self.sensitivity_settings = custom_settings["sensitivity_settings"].Clone()
        custom_settings.RemoveValue("response_function_settings")
        custom_settings.RemoveValue("sensitivity_settings")
        # Construct the base solver.
        super(PotentialAdjointSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "Construction finished")

    def AddVariables(self):
        super(PotentialAdjointSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)

        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL, self.main_model_part)

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        if self.response_function_settings["response_type"].GetString() == "adjoint_lift_jump_coordinates":
            self.response_function = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointLiftJumpCoordinatesResponseFunction(self.main_model_part, self.response_function_settings)
        else:
            raise Exception("invalid response_type: " + self.response_function_settings["response_type"].GetString())

        self.adjoint_postprocess=KratosMultiphysics.SensitivityBuilder(self.sensitivity_settings,self.main_model_part, self.response_function)
        self.adjoint_postprocess.Initialize()

        scheme = KratosMultiphysics.ResidualBasedAdjointStaticScheme(self.response_function)
        move_mesh_flag = False #USER SHOULD NOT CHANGE THIS

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.incompressible_solution_stratety = KratosMultiphysics.ResidualBasedLinearStrategy(
            self.main_model_part,
            scheme,
            self.linear_solver,
            builder_and_solver,
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["calculate_solution_norm"].GetBool(),
            move_mesh_flag)

        (self.incompressible_solution_stratety).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.incompressible_solution_stratety.Check()

        self.response_function.Initialize()

        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "Finished initialization.")
    def PrepareModelPart(self):
        super(PotentialAdjointSolver, self).PrepareModelPart()
       # defines how the primal elements should be replaced with their adjoint counterparts
        replacement_settings = KratosMultiphysics.Parameters("""
            {
                "element_name_table" :
                {
                    "IncompressiblePotentialFlowElement2D3N" : "AdjointIncompressiblePotentialFlowElement2D3N"
                },
                "condition_name_table" :
                {
                    "PotentialWallCondition2D2N"             : "AdjointIncompressiblePotentialWallCondition2D2N"
                }
            }
        """)

        KratosMultiphysics.StructuralMechanicsApplication.ReplaceMultipleElementsAndConditionsProcess(self.main_model_part, replacement_settings).Execute()
        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "ModelPart prepared for Solver.")

    def InitializeSolutionStep(self):
        super(PotentialAdjointSolver, self).InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(PotentialAdjointSolver, self).FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.adjoint_postprocess.UpdateSensitivities()

    def SolveSolutionStep(self):
        super(PotentialAdjointSolver, self).SolveSolutionStep()

