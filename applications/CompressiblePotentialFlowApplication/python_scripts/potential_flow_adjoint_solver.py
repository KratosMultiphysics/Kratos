from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()
from python_solver import PythonSolver
from potential_flow_solver import LaplacianSolver
# import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def CreateSolver(model, custom_settings):
    return PotentialAdjointSolver(model, custom_settings)

class PotentialAdjointSolver(LaplacianSolver):
    def __init__(self, model, custom_settings):
        adjoint_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_settings" : {
                "scheme_type": "adjoint_structural"
            },
            "element_replace_settings" : {
                "element_name":"IncompressibleAdjointPotentialFlowElement",
                "condition_name": "IncompressibleAdjointPotentialWallCondition"
            }
        }
        """)

        self.validate_and_transfer_matching_settings(custom_settings, adjoint_settings)
        self.scheme_settings = adjoint_settings["scheme_settings"]
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.DISTANCE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CompressiblePotentialFlowApplication.COORDINATES_SENSITIVITY)

        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "Variables ADDED")
        
    def AddDofs(self):
        super(PotentialAdjointSolver, self).AddDofs()
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.CompressiblePotentialFlowApplication.ADJOINT_AUXILIARY_VELOCITY_POTENTIAL, self.main_model_part)
        
    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        if self.response_function_settings["response_type"].GetString() == "adjoint_lift":
            self.response_function = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointLiftResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_lift_coordinates":
            self.response_function = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointLiftCoordinatesResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_lift_coordinates_global":
            self.response_function = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointLiftGlobalCoordinatesResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_potential_coordinates":
            self.response_function = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointPotentialCoordinatesResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_lift_jump_coordinates":
            self.response_function = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointLiftJumpCoordinatesResponseFunction(self.main_model_part, self.response_function_settings)
        else:
            raise Exception("invalid response_type: " + self.response_function_settings["response_type"].GetString())

        self.adjoint_postprocess = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointPostprocess(self.main_model_part, self.response_function, self.sensitivity_settings)
        self.adjoint_postprocess.Initialize()

        scheme = KratosMultiphysics.CompressiblePotentialFlowApplication.AdjointPotentialStaticScheme(self.scheme_settings, self.response_function)
        move_mesh_flag = False #USER SHOULD NOT CHANGE THIS
    
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(
            self.main_model_part, 
            scheme, 
            self.linear_solver,
            builder_and_solver,
            self.settings["compute_reactions"].GetBool(), 
            self.settings["reform_dofs_at_each_step"].GetBool(), 
            self.settings["calculate_solution_norm"].GetBool(), 
            move_mesh_flag)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Check()

        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "Finished initialization.")
    def PrepareModelPart(self):
        super(PotentialAdjointSolver, self).PrepareModelPart()
        # TODO Why does replacement need to happen after reading materials?
        KratosMultiphysics.CompressiblePotentialFlowApplication.ReplaceElementsAndConditionsAdjointProcess(self.main_model_part).Execute()
        self.print_on_rank_zero("::[PotentialAdjointSolver]:: ", "ModelPart prepared for Solver.")

    def InitializeSolutionStep(self):
        super(PotentialAdjointSolver, self).InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(PotentialAdjointSolver, self).FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        super(PotentialAdjointSolver, self).SolveSolutionStep()
        '''after adjoint solution,adjoint_postprocess calculate sensitivities'''
        self.adjoint_postprocess.UpdateSensitivities() # TODO call postprocess here or in FinalizeSolutionStep ?

