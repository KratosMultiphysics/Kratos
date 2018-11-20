from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import structural_mechanics_solver

def CreateSolver(model, custom_settings):
    return StructuralMechanicsAdjointStaticSolver(model, custom_settings)

class StructuralMechanicsAdjointStaticSolver(structural_mechanics_solver.MechanicalSolver):

    def __init__(self, model, custom_settings):

        adjoint_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_settings" : {
                "scheme_type": "adjoint_structural"
            }
        }
        """)

        self.validate_and_transfer_matching_settings(custom_settings, adjoint_settings)
        self.scheme_settings = adjoint_settings["scheme_settings"]

        self.response_function_settings = custom_settings["response_function_settings"].Clone()
        self.sensitivity_settings = custom_settings["sensitivity_settings"].Clone()
        custom_settings.RemoveValue("response_function_settings")
        custom_settings.RemoveValue("sensitivity_settings")
        # Construct the base solver.
        super(StructuralMechanicsAdjointStaticSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "Construction finished")

    def AddVariables(self):
        super(StructuralMechanicsAdjointStaticSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
        # TODO evaluate if these variables should be historical
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "Variables ADDED")

    def PrepareModelPart(self):
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]!= 3):
            raise Exception("there are currently only 3D adjoint elements available")
        super(StructuralMechanicsAdjointStaticSolver, self).PrepareModelPart()
        # TODO Why does replacement need to happen after reading materials?
        StructuralMechanicsApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(self.main_model_part).Execute()
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "ModelPart prepared for Solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z, self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z, self.main_model_part)
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "DOF's ADDED.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        if self.response_function_settings["response_type"].GetString() == "adjoint_local_stress":
            self.response_function = StructuralMechanicsApplication.AdjointLocalStressResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_nodal_displacement":
            self.response_function = StructuralMechanicsApplication.AdjointNodalDisplacementResponseFunction(self.main_model_part, self.response_function_settings)
        elif self.response_function_settings["response_type"].GetString() == "adjoint_linear_strain_energy":
            self.response_function = StructuralMechanicsApplication.AdjointLinearStrainEnergyResponseFunction(self.main_model_part, self.response_function_settings)
        else:
            raise Exception("invalid response_type: " + self.response_function_settings["response_type"].GetString())

        self.adjoint_postprocess = StructuralMechanicsApplication.AdjointPostprocess(self.main_model_part, self.response_function, self.sensitivity_settings)
        self.adjoint_postprocess.Initialize()

        super(StructuralMechanicsAdjointStaticSolver, self).Initialize()

        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "Finished initialization.")

    def Solve(self):
        if self.response_function_settings["response_type"].GetString() == "adjoint_linear_strain_energy":
            self._SolveSolutionStepSpecialLinearStrainEnergy()
        else:
            super(StructuralMechanicsAdjointStaticSolver, self).Solve()


    def InitializeSolutionStep(self):
        super(StructuralMechanicsAdjointStaticSolver, self).InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(StructuralMechanicsAdjointStaticSolver, self).FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.adjoint_postprocess.UpdateSensitivities()

    def SolveSolutionStep(self):
        if self.response_function_settings["response_type"].GetString() == "adjoint_linear_strain_energy":
            self._SolveSolutionStepSpecialLinearStrainEnergy()
        else:
            super(StructuralMechanicsAdjointStaticSolver, self).SolveSolutionStep()
        #after adjoint solution, calculate sensitivities
        #self.response_function.UpdateSensitivities()

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
        self.scheme_settings.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        return StructuralMechanicsApplication.AdjointStructuralStaticScheme(self.scheme_settings, self.response_function)
