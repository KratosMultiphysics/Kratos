from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_mechanics_solver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return AdjointStructuralSolver(main_model_part, custom_settings)

class AdjointStructuralSolver(structural_mechanics_solver.MechanicalSolver):

    def __init__(self, main_model_part, custom_settings):

        adjoint_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_settings" : {
                "scheme_type": "structural"
            }
        }
        """)

        self.validate_and_transfer_matching_settings(custom_settings, adjoint_settings)
        self.scheme_settings = adjoint_settings["scheme_settings"]

        self.response_function_settings = custom_settings["response_function_settings"].Clone()
        custom_settings.RemoveValue("response_function_settings")
        # Construct the base solver.
        super(AdjointStructuralSolver, self).__init__(main_model_part, custom_settings)
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "Construction finished")

    def AddVariables(self):
        super(AdjointStructuralSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
        # TODO evaluate if these variables should be historical
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "Variables ADDED")

    def PrepareModelPartForSolver(self):
        # here we replace the dummy elements we read with proper elements
        self.settings.AddEmptyValue("element_replace_settings")
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                {
                "add_string": "Adjoint",
                "add_before_in_element_name": "Element",
                "add_before_in_condition_name": "Condition",
                "elements_conditions_to_ignore": "ShapeOptimizationCondition",
                "from_primal_to_adjoint": true
                }
                """)

        elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            raise Exception("there is currently no 2D adjoint element")
        else:
            raise Exception("domain size is not 2 or 3")

        StructuralMechanicsApplication.ReplaceElementsAndConditionsForAdjointProblemProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
        super(AdjointStructuralSolver, self).PrepareModelPartForSolver()
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "ModelPart prepared for Solver.")

    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X)
            node.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y)
            node.AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z)
            if self.settings["rotation_dofs"].GetBool():
                node.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X)
                node.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y)
                node.AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z)
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "DOF's ADDED.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.response_function_settings["response_type"].GetString() == "adjoint_local_stress":
            if (domain_size == 2):
                raise Exception("Currently only availible for 3D. Your choice is 2D")
            elif (domain_size == 3):
                self.response_function = StructuralMechanicsApplication.AdjointLocalStressResponseFunction(self.main_model_part, self.response_function_settings)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        elif self.response_function_settings["response_type"].GetString() == "adjoint_nodal_displacement":
            if (domain_size == 2):
                raise Exception("Currently only availible for 3D. Your choice is 2D")
            elif (domain_size == 3):
                self.response_function = StructuralMechanicsApplication.AdjointNodalDisplacementResponseFunction(self.main_model_part, self.response_function_settings)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        elif self.response_function_settings["response_type"].GetString() == "adjoint_strain_energy":
            if (domain_size == 2):
                raise Exception("Currently only availible for 3D. Your choice is 2D")
            elif (domain_size == 3):
                self.response_function = StructuralMechanicsApplication.AdjointStrainEnergyResponseFunction(self.main_model_part, self.response_function_settings)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        else:
            raise Exception("invalid response_type: " + self.response_function_settings["response_type"].GetString())

        super(AdjointStructuralSolver, self).Initialize()

        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "Finished initialization.")

    def Solve(self):
        if self.response_function_settings["response_type"].GetString() == "adjoint_strain_energy":
            self._SolveSpecialStrainEnergy()
        else:
            if self.settings["clear_storage"].GetBool():
                self.Clear()
            mechanical_solution_strategy = self.get_mechanical_solution_strategy()
            mechanical_solution_strategy.Solve()

    def _SolveSpecialStrainEnergy(self):

        self.response_function.Initialize()

        for node in self.main_model_part.Nodes:
            adjoint_displacement = 0.5 * node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            node.SetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT, adjoint_displacement )
            if self.settings["rotation_dofs"].GetBool():
                adjoint_rotation = 0.5 * node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
                node.SetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_ROTATION, adjoint_rotation )

        self.response_function.FinalizeSolutionStep()

    def _create_mechanical_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        builder_and_solver = self.get_builder_and_solver()
        calculate_reaction_flag = False
        reform_dof_set_at_each_step = False
        calculate_norm_dx_flag = False
        move_mesh_flag = False

        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              mechanical_scheme,
                                                              linear_solver,
                                                              builder_and_solver,
                                                              calculate_reaction_flag,
                                                              reform_dof_set_at_each_step,
                                                              calculate_norm_dx_flag,
                                                              move_mesh_flag)

    def _create_solution_scheme(self):
        return StructuralMechanicsApplication.AdjointStructuralScheme(self.scheme_settings, self.response_function)
