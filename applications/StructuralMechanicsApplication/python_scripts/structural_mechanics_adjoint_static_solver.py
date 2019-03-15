from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_mechanics_solver

# Import Kratos Utility
import restart_utility

def CreateSolver(model, custom_settings):
    return StructuralMechanicsAdjointStaticSolver(model, custom_settings)

class StructuralMechanicsAdjointStaticSolver(structural_mechanics_solver.MechanicalSolver):

    def __init__(self, model, custom_settings):

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
        elif self.response_function_settings["response_type"].GetString() == "adjoint_non_linear_strain_energy":
            self.response_function = StructuralMechanicsApplication.AdjointNonlinearStrainEnergyResponseFunction(self.main_model_part, self.response_function_settings)
        else:
            raise Exception("invalid response_type: " + self.response_function_settings["response_type"].GetString())

        self.adjoint_postprocess = StructuralMechanicsApplication.AdjointPostprocess(self.main_model_part, self.response_function, self.sensitivity_settings)
        self.adjoint_postprocess.Initialize()

        super(StructuralMechanicsAdjointStaticSolver, self).Initialize()
        self.response_function.Initialize()

        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "Finished initialization.")

    def InitializeSolutionStep(self):
        super(StructuralMechanicsAdjointStaticSolver, self).InitializeSolutionStep()

        # TODO Armin: hdf5 is reading the displacement but not updating the coordinates
        # TODO Mahmoud: check why KratosMultiphysics.VariableUtils().UpdateInitialToCurrentConfiguration update the nodal position here?
        for node in self.main_model_part.Nodes:
            node.X = node.X0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)[0]
            node.Y = node.Y0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)[1]
            node.Z = node.Z0 + node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)[2]

        # Replacing pointers to primal elements with pointers to the serialized elements
        self._load_serialized_model()
        StructuralMechanicsApplication.ReplaceElementsWithSerializedElementsProcess(self.main_model_part, self.loaded_model_part).Execute()
        self.print_on_rank_zero("::[AdjointMechanicalSolver]:: ", "replace primal elements with serialized elements")

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

    #TODO Mahmoud: check if this could be done directly from json file instead of doing it here
    def _load_serialized_model(self):
        loaded_model = KratosMultiphysics.Model()
        model_part_name = self.settings["model_part_name"].GetString()
        self.loaded_model_part = loaded_model.CreateModelPart(model_part_name)
        restart_parameters_load = KratosMultiphysics.Parameters("""
        {
            "input_filename"                 : "test_restart_file",
            "restart_load_file_label"        : "",
            "serializer_trace"               : "no_trace",
            "load_restart_files_from_folder" : false
        }
        """)
        restart_parameters_load["restart_load_file_label"].SetString( str(round(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME], 3) ))
        rest_utility_load = restart_utility.RestartUtility(self.loaded_model_part, restart_parameters_load)
        rest_utility_load.LoadRestart()