# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return StructuralMechanicsAdjointDynamicSolver(model, custom_settings)

class StructuralMechanicsAdjointDynamicSolver(MechanicalSolver):

    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "response_function_settings" : {},
            "sensitivity_settings" : {},
            "time_integration_method" : "implicit",
            "scheme_type"             : "bossak",
            "damp_factor_m"           :-0.3,
            "clear_storage": false
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.AUX_ADJOINT_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_VECTOR_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_VECTOR_3)
        #self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY)
        if self.settings["rotation_dofs"].GetBool():
            #raise RuntimeError("Adding rotation dofs")
            print("adding rotation dofs in adjoint solver")
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ANGULAR_ACCELERATION)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.AUX_ADJOINT_ROTATION_VECTOR_1)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION_VECTOR_2)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION_VECTOR_3)
        # TODO evaluate if these variables should be historical
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.TEMPERATURE_SENSITIVITY)
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Variables ADDED")

    def PrepareModelPart(self):
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]!= 3):
            raise Exception("there are currently only 3D adjoint elements available")
        super().PrepareModelPart()
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
                    "TrussElement3D2N"               : "AdjointFiniteDifferenceTrussElement3D2N",
                    "TotalLagrangianElement2D3N"     : "TotalLagrangianAdjointElement2D3N",
                    "TotalLagrangianElement2D4N"     : "TotalLagrangianAdjointElement2D4N",
                    "TotalLagrangianElement2D6N"     : "TotalLagrangianAdjointElement2D6N",
                    "TotalLagrangianElement3D4N"     : "TotalLagrangianAdjointElement3D4N",
                    "TotalLagrangianElement3D8N"     : "TotalLagrangianAdjointElement3D8N",
                    "SmallDisplacementElement3D4N"   : "AdjointFiniteDifferencingSmallDisplacementElement3D4N",
                    "SmallDisplacementElement3D6N"   : "AdjointFiniteDifferencingSmallDisplacementElement3D6N",
                    "SmallDisplacementElement3D8N"   : "AdjointFiniteDifferencingSmallDisplacementElement3D8N",
                    "SpringDamperElement3D"          : "AdjointFiniteDifferenceSpringDamperElement3D2N",
                    "SpringDamperElement3D2N"        : "AdjointFiniteDifferenceSpringDamperElement3D2N" 
                },
                "condition_name_table" :
                {
                    "PointLoadCondition2D1N"         : "AdjointSemiAnalyticPointLoadCondition2D1N",
                    "PointLoadCondition3D1N"         : "AdjointSemiAnalyticPointLoadCondition3D1N",
                    "SurfaceLoadCondition3D3N"       : "AdjointSemiAnalyticSurfaceLoadCondition3D3N",
                    "SurfaceLoadCondition3D4N"       : "AdjointSemiAnalyticSurfaceLoadCondition3D4N",
                    "SmallDisplacementSurfaceLoadCondition3D3N" : "AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D3N",
                    "SmallDisplacementSurfaceLoadCondition3D4N" : "AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D4N",
                    "LineLoadCondition3D2N"                     : "AdjointSemiAnalyticLineLoadCondition3D2N",
                    "SmallDisplacementLineLoadCondition3D2N"    : "AdjointSemiAnalyticSmallDisplacementLineLoadCondition3D2N"
                },
                "ignore_conditions" : [
                    "SurfaceCondition3D3N",
                    "SurfaceCondition3D4N",
                    "PointCondition3D1N"
                ]
            }
        """) # TODO remove "Condition3D" after issue#4439 is resolved; remove SpringDamperElement3D2N, it is deprecated

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
        #raise Exception("invalid response_type: " + response_type)
        self.sensitivity_builder = KratosMultiphysics.SensitivityBuilder(self.settings["sensitivity_settings"], self.main_model_part, self.response_function)
        self.sensitivity_builder.Initialize()

        super().Initialize()
        self.response_function.Initialize()

        # set the strategy to only build LHS once, and build RHS for all other solves
        self._GetSolutionStrategy().SetRebuildLevel(1)

        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Finished initialization.")

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        
        #print("intialized mech solver")
        self.response_function.InitializeSolutionStep()
        #print("initialized adj dyn solver")

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.sensitivity_builder.UpdateSensitivities()

    def SolveSolutionStep(self):
        # if self.settings["response_function_settings"]["response_type"].GetString() == "adjoint_linear_strain_energy":
        #     self._SolveSolutionStepSpecialLinearStrainEnergy()
        #     return True
        # else:
        #     return super().SolveSolutionStep()
        return super().SolveSolutionStep()

    # def _SolveSolutionStepSpecialLinearStrainEnergy(self):
    #     for node in self.main_model_part.Nodes:
    #         adjoint_displacement = 0.5 * node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
    #         node.SetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT, adjoint_displacement )
    #         if self.settings["rotation_dofs"].GetBool():
    #             adjoint_rotation = 0.5 * node.GetSolutionStepValue(KratosMultiphysics.ROTATION)
    #             node.SetSolutionStepValue(StructuralMechanicsApplication.ADJOINT_ROTATION, adjoint_rotation )

    def _CreateSolutionStrategy(self):
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

    def _CreateScheme(self):
        settings = KratosMultiphysics.Parameters("""{
            "name"         : "adjoint_bossak",
            "scheme_type"  : "bossak",
            "alpha_bossak" : -0.3
        }""")
        return KratosMultiphysics.ResidualBasedAdjointBossakScheme(settings, self.response_function)
    
    def _ComputeDeltaTime(self):
        """This function returns the delta time
        Note that the analysis does not support the automatic time stepping
        Also note that as it requires to go backwards in time, here we return minus the user time step
        """
        # if self.settings["time_stepping"]["automatic_time_step"].GetBool():
        #     raise Exception("Automatic time stepping is not supported by adjoint fluid solver.")

        delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        if delta_time > 0.0:
            # Expected positive time_step is reversed to advance backwards in time
            delta_time *= -1.0
        else:
            # TODO: Remove this check after the backwards compatibility period
            # In order to keep backwards compatibility, we throw a warning if a negative time_step is provided
            KratosMultiphysics.Logger.PrintWarning("Setting a negative 'time_step' is no longer needed by 'StructuralMechanicsAdjointDynamicSolver'.")

        return delta_time
    
    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt
        new_step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] - 1

        for mp_name in self.model.GetModelPartNames():
            mp : KratosMultiphysics.ModelPart = self.model.GetModelPart(mp_name)
            mp.ProcessInfo[KratosMultiphysics.STEP] = new_step

        self.main_model_part.CloneTimeStep(new_time)
        return new_time
