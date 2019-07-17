from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import structural_mechanics_solver

def CreateSolver(model, custom_settings):
    return StructuralMechanicsDirectSensitivitySolver(model, custom_settings)

class StructuralMechanicsDirectSensitivitySolver(structural_mechanics_solver.MechanicalSolver):

    def __init__(self, model, custom_settings):        
        direct_settings = KratosMultiphysics.Parameters("""
        {
            "scheme_settings" : {
                "scheme_type": "direct_structural"
            }
        }
        """)

        self.validate_and_transfer_matching_settings(custom_settings, direct_settings)
        self.scheme_settings = direct_settings["scheme_settings"]

        self.direct_settings = custom_settings["variable_settings"].Clone()            
        custom_settings.RemoveValue("variable_settings")        
        # Construct the base solver.
        super(StructuralMechanicsDirectSensitivitySolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[DirectSensitivitySolver]:: ", "Construction finished")  

    def AddVariables(self):
        super(StructuralMechanicsDirectSensitivitySolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
        # TODO evaluate if these variables should be historical
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD_SENSITIVITY)
        self.print_on_rank_zero("::[DirectSensitivitySolver]:: ", "Variables ADDED")

    def PrepareModelPart(self):
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]!= 3):
            raise Exception("there are currently only 3D direct elements available")
        super(StructuralMechanicsDirectSensitivitySolver, self).PrepareModelPart()

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
                    "ShellThinElement3D3N"            : "AdjointFiniteDifferencingShellThinElement3D3N",                    
                    "CrLinearBeamElement3D2N"         : "AdjointFiniteDifferenceCrBeamElementLinear3D2N",
                    "TrussLinearElement3D2N"          : "AdjointFiniteDifferenceTrussLinearElement3D2N",
                    "TrussElement3D2N"                : "AdjointFiniteDifferenceTrussElement3D2N"
                },
                "condition_name_table" :
                {
                    "PointLoadCondition2D1N"          : "AdjointSemiAnalyticPointLoadCondition2D1N",
                    "PointLoadCondition3D1N"          : "AdjointSemiAnalyticPointLoadCondition3D1N",
                    "PointMomentCondition3D1N"        : "AdjointSemiAnalyticPointMomentCondition3D1N"
                }
            }
        """)

        StructuralMechanicsApplication.ReplaceMultipleElementsAndConditionsProcess(self.main_model_part, replacement_settings).Execute()
        process_info.SetValue(StructuralMechanicsApplication.IS_ADJOINT, True)

        self.print_on_rank_zero("::[DirectSensitivitySolver]:: ", "ModelPart prepared for Solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z, self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z, self.main_model_part)
        self.print_on_rank_zero("::[DirectSensitivitySolver]:: ", "DOF's ADDED.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """ 
        
        # Initialize the design variable 
        if  self.direct_settings["variable_type"].GetString() == "element_data_type":
            self.variable = StructuralMechanicsApplication.DirectSensitivityElementDataVariable(self.main_model_part, self.direct_settings)
        elif  self.direct_settings["variable_type"].GetString() == "nodal_coordinates_type":
            self.variable = StructuralMechanicsApplication.DirectSensitivityNodalCoordinatesVariable(self.main_model_part, self.direct_settings)
        elif  self.direct_settings["variable_type"].GetString() == "point_load_type":
            self.variable = StructuralMechanicsApplication.DirectSensitivityPointLoadVariable(self.main_model_part, self.direct_settings) 
        else:
            raise Exception("invalid variable_type: " + self.direct_settings["variable_type"].GetString())       

        super(StructuralMechanicsDirectSensitivitySolver, self).Initialize()

        self.print_on_rank_zero("::[DirectSensitivitySolver]:: ", "Finished initialization.")

              
    def InitializeSolutionStep(self):
        super(StructuralMechanicsDirectSensitivitySolver, self).InitializeSolutionStep()
        self.variable.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(StructuralMechanicsDirectSensitivitySolver, self).FinalizeSolutionStep()
        self.variable.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        super(StructuralMechanicsDirectSensitivitySolver, self).SolveSolutionStep()         
        

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
        return StructuralMechanicsApplication.DirectStructuralStaticScheme(self.scheme_settings, self.variable)
