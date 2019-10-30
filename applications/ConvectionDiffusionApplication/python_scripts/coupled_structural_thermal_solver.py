from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA
import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledThermoMechanicalSolver(main_model_part, custom_settings)

class CoupledThermoMechanicalSolver(PythonSolver):

    def __init__(self, model, custom_settings):

        super(CoupledThermoMechanicalSolver, self).__init__(model, custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "ThermoMechanicallyCoupled",
            "domain_size" : -1,
            "echo_level": 0,
            "structural_solver_settings": {
                "solver_type": "Static",
                "model_part_name"                 : "Structure",
                "domain_size"                     : 2,
                "echo_level"                      : 1,
                "analysis_type"                   : "non_linear",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                },
            "material_import_settings"        : {
                "materials_filename" : "StructuralMaterials.json"
            }
            },
            "thermal_solver_settings": {
                "solver_type": "Transient",
                "analysis_type": "linear",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                        "materials_filename": "ThermalMaterials.json"
                }
            }
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural
        self.structural_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structural_solver_settings"],"OpenMP")

        from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")

    def AddVariables(self):
        # Import the structural and thermal solver variables. Then merge them to have them in both structural and thermal solvers.
        self.structural_solver.AddVariables()
        self.thermal_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.structural_solver.main_model_part, self.thermal_solver.main_model_part)

    def ImportModelPart(self):
        # Call the structural solver to import the model part from the mdpa
        self.structural_solver.ImportModelPart()

        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

        # Here the structural model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff2D",
                                      "ThermalFace2D2N")
        else:
            modeler.GenerateModelPart(self.structural_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "EulerianConvDiff3D",
                                      "ThermalFace3D3N")

        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

    def PrepareModelPart(self):
        self.structural_solver.PrepareModelPart()
        self.thermal_solver.PrepareModelPart()

    def AddDofs(self):
        self.structural_solver.AddDofs()
        self.thermal_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.structural_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.structural_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.structural_solver.Initialize()
        self.thermal_solver.Initialize()

    def Clear(self):
        (self.structural_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.structural_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.structural_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.structural_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        self.structural_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.structural_solver.Predict()
        self.thermal_solver.Predict()

    def SolveSolutionStep(self):
        fluid_is_converged = self.structural_solver.SolveSolutionStep()
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()

        return (fluid_is_converged and thermal_is_converged)

    def FinalizeSolutionStep(self):
        self.structural_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

