import sys

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.GeoMechanicsApplication as KratosGMA

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

# Import base class file
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_solver import GeoMechanicalSolver as GeoSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledThermalPressureSolver(main_model_part, custom_settings)

class CoupledThermalPressureSolver(GeoSolver):

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "ThermalPressureCoupled",
            "model_part_name": "PorousDomain",
            "domain_size" : -1,
            "echo_level": 0,
            "reduction_factor": 0.5,
            "increase_factor": 2.0,
            "min_iterations": 6,
            "max_iterations": 15,
            "number_cycles": 5,
            "solution_type": "quasi_static",
            "reset_displacements": false,
            "start_time": 0.0,
            "rebuild_level": 2,
            "time_stepping": {
                "time_step" : 1.0
            },
            "model_import_settings": {
                    "input_type"    : "mdpa",
                    "input_filename": "unknown_name"
            },
            "pressure_solver_settings": {
                "solver_type"                     : "Pw",
                "model_part_name"                 : "PressureModelPart",
                "domain_size"                     : 2,
                "echo_level"                      : 1,
                "max_iterations"                  : 15,
                "model_import_settings": {
                    "input_type"    : "mdpa",
                    "input_filename": "unknown_name"
                },
                "material_import_settings": {
                    "materials_filename" : "MaterialParameters.json"
                }
            },
            "thermal_solver_settings": {
                "solver_type"                     : "T",
                "model_part_name"                 : "ThermalModelPart",
                "domain_size"                     : 2,
                "echo_level"                      : 1,
                "max_iterations"                  : 15,
                "model_import_settings": {
                    "input_type"   : "use_input_model_part"
                },
                "material_import_settings": {
                    "materials_filename" : "MaterialParameters.json"
                }
            },
            "time_integration_method": "implicit"
        }
        """)

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):

        super(CoupledThermalPressureSolver, self).__init__(model, custom_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        from KratosMultiphysics.GeoMechanicsApplication import geomechanics_solvers_wrapper
        self.pressure_solver = geomechanics_solvers_wrapper.CreateSolverByParameters(self.model, self.settings["pressure_solver_settings"], "OpenMP")
        self.thermal_solver = geomechanics_solvers_wrapper.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"], "OpenMP")

    def AddVariables(self):
        # Import the structural and thermal solver variables. Then merge them to have them in both pressure and thermal solvers.
        self.pressure_solver.AddVariables()
        self.thermal_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.pressure_solver.main_model_part, self.thermal_solver.main_model_part)

    def ImportModelPart(self):
        # Call the structural solver to import the model part from the mdpa
        self.pressure_solver.ImportModelPart()

    def PrepareModelPart(self):
        self.pressure_solver.PrepareModelPart()
        

        # Save the convection diffusion settings
        #convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

        # Here the structural model part is cloned to be thermal model part so that the nodes are shared
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        if self.domain_size == 2:
            modeler.GenerateModelPart(self.pressure_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "TransientThermalElement2D3N",
                                      "TNormalFluxCondition2D2N")
        else:
            modeler.GenerateModelPart(self.pressure_solver.main_model_part,
                                      self.thermal_solver.main_model_part,
                                      "TransientThermalElement3D8N",
                                      "TNormalFluxCondition2D3N")

        # Set the saved convection diffusion settings to the new thermal model part
        #self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

        self.thermal_solver.PrepareModelPart()


    def AddDofs(self):
        self.pressure_solver.AddDofs()
        self.thermal_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.pressure_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def GetMinimumBufferSize(self):
        buffer_size_pressure = self.pressure_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_pressure, buffer_size_thermal)

    def Initialize(self):
        self.pressure_solver.Initialize()
        self.thermal_solver.Initialize()

    def Clear(self):
        (self.pressure_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.pressure_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.pressure_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.pressure_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        self.pressure_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.pressure_solver.Predict()
        self.thermal_solver.Predict()

    def SolveSolutionStep(self):
        pressure_is_converged = self.pressure_solver.SolveSolutionStep()
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()

        return (pressure_is_converged and thermal_is_converged)

    def FinalizeSolutionStep(self):
        self.pressure_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()