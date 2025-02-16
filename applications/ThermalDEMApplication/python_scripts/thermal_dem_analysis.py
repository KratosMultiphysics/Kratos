# Imports
from   KratosMultiphysics import *
from   KratosMultiphysics.DEMApplication import *
from   KratosMultiphysics.ThermalDEMApplication import *
from   KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from   KratosMultiphysics.ThermalDEMApplication.thermal_dem_io import ThermalDEMIo
import KratosMultiphysics.ThermalDEMApplication.default_input_settings as DefaultSettings
from   importlib import import_module

# Auxiliary functions
def SetDefaultBoolParameterIfNotExists(default_parameters, key):
    if not key in default_parameters.keys():
        default_parameters.AddEmptyValue(key)
        default_parameters[key].SetBool(False)

# Analysis class
class ThermalDEMAnalysis(DEMAnalysisStage):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

    def GetDefaultInputParameters(self):
        # Get default values of mechanical parameters
        dem_parameters = super().GetDefaultInputParameters()

        # Add default values of thermal parameters
        thermal_settings = DefaultSettings.GetDefaultInputSettings()
        dem_parameters.AddValue("thermal_settings", thermal_settings)

        # Add default values of post options
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostTemperature")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostHeatFlux")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphParticleTempAll")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphParticleTempMin")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphParticleTempMax")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphParticleTempAvg")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphParticleTempAvgVol")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphParticleTempDev")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphMechanicalEnergy")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphDissipatedEnergy")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphThermalEnergy")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphHeatFluxContributions")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphHeatGenerationValues")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostGraphHeatGenerationContributions")
        SetDefaultBoolParameterIfNotExists(dem_parameters, "PostHeatMapGeneration")

        return dem_parameters

    def SetGraphicalOutput(self):
        self.demio = ThermalDEMIo(self.model, self.DEM_parameters, self.post_path, self.all_model_parts)
        if self.DEM_parameters["post_vtk_option"].GetBool():
            import KratosMultiphysics.DEMApplication.dem_vtk_output as dem_vtk_output
            self.vtk_output = dem_vtk_output.VtkOutput(self.main_path, self.problem_name, self.spheres_model_part, self.rigid_face_model_part)
    
    def _CreateSolver(self):
        def SetSolverStrategy():
            strategy_file_name = self.DEM_parameters["solver_settings"]["strategy"].GetString()
            if (strategy_file_name != "thermal_sphere_strategy"):
                raise Exception('ThermalDEM', 'Strategy \'' + strategy_file_name + '\' is not implemented.')
            else:
                imported_module = import_module("KratosMultiphysics.ThermalDEMApplication" + "." + strategy_file_name)
                return imported_module

        # Create cpp strategy for thermal solver
        return SetSolverStrategy().ExplicitStrategy(self.all_model_parts,
                                                    self.creator_destructor,
                                                    self.dem_fem_search,
                                                    self.DEM_parameters,
                                                    self.procedures)

    def GetRVEUtility(self):
        if self.DEM_parameters["Dimension"].GetInt() == 2:
            return RVEWallBoundaryThermal2D(self.rve_evaluation_frequency, self.rve_write_frequency, self.rve_consolidation_stop_criterion, self.rve_consolidation_limit_value, self.rve_inner_volume_offset)
        else:
            raise Exception('Error: The selected RVE utility is not implemented')

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        
        # Update time to print in all model parts (if it has not already been updated in super class)
        if not self.DEM_parameters["ContactMeshOption"].GetBool() and self._GetSolver().write_graph:
            self.UpdateIsTimeToPrintInModelParts(self.IsTimeToPrintPostProcess())

# Initializer
if __name__ == "__main__":
    model = KratosMultiphysics.Model()
    with open("ProjectParametersDEM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    ThermalDEMAnalysis(model, parameters).Run()
