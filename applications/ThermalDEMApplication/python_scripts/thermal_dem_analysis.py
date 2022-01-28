from   KratosMultiphysics import *
from   KratosMultiphysics.DEMApplication import *
from   KratosMultiphysics.ThermalDEMApplication import *
from   KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from   KratosMultiphysics.ThermalDEMApplication.thermal_dem_io import ThermalDEMIo
from   importlib import import_module

class ThermalDEMAnalysis(DEMAnalysisStage):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

    def _CreateSolver(self):
        def SetSolverStrategy():
            strategy_file_name = self.DEM_parameters["solver_settings"]["strategy"].GetString()
            if (strategy_file_name != "thermal_sphere_strategy"):
                raise Exception('ThermalDEM', 'Strategy \'' + strategy_file_name + '\' is not implemented.')
            else:
                imported_module = import_module("KratosMultiphysics.ThermalDEMApplication" + "." + strategy_file_name)
                return imported_module

        return SetSolverStrategy().ExplicitStrategy(self.all_model_parts,
                                                    self.creator_destructor,
                                                    self.dem_fem_search,
                                                    self.DEM_parameters,
                                                    self.procedures)

    def SetGraphicalOutput(self):
        self.demio = ThermalDEMIo(self.model, self.DEM_parameters, self.post_path, self.all_model_parts)
        if self.DEM_parameters["post_vtk_option"].GetBool():
            import KratosMultiphysics.DEMApplication.dem_vtk_output as dem_vtk_output
            self.vtk_output = dem_vtk_output.VtkOutput(self.main_path, self.problem_name, self.spheres_model_part, self.rigid_face_model_part)

if __name__ == "__main__":
    model = KratosMultiphysics.Model()
    with open("ProjectParametersDEM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    ThermalDEMAnalysis(model, parameters).Run()
