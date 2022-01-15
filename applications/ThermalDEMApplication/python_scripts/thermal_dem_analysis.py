from   importlib import import_module
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
from   KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

class ThermalDEMAnalysis(DEMAnalysisStage):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

    def _CreateSolver(self):
        def SetSolverStrategy():
            strategy_file_name = self.DEM_parameters["solver_settings"]["strategy"].GetString()
            imported_module = import_module("KratosMultiphysics.ThermalDEMApplication" + "." + strategy_file_name)
            return imported_module

        return SetSolverStrategy().ExplicitStrategy(self.all_model_parts,
                                                    self.creator_destructor,
                                                    self.dem_fem_search,
                                                    self.DEM_parameters,
                                                    self.procedures)
