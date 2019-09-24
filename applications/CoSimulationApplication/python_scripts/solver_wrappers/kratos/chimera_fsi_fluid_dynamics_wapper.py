from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Other imports
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.ChimeraApplication.fluid_chimera_analysis import FluidChimeraAnalysis

def Create(cosim_solver_settings, level):
    return ChimeraFsiFluidDynamicsWrapper(cosim_solver_settings, level)

class ChimeraFsiFluidDynamicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    def _CreateAnalysisStage(self):
        return FluidChimeraAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        self._analysis_stage.FinalizeSolutionStep()
        self._analysis_stage.InitializeSolutionStep()
        super(ChimeraFsiFluidDynamicsWrapper, self).SolveSolutionStep()

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__
