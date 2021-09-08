# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Additional imports
import pdb

# Importing PfemFluidDynamics
if not CheckIfApplicationsAvailable("PfemFluidDynamicsApplication"):
    raise ImportError("The PfemFluidDynamicsApplication is not available!")
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis import PfemFluidDynamicsAnalysis

def Create(settings, model, solver_name):
    return PfemFluidDynamicsWrapper(settings, model, solver_name)

class PfemFluidDynamicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the PfemFluidDynamicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return PfemFluidDynamicsAnalysis(self.model, self.project_parameters)

    def Initialize(self):
        super().Initialize()

        # save list of model parts to be fixed during PFEM solving step
        self.list_of_fix_free_model_parts = self.settings["solver_wrapper_settings"]["fix_free_vel_model_part"].GetStringArray()


    def SolveSolutionStep(self):
        ##### Fix the velocity on the pfem interface nodes only
        pdb.set_trace()
        for fix_model_part in self.list_of_fix_free_model_parts:
            fix_nodes_model_part = KM.PfemFluidDynamicsApplication.FixFreeVelocityOnNodesProcess(self._analysis_stage._GetSolver().model[fix_model_part], 0)
            fix_nodes_model_part.Execute()

        # solve PFEM
        super().SolveSolutionStep()

        ##### Free the velocity on the pfem interface nodes only
        for free_model_part in self.list_of_fix_free_model_parts:
            free_nodes_model_part = KM.PfemFluidDynamicsApplication.FixFreeVelocityOnNodesProcess(self._analysis_stage._GetSolver().model[free_model_part], 1)
            free_nodes_model_part.Execute()
