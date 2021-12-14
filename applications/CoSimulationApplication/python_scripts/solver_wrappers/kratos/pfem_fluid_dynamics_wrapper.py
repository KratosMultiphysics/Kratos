# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

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
        # Fix the velocity on the pfem interface nodes only
        for fix_model_part_name in self.list_of_fix_free_model_parts:
            fix_model_part = self._analysis_stage._GetSolver().model[fix_model_part_name]
            fix_nodes_model_part = KM.PfemFluidDynamicsApplication.PFEMFixFreeVelocityOnNodesProcess(fix_model_part, True)
            fix_nodes_model_part.Execute()
            self._PrintFixationStatus(fix_model_part_name, True)

        # solve PFEM
        super().SolveSolutionStep()

        # Free the velocity on the pfem interface nodes only
        for free_model_part_name in self.list_of_fix_free_model_parts:
            free_model_part = self._analysis_stage._GetSolver().model[free_model_part_name]
            free_nodes_model_part = KM.PfemFluidDynamicsApplication.PFEMFixFreeVelocityOnNodesProcess(free_model_part, False)
            free_nodes_model_part.Execute()
            self._PrintFixationStatus(free_model_part_name, False)

    def _PrintFixationStatus(self, model_part_name, fix_free_status):
        if self.echo_level > 0:
            if fix_free_status:
                result_msg = " FIXED VELOCITY VALUES FOR THE NODES IN MODEL PART: " + model_part_name
            else:
                result_msg = " FREED VELOCITY VALUES FOR THE NODES IN MODEL PART: " + model_part_name

            cs_tools.cs_print_info(self._ClassName(), result_msg)
