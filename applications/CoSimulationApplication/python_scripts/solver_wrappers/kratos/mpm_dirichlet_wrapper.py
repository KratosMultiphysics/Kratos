# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing MPM
if not CheckIfApplicationsAvailable("MPMApplication"):
    raise ImportError("The MPMApplication is not available!")
import KratosMultiphysics.MPMApplication as KPM
from KratosMultiphysics.MPMApplication.mpm_analysis import MPMAnalysis

# Other imports
import math

def Create(settings, model, solver_name):
    return MPMDirichletWrapper(settings, model, solver_name)

class MPMDirichletWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the MPMApplication of Kratos."""
    """It is designed for the Dirichlet Interface in the MPMApplication"""

    def _CreateAnalysisStage(self):
        return MPMAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        coupling_model_part = self.model.GetModelPart("MPM_Coupling_Dirichlet_Interface")
        model_part_name = self.settings["solver_wrapper_settings"]["interface_model_part_name"].GetString()
        model_part = self.model.GetModelPart(model_part_name)

        ## Transfer information from coupling_mp to mp
        for coupling_node in coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id

            ## IMPOSED DISPLACEMENT
            total_displacement = coupling_node.GetSolutionStepValue(KM.DISPLACEMENT,0)
            old_displacement = model_part.GetCondition(coupling_id).CalculateOnIntegrationPoints(KPM.MPC_DISPLACEMENT, model_part.ProcessInfo)[0]
            incremental_displacement = total_displacement - old_displacement
            model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KPM.MPC_IMPOSED_DISPLACEMENT, [incremental_displacement], model_part.ProcessInfo)

            ## ADD VELOCITY
            current_velocity = coupling_node.GetSolutionStepValue(KM.VELOCITY,0)
            model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KPM.MPC_VELOCITY, [current_velocity], model_part.ProcessInfo)

            ## ADD NORMAL
            normal = coupling_node.GetSolutionStepValue(KM.NORMAL,0)
            # Check and see whether the normal is not zero
            norm_normal = math.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
            if norm_normal > 1.e-10:
                model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KPM.MPC_NORMAL, [normal], model_part.ProcessInfo)

        super().SolveSolutionStep()

        ### Get contact force from mp to coupling_mp
        for mpc in model_part.Conditions:
            if (mpc.Is(KM.INTERFACE)):
                coupling_id   = mpc.Id
                contact_force = mpc.CalculateOnIntegrationPoints(KPM.MPC_CONTACT_FORCE, model_part.ProcessInfo)[0]
                coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KM.CONTACT_FORCE,0,contact_force)
