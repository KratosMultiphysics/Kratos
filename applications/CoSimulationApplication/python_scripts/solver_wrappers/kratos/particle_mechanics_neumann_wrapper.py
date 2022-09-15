# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing ParticleMechanics
if not CheckIfApplicationsAvailable("ParticleMechanicsApplication"):
    raise ImportError("The ParticleMechanicsApplication is not available!")
import KratosMultiphysics.ParticleMechanicsApplication as KPM
from KratosMultiphysics.ParticleMechanicsApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

def Create(settings, model, solver_name):
    return ParticleMechanicsNeumannWrapper(settings, model, solver_name)

class ParticleMechanicsNeumannWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the ParticleMechanicsApplication of Kratos"""
    """It is designed for the Neumann Interface in the ParticleMechanicsApplication"""

    def _CreateAnalysisStage(self):
        return ParticleMechanicsAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        coupling_model_part = self.model.GetModelPart("MPM_Coupling_Neumann_Interface")
        model_part_name = self.settings["solver_wrapper_settings"]["interface_model_part_name"].GetString()
        model_part = self.model.GetModelPart(model_part_name)
        
        ## Transfer information from coupling_mp to mp
        for coupling_node in coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id

            ## IMPOSED Point load
            point_load = coupling_node.GetSolutionStepValue(KM.CONTACT_FORCE)
            model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KPM.POINT_LOAD, [point_load], model_part.ProcessInfo)
        
        super().SolveSolutionStep()

        ### Save displacement of mpc in coupling node
        for mpc in model_part.Conditions:
            if (mpc.Is(KM.INTERFACE)):
                coupling_id   = mpc.Id

                # Update displacement
                delta_x = mpc.CalculateOnIntegrationPoints(KPM.MPC_DELTA_DISPLACEMENT, model_part.ProcessInfo)[0]
                displacement = coupling_model_part.GetNode(coupling_id).GetSolutionStepValue(KM.DISPLACEMENT)
                displacement += delta_x
                coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KM.DISPLACEMENT,0,displacement)

                # Update velocity
                velocity = mpc.CalculateOnIntegrationPoints(KPM.MPC_VELOCITY, model_part.ProcessInfo)[0]
                coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KM.VELOCITY,0,velocity)
