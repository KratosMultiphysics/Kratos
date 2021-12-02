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
import numpy as np

# Other imports
import math

def Create(settings, model, solver_name):
    return ParticleMechanicsNeumannWrapper(settings, model, solver_name)

class ParticleMechanicsNeumannWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the ParticleMechanicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return ParticleMechanicsAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        coupling_model_part = self.model.GetModelPart("MPM_Coupling_Neumann_Interface")
        model_part_name = self.project_parameters["coupling_settings"]["interface_model_part_name"].GetString() # TODO this should be specified in "solver_wrapper_settings" in teh cosim-json
        model_part = self.model.GetModelPart(model_part_name)
        
        mpm_model_part = model_part.GetRootModelPart()
        model_part_discretized_line_load = mpm_model_part.GetSubModelPart("inner_discretized_load_particles")
        
        # calculate_mpm=False
        ## Transfer information from coupling_mp to mp
        for coupling_node in coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id

            ## IMPOSED Point load
            point_load = coupling_node.GetSolutionStepValue(KM.CONTACT_FORCE)

            counter = 0
            corresponding_conditions =[]
            for point_cond in model_part_discretized_line_load.Conditions:
                corresponding_condition_id = point_cond.CalculateOnIntegrationPoints(KPM.MPC_CORRESPONDING_CONDITION_ID, model_part.ProcessInfo)[0]
                if (corresponding_condition_id == coupling_id):
                    counter +=1
                    corresponding_conditions.append(point_cond)
            
            point_load = point_load /(counter + 1)
            
            model_part.GetCondition(coupling_id).SetValuesOnIntegrationPoints(KPM.POINT_LOAD, [point_load], model_part.ProcessInfo)

            for point_cond in corresponding_conditions:
                point_cond.SetValuesOnIntegrationPoints(KPM.POINT_LOAD, [point_load], model_part.ProcessInfo)

        
        super().SolveSolutionStep()

        ### Save displacement of mpc in coupling node
        for mpc in model_part.Conditions:
            if (mpc.Is(KM.INTERFACE)):
                coupling_id   = mpc.Id
                delta_x = mpc.CalculateOnIntegrationPoints(KPM.MPC_DISPLACEMENT, model_part.ProcessInfo)[0]

                node = coupling_model_part.GetNode(coupling_id)
                du = (node.X-node.X0) + delta_x[0]
                dw = (node.Y-node.Y0) + delta_x[1]
                dz = (node.Z -node.Z0) + delta_x[2]
                displacement = [du,dw,dz]
                coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KM.DISPLACEMENT,0,displacement)

                velocity = mpc.CalculateOnIntegrationPoints(KPM.MPC_VELOCITY, model_part.ProcessInfo)[0]
                coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KM.VELOCITY,0,velocity)
