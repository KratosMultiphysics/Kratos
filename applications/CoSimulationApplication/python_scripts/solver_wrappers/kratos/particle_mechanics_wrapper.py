from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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

# Other imports
import math

def Create(settings, solver_name):
    return ParticleMechanicsWrapper(settings, solver_name)

class ParticleMechanicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the ParticleMechanicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return ParticleMechanicsAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        coupling_model_part = self.model.GetModelPart("MPM_Coupling_Interface")
        model_part_name = self.project_parameters["coupling_settings"]["interface_model_part_name"].GetString() # TODO this should be specified in "solver_wrapper_settings" in teh cosim-json
        model_part = self.model.GetModelPart(model_part_name)

        ## Transfer information from coupling_mp to mp
        for coupling_node in coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id

            ## IMPOSED DISPLACEMENT
            total_displacement = coupling_node.GetSolutionStepValue(KM.DISPLACEMENT,0)
            old_displacement = model_part.GetCondition(coupling_id).GetValue(KPM.MPC_DISPLACEMENT)
            incremental_displacement = total_displacement - old_displacement
            model_part.GetCondition(coupling_id).SetValue(KPM.MPC_IMPOSED_DISPLACEMENT,incremental_displacement)

            ## ADD VELOCITY
            current_velocity = coupling_node.GetSolutionStepValue(KM.VELOCITY,0)
            model_part.GetCondition(coupling_id).SetValue(KPM.MPC_VELOCITY, current_velocity)

            ## ADD NORMAL
            normal = coupling_node.GetSolutionStepValue(KM.NORMAL,0)
            # Check and see whether the normal is not zero
            norm_normal = math.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
            if norm_normal > 1.e-10:
                model_part.GetCondition(coupling_id).SetValue(KPM.MPC_NORMAL, normal)

        super(ParticleMechanicsWrapper, self).SolveSolutionStep()

        ### Get contact force from mp to coupling_mp
        for mpc in model_part.Conditions:
            if (mpc.Is(KM.INTERFACE)):
                coupling_id   = mpc.Id
                contact_force = mpc.GetValue(KPM.MPC_CONTACT_FORCE)
                coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KM.CONTACT_FORCE,0,contact_force)
