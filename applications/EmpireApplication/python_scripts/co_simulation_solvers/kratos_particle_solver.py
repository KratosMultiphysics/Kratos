from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from kratos_base_field_solver import KratosBaseFieldSolver

# Other imports
from KratosMultiphysics.ParticleMechanicsApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

import math

def CreateSolver(cosim_solver_settings, level):
    return KratosParticleSolver(cosim_solver_settings, level)

class KratosParticleSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return ParticleMechanicsAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        coupling_model_part = self.model.GetModelPart("MPM_Coupling_Interface")
        model_part_name = self.project_parameters["coupling_settings"]["interface_model_part_name"].GetString()
        model_part = self.model.GetModelPart(model_part_name)

        ## Transfer information from coupling_mp to mp
        for coupling_node in coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id

            ## IMPOSED DISPLACEMENT
            total_displacement = coupling_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0)
            old_displacement = model_part.GetCondition(coupling_id).GetValue(KratosParticle.MPC_DISPLACEMENT)
            incremental_displacement = total_displacement - old_displacement
            model_part.GetCondition(coupling_id).SetValue(KratosParticle.MPC_IMPOSED_DISPLACEMENT,incremental_displacement)

            ## ADD VELOCITY
            current_velocity = coupling_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
            model_part.GetCondition(coupling_id).SetValue(KratosParticle.MPC_VELOCITY, current_velocity)

            ## ADD NORMAL
            normal = coupling_node.GetSolutionStepValue(KratosMultiphysics.NORMAL,0)
            # Check and see whether the normal is not zero
            norm_normal = math.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])
            if norm_normal > 1.e-10:
                model_part.GetCondition(coupling_id).SetValue(KratosParticle.MPC_NORMAL, normal)

        super(KratosParticleSolver, self).SolveSolutionStep()

        ### Get contact force from mp to coupling_mp
        for mpc in model_part.Conditions:
            if (mpc.Is(KratosMultiphysics.INTERFACE)):
                coupling_id   = mpc.Id
                contact_force = mpc.GetValue(KratosParticle.MPC_CONTACT_FORCE)
                coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE,0,contact_force)


    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__

