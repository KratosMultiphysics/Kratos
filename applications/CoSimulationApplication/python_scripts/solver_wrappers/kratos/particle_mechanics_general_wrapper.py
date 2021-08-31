# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
import KratosMultiphysics.ParticleMechanicsApplication as KPM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing ParticleMechanics
if not CheckIfApplicationsAvailable("ParticleMechanicsApplication"):
    raise ImportError("The ParticleMechanicsApplication is not available!")
from KratosMultiphysics.ParticleMechanicsApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

def Create(settings, model, solver_name):
    return ParticleMechanicsWrapper(settings, model, solver_name)

class ParticleMechanicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the ParticleMechanicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return ParticleMechanicsAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        mpm_background_grid_model_part = self.model.GetModelPart("Background_Grid")

        coupling_model_part = mpm_background_grid_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto3")
        
        ## Transfer information from coupling_mp to mp
        for coupling_node in coupling_model_part.Nodes:
            
            print(coupling_node.GetSolutionStepValue(KPM.IMPOSED_DISPLACEMENT_Y))
            # print(coupling_node.GetSolutionStepValue(KM.DISPLACEMENT))

        super().SolveSolutionStep()
