# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
import KratosMultiphysics as KM
# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing StructuralMechanics
if not CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    raise ImportError("The StructuralMechanicsApplication is not available!")
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

def Create(settings, model, solver_name):
    return StructuralMechanicsWrapper(settings, model, solver_name)

class StructuralMechanicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the StructuralMechanicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return StructuralMechanicsAnalysis(self.model, self.project_parameters)

    def SolveSolutionStep(self):
        

        super().SolveSolutionStep()
        coupling_model_part = self.model.GetModelPart("Structure.LineLoad2D_Load_on_lines_Auto1")
        
        ## Transfer information from coupling_mp to mp
        for coupling_node in coupling_model_part.Nodes:
            print(coupling_node.GetSolutionStepValue(KM.DISPLACEMENT_Y))
            
