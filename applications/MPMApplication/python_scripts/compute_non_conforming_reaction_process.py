# Importing the Kratos Library
import KratosMultiphysics
# Import applications
import KratosMultiphysics.MPMApplication as KratosMPM
# Import base class file
from KratosMultiphysics.MPMApplication.compute_reaction_process import ComputeReactionProcess

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeNonConformingReactionProcess(model, settings["Parameters"])

class ComputeNonConformingReactionProcess(ComputeReactionProcess):
    """
    The specific implementation for the output of reaction on non-conforming boundaries
    defined through material point conditions
    """
    def _GetFileHeader(self):
        header  = f'# Non-conforming reaction for model part {self.params["model_part_name"].GetString()}\n'
        header +=  '# Time Fx Fy Fz\n'
        return header

    def _GetReaction(self):
        return KratosMPM.ReactionUtilities.CalculateNonConformingReaction(self.model_part)
