import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.compute_reaction_process import ComputeReactionProcess

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeGridConformingReactionProcess(model, settings["Parameters"])

class ComputeGridConformingReactionProcess(ComputeReactionProcess):
    """
    The specific implementation for the output of reaction on grid conforming boundaries
    """
    def _GetFileHeader(self):
        header  = f'# Grid conforming reaction for model part {self.model_part_name}\n'
        header +=  '# Time Fx Fy Fz\n'
        return header

    def _GetReaction(self):
        return KratosMPM.ReactionUtilities.CalculateGridConformingReaction(self.model_part)
