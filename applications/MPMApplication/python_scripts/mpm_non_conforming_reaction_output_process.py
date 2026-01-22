import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.mpm_reaction_output_process import MPMReactionOutputProcess

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMNonConformingReactionOutputProcess(model, settings["Parameters"])

class MPMNonConformingReactionOutputProcess(MPMReactionOutputProcess):
    """
    The specific implementation for the output of reaction on non-conforming boundaries
    defined through material point conditions
    """
    def _GetFileHeader(self):
        header  = f'# Non-conforming reaction for model part {self.model_part_name}\n'
        header +=  '# Time Fx Fy Fz\n'
        return header

    def _GetReaction(self):
        return KratosMPM.ReactionUtilities.CalculateNonConformingReaction(self.model_part)
