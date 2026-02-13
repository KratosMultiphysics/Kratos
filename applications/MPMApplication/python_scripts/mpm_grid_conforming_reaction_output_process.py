import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.mpm_reaction_output_process import MPMReactionOutputProcess

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.OutputProcess:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMGridConformingReactionOutputProcess(model, settings["Parameters"])

class MPMGridConformingReactionOutputProcess(MPMReactionOutputProcess):
    """
    The specific implementation for the output of reaction on grid conforming boundaries
    """
    def __init__(self, model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        super().__init__()
        if not self.model_part_name.startswith('Background_Grid'):
            err_msg  = f"'{self.__class__.__name__}' prints the sum of the reaction of a given SubModelPart of the background grid.\n"
            err_msg += "Model part name must start with 'Background_Grid'. Input model part name: '{self.model_part_name}'"
            raise Exception(err_msg)

    def _GetFileHeader(self):
        header  = f'# Grid conforming reaction for model part {self.model_part_name}\n'
        header +=  '# Time Fx Fy Fz\n'
        return header

    def _GetReaction(self):
        return KratosMPM.ReactionUtilities.CalculateGridConformingReaction(self.model_part)
