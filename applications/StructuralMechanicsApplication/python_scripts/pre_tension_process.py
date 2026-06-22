# --- Core Imports ---
import KratosMultiphysics

# --- Structural Imports ---
from KratosMultiphysics.StructuralMechanicsApplication import InsertPreTensionOperation


class PreTensionProcess(KratosMultiphysics.Process):
    """ @see @ref pre_tensioning "Pre-Tensioning"
        @ingroup pre_tensioning
    """

    def __init__(
        self,
        model: KratosMultiphysics.Model,
        parameters: KratosMultiphysics.Parameters) -> None:
            super().__init__()
            parameters.AddMissingParameters(self.GetDefaultParameters())
            self.__operation: InsertPreTensionOperation = self.__MakeOperation(model, parameters)

            # The operation needs to execute upon construction.
            # Why? Because it should be a "Modeler", but modeler's
            # don't have access to DoFs yet. And also because other
            # processes immediately check for existing sub model parts
            # instead of giving a chance to other processes do their
            # jobs. So here we are.
            self.__operation.Execute()


    def ExecuteBeforeSolutionLoop(self) -> None:
        #self.__operation.Execute()
        pass


    def GetDefaultParameters(self) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters(R"""{
            "model_part_name" : "",
            "verbosity" : 1
        }""")


    @classmethod
    def __MakeOperation(
        cls,
        model: KratosMultiphysics.Model,
        parameters: KratosMultiphysics.Parameters) -> InsertPreTensionOperation:
            return InsertPreTensionOperation(model, parameters)


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    return PreTensionProcess(
        model,
        parameters["Parameters"])
