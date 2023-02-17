import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion

class TransformationTechnique(Kratos.Process):
    def __init__(self):
        super().__init__()

    def TransformSensitivity(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        raise NotImplementedError("Calling base class TransformationTechnique::TransformSensitivity. Please implement it in the derrived class")

    def TransformUpdate(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        raise NotImplementedError("Calling base class TransformationTechnique::TransformUpdate. Please implement it in the derrived class")