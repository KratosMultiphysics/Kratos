import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.transformation_techniques.transformation_technique import TransformationTechnique
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolderUnion

class VertexMorphingTransformationTechnique(TransformationTechnique):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.parameters = parameters
        self.model = model
        self.__mappers = {}

    def InitializeSolutionStep(self):
        # update the mappers
        for __mapper in self.__mappers.values():
            __mapper.Update()

    def TransformSensitivity(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        self.__GetMapper(container_variable_data_holder).Map(container_variable_data_holder.Clone(), container_variable_data_holder)

    def TransformUpdate(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        self.__GetMapper(container_variable_data_holder).InverseMap(container_variable_data_holder, container_variable_data_holder.Clone())

    def __GetMapper(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        if container_variable_data_holder.GetContainer() not in self.__mappers.keys():
            if isinstance(container_variable_data_holder.GetContainer(), Kratos.ModelPart.Nodes):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.VertexMorphingNodalContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.parameters.Clone())
            elif isinstance(container_variable_data_holder.GetContainer(), Kratos.ModelPart.Conditions):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.VertexMorphingConditionContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.parameters.Clone())
            elif isinstance(container_variable_data_holder.GetContainer(), Kratos.ModelPart.Elements):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.VertexMorphingElementContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.parameters.Clone())

            self.__mappers[container_variable_data_holder.GetContainer()].Update()

        return self.__mappers[container_variable_data_holder.GetContainer()]