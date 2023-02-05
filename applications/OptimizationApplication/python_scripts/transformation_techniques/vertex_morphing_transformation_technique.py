from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.transformation_techniques.transformation_technique import TransformationTechnique
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolderUnion
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerTypes

class VertexMorphingTransformationTechnique(TransformationTechnique):
    __MapperTypes = Union[
        KratosOA.VertexMorphingNodalContainerVariableDataMapper,
        KratosOA.VertexMorphingConditionContainerVariableDataMapper,
        KratosOA.VertexMorphingElementContainerVariableDataMapper
    ]

    def __init__(self, _: Kratos.Model, parameters: Kratos.Parameters, __: OptimizationInfo):
        self.__parameters = parameters
        self.__mappers: 'dict[ContainerTypes, VertexMorphingTransformationTechnique.__MapperTypes]' = {}

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
            if isinstance(container_variable_data_holder.GetContainer(), Kratos.NodesArray):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.VertexMorphingNodalContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.__parameters.Clone())
            elif isinstance(container_variable_data_holder.GetContainer(), Kratos.ConditionsArray):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.VertexMorphingConditionContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.__parameters.Clone())
            elif isinstance(container_variable_data_holder.GetContainer(), Kratos.ElementsArray):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.VertexMorphingElementContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.__parameters.Clone())

            self.__mappers[container_variable_data_holder.GetContainer()].Update()

        return self.__mappers[container_variable_data_holder.GetContainer()]