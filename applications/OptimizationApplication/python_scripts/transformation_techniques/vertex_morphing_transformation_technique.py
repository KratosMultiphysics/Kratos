from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.transformation_techniques.transformation_technique import TransformationTechnique
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerTypes

class VertexMorphingTransformationTechnique(TransformationTechnique):
    __MapperTypes = Union[
        KratosOA.Mappers.VertexMorphingNodalContainerVariableDataMapper,
        KratosOA.Mappers.VertexMorphingConditionContainerVariableDataMapper,
        KratosOA.Mappers.VertexMorphingElementContainerVariableDataMapper
    ]

    def __init__(self, _: Kratos.Model, parameters: Kratos.Parameters, __: OptimizationInfo):
        super().__init__()

        self.__parameters = parameters
        self.__mappers: 'dict[ContainerTypes, VertexMorphingTransformationTechnique.__MapperTypes]' = {}

    def ExecuteInitializeSolutionStep(self):
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
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.Mappers.VertexMorphingNodalContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.__parameters.Clone())
            elif isinstance(container_variable_data_holder.GetContainer(), Kratos.ConditionsArray):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.Mappers.VertexMorphingConditionContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.__parameters.Clone())
            elif isinstance(container_variable_data_holder.GetContainer(), Kratos.ElementsArray):
                self.__mappers[container_variable_data_holder.GetContainer()] = KratosOA.Mappers.VertexMorphingElementContainerVariableDataMapper(container_variable_data_holder.GetModelPart(), container_variable_data_holder.GetModelPart(), self.__parameters.Clone())

            self.__mappers[container_variable_data_holder.GetContainer()].Update()

        return self.__mappers[container_variable_data_holder.GetContainer()]