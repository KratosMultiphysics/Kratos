__all__ = [
    "ContainerProxy"
]

# --- Core Imports ---
import KratosMultiphysics

# --- STD Imports ---
import typing


##! @addtogroup KratosCore
##! @{


class ContainerProxy:
    """ @brief Wrapper class representing an array of @ref Node, @ref Element, or @ref Condition s.
        @classname ContainerProxy
        @details @ref ContainerProxy provides a uniform interface for the historical/non-historical node,
                 element, and condition containers of a @ref ModelPart. @ref ContainerProxy can be iterated
                 on and yields @ref DynamicEntityProxy instances that provide a generic interface for each
                 entity type (historical/non-historical nodes, elements, and conditions).
    """

    def __init__(self,
                 data_location: KratosMultiphysics.Globals.DataLocation,
                 model_part: KratosMultiphysics.ModelPart):
        self.__data_location = data_location
        self.__container: typing.Union[KratosMultiphysics.Node,KratosMultiphysics.Element,KratosMultiphysics.Condition]
        if data_location in (KratosMultiphysics.Globals.DataLocation.NodeHistorical, KratosMultiphysics.Globals.DataLocation.NodeNonHistorical):
            self.__container = model_part.Nodes
        elif data_location == KratosMultiphysics.Globals.DataLocation.Element:
            self.__container = model_part.Elements
        elif data_location == KratosMultiphysics.Globals.DataLocation.Condition:
            self.__container = model_part.Conditions
        else:
            raise ValueError(f"Invalid data location '{data_location}' (expecting 'NodeHistorical', 'NodeNonHistorical', 'Element', or 'Condition')")


    def __iter__(self) -> "Iterator":
        return ContainerProxy.Iterator(self.__data_location, self.__container)


    def __len__(self) -> int:
        return len(self.__container)


    def __contains__(self, entity: typing.Union[KratosMultiphysics.Node,KratosMultiphysics.Element,KratosMultiphysics.Condition]) -> bool:
        return entity in self.__container


    class Iterator:

        def __init__(self,
                     data_location: KratosMultiphysics.Globals.DataLocation,
                     container: typing.Union[KratosMultiphysics.NodesArray,KratosMultiphysics.ElementsArray,KratosMultiphysics.ConditionsArray]):
            self.__data_location = data_location
            self.__wrapped = container.__iter__()


        def __next__(self) -> KratosMultiphysics.EntityProxy:
            return KratosMultiphysics.EntityProxy(self.__data_location, self.__wrapped.__next__())


##! @}
