//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "utilities/parallel_utilities.h"

// Include base h
#include "vertex_morphing_container_variable_data_mapper.h"

namespace Kratos {

template<class TContainerType>
VertexMorphingContainerVariableDataMapper<TContainerType>::VertexMorphingContainerVariableDataMapper(
    Model& rModel,
    Parameters Params)
    : mrModel(rModel)
{

}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::Update()
{

}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::Map(
    ContainerVariableDataHolderType& rOutputDataContainer,
    const ContainerVariableDataHolderType& rInputDataContainer)
{

}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::InverseMap(
    ContainerVariableDataHolderType& rOutputDataContainer,
    const ContainerVariableDataHolderType& rInputDataContainer)
{
}

// template instantiations
template class VertexMorphingContainerVariableDataMapper<ModelPart::NodesContainerType>;
template class VertexMorphingContainerVariableDataMapper<ModelPart::ElementsContainerType>;
template class VertexMorphingContainerVariableDataMapper<ModelPart::ConditionsContainerType>;

} // namespace Kratos