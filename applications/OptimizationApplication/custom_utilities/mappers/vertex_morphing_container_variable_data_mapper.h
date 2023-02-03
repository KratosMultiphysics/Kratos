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

#pragma once

// System includes
#include <variant>
#include <functional>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"

// Application includes
#include "custom_utilities/mappers/container_variable_data_mapper.h"
#include "custom_utilities/mappers/entity_point.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) VertexMorphingContainerVariableDataMapper : public ContainerVariableDataMapper<TContainerType>
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = ContainerVariableDataMapper<TContainerType>;

    using ContainerVariableDataHolderType = ContainerVariableDataHolderBase<TContainerType>;

    using TEntityType = typename TContainerType::data_type;

    using EntityPointType = EntityPoint<TEntityType>;

    using EntityPointTypePointer = typename EntityPointType::Pointer;

    using EntityPointVector = std::vector<EntityPointTypePointer>;

    using EntityPointIterator = typename std::vector<EntityPointTypePointer>::iterator;

    using DoubleVectorIterator = std::vector<double>::iterator;

    using  BucketType = Bucket<3, EntityPointType, EntityPointVector, EntityPointTypePointer, EntityPointIterator, DoubleVectorIterator>;

    using KDTree = Tree<KDTreePartition<BucketType>>;

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(VertexMorphingContainerVariableDataMapper);

    ///@}
    ///@name LifeCycle
    ///@{

    VertexMorphingContainerVariableDataMapper(
        Model& rModel,
        Parameters Params);

    ~VertexMorphingContainerVariableDataMapper() override = default;

    ///@}
    ///@name Public operations

    void Update() override;

    void Map(
        ContainerVariableDataHolderType& rOutputDataContainer,
        const ContainerVariableDataHolderType& rInputDataContainer) override;

    void InverseMap(
        ContainerVariableDataHolderType& rOutputDataContainer,
        const ContainerVariableDataHolderType& rInputDataContainer) override;

    ///@}
private:
    ///@name Private member variables
    ///@{

    Model& mrModel;

    ///@}
};

///@}
} // namespace Kratos