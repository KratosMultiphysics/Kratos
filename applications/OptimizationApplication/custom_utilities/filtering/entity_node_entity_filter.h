//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "expression/container_expression.h"

// Application includes
#include "custom_utilities/filtering/filter.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) EntityNodeEntityFilter: public Filter<TContainerType>
{
public:
    ///@name Type definitions
    ///@{

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(EntityNodeEntityFilter);

    ///@}
    ///@name LifeCycle
    ///@{

    EntityNodeEntityFilter(ModelPart& rModelPart);

    ~EntityNodeEntityFilter() override = default;

    ///@}
    ///@name Public operations

    void Update() override;

    ContainerExpression<TContainerType> FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const override;

    ContainerExpression<TContainerType> FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const override;

    std::string Info() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart& mrModelPart;

    ContainerExpression<ModelPart::NodesContainerType> mNeighbourEntities;

    ContainerExpression<TContainerType> mEntityDomainSize;

    ///@}
};

///@}
} // namespace Kratos