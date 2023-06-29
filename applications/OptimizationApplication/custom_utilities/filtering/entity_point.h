//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "geometries/point.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

template<class TEntityType>
class KRATOS_API(OPTIMIZATION_APPLICATION) EntityPoint : public Point
{
public:
    ///@name Type definitions
    ///@{

    using EntityType = TEntityType;

    /// Pointer definition of ContainerData
    KRATOS_CLASS_POINTER_DEFINITION(EntityPoint);

    ///@}
    ///@name Life cycle
    ///@{

    EntityPoint() = default;

    EntityPoint(
        const TEntityType& rEntity,
        const IndexType Id);

    IndexType Id() const;

    const TEntityType& GetEntity() const;

    ///@}
private:
    ///@name Private member variables
    ///@{

    const IndexType mId = 0;

    const TEntityType* mpEntity;

    ///@}
    ///@name Private static methods
    ///@{

    static Point GetPoint(const TEntityType& rEntity);

    ///@}
};

///@}
} // namespace Kratos