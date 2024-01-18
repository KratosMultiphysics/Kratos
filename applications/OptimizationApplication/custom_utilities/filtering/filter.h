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
#include "expression/container_expression.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) Filter
{
public:
    ///@name Type definitions
    ///@{

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(Filter);

    ///@}
    ///@name LifeCycle
    ///@{

    Filter() = default;

    virtual ~Filter() = default;

    ///@}
    ///@name Public operations

    virtual void Update()
    {
        KRATOS_ERROR << "Calling base class Filter::Update method. Please define it in derrived class.";
    }

    virtual ContainerExpression<TContainerType> FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
    {
        KRATOS_ERROR << "Calling base class Filter::FilterField method. Please define it in derrived class.";
    }

    virtual ContainerExpression<TContainerType> FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const
    {
        KRATOS_ERROR << "Calling base class Filter::FilterIntegratedField method. Please define it in derrived class.";
    }

    virtual std::string Info() const
    {
        KRATOS_ERROR << "Calling base class Filter::Info method. Please define it in derrived class.";
    }

    ///@}
};

///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const Filter<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}

///@}
} // namespace Kratos