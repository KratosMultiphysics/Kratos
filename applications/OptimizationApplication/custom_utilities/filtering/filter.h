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

    /**
     * @brief Check method used to check the inputs for the filter.
     */
    virtual void Check()
    {
        KRATOS_ERROR << "Calling base class Filter::Check method. Please define it in derived class.";
    }

    virtual void Initialize()
    {
        KRATOS_ERROR << "Calling base class Filter::Initialize method. Please define it in derived class.";
    }

    virtual void Finalize()
    {
        KRATOS_ERROR << "Calling base class Filter::Finalize method. Please define it in derived class.";
    }

    /**
     * @brief Informs filter to update its internal data because the physical space design has changed.
     *
     */
    virtual void Update()
    {
        KRATOS_ERROR << "Calling base class Filter::Update method. Please define it in derived class.";
    }

    /**
     * @brief Returns the filtered field.
     */
    virtual ContainerExpression<TContainerType> FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
    {
        KRATOS_ERROR << "Calling base class Filter::FilterField method. Please define it in derived class.";
    }

    /**
     * @brief Returns a filtered field which is already an integrated field.
     */
    virtual ContainerExpression<TContainerType> FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const
    {
        KRATOS_ERROR << "Calling base class Filter::FilterIntegratedField method. Please define it in derived class.";
    }

    virtual std::string Info() const
    {
        KRATOS_ERROR << "Calling base class Filter::Info method. Please define it in derived class.";
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