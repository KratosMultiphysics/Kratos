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
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "container_variable_data_holder_base.h"

namespace Kratos {

///@name Kratos Classes
///@{

template <class TContainerType, class TContainerDataIO>
class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableDataHolder : public ContainerVariableDataHolderBase<TContainerType> {
public:
    ///@name Type definitions
    ///@{

    using BaseType = ContainerVariableDataHolderBase<TContainerType>;

    using IndexType = std::size_t;

    /// Pointer definition of ContainerVariableDataHolder
    KRATOS_CLASS_POINTER_DEFINITION(ContainerVariableDataHolder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ContainerVariableDataHolder(ModelPart& rModelPart)
        : BaseType(rModelPart)
    {
    }

    /// Copy constructor
    ContainerVariableDataHolder(const ContainerVariableDataHolder& rOther)
        : BaseType(rOther)
    {
    }

    /// Copy constructor with base class used to transfer data between compatible data containers
    ContainerVariableDataHolder(const BaseType& rOther)
        : BaseType(rOther)
    {
    }

    ~ContainerVariableDataHolder() override = default;

    ///@}
    ///@name Public operations
    ///@{

    ContainerVariableDataHolder::Pointer Clone() const;

    ContainerVariableDataHolder::Pointer CloneWithDataInitializedToZero() const;

    template <class TDataType>
    void ReadDataFromContainerVariable(
        const Variable<TDataType>& rVariable);

    template <class TDataType>
    void AssignDataToContainerVariable(
        const Variable<TDataType>& rVariable);

    template <class TDataType>
    void SetDataForContainerVariable(
        const Variable<TDataType>& rVariable,
        const TDataType& rValue);

    template <class TDataType>
    void SetDataForContainerVariableToZero(
        const Variable<TDataType>& rVariable);

    ///@}
    ///@name Operators
    ///@{

    ContainerVariableDataHolder operator+(const ContainerVariableDataHolder& rOther) const;

    ContainerVariableDataHolder& operator+=(const ContainerVariableDataHolder& rOther);

    ContainerVariableDataHolder operator+(const double Value) const;

    ContainerVariableDataHolder& operator+=(const double Value);

    ContainerVariableDataHolder operator-(const ContainerVariableDataHolder& rOther) const;

    ContainerVariableDataHolder& operator-=(const ContainerVariableDataHolder& rOther);

    ContainerVariableDataHolder operator-(const double Value) const;

    ContainerVariableDataHolder& operator-=(const double Value);

    ContainerVariableDataHolder operator*(const ContainerVariableDataHolder& rOther) const;

    ContainerVariableDataHolder& operator*=(const ContainerVariableDataHolder& rOther);

    ContainerVariableDataHolder operator*(const double Value) const;

    ContainerVariableDataHolder& operator*=(const double Value);

    ContainerVariableDataHolder operator/(const ContainerVariableDataHolder& rOther) const;

    ContainerVariableDataHolder& operator/=(const ContainerVariableDataHolder& rOther);

    ContainerVariableDataHolder operator/(const double Value) const;

    ContainerVariableDataHolder& operator/=(const double Value);

    ContainerVariableDataHolder operator^(const ContainerVariableDataHolder& rOther) const;

    ContainerVariableDataHolder& operator^=(const ContainerVariableDataHolder& rOther);

    ContainerVariableDataHolder operator^(const double Value) const;

    ContainerVariableDataHolder& operator^=(const double Value);

    ContainerVariableDataHolder& operator=(const ContainerVariableDataHolder& rOther);

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}
};

///@}
/// output stream function
template<class TContainerType, class TContainerDataIO>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos
