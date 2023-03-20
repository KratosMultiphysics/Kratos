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
#include "container_variable_data.h"

namespace Kratos {

///@name Kratos Classes
///@{

template <class TContainerType, class TContainerDataIO>
class KRATOS_API(OPTIMIZATION_APPLICATION) SpecializedContainerVariableData : public ContainerVariableData<TContainerType> {
public:
    ///@name Type definitions
    ///@{

    using BaseType = ContainerVariableData<TContainerType>;

    using IndexType = std::size_t;

    /// Pointer definition of SpecializedContainerVariableData
    KRATOS_CLASS_POINTER_DEFINITION(SpecializedContainerVariableData);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    SpecializedContainerVariableData(ModelPart& rModelPart)
        : BaseType(rModelPart)
    {
    }

    /// Copy constructor with base class used to transfer data between compatible data containers
    SpecializedContainerVariableData(const BaseType& rOther)
        : BaseType(rOther)
    {
    }

    SpecializedContainerVariableData& operator=(const SpecializedContainerVariableData& rOther);

    ///@}
    ///@name Public operations
    ///@{

    SpecializedContainerVariableData::Pointer Clone() const;

    SpecializedContainerVariableData::Pointer CloneWithDataInitializedToZero() const;

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

    SpecializedContainerVariableData operator+(const SpecializedContainerVariableData& rOther) const;

    SpecializedContainerVariableData& operator+=(const SpecializedContainerVariableData& rOther);

    SpecializedContainerVariableData operator+(const double Value) const;

    SpecializedContainerVariableData& operator+=(const double Value);

    SpecializedContainerVariableData operator-(const SpecializedContainerVariableData& rOther) const;

    SpecializedContainerVariableData& operator-=(const SpecializedContainerVariableData& rOther);

    SpecializedContainerVariableData operator-(const double Value) const;

    SpecializedContainerVariableData& operator-=(const double Value);

    SpecializedContainerVariableData operator*(const SpecializedContainerVariableData& rOther) const;

    SpecializedContainerVariableData& operator*=(const SpecializedContainerVariableData& rOther);

    SpecializedContainerVariableData operator*(const double Value) const;

    SpecializedContainerVariableData& operator*=(const double Value);

    SpecializedContainerVariableData operator/(const SpecializedContainerVariableData& rOther) const;

    SpecializedContainerVariableData& operator/=(const SpecializedContainerVariableData& rOther);

    SpecializedContainerVariableData operator/(const double Value) const;

    SpecializedContainerVariableData& operator/=(const double Value);

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
    const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos
