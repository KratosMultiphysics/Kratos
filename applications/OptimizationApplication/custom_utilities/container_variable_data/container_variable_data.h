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
#include <cmath>
#include <string>
#include <variant>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/container_variable_data/expressions.h"

namespace Kratos {

///@name Kratos Classes
///@{

template <class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ContainerVariableData {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(ContainerVariableData);

    ///@}
    ///@name Life cycle
    ///#{

    virtual ~ContainerVariableData() = default;

    ///@}
    ///@name Public operations
    ///@{

    void CopyDataFrom(const ContainerVariableData<TContainerType>& rOther);

    void SetDataToZero();

    ///@}
    ///@name Input and output
    ///@{

    void SetExpression(Expression::Pointer pExpression);

    const Expression& GetExpression() const;

    const Expression::Pointer pGetExpression() const;

    IndexType GetDataDimension() const;

    ModelPart& GetModelPart();

    const ModelPart& GetModelPart() const;

    TContainerType& GetContainer();

    const TContainerType& GetContainer() const;

    virtual std::string Info() const;

    std::string PrintData() const;

    ///@}
protected:
    ///@name Life cycle
    ///@{

    /// Constructor with the model part
    ContainerVariableData(ModelPart& rModelPart);

    /// Copy constructor
    ContainerVariableData(const ContainerVariableData& rOther);

    ///@}
    ///@name Protected member variables
    ///@{

    Expression::Pointer mpExpression = nullptr;

    ModelPart* const mpModelPart;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ContainerVariableData<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}

} // namespace Kratos