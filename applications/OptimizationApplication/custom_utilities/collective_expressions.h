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
#include <vector>
#include <variant>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "containers/container_expression/specialized_container_expression.h"
#include "containers/container_expression/container_data_io.h"

// Application includes
#include "container_properties_data_io.h"

namespace Kratos {

///@name Kratos Classes
///@{
/**
 * @brief Construct a new CollectiveExpressions instance
 * @details This constructs a CollectiveExpressions instance which can hold list of ContainerExpresions of different types.
 *          The list within the instance can be treated as a single value, therefore all the binary operations present in the ContainerExpressions
 *          are also available with this container. This uses std::variant to store different types of containers. Since these containers memeory footprint
 *          is the same, the memory footprint of the std::variant will be almost equal to the memory footprint of a single container for all the container types.
 *
 */
class KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpressions {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using HistoricalExpressionPointer =  SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>>::Pointer;

    using NodalExpressionPointer =  SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer;

    using ConditionExpressionPointer =  SpecializedContainerExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer;

    using ElementExpressionPointer =  SpecializedContainerExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer;

    using ConditionPropertiesExpressionPointer =  SpecializedContainerExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::Properties>>::Pointer;

    using ElementPropertiesExpressionPointer =  SpecializedContainerExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::Properties>>::Pointer;

    using CollectiveExpressionType = std::variant<
            HistoricalExpressionPointer,
            NodalExpressionPointer,
            ConditionExpressionPointer,
            ElementExpressionPointer,
            ConditionPropertiesExpressionPointer,
            ElementPropertiesExpressionPointer>;

    KRATOS_CLASS_POINTER_DEFINITION(CollectiveExpressions);

    ///@}
    ///@name Life cycle
    ///#{

    // Default constructor
    CollectiveExpressions() noexcept = default;


    // Constructor with list
    CollectiveExpressions(const std::vector<CollectiveExpressionType>& rExpressionPointersList);

    // Copy consttructor
    CollectiveExpressions(const CollectiveExpressions& rOther);

    // destructor
    ~CollectiveExpressions() = default;

    ///@}
    ///@name Public operations
    ///@{

    CollectiveExpressions Clone() const;

    void SetToZero();

    void Add(const CollectiveExpressionType& pExpression);

    void Add(const CollectiveExpressions& rCollectiveExpression);

    void Clear();

    std::vector<CollectiveExpressionType> GetContainerExpressions();

    std::vector<CollectiveExpressionType> GetContainerExpressions() const;

    bool IsCompatibleWith(const CollectiveExpressions& rOther) const;

    ///@}
    ///@name Public operators
    ///@{

    CollectiveExpressions operator+(const CollectiveExpressions& rOther) const;

    CollectiveExpressions& operator+=(const CollectiveExpressions& rOther);

    CollectiveExpressions operator+(const double Value) const;

    CollectiveExpressions& operator+=(const double Value);

    CollectiveExpressions operator-(const CollectiveExpressions& rOther) const;

    CollectiveExpressions& operator-=(const CollectiveExpressions& rOther);

    CollectiveExpressions operator-(const double Value) const;

    CollectiveExpressions& operator-=(const double Value);

    CollectiveExpressions operator*(const CollectiveExpressions& rOther) const;

    CollectiveExpressions& operator*=(const CollectiveExpressions& rOther);

    CollectiveExpressions operator*(const double Value) const;

    CollectiveExpressions& operator*=(const double Value);

    CollectiveExpressions operator/(const CollectiveExpressions& rOther) const;

    CollectiveExpressions& operator/=(const CollectiveExpressions& rOther);

    CollectiveExpressions operator/(const double Value) const;

    CollectiveExpressions& operator/=(const double Value);

    CollectiveExpressions Pow(const CollectiveExpressions& rOther) const;

    CollectiveExpressions Pow(const double Value) const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const;

    ///@}
private:
    ///@name Private variables
    ///@{

    std::vector<CollectiveExpressionType> mExpressionPointersList;

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const CollectiveExpressions& rThis)
{
    return rOStream << rThis.Info();
}

///@}

} // namespace Kratos