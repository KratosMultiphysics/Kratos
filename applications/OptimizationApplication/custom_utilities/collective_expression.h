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
#include "expression/container_expression.h"
#include "expression/traits.h"
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "collective_expression_arithmetic_operators.h"

namespace Kratos {

///@name Kratos Classes
///@{
/**
 * @brief Construct a new CollectiveExpression instance
 * @details This constructs a CollectiveExpression instance which can hold list of ContainerExpresions of different types.
 *          The list within the instance can be treated as a single value, therefore all the binary operations present in the ContainerExpressions
 *          are also available with this container. This uses std::variant to store different types of containers. Since these containers memeory footprint
 *          is the same, the memory footprint of the std::variant will be almost equal to the memory footprint of a single container for all the container types.
 *
 */
class KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using CollectiveExpressionType = std::variant<
            ContainerExpression<ModelPart::NodesContainerType>::Pointer,
            ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
            ContainerExpression<ModelPart::ElementsContainerType>::Pointer>;

    KRATOS_CLASS_POINTER_DEFINITION(CollectiveExpression);

    ///@}
    ///@name Life cycle
    ///#{

    // Default constructor
    CollectiveExpression() noexcept = default;


    /**
     * @brief Construct a list of container expressions
     *
     * @param rExpressionPointersList       List of container expressions
     */
    CollectiveExpression(const std::vector<CollectiveExpressionType>& rExpressionPointersList);

    // Copy consttructor
    CollectiveExpression(const CollectiveExpression& rOther);

    /// Assignment operator
    CollectiveExpression& operator=(const CollectiveExpression& rOther);

    // destructor
    ~CollectiveExpression() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clones this collective expressions by cloning all the ContainerExpressions in the instance.
     *
     * This clones current instance of CollectiveExpression by cloning all the ContainerExpressions. Cloning
     * a single ContainerExpression is a light weight operation, hence this operation also a light weight operation.
     *
     * @return CollectiveExpression        Cloned collective expressions which has all the clones of ContainerExpressions.
     */
    CollectiveExpression Clone() const;

    /**
     * @brief Add a ContainerExpression to the current instance.
     *
     * @param pExpression   The ContainerExpression to be added.
     */
    void Add(const CollectiveExpressionType& pExpression);

    /**
     * @brief Add a CollectiveExpression to the current instance.
     *
     * This method appends the current instance's list of ContainerExpressions with the new ContainerExpressions
     * in @ref rCollectiveExpression.
     *
     * @param rCollectiveExpression     The CollectiveExpression to be appended to the current CollectiveExpression.
     */
    void Add(const CollectiveExpression& rCollectiveExpression);

    /**
     * @brief Clear list of all the ContainerExpressions.
     *
     * This does not destroy the corresponding ContainerExpression (uses SmartPointers). This removes
     * all the ContainerExpression objects from the current list.
     *
     */
    void Clear();

    /**
     * @brief Get the Collective Flattened Data Size.
     *
     * This method returns the total number of double values present in all the ContainerExpressions.
     *
     * @return IndexType        Total number of double values present in all the ContainerExpressions.
     */
    IndexType GetCollectiveFlattenedDataSize() const;

    std::vector<CollectiveExpressionType> GetContainerExpressions();

    std::vector<CollectiveExpressionType> GetContainerExpressions() const;

    bool IsCompatibleWith(const CollectiveExpression& rOther) const;

    ///@}
    ///@name Public operators
    ///@{

    CollectiveExpression& operator+=(const CollectiveExpression& rOther);

    CollectiveExpression& operator+=(const double Value);

    CollectiveExpression& operator-=(const CollectiveExpression& rOther);

    CollectiveExpression& operator-=(const double Value);

    CollectiveExpression& operator*=(const CollectiveExpression& rOther);

    CollectiveExpression& operator*=(const double Value);

    CollectiveExpression& operator/=(const CollectiveExpression& rOther);

    CollectiveExpression& operator/=(const double Value);

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
    const CollectiveExpression& rThis)
{
    return rOStream << rThis.Info();
}

///@}

} // namespace Kratos