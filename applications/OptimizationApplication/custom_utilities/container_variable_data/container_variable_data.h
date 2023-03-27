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
#include <optional>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/container_variable_data/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Container variable data holder
 *
 * Instances of this class are used to hold any type of data in a container
 * of TContainerType. TContainerType can be a container of nodes, conditions or elements.
 * The data is stored as an Expression. The expression can be one of the followings:
 *
 *      1. A literal expression is a concrete value expression. Followints are the types:
 *              LiteralExpression<double>: This expression holds a single double value for
 *                                       all the entities of the TContainerType. This is light weight.
 *              LiteralExpression<array_1d<double, 3>>: This expression holds a single array3 value for
 *                                       all the entities of the TContainerType. This is light weight.
 *              LiteralVectorExpression: This expression can hold double, array3 different values for
 *                                       all the entities of the TContainerType. In this case, the dimension
 *                                       of the original data is stored, and then higher dimensional entity data
 *                                       is flattened out to a single vector. This will occupy almost the same memory as
 *                                       the variable data is stored in nodes/conditions/elements.
 *      2. A binary expression: These expressions does not hold any values, hence they are light weight.
 *                              They are used to store keep track of the operations carried on the LiteralExpressions.
 *
 * Literal expressions are created on the following cases:
 *      1. When data is set using the SpecializedContainerVariableData::ReadData. Here, a LiteralVectorExpression is created.
 *      2. When data is reset to one value using either ContainerVariableData::SetDataToZero or SpecializedContainerVariableData::SetData
 *         or SpecializedContainerVariableData::SetZero. Here a LiteralExpression<double> or LiteralExpression<array_1d<double, 3>> is created.
 *      3. When a SpecializedContainerVariableData is used with "+", "-", "*", "/", "Pow" operators with double values in right operand.
 *
 * BinaryExpressions are created on the followin cases:
 *      1. When a SpecializedContainerVariableData is operated with "+", "-", "*", "/", "Pow".
 *
 * ContainerVariableData only holds double vector if any nodal, condition or element variable data needs to be stored for future calculations
 * where the variable can be released to store new data. Hence same variable can be used to store different data in the same container.
 *
 * When operatos such as "+", "-", "*", "/", "Pow" are used on these contaienrs in python or c++, it creates an expression. They are not
 * evaluated and the results are not stored on vectors hence avoiding unnecessary computations and unnecessary memory allocations and de-allocations.
 *
 * This class's constructors are protected, hence no objects of this class can be created. This class is used as the common
 * interface to transfer data between compatible data containers such as:
 *      1. Nodal historical and non-historical.
 *      2. Element and element properties.
 *      3. Condition and condition properties.
 *
 * This class is optimized and compatible with OpenMP and MPI.
 *
 * @tparam TContainerType
 */
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

    /**
     * @brief Copies the data from another same type container variable data.
     *
     * This method is used to copy data from another variable data. The model
     * parts should be matching to successfully copy data. This does not
     * copy the data, it only copies the expression pointer, hence this
     * operation is also a light weight operation.
     *
     * @param rOther        Other container variable data
     */
    void CopyDataFrom(const ContainerVariableData<TContainerType>& rOther);

    /**
     * @brief Set the Data To Zero container.
     *
     * This replaces the expression with a LiteralExpression<double> with 0.0.
     * This is also a light weight operation since it only creates on literal with
     * one 0.0 value.
     *
     */
    void SetDataToZero();

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Set the Expression of the container data
     *
     * @param pExpression       Expression to be used in this container variable data
     */
    void SetExpression(Expression::Pointer pExpression);

    /**
     * @brief Get the Expression
     *
     * @return const Expression&    Returns the reference of the expression
     */
    const Expression& GetExpression() const;

    /**
     * @brief Get the expression pointer
     *
     * @return const Expression::Pointer Returns the pointer of the expression
     */
    const Expression::Pointer pGetExpression() const;

    /**
     * @brief Get the shape of the expression data
     *
     * This returns the dimension of the underlying expression in each element of the vector.
     *
     * @return const std::vector<IndexType>
     */
    const std::vector<IndexType> GetShape() const;

    /**
     * @brief Get the Local Size of the data
     *
     * This returns the local size which is the sum of available indices in each dimension.
     *
     * @return IndexType
     */
    IndexType GetLocalSize() const;

    /**
     * @brief Get the Model Part used in the container data
     *
     * @return ModelPart&       Model part
     */
    ModelPart& GetModelPart();

    /**
     * @brief Get the Model Part used in the container
     *
     * @return const ModelPart& Model part
     */
    const ModelPart& GetModelPart() const;

    /**
     * @brief Get the Container of the model part
     *
     * This returns the container of the model part on which this
     * container object is responsible for. It always returns the
     * local mesh container.
     *
     * @return TContainerType&      Container of the model part
     */
    TContainerType& GetContainer();

    /**
     * @brief Get the Container of the model part
     *
     * This returns the container of the model part on which this
     * container object is responsible for. It always returns the
     * local mesh container.
     *
     * @return const TContainerType&  Container of the model part
     */
    const TContainerType& GetContainer() const;

    /**
     * @brief Get the info string
     *
     * @return std::string
     */
    virtual std::string Info() const;

    /**
     * @brief Prints containing data
     *
     * @return std::string
     */
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

    std::optional<Expression::Pointer> mpExpression;

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