//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
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
#include "containers/container_expression/expressions/expression.h"

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
 *      1. A literal expression is a concrete value expression. Followints are the eg. types:
 *              LiteralExpression<double>: This expression holds a single double value for
 *                                       all the entities of the TContainerType. This is light weight.
 *              LiteralExpression<array_1d<double, 3>>: This expression holds a single array3 value for
 *                                       all the entities of the TContainerType. This is light weight.
 *              LiteralFlatExpression: This expression can hold double, array3 different values for
 *                                       all the entities of the TContainerType. In this case, the dimension
 *                                       of the original data is stored, and then higher dimensional entity data
 *                                       is flattened out to a single vector. This will occupy almost the same memory as
 *                                       the variable data is stored in nodes/conditions/elements.
 *      2. A binary expression: These expressions does not hold any values, hence they are light weight.
 *                              They are used to store keep track of the operations carried on the LiteralExpressions.
 *
 * Literal expressions are created on the following cases:
 *      1. When data is set using the SpecializedContainerExpression::Read. Here, a LiteralFlatExpression is created.
 *      2. When data is reset to one value using either ContainerExpression::SetDataToZero or SpecializedContainerExpression::SetData
 *         or SpecializedContainerExpression::SetZero. Here a LiteralExpression<double> or LiteralExpression<array_1d<double, 3>> is created.
 *      3. When a SpecializedContainerExpression is used with "+", "-", "*", "/", "Pow" operators with double values in right operand.
 *
 * BinaryExpressions are created on the followin cases:
 *      1. When a SpecializedContainerExpression is operated with "+", "-", "*", "/", "Pow".
 *
 * ContainerExpression only holds double vector if any nodal, condition or element variable data needs to be stored for future calculations
 * where the variable can be released to store new data. Hence same variable can be used to store different data in the same container.
 *
 * When operatos such as "+", "-", "*", "/", "Pow" are used on these containers in python or c++, it creates an expression. They are not
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
class KRATOS_API(KRATOS_CORE) ContainerExpression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(ContainerExpression);

    ///@}
    ///@name Life cycle
    ///#{

    virtual ~ContainerExpression() = default;

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
    void CopyFrom(const ContainerExpression<TContainerType>& rOther);

    /**
     * @brief Reads data from c like interface
     *
     * This method can read data from a c-interface where @ref pBegin is the
     * starting pointer to a contiguous array having space and data
     * to all @ref NumberOfEntities where each entity having data for
     * shape @ref rShape.
     *
     * Eg: 1) If NumberOfEntities = 10, and rShape = {4}, then
     *        pBegin should point to starting double value's pointer
     *        where it has 40 doubles in a contiguous manner.
     *     2) iF NumberOfEntities = 10 and rShape = (4, 3), then
     *        pBegin should point to starting double value's pointer
     *        where it has 120 doubles in a contiguous manner.
     *        The matrix within these containers are using row first notation.
     *              [
     *                  [1, 2, 3]
     *                  [3, 4, 5]
     *              ] = {1, 2, 3, 3, 4, 5}
     *
     * @param pBegin            Starting pointer to the data.
     * @param NumberOfEntities  Number of entities present in data.
     * @param pShapeBegin       Starting  point of the shape of data in each entity.
     * @param ShapeSize         Size of the shape.
     */
    void Read(
        double const* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    /**
     * @brief Move data from pBegin array to internal structure.
     *
     * @warning This instance takes ownership of the passed array.
     *
     * @param pBegin            Starting pointer to the data.
     * @param NumberOfEntities  Number of entities present in data.
     * @param pShapeBegin       Starting  point of the shape of data in each entity.
     * @param ShapeSize         Size of the shape.
     */
    void MoveFrom(
        double* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize);

    /**
     * @brief Assign the data in the expression to c-like interfaces
     *
     * This method can assign data to a c-interface where @ref pBegin is the
     * starting pointer to a contiguous array having space and data
     * to all @ref NumberOfEntities where each entity having data for
     * shape @ref rShape. This is checking whether the internal @ref NumberOfEntities
     * and internal @ref Expression shape is matching with the user input.
     *
     * Eg: 1) If NumberOfEntities = 10, and rShape = {4}, then
     *        pBegin should point to starting double value's pointer
     *        where it has 40 doubles in a contiguous manner.
     *     2) iF NumberOfEntities = 10 and rShape = (4, 3), then
     *        pBegin should point to starting double value's pointer
     *        where it has 120 doubles in a contiguous manner.
     *        The matrix within these containers are using row first notation.
     *              [
     *                  [1, 2, 3]
     *                  [3, 4, 5]
     *              ] = {1, 2, 3, 3, 4, 5}
     *
     * @param pBegin            Starting pointer to the data.
     * @param NumberOfEntities  Number of entities present in data.
     * @param rShape            Shape of data in each entity.
     */
    void Evaluate(
        double* pBegin,
        const int NumberOfEntities,
        int const* pShapeBegin,
        const int ShapeSize) const;

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
    IndexType GetFlattenedSize() const;

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
    ContainerExpression(ModelPart& rModelPart);

    /// Copy constructor
    ContainerExpression(const ContainerExpression& rOther);

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
    const ContainerExpression<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}

} // namespace Kratos