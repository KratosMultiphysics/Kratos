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
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "containers/container_expression/container_expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Construct a new specialized @ref ContainerExpression.
 *
 * This class is used to represent objects which are the specializations
 * of ContainerExpression. The base class is @ref ContainerExpression<TContainerType>.
 * Hence this class can be specialized to read and write data to different entity containers
 * for each container type.
 *
 * This class does not have any member variables, instead it uses expressions to represent the
 * data. Therefore, all the operations such as "+", "-", "*", etc are light weight because
 * they do not compute the value in place, they just create an expression which keeps track fo the expressions.
 * These expressions are evaluated only if @ref Evaluate. This also does not create additional
 * vectors to hold the resultant value of the expression. It uses the model parts respective containers entity input/output
 * method specified to write evaluated entity resultant values to model part entities. Hence, these
 * SpecializedContainerExpression can be easily visualized using common variable data visualization methods.
 *
 * Copy constructor is introduced with the base class type because, this allows copying data
 * between compatible @ref SpecializedContainerExpression containers. This copy is also light weight since
 * this only copies the pointers, not the data itself.
 *
 * This class can take advantage of OpenMP and MPI.
 *
 * @tparam TContainerType           Container type.
 * @tparam TContainerDataIO         Container entity input/output type.
 * @tparam TMeshType                Mesh type, should be Local, Ghost or Interface.
 */
template <class TContainerType, class TContainerDataIO, class TMeshType = MeshType::Local>
class SpecializedContainerExpression : public ContainerExpression<TContainerType, TMeshType> {
public:
    ///@name Type definitions
    ///@{

    using BaseType = ContainerExpression<TContainerType, TMeshType>;

    using IndexType = std::size_t;

    /// Pointer definition of SpecializedContainerExpression
    KRATOS_CLASS_POINTER_DEFINITION(SpecializedContainerExpression);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    SpecializedContainerExpression(ModelPart& rModelPart)
        : BaseType(rModelPart)
    {
    }

    /// Copy constructor with base class used to transfer data between compatible data containers
    SpecializedContainerExpression(const BaseType& rOther)
        : BaseType(rOther)
    {
    }

    SpecializedContainerExpression& operator=(const SpecializedContainerExpression& rOther);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clones the existing data container.
     *
     * This clones existing specialized data container. This is light weight operation
     * since this just clones the expression pointer. No data copying for the underlying
     * data in expression is done.
     *
     * @return SpecializedContainerExpression::Pointer
     */
    SpecializedContainerExpression::Pointer Clone() const;

    using BaseType::Read;

    using BaseType::MoveFrom;

    /**
     * @brief Read data from the entities of the container type in model part
     *
     * This creates a LiteralVectorExpression with a vecor of size number_of_entities * dimension_of_variable
     * And this vector is populated with the rVariable values from entities in the container of the model part.
     *
     * This is compatible with OpenMP and MPI.
     *
     * @tparam TDataType
     * @param rVariable         Variable of the entities in the container
     */
    template <class TDataType>
    void Read(
        const Variable<TDataType>& rVariable);

    using BaseType::Evaluate;

    /**
     * @brief Assigns data of the expression to the entities of the container type in model part
     *
     * This method evaluates the expression in this container object and assigns value to the
     * respective variable of the entities in the container. This method does not allocate additional
     * memory of a vector to store the evaluated result of the expression hence this is memory efficient.
     *
     * The given variable dimension and the dimension of the expression in the container data should match.
     * Otherwise, an error is thrown. This is an expensive operation since this evaluates the expression for all the entities and stores them
     * in their respective container.
     *
     * This is compatible with OpenMP and MPI.
     *
     * @tparam TDataType
     * @param rVariable         Variable the data to be written to.
     */
    template <class TDataType>
    void Evaluate(
        const Variable<TDataType>& rVariable);

    /**
     * @brief Set the Data For Container Variable to given value.
     *
     * This replaces the existing expression with LiteralExpression<double>
     * or LiteralArra3Expression. This is a light weight operation.
     *
     * @tparam TDataType
     * @param rValue            Value the container data is set to.
     */
    template <class TDataType>
    void SetData(
        const TDataType& rValue);

    /**
     * @brief Set the Data To Variable Zero Value
     *
     * This replaces the existing expression with LiteralExpression<double>
     * or LiteralArra3Expression. This is a light weight operation.
     *
     * @tparam TDataType
     * @param rVariable         Variable which is used to get the zero value to set to.
     */
    template <class TDataType>
    void SetZero(
        const Variable<TDataType>& rVariable);

    ///@}
    ///@name Operators
    ///@{

    SpecializedContainerExpression operator+(const SpecializedContainerExpression& rOther) const;

    SpecializedContainerExpression& operator+=(const SpecializedContainerExpression& rOther);

    SpecializedContainerExpression operator+(const double Value) const;

    SpecializedContainerExpression& operator+=(const double Value);

    SpecializedContainerExpression operator-(const SpecializedContainerExpression& rOther) const;

    SpecializedContainerExpression& operator-=(const SpecializedContainerExpression& rOther);

    SpecializedContainerExpression operator-(const double Value) const;

    SpecializedContainerExpression& operator-=(const double Value);

    SpecializedContainerExpression operator*(const SpecializedContainerExpression& rOther) const;

    SpecializedContainerExpression& operator*=(const SpecializedContainerExpression& rOther);

    SpecializedContainerExpression operator*(const double Value) const;

    SpecializedContainerExpression& operator*=(const double Value);

    SpecializedContainerExpression operator/(const SpecializedContainerExpression& rOther) const;

    SpecializedContainerExpression& operator/=(const SpecializedContainerExpression& rOther);

    SpecializedContainerExpression operator/(const double Value) const;

    SpecializedContainerExpression& operator/=(const double Value);

    SpecializedContainerExpression Pow(const SpecializedContainerExpression& rOther) const;

    SpecializedContainerExpression Pow(const double Value) const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}
};

///@}
/// output stream function
template <class TContainerType, class TContainerDataIO, class TMeshType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rThis)
{
    return rOStream << rThis.Info();
}

} // namespace Kratos

// includes the implementations of here
// This is kept headers only because, this can be customized
// in applications if required.
#include "specialized_container_expression_impl.h"
