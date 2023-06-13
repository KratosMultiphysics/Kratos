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
#include "expression/container_expression.h"
#include "expression/view_operators.h"

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
template <class TContainerType, class TContainerDataIO, MeshType TMeshType = MeshType::Local>
class  SpecializedContainerExpression : public ContainerExpression<TContainerType, TMeshType> {
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

    /**
     * @brief Returns a slice of the current expression. Slicing is based on the entity values.
     *
     * @details This method returns a new sliced @ref SpecializedContainerExpression.
     *          This slicing is done on each entitiy's data array, and not on the flattened
     *          expression. Following is an example:
     *
     *          Assume a SpecializedContainerExpression with an expression of shape [5] and 2 entities with
     *          following data values in the flattened representation.
     *
     *          data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
     *                  <---- 1 ----> <----- 2 ----->
     *
     *          Data for entity 1 is represented with <--1-->.
     *
     *          If the Offset is 1, and Stride is 3 then, the @ref SpecializedContainerExpression output
     *          of this method will produce the lazy expression which will give following data when SpecializedContainerExpression::Evaluate
     *          is called.
     *
     *          output_data = [2, 3, 4, 7, 8, 9]
     *          output containers shape = [3] = equal to Stride.
     *
     *          Slicing will always create one dimensional vector even if the input has more than one dimension.
     *          @see Reshape to reshape the one dimensional vector to desired shape if required.
     *
     *          This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
     *
     * @param Offset                            Offset of entity value to be considered to start slicing.
     * @param Stride                            Length of components from the offset in the entity value.
     * @return SpecializedContainerExpression   Container expresssion having the sliced lazy expression.
     */
    SpecializedContainerExpression Slice(
        const IndexType Offset,
        const IndexType Stride) const;

    /**
     * @brief Returns current expression reshaped to specified shape.
     *
     * @details This method returns a new reshaped @ref SpecializedContainerExpression.
     *          This reshaping is done on each entitiy's data array, and not on the flattened
     *          expression. Following is an example:
     *
     *          Assume a SpecializedContainerExpression with an expression of shape [2, 3] and 2 entities with
     *          following data values in the flattened representation.
     *
     *          data = [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [10, 11, 12]]]
     *                  <-------- 1 --------->  <----------- 2 ----------->
     *
     *          If the rShape = [3, 2] then the returned @ref SpecializedContainerExpression will have an lazy
     *          expression which will return following output data when Evaluate is called.
     *
     *          output_data = [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]
     *          output containers shape = [3, 2]
     *
     *          This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
     *
     * @param rShape                            New shape to used to reshape the existing expression.
     * @return SpecializedContainerExpression   New container expression with data reshaped to @ref rShape.
     */
    SpecializedContainerExpression Reshape(const std::vector<IndexType>& rShape) const;

    /**
     * @brief Returns current expression reshaped to specified shape.
     *
     * @details This method is same as the previous Reshape method. It accepts start and end iterators
     *          of the shape.
     *
     * @see Reshape(const IndexType Offset, const IndexType Stride)
     * @param Begin                             Start iterator of the shape iterator.
     * @param End                               End iterator of the shape iterator.
     * @return SpecializedContainerExpression   New container expression with data reshaped to @ref rShape.
     */
    template<class TIteratorType>
    SpecializedContainerExpression Reshape(
        TIteratorType Begin,
        TIteratorType End) const
    {
        SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
        result.mpExpression = Kratos::Reshape(*this->mpExpression, Begin, End);
        return result;
    }

    /**
     * @brief Returns container expression which combines current and all the other container expressions provided.
     *
     * @details     This method provides a combined container expression as explained in the following example.
     *              All the given container expressions in @ref rListOfOthers should have the same number of
     *              entities as in the current container expression. The combination is done in the following order
     *              (The values won't be evaluated unless Evaluate is called.).
     *                  1. First entity values of current expression put to the entity values of the output expression.
     *                  2. Then entity values of each expression of @ref rListOfOthers is put one after other.
     *
     *              Assume current exression has following data with item shape [2] fand with 2 entities
     *
     *              data = [1, 2, 3, 4]
     *
     *              If the @ref rListOfOthers has following Container expressions
     *                  rListOfOthers[0] = data{5, 6} with 2 entities, and shape = []
     *                  rListOfOthers[1] = data{7, 8, 9, 10} with 2 entities, and shape = [2]
     *
     *              Then the resulting container expression will have shape of [5] with 2 entities.
     *              Following is the output data array which can be obtained after calling Evaluate method
     *
     *                  output_data = [1, 2, 5, 7, 8, 3, 4, 6, 9, 10]
     *
     *              This creates a lazy expression, hence it has a constant cost complexity irrespective of the data size.
     *
     * @param rListOfOthers                         List of other container expressions to be combined with the existing one.
     * @return SpecializedContainerExpression       New container expression with data combined.
     */
    SpecializedContainerExpression Comb(const std::vector<typename BaseType::Pointer>& rListOfOthers) const;

    /**
     * @brief Returns container expression which combines current rOther expression.
     *
     * @details @see Com(const std::vector<typename BaseType::Pointer>& rListOfOthers) method for details.
     *
     * @param rOther                                Other contaienr expression to be combined.
     * @return SpecializedContainerExpression       New container expression with data combined.
     */
    SpecializedContainerExpression Comb(const BaseType& rOther) const;

    /**
     * @brief Returns container expression which combines current and other expressions given by Begin and end.
     *
     * @param Begin                             Begining of the container expressions list.
     * @param End                               End of the container expressions list.
     * @return SpecializedContainerExpression       New container expression with data combined.
     */
    template<class TIteratorType>
    SpecializedContainerExpression Comb(
        TIteratorType Begin,
        TIteratorType End) const
    {
        SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
        std::vector<Expression::ConstPointer> expressions;
        expressions.push_back(this->pGetExpression());
        for (auto itr = Begin; itr != End; ++itr) {
            expressions.push_back((*itr)->pGetExpression());
        }
        result.mpExpression = Kratos::Comb(expressions.begin(), expressions.end());
        return result;
    }

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

    SpecializedContainerExpression Power(const SpecializedContainerExpression& rOther) const;

    SpecializedContainerExpression Power(const double Value) const;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}
};

///@}
/// output stream function
template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
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
