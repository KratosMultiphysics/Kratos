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

/**
 * @brief Construct a new specialized container variable data
 *
 * This class is used to represent objects which are the specializations
 * of ContainerVariableData. The base class ifof type ContainerVariableData<TContainerType>.
 * Hence this class can be specialized to read and write data to different entity containers
 * for each container type.
 *
 * This class does not have any member variables, hence it uses expressions to represent the
 * data. Therefore, all the operations such as "+", "-", "*", etc are light weight because
 * they do not compute the value in place, they just create an expression which keeps track fo the expressions.
 * These expressions are evaluated only if AssignData. This also does not create additional
 * vectors to hold the resultant value of the expression. It uses the model parts respective containers entity input/output
 * method specified to write evaluated entity resultant values to model part entities. Hence, these
 * SpecializedContainerVariableData can be easily visualized using common variable data visualization methods.
 *
 * Copy constructor is introduced with the base class type becaues, this allows copying data
 * between compatible SpecializedContainerVariableData containers. This copy is also light weight since
 * this only copies the pointers, not the data itself.
 *
 * This class is optimized and compatible with OpenMP and MPI.
 *
 * @tparam TContainerType           Container type
 * @tparam TContainerDataIO         Container entity input/output type
 */
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

    /**
     * @brief Clones the existing data container.
     *
     * This clones existing specialized data container. This is light weight operation
     * since this just clones the expression pointer. No data copying for the underlying
     * data in expression is done.
     *
     * @return SpecializedContainerVariableData::Pointer
     */
    SpecializedContainerVariableData::Pointer Clone() const;

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
    void ReadData(
        const Variable<TDataType>& rVariable);

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
    void AssignData(
        const Variable<TDataType>& rVariable);

    /**
     * @brief Set the Data For Container Variable to given value.
     *
     * This replaces the existing expression with LiteralDoubleExpression
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
     * This replaces the existing expression with LiteralDoubleExpression
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
