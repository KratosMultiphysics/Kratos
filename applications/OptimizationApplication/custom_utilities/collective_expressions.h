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
#include "expression/specialized_container_expression.h"
#include "expression/container_data_io.h"

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

    using VariableTypes = std::variant<
            const Variable<double>*,
            const Variable<array_1d<double, 3>>*,
            const Variable<array_1d<double, 4>>*,
            const Variable<array_1d<double, 6>>*,
            const Variable<array_1d<double, 9>>*,
            const Variable<Vector>*,
            const Variable<Matrix>*>;

    KRATOS_CLASS_POINTER_DEFINITION(CollectiveExpressions);

    ///@}
    ///@name Life cycle
    ///#{

    // Default constructor
    CollectiveExpressions() noexcept = default;


    /**
     * @brief Construct a list of container expressions
     *
     * @param rExpressionPointersList       List of container expressions
     */
    CollectiveExpressions(const std::vector<CollectiveExpressionType>& rExpressionPointersList);

    // Copy consttructor
    CollectiveExpressions(const CollectiveExpressions& rOther);

    // destructor
    ~CollectiveExpressions() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clones this collective expressions by cloning all the ContainerExpressions in the instance.
     *
     * This clones current instance of CollectiveExpressions by cloning all the ContainerExpressions. Cloning
     * a single ContainerExpression is a light weight operation, hence this operation also a light weight operation.
     *
     * @return CollectiveExpressions        Cloned collective expressions which has all the clones of ContainerExpressions.
     */
    CollectiveExpressions Clone() const;

    /**
     * @brief Set the To all the ContainerExpression objects.
     *
     */
    void SetToZero();

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
    void Add(const CollectiveExpressions& rCollectiveExpression);

    /**
     * @brief Clear list of all the ContainerExpressions.
     *
     * This does not destroy the corresponding ContainerExpression (uses SmartPointers). This removes
     * all the ContainerExpression objects from the current list.
     *
     */
    void Clear();

    /**
     * @brief Read and copy data from pBegin array to internal structure.
     *
     * This method read and copy data from @ref pBegin to internal structures of containing ContainerExpressions. The shape of the each
     * data set for each ContainerExpression is provided by the flattened matrix alike in pListShapeBegin. This matrix alike will have different
     * columns in each row. Each row corresponds to the ContainerExpression index in the list, and columns in each row represents
     * number of dimensions for that ContainerExpression. The value in column represents the size of each dimension. @ref NumberOfEntities
     * vector represents how much data for entities is present for each ContainerExpression in @ref pBegin.
     *
     * @param pBegin                Starting pointer to the data to read from.
     * @param NumberOfEntities      Number of entities for each container.
     * @param pListShapeBegin       Matrix of shape sizes. Rows for each container, columns for size of each dimension.
     * @param ShapeSizes            Vector containing number of dimensions for each container.
     * @param NumberOfContainers    Number of containers.
     */
    void Read(
        double const* pBegin,
        int const* NumberOfEntities,
        int const** pListShapeBegin,
        int const* ShapeSizes,
        const int NumberOfContainers);

    /**
     * @brief Read data from the ContaineExpressions
     *
     * This method reads data from each ContainerExpression corresponding to @ref rVariable.
     *
     * @param rVariable     Variable to be used to read data from each container in ContainerExpressions.
     */
    void Read(const VariableTypes& rVariable);

    /**
     * @brief Read data from the ContaineExpressions
     *
     * This method reads data from each ContainerExpression corresponding to each variable in @ref rVariables list.
     *
     * @param rVariables     Variables to be used to read data from each container in ContainerExpressions.
     */
    void Read(const std::vector<VariableTypes>& rVariables);

    /**
     * @brief Move data from pBegin array to internal structure.
     *
     * This method moves data from @ref pBegin to internal structures of containing ContainerExpressions. The shape of the each
     * data set for each ContainerExpression is provided by the flattened matrix alike in pListShapeBegin. This matrix alike will have different
     * columns in each row. Each row corresponds to the ContainerExpression index in the list, and columns in each row represents
     * number of dimensions for that ContainerExpression. The value in column represents the size of each dimension. @ref NumberOfEntities
     * vector represents how much data for entities is present for each ContainerExpression in @ref pBegin.
     *
     * @warning This instance does not take the ownership of the passed array. Hence, life time is not managed by this instance.
     *
     * @param pBegin                Starting pointer to the data to read from.
     * @param NumberOfEntities      Number of entities for each container.
     * @param pListShapeBegin       Matrix of shape sizes. Rows for each container, columns for size of each dimension.
     * @param ShapeSizes            Vector containing number of dimensions for each container.
     * @param NumberOfContainers    Number of containers.
     */
    void MoveFrom(
        double* pBegin,
        int const* NumberOfEntities,
        int const** pListShapeBegin,
        int const* ShapeSizes,
        const int NumberOfContainers);

    /**
     * @brief Assign the data in the expression to c-like interfaces.
     *
     * This method can assign data in all the ContainerExpressions to a c-interface where @ref pBegin is the
     * starting pointer to a contiguous array having @ref Size elements where @ref Size corresponds
     * to total number of double values present in each ContainerExpression.
     *
     * @param pBegin            Starting pointer to the data.
     * @param Size              Size of the raw double data vector.
     */
    void Evaluate(
        double* pBegin,
        const int Size) const;

    /**
     * @brief Assigns evaluated expressions to the rVariable in all the ContainerExpressions.
     *
     * @param rVariable     Variable to be used to assign evaluated data in each ContainerExpression.
     */
    void Evaluate(const VariableTypes& rVariable);

    /**
     * @brief Assigns evaluated expressions to repective variables in all the ContainerExpressions.
     *
     * @param rVariables    List of variables each corresponding to a variable which is used in assignment.
     */
    void Evaluate(const std::vector<VariableTypes>& rVariables);

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