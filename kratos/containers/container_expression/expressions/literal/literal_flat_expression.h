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
#include <vector>
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/container_expression/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Expression to hold a literal with a flattened data structure
 *
 * This expression holds data of a literal with rShape by
 * flattening the data structure to a vector. This allocates
 * a Vector with large sizes.
 *
 */
class KRATOS_API(KRATOS_CORE) LiteralFlatExpression : public Expression {
public:
    ///@name Type definition
    ///@{

    using Pointer = Kratos::intrusive_ptr<LiteralFlatExpression>;

    ///@}
    ///@name Life cycle
    ///@{

    LiteralFlatExpression(
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    LiteralFlatExpression(
        double* pDataBegin,
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    LiteralFlatExpression(const LiteralFlatExpression& rOther) = delete;

    ~LiteralFlatExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief A specific create method is used in here to return a pointer to the LiteralFlatExpression
     *
     * This method is used to return the derrived class LiteralFlatExpression pointer
     * so that additional SetData method can be accessed to initialize
     * the vector data.
     *
     * @param NumberOfEntities                  Number of entities in the container.
     * @param rShape                            Shape of the data in each entitiy.
     * @return LiteralFlatExpression::Pointer   Returns an intrusive pointer to LiteralFlatExpression.
     */
    static LiteralFlatExpression::Pointer Create(
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    static LiteralFlatExpression::Pointer Create(
        double* pDataBegin,
        const IndexType NumberOfEntities,
        const std::vector<IndexType>& rShape);

    void SetData(
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex,
        const double Value);

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

    ///@}
protected:
    ///@name Protected classes
    ///@{

    class Data
    {
    public:
        ///@}
        ///@name Life cycle
        ///@{

        /**
         * @brief Construct a new Data object where the array memory is managed by the object.
         *
         * @param Size      Size of the allocated array.
         */
        Data(const IndexType Size): mpBegin(new double[Size]), mIsManaged(true) {}

        /**
         * @brief Construct a new Data object, where the underlying array memory is not managed by the object.
         *
         * @param pBegin    Pointer to the memory array.
         */
        Data(double* pBegin): mpBegin(pBegin), mIsManaged(false) {}

        ~Data() { if (mIsManaged) { delete[] mpBegin; } }

        ///@}
        ///@name Operators
        ///@{

        double& operator[](const IndexType Index) noexcept { return mpBegin[Index]; }

        double operator[](const IndexType Index) const noexcept { return mpBegin[Index]; }

        ///@}
    private:
        ///@name Private member variables
        ///@{

        double* mpBegin;

        bool mIsManaged;

        ///@}
    };

    ///@}
    ///@name Protected member variables
    ///@{

    const std::vector<IndexType> mShape;

    Data mData;

    ///@}
};

class LiteralScalarFlatExpression : public LiteralFlatExpression
{
public:

    using LiteralFlatExpression::LiteralFlatExpression;

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

};

class LiteralNonScalarFlatExpression : public LiteralFlatExpression
{
public:

    using LiteralFlatExpression::LiteralFlatExpression;

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

};

} // namespace Kratos