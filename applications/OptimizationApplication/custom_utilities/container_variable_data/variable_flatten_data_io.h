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
#include <vector>

// Project includes
#include "includes/define.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

// Dummy generic class to be specialized later
template <class TDataType>
class VariableFlattenDataIO {};

template<>
class VariableFlattenDataIO<double>
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using TDataType = double;

    ///@}
    ///@name Life cycle
    ///@{

    VariableFlattenDataIO(
        const TDataType SampleValue,
        const IndexType DomainSize)
    {
    }

    VariableFlattenDataIO(const std::vector<IndexType>& rShape)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    void Assign(
        TDataType& rOutput,
        double const* pValueBegin) const
    {
        rOutput = *pValueBegin;
    }

    void Read(
        double* pValueBegin,
        const TDataType Value) const
    {
        *pValueBegin = Value;
    }

    const std::vector<IndexType> GetShape() const
    {
        return {1};
    }

    ///@}
};

template<>
class VariableFlattenDataIO<array_1d<double, 3>>
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using TDataType = array_1d<double, 3>;

    ///@}
    ///@name Life cycle
    ///@{

    VariableFlattenDataIO(
        const TDataType& SampleValue,
        const IndexType DomainSize)
        : mLocalSize(DomainSize)
    {
        KRATOS_ERROR_IF_NOT(mLocalSize == 2 || mLocalSize == 3)
            << "Only supports 2 or 3 local size. [ local size = " << mLocalSize
            << " ].\n";
    }

    VariableFlattenDataIO(const std::vector<IndexType>& rShape)
    {
        KRATOS_ERROR_IF_NOT(rShape.size() == 1)
            << "Trying to use variable flatten data io with unsupported data\n";

        mLocalSize = rShape[0];

        KRATOS_ERROR_IF_NOT(mLocalSize == 2 || mLocalSize == 3)
            << "Only supports 2 or 3 local size. [ local size = " << mLocalSize
            << " ].\n";
    }

    ///@}
    ///@name Public operations
    ///@{

    void Assign(
        TDataType& rOutput,
        double const* pValueBegin) const
    {
        rOutput[0] = *(pValueBegin++);
        rOutput[1] = *(pValueBegin++);
        if (mLocalSize == 3) rOutput[2] = *(pValueBegin++);

    }

    void Read(
        double* pValueBegin,
        const TDataType& Value) const
    {
        *(pValueBegin++) = Value[0];
        *(pValueBegin++) = Value[1];
        if (mLocalSize == 3) *(pValueBegin++) = Value[2];
    }

    const std::vector<IndexType> GetShape() const
    {
        return {mLocalSize};
    }

    ///@}
private:
    ///@name Private members
    ///@{

    IndexType mLocalSize;

    ///@}
};

} // namespace Kratos