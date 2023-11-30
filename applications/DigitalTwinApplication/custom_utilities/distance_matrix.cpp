//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes

// Application includes

// Include base h
#include "distance_matrix.h"

namespace Kratos {

DistanceMatrix::DistanceMatrix(const std::vector<double>& rCompressedDistancesMatrix)
    : mCompressedDistanceMatrix(rCompressedDistancesMatrix)
{
    mNumberOfItems = ((1.0 + std::sqrt(8.0*mCompressedDistanceMatrix.size() + 1.0)) / 2);

    KRATOS_ERROR_IF_NOT(mNumberOfItems * (mNumberOfItems - 1) == 2 * mCompressedDistanceMatrix.size())
        << "Invalid size in compressed distance matrix [ rCompressedDistanceMatrix.size() = "
        << rCompressedDistancesMatrix.size() << ", derived number of entities = " << mNumberOfItems << " ].\n";
}

double DistanceMatrix::GetDistance(
    const IndexType Index1,
    const IndexType Index2) const
{
    IndexType local_1{Index1}, local_2{Index2};
    if (Index1 > Index2) {
        local_1 = Index2;
        local_2 = Index1;
    }

    return mCompressedDistanceMatrix[mNumberOfItems * local_1 + local_2 - static_cast<IndexType>((local_1 + 2) * (local_1 + 1) / 2)];
}

} /* namespace Kratos.*/