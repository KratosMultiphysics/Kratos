// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Gennady Markelov
//

#pragma once

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <algorithm>
#include <cstdlib>
#include <vector>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GenericUtilities
{
public:
    // Implementation based on Raymond Chen's article "Applying a permutation to a vector, part 1"
    // (see https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095)
    template <typename VectorType, typename IndexSequenceType>
    static VectorType PermutedVector(const VectorType& rVector, const IndexSequenceType& rIndices)
    {
        auto result = VectorType(rVector.size());
        std::ranges::transform(rIndices, result.begin(),
                               [&rVector](auto Index) { return rVector[Index]; });
        return result;
    }

    template <typename MatrixType, typename IndexSequenceType>
    static MatrixType MatrixWithPermutedColumns(const MatrixType& rMatrix, const IndexSequenceType& rIndices)
    {
        auto result = MatrixType(rMatrix.size1(), rMatrix.size2());
        for (auto i = std::size_t{0}; i < rMatrix.size1(); ++i) {
            for (auto j = std::size_t{0}; j < rMatrix.size2(); ++j) {
                result(i, j) = rMatrix(i, rIndices[j]);
            }
        }
        return result;
    }

    template <class ContainerType, class OutputIt>
    static void GetIdsFromEntityContents(const ContainerType& rContainer, OutputIt Out)
    {
        std::ranges::transform(rContainer, Out, [](const auto& rItem) { return rItem.Id(); });
    }

    template <class ContainerType>
    static auto GetIdsFromEntityContents(const ContainerType& rContainer)
    {
        std::vector<std::size_t> result;
        result.reserve(rContainer.size());
        GetIdsFromEntityContents(rContainer, std::back_inserter(result));
        return result;
    }

}; // class GenericUtilities

} // namespace Kratos
