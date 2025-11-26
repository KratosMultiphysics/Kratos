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

    template <class ContainerType, class IdsType>
    static void CollectIdsFromEntity(const ContainerType& rContainer, IdsType& rIds)
    {
        for (const auto& r_item : rContainer) {
            const auto id = r_item.Id();

            if constexpr (requires { rIds.insert(id); }) {
                // Works for set/unordered_set
                rIds.insert(id);
            } else if constexpr (requires { rIds.push_back(id); }) {
                // Works for vector/deque
                rIds.push_back(id);
            } else {
                static_assert(!std::is_same_v<IdsType, IdsType>,
                              "OutputContainerType does not support insert or push_back.");
            }
        }
    }

    template <class ContainerType>
    static auto CollectIdsFromEntity(const ContainerType& rContainer)
    {
        using ItemRef = decltype(*std::begin(rContainer));
        using Item    = std::remove_cv_t<std::remove_reference_t<ItemRef>>;
        using IdType  = decltype(std::declval<Item>().Id());
        std::vector<IdType> result;
        CollectIdsFromEntity(rContainer, result);
        return result;
    }

}; // class GenericUtilities

} // namespace Kratos
