// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include "includes/ublas_interface.h"

namespace Kratos
{

class UBlasUtils
{
public:
    static Vector MakeVector(const std::initializer_list<double>& values);

    template <typename InputIt>
    static Matrix MakeDiagonalMatrix(InputIt BeginOfDiagonalEntries,
                                     InputIt EndOfDiagonalEntries)
    {
        const auto size = static_cast<std::size_t>(EndOfDiagonalEntries - BeginOfDiagonalEntries);
        Matrix result{size, size, 0.0};
        for (auto i = std::size_t{0}; i < size; ++i) {
            result(i, i) = *BeginOfDiagonalEntries;
            ++BeginOfDiagonalEntries;
        }
        return result;
    }
};

}
