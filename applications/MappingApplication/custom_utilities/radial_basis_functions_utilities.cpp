//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:    Juan I. Camarotti

// System includes

// External includes

// Project includes
#include "radial_basis_functions_utilities.h"

namespace Kratos
{
namespace RadialBasisFunctionsUtilities
{
    // Extracted from: Review of coupling methods for non-matching meshes (https://doi.org/10.1016/j.cma.2006.03.017)
    double CalculateWendlandC2SupportRadius(const Matrix& rPoints, const double k = 2.5)
    {
        const SizeType n_points = rPoints.size1();
        KRATOS_ERROR_IF(n_points < 2) << "At least two points are required to estimate spacing." << std::endl;

        double total_distance = 0.0;
        SizeType count = 0;

        for (SizeType i = 0; i < n_points; ++i) {
            for (SizeType j = i + 1; j < n_points; ++j) {
                const double distance = norm_2(row(rPoints, i) - row(rPoints, j));
                total_distance += distance;
                ++count;
            }
        }

        const double average_spacing = total_distance / count;
        return k * average_spacing; // Support radius for Wendland C2
    }

} // namespace RadialBasisFunctionsUtilities
} // namespace Kratos.