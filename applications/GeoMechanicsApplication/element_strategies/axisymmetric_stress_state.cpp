// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "axisymmetric_stress_state.h"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

void AxisymmetricStressState::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector& Np, const Geometry<Node>& rGeometry)
{
    KRATOS_TRY

    const double radius = GeoElementUtilities::CalculateRadius(Np, rGeometry);

    const SizeType dimension = rGeometry.WorkingSpaceDimension();
    const SizeType number_of_nodes = rGeometry.size();

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const IndexType index = dimension * i;

        rB(INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X) = GradNpT(i, INDEX_X);
        rB(INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y) = GradNpT(i, INDEX_Y);
        rB(INDEX_2D_PLANE_STRAIN_ZZ, index + INDEX_X) = Np[i] / radius;
        rB(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X) = GradNpT(i, INDEX_Y);
        rB(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y) = GradNpT(i, INDEX_X);
    }

    KRATOS_CATCH("")

}

} // namespace Kratos
