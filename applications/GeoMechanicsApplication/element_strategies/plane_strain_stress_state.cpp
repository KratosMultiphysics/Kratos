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
#include "plane_strain_stress_state.h"
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

void PlaneStrainStressState::CalculateBMatrix(Matrix& rB,
                                              const Matrix& GradNpT,
                                              const Vector& Np,
                                              const Geometry<Node>& rGeometry)
{
    const SizeType dimension       = rGeometry.WorkingSpaceDimension();
    const SizeType number_of_nodes = rGeometry.size();

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const std::size_t index = dimension * i;

        rB(INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X) = GradNpT(i, INDEX_X);
        rB(INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y) = GradNpT(i, INDEX_Y);
        rB(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X) = GradNpT(i, INDEX_Y);
        rB(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y) = GradNpT(i, INDEX_X);
    }
}

double PlaneStrainStressState::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointsArrayType& IntegrationPoints,
                                                             unsigned int PointNumber,
                                                             double detJ,
                                                             const Geometry<Node>& rGeometry)
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}

} // namespace Kratos
