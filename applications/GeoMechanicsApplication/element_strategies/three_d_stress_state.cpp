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
#include "three_d_stress_state.h"
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

void ThreeDStressState::CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector& Np, const Geometry<Node>& rGeometry)
{
    const SizeType dimension       = rGeometry.WorkingSpaceDimension();
    const SizeType number_of_nodes = rGeometry.size();

    KRATOS_ERROR_IF(dimension != 3) << "ThreeDStressState::CalculateBMatrix called for a "
                                    << dimension << "D element." << std::endl;

    SizeType index;
    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        index = dimension * i;

        rB(INDEX_3D_XX, index + INDEX_X) = GradNpT(i, INDEX_X);
        rB(INDEX_3D_YY, index + INDEX_Y) = GradNpT(i, INDEX_Y);
        rB(INDEX_3D_ZZ, index + INDEX_Z) = GradNpT(i, INDEX_Z);
        rB(INDEX_3D_XY, index + INDEX_X) = GradNpT(i, INDEX_Y);
        rB(INDEX_3D_XY, index + INDEX_Y) = GradNpT(i, INDEX_X);
        rB(INDEX_3D_YZ, index + INDEX_Y) = GradNpT(i, INDEX_Z);
        rB(INDEX_3D_YZ, index + INDEX_Z) = GradNpT(i, INDEX_Y);
        rB(INDEX_3D_XZ, index + INDEX_X) = GradNpT(i, INDEX_Z);
        rB(INDEX_3D_XZ, index + INDEX_Z) = GradNpT(i, INDEX_X);
    }
}

double ThreeDStressState::CalculateIntegrationCoefficient(Geometry<Node>::IntegrationPointsArrayType& IntegrationPoints,
                                                        unsigned int PointNumber,
                                                        double detJ,
                                                        const Geometry<Node>& rGeometry)
{
    return IntegrationPoints[PointNumber].Weight() * detJ;
}
} // namespace Kratos
