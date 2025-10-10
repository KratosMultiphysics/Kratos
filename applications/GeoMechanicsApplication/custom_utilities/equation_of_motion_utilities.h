// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

// Project includes

// Application includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoEquationOfMotionUtilities
{
public:
    static Matrix CalculateMassMatrix(std::size_t                dimension,
                                      std::size_t                number_U_nodes,
                                      std::size_t                NumberIntegrationPoints,
                                      const Matrix&              Nu_container,
                                      const std::vector<double>& rSolidDensities,
                                      const std::vector<double>& rIntegrationCoefficients);

    static Vector CalculateDetJsInitialConfiguration(const Geometry<Node>& rGeom,
                                                     const GeometryData::IntegrationMethod IntegrationMethod);

    static Matrix CalculateDampingMatrix(double        RayleighAlpha,
                                         double        RayleighBeta,
                                         const Matrix& rMassMatrix,
                                         const Matrix& rStiffnessMatrix);

    template<typename TMatrix>
    static void CalculateStiffnessMatrixGPoint(TMatrix&       rStiffness,
                                               const Matrix& rB,
                                               const Matrix& rConstitutiveMatrix,
                                               double        IntegrationCoefficient)
    {
        const std::size_t strain_size1   = rConstitutiveMatrix.size1(); // rows
        const std::size_t number_of_dofs = rB.size2();

        Matrix CB(strain_size1, number_of_dofs);
        noalias(CB) = prod(rConstitutiveMatrix, rB);

        Matrix transposed_B(number_of_dofs, strain_size1);
        noalias(transposed_B) = trans(rB);

        noalias(rStiffness) = prod(transposed_B, CB);

        rStiffness *= IntegrationCoefficient;
    }

    static Matrix CalculateStiffnessMatrix(const std::vector<Matrix>& rBs,
                                           const std::vector<Matrix>& rConstitutiveMatrices,
                                           const std::vector<double>& rIntegrationCoefficients);

    static Vector CalculateInternalForceVector(const std::vector<Matrix>& rBs,
                                               const std::vector<Vector>& rStressVectors,
                                               const std::vector<double>& rIntegrationCoefficients);

}; /* Class GeoEquationOfMotionUtilities*/
} /* namespace Kratos.*/
