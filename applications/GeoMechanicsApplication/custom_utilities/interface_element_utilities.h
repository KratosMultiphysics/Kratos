// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "containers/array_1d.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) InterfaceElementUtilities
{
public:
    static void CalculateNuMatrix(BoundedMatrix<double, 2, 4>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint);
    static void CalculateNuMatrix(BoundedMatrix<double, 2, 8>& rNu, const Matrix& Ncontainer, const unsigned int& GPoint);
    static void CalculateNuMatrix(BoundedMatrix<double, 3, 12>& rNu,
                                  const Matrix&                 Ncontainer,
                                  const unsigned int&           GPoint);
    static void CalculateNuMatrix(BoundedMatrix<double, 3, 18>& rNu,
                                  const Matrix&                 Ncontainer,
                                  const unsigned int&           GPoint);
    static void CalculateNuMatrix(BoundedMatrix<double, 3, 24>& rNu,
                                  const Matrix&                 Ncontainer,
                                  const unsigned int&           GPoint);

    static void FillPermeabilityMatrix(BoundedMatrix<double, 2, 2>& rPermeabilityMatrix,
                                       const double&                JointWidth,
                                       const double&                Transversal_Permeability);
    static void FillPermeabilityMatrix(BoundedMatrix<double, 3, 3>& rPermeabilityMatrix,
                                       const double&                JointWidth,
                                       const double&                Transversal_Permeability);

    static void CalculateVoigtVector(array_1d<double, 2>& rVoigtVector);
    static void CalculateVoigtVector(array_1d<double, 3>& rVoigtVector);

    static void CalculateLinkPermeabilityMatrix(BoundedMatrix<double, 2, 2>& rPermeabilityMatrix,
                                                const double&                JointWidth);
    static void CalculateLinkPermeabilityMatrix(BoundedMatrix<double, 3, 3>& rPermeabilityMatrix,
                                                const double&                JointWidth);

    static Matrix Calculate2DRotationMatrix(const Geometry<Node>&      rGeometry,
                                            const array_1d<double, 3>& rLocalCoordinate);

private:
    static std::vector<array_1d<double, 3>> CalculateMidPoints(const Geometry<Node>& rGeometry);
    static array_1d<double, 3> CalculateTangentialVector(const Matrix& rShapeFunctionGradients,
                                                         const std::vector<array_1d<double, 3>>& rLocations);
};

} // namespace Kratos

