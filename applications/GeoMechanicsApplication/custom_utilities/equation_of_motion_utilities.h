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

class GeoEquationOfMotionUtilities
{
public:
    static Matrix CalculateMassMatrix(std::size_t                dimension,
                                      std::size_t                number_U_nodes,
                                      std::size_t                NumberIntegrationPoints,
                                      const Matrix&              Nu_container,
                                      const Vector&              rSolidDensities,
                                      const std::vector<double>& rIntegrationCoefficients);

    static Vector CalculateDetJsInitialConfiguration(const Geometry<Node>& rGeom,
                                                     const GeometryData::IntegrationMethod IntegrationMethod);

    static Matrix CalculateDampingMatrix(double        RayleighAlpha,
                                         double        RayleighBeta,
                                         const Matrix& rMassMatrix,
                                         const Matrix& rStiffnessMatrix);

}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
