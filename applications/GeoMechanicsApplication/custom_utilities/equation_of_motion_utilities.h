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
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

class GeoEquationOfMotionUtilities
{
public:
    static Matrix CalculateMassMatrix(SizeType                   dimension,
                                      SizeType                   number_U_nodes,
                                      SizeType                   NumberIntegrationPoints,
                                      const Matrix&              Nu_container,
                                      const Vector&              rSolidDensities,
                                      const std::vector<double>& rIntegrationCoefficients);

    static Vector CalculateDetJsInitialConfiguration(const Geometry<Node>& rGeom,
                                                     const GeometryData::IntegrationMethod IntegrationMethod);

}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
