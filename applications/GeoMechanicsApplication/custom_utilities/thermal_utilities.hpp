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
//                   John van Esch
//                   

#pragma once

// System includes
#include "includes/node.h"
#include "geometries/geometry.h"

// External includes

// Project includes

#pragma once 

namespace Kratos {
class ThermalUtilities {

public:

    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;

    // ============================================================================================
    // ============================================================================================
    static double CalculateWaterDensityOnIntegrationPoints(
        const Vector& N,
    	const GeometryType& rGeom)
    {
        const unsigned int TNumNodes = rGeom.PointsNumber();
        Vector TemperatureVector;
        TemperatureVector.resize(TNumNodes);
        //Nodal Variables
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            TemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
        }
    
        double temp = inner_prod(N, TemperatureVector);
        double density = +9.998396e+2 + 6.764771e-2 * temp - 8.993699e-3 * std::pow(temp, 2)
                         + 9.143518e-5 * std::pow(temp, 3) - 8.907391e-7 * std::pow(temp, 4)
                         + 5.291959e-9 * std::pow(temp, 5) - 1.359813e-11 * std::pow(temp, 6);
        return density;
    }
    
    // ============================================================================================
    // ============================================================================================
    static double CalculateWaterViscosityOnIntegrationPoints(
        const Vector& N,
    	const GeometryType& rGeom)
    {
        const unsigned int TNumNodes = rGeom.PointsNumber();
        Vector TemperatureVector;
        TemperatureVector.resize(TNumNodes);
        //Nodal Variables
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            TemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
        }
    
        double temp = inner_prod(N, TemperatureVector);
        double c1 = 247.8 / (temp + 133.0);
        double viscosity = 2.4318e-5 * std::pow(10.0, c1);
        return viscosity;
    }

}; // class ThermalUtilities
}

