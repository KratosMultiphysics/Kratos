//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mohamed Nabi
//                   



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
    template<unsigned int TDim, unsigned int TNumNodes>
    static double CalculateWaterDensityOnIntegrationPoints(
        const Vector& N,
    	const GeometryType& rGeom)
    {
        array_1d<double, TNumNodes> TemperatureVector;
        //Nodal Variables
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            TemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
        }
    
        double temp = inner_prod(N, TemperatureVector);
        return
            + 9.998396e+2 + 6.764771e-2 * temp - 8.993699e-3 * std::pow(temp, 2)
            + 9.143518e-5 * std::pow(temp, 3) - 8.907391e-7 * std::pow(temp, 4)
            + 5.291959e-9 * std::pow(temp, 5) - 1.359813e-11 * std::pow(temp, 6);
    }
    
    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    static double CalculateWaterViscosityOnIntegrationPoints(
        const Vector& N,
    	const GeometryType& rGeom)
    {
        array_1d<double, TNumNodes> TemperatureVector;
        //Nodal Variables
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            TemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
        }
    
        double temp = inner_prod(N, TemperatureVector);
        double c1 = 247.8 / (temp + 133.0);
        return 1.0 / (2.4318e-5 * std::pow(10.0, c1));
    }
	
	

}; // class StructuralMechanicsElementUtilities.
}

