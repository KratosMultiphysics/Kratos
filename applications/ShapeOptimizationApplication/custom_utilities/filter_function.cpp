// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Aditya Ghantasala, https://github.com/adityaghantasala
//
// ==============================================================================


// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "filter_function.h"

// ==============================================================================

namespace Kratos
{

FilterFunction::FilterFunction(const std::string FilterFunctionType)
{
    // Set type of weighting function

    // Type 1: Gaussian function
    if (FilterFunctionType == "gaussian")
        mFilterFunctional =  [](double radius, double distance) {return std::max(0.0, exp(-(distance*distance) / (2 * radius * radius / 9.0)));};

    // Type 2: Linear function
    else if (FilterFunctionType == "linear")
        mFilterFunctional =  [](double radius, double distance) {return std::max(0.0, (radius - distance) / radius);};

    // Type 3: Constant function
    else if (FilterFunctionType == "constant")
        mFilterFunctional = [](double radius, double distance) {return 1.0;};

    // Type 4: Cosine function
    else if (FilterFunctionType == "cosine")
        mFilterFunctional = [](double radius, double distance) {return std::max(0.0, 1-0.5*(1-std::cos(Globals::Pi/radius*distance)));};

    // Type 5: Quartic function
    else if (FilterFunctionType == "quartic")
        mFilterFunctional = [](double radius, double distance) {return std::max(0.0, (pow(distance-radius,4.0)/pow(radius,4.0)));};

    // Type 6: Green's function
    else if (FilterFunctionType == "green")
        mFilterFunctional = [](double radius, double distance) {return (1.0/((4*Kratos::Globals::Pi*distance)/(radius*radius)+1)) * exp(-distance/radius);};  

    // Throw error message in case of wrong specification
    else
        KRATOS_ERROR << "Specified kernel function of type : "<< FilterFunctionType << " is not recognized. \n \t Options are: constant, linear , gaussian, cosine, quartic and green." << std::endl;
}

double FilterFunction::ComputeWeight(const Array3DType& ICoord, const Array3DType& JCoord, const double Radius) const
{
    KRATOS_TRY;

    // Compute distance vector
    const double distance = GetDistance(ICoord, JCoord);

    // Depending on which weighting function is chosen, compute weight
    return mFilterFunctional(Radius, distance);

    KRATOS_CATCH("");
}

} // namespace Kratos.
