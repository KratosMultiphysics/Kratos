//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel
//                   Aditya Ghantasala
//                   Suneth Warnakulasuriya
//

// System includes
#include <string>

// Project includes
#include "mapping_filter_functions.h"

namespace Kratos
{

MappingFilterFunction::MappingFilterFunction(const std::string& rFilterFunctionType)
{
    // Set type of weighting function

    // Type 1: Gaussian function
    if (rFilterFunctionType == "gaussian")
        mMappingFilterFunctional =  [](double radius, double distance) {return std::max(0.0, exp(-(distance*distance) / (2 * radius * radius / 9.0)));};

    // Type 2: Linear function
    else if (rFilterFunctionType == "linear")
        mMappingFilterFunctional =  [](double radius, double distance) {return std::max(0.0, (radius - distance) / radius);};

    // Type 3: Constant function
    else if (rFilterFunctionType == "constant")
        mMappingFilterFunctional = [](double radius, double distance) {return 1.0;};

    // Type 4: Cosine function
    else if (rFilterFunctionType == "cosine")
        mMappingFilterFunctional = [](double radius, double distance) {return std::max(0.0, 1-0.5*(1-std::cos(Globals::Pi/radius*distance)));};

    // Type 5: Quartic function
    else if (rFilterFunctionType == "quartic")
        mMappingFilterFunctional = [](double radius, double distance) {return std::max(0.0, (pow(distance-radius,4.0)/pow(radius,4.0)));};

    // Throw error message in case of wrong specification
    else
        KRATOS_ERROR << "Specified kernel function of type : "
                     << rFilterFunctionType << " is not recognized. \n \t Options are:"
                     << "\n\tconstant"
                     << "\n\tlinear"
                     << "\n\tgaussian"
                     << "\n\tcosine"
                     << "\n\tquartic.\n";
}

double MappingFilterFunction::GetDistance(
    const Array3DType& ICoord,
    const Array3DType& JCoord) const
{
    const Array3DType dist_vector = ICoord - JCoord;
    return sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
}

double MappingFilterFunction::ComputeWeight(
    const Array3DType& ICoord,
    const Array3DType& JCoord,
    const double Radius) const
{
    KRATOS_TRY;

    // Compute distance vector
    const double distance = GetDistance(ICoord, JCoord);

    // Depending on which weighting function is chosen, compute weight
    return mMappingFilterFunctional(Radius, distance);

    KRATOS_CATCH("");
}

} // namespace Kratos.
