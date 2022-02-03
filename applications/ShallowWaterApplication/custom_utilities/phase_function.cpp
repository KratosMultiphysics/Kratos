//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "phase_function.h"

namespace Kratos
{

double PhaseFunction::InverseHeight(const double Height, const double Epsilon)
{
    const double h4 = std::pow(Height, 4);
    const double epsilon4 = std::pow(Epsilon, 4);
    return std::sqrt(2) * std::max(Height, .0) / std::sqrt(h4 + std::max(h4, epsilon4));
}

double PhaseFunction::WetFraction(const double Height, const double Epsilon)
{
    return Height * InverseHeight(Height, Epsilon);
}

}  // namespace Kratos.
