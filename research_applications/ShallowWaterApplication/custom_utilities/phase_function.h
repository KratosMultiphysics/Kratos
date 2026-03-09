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

#ifndef KRATOS_PHASE_FUNCTION_H_INCLUDED
#define KRATOS_PHASE_FUNCTION_H_INCLUDED


// System includes
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @ingroup ShallowWaterApplication
 * @class PhaseFunction
 * @brief This class is a wrapper of useful utilities for shallow water computations
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) PhaseFunction
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PhaseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static double InverseHeight(const double Height, const double Epsilon);

    static double WetFraction(double Height, double Epsilon);

    ///@}

}; // Class PhaseFunction

///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_PHASE_FUNCTION_H_INCLUDED  defined
