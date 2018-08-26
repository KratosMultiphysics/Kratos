//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#if !defined(KRATOS_STRESS_RESPONSE_DEFINITIONS_H_INCLUDED )
#define  KRATOS_STRESS_RESPONSE_DEFINITIONS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
    enum class TracedStressType
    {
        FX,
        FY,
        FZ,
        MX,
        MY,
        MZ,
        FXX,
        FXY,
        FXZ,
        FYX,
        FYY,
        FYZ,
        FZX,
        FZY,
        FZZ,
        MXX,
        MXY,
        MXZ,
        MYX,
        MYY,
        MYZ,
        MZX,
        MZY,
        MZZ,
        PK2
    };

    enum class StressTreatment
    {
        Mean,
        Node,
        GaussPoint
    };

    namespace StressResponseDefinitions
    {

        TracedStressType ConvertStringToTracedStressType(const std::string& Str);

        StressTreatment ConvertStringToStressTreatment(const std::string& Str);

    } // namespace StressResponseDefinitions.


}  // namespace Kratos.

#endif // KRATOS_STRESS_RESPONSE_DEFINITIONS_H_INCLUDED  defined


