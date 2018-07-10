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

#if !defined(KRATOS_RESPONSE_DATA_H_INCLUDED )
#define  KRATOS_RESPONSE_DATA_H_INCLUDED


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
        MZZ
    };

    enum class StressTreatment
    {
        Mean,
        Node,
        GaussPoint,
        StressTreatmentNotAvailable
    };

    namespace ResponseData
    {
        
        TracedStressType ConvertStressType(const std::string& Str);

        StressTreatment ConvertStressTreatment(const std::string& Str);

    } // namespace ResponseData.


}  // namespace Kratos.

#endif // KRATOS_RESPONSE_DATA_H_INCLUDED  defined


