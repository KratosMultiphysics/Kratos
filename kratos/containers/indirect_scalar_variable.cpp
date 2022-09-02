//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes

// Application includes

// Include base h
#include "indirect_scalar_variable.h"

namespace Kratos
{

#ifdef KRATOS_SMP_CXX11
    thread_local IndirectScalarVariable::mDefaultValue = 0.0;
#else
    double IndirectScalarVariable::mDefaultValue = 0.0;
#endif

} // namespace Kratos
