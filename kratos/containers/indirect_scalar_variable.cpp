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

#if defined(KRATOS_SMP_OPENMP) && defined(KRATOS_COMPILED_IN_WINDOWS)
    std::vector<double> IndirectScalarVariable::mDefaultValues = std::vector<double>(OpenMPUtils::GetNumThreads(), 0.0);
#elif defined(KRATOS_SMP_CXX11)
    thread_local double IndirectScalarVariable::mDefaultValue = 0.0;
#else
    double IndirectScalarVariable::mDefaultValue = 0.0;
#endif

} // namespace Kratos
