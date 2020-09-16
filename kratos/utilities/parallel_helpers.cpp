//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//


// System includes

// External includes

// Project includes
#include "includes/kernel.h"
#include "utilities/parallel_helpers.h"


namespace Kratos {

int ParallelHelpers::GetNumThreads()
{
    return Kernel::GetNumThreads();
}

void ParallelHelpers::SetNumThreads(const int NumThreads)
{
    return Kernel::SetNumThreads(NumThreads);
}

}  // namespace Kratos.
