//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// Project includes
#include "linear_solvers/linear_solver.h"
#include "spaces/eigen_space.h"
#include "spaces/ublas_space.h"

namespace Kratos {


KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) LinearSolver<
    TEigenSparseSpace<double>,
    TUblasDenseSpace<double>
>;


KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) LinearSolver<
    TEigenSparseSpace<float>,
    TUblasDenseSpace<double>
>;


} // namespace Kratos
