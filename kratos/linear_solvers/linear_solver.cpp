//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Project Includes
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"


namespace Kratos {


template class LinearSolver<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;

template class LinearSolver<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos
