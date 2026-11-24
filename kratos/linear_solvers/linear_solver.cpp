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

// Project includes
#include "linear_solvers/linear_solver.h"

namespace Kratos {

template class LinearSolver<
    TDefaultSparseSpace<double>,
    TDefaultDenseSpace<double>
>;

template class LinearSolver<
    TDefaultSparseSpace<float>,
    TDefaultDenseSpace<double>
>;

} // namespace Kratos
