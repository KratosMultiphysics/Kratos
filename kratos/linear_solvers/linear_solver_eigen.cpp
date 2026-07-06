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
#include "linear_solvers/linear_solver_eigen.h"


namespace Kratos {


template class LinearSolver<
    TEigenSparseSpace<double>,
    TUblasDenseSpace<double>
>;


template class LinearSolver<
    TEigenSparseSpace<float>,
    TUblasDenseSpace<double>
>;


} // namespace Kratos
