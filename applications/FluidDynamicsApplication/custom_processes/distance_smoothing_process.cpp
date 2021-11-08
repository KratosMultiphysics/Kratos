//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

// Project includes

// Application includes
#include "distance_smoothing_process.h"

namespace Kratos
{

typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

template class DistanceSmoothingProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>;
template class DistanceSmoothingProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>;

}; // namespace Kratos
