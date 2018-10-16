//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    //typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    //typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


}





}  // namespace Python.

} // Namespace Kratos

