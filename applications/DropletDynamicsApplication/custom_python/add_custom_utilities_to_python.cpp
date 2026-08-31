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


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "spaces/default_spaces.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef TDefaultSparseSpace<double> SparseSpaceType;
    typedef TDefaultDenseSpace<double> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

}

} // namespace Python.
} // Namespace Kratos
