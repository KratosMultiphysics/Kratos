//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "add_convergence_accelerators_to_python.h"
#include "add_convergence_accelerator_to_python.h"
#include "includes/define.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "spaces/ublas_space.h"

namespace Kratos
{

namespace Python
{

void AddConvergenceAcceleratorsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

    // Convergence accelerator base class
    AddBaseConvergenceAcceleratorToPython<SparseSpaceType, LocalSpaceType>(m, "ConvergenceAccelerator");

}

}  // namespace Python.

} // Namespace Kratos
