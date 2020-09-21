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
#include "python/add_convergence_accelerator_to_python.h"
#include "includes/define.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "trilinos_space.h"

namespace Kratos {
namespace Python {

void AddTrilinosConvergenceAcceleratorsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType  = UblasSpace<double, Matrix, Vector>;

    // Convergence accelerator base class
    AddBaseConvergenceAcceleratorToPython<TrilinosSparseSpaceType, TrilinosLocalSpaceType>(m, "TrilinosConvergenceAccelerator");

}

}  // namespace Python.
} // Namespace Kratos
