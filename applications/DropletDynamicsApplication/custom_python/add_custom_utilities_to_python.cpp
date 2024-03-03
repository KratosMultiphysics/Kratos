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
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/find_nodal_h_process_max.h"
#include "custom_utilities/find_conservative_elements.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    
    // Add FindNodalHProcessMax to Python
    py::class_<FindNodalHProcessMax<true>, FindNodalHProcessMax<true>::Pointer, Process>(m, "FindNodalHProcessMax")
        .def(py::init<ModelPart&>())
        .def("Execute", &FindNodalHProcessMax<true>::Execute);

    py::class_<FindNodalHProcessMax<false>, FindNodalHProcessMax<false>::Pointer, Process>(m, "FindNonHistoricalNodalHProcessMax")
        .def(py::init<ModelPart&>())
        .def("Execute", &FindNodalHProcessMax<false>::Execute);

    // Add FindConservativeElementsProcess to Python
    py::class_<FindConservativeElementsProcess<FindConservativeElementsSettings::SaveAsHistoricalVariable>, FindConservativeElementsProcess<FindConservativeElementsSettings::SaveAsHistoricalVariable>::Pointer, Process>(m,"FindNodalNighbersProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<FindConservativeElementsProcess<FindConservativeElementsSettings::SaveAsNonHistoricalVariable>, FindConservativeElementsProcess<FindConservativeElementsSettings::SaveAsNonHistoricalVariable>::Pointer, Process>(m,"FindNodalNighbersNonHistoricalProcess")
    .def(py::init<ModelPart&>())
    ;

}

} // namespace Python.
} // Namespace Kratos
