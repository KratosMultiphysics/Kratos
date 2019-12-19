//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Saransh Saxena
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_processes/metrics_temperature_gradient_process.h"


namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /*typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;*/

    py::class_<MetricsTemperatureGradientProcess<2>, MetricsTemperatureGradientProcess<2>::Pointer, Process>(m, "MetricsTemperatureGradientProcess2D")
    .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<MetricsTemperatureGradientProcess<3>, MetricsTemperatureGradientProcess<3>::Pointer, Process>(m, "MetricsTemperatureGradientProcess3D")
    .def(py::init<ModelPart&, Parameters>())
;

}

} // namespace Python.
} // Namespace Kratos