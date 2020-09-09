//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_processes/weighted_divergence_calculation_process.h"
#include "custom_processes/metrics_divergencefree_process.h"
#include "custom_processes/calculate_divergence_process.h"


namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<WeightedDivergenceCalculationProcess, WeightedDivergenceCalculationProcess::Pointer, Process >
        (m, "WeightedDivergenceCalculationProcess")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<MetricDivergenceFreeProcess<2>, MetricDivergenceFreeProcess<2>::Pointer, Process>(m, "MetricDivergenceFreeProcess2D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<MetricDivergenceFreeProcess<3>, MetricDivergenceFreeProcess<3>::Pointer, Process>(m, "MetricDivergenceFreeProcess3D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<CalculateDivergenceProcess, CalculateDivergenceProcess::Pointer, Process >
        (m, "DivergenceProcess")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, Parameters>())
    ;

}

} // namespace Python.
} // Namespace Kratos
