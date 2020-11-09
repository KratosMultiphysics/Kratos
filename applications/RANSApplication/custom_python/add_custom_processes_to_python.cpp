//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"

// Application includes
#include "custom_processes/rans_formulation_process.h"
#include "custom_processes/rans_wall_function_update_process.h"
#include "custom_processes/rans_k_turbulent_intensity_inlet_process.h"
#include "custom_processes/rans_nut_y_plus_wall_function_update_process.h"
#include "custom_processes/rans_epsilon_turbulent_mixing_length_inlet_process.h"
#include "custom_processes/rans_nut_k_epsilon_update_process.h"
#include "custom_processes/rans_omega_turbulent_mixing_length_inlet_process.h"
#include "custom_processes/rans_nut_k_omega_update_process.h"
#include "custom_processes/rans_wall_distance_calculation_process.h"
#include "custom_processes/rans_nut_k_omega_sst_update_process.h"
#include "custom_processes/rans_apply_exact_nodal_periodic_condition_process.h"
#include "custom_processes/rans_apply_flag_to_skin_process.h"
#include "custom_processes/rans_clip_scalar_variable_process.h"
#include "custom_processes/rans_line_output_process.h"
#include "custom_processes/rans_compute_reactions_process.h"

// Include base h
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RansKTurbulentIntensityInletProcess, RansKTurbulentIntensityInletProcess::Pointer, Process>(m, "RansKTurbulentIntensityInletProcess")
        .def(py::init<Model&, Parameters&>());

    // k-epsilon specific processes
    py::class_<RansEpsilonTurbulentMixingLengthInletProcess, RansEpsilonTurbulentMixingLengthInletProcess::Pointer, Process>(m, "RansEpsilonTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    // k-omega specific processes
    py::class_<RansOmegaTurbulentMixingLengthInletProcess, RansOmegaTurbulentMixingLengthInletProcess::Pointer, Process>(m, "RansOmegaTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    // misc. processes
    py::class_<RansApplyExactNodalPeriodicConditionProcess, RansApplyExactNodalPeriodicConditionProcess::Pointer, Process>(m, "RansApplyExactNodalPeriodicConditionProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansApplyFlagToSkinProcess, RansApplyFlagToSkinProcess::Pointer, Process>(m, "RansApplyFlagToSkinProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansLineOutputProcess, RansLineOutputProcess::Pointer, Process>(m, "RansLineOutputProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansComputeReactionsProcess, RansComputeReactionsProcess::Pointer, Process>(m, "RansComputeReactionsProcess")
        .def(py::init<Model&, Parameters&>());

    // adding RansFormulationProcesses
    py::class_<RansFormulationProcess, RansFormulationProcess::Pointer, Process>(m, "RansFormulationProcess")
        .def(py::init<>())
        .def("ExecuteBeforeCouplingSolveStep", &RansFormulationProcess::ExecuteBeforeCouplingSolveStep)
        .def("ExecuteAfterCouplingSolveStep", &RansFormulationProcess::ExecuteAfterCouplingSolveStep);

    py::class_<RansClipScalarVariableProcess, RansClipScalarVariableProcess::Pointer, RansFormulationProcess>(m, "RansClipScalarVariableProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansNutKEpsilonUpdateProcess, RansNutKEpsilonUpdateProcess::Pointer, RansFormulationProcess>(m, "RansNutKEpsilonUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    py::class_<RansNutKOmegaSSTUpdateProcess, RansNutKOmegaSSTUpdateProcess::Pointer, RansFormulationProcess>(m, "RansNutKOmegaSSTUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const double, const int>());

    py::class_<RansNutKOmegaUpdateProcess, RansNutKOmegaUpdateProcess::Pointer, RansFormulationProcess>(m, "RansNutKOmegaUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const int>());

    py::class_<RansNutYPlusWallFunctionUpdateProcess, RansNutYPlusWallFunctionUpdateProcess::Pointer, RansFormulationProcess>(m, "RansNutYPlusWallFunctionUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    py::class_<RansWallFunctionUpdateProcess, RansWallFunctionUpdateProcess::Pointer, RansFormulationProcess>(m, "RansWallFunctionUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    py::class_<RansWallDistanceCalculationProcess, RansWallDistanceCalculationProcess::Pointer, RansFormulationProcess>(m, "RansWallDistanceCalculationProcess")
        .def(py::init<Model&, Parameters&>());

}
} // namespace Python
} // namespace Kratos
