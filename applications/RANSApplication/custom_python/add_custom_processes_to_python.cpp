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

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

    using RansWallFunctionUpdateProcessType = RansWallFunctionUpdateProcess;
    py::class_<RansWallFunctionUpdateProcessType, RansWallFunctionUpdateProcessType::Pointer, Process>(m, "RansWallFunctionUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    using RansKTurbulentIntensityInletProcessType = RansKTurbulentIntensityInletProcess;
    py::class_<RansKTurbulentIntensityInletProcessType, RansKTurbulentIntensityInletProcessType::Pointer, Process>(m, "RansKTurbulentIntensityInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutYPlusWallFunctionUpdateProcessType = RansNutYPlusWallFunctionUpdateProcess;
    py::class_<RansNutYPlusWallFunctionUpdateProcessType, RansNutYPlusWallFunctionUpdateProcessType::Pointer, Process>(m, "RansNutYPlusWallFunctionUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    // k-epsilon specific processes
    using RansEpsilonTurbulentMixingLengthInletProcessType = RansEpsilonTurbulentMixingLengthInletProcess;
    py::class_<RansEpsilonTurbulentMixingLengthInletProcessType, RansEpsilonTurbulentMixingLengthInletProcessType::Pointer, Process>(m, "RansEpsilonTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutKEpsilonUpdateProcessType = RansNutKEpsilonUpdateProcess;
    py::class_<RansNutKEpsilonUpdateProcessType, RansNutKEpsilonUpdateProcessType::Pointer, Process>(m, "RansNutKEpsilonUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    // k-omega specific processes
    using RansOmegaTurbulentMixingLengthInletProcessType = RansOmegaTurbulentMixingLengthInletProcess;
    py::class_<RansOmegaTurbulentMixingLengthInletProcessType, RansOmegaTurbulentMixingLengthInletProcessType::Pointer, Process>(m, "RansOmegaTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutKOmegaUpdateProcessType = RansNutKOmegaUpdateProcess;
    py::class_<RansNutKOmegaUpdateProcessType, RansNutKOmegaUpdateProcessType::Pointer, Process>(m, "RansNutKOmegaUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const int>());

    // k-omega-sst specific processes
    using RansWallDistanceCalculationProcessType = RansWallDistanceCalculationProcess<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<RansWallDistanceCalculationProcessType, RansWallDistanceCalculationProcessType::Pointer, Process>(m, "RansWallDistanceCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutKOmegaSSTUpdateProcessType = RansNutKOmegaSSTUpdateProcess;
    py::class_<RansNutKOmegaSSTUpdateProcessType, RansNutKOmegaSSTUpdateProcessType::Pointer, Process>(m, "RansNutKOmegaSSTUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const double, const int>());

    // misc. processes
    using RansApplyExactNodalPeriodicConditionProcessType = RansApplyExactNodalPeriodicConditionProcess;
    py::class_<RansApplyExactNodalPeriodicConditionProcessType, RansApplyExactNodalPeriodicConditionProcessType::Pointer, Process>(m, "RansApplyExactNodalPeriodicConditionProcess")
        .def(py::init<Model&, Parameters&>());

    using RansApplyFlagToSkinProcessType = RansApplyFlagToSkinProcess;
    py::class_<RansApplyFlagToSkinProcessType, RansApplyFlagToSkinProcessType::Pointer, Process>(m, "RansApplyFlagToSkinProcess")
        .def(py::init<Model&, Parameters&>());

    using RansClipScalarVariableProcessType = RansClipScalarVariableProcess;
    py::class_<RansClipScalarVariableProcessType, RansClipScalarVariableProcessType::Pointer, Process>(m, "RansClipScalarVariableProcess")
        .def(py::init<Model&, Parameters&>());

    using RansLineOutputProcessType = RansLineOutputProcess;
    py::class_<RansLineOutputProcessType, RansLineOutputProcessType::Pointer, Process>(m, "RansLineOutputProcess")
        .def(py::init<Model&, Parameters&>());

    using RansComputeReactionsProcessType = RansComputeReactionsProcess;
    py::class_<RansComputeReactionsProcessType, RansComputeReactionsProcessType::Pointer, Process>(m, "RansComputeReactionsProcess")
        .def(py::init<Model&, Parameters&>());
}
} // namespace Python
} // namespace Kratos
