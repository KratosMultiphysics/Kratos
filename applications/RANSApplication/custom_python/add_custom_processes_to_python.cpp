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
#include "pybind11/stl.h"

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
#include "custom_processes/rans_epsilon_turbulent_mixing_length_inlet_process.h"
#include "custom_processes/rans_omega_turbulent_mixing_length_inlet_process.h"
#include "custom_processes/rans_wall_distance_calculation_process.h"
#include "custom_processes/rans_apply_exact_nodal_periodic_condition_process.h"
#include "custom_processes/rans_apply_flag_to_skin_process.h"
#include "custom_processes/rans_clip_scalar_variable_process.h"
#include "custom_processes/rans_line_output_process.h"
#include "custom_processes/rans_nut_nodal_update_process.h"
#include "custom_processes/rans_compute_reactions_process.h"
#include "custom_processes/rans_variable_data_transfer_process.h"
#include "custom_processes/rans_initialize_bossak_previous_step_variable_derivatives_process.h"
#include "custom_processes/rans_omega_viscous_log_wall_process.h"
#include "custom_processes/rans_omega_viscous_log_binomial_wall_process.h"
#include "custom_processes/rans_wall_properties_update_process.h"
#include "custom_processes/rans_compute_y_plus_process.h"
#include "custom_processes/rans_vtk_output_process.h"
#include "custom_processes/rans_omega_automatic_inlet_process.h"
#include "custom_processes/rans_smooth_clip_scalar_variable_process.h"

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

    py::class_<RansInitializeBossakPreviousStepVariableDerivatives, RansInitializeBossakPreviousStepVariableDerivatives::Pointer, Process>(m, "RansInitializeBossakPreviousStepVariableDerivatives")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansOmegaAutomaticInletProcess, RansOmegaAutomaticInletProcess::Pointer, Process>(m, "RansOmegaAutomaticInletProcess")
        .def(py::init<Model&, Parameters&>());

    // adding RansFormulationProcesses
    py::class_<RansFormulationProcess, RansFormulationProcess::Pointer, Process>(m, "RansFormulationProcess")
        .def(py::init<>())
        .def("ExecuteBeforeCouplingSolveStep", &RansFormulationProcess::ExecuteBeforeCouplingSolveStep)
        .def("ExecuteAfterCouplingSolveStep", &RansFormulationProcess::ExecuteAfterCouplingSolveStep);

    py::class_<RansClipScalarVariableProcess, RansClipScalarVariableProcess::Pointer, RansFormulationProcess>(m, "RansClipScalarVariableProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansNutNodalUpdateProcess, RansNutNodalUpdateProcess::Pointer, RansFormulationProcess>(m, "RansNutNodalUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const int>());

    py::class_<RansWallFunctionUpdateProcess, RansWallFunctionUpdateProcess::Pointer, RansFormulationProcess>(m, "RansWallFunctionUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const int>());

    py::class_<RansWallDistanceCalculationProcess, RansWallDistanceCalculationProcess::Pointer, RansFormulationProcess>(m, "RansWallDistanceCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansComputeReactionsProcess, RansComputeReactionsProcess::Pointer, RansFormulationProcess>(m, "RansComputeReactionsProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const std::vector<std::string>&, const int>());

    py::class_<RansVariableDataTransferProcess, RansVariableDataTransferProcess::Pointer, RansFormulationProcess>(m, "RansVariableDataTransferProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const std::string&, const std::vector<std::string>&, const std::vector<std::tuple<const std::string, const bool, const int, const std::string, const bool, const int>>&, const int>())
        .def(py::init<Model&, Model&, const std::string&, const std::string&, const std::vector<std::string>&, const std::vector<std::tuple<const std::string, const bool, const int, const std::string, const bool, const int>>&, const int>());

    py::class_<RansWallPropertiesUpdateProcess, RansWallPropertiesUpdateProcess::Pointer, RansFormulationProcess>(m, "RansWallPropertiesUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const bool, const bool, const std::vector<std::string>&, const int>());

    py::class_<RansOmegaViscousLogWallProcess, RansOmegaViscousLogWallProcess::Pointer, RansFormulationProcess>(m, "RansOmegaViscousLogWallProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansOmegaViscousLogBinomialWallProcess, RansOmegaViscousLogBinomialWallProcess::Pointer, RansFormulationProcess>(m, "RansOmegaViscousLogBinomialWallProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansComputeYPlusProcess, RansComputeYPlusProcess::Pointer, RansFormulationProcess>(m, "RansComputeYPlusProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansVTKOutputProcess, RansVTKOutputProcess::Pointer, RansFormulationProcess>(m, "RansVTKOutputProcess")
        .def(py::init<Model&, Parameters&>());

    py::class_<RansSmoothClipScalarVariableProcess, RansSmoothClipScalarVariableProcess::Pointer, RansFormulationProcess>(m, "RansSmoothClipScalarVariableProcess")
        .def(py::init<Model&, Parameters&>());
}
} // namespace Python
} // namespace Kratos
