// System includes

#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

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
#include "custom_processes/rans_epsilon_turbulent_mixing_inlet_process.h"
#include "custom_processes/rans_nut_k_epsilon_high_re_update_process.h"
#include "custom_processes/rans_omega_turbulent_mixing_inlet_process.h"
#include "custom_processes/rans_nut_k_omega_update_process.h"
#include "custom_processes/rans_wall_distance_calculation_process.h"
#include "custom_processes/rans_nut_k_omega_sst_update_process.h"

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
    py::class_<RansWallFunctionUpdateProcessType, RansWallFunctionUpdateProcessType::Pointer, Process>(
        m, "RansWallFunctionUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    using RansKTurbulentIntensityInletProcessType = RansKTurbulentIntensityInletProcess;
    py::class_<RansKTurbulentIntensityInletProcessType, RansKTurbulentIntensityInletProcessType::Pointer, Process>(m, "RansKTurbulentIntensityInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutYPlusWallFunctionUpdateProcessType = RansNutYPlusWallFunctionUpdateProcess;
    py::class_<RansNutYPlusWallFunctionUpdateProcessType, RansNutYPlusWallFunctionUpdateProcessType::Pointer, Process>(
        m, "RansNutYPlusWallFunctionUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    // k-epsilon specific processes
    using RansEpsilonTurbulentMixingLengthInletProcessType = RansEpsilonTurbulentMixingLengthInletProcess;
    py::class_<RansEpsilonTurbulentMixingLengthInletProcessType, RansEpsilonTurbulentMixingLengthInletProcessType::Pointer, Process>(m, "RansEpsilonTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutKEpsilonHighReUpdateProcessType = RansNutKEpsilonHighReUpdateProcess;
    py::class_<RansNutKEpsilonHighReUpdateProcessType, RansNutKEpsilonHighReUpdateProcessType::Pointer, Process>(
        m, "RansNutKEpsilonHighReUpdateProcess")
        .def(py::init<Model&, Parameters&>())
        .def(py::init<Model&, const std::string&, const double, const double, const int>());

    // k-omega specific processes
    using RansOmegaTurbulentMixingLengthInletProcessType = RansOmegaTurbulentMixingLengthInletProcess;
    py::class_<RansOmegaTurbulentMixingLengthInletProcessType, RansOmegaTurbulentMixingLengthInletProcessType::Pointer, Process>(m, "RansOmegaTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutKOmegaUpdateProcessType = RansNutKOmegaUpdateProcess;
    py::class_<RansNutKOmegaUpdateProcessType, RansNutKOmegaUpdateProcessType::Pointer, Process>(
        m, "RansNutKOmegaUpdateProcess")
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

}

} // namespace Python.
} // Namespace Kratos
