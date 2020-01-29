// System includes

#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

#include "custom_python/add_custom_auxiliary_processes_to_python.h"

// RANS auxiliary processes
#include "custom_processes/auxiliary_processes/rans_apply_exact_nodal_periodic_condition_process.h"
#include "custom_processes/auxiliary_processes/rans_apply_flag_process.h"
#include "custom_processes/auxiliary_processes/rans_check_scalar_bounds_process.h"
#include "custom_processes/auxiliary_processes/rans_check_vector_bounds_process.h"
#include "custom_processes/auxiliary_processes/rans_clip_scalar_variable_process.h"
#include "custom_processes/auxiliary_processes/rans_epsilon_turbulent_mixing_inlet_process.h"
#include "custom_processes/auxiliary_processes/rans_find_condition_parent_process.h"
#include "custom_processes/auxiliary_processes/rans_k_turbulent_intensity_inlet_process.h"
#include "custom_processes/auxiliary_processes/rans_line_output_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_low_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_y_plus_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_scalar_cell_center_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_vector_align_process.h"
#include "custom_processes/auxiliary_processes/rans_vector_cell_center_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_wall_distance_calculation_process.h"

// RANS sensitivity processes
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_velocity_sensitivities_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_y_plus_wall_function_sensitivities_process.h"

namespace Kratos
{
namespace Python
{
void AddCustomAuxiliaryProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

    using RansCheckScalarBoundsProcessType = RansCheckScalarBoundsProcess;
    py::class_<RansCheckScalarBoundsProcessType, RansCheckScalarBoundsProcessType::Pointer, Process>(
        m, "RansCheckScalarBoundsProcess")
        .def(py::init<Model&, Parameters&>());

    using RansCheckVectorBoundsProcessType = RansCheckVectorBoundsProcess;
    py::class_<RansCheckVectorBoundsProcessType, RansCheckVectorBoundsProcessType::Pointer, Process>(
        m, "RansCheckVectorBoundsProcess")
        .def(py::init<Model&, Parameters&>());

    using RansScalarCellCenterAveragingProcessType = RansScalarCellCenterAveragingProcess;
    py::class_<RansScalarCellCenterAveragingProcessType, RansScalarCellCenterAveragingProcessType::Pointer, Process>(
        m, "RansScalarCellCenterAveragingProcess")
        .def(py::init<Model&, Parameters&>());

    using RansVectorCellCenterAveragingProcessType = RansVectorCellCenterAveragingProcess;
    py::class_<RansVectorCellCenterAveragingProcessType, RansVectorCellCenterAveragingProcessType::Pointer, Process>(
        m, "RansVectorCellCenterAveragingProcess")
        .def(py::init<Model&, Parameters&>());

    using RansApplyFlagProcessType = RansApplyFlagProcess;
    py::class_<RansApplyFlagProcessType, RansApplyFlagProcessType::Pointer, Process>(
        m, "RansApplyFlagProcess")
        .def(py::init<Model&, Parameters&>());

    using RansFindConditionParentProcessType = RansFindConditionParentProcess;
    py::class_<RansFindConditionParentProcessType, RansFindConditionParentProcessType::Pointer, Process>(
        m, "RansFindConditionParentProcess")
        .def(py::init<Model&, Parameters&>());

    using RansVectorAlignProcessType = RansVectorAlignProcess;
    py::class_<RansVectorAlignProcessType, RansVectorAlignProcessType::Pointer, Process>(
        m, "RansVectorAlignProcess")
        .def(py::init<Model&, Parameters&>());

    using RansWallDistanceCalculationProcessType =
        RansWallDistanceCalculationProcess<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<RansWallDistanceCalculationProcessType, RansWallDistanceCalculationProcessType::Pointer, Process>(
        m, "RansWallDistanceCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    using RansLogarithmicYPlusCalculationProcessType = RansLogarithmicYPlusCalculationProcess;
    py::class_<RansLogarithmicYPlusCalculationProcessType, RansLogarithmicYPlusCalculationProcessType::Pointer, Process>(
        m, "RansLogarithmicYPlusCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutKEpsilonHighReCalculationProcessType = RansNutKEpsilonHighReCalculationProcess;
    py::class_<RansNutKEpsilonHighReCalculationProcessType, RansNutKEpsilonHighReCalculationProcessType::Pointer, Process>(
        m, "RansNutKEpsilonHighReCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    using RansKTurbulentIntensityInletProcessType = RansKTurbulentIntensityInletProcess;
    py::class_<RansKTurbulentIntensityInletProcessType, RansKTurbulentIntensityInletProcessType::Pointer, Process>(
        m, "RansKTurbulentIntensityInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansEpsilonTurbulentMixingLengthInletProcessType =
        RansEpsilonTurbulentMixingLengthInletProcess;
    py::class_<RansEpsilonTurbulentMixingLengthInletProcessType,
               RansEpsilonTurbulentMixingLengthInletProcessType::Pointer, Process>(
        m, "RansEpsilonTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    using RansClipScalarVariableProcessType = RansClipScalarVariableProcess;
    py::class_<RansClipScalarVariableProcessType, RansClipScalarVariableProcessType::Pointer, Process>(
        m, "RansClipScalarVariableProcess")
        .def(py::init<Model&, Parameters&>());

    using RansApplyExactNodalPeriodicConditionProcessType =
        RansApplyExactNodalPeriodicConditionProcess;
    py::class_<RansApplyExactNodalPeriodicConditionProcessType,
               RansApplyExactNodalPeriodicConditionProcessType::Pointer, Process>(
        m, "RansApplyExactNodalPeriodicConditionProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutYPlusWallFunctionProcessType = RansNutYPlusWallFunctionProcess;
    py::class_<RansNutYPlusWallFunctionProcessType, RansNutYPlusWallFunctionProcessType::Pointer, Process>(
        m, "RansNutYPlusWallFunctionProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutLowReCalculationProcessType = RansNutLowReCalculationProcess;
    py::class_<RansNutLowReCalculationProcessType, RansNutLowReCalculationProcessType::Pointer, Process>(
        m, "RansNutLowReCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    using RansLineOutputProcessType = RansLineOutputProcess;
    py::class_<RansLineOutputProcessType, RansLineOutputProcessType::Pointer, Process>(
        m, "RansLineOutputProcess")
        .def(py::init<Model&, Parameters&>());

    using RansLogarithmicYPlusVelocitySensitivitiesProcessType = RansLogarithmicYPlusVelocitySensitivitiesProcess;
    py::class_<RansLogarithmicYPlusVelocitySensitivitiesProcessType, RansLogarithmicYPlusVelocitySensitivitiesProcessType::Pointer, Process>(
        m, "RansLogarithmicYPlusVelocitySensitivitiesProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutKEpsilonHighReSensitivitiesProcessType = RansNutKEpsilonHighReSensitivitiesProcess;
    py::class_<RansNutKEpsilonHighReSensitivitiesProcessType, RansNutKEpsilonHighReSensitivitiesProcessType::Pointer, Process>(
        m, "RansNutKEpsilonHighReSensitivitiesProcess")
        .def(py::init<Model&, Parameters&>());

    using RansNutYPlusWallFunctionSensitivitiesProcessType = RansNutYPlusWallFunctionSensitivitiesProcess;
    py::class_<RansNutYPlusWallFunctionSensitivitiesProcessType, RansNutYPlusWallFunctionSensitivitiesProcessType::Pointer, Process>(
        m, "RansNutYPlusWallFunctionSensitivitiesProcess")
        .def(py::init<Model&, Parameters&>());
}

} // namespace Python.
} // Namespace Kratos
