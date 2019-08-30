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

// Application includes
#include "custom_processes/solving_strategies/k_epsilon_co_solving_process.h"
#include "custom_processes/solving_strategies/scalar_co_solving_process.h"


// RANS auxiliary processes
#include "custom_processes/auxiliary_processes/rans_check_scalar_bounds_process.h"
#include "custom_processes/auxiliary_processes/rans_epsilon_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_apply_exact_nodal_periodic_condition_process.h"
#include "custom_processes/auxiliary_processes/rans_apply_flag_process.h"
#include "custom_processes/auxiliary_processes/rans_clip_scalar_variable_process.h"
#include "custom_processes/auxiliary_processes/rans_epsilon_turbulent_mixing_inlet_process.h"
#include "custom_processes/auxiliary_processes/rans_epsilon_wall_friction_velocity_process.h"
#include "custom_processes/auxiliary_processes/rans_find_condition_parent_process.h"
#include "custom_processes/auxiliary_processes/rans_k_turbulent_intensity_inlet_process.h"
#include "custom_processes/auxiliary_processes/rans_k_wall_friction_velocity_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_low_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_y_plus_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_scalar_cell_center_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_vector_align_process.h"
#include "custom_processes/auxiliary_processes/rans_vector_cell_center_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_wall_distance_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_check_scalar_condition_bounds_process.h"
#include "custom_processes/auxiliary_processes/rans_check_vector_bounds_process.h"

namespace Kratos
{
namespace Python
{
void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    // Adding solving strategies
    typedef ScalarCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> ScalarCoSolvingProcessType;
    py::class_<ScalarCoSolvingProcessType, ScalarCoSolvingProcessType::Pointer, Process>(
        m, "ScalarCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&, Variable<double>&>())
        .def("AddStrategy", &ScalarCoSolvingProcessType::AddStrategy)
        .def("AddAuxiliaryProcess", &ScalarCoSolvingProcessType::AddAuxiliaryProcess);

    typedef KEpsilonCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> KEpsilonCoSolvingProcessType;
    py::class_<KEpsilonCoSolvingProcessType, KEpsilonCoSolvingProcessType::Pointer, ScalarCoSolvingProcessType, Process>(
        m, "KEpsilonCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&>());

    // Adding auxiliary processes
    typedef RansNutKWallFunctionProcess RansNutKWallFunctionProcessType;
    py::class_<RansNutKWallFunctionProcessType, RansNutKWallFunctionProcessType::Pointer, Process>(
        m, "RansNutKWallFunctionProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansEpsilonWallFunctionProcess RansEpsilonWallFunctionProcessType;
    py::class_<RansEpsilonWallFunctionProcessType, RansEpsilonWallFunctionProcessType::Pointer, Process>(
        m, "RansEpsilonWallFunctionProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansCheckScalarBoundsProcess RansCheckScalarBoundsProcessType;
    py::class_<RansCheckScalarBoundsProcessType, RansCheckScalarBoundsProcessType::Pointer, Process>(
        m, "RansCheckScalarBoundsProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansCheckVectorBoundsProcess RansCheckVectorBoundsProcessType;
    py::class_<RansCheckVectorBoundsProcessType, RansCheckVectorBoundsProcessType::Pointer, Process>(
        m, "RansCheckVectorBoundsProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansCheckScalarConditionBoundsProcess RansCheckScalarConditionBoundsProcessType;
    py::class_<RansCheckScalarConditionBoundsProcessType, RansCheckScalarConditionBoundsProcessType::Pointer, Process>(
        m, "RansCheckScalarConditionBoundsProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansScalarCellCenterAveragingProcess RansScalarCellCenterAveragingProcessType;
    py::class_<RansScalarCellCenterAveragingProcessType, RansScalarCellCenterAveragingProcessType::Pointer, Process>(
        m, "RansScalarCellCenterAveragingProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansVectorCellCenterAveragingProcess RansVectorCellCenterAveragingProcessType;
    py::class_<RansVectorCellCenterAveragingProcessType, RansVectorCellCenterAveragingProcessType::Pointer, Process>(
        m, "RansVectorCellCenterAveragingProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansApplyFlagProcess RansApplyFlagProcessType;
    py::class_<RansApplyFlagProcessType, RansApplyFlagProcessType::Pointer, Process>(
        m, "RansApplyFlagProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansFindConditionParentProcess RansFindConditionParentProcessType;
    py::class_<RansFindConditionParentProcessType, RansFindConditionParentProcessType::Pointer, Process>(
        m, "RansFindConditionParentProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansVectorAlignProcess RansVectorAlignProcessType;
    py::class_<RansVectorAlignProcessType, RansVectorAlignProcessType::Pointer, Process>(
        m, "RansVectorAlignProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansWallDistanceCalculationProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> RansWallDistanceCalculationProcessType;
    py::class_<RansWallDistanceCalculationProcessType, RansWallDistanceCalculationProcessType::Pointer, Process>(
        m, "RansWallDistanceCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansLogarithmicYPlusCalculationProcess RansLogarithmicYPlusCalculationProcessType;
    py::class_<RansLogarithmicYPlusCalculationProcessType, RansLogarithmicYPlusCalculationProcessType::Pointer, Process>(
        m, "RansLogarithmicYPlusCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansNutKEpsilonHighReCalculationProcess RansNutKEpsilonHighReCalculationProcessType;
    py::class_<RansNutKEpsilonHighReCalculationProcessType, RansNutKEpsilonHighReCalculationProcessType::Pointer, Process>(
        m, "RansNutKEpsilonHighReCalculationProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansKTurbulentIntensityInletProcess RansKTurbulentIntensityInletProcessType;
    py::class_<RansKTurbulentIntensityInletProcessType, RansKTurbulentIntensityInletProcessType::Pointer, Process>(
        m, "RansKTurbulentIntensityInletProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansEpsilonTurbulentMixingLengthInletProcess RansEpsilonTurbulentMixingLengthInletProcessType;
    py::class_<RansEpsilonTurbulentMixingLengthInletProcessType,
               RansEpsilonTurbulentMixingLengthInletProcessType::Pointer, Process>(
        m, "RansEpsilonTurbulentMixingLengthInletProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansKWallFrictionVelocityProcess RansKWallFrictionVelocityProcessType;
    py::class_<RansKWallFrictionVelocityProcessType, RansKWallFrictionVelocityProcessType::Pointer, Process>(
        m, "RansKWallFrictionVelocityProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansEpsilonWallFrictionVelocityProcess RansEpsilonWallFrictionVelocityProcessType;
    py::class_<RansEpsilonWallFrictionVelocityProcessType, RansEpsilonWallFrictionVelocityProcessType::Pointer, Process>(
        m, "RansEpsilonWallFrictionVelocityProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansClipScalarVariableProcess RansClipScalarVariableProcessType;
    py::class_<RansClipScalarVariableProcessType, RansClipScalarVariableProcessType::Pointer, Process>(
        m, "RansClipScalarVariableProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansApplyExactNodalPeriodicConditionProcess RansApplyExactNodalPeriodicConditionProcessType;
    py::class_<RansApplyExactNodalPeriodicConditionProcessType,
               RansApplyExactNodalPeriodicConditionProcessType::Pointer, Process>(
        m, "RansApplyExactNodalPeriodicConditionProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansNutYPlusWallFunctionProcess RansNutYPlusWallFunctionProcessType;
    py::class_<RansNutYPlusWallFunctionProcessType, RansNutYPlusWallFunctionProcessType::Pointer, Process>(
        m, "RansNutYPlusWallFunctionProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansNutLowReCalculationProcess RansNutLowReCalculationProcessType;
    py::class_<RansNutLowReCalculationProcessType, RansNutLowReCalculationProcessType::Pointer, Process>(
        m, "RansNutLowReCalculationProcess")
        .def(py::init<Model&, Parameters&>());
}

} // namespace Python.
} // Namespace Kratos
