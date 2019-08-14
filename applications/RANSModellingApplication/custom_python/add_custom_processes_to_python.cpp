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
// #include "custom_processes/solving_strategies/k_epsilon_steady_co_solving_process.h"
#include "custom_processes/solving_strategies/scalar_co_solving_process.h"

// RANS Y Plus models
#include "custom_processes/y_plus_model_processes/rans_logarithmic_y_plus_model_process.h"
#include "custom_processes/y_plus_model_processes/rans_tke_y_plus_model_process.h"

// RANS initialization processes
#include "custom_processes/epsilon_turbulent_mixing_length_evaluation_process.h"
#include "custom_processes/k_epsilon_evaluation_utau_process.h"
#include "custom_processes/k_turbulent_intensity_evaluation_process.h"

// RANS wall processes
#include "custom_processes/wall_processes/rans_exact_wall_distance_calculation_process.h"
#include "custom_processes/wall_processes/rans_wall_velocity_calculation_process.h"

// RANS auxiliary processes
#include "custom_processes/auxiliary_processes/rans_check_scalar_bounds_process.h"
#include "custom_processes/auxiliary_processes/rans_epsilon_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_scalar_neighbour_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_y_plus_wall_distance_calculation_process.h"

#include "custom_processes/auxiliary_processes/rans_apply_exact_nodal_periodic_condition_process.h"
#include "custom_processes/auxiliary_processes/rans_apply_flag_process.h"
#include "custom_processes/auxiliary_processes/rans_clip_scalar_variable_by_neighbour_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_clip_scalar_variable_process.h"
#include "custom_processes/auxiliary_processes/rans_epsilon_turbulent_mixing_inlet_process.h"
#include "custom_processes/auxiliary_processes/rans_epsilon_wall_friction_velocity_process.h"
#include "custom_processes/auxiliary_processes/rans_find_condition_parent_process.h"
#include "custom_processes/auxiliary_processes/rans_k_turbulent_intensity_inlet_process.h"
#include "custom_processes/auxiliary_processes/rans_k_wall_friction_velocity_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_low_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_y_plus_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_scalar_cell_center_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_vector_align_process.h"
#include "custom_processes/auxiliary_processes/rans_vector_cell_center_averaging_process.h"
#include "custom_processes/auxiliary_processes/rans_wall_distance_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_y_plus_k_calculation_process.h"
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

    // typedef KEpsilonSteadyCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> KEpsilonSteadyCoSolvingProcessType;
    // py::class_<KEpsilonSteadyCoSolvingProcessType, KEpsilonSteadyCoSolvingProcessType::Pointer, ScalarCoSolvingProcessType, Process>(
    //     m, "KEpsilonSteadyCoSolvingProcess")
    //     .def(py::init<ModelPart&, Parameters&, Process&, Process&>());

    // Adding y_plus calculation models
    py::class_<RansLogarithmicYPlusModelProcess, RansLogarithmicYPlusModelProcess::Pointer, Process>(
        m, "RansLogarithmicYPlusModelProcess")
        .def(py::init<ModelPart&, Parameters&>());
    py::class_<RansTKEYPlusModelProcess, RansTKEYPlusModelProcess::Pointer, Process>(
        m, "RansTKEYPlusModelProcess")
        .def(py::init<ModelPart&, Parameters&>());

    // Adding initialization processes
    py::class_<RansKEpsilonEvaluationUtauProcess, RansKEpsilonEvaluationUtauProcess::Pointer, Process>(
        m, "RansKEpsilonEvaluationUtauProcess")
        .def(py::init<Model&, Parameters&, Process&>());

    py::class_<RansKTurbulentIntensityEvaluationProcess, RansKTurbulentIntensityEvaluationProcess::Pointer, Process>(
        m, "RansKTurbulentIntensityEvaluationProcess")
        .def(py::init<ModelPart&, Parameters&, const bool>());

    py::class_<RansEpsilonTurbulentMixingLengthEvaluationProcess,
               RansEpsilonTurbulentMixingLengthEvaluationProcess::Pointer, Process>(
        m, "RansEpsilonTurbulentMixingLengthEvaluationProcess")
        .def(py::init<ModelPart&, Parameters&, const bool>());

    // Adding wall distance calculation processes
    typedef RansExactWallDistanceCalculationProcess<2, SparseSpaceType, LocalSpaceType, LinearSolverType> RansExactWallDistanceCalculationProcessType2D;
    py::class_<RansExactWallDistanceCalculationProcessType2D, RansExactWallDistanceCalculationProcessType2D::Pointer, Process>(
        m, "RansExactWallDistanceCalculationProcess2D")
        .def(py::init<ModelPart&, Parameters&>());

    typedef RansExactWallDistanceCalculationProcess<3, SparseSpaceType, LocalSpaceType, LinearSolverType> RansExactWallDistanceCalculationProcessType3D;
    py::class_<RansExactWallDistanceCalculationProcessType3D, RansExactWallDistanceCalculationProcessType3D::Pointer, Process>(
        m, "RansExactWallDistanceCalculationProcess3D")
        .def(py::init<ModelPart&, Parameters&>());

    // Adding wall velocity calculation processes
    typedef RansWallVelocityCalculationProcess RansWallVelocityCalculationProcessType;
    py::class_<RansWallVelocityCalculationProcessType, RansWallVelocityCalculationProcessType::Pointer, Process>(
        m, "RansWallVelocityCalculationProcess")
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

    typedef RansScalarNeighbourAveragingProcess RansScalarNeighbourAveragingProcessType;
    py::class_<RansScalarNeighbourAveragingProcessType, RansScalarNeighbourAveragingProcessType::Pointer, Process>(
        m, "RansScalarNeighbourAveragingProcess")
        .def(py::init<ModelPart&, Parameters&>());

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

    typedef RansYPlusWallDistanceCalculationProcess RansYPlusWallDistanceCalculationProcessType;
    py::class_<RansYPlusWallDistanceCalculationProcessType, RansYPlusWallDistanceCalculationProcessType::Pointer, Process>(
        m, "RansYPlusWallDistanceCalculationProcess")
        .def(py::init<ModelPart&, Parameters&>());

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

    typedef RansNutHighReCalculationProcess RansNutHighReCalculationProcessType;
    py::class_<RansNutHighReCalculationProcessType, RansNutHighReCalculationProcessType::Pointer, Process>(
        m, "RansNutHighReCalculationProcess")
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

    typedef RansClipScalarVariableByNeighbourAveragingProcess RansClipScalarVariableByNeighbourAveragingProcessType;
    py::class_<RansClipScalarVariableByNeighbourAveragingProcessType,
               RansClipScalarVariableByNeighbourAveragingProcessType::Pointer, Process>(
        m, "RansClipScalarVariableByNeighbourAveragingProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansApplyExactNodalPeriodicConditionProcess RansApplyExactNodalPeriodicConditionProcessType;
    py::class_<RansApplyExactNodalPeriodicConditionProcessType,
               RansApplyExactNodalPeriodicConditionProcessType::Pointer, Process>(
        m, "RansApplyExactNodalPeriodicConditionProcess")
        .def(py::init<Model&, Parameters&>());

    typedef RansYPlusKCalculationProcess RansYPlusKCalculationProcessType;
    py::class_<RansYPlusKCalculationProcessType, RansYPlusKCalculationProcessType::Pointer, Process>(
        m, "RansYPlusKCalculationProcess")
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
