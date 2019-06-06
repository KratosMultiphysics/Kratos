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
#include "custom_processes/solving_strategies/k_epsilon_steady_co_solving_process.h"
#include "custom_processes/solving_strategies/scalar_co_solving_process.h"

// RANS Y Plus models
#include "custom_processes/y_plus_model_processes/rans_logarithmic_y_plus_model_process.h"
#include "custom_processes/y_plus_model_processes/rans_tke_y_plus_model_process.h"

// RANS initialization processes
#include "custom_processes/k_epsilon_evaluation_utau_process.h"

// RANS wall processes
#include "custom_processes/wall_processes/rans_exact_wall_distance_calculation_process.h"
#include "custom_processes/wall_processes/rans_wall_velocity_calculation_process.h"

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
        .def("AddStrategy", &ScalarCoSolvingProcessType::AddStrategy);

    typedef KEpsilonCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> KEpsilonCoSolvingProcessType;
    py::class_<KEpsilonCoSolvingProcessType, KEpsilonCoSolvingProcessType::Pointer, ScalarCoSolvingProcessType, Process>(
        m, "KEpsilonCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&, Process&>());

    typedef KEpsilonSteadyCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType> KEpsilonSteadyCoSolvingProcessType;
    py::class_<KEpsilonSteadyCoSolvingProcessType, KEpsilonSteadyCoSolvingProcessType::Pointer, ScalarCoSolvingProcessType, Process>(
        m, "KEpsilonSteadyCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters&, Process&>());

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
        .def(py::init<ModelPart&>());
}

} // namespace Python.
} // Namespace Kratos
