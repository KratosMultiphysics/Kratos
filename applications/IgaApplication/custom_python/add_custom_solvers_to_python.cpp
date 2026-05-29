//
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Manuel Messmer
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_solvers_to_python.h"
#include "custom_solvers/additive_schwarz_preconditioner.h"
#include "custom_solvers/weak_dof_scaling_solver.h"
#include "factories/standard_linear_solver_factory.h"
#include "factories/standard_preconditioner_factory.h"

#include "iga_application_variables.h"


namespace Kratos {
namespace Python {

void AddCustomSolversToPython(
    pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;
    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
    typedef AdditiveSchwarzPreconditioner<SpaceType,  LocalSpaceType> AdditiveSchwarzPreconditionerType;
    typedef WeakDofScalingSolver<SpaceType, LocalSpaceType> WeakDofScalingSolverType;

    py::class_<AdditiveSchwarzPreconditionerType, AdditiveSchwarzPreconditionerType::Pointer, PreconditionerType>(m,"AdditiveSchwarzPreconditioner")
        .def(py::init<>() )
        .def("__str__", PrintObject<AdditiveSchwarzPreconditionerType>)
    ;

    py::class_<WeakDofScalingSolverType, WeakDofScalingSolverType::Pointer, LinearSolverType>(m,"WeakDofScalingSolver")
        .def(py::init<Parameters>())
        .def("__str__", PrintObject<WeakDofScalingSolverType>)
    ;

    static auto AdditiveSchwarzPreconditionerFactory = StandardPreconditionerFactory<SpaceType,LocalSpaceType,AdditiveSchwarzPreconditionerType>();
    static auto WeakDofScalingSolverFactory = StandardLinearSolverFactory<SpaceType, LocalSpaceType, WeakDofScalingSolverType>();
    KRATOS_REGISTER_PRECONDITIONER("additive_schwarz", AdditiveSchwarzPreconditionerFactory);
    KRATOS_REGISTER_LINEAR_SOLVER("weak_dof_scaling", WeakDofScalingSolverFactory);
}

} // namespace Python
} // Namespace Kratos
