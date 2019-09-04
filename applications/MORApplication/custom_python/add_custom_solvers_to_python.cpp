/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Matthias Ebert
*/

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_solvers_to_python.h"
#include "linear_solvers/linear_solver.h"
#include "custom_solvers/gen_eigensystem_solver.h"

#include "spaces/ublas_space.h"

//#include "factories/standard_linear_solver_factory.h"
//#include "eigen_solvers_application.h"

namespace Kratos {
namespace Python {


void AddCustomSolversToPython(pybind11::module& m)
{
    namespace py = pybind11;

//TODO: check what's really needed

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef Preconditioner<SparseSpaceType, LocalSpaceType> PreconditionerType;
    typedef Reorderer<SparseSpaceType, LocalSpaceType> ReordererType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    typedef GenEigensystemSolver< SparseSpaceType, LocalSpaceType, PreconditionerType, ReordererType> EigenSystemSolverType;


    py::class_<EigenSystemSolverType, typename EigenSystemSolverType::Pointer, LinearSolverType >
        (m, "GenEigensystemSolver")
        .def(py::init<Parameters>())
    ;


// version for broke down complex version without LinSolver
/*

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef GenEigensystemSolver< SparseSpaceType, LocalSpaceType> EigenSystemSolverType;

    py::class_<EigenSystemSolverType, typename EigenSystemSolverType::Pointer >
        (m, "GenEigensystemSolver")
        .def(py::init<Parameters>())
    ;

*/

}

/*
void AddCustomSolversToPython(pybind11::module& m)
{
    //using complex = std::complex<double>;


    // --- eigensystem solver

    register_gen_eigensystem_solver<SparseLU<double>>(m, "GenEigensystemSolver");
 
}
*/

} // namespace Python
} // namespace Kratos
