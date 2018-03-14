// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "custom_linear_solvers/mixedulm_linear_solver.h"

namespace Kratos
{

namespace Python
{

void  AddCustomLinearSolversToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef IterativeSolver<SpaceType,  LocalSpaceType> IterativeSolverType;

    typedef MixedULMLinearSolver<SpaceType,  LocalSpaceType> MixedULMLinearSolverType;

    using namespace boost::python;

    class_<MixedULMLinearSolverType, MixedULMLinearSolverType::Pointer, bases<IterativeSolverType>, boost::noncopyable >("MixedULMLinearSolver",init<LinearSolverType::Pointer>())
    .def(init<LinearSolverType::Pointer ,double, const std::size_t >())
    .def(init<LinearSolverType::Pointer, Parameters>())
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

