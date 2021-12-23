// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "custom_linear_solvers/mixedulm_linear_solver.h"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void  AddCustomLinearSolversToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef IterativeSolver<SpaceType,  LocalSpaceType> IterativeSolverType;

    typedef MixedULMLinearSolver<SpaceType,  LocalSpaceType> MixedULMLinearSolverType;

    py::class_<MixedULMLinearSolverType, typename MixedULMLinearSolverType::Pointer, IterativeSolverType>(m, "MixedULMLinearSolver")
    .def(py::init<LinearSolverType::Pointer >())
    .def(py::init<LinearSolverType::Pointer,double, const std::size_t >())
    .def(py::init<LinearSolverType::Pointer, Parameters>())
    ;
}

}  // namespace Python.

} // Namespace Kratos

