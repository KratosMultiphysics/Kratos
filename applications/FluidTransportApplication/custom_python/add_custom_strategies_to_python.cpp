//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/kratos_parameters.h"

//schemes
#include "custom_strategies/schemes/generalized_newmark_GN11_scheme.hpp"
#include "custom_strategies/schemes/explicit_forward_euler_scheme.hpp"

//strategies
//#include "solving_strategies/strategies/solving_strategy.h"

//builders and solvers

//linear solvers
//#include "linear_solvers/linear_solver.h"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    //typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    //typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    //typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    //typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;

    typedef GeneralizedNewmarkGN11Scheme< SparseSpaceType, LocalSpaceType >  GeneralizedNewmarkGN11SchemeType;
    typedef ExplicitForwardEulerScheme< SparseSpaceType, LocalSpaceType >  ExplicitForwardEulerSchemeType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    py::class_< GeneralizedNewmarkGN11SchemeType, typename GeneralizedNewmarkGN11SchemeType::Pointer, BaseSchemeType >
    (m, "GeneralizedNewmarkGN11Scheme")
    .def(py::init< double >());

    py::class_< ExplicitForwardEulerSchemeType, typename ExplicitForwardEulerSchemeType::Pointer, BaseSchemeType >
    (m, "ExplicitForwardEulerScheme")
    .def(py::init< double >());

}

}  // namespace Python.
} // Namespace Kratos
