// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Aron Noordam
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/kratos_parameters.h"

//linear solvers
#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "custom_solvers/pre_factorized_skyline_lu_factorization_solver.h"


namespace Kratos {
namespace Python {

void AddCustomSolversToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using ReordererType = Reorderer<SparseSpaceType, LocalSpaceType >;
    using DirectSolverType = DirectSolver<SparseSpaceType, LocalSpaceType, ReordererType >;
    using PreFactorizedSkylineLUFactorizationSolverType = PreFactorizedSkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType, ReordererType >;


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    py::class_<PreFactorizedSkylineLUFactorizationSolverType, PreFactorizedSkylineLUFactorizationSolverType::Pointer, DirectSolverType>(m, "PreFactorizedSkylineLUFactorizationSolver")
        .def(py::init< >())
        .def(py::init<Parameters>())
        .def("__str__", PrintObject<PreFactorizedSkylineLUFactorizationSolverType>)
        ;


}

} // Namespace Python.
} // Namespace Kratos
