//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "containers/flags.h"
#include "spaces/ublas_space.h"
#include "boost/numeric/ublas/matrix.hpp"

//---strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/MPM_strategy.h"

//---convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//---schemes
#include "custom_strategies/schemes/mpm_residual_based_bossak_scheme.hpp"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

//---builders and solvers
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

//---linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos{
namespace Python{

    namespace py = pybind11;

    void AddCustomStrategiesToPython(pybind11::module& m)
    {
        typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        //base types
        typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
        typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
        typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

        //custom strategy types
        typedef MPMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType,2> MPMStrategyType2D;
        typedef MPMStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType,3> MPMStrategyType3D;

        //custom scheme types
        typedef MPMResidualBasedBossakScheme< SparseSpaceType, LocalSpaceType >  MPMResidualBasedBossakSchemeType;

        // MPM Residual Based Bossak Scheme Type
        py::class_< MPMResidualBasedBossakSchemeType,typename MPMResidualBasedBossakSchemeType::Pointer, BaseSchemeType >(m,"MPMResidualBasedBossakScheme")
            .def(py::init < ModelPart&, unsigned int, unsigned int, double, double, bool>())
            .def("Initialize", &MPMResidualBasedBossakSchemeType::Initialize)
            ;

        // Strategy Type
        py::class_< MPMStrategyType2D,typename MPMStrategyType2D::Pointer, BaseSolvingStrategyType >(m,"MPM2D")
            .def(py::init< ModelPart&, ModelPart&, ModelPart&, LinearSolverType::Pointer, std::string, int, bool, bool, bool, bool>() )
            .def( "SetEchoLevel", &MPMStrategyType2D::SetEchoLevel)
            ;

        py::class_< MPMStrategyType3D,typename MPMStrategyType3D::Pointer, BaseSolvingStrategyType >(m,"MPM3D")
            .def(py::init< ModelPart&, ModelPart&, ModelPart&, LinearSolverType::Pointer, std::string, int, bool, bool, bool, bool>() )
            .def( "SetEchoLevel", &MPMStrategyType3D::SetEchoLevel)
            ;
    }

}  // namespace Python.
} // Namespace Kratos

