//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"
#include "custom_strategies/custom_builders_and_solvers/block_residual_based_builder_and_solver.hpp"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergence_criteria/residual_criteria.hpp"
#include "custom_strategies/custom_convergence_criteria/displacement_criteria.hpp"

//schemes
#include "custom_strategies/custom_schemes/residual_based_static_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_newmark_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_bossak_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_rotation_bossak_scheme.hpp"


//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    //base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaBaseType;

    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;


    //custom types
    typedef ResidualBasedStaticScheme< SparseSpaceType, LocalSpaceType > ResidualBasedStaticSchemeType;
    typedef ResidualBasedNewmarkScheme< SparseSpaceType, LocalSpaceType > ResidualBasedNewmarkSchemeType;
    typedef ResidualBasedBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakSchemeType;

    typedef ResidualBasedRotationBossakScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedRotationBossakSchemeType;


    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // class_< TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
    // 	      bases< BaseSolvingStrategyType >,  boost::noncopyable >
    // ("TestStrategy",
    //  init<ModelPart&, LinearSolverType::Pointer, int, int, bool >() )
    // .def("MoveNodes",&TestStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
    // ;


    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    // Residual Based Builder and Solver
    typedef ResidualBasedBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBuilderAndSolverType;

    class_< ResidualBasedBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("ResidualBasedBuilderAndSolver", init< LinearSolverType::Pointer > ());


    // Block Residual Based Builder and Solver
    typedef BlockResidualBasedBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BlockResidualBasedBuilderAndSolverType;
    class_< BlockResidualBasedBuilderAndSolverType, bases<BuilderAndSolverType>, boost::noncopyable > ("BlockResidualBasedBuilderAndSolver", init< LinearSolverType::Pointer > ());


    //********************************************************************
    //*************************SHCHEME CLASSES****************************
    //********************************************************************

    class_< ResidualBasedStaticSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedStaticScheme", init< >() )

            .def("Initialize", &ResidualBasedStaticScheme<SparseSpaceType, LocalSpaceType>::Initialize)

            ;

    class_< ResidualBasedNewmarkSchemeType,
            bases< BaseSchemeType >, boost::noncopyable >
            (
                "ResidualBasedNewmarkScheme", init< double >() )

            .def("Initialize", &ResidualBasedNewmarkScheme<SparseSpaceType, LocalSpaceType>::Initialize)

            ;



    class_< ResidualBasedBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedBossakScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;


    class_< ResidualBasedRotationBossakSchemeType,
            bases< BaseSchemeType >,  boost::noncopyable >
            (
                "ResidualBasedRotationBossakScheme", init< double , double >() )

            .def("Initialize", &ResidualBasedRotationBossakScheme<SparseSpaceType, LocalSpaceType>::Initialize)
            ;


    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************


    class_< DisplacementConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ,
            bases< ConvergenceCriteriaBaseType >, boost::noncopyable >
            (
                "DisplacementConvergenceCriteria", init<double, double >()
            );



    class_< ResidualConvergenceCriteria< SparseSpaceType,  LocalSpaceType > ,
            bases< ConvergenceCriteriaBaseType >, boost::noncopyable >
            (
                "ResidualConvergenceCriteria", init<double, double >()
            );


}

}  // namespace Python.

} // Namespace Kratos

