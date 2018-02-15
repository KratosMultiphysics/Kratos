//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

/* External includes */

/* Project includes */
#include "testing/testing.h"

/* Utility includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

/* Element include */
#include "tests/test_element.h"

// Linear solvers
#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

// The most basic scheme (static)
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

// The most extended criteria (residual criteria)
#include "solving_strategies/convergencecriterias/residual_criteria.h"

// The most builder and solver (the block builder and solver)
#include "solving_strategies/convergencecriterias/residual_criteria.h"

// The strategies to test
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos 
{
    namespace Testing 
    {
        /// Tests
        // NOTE: The strategies test many things simulataneously
        // TODO: Create test for the other components
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        
        // The direct solver
        typedef Reorderer<SparseSpaceType,  LocalSpaceType > ReordererType;
        typedef DirectSolver<SparseSpaceType,  LocalSpaceType, ReordererType > DirectSolverType;
        typedef LinearSolver<SparseSpaceType,LocalSpaceType> LinearSolverType;
        typedef SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;
        
        // The convergence criteria type
        typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
        typedef ResidualCriteria< SparseSpaceType, LocalSpaceType > ResidualCriteriaType;
        
        // The builder ans solver type
        typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedBlockBuilderAndSolverType;
        
        // The time scheme
        typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
        typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeType;
        
        // The strategies
        typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
        typedef ResidualBasedLinearStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedLinearStrategyType;
        typedef ResidualBasedNewtonRaphsonStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > ResidualBasedNewtonRaphsonStrategyType;
        
        // Dof arrays
        typedef PointerVectorSet<Dof<double>, SetIdentityFunction<Dof<double>>, std::less<SetIdentityFunction<Dof<double>>::result_type>, std::equal_to<SetIdentityFunction<Dof<double>>::result_type>, Dof<double>* > DofsArrayType;
        
        typedef TestElement::ResidualType ResidualType;
        
        static inline void BasicTestStrategyDisplacement(
            ModelPart& ModelPart,
            const ResidualType ThisResidualType
            )
        {
            ModelPart.SetBufferSize(3);
            
            ModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            ModelPart.AddNodalSolutionStepVariable(VELOCITY);
            ModelPart.AddNodalSolutionStepVariable(ACCELERATION);
            
            NodeType::Pointer pnode = ModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            std::vector<NodeType::Pointer> geom(1);
            geom[0] = pnode;
            GeometryType::Pointer pgeom = Kratos::make_shared<GeometryType>(geom);
            Element::Pointer pelem = Kratos::make_shared<TestElement>(1, pgeom, ThisResidualType);
            ModelPart.AddElement(pelem);
            
            pnode->AddDof(DISPLACEMENT_X);
            pnode->AddDof(DISPLACEMENT_Y);
            pnode->AddDof(DISPLACEMENT_Z);
            
            // Set initial solution
            (pnode->FastGetSolutionStepValue(DISPLACEMENT)).clear();
            (pnode->FastGetSolutionStepValue(DISPLACEMENT, 1)).clear();
            (pnode->FastGetSolutionStepValue(DISPLACEMENT, 2)).clear();
        }
     
        /** 
         * Checks if the Linear strategy performs correctly the resolution of the system
         */
        
        KRATOS_TEST_CASE_IN_SUITE(DisplacementLinearStrategy, KratosCoreStrategiesFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            ModelPart model_part("Main");
            
            typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            typename LinearSolverType::Pointer psolver = typename LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            typename BuilderAndSolverType::Pointer pbuildandsolve = typename BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );
            
            typename SolvingStrategyType::Pointer pstrategy = typename SolvingStrategyType::Pointer( new ResidualBasedLinearStrategyType(model_part, pscheme, psolver, pbuildandsolve));

            BasicTestStrategyDisplacement(model_part, ResidualType::LINEAR);
            
            NodeType::Pointer pnode = model_part.pGetNode(1);
            
            double time = 0.0;
            const double delta_time = 1.0e-4;
            const unsigned int number_iterations = 1;
            for (unsigned int iter = 0; iter < number_iterations; ++iter) {
                time += delta_time;
       
                model_part.CloneTimeStep(time);
               
                pstrategy->SetEchoLevel(0);
                pstrategy->Solve();

                const array_1d<double, 3>& solution = pnode->FastGetSolutionStepValue(DISPLACEMENT);
           
//                 // Debug
//                 KRATOS_WATCH(solution)
                
                array_1d<double, 3> ref_solution;
                ref_solution[0] = 1.0;
                ref_solution[1] = 1.0;
                ref_solution[2] = 1.0;
                
                KRATOS_CHECK_LESS_EQUAL(norm_2(solution - ref_solution), tolerance);
            }
        }

        /** 
         * Checks if the Newton Rapshon strategy performs correctly the resolution of the system
         */
        
        KRATOS_TEST_CASE_IN_SUITE(DisplacementNRStrategy, KratosCoreStrategiesFastSuite) 
        {
            constexpr double tolerance = 1e-6;
            
            ModelPart model_part("Main");
            
            typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticSchemeType() );
            typename LinearSolverType::Pointer psolver = typename LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
            typename ConvergenceCriteriaType::Pointer pcriteria = typename ConvergenceCriteriaType::Pointer( new ResidualCriteriaType(1.0e-4, 1.0e-9) );
            typename BuilderAndSolverType::Pointer pbuildandsolve = typename BuilderAndSolverType::Pointer( new ResidualBasedBlockBuilderAndSolverType(psolver) );
            
            typename SolvingStrategyType::Pointer pstrategy = typename SolvingStrategyType::Pointer( new ResidualBasedNewtonRaphsonStrategyType(model_part, pscheme, psolver, pcriteria, pbuildandsolve));

            BasicTestStrategyDisplacement(model_part, ResidualType::NON_LINEAR);
            
            NodeType::Pointer pnode = model_part.pGetNode(1);
            
            double time = 0.0;
            const double delta_time = 1.0e-4;
            const unsigned int number_iterations = 1;
            for (unsigned int iter = 0; iter < number_iterations; ++iter) {
                time += delta_time;
       
                model_part.CloneTimeStep(time);
               
                array_1d<double, 3> init_vector;
                init_vector[0] = 0.5;
                init_vector[1] = 0.5;
                init_vector[2] = 0.5;
                pnode->FastGetSolutionStepValue(DISPLACEMENT) = init_vector;
                
                pcriteria->SetEchoLevel(0);
                pstrategy->SetEchoLevel(0);
                pstrategy->Solve();

                const array_1d<double, 3>& solution = pnode->FastGetSolutionStepValue(DISPLACEMENT);
           
//                 // Debug
//                 KRATOS_WATCH(solution)
                
                array_1d<double, 3> ref_solution;
                ref_solution[0] = 1.0;
                ref_solution[1] = 1.0;
                ref_solution[2] = 1.0;
                
                KRATOS_CHECK_LESS_EQUAL(norm_2(solution - ref_solution), tolerance);
            }
        }
        
    } // namespace Testing
}  // namespace Kratos.

