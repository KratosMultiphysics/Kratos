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

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"


// Utility includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/residual_based_newmark_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"

namespace Kratos 
{
    namespace Testing 
    {
        /// Tests
       
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        
        typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
        typedef PointerVectorSet<Dof<double>, SetIdentityFunction<Dof<double>>, std::less<SetIdentityFunction<Dof<double>>::result_type>, std::equal_to<SetIdentityFunction<Dof<double>>::result_type>, Dof<double>* > DofsArrayType;
        
        static inline DofsArrayType BasicTestSchemeDisplacement(
            ModelPart& ModelPart,
            SchemeType::Pointer pScheme,
            std::vector< Dof<double>::Pointer >& DoF,
            const double DeltaTime
            )
        {
            ModelPart.SetBufferSize(3);
            
            ModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
            ModelPart.AddNodalSolutionStepVariable(VELOCITY);
            ModelPart.AddNodalSolutionStepVariable(ACCELERATION);
            
            auto pnode = ModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            
            pnode->AddDof(DISPLACEMENT_X);
            pnode->AddDof(DISPLACEMENT_Y);
            pnode->AddDof(DISPLACEMENT_Z);
            
            ModelPart.GetProcessInfo()[DELTA_TIME] = DeltaTime;
            
            DoF.reserve(3);
            DoF.push_back(pnode->pGetDof(DISPLACEMENT_X));
            DoF.push_back(pnode->pGetDof(DISPLACEMENT_Y));
            DoF.push_back(pnode->pGetDof(DISPLACEMENT_Z));
            
            // Set initial solution
            const array_1d<double, 3> zero_vector = ZeroVector(3);
            pnode->FastGetSolutionStepValue(DISPLACEMENT) = zero_vector;
            pnode->FastGetSolutionStepValue(DISPLACEMENT, 1) = zero_vector;
            pnode->FastGetSolutionStepValue(DISPLACEMENT, 2) = zero_vector;
            pnode->FastGetSolutionStepValue(VELOCITY) = zero_vector;
            pnode->FastGetSolutionStepValue(VELOCITY, 1) = zero_vector;
            pnode->FastGetSolutionStepValue(VELOCITY, 2) = zero_vector;
            pnode->FastGetSolutionStepValue(ACCELERATION) = zero_vector;
            pnode->FastGetSolutionStepValue(ACCELERATION, 1) = zero_vector;
            pnode->FastGetSolutionStepValue(ACCELERATION, 2) = zero_vector;
            
            DofsArrayType Doftemp;
            Doftemp.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
            {
                Doftemp.push_back( it->get() );
            }
            Doftemp.Sort();
            
			CompressedMatrix A(boost::numeric::ublas::zero_matrix<double>(3, 3));
            Vector Dx = ZeroVector(3);
            Vector b = ZeroVector(3);
            
            pScheme->Initialize(ModelPart);
            
            return Doftemp;
        }
     
        /** 
         * Checks if the Newmark scheme performs correctly the integration
         */
        
        KRATOS_TEST_CASE_IN_SUITE(DisplacementNewmarkSchemeTest, KratosCoreFastSuite)
        {
            Model current_model;

            constexpr double tolerance = 1e-6;
            
            ModelPart& model_part = current_model.CreateModelPart("Main");
            
            typedef ResidualBasedNewmarkDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedNewmarkDisplacementSchemeType;
            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedNewmarkDisplacementSchemeType() );
            
            const double delta_time = 1.0e-4;

            std::vector< Dof<double>::Pointer > DoF;
            DofsArrayType Doftemp = BasicTestSchemeDisplacement(model_part, pscheme, DoF, delta_time);
            
            CompressedMatrix A = boost::numeric::ublas::zero_matrix<double>(3, 3);
            Vector Dx = ZeroVector(3);
            Vector b = ZeroVector(3);
            
            double time = 0;
            
            Node<3>::Pointer pnode = model_part.pGetNode(1);
            
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - delta_time);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * delta_time);
            pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - delta_time);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * delta_time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - delta_time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * delta_time);
            
            pscheme->Initialize(model_part);
            
            const unsigned int number_iterations = 10;
            for (unsigned int iter = 0; iter < number_iterations; ++iter)
            {
                time += delta_time;
                
                model_part.CloneTimeStep(time);
                
                Dx[0] = std::cos(time) - std::cos(time - delta_time);
                
                pscheme->InitializeSolutionStep(model_part, A, Dx, b);
                pscheme->Update(model_part, Doftemp, A, Dx, b);
                
                const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);
                
//                 // Debug
//                 std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;
                
                KRATOS_CHECK_LESS_EQUAL(std::abs(x - std::cos(time)), tolerance);
                KRATOS_CHECK_LESS_EQUAL(std::abs(v + std::sin(time)), tolerance);
                KRATOS_CHECK_LESS_EQUAL(std::abs(a + std::cos(time)), tolerance);
            }
        }
        
        /** 
         * Checks if the Bossak scheme performs correctly the integration
         */
        
        KRATOS_TEST_CASE_IN_SUITE(DisplacementBossakSchemeTest, KratosCoreFastSuite)
        {
            Model current_model;

            constexpr double tolerance = 1e-6;
            
            ModelPart& model_part = current_model.CreateModelPart("Main");
            
            typedef ResidualBasedBossakDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementSchemeType;
            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedBossakDisplacementSchemeType() );
            
            const double delta_time = 1.0e-4;

            std::vector< Dof<double>::Pointer > DoF;
            DofsArrayType Doftemp = BasicTestSchemeDisplacement(model_part, pscheme, DoF, delta_time);
            
            CompressedMatrix A = boost::numeric::ublas::zero_matrix<double>(3, 3);
            Vector Dx = ZeroVector(3);
            Vector b = ZeroVector(3);
            
            double time = 0;
            
            Node<3>::Pointer pnode = model_part.pGetNode(1);
            
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - delta_time);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * delta_time);
            pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - delta_time);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * delta_time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - delta_time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * delta_time);
            
            pscheme->Initialize(model_part);
            
            const unsigned int number_iterations = 10;
            for (unsigned int iter = 0; iter < number_iterations; ++iter)
            {
                time += delta_time;
                
                model_part.CloneTimeStep(time);
                
                Dx[0] = std::cos(time) - std::cos(time - delta_time);
                
                pscheme->InitializeSolutionStep(model_part, A, Dx, b);
                pscheme->Update(model_part, Doftemp, A, Dx, b);
                
                const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);
                
//                 // Debug
//                 std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;
                
                KRATOS_CHECK_LESS_EQUAL(std::abs(x - std::cos(time)), tolerance);
                KRATOS_CHECK_LESS_EQUAL(std::abs(v + std::sin(time)), tolerance);
                KRATOS_CHECK_LESS_EQUAL(std::abs(a + std::cos(time)), tolerance);
            }
        }
        
        /** 
         * Checks if the BDF2 scheme performs correctly the integration
         */
        
        KRATOS_TEST_CASE_IN_SUITE(DisplacementBDF2SchemeTest, KratosCoreFastSuite)
        {
            Model current_model;
            
            constexpr double tolerance = 1e-6;
            
            ModelPart& model_part = current_model.CreateModelPart("Main");
            
            typedef ResidualBasedBDFDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBDFDisplacementSchemeType;
            SchemeType::Pointer pscheme = SchemeType::Pointer( new ResidualBasedBDFDisplacementSchemeType(2) );
            
            const double delta_time = 1.0e-4;

            std::vector< Dof<double>::Pointer > DoF;
            DofsArrayType Doftemp = BasicTestSchemeDisplacement(model_part, pscheme, DoF, delta_time);
            
            CompressedMatrix A = boost::numeric::ublas::zero_matrix<double>(3, 3);
            Vector Dx = ZeroVector(3);
            Vector b = ZeroVector(3);
            
            double time = 0;
            
            Node<3>::Pointer pnode = model_part.pGetNode(1);
            
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X) = std::cos(time);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = std::cos(time - delta_time);
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X, 2) = std::cos(time - 2.0 * delta_time);
            pnode->FastGetSolutionStepValue(VELOCITY_X) = -std::sin(time);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 1) = -std::sin(time - delta_time);
            pnode->FastGetSolutionStepValue(VELOCITY_X, 2) = -std::sin(time - 2.0 * delta_time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X) = -std::cos(time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 1) = -std::cos(time - delta_time);
            pnode->FastGetSolutionStepValue(ACCELERATION_X, 2) = -std::cos(time - 2.0 * delta_time);
            
            pscheme->Initialize(model_part);
            
            const unsigned int number_iterations = 10;
            for (unsigned int iter = 0; iter < number_iterations; ++iter)
            {
                time += delta_time;
                
                model_part.CloneTimeStep(time);
                
                Dx[0] = std::cos(time) - std::cos(time - delta_time);
                
                pscheme->InitializeSolutionStep(model_part, A, Dx, b);
                pscheme->Update(model_part, Doftemp, A, Dx, b);
                
                const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
                const double v = pnode->FastGetSolutionStepValue(VELOCITY_X);
                const double a = pnode->FastGetSolutionStepValue(ACCELERATION_X);
                
//                 // Debug
//                 std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;
                
                KRATOS_CHECK_LESS_EQUAL(std::abs(x - std::cos(time)), tolerance);
                KRATOS_CHECK_LESS_EQUAL(std::abs(v + std::sin(time)), tolerance);
                KRATOS_CHECK_LESS_EQUAL(std::abs(a + std::cos(time)), tolerance);
            }
        }
        
    } // namespace Testing
}  // namespace Kratos.

