// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

// System includes

/* External includes */

#include "dgeoapplication.h"

/* Utility includes */
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

// The most basic scheme (static)
#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"

// The most builder and solver (the block builder and solver)
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// The strategies to test
#include <custom_processes/apply_component_table_process.hpp>
#include <custom_processes/apply_constant_hydrostatic_pressure_process.hpp>
#include <linear_solvers/skyline_lu_factorization_solver.h>

#include <solving_strategies/convergencecriterias/mixed_generic_criteria.h>
#include <solving_strategies/strategies/implicit_solving_strategy.h>
#include <solving_strategies/strategies/residualbased_newton_raphson_strategy.h>

#include "custom_strategies/strategies/geo_mechanics_newton_raphson_erosion_process_strategy.hpp"

namespace Kratos
{
    class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoFlow: public KratosGeoApplication
    {
    public:
        KratosGeoFlow();
        ~KratosGeoFlow(){};

        int execute_application(std::string workingDirectory, std::string parameterName,
                                double minCriticalHead, double maxCriticalHead, double stepCriticalHead,
                                std::string criticalHeadBoundaryModelPartName,
                                std::function<void(char*)> logCallback,
                                std::function<void(double)> reportProgress,
                                std::function<void(char*)> reportTextualProgress,
                                std::function<bool()> shouldCancel);

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        // The direct solver
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType> SkylineLUFactorizationSolverType;

        // The convergence criteria type
        typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
        typedef MixedGenericCriteria<SparseSpaceType, LocalSpaceType> MixedGenericCriteriaType;
        typedef typename MixedGenericCriteriaType::ConvergenceVariableListType ConvergenceVariableListType;

        // The builder ans solver type
        typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverType;

        // The time scheme
        typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
        typedef BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> BackwardEulerQuasistaticPwSchemeType;

        // The strategies
        typedef ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>
            ImplicitSolvingStrategyType;

        typedef ResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>
            ResidualBasedNewtonRaphsonStrategyType;

        typedef GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>
            GeoMechanicsNewtonRaphsonErosionProcessStrategyType;

        // Dof arrays
        typedef SetIdentityFunction<Dof<double>> result_type;
        typedef PointerVectorSet<Dof<double>, SetIdentityFunction<Dof<double>>, std::less<result_type>,
                                 std::equal_to<result_type>, Dof<double> *>
            DofsArrayType;

        ConvergenceCriteriaType::Pointer convergence_criteria() override;
        LinearSolverType::Pointer solver_settings() override;
        ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer strategy_settings(ModelPart &model_part);

    private:

        shared_ptr<Process> FindRiverBoundaryByName(std::string criticalHeadBoundaryModelPartName,
                                                    std::vector<std::shared_ptr<Process>> processes);

        
    };
}