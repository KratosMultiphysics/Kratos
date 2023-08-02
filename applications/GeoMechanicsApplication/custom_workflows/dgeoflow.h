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

#include <geo_mechanics_application.h>
#include "includes/kernel.h"

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
    class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosExecute
    {
    public:
        KratosExecute();
        ~KratosExecute(){};

        int execute_flow_analysis(std::string workingDirectory, std::string parameterName,
                                  double minCriticalHead, double maxCriticalHead, double stepCriticalHead,
                                  std::string criticalHeadBoundaryModelPartName,
                                  std::function<void(char*)> logCallback,
                                  std::function<void(double)> reportProgress,
                                  std::function<void(char*)> reportTextualProgress,
                                  std::function<bool()> shouldCancel);

        typedef Node NodeType;
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

        ConvergenceCriteriaType::Pointer setup_criteria_dgeoflow();
        LinearSolverType::Pointer setup_solver_dgeoflow();
        GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer setup_strategy_dgeoflow(ModelPart &model_part);
        void parseMaterial(Model &model, std::string filepath);

        Parameters openProjectParamsFile(std::string filepath);
        std::vector<std::shared_ptr<Process>> parseProcess(ModelPart &model_part, Parameters projFile);

    private:
        // Initial Setup
        Model current_model;
        Kernel kernel;
        KratosGeoMechanicsApplication::Pointer geoApp;
        
        void ResetModelParts();

        int echoLevel = 1;

        int GetEchoLevel();

        void SetEchoLevel(int level);

        shared_ptr<Process> FindRiverBoundaryByName(std::string criticalHeadBoundaryModelPartName,
                                                    std::vector<std::shared_ptr<Process>> processes);

        shared_ptr<Process> FindRiverBoundaryAutomatically(KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer p_solving_strategy,
                                                           std::vector<std::shared_ptr<Process>> processes);

        int mainExecution(ModelPart &model_part,
                          std::vector<std::shared_ptr<Process>> processes,
                          KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer p_solving_strategy,
                          double time, double delta_time, double number_iterations);
    };
}