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
    struct CriticalHeadInfo;
    struct CallBackFunctions;

    class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosExecute
    {
    public:
        KratosExecute();

        int ExecuteFlowAnalysis(std::string_view                       rWorkingDirectory,
                                const std::string&                       rProjectParamsFileName,
                                double                                   minCriticalHead,
                                double                                   maxCriticalHead,
                                double                                   stepCriticalHead,
                                std::string_view                       rCriticalHeadBoundaryModelPartName,
                                const std::function<void(const char*)>&  rLogCallback,
                                const std::function<void(double)>&       rReportProgress,
                                const std::function<void(const char*)>&  rReportTextualProgress,
                                const std::function<bool()>&             rShouldCancel);

        void ExecuteWithoutPiping(ModelPart& model_part,
                                  const Kratos::Parameters& gid_output_settings) const;

        int ExecuteWithPiping(ModelPart& model_part,
                              const Kratos::Parameters& gid_output_settings,
                              const CriticalHeadInfo& criticalHeadInfo,
                              LoggerOutput::Pointer p_output,
                              const CallBackFunctions& rCallBackFunctions);

        void AddNodalSolutionStepVariables(ModelPart& model_part) const;

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

        typedef GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>
            GeoMechanicsNewtonRaphsonErosionProcessStrategyType;

        int FindCriticalHead(ModelPart& model_part,
                               const Kratos::Parameters& gid_output_settings,
                               const CriticalHeadInfo& criticalHeadInfo,
                               LoggerOutput::Pointer p_output,
                               const shared_ptr<Process>& RiverBoundary,
                               const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer& p_solving_strategy,
                               const CallBackFunctions& rCallBackFunctions);

        void HandleCriticalHeadFound(const CriticalHeadInfo& criticalHeadInfo);

        void HandleCancellationAndReset(const std::function<void(const char*)>& rLogCallback,
                                    LoggerOutput::Pointer p_output);

        // Dof arrays
        typedef SetIdentityFunction<Dof<double>> result_type;
        typedef PointerVectorSet<Dof<double>, SetIdentityFunction<Dof<double>>, std::less<result_type>,
                                 std::equal_to<result_type>, Dof<double> *>
            DofsArrayType;

        static ConvergenceCriteriaType::Pointer setup_criteria_dgeoflow();
        static LinearSolverType::Pointer setup_solver_dgeoflow();
        static GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer setup_strategy_dgeoflow(ModelPart &model_part);
        void parseProcess(ModelPart &model_part, Parameters projFile);

    private:
        // Initial Setup
        Model current_model;
        Kernel kernel;
        KratosGeoMechanicsApplication::Pointer geoApp;
        std::string mWorkingDirectory;
        std::string mCriticalHeadBoundaryModelPartName;
        bool pipingSuccess = false;
        double criticalHead;
        double currentHead;
        bool exitLoop = false;
        std::vector<std::shared_ptr<Process>> mProcesses;

        void ResetModelParts();

        int echoLevel = 1;

        [[nodiscard]] int GetEchoLevel() const;

        void SetEchoLevel(int level);

        shared_ptr<Process> FindRiverBoundaryByName(const std::string& rCriticalHeadBoundaryModelPartName) const;

        shared_ptr<Process> FindRiverBoundaryAutomatically(const KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer& rpSolvingStrategy) const;

        int MainExecution(ModelPart&                                                          rModelPart,
                          const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer& rpSolvingStrategy,
                          double                                                              Time,
                          double                                                              DeltaTime,
                          unsigned int                                                        NumberOfIterations) const;
    };
}