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

#include "includes/kernel.h"
#include <geo_mechanics_application.h>

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

    using NodeType        = Node;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

    // The direct solver
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType>;

    // The convergence criteria type
    using ConvergenceCriteriaType  = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;
    using MixedGenericCriteriaType = MixedGenericCriteria<SparseSpaceType, LocalSpaceType>;
    using ConvergenceVariableListType = typename MixedGenericCriteriaType::ConvergenceVariableListType;

    using GeoMechanicsNewtonRaphsonErosionProcessStrategyType =
        GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    static ConvergenceCriteriaType::Pointer setup_criteria_dgeoflow();
    static LinearSolverType::Pointer        setup_solver_dgeoflow();
    static GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer setup_strategy_dgeoflow(ModelPart& rModelPart);
    void ParseProcesses(ModelPart& rModelPart, Parameters projFile);

    struct CriticalHeadInfo {
        double minCriticalHead  = 0.0;
        double maxCriticalHead  = 0.0;
        double stepCriticalHead = 0.0;

        CriticalHeadInfo(double minCriticalHead, double maxCriticalHead, double stepCriticalHead)
            : minCriticalHead(minCriticalHead), maxCriticalHead(maxCriticalHead), stepCriticalHead(stepCriticalHead)
        {
        }
    };

    struct CallBackFunctions {
        std::function<void(const char*)> LogCallback;
        std::function<void(double)>      ReportProgress;
        std::function<void(const char*)> ReportTextualProgress;
        std::function<bool()>            ShouldCancel;

        CallBackFunctions(std::function<void(const char*)> LogCallback,
                          std::function<void(double)>      ReportProgress,
                          std::function<void(const char*)> ReportTextualProgress,
                          std::function<bool()>            ShouldCancel)
            : LogCallback(std::move(LogCallback)),
              ReportProgress(std::move(ReportProgress)),
              ReportTextualProgress(std::move(ReportTextualProgress)),
              ShouldCancel(std::move(ShouldCancel))
        {
        }
    };

    int ExecuteFlowAnalysis(std::string_view         WorkingDirectory,
                            const std::string&       rProjectParamsFileName,
                            const CriticalHeadInfo&  rCriticalHeadInfo,
                            std::string_view         CriticalHeadBoundaryModelPartName,
                            const CallBackFunctions& rCallBackFunctions);

    void ExecuteWithoutPiping(ModelPart&                rModelPart,
                              const Kratos::Parameters& rGidOutputSettings,
                              const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer pSolvingStrategy) const;

    int ExecuteWithPiping(ModelPart&                rModelPart,
                          const Kratos::Parameters& rGidOutputSettings,
                          const CriticalHeadInfo&   rCriticalHeadInfo,
                          LoggerOutput::Pointer     pOutput,
                          const CallBackFunctions&  rCallBackFunctions,
                          const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer pSolvingStrategy);

    void WriteCriticalHeadResultToFile() const;

    void AddNodalSolutionStepVariables(ModelPart& rModelPart) const;

    int FindCriticalHead(ModelPart&                 rModelPart,
                         const Kratos::Parameters&  rGidOutputSettings,
                         const CriticalHeadInfo&    rCriticalHeadInfo,
                         LoggerOutput::Pointer      pOutput,
                         const shared_ptr<Process>& pRiverBoundary,
                         const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer pSolvingStrategy,
                         const CallBackFunctions& rCallBackFunctions);

    void HandleCriticalHeadFound(const CriticalHeadInfo& rCriticalHeadInfo);

    void HandleCleanUp(const CallBackFunctions& rCallBackFunctions, LoggerOutput::Pointer pOutput);

private:
    // Initial Setup
    Model                                  mCurrentModel;
    Kernel                                 mKernel;
    KratosGeoMechanicsApplication::Pointer mpGeoApp;
    std::string                            mWorkingDirectory;
    std::string                            mCriticalHeadBoundaryModelPartName;
    bool                                   mPipingSuccess = false;
    double                                 mCriticalHead  = 0.0;
    double                                 mCurrentHead   = 0.0;
    std::vector<std::shared_ptr<Process>>  mProcesses;
    int                                    mEchoLevel = 1;

    void ResetModelParts();

    [[nodiscard]] int GetEchoLevel() const;

    void SetEchoLevel(int level);

    shared_ptr<Process> FindRiverBoundaryByName(const std::string& CriticalHeadBoundaryModelPartName) const;

    shared_ptr<Process> FindRiverBoundaryAutomatically(
        const KratosExecute::GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer rpSolvingStrategy) const;

    int MainExecution(ModelPart& rModelPart,
                      const GeoMechanicsNewtonRaphsonErosionProcessStrategyType::Pointer rpSolvingStrategy,
                      double       Time,
                      double       DeltaTime,
                      unsigned int NumberOfIterations) const;

    bool AreExceedingMaxCriticalHead(double CurrentHead, double MaxCriticalHead) const;
};
} // namespace Kratos