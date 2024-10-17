// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include <filesystem>
#include <functional>
#include <string>

#include "includes/kernel.h"
#include "includes/kratos_export_api.h"

#include "custom_utilities/process_factory.hpp"
#include "geo_mechanics_application.h"
#include "linear_solvers_application.h"
#include "structural_mechanics_application.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

class InputUtility;
class ProcessInfoParser;
class Process;
class TimeLoopExecutorInterface;
class TimeIncrementor;
class StrategyWrapper;

class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoSettlement
{
public:
    KratosGeoSettlement(std::unique_ptr<InputUtility>              pInputUtility,
                        std::unique_ptr<ProcessInfoParser>         pProcessInfoParser,
                        std::unique_ptr<TimeLoopExecutorInterface> pTimeLoopExecutorInterface);

    ~KratosGeoSettlement();

    int RunStage(const std::filesystem::path&            rWorkingDirectory,
                 const std::filesystem::path&            rProjectParametersFile,
                 const std::function<void(const char*)>& rLogCallback,
                 const std::function<void(double)>&      rReportProgress,
                 const std::function<void(const char*)>& rReportTextualProgress,
                 const std::function<bool()>&            rShouldCancel);

    const InputUtility* GetInterfaceInputUtility() const;

private:
    ModelPart& AddNewModelPart(const std::string& rModelPartName);
    void       PrepareModelPart(const Parameters& rSolverSettings);

    ModelPart& GetMainModelPart();
    ModelPart& GetComputationalModelPart();

    static void                           AddNodalSolutionStepVariablesTo(ModelPart& rModelPart);
    static void                           AddDegreesOfFreedomTo(ModelPart& rModelPart);
    void                                  InitializeProcessFactory();
    std::vector<std::shared_ptr<Process>> GetProcesses(const Parameters& project_parameters) const;
    static std::unique_ptr<TimeIncrementor> MakeTimeIncrementor(const Parameters& rProjectParameters);
    std::shared_ptr<StrategyWrapper> MakeStrategyWrapper(const Parameters& rProjectParameters,
                                                         const std::filesystem::path& rWorkingDirectory);
    LoggerOutput::Pointer            CreateLoggingOutput(std::stringstream& rKratosLogBuffer) const;
    void FlushLoggingOutput(const std::function<void(const char*)>& rLogCallback,
                            LoggerOutput::Pointer                   pLoggerOutput,
                            const std::stringstream&                rKratosLogBuffer) const;

    template <typename TVariableType>
    void RestoreValuesOfNodalVariable(const TVariableType& rVariable, Node::IndexType SourceIndex, Node::IndexType DestinationIndex)
    {
        if (!GetComputationalModelPart().HasNodalSolutionStepVariable(rVariable)) return;

        VariableUtils{}.SetHistoricalVariableToZero(rVariable, GetComputationalModelPart().Nodes());

        block_for_each(GetComputationalModelPart().Nodes(),
                       [&rVariable, SourceIndex, DestinationIndex](auto& node) {
            node.GetSolutionStepValue(rVariable, DestinationIndex) =
                node.GetSolutionStepValue(rVariable, SourceIndex);
        });
    }

    template <typename ProcessType>
    std::function<ProcessFactory::ProductType(const Parameters&)> MakeCreatorFor()
    {
        return [&model = mModel](const Parameters& rProcessSettings) {
            auto& model_part = model.GetModelPart(rProcessSettings["model_part_name"].GetString());
            return std::make_unique<ProcessType>(model_part, rProcessSettings);
        };
    }

    Kernel                                        mKernel;
    Model                                         mModel;
    std::string                                   mModelPartName;
    KratosGeoMechanicsApplication::Pointer        mpGeoApp;
    KratosLinearSolversApplication::Pointer       mpLinearSolversApp;
    KratosStructuralMechanicsApplication::Pointer mpStructuralMechanicsApp;
    std::unique_ptr<ProcessFactory>            mProcessFactory = std::make_unique<ProcessFactory>();
    std::unique_ptr<InputUtility>              mpInputUtility;
    std::unique_ptr<ProcessInfoParser>         mpProcessInfoParser;
    std::unique_ptr<TimeLoopExecutorInterface> mpTimeLoopExecutor;
    const std::string mComputationalSubModelPartName{"settlement_computational_model_part"};
};

} // namespace Kratos
