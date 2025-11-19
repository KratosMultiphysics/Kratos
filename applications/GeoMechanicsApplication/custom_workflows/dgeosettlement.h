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

#include "custom_utilities/node_utilities.h"
#include "custom_utilities/process_factory.hpp"

namespace Kratos
{

class InputUtility;
class ProcessInfoParser;
class Process;
class TimeLoopExecutorInterface;
class TimeIncrementor;
class StrategyWrapper;
class KratosGeoMechanicsApplication;
class KratosStructuralMechanicsApplication;
class KratosLinearSolversApplication;

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
                 const std::function<void(double)>&      rProgressDelegate,
                 const std::function<void(const char*)>& rReportTextualProgress,
                 const std::function<bool()>&            rShouldCancel);

    const InputUtility* GetInterfaceInputUtility() const;

private:
    ModelPart& AddNewModelPart(const std::string& rModelPartName);
    void       PrepareModelPart(const Parameters& rSolverSettings);
    void AddProcessesSubModelPartList(const Parameters& rProjectParameters, Parameters& rSolverSettings);

    ModelPart& GetMainModelPart();
    ModelPart& GetComputationalModelPart();

    static void                           AddNodalSolutionStepVariablesTo(ModelPart& rModelPart);
    static void                           AddDegreesOfFreedomTo(ModelPart& rModelPart);
    void                                  InitializeProcessFactory();
    std::vector<std::shared_ptr<Process>> GetProcesses(const Parameters& project_parameters) const;
    static std::unique_ptr<TimeIncrementor> MakeTimeIncrementor(const Parameters& rProjectParameters);
    std::shared_ptr<StrategyWrapper> MakeStrategyWrapper(const Parameters& rProjectParameters,
                                                         const std::filesystem::path& rWorkingDirectory);
    static LoggerOutput::Pointer     CreateLoggingOutput(std::stringstream& rKratosLogBuffer);
    static void FlushLoggingOutput(const std::function<void(const char*)>& rLogCallback,
                                   LoggerOutput::Pointer                   pLoggerOutput,
                                   const std::stringstream&                rKratosLogBuffer);

    template <typename TVariableType>
    void ResetValuesOfNodalVariable(const TVariableType& rVariable)
    {
        if (!GetComputationalModelPart().HasNodalSolutionStepVariable(rVariable)) return;

        NodeUtilities::AssignUpdatedVectorVariableToNodes(GetComputationalModelPart().Nodes(),
                                                          rVariable, rVariable.Zero(), 0);
        NodeUtilities::AssignUpdatedVectorVariableToNodes(GetComputationalModelPart().Nodes(),
                                                          rVariable, rVariable.Zero(), 1);
    }

    template <typename ProcessType>
    std::function<ProcessFactory::ProductType(const Parameters&)> MakeCreatorFor()
    {
        return [&model = mModel](const Parameters& rProcessSettings) {
            auto& model_part = model.GetModelPart(rProcessSettings["model_part_name"].GetString());
            return std::make_unique<ProcessType>(model_part, rProcessSettings);
        };
    }

    template <typename ProcessType>
    std::function<ProcessFactory::ProductType(const Parameters&)> MakeCreatorWithModelFor()
    {
        return [&model = mModel](const Parameters& rProcessSettings) {
            return std::make_unique<ProcessType>(model, rProcessSettings);
        };
    }

    Kernel                                                mKernel;
    Model                                                 mModel;
    std::string                                           mModelPartName;
    std::shared_ptr<KratosGeoMechanicsApplication>        mpGeoApp;
    std::shared_ptr<KratosLinearSolversApplication>       mpLinearSolversApp;
    std::shared_ptr<KratosStructuralMechanicsApplication> mpStructuralMechanicsApp;
    std::unique_ptr<ProcessFactory>            mProcessFactory = std::make_unique<ProcessFactory>();
    std::unique_ptr<InputUtility>              mpInputUtility;
    std::unique_ptr<ProcessInfoParser>         mpProcessInfoParser;
    std::unique_ptr<TimeLoopExecutorInterface> mpTimeLoopExecutor;
    const std::string mComputationalSubModelPartName{"settlement_computational_model_part"};
};

} // namespace Kratos
