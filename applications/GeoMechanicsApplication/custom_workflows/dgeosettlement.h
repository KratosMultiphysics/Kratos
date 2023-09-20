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

#include "geo_mechanics_application.h"

namespace Kratos
{

class ProcessFactory;
class InputUtility;
class ProcessInfoParser;

class KRATOS_API(GEO_MECHANICS_APPLICATION) KratosGeoSettlement
{
public:
    explicit KratosGeoSettlement(std::unique_ptr<InputUtility> pInputUtility,
                                 std::unique_ptr<ProcessInfoParser> pProcessInfoParser);
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
    static void AddNodalSolutionStepVariablesTo(ModelPart& rModelPart);
    static void AddDegreesOfFreedomTo(ModelPart& rModelPart);
    void InitializeProcessFactory();

    Kernel mKernel;
    Model mModel;
    std::string mModelPartName;
    KratosGeoMechanicsApplication::Pointer mpGeoApp;
    std::unique_ptr<ProcessFactory> mProcessFactory;
    std::unique_ptr<InputUtility> mpInputUtility;
    std::unique_ptr<ProcessInfoParser> mpProcessInfoParser;
};

}
