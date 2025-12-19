// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/kratos_export_api.h"
#include "includes/smart_pointers.h"
#include "includes/table.h"
#include "processes/process.h"

#include <string>

namespace Kratos
{

class ModelPart;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyComponentTableProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyComponentTableProcess);

    using TableType = Table<double, double>;

    ApplyComponentTableProcess(ModelPart& rModelPart, Parameters ProcessSettings);
    ApplyComponentTableProcess(const ApplyComponentTableProcess&)            = delete;
    ApplyComponentTableProcess& operator=(const ApplyComponentTableProcess&) = delete;
    ~ApplyComponentTableProcess() override                                   = default;

    void                      ExecuteInitialize() override;
    void                      ExecuteInitializeSolutionStep() override;
    void                      ExecuteFinalize() override;
    [[nodiscard]] std::string Info() const override;

private:
    ModelPart&         mrModelPart;
    std::string        mVariableName;
    bool               mIsFixed;
    bool               mIsFixedProvided;
    double             mInitialValue;
    TableType::Pointer mpTable;
    double             mTimeUnitConverter;
};

} // namespace Kratos