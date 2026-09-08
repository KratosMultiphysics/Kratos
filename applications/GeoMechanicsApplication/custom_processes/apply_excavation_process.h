// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Lorenzo Gracia,
//                   Aron Noordam,
//                   Vahid Galavi,
//                   Marjan Fathian
//
#pragma once

#include "processes/process.h"

#include <string>

namespace Kratos
{
class Model;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyExcavationProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyExcavationProcess);

    ApplyExcavationProcess(Model& rModel, const Parameters& rProcessSettings);

    ~ApplyExcavationProcess() override = default;

    ApplyExcavationProcess(const ApplyExcavationProcess&)            = delete;
    ApplyExcavationProcess& operator=(const ApplyExcavationProcess&) = delete;

    void                      ExecuteInitialize() override;
    [[nodiscard]] std::string Info() const override;

private:
    std::vector<std::reference_wrapper<ModelPart>> mrModelParts;
    bool                                           mDeactivateSoilPart;
};

} // namespace Kratos