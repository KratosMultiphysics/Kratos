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

#pragma once

#include "includes/kratos_export_api.h"
#include "processes/process.h"

namespace Kratos
{

class Model;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplySeepageBoundaryProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplySeepageBoundaryProcess);

    ApplySeepageBoundaryProcess() = default;
    ApplySeepageBoundaryProcess(Model& rModel, const Parameters& rProcessSettings);

    [[nodiscard]] std::string Info() const override;
};

} // namespace Kratos