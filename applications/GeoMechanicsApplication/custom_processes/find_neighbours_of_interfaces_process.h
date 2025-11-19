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

#include "includes/kratos_export_api.h"
#include "processes/process.h"

#include <functional>
#include <vector>

namespace Kratos
{

class Model;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) FindNeighboursOfInterfacesProcess : public Process
{
public:
    FindNeighboursOfInterfacesProcess(Model& rModel, const Parameters& rProcessSettings);
    ~FindNeighboursOfInterfacesProcess() override;
    void ExecuteInitialize() override;

    [[nodiscard]] std::string Info() const override;

private:
    std::vector<std::reference_wrapper<ModelPart>> mrModelParts;
    ModelPart&                                     mrMainModelPart;
};

} // namespace Kratos