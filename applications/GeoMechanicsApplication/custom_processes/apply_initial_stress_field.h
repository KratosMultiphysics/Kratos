// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors: Richard Faasse
//

#pragma once

#include "includes/kratos_export_api.h"
#include "processes/process.h"

namespace Kratos
{

class ModelPart;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyInitialStressField : public Process
{
public:
    ApplyInitialStressField(ModelPart& rModelPart, const Parameters&);

private:
    ModelPart& mrModelPart;
};

} // namespace Kratos
