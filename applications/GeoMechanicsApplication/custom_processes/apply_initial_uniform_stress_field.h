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

///
/// @brief This process applies a uniform stress field to all elements in a model part.
/// The elements in the modelpart need to be able to calculate and set the CAUCHY_STRESS_VECTOR variable.
/// The Parameters object should contain a "value" field, which is a vector of size 6 representing the stress components.
///
class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyInitialUniformStressField : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyInitialUniformStressField);

    ApplyInitialUniformStressField(ModelPart& rModelPart, const Parameters& rParameters);
    void ExecuteInitialize() override;

private:
    ModelPart& mrModelPart;
    Vector     mImposedStressVector;
};

} // namespace Kratos
