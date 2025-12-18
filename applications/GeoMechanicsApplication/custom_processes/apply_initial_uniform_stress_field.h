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

///
/// @brief This process applies an initial uniform stress field to all elements in a model part.
/// The elements in the model part need to be able to calculate and set the CAUCHY_STRESS_VECTOR variable.
/// The Parameters object should contain a "value" field, which is a vector representing the stress components.
/// The vector should have a length equal to the strain size (e.g. 4 for plane strain and axisymmetric cases, 6 for 3D).
/// Note that this means that if you want to apply a uniform stress field to
/// elements with different strain sizes, you will need to apply the process multiple times with separate model parts.
///
class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyInitialUniformStressField : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyInitialUniformStressField);

    ApplyInitialUniformStressField(ModelPart& rModelPart, const Parameters& rParameters);
    void                      ExecuteInitialize() override;
    [[nodiscard]] std::string Info() const override;

private:
    ModelPart& mrModelPart;
    Vector     mImposedStressVector;
};

} // namespace Kratos
