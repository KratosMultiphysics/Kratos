// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "processes/process.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

class ModelPart;
class Parameters;
class Element;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ResetDisplacementProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ResetDisplacementProcess);
    ResetDisplacementProcess(ModelPart& rModelPart, const Parameters&);
    ~ResetDisplacementProcess() override = default;

    ResetDisplacementProcess(const ResetDisplacementProcess&)            = delete;
    ResetDisplacementProcess& operator=(const ResetDisplacementProcess&) = delete;

    void ExecuteInitialize() override;
    int  Check() override;

private:
    ModelPart& mrModelPart;

    static void CheckRetrievedElementData(const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws,
                                          const std::vector<Vector>& rStressesOnIntegrationPoints,
                                          const Element&             element);
};

} // namespace Kratos