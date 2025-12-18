// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) SetAbsorbingBoundaryParametersProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SetAbsorbingBoundaryParametersProcess);

    SetAbsorbingBoundaryParametersProcess(ModelPart& model_part, Parameters rParameters);
    SetAbsorbingBoundaryParametersProcess(const SetAbsorbingBoundaryParametersProcess&) = delete;
    SetAbsorbingBoundaryParametersProcess& operator=(const SetAbsorbingBoundaryParametersProcess&) = delete;
    ~SetAbsorbingBoundaryParametersProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    std::string Info() const override;

private:
    /// Member Variables
    ModelPart& mrModelPart;
    Vector     mAbsorbingFactors;
    double     mVirtualThickness;
};

} // namespace Kratos