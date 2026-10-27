// License: BSD License (kratos/license.txt)
#pragma once
#include "containers/model.h"
#include "processes/process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos {
/** @brief Build a steady Q1 Stokes Gap-SBM layer without adding DOFs.
 * Run once after conservative activation and fluid material setup. Background
 * and inner meshes remain independent; all gap elements reuse active parents.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) CreateFluidPatchGapCouplingConditionsProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CreateFluidPatchGapCouplingConditionsProcess);
    CreateFluidPatchGapCouplingConditionsProcess(Model& rModel, Parameters Settings);
    void Execute() override;
    const Parameters GetDefaultParameters() const override;
private:
    Model& mrModel;
    Parameters mSettings;
};
} // namespace Kratos
