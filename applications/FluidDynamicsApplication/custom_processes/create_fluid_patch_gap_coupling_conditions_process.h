//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#pragma once
#include "containers/model.h"
#include "processes/process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos {
/** @brief Build a Q1 Stokes Gap-SBM layer without adding DOFs.
 * Run once after conservative activation and fluid material setup. Background
 * and inner meshes remain independent; all gap elements reuse active parents.
 * allow_transient opts into BDF2 with a fixed background. For motion, discard
 * this process' owned gap entities and rebuild after updating activation.
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
