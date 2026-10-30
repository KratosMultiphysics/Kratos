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
 * The stationary two_phase_stokes formulation instead takes the water and air
 * complete-element patches as background and inner, respectively. Their open
 * horizontal interface partitions must match and the required
 * interface_height setting locates it.
 * two_phase_qsvms reuses time-integrated QSVMS with zero BDF coefficients,
 * zero dynamic tau, fixed meshes and ASGS for a stationary solve.
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
