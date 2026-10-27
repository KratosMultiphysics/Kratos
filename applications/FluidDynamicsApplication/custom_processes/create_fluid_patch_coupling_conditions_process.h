// License: BSD License (kratos/license.txt)

#pragma once

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/**
 * @brief Create P1/P1 or Q1/Q1 fluid interface conditions between two independent meshes.
 * @details Extracts active boundary edges and intersects collinear edge pairs.
 * Nonmatching edge partitions are integrated on their common subsegments.
 * Patches and destination must belong to the same root ModelPart; no nodes
 * are merged or created. The destination must be empty. Execute once after
 * mesh creation and element activation; moving interfaces are not updated.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION)
    CreateFluidPatchCouplingConditionsProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CreateFluidPatchCouplingConditionsProcess);

    CreateFluidPatchCouplingConditionsProcess(Model& rModel, Parameters Settings);
    void Execute() override;
    const Parameters GetDefaultParameters() const override;
    std::string Info() const override { return "CreateFluidPatchCouplingConditionsProcess"; }

private:
    Model& mrModel;
    Parameters mSettings;
};

} // namespace Kratos
