// BSD License, Kratos default license: kratos/license.txt
#pragma once

#include "includes/model_part.h"

namespace Kratos {
/// Materialize air-side scalar Gap-SBM entities in an auxiliary model part.
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) LaplacianGapGeometryUtility
{
public:
    static void Create(ModelPart& rPhase, ModelPart& rComputing, const Matrix& rContour,
                       Properties::Pointer pConditionProperties, double Tolerance);
};
}
