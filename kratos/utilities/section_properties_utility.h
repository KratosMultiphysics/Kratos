#pragma once

#include "includes/model_part.h"

namespace Kratos
{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SectionPropertiesUtility
{
public:

    static void InterpretSections(ModelPart& rModelPart);

private:

    static void ComputeRectangularSection(
        Properties& rProperties,
        const Vector& rParameters);
};

} // namespace Kratos
