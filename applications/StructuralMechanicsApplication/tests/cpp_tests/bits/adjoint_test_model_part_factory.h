#if !defined(KRATOS_ADJOINT_CPP_TEST_MODEL_PART_FACTORY_H_INCLUDED)
#define KRATOS_ADJOINT_CPP_TEST_MODEL_PART_FACTORY_H_INCLUDED

#include "includes/model_part.h"

namespace Kratos
{
ModelPart& CreateStructuralMechanicsAdjointTestModelPart(ModelPart* pPrimalModelPart);
}

#endif