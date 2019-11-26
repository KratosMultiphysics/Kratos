#if !defined(KRATOS_STRUCTURAL_MECHANICS_CPP_TESTS_BITS_TEST_MODEL_PART_FACTORY_H_INCLUDED)
#define KRATOS_STRUCTURAL_MECHANICS_CPP_TESTS_BITS_TEST_MODEL_PART_FACTORY_H_INCLUDED

#include "containers/model.h"
#include "includes/constitutive_law.h"
#include "includes/model_part.h"
#include <functional>

namespace Kratos
{
ModelPart& CreateStructuralMechanicsTestModelPart(Model* pModel,
                                                  const Element& rElementPrototype,
                                                  const ConstitutiveLaw& rCLPrototype,
                                                  std::function<void(ModelPart*)> CustomBC = {});
}

#endif