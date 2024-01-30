//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/expect.h"
#include "utilities/entities_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(EntityIdentifierHasPrototype, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("MainModelPart");
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_tri_geom = r_model_part.CreateNewGeometry("Triangle2D3", 1, {{1, 2, 3}});
    auto entity_identifier = EntitiesUtilities::EntitityIdentifier<Element>("Element2D3N");
    KRATOS_EXPECT_TRUE(entity_identifier.HasPrototypeEntity(*p_tri_geom));
}

}  // namespace Kratos::Testing.
