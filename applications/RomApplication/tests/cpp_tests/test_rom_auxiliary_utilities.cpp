//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

// System includes
#include <algorithm>

/* External includes */

/* Project includes */
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_utilities/rom_auxiliary_utilities.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(RomAuxiliaryUtilitiesGetHRomElementChildrenIds, RomApplicationFastSuite)
{
    // Two triangles sharing the edge 2-3, with a line condition on each boundary edge
    //  4 ---- 3
    //  | \ e2 |
    //  |  \   |
    //  | e1 \ |
    //  1 ---- 2
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    auto p_prop = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    r_model_part.CreateNewElement("Element2D3N", 1, {1, 2, 4}, p_prop);
    r_model_part.CreateNewElement("Element2D3N", 2, {2, 3, 4}, p_prop);
    r_model_part.CreateNewCondition("LineCondition2D2N", 1, std::vector<IndexType>{1, 2}, p_prop); // child of element 1
    r_model_part.CreateNewCondition("LineCondition2D2N", 2, std::vector<IndexType>{2, 3}, p_prop); // child of element 2
    r_model_part.CreateNewCondition("LineCondition2D2N", 3, std::vector<IndexType>{3, 4}, p_prop); // child of element 2
    r_model_part.CreateNewCondition("LineCondition2D2N", 4, std::vector<IndexType>{4, 1}, p_prop); // child of element 1

    // Only element 1 in the HROM: its two conditions are returned
    std::map<std::string, std::map<IndexType, double>> hrom_weights;
    hrom_weights["Elements"][1] = 1.5;
    hrom_weights["Conditions"];
    auto children_ids = RomAuxiliaryUtilities::GetHRomElementChildrenIds(r_model_part, hrom_weights);
    std::sort(children_ids.begin(), children_ids.end());
    KRATOS_EXPECT_EQ(children_ids, std::vector<IndexType>({1, 4}));

    // Conditions already in the HROM are not returned again
    hrom_weights["Conditions"][4] = 0.0;
    children_ids = RomAuxiliaryUtilities::GetHRomElementChildrenIds(r_model_part, hrom_weights);
    KRATOS_EXPECT_EQ(children_ids, std::vector<IndexType>({1}));

    // Only element 2 in the HROM: condition 1 shares node 2 with it but is not one of its faces
    hrom_weights["Conditions"].clear();
    hrom_weights["Elements"].clear();
    hrom_weights["Elements"][2] = 0.0;
    children_ids = RomAuxiliaryUtilities::GetHRomElementChildrenIds(r_model_part, hrom_weights);
    std::sort(children_ids.begin(), children_ids.end());
    KRATOS_EXPECT_EQ(children_ids, std::vector<IndexType>({2, 3}));

    // No HROM elements, no children
    hrom_weights["Elements"].clear();
    children_ids = RomAuxiliaryUtilities::GetHRomElementChildrenIds(r_model_part, hrom_weights);
    KRATOS_EXPECT_EQ(children_ids.size(), 0);
}

} // namespace Kratos::Testing
