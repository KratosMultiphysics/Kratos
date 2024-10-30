// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "contact_structural_mechanics_application_variables.h"
#include "tests/cpp_tests/contact_structural_mechanics_fast_suite.h"

// Processes
#include "custom_processes/assign_parent_element_conditions_process.h"

namespace Kratos::Testing 
{

/** 
* Checks the correct work of the AssignParentElementConditionsProcess
*/
KRATOS_TEST_CASE_IN_SUITE(AssignParentElementConditionsProcess1, KratosContactStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 3);
    
    Properties::Pointer p_prop = r_model_part.CreateNewProperties(0);

    // First we create the nodes
    r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);

    // Now we create the elements and conditions
    auto p_element = r_model_part.CreateNewElement("Element2D4N", 1, {{1,2,3,4}}, p_prop);
    auto p_cond1 = r_model_part.CreateNewCondition("LineCondition2D2N", 1, {{1,2}}, p_prop);
    auto p_cond2 = r_model_part.CreateNewCondition("LineCondition2D2N", 2, {{2,3}}, p_prop);
    auto p_cond3 = r_model_part.CreateNewCondition("LineCondition2D2N", 3, {{3,4}}, p_prop);
    auto p_cond4 = r_model_part.CreateNewCondition("LineCondition2D2N", 4, {{4,1}}, p_prop);

    // Call the process
    AssignParentElementConditionsProcess(r_model_part, r_model_part).Execute();

    // Check the parent element
    KRATOS_EXPECT_EQ(p_cond1->GetValue(PARENT_ELEMENT), p_element);
    KRATOS_EXPECT_EQ(p_cond2->GetValue(PARENT_ELEMENT), p_element);
    KRATOS_EXPECT_EQ(p_cond3->GetValue(PARENT_ELEMENT), p_element);
    KRATOS_EXPECT_EQ(p_cond4->GetValue(PARENT_ELEMENT), p_element);
}

}  // namespace Kratos::Testing.
