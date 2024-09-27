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
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/contact_utilities.h"
#include "tests/test_utilities/cpp_tests_utilities.h"
#include "tests/cpp_tests/contact_structural_mechanics_fast_suite.h"

namespace Kratos::Testing
{
/**
* Checks the correct work of the contact utilities
*/
KRATOS_TEST_CASE_IN_SUITE(CheckModelPartHasRotationDoF, KratosContactStructuralMechanicsFastSuite)
{
    Model current_model;

    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    r_model_part.SetBufferSize(2);

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(ROTATION);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(DOMAIN_SIZE, 2);

    CppTestsUtilities::Create2DGeometry(r_model_part);

    VariableUtils().AddDof(DISPLACEMENT_X, r_model_part);
    VariableUtils().AddDof(DISPLACEMENT_Y, r_model_part);
    VariableUtils().AddDof(DISPLACEMENT_Z, r_model_part);

    // Check there is not rotation
    KRATOS_EXPECT_FALSE(ContactUtilities::CheckModelPartHasRotationDoF(r_model_part));

    VariableUtils().AddDof(ROTATION_X, r_model_part);
    VariableUtils().AddDof(ROTATION_Y, r_model_part);
    VariableUtils().AddDof(ROTATION_Z, r_model_part);

    // Check there is not rotation
    KRATOS_EXPECT_TRUE(ContactUtilities::CheckModelPartHasRotationDoF(r_model_part));
}
}  // namespace Kratos::Testing.
