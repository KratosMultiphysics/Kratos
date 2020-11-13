// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:	   BSD License
//				   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/variable_utils.h"
#include "utilities/cpp_tests_utilities.h"
#include "custom_utilities/contact_utilities.h"

namespace Kratos
{
    namespace Testing
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
            KRATOS_CHECK_IS_FALSE(ContactUtilities::CheckModelPartHasRotationDoF(r_model_part));

            VariableUtils().AddDof(ROTATION_X, r_model_part);
            VariableUtils().AddDof(ROTATION_Y, r_model_part);
            VariableUtils().AddDof(ROTATION_Z, r_model_part);

            // Check there is not rotation
            KRATOS_CHECK(ContactUtilities::CheckModelPartHasRotationDoF(r_model_part));
        }

    } // namespace Testing
}  // namespace Kratos.
