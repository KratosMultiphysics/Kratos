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

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "processes/compute_bdfcoefficients_process.h"

namespace Kratos {
  namespace Testing {

	KRATOS_TEST_CASE_IN_SUITE(ComputeBDFCoefficientsProcess, KratosCoreFastSuite)
	{
        const double tolerance = 1.0e-8;

        // Create test model part
        Model Model;
        auto &r_model_part = Model.CreateModelPart("ModelPart");

        // Set fake delta times to test
        r_model_part.CloneTimeStep(0.1);
        r_model_part.CloneTimeStep(0.3);

        // Set first order BDF coefficients process
        ComputeBDFCoefficientsProcess bdf_coefs_proc_1(r_model_part, 1);

        // Check initialization
        const auto &r_bdf_coefs_1 = r_model_part.GetProcessInfo()[BDF_COEFFICIENTS];
        KRATOS_CHECK_NEAR(r_bdf_coefs_1[0], 0.0, tolerance);
        KRATOS_CHECK_NEAR(r_bdf_coefs_1[1], 0.0, tolerance);

        // Check 1st order results
        bdf_coefs_proc_1.Execute();
        KRATOS_CHECK_NEAR(r_bdf_coefs_1[0],  5.0, tolerance);
        KRATOS_CHECK_NEAR(r_bdf_coefs_1[1], -5.0, tolerance);

        // Set second order BDF coefficients process
        ComputeBDFCoefficientsProcess bdf_coefs_proc_2(r_model_part, 2);

        // Check 2nd order results
        bdf_coefs_proc_2.Execute();
        const auto &r_bdf_coefs_2 = r_model_part.GetProcessInfo()[BDF_COEFFICIENTS];
        KRATOS_CHECK_NEAR(r_bdf_coefs_2[0], 25.0/3.0, tolerance);
        KRATOS_CHECK_NEAR(r_bdf_coefs_2[1], -15.0, tolerance);
        KRATOS_CHECK_NEAR(r_bdf_coefs_2[2], 20.0/3.0, tolerance);
    }

}
}  // namespace Kratos.
