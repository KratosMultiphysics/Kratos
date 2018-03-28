//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "includes/constitutive_law.h"

namespace Kratos 
{
    namespace Testing 
    {
        /**
        * Checks the correct work of the Has methods
        */

        KRATOS_TEST_CASE_IN_SUITE(TestHasMethods, KratosConstitutiveLawFastSuite)
        {
            ConstitutiveLaw this_cl = ConstitutiveLaw();

            KRATOS_CHECK_IS_FALSE(this_cl.Has(IS_RESTARTED));                 // Bool
            KRATOS_CHECK_IS_FALSE(this_cl.Has(DOMAIN_SIZE));                  // Integer
            KRATOS_CHECK_IS_FALSE(this_cl.Has(NODAL_H));                      // Double
            KRATOS_CHECK_IS_FALSE(this_cl.Has(DISPLACEMENT));                 // Array 1D
            KRATOS_CHECK_IS_FALSE(this_cl.Has(INITIAL_STRAIN));               // Vector
            KRATOS_CHECK_IS_FALSE(this_cl.Has(GREEN_LAGRANGE_STRAIN_TENSOR)); // Matrix
        }

        /**
        * Checks the correct work of the GetStrainMeasure method
        */

        KRATOS_TEST_CASE_IN_SUITE(TestGetStrainMeasureMethod, KratosConstitutiveLawFastSuite)
        {
            ConstitutiveLaw this_cl = ConstitutiveLaw();

            KRATOS_CHECK_EQUAL(this_cl.GetStrainMeasure(), ConstitutiveLaw::StrainMeasure_Infinitesimal);
        }

        /**
        * Checks the correct work of the ValidateInput method
        */

        KRATOS_TEST_CASE_IN_SUITE(TestValidateInputMethod, KratosConstitutiveLawFastSuite)
        {
            ConstitutiveLaw this_cl = ConstitutiveLaw();
            Properties prop = Properties();

            KRATOS_CHECK_IS_FALSE(this_cl.ValidateInput(prop));
        }

        /**
        * Checks the correct work of the IsIncremental method
        */

        KRATOS_TEST_CASE_IN_SUITE(TestIsIncrementalMethod, KratosConstitutiveLawFastSuite)
        {
            ConstitutiveLaw this_cl = ConstitutiveLaw();

            KRATOS_CHECK_IS_FALSE(this_cl.IsIncremental());
        }

        /**
        * Checks the correct work of the GetStressMeasure method
        */

        KRATOS_TEST_CASE_IN_SUITE(TestGetStressMeasureMethod, KratosConstitutiveLawFastSuite)
        {
            ConstitutiveLaw this_cl = ConstitutiveLaw();

            KRATOS_CHECK_EQUAL(this_cl.GetStressMeasure(), ConstitutiveLaw::StressMeasure_PK1);
        }

    } // namespace Testing
}  // namespace Kratos.
