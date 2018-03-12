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

    } // namespace Testing
}  // namespace Kratos.
