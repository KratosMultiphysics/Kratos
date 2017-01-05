//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

// Project includes
#include "testing/testing.h"
#include "utilities/brent_iteration.h"

namespace Kratos {
    namespace Testing {

        namespace Internals {
            double Brent_Test_Function1(double x)
            {
                return (x-2.0)*(x+1.0); // Roots at 2.0 and -1.0
            }
        }

        KRATOS_TEST_CASE_IN_SUITE(BrentIteration, KratosCoreFastSuite)
        {
            double iterTol = 1e-5;
            double maxIter = 10;
            double checkTol = 1e-3;

            // Standard operation
            KRATOS_CHECK_NEAR(BrentIteration::FindRoot(Internals::Brent_Test_Function1,1.0,3.0,iterTol,maxIter),2.0,checkTol);

            // Passing solution as input
            KRATOS_CHECK_NEAR(BrentIteration::FindRoot(Internals::Brent_Test_Function1,2.0,3.0,iterTol,maxIter),2.0,checkTol);
            KRATOS_CHECK_NEAR(BrentIteration::FindRoot(Internals::Brent_Test_Function1,1.0,2.0,iterTol,maxIter),2.0,checkTol);

            // Wrong input: initial guesses on the same side of root
            try {
                BrentIteration::FindRoot(Internals::Brent_Test_Function1,2.5,3.0,iterTol,maxIter);
                KRATOS_THROW_ERROR(Exception,"Test Failed: BrentIteration did not throw expected error","");
            }
            catch (Exception& e) {
                KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(std::string(e.what()), "Error in BrentIteration::FindRoot: The images of both initial guesses have the same sign.");
            }

        }
    }
}


