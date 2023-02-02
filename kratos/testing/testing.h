//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#if !defined(KRATOS_TESTING_H_INCLUDED )
#define  KRATOS_TESTING_H_INCLUDED

#include "includes/kernel.h"
#include "includes/expects.h"  // It is almost always necessary. includes the exception
#include "includes/data_communicator.h"

#include "testing/test_suite.h" // This includes the test_case.h which includes tester.h
#include "testing/test_skipped_exception.h" // Macros and exception class used to skip tests.

namespace Kratos {
namespace Testing {

class KernelTest : public ::testing::Test 
{
    // We chose to initialize the mKernel in the constructor and not in the
    // fixture SetUp because this way we can perfectly control the life
    // cycle of the Kernel (as oposite of using a unique_ptr, for example) 
	protected:
        KernelTest() : mKernel() {}
        ~KernelTest() {}

	    Kratos::Kernel mKernel;
};

DataCommunicator& GetDefaultDataCommunicator();

}
}

#endif // KRATOS_TEST_SUITE_H_INCLUDED  defined
