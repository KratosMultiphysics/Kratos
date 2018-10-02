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

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "includes/amatrix_interface.h"


namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(AMatrixSum, KratosCoreFastSuite)
		{
			DenseVector<double> a{ 1.0, 2.0, 3.1 };
			KRATOS_CHECK_EQUAL(sum(a), 6.1);
		}

    }
}  // namespace Kratos.
