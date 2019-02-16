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

#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it

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
		KRATOS_TEST_CASE_IN_SUITE(AMatrixSwapOfVectors, KratosCoreFastSuite)
		{
			DenseVector<int> a{ 1, 2, 3 };
			DenseVector<int> b;
			b.swap(a);
			KRATOS_CHECK_EQUAL(b[0], 1);
			KRATOS_CHECK_EQUAL(b[1], 2);
			KRATOS_CHECK_EQUAL(b[2], 3);
			KRATOS_CHECK_EQUAL(b.size(), 3);
			KRATOS_CHECK_EQUAL(a.size(), 0);
		}

    }
}  // namespace Kratos.
#endif // ifdef KRATOS_USE_AMATRIX
