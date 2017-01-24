//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_4.h"
#include "tests/geometries/test_geometry.h"

namespace Kratos {
	namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(TestTetrahedra3D4, KratosCoreFastSuite)
		{
			Geometry<Point<3, double>> geom;

			// KRATOS_CHECK_EXCEPTION(geom.Area());
			// KRATOS_CHECK_EXCEPTION(geom.Volume());
		}

	}
}  // namespace Kratos.
