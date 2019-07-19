//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_TEST_UTILITY_H_INCLUDED )
#define  KRATOS_NURBS_TEST_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

class NurbsTestUtility
{
public:     // static methods
    static array_1d<double, 3> Point(const double X, const double Y,
        const double Z)
    {
        array_1d<double, 3> point;

        point[0] = X;
        point[1] = Y;
        point[2] = Z;

        return point;
    }

    static void ArrayAlmostEqual(const array_1d<double, 3>& rActual,
        const array_1d<double, 3>& rExpected, const double Tolerance = 1e-5)
    {
        for (size_t i = 0; i < 3; i++)
        {
            KRATOS_CHECK_NEAR(rActual[i], rExpected[i], Tolerance);
        }
    }
}; // class NurbsTestUtility

} // namespace Testing
} // namespace Kratos

#endif // KRATOS_NURBS_TEST_UTILITY_H_INCLUDED defined