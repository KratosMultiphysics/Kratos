//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

#if !defined(IGA_SURFACE_UTILITIES_H_INCLUDED)
#define IGA_SURFACE_UTILITIES_H_INCLUDED

#include "includes/model_part.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

#include "utilities/math_utils.h"

namespace Kratos
{

namespace IgaSurfaceUtilties
{

    static void CalculateBaseVectors(
        const Element::GeometryType& rGeometry,
        const Matrix& rDN_De,
        array_1d<double, 3>& r_g1,
        array_1d<double, 3>& r_g2,
        array_1d<double, 3>& r_g3)
    {
        Matrix jacobian = ZeroMatrix(3, 2);
        IgaGeometryUtilities::CalculateJacobian(rGeometry, rDN_De, 3, 2, jacobian);

        //basis vectors g1 and g2
        if (r_g1.size() != 3)
            r_g1.resize(3, false);
        r_g1 = ZeroVector(3);
        if (r_g2.size() != 3)
            r_g2.resize(3, false);
        r_g2 = ZeroVector(3);
        if (r_g3.size() != 3)
            r_g3.resize(3, false);
        r_g3 = ZeroVector(3);

        r_g1[0] = jacobian(0, 0);
        r_g2[0] = jacobian(0, 1);
        r_g1[1] = jacobian(1, 0);
        r_g2[1] = jacobian(1, 1);
        r_g1[2] = jacobian(2, 0);
        r_g2[2] = jacobian(2, 1);

        MathUtils<double>::CrossProduct(r_g3, r_g1, r_g2);
    }

    static void CalculateBaseVector(
        const Element::GeometryType& rGeometry,
        const Matrix& rDN_De,
        array_1d<double, 3>& r_g3)
    {
        array_1d<double, 3> g1 = ZeroVector(3);
        array_1d<double, 3> g2 = ZeroVector(3);

        CalculateBaseVectors(rGeometry, rDN_De, g1, g2, r_g3);
    }
}

} // namespace Kratos

#endif // IGA_SURFACE_UTILITIES_H_INCLUDED
