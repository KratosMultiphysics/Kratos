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

#if !defined(IGA_CURVE_ON_SURFACE_UTILITIES_H_INCLUDED)
#define IGA_CURVE_ON_SURFACE_UTILITIES_H_INCLUDED

#include "includes/model_part.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

namespace Kratos
{

    namespace IgaCurveOnSurfaceUtilities
    {
        /* returns the tangtial vector t2 oof the parameter curve
        *  on the surface.
        */
        static void CalculateTangent(
            const Element::GeometryType& rGeometry,
            const Matrix& rDN_De,
            const array_1d<double, 2>& rTangents,
            array_1d<double, 3>& rTangentVector)
        {
            Matrix Jacobian = ZeroMatrix(3, 2);
            IgaGeometryUtilities::CalculateJacobian(rGeometry, rDN_De, 3, 2, Jacobian);

            //basis vectors g1 and g2
            array_1d<double, 3> g1;
            array_1d<double, 3> g2;

            g1[0] = Jacobian(0, 0);
            g2[0] = Jacobian(0, 1);
            g1[1] = Jacobian(1, 0);
            g2[1] = Jacobian(1, 1);
            g1[2] = Jacobian(2, 0);
            g2[2] = Jacobian(2, 1);

            rTangentVector = g1 * rTangents[0] + g2 * rTangents[1];
        }
    }

} // namespace Kratos

#endif // IGA_CURVE_ON_SURFACE_UTILITIES_H_INCLUDED
