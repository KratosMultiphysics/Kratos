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
#include "includes/node.h"
#include "custom_utilities/node_surface_geometry_3d.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

#include "custom_utilities/anurbs.h"

#include "utilities/math_utils.h"

namespace Kratos
{

namespace IgaSurfaceUtilities
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

    static void CalculateInitialBaseVectors(
        const Element::GeometryType& rGeometry,
        const Matrix& rDN_De,
        array_1d<double, 3>& r_g1,
        array_1d<double, 3>& r_g2,
        array_1d<double, 3>& r_g3)
    {
        Matrix jacobian = ZeroMatrix(3, 2);
        IgaGeometryUtilities::CalculateInitialJacobian(rGeometry, rDN_De, 3, 2, jacobian);

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

    static bool ProjectNodeOnSurface(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
        const array_1d<double, 2>& rLocalCoordinates,
        const array_1d<double, 3>& rPoint,
        const double Accuracy,
        const double Tolerance,
        const double NumberOfIterations,
        array_1d<double, 2>& rNewLocalCoordinates
    )
    {
        return true;
        //return IgaSurfaceUtilities::NewtonRaphson(
        //    pSurface,
        //    rPoint,
        //    rLocalCoordinates[0],
        //    rLocalCoordinates[1],
        //    rNewLocalCoordinates[0],
        //    rNewLocalCoordinates[1],
        //    Accuracy,
        //    Tolerance,
        //    NumberOfIterations
        //);
    }

    static bool NewtonRaphson(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
        const array_1d<double, 3>& rPoint,
        const double InitialU,
        const double InitialV,
        double& rU,
        double& rV,
        const double Accuracy,
        const double Tolerance,
        const double MaxIterations
    )
    {
        rU = InitialU;
        rV = InitialV;

        //ScalarType eps1 = 0.00001;
        //ScalarType eps2 = 0.000005;

        //ScalarType minu = m_surface->DomainU().T0();
        //ScalarType maxu = m_surface->DomainU().T1();
        //ScalarType minv = m_surface->DomainV().T0();
        //ScalarType maxv = m_surface->DomainV().T1();

        for (int i = 0; i < MaxIterations; i++) {
            const auto s = pSurface->DerivativesAt(rU, rV, 2);

            const array_1d<double, 3> distance = s[0] - rPoint;

            const double c1v = norm_2(distance);

            if (c1v < Accuracy) {
                return true;
            }

            double s1_l = norm_2(s[1]);
            double s2_l = norm_2(s[2]);

            double c2an = std::abs(inner_prod(s[1], distance));
            double c2ad = s1_l * c1v;

            double c2bn = std::abs(inner_prod(s[2], distance));
            double c2bd = s2_l * c1v;

            double c2av = c2an / c2ad;
            double c2bv = c2bn / c2bd;

            if (c2av < Tolerance && c2bv < Tolerance) {
                return true;
            }

            double f = inner_prod(s[1], distance);
            double g = inner_prod(s[2], distance);

            double J_00 = inner_prod(s[1], s[1]) + inner_prod(s[3], distance);
            double J_01 = inner_prod(s[1], s[2]) + inner_prod(s[4], distance);
            double J_11 = inner_prod(s[2], s[2]) + inner_prod(s[5], distance);

            double det_J = J_00 * J_11 - J_01 * J_01;

            double d_u = (g * J_01 - f * J_11) / det_J;
            double d_v = (f * J_01 - g * J_00) / det_J;

            rU += d_u;
            rV += d_v;
        }

        return false;
    }
}

} // namespace Kratos

#endif // IGA_SURFACE_UTILITIES_H_INCLUDED
