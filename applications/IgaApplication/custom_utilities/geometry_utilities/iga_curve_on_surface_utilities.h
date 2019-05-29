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
            Matrix jacobian = ZeroMatrix(3, 2);
            IgaGeometryUtilities::CalculateJacobian(rGeometry, rDN_De, 3, 2, jacobian);

            //basis vectors g1 and g2
            array_1d<double, 3> g1;
            array_1d<double, 3> g2;

            g1[0] = jacobian(0, 0);
            g2[0] = jacobian(0, 1);
            g1[1] = jacobian(1, 0);
            g2[1] = jacobian(1, 1);
            g1[2] = jacobian(2, 0);
            g2[2] = jacobian(2, 1);

            rTangentVector = g1 * rTangents[0] + g2 * rTangents[1];
        }

        static void CalculateNormal(
            const Element::GeometryType& rGeometry,
            const Matrix& rDN_De,
            const array_1d<double, 2>& rTangents,
            array_1d<double, 3>& rNormalVector)
        {
            Matrix jacobian = ZeroMatrix(3, 2);
            IgaGeometryUtilities::CalculateJacobian(rGeometry, rDN_De, 3, 2, jacobian);

            //basis vectors g1 and g2
            array_1d<double, 3> g1;
            array_1d<double, 3> g2;

            g1[0] = jacobian(0, 0);
            g2[0] = jacobian(0, 1);
            g1[1] = jacobian(1, 0);
            g2[1] = jacobian(1, 1);
            g1[2] = jacobian(2, 0);
            g2[2] = jacobian(2, 1);

            rNormalVector = g2 * rTangents[0] - g1 * rTangents[1];
        }

        static void CalculateVariationRotation(
            const Element::GeometryType& rGeometry,
            const Matrix& rDN_De,
            const array_1d<double, 2>& rTangents,
            const array_1d<double, 3>& rG10,
            const array_1d<double, 3>& rG20,
            const array_1d<double, 3>& rG30,
            const array_1d<double, 3>& rg1,
            const array_1d<double, 3>& rg2,
            const array_1d<double, 3>& rg3,
            array_1d<double, 2>& rPhi,
            Vector& rPhi_r,
            Matrix& rPhi_rs)
        {
            const int number_of_control_points = rGeometry.size();
            const int mat_size = number_of_control_points * 3;

            if (rPhi_r.size() != mat_size)
                rPhi_r.resize(mat_size, false);
            rPhi_r = ZeroVector(mat_size);
            if ((rPhi_rs.size1() != mat_size) && (rPhi_rs.size2() != mat_size))
                rPhi_rs.resize(mat_size, mat_size, false);
            rPhi_rs = ZeroMatrix(mat_size, mat_size);

            // T1 normal to trim, T2 tangential to trim
            array_1d<double, 3> T2 = rTangents(0)*rG10 + rTangents(1)*rG20;
            array_1d<double, 3> T1 = ZeroVector(3);
            MathUtils<double>::CrossProduct(T1, T2, rG30);

            T2 = T2 / norm_2(T2);
            T1 = T1 / norm_2(T1);

            // computation of the a3 displacement
            array_1d<double, 3> w = rg3 - rG30;
            array_1d<double, 3> sinus_omega_vector;
            MathUtils<double>::CrossProduct(sinus_omega_vector, rG30, w);

            array_1d<double, 2> sinus_omega;
            sinus_omega(0) = inner_prod(sinus_omega_vector, T2);
            sinus_omega(1) = inner_prod(sinus_omega_vector, T1);

            if (sinus_omega(0) > 1.0)
                sinus_omega(0) = 0.999999;
            if (sinus_omega(1) > 1.0)
                sinus_omega(1) = 0.999999;
            rPhi(0) = asin(sinus_omega(0));
            rPhi(1) = asin(sinus_omega(1));

            //variation of the a3
            array_1d<double, 3> t3 = rg3 / norm_2(rg3);
            array_1d<double, 3> tilde_t3; //g3
            MathUtils<double>::CrossProduct(tilde_t3, rg1, rg2);
            double Length_t3 = norm_2(tilde_t3);

            for (unsigned int n = 0; n < number_of_control_points; n++)
            {
                for (unsigned int i = 0; i < 3; i++)
                {
                    //variations of the basis vectors
                    array_1d<double, 3> a1_r = ZeroVector(3);
                    array_1d<double, 3> a2_r = ZeroVector(3);

                    a1_r(i) = rDN_De(n, 0);
                    a2_r(i) = rDN_De(n, 1);

                    //variation of the non normalized local vector
                    Vector a1_rxg2 = ZeroVector(3);
                    MathUtils<double>::CrossProduct(a1_rxg2, a1_r, rg2);
                    Vector g1xa2_r = ZeroVector(3);
                    MathUtils<double>::CrossProduct(g1xa2_r, rg1, a2_r);
                    array_1d<double, 3> tilde_3_r = a1_rxg2 + g1xa2_r;
                    double line_t3_r = inner_prod(t3, tilde_3_r);
                    array_1d<double, 3> t3_r = tilde_3_r / Length_t3 - line_t3_r * t3 / Length_t3;
                    array_1d<double, 3> sinus_omega_r = ZeroVector(3);
                    MathUtils<double>::CrossProduct(sinus_omega_r, rG30, t3_r);
                    rPhi_r(n * 3 + i) = 1.0 / sqrt(1.0 - pow(sinus_omega(0), 2))*inner_prod(sinus_omega_r, T2);
                    // if needed at some point:
                    //Phi_r_2(i * 3 + j) = 1.0 / sqrt(1.0 - pow(SinusOmega(1), 2))*inner_prod(SinusOmega_r, T1);
                }
            }
            for (unsigned int n = 0; n < number_of_control_points; n++)
            {
                for (unsigned int i = 0; i < 3; i++)
                {
                    //variations of the basis vectors
                    array_1d<double, 3> a1_r_n = ZeroVector(3);
                    array_1d<double, 3> a2_r_n = ZeroVector(3);

                    a1_r_n(i) = rDN_De(n, 0);
                    a2_r_n(i) = rDN_De(n, 1);

                    //variation of the non normalized local vector
                    Vector a1_r_nxg2 = ZeroVector(3);
                    MathUtils<double>::CrossProduct(a1_r_nxg2, a1_r_n, rg2);
                    Vector g1xa2_r_n = ZeroVector(3);
                    MathUtils<double>::CrossProduct(g1xa2_r_n, rg1, a2_r_n);
                    array_1d<double, 3> tilde_3_r_n = a1_r_nxg2 + g1xa2_r_n;
                    double line_t3_r_n = inner_prod(t3, tilde_3_r_n);
                    array_1d<double, 3> t3_r_n = tilde_3_r_n / Length_t3 - line_t3_r_n * t3 / Length_t3;
                    array_1d<double, 3> sinus_omega_r_n = ZeroVector(3);
                    MathUtils<double>::CrossProduct(sinus_omega_r_n, rG30, t3_r_n);

                    for (unsigned int m = 0; m < number_of_control_points; m++)
                    {
                        for (unsigned int j = 0; j < 3; j++)
                        {
                            //variations of the basis vectors
                            array_1d<double, 3> a1_r_m = ZeroVector(3);
                            array_1d<double, 3> a2_r_m = ZeroVector(3);

                            a1_r_m(j) = rDN_De(m, 0);
                            a2_r_m(j) = rDN_De(m, 1);

                            Vector a1_r_mxg2 = ZeroVector(3);
                            MathUtils<double>::CrossProduct(a1_r_mxg2, a1_r_m, rg2);
                            Vector g1xa2_r_m = ZeroVector(3);
                            MathUtils<double>::CrossProduct(g1xa2_r_m, rg1, a2_r_m);
                            array_1d<double, 3> tilde_3_r_m = a1_r_mxg2 + g1xa2_r_m;
                            double line_t3_r_m = inner_prod(t3, tilde_3_r_m);
                            array_1d<double, 3> t3_r_m = tilde_3_r_m / Length_t3 - line_t3_r_m * t3 / Length_t3;
                            array_1d<double, 3> sinus_omega_r_m;
                            MathUtils<double>::CrossProduct(sinus_omega_r_m, rG30, t3_r_m);

                            Vector a1_r_nxa2_r_m = ZeroVector(3);
                            MathUtils<double>::CrossProduct(a1_r_nxa2_r_m, a1_r_n, a2_r_m);
                            Vector a1_r_mxa2_r_n = ZeroVector(3);
                            MathUtils<double>::CrossProduct(a1_r_mxa2_r_n, a1_r_m, a2_r_n);
                            array_1d<double, 3> tilde_t3_rs = a1_r_nxa2_r_m + a1_r_mxa2_r_n;
                            double line_t3_rs = inner_prod(t3_r_m, tilde_3_r_n) + inner_prod(t3, tilde_t3_rs);
                            array_1d<double, 3> t3_rs = (tilde_t3_rs*Length_t3 - line_t3_r_m * tilde_3_r_n) / pow(Length_t3, 2)
                                - line_t3_rs * t3 / Length_t3 - line_t3_r_n * (t3_r_m * Length_t3 - line_t3_r_m * t3) / pow(Length_t3, 2);
                            array_1d<double, 3> sinus_omega_rs;
                            MathUtils<double>::CrossProduct(sinus_omega_rs, rG30, t3_rs);

                            rPhi_rs(n * 3 + i, m * 3 + j) = inner_prod(sinus_omega_rs, T2) / sqrt(1.0 - pow(sinus_omega(0), 2))
                                + inner_prod(sinus_omega_r_m, T2)*inner_prod(sinus_omega_r_n, T2)*sinus_omega(0) / pow(1.0
                                    - pow(sinus_omega(0), 2), 1.5);
                        }
                    }
                }
            }
        }
    }
} // namespace Kratos

#endif // IGA_CURVE_ON_SURFACE_UTILITIES_H_INCLUDED
