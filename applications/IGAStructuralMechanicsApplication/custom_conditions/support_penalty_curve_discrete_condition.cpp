//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/support_penalty_curve_discrete_condition.h"
#include "utilities/math_utils.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

#include "custom_utilities/iga_flags.h"

namespace Kratos
{
    void SupportPenaltyCurveDiscreteCondition::Initialize()
    {
        KRATOS_TRY

        //Constitutive Law initialisation
        CurveBaseDiscreteCondition::Initialize();

        Vector g1 = ZeroVector(3);
        Vector g2 = ZeroVector(3);
        Vector g3 = ZeroVector(3);

        GetBaseVectorsSurface(GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES), g1, g2, g3);

        m_g10 = g1;
        m_g20 = g2;
        m_g30 = g3;

        KRATOS_CATCH("")
    }
//************************************************************************************
//************************************************************************************
    void SupportPenaltyCurveDiscreteCondition::CalculateRotation(const Matrix& ShapeFunctionDerivatives,
        Vector& Phi_r, Matrix& Phi_rs, array_1d<double, 2>& Phi)
    {
        KRATOS_TRY
        int number_of_points = ShapeFunctionDerivatives.size1();

        Vector Tangents = this->GetValue(TANGENTS);

        //basis vectors g1, g2 and g3 in current configuration
        Vector g1, g2, g3;
        GetBaseVectorsSurface(ShapeFunctionDerivatives, g1, g2, g3);

        // T1 normal to trim, T2 tangential to trim
        array_1d<double, 3> T2 = Tangents(0)*m_g10 + Tangents(1)*m_g20;
        array_1d<double, 3> T1 = ZeroVector(3);
        MathUtils<double>::CrossProduct(T1, T2, m_g30);

        T2 = T2 / norm_2(T2);
        T1 = T1 / norm_2(T1);

        // computation of the a3 displacement
        array_1d<double, 3> w = g3 - m_g30;
        array_1d<double, 3> SinusOmegaVector;
        MathUtils<double>::CrossProduct(SinusOmegaVector, m_g30, w);

        array_1d<double, 2> SinusOmega;
        SinusOmega(0) = inner_prod(SinusOmegaVector, T2);
        SinusOmega(1) = inner_prod(SinusOmegaVector, T1);

        if (SinusOmega(0) > 1.0)
            SinusOmega(0) = 0.999999;
        if (SinusOmega(1) > 1.0)
            SinusOmega(1) = 0.999999;
        Phi(0) = asin(SinusOmega(0));
        Phi(1) = asin(SinusOmega(1));

        //variation of the a3
        array_1d<double, 3> t3 = g3/norm_2(g3);
        array_1d<double, 3> tilde_t3; //g3
        MathUtils<double>::CrossProduct(tilde_t3, g1, g2);
        double Length_t3 = norm_2(tilde_t3);

        for (unsigned int n = 0; n < number_of_points; n++)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                //variations of the basis vectors
                array_1d<double, 3> a1_r = ZeroVector(3);
                array_1d<double, 3> a2_r = ZeroVector(3);

                a1_r(i) = ShapeFunctionDerivatives(n, 0);
                a2_r(i) = ShapeFunctionDerivatives(n, 1);

                //variation of the non normalized local vector
                Vector a1_rxg2 = ZeroVector(3);
                MathUtils<double>::CrossProduct(a1_rxg2, a1_r, g2);
                Vector g1xa2_r = ZeroVector(3);
                MathUtils<double>::CrossProduct(g1xa2_r, g1, a2_r);
                array_1d<double, 3> tilde_3_r = a1_rxg2 + g1xa2_r;
                double line_t3_r = inner_prod(t3, tilde_3_r);
                array_1d<double, 3> t3_r = tilde_3_r / Length_t3 - line_t3_r * t3 / Length_t3;
                array_1d<double, 3> SinusOmega_r = ZeroVector(3);
                MathUtils<double>::CrossProduct(SinusOmega_r, m_g30, t3_r);
                Phi_r(n * 3 + i) = 1.0 / sqrt(1.0 - pow(SinusOmega(0), 2))*inner_prod(SinusOmega_r, T2);
                // if needed at some point:
                //Phi_r_2(i * 3 + j) = 1.0 / sqrt(1.0 - pow(SinusOmega(1), 2))*inner_prod(SinusOmega_r, T1);
            }
        }
        for (unsigned int n = 0; n < number_of_points; n++)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                //variations of the basis vectors
                array_1d<double, 3> a1_r_n = ZeroVector(3);
                array_1d<double, 3> a2_r_n = ZeroVector(3);

                a1_r_n(i) = ShapeFunctionDerivatives(n, 0);
                a2_r_n(i) = ShapeFunctionDerivatives(n, 1);

                //variation of the non normalized local vector
                Vector a1_r_nxg2 = ZeroVector(3);
                MathUtils<double>::CrossProduct(a1_r_nxg2, a1_r_n, g2);
                Vector g1xa2_r_n = ZeroVector(3);
                MathUtils<double>::CrossProduct(g1xa2_r_n, g1, a2_r_n);
                array_1d<double, 3> tilde_3_r_n = a1_r_nxg2 + g1xa2_r_n;
                double line_t3_r_n = inner_prod(t3, tilde_3_r_n);
                array_1d<double, 3> t3_r_n = tilde_3_r_n / Length_t3 - line_t3_r_n * t3 / Length_t3;
                array_1d<double, 3> SinusOmega_r_n = ZeroVector(3);
                MathUtils<double>::CrossProduct(SinusOmega_r_n, m_g30, t3_r_n);

                for (unsigned int m = 0; m < number_of_points; m++)
                {
                    for (unsigned int j = 0; j < 3; j++)
                    {
                        //variations of the basis vectors
                        array_1d<double, 3> a1_r_m = ZeroVector(3);
                        array_1d<double, 3> a2_r_m = ZeroVector(3);

                        a1_r_m(j) = ShapeFunctionDerivatives(m, 0);
                        a2_r_m(j) = ShapeFunctionDerivatives(m, 1);

                        Vector a1_r_mxg2 = ZeroVector(3);
                        MathUtils<double>::CrossProduct(a1_r_mxg2, a1_r_m, g2);
                        Vector g1xa2_r_m = ZeroVector(3);
                        MathUtils<double>::CrossProduct(g1xa2_r_m, g1, a2_r_m);
                        array_1d<double, 3> tilde_3_r_m = a1_r_mxg2 + g1xa2_r_m;
                        double line_t3_r_m = inner_prod(t3, tilde_3_r_m);
                        array_1d<double, 3> t3_r_m = tilde_3_r_m / Length_t3 - line_t3_r_m * t3 / Length_t3;
                        array_1d<double, 3> SinusOmega_r_m;
                        MathUtils<double>::CrossProduct(SinusOmega_r_m, m_g30, t3_r_m);

                        Vector a1_r_nxa2_r_m = ZeroVector(3);
                        MathUtils<double>::CrossProduct(a1_r_nxa2_r_m, a1_r_n, a2_r_m);
                        Vector a1_r_mxa2_r_n = ZeroVector(3);
                        MathUtils<double>::CrossProduct(a1_r_mxa2_r_n, a1_r_m, a2_r_n);
                        array_1d<double, 3> tilde_t3_rs = a1_r_nxa2_r_m + a1_r_mxa2_r_n;
                        double line_t3_rs = inner_prod(t3_r_m, tilde_3_r_n) + inner_prod(t3, tilde_t3_rs);
                        array_1d<double, 3> t3_rs = (tilde_t3_rs*Length_t3 - line_t3_r_m * tilde_3_r_n) / pow(Length_t3, 2)
                            - line_t3_rs * t3 / Length_t3 - line_t3_r_n * (t3_r_m * Length_t3 - line_t3_r_m * t3) / pow(Length_t3, 2);
                        array_1d<double, 3> SinusOmega_rs;
                        MathUtils<double>::CrossProduct(SinusOmega_rs, m_g30, t3_rs);

                        Phi_rs(n * 3 + i, m * 3 + j) = inner_prod(SinusOmega_rs, T2) / sqrt(1.0 - pow(SinusOmega(0), 2))
                            + inner_prod(SinusOmega_r_m, T2)*inner_prod(SinusOmega_r_n, T2)*SinusOmega(0) / pow(1.0
                                - pow(SinusOmega(0), 2), 1.5);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }



//************************************************************************************
//************************************************************************************
    void SupportPenaltyCurveDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY

        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        //resizing the system in case it does not have the right size
        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

        //resizing as needed the RHS
        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size, false);
        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

        //Read in Data
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);

        double Penalty = GetProperties()[PENALTY_FACTOR];


        //ROTATIONS
        if (Has(ROTATION))
        {
            const array_1d<double, 3>& rotation = this->GetValue(ROTATION);

            Vector Phi_r = ZeroVector(mat_size);
            Matrix Phi_rs = ZeroMatrix(mat_size, mat_size);
            array_1d<double, 2> Phi = ZeroVector(2);

            CalculateRotation(DN_De, Phi_r, Phi_rs, Phi);

            for (unsigned int n = 0; n < mat_size; n++)
            {
                for (unsigned int m = 0; m < mat_size; m++)
                {
                    rLeftHandSideMatrix(n, m) += (Phi_r(n)*Phi_r(m) + Phi(0)*Phi_rs(n, m));
                }
                rRightHandSideVector(n) -= Phi(0)*Phi_r(n);
            }
        }

        if (Has(DISPLACEMENT))
        {
            // DISPLACEMENTS
            Matrix Stiffness = ZeroMatrix(3, mat_size);
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                if (Is(IGAFlags::FIX_DISPLACEMENT_X))
                {
                    Stiffness(0, 3 * i) = N[i];
                }
                if (Is(IGAFlags::FIX_DISPLACEMENT_Y))
                {
                    Stiffness(1, 3 * i + 1) = N[i];
                }
                if (Is(IGAFlags::FIX_DISPLACEMENT_Z))
                {
                    Stiffness(2, 3 * i + 2) = N[i];
                }
            }
            const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);

            Vector TDisplacements(mat_size);
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                int index = 3 * i;
                TDisplacements[index] = (disp[0] - displacement[0]);
                TDisplacements[index + 1] = (disp[1] - displacement[1]);
                TDisplacements[index + 2] = (disp[2] - displacement[2]);
            }

            noalias(rLeftHandSideMatrix) += prod(trans(Stiffness), Stiffness);
            noalias(rRightHandSideVector) -= prod(prod(trans(Stiffness), Stiffness), TDisplacements);
        }
        //noalias(rRightHandSideVector) -= prod(prod(trans(Stiffness), Stiffness), TDisplacements);

        Vector t2 = ZeroVector(3);
        CalculateBaseVector(t2, DN_De);

        double weighting = integration_weight * norm_2(t2);

        //MAPPING
        rLeftHandSideMatrix  *= weighting * Penalty;
        rRightHandSideVector *= weighting * Penalty;

        KRATOS_CATCH("")
    }
} // Namespace Kratos


