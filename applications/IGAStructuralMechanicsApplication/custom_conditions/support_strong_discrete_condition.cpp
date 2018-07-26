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
#include "custom_conditions/support_strong_discrete_condition.h"
#include "utilities/math_utils.h"


#include "custom_utilities/iga_flags.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    void SupportStrongDiscreteCondition::GetBaseVectorsSurface(
        const Matrix& DN_De,
        Vector& g1,
        Vector& g2,
        Vector& g3)
    {
        Matrix J = ZeroMatrix(3, 2);
        CalculateJacobian(DN_De, J);

        //basis vectors g1 and g2
        if (g1.size() != 3)
            g1.resize(3, false);
        g1 = ZeroVector(3);
        if (g2.size() != 3)
            g2.resize(3, false);
        g2 = ZeroVector(3);
        if (g3.size() != 3)
            g3.resize(3, false);
        g3 = ZeroVector(3);

        g1[0] = J(0, 0);
        g2[0] = J(0, 1);
        g1[1] = J(1, 0);
        g2[1] = J(1, 1);
        g1[2] = J(2, 0);
        g2[2] = J(2, 1);

        MathUtils<double>::CrossProduct(g3, g1, g2);
    }
//************************************************************************************
//************************************************************************************
    void SupportStrongDiscreteCondition::CalculateRotation(const Matrix& ShapeFunctionDerivatives,
        Vector& Phi_r, array_1d<double, 2>& Phi)
    {
        KRATOS_TRY
        int number_of_points = ShapeFunctionDerivatives.size1();
        std::cout << "check 1" << std::endl;
        Vector Tangents = this->GetValue(TANGENTS);
        std::cout << "check 2" << std::endl;

        //basis vectors g1, g2 and g3 in current configuration
        Vector g1, g2, g3;
        GetBaseVectorsSurface(ShapeFunctionDerivatives, g1, g2, g3);
        std::cout << "check 3" << std::endl;

        // t1 normal to trim, t2 tangential to trim
        array_1d<double, 3> t2 = Tangents(0)*g1 + Tangents(1)*g2;
        array_1d<double, 3> t1 = ZeroVector(3);
        MathUtils<double>::CrossProduct(t1, t2, g3);

        t2 = t2 / norm_2(t2);
        t1 = t1 / norm_2(t1);

        // computation of the a3 displacement
        array_1d<double, 3> w = g3 - g3;
        array_1d<double, 3> SinusOmegaVector;
        MathUtils<double>::CrossProduct(SinusOmegaVector, g3, w);

        array_1d<double, 2> SinusOmega;
        SinusOmega(0) = inner_prod(SinusOmegaVector, t2);
        SinusOmega(1) = inner_prod(SinusOmegaVector, t1);

        if (SinusOmega(0) > 1.0)
            SinusOmega(0) = 0.999999;
        if (SinusOmega(1) > 1.0)
            SinusOmega(1) = 0.999999;
        Phi(0) = asin(SinusOmega(0));
        Phi(1) = asin(SinusOmega(1));

        //variation of the a3
        array_1d<double, 3> t3 = g3;
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
                MathUtils<double>::CrossProduct(SinusOmega_r, g3, t3_r);
                Phi_r(n * 3 + i) = 1.0 / sqrt(1.0 - pow(SinusOmega(0), 2))*inner_prod(SinusOmega_r, t2);
                // if needed at some point:
                //Phi_r_2(i * 3 + j) = 1.0 / sqrt(1.0 - pow(SinusOmega(1), 2))*inner_prod(SinusOmega_r, t1);
            }
        }
        KRATOS_CATCH("")
    }

    void SupportStrongDiscreteCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        const int number_of_control_points = GetGeometry().size();
        if (Is(IGAFlags::FIX_DISPLACEMENT_X))
        {
            std::cout << "des net" << std::endl;
            KRATOS_WATCH(IGAFlags::FIX_DISPLACEMENT_X)
        }
        if (IGAFlags::FIX_DISPLACEMENT_Y)
        {
            std::cout << "des net" << std::endl;
            KRATOS_WATCH(IGAFlags::FIX_DISPLACEMENT_Y)
        }
        if (IGAFlags::FIX_DISPLACEMENT_Z)
        {
            std::cout << "des net" << std::endl;
            KRATOS_WATCH(IGAFlags::FIX_DISPLACEMENT_Z)
        }


        if (Has(DISPLACEMENT))
        {
            const array_1d<double, 3>& displacement = this->GetValue(DISPLACEMENT);
            const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);

            KRATOS_WATCH(displacement)
            
            // DISPLACEMENTS
            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                if (N[i] > 1e-7)
                {
                    if (Is(IGAFlags::FIX_DISPLACEMENT_X))
                    {
                        GetGeometry()[i].SetValue(DISPLACEMENT_X, displacement[0]);
                        GetGeometry()[i].Fix(DISPLACEMENT_X);
                    }
                    if (Is(IGAFlags::FIX_DISPLACEMENT_Y))
                    {
                        GetGeometry()[i].SetValue(DISPLACEMENT_Y, displacement[1]);
                        GetGeometry()[i].Fix(DISPLACEMENT_Y);
                    }
                    if (Is(IGAFlags::FIX_DISPLACEMENT_Z))
                    {
                        GetGeometry()[i].SetValue(DISPLACEMENT_Z, displacement[2]);
                        GetGeometry()[i].Fix(DISPLACEMENT_Z);
                    }
                }
            }

            if (Has(ROTATION))
            {
                const array_1d<double, 3>& rotation = this->GetValue(ROTATION);
                const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

                for (unsigned int i = 0; i < number_of_control_points; i++)
                {
                    if (DN_De(i, 0) > 1e-7 || DN_De(i, 1) > 1e-7)
                    {
                        if (Is(IGAFlags::FIX_DISPLACEMENT_X))
                        {
                            GetGeometry()[i].SetValue(DISPLACEMENT_X, displacement[0]);
                            GetGeometry()[i].Fix(DISPLACEMENT_X);
                        }

                        if (Is(IGAFlags::FIX_DISPLACEMENT_Y))
                        {
                            GetGeometry()[i].SetValue(DISPLACEMENT_Y, displacement[1]);
                            GetGeometry()[i].Fix(DISPLACEMENT_Y);
                        }
                        if (Is(IGAFlags::FIX_DISPLACEMENT_Z))
                        {
                            GetGeometry()[i].SetValue(DISPLACEMENT_Z, displacement[2]);
                            GetGeometry()[i].Fix(DISPLACEMENT_Z);
                        }
                    }
                }
            }
        }
    }

//************************************************************************************
//************************************************************************************
    void SupportStrongDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
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

    }
} // Namespace Kratos


