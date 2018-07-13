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
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"

// External includes

// Project includes
#include "custom_elements/membrane_discrete_element.h"

// Application includes
#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    void MembraneDiscreteElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        // definition of problem size
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
        }
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
        }

        //reading in of integration weight, shape function values and shape function derivatives
        const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables(3);
        CalculateConstitutiveVariables(actual_metric, constitutive_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

        Vector prestress_tensor = ZeroVector(3);
        CalculatePresstressTensor(prestress_tensor, actual_metric);
        constitutive_variables.S += prestress_tensor;

        // calculate B MATRICES
        Matrix BMembrane = ZeroMatrix(3, mat_size);
        CalculateBMembrane(BMembrane, actual_metric);

        // Nonlinear Deformation
        SecondVariations SecondVariationsStrain = SecondVariations(mat_size);
        CalculateSecondVariationStrainMembrane(SecondVariationsStrain, actual_metric);

        double weighting = integration_weight * mInitialMetric.detJ * GetProperties()[THICKNESS];

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            //adding membrane contributions to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, 
                BMembrane, 
                constitutive_variables.D,
                weighting);

            //adding non-linear-contribution to Stiffness-Matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                SecondVariationsStrain,
                constitutive_variables.S,
                weighting);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            // operation performed: rRightHandSideVector -= Weight*IntForce
            noalias(rRightHandSideVector) -= weighting * prod(trans(BMembrane), constitutive_variables.S);
        }
        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    void MembraneDiscreteElement::Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rVariable == VON_MISES_STRESS)
        {
            Vector stresses = ZeroVector(3);
            CalculateStresses(stresses, rCurrentProcessInfo);

            Vector principal_forces = ZeroVector(3);
            // Principal normal forces
            principal_forces[0] = 0.5*(stresses[0] + stresses[1]) + sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));
            principal_forces[1] = 0.5*(stresses[0] + stresses[1]) - sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));
            principal_forces[2] = sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));

            double vMises = sqrt(pow(principal_forces[0], 2) + pow(principal_forces[1], 2) - principal_forces[0] * principal_forces[1] + 3 * pow(principal_forces[2], 2));

            rOutput = vMises;
        }
        else
        {
            //SurfaceBaseDiscreteElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
        }
    }

    void MembraneDiscreteElement::Calculate(
        const Variable<Vector>& rVariable,
        Vector& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rVariable == STRESSES)
        {
            Vector stresses = ZeroVector(3);

            CalculateStresses(stresses, rCurrentProcessInfo);

            rValues = stresses;
        }
        else if (rVariable == PRINCIPAL_STRESSES)
        {
            Vector stresses = ZeroVector(3);
            CalculateStresses(stresses, rCurrentProcessInfo);

            Vector principal_stresses = ZeroVector(3);
            // Principal normal forces
            principal_stresses[0] = 0.5*(stresses[0] + stresses[1]) + sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));
            principal_stresses[1] = 0.5*(stresses[0] + stresses[1]) - sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));
            principal_stresses[2] = sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));

            rValues = principal_stresses;
        }
        else if (rVariable == PRINCIPAL_FORCES)
        {
            Vector stresses = ZeroVector(3);
            CalculateStresses(stresses, rCurrentProcessInfo);

            Vector principal_stresses = ZeroVector(3);
            // Principal normal forces
            principal_stresses[0] = 0.5*(stresses[0] + stresses[1]) + sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));
            principal_stresses[1] = 0.5*(stresses[0] + stresses[1]) - sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));
            principal_stresses[2] = sqrt(pow((stresses[0] - stresses[1]) / 2, 2) + pow(stresses[2], 2));

            const double thickness = GetValue(THICKNESS);

            rValues = principal_stresses * thickness;
        }
        else
        {
            SurfaceBaseDiscreteElement::Calculate(rVariable, rValues, rCurrentProcessInfo);
        }
    }

    void MembraneDiscreteElement::CalculateStresses(
        Vector& rStresses,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rStresses.size() != 5)
            rStresses.resize(5);
        noalias(rStresses) = ZeroVector(5); //resetting LHS
        

        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        MetricVariables actual_metric(3);
        CalculateMetric(actual_metric);
        ConstitutiveVariables constitutive_variables(3);
        CalculateConstitutiveVariables(actual_metric, constitutive_variables, Values, ConstitutiveLaw::StressMeasure_PK2);

        double detF = actual_metric.dA / mInitialMetric.dA;

        Vector prestress_tensor = ZeroVector(3);
        CalculatePresstressTensor(prestress_tensor, actual_metric);
        constitutive_variables.S += prestress_tensor;

        // PK2 normal force in local cartesian e1, e2
        rStresses[0] = constitutive_variables.S[0];
        rStresses[1] = constitutive_variables.S[1];
        rStresses[2] = constitutive_variables.S[2];
    }

    void MembraneDiscreteElement::CalculatePresstressTensor(
        Vector& rPrestressTensor,
        MetricVariables& rMetric
    )
    {
        rPrestressTensor.resize(3);

        const Vector prestress_variable = this->GetValue(MEMBRANE_PRESTRESS_TENSOR_PK2);
        const double thickness = this->GetValue(THICKNESS);

        rPrestressTensor = prestress_variable;
    }

    void MembraneDiscreteElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;
        const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        const double thickness = this->GetProperties().GetValue(THICKNESS);
        const double density = this->GetProperties().GetValue(DENSITY);
        const double mass = thickness * density * mInitialMetric.dA * integration_weight;

        const int number_of_control_points = N.size();
        const int mat_size = 3 * number_of_control_points;

        if (rMassMatrix.size1() != mat_size)
            rMassMatrix.resize(mat_size, mat_size, false);
        rMassMatrix = ZeroMatrix(mat_size, mat_size);

        for (int r = 0; r<number_of_control_points; r++)
        {
            for (int s = 0; s<number_of_control_points; s++)
            {
                rMassMatrix(3 * s,     3 * r)     = N(s)*N(r) * mass;
                rMassMatrix(3 * s + 1, 3 * r + 1) = rMassMatrix(3 * s, 3 * r);
                rMassMatrix(3 * s + 2, 3 * r + 2) = rMassMatrix(3 * s, 3 * r);
            }
        }
        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************
    void MembraneDiscreteElement::CalculateMetric(
        MetricVariables& rMetric
    )
    {
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        Jacobian(DN_De, rMetric.J);


        rMetric.g1[0] = rMetric.J(0, 0);
        rMetric.g2[0] = rMetric.J(0, 1);
        rMetric.g1[1] = rMetric.J(1, 0);
        rMetric.g2[1] = rMetric.J(1, 1);
        rMetric.g1[2] = rMetric.J(2, 0);
        rMetric.g2[2] = rMetric.J(2, 1);

        //basis vector g3
        MathUtils<double>::CrossProduct(rMetric.g3, rMetric.g1, rMetric.g2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.g3);
        //normal vector _n
        Vector n = rMetric.g3 / rMetric.dA;


        //GetCovariantMetric
        rMetric.gab[0] = pow(rMetric.g1[0], 2) + pow(rMetric.g1[1], 2) + pow(rMetric.g1[2], 2);
        rMetric.gab[1] = pow(rMetric.g2[0], 2) + pow(rMetric.g2[1], 2) + pow(rMetric.g2[2], 2);
        rMetric.gab[2] = rMetric.g1[0] * rMetric.g2[0] + rMetric.g1[1] * rMetric.g2[1] + rMetric.g1[2] * rMetric.g2[2];

        Hessian(rMetric.H, DDN_DDe);

        rMetric.curvature[0] = rMetric.H(0, 0)*n[0] + rMetric.H(1, 0)*n[1] + rMetric.H(2, 0)*n[2];
        rMetric.curvature[1] = rMetric.H(0, 1)*n[0] + rMetric.H(1, 1)*n[1] + rMetric.H(2, 1)*n[2];
        rMetric.curvature[2] = rMetric.H(0, 2)*n[0] + rMetric.H(1, 2)*n[1] + rMetric.H(2, 2)*n[2];


        //contravariant rMetric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.gab_con[0] = invdetGab*rMetric.gab[1];
        rMetric.gab_con[2] = -invdetGab*rMetric.gab[2];
        rMetric.gab_con[1] = invdetGab*rMetric.gab[0];


        array_1d<double, 3> g_con_1 = rMetric.g1*rMetric.gab_con[0] + rMetric.g2*rMetric.gab_con[2];
        array_1d<double, 3> g_con_2 = rMetric.g1*rMetric.gab_con[2] + rMetric.g2*rMetric.gab_con[1];


        //local cartesian coordinates
        double lg1 = norm_2(rMetric.g1);
        array_1d<double, 3> e1 = rMetric.g1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;

        //Matrix T_G_E = ZeroMatrix(3, 3);
        //Transformation matrix T from contravariant to local cartesian basis
        double eG11 = inner_prod(e1, rMetric.g1);
        double eG12 = inner_prod(e1, rMetric.g2);
        double eG21 = inner_prod(e2, rMetric.g1);
        double eG22 = inner_prod(e2, rMetric.g2);

        rMetric.Q = ZeroMatrix(3, 3);
        rMetric.Q(0, 0) = eG11*eG11;
        rMetric.Q(0, 1) = eG12*eG12;
        rMetric.Q(0, 2) = 2.0*eG11*eG12;
        rMetric.Q(1, 0) = eG21*eG21;
        rMetric.Q(1, 1) = eG22*eG22;
        rMetric.Q(1, 2) = 2.0*eG21*eG22;
        rMetric.Q(2, 0) = 2.0*eG11*eG21;
        rMetric.Q(2, 1) = 2.0*eG12*eG22;
        rMetric.Q(2, 2) = 2.0*eG11*eG22 + eG12*eG21;

        rMetric.T = ZeroMatrix(3, 3);
        rMetric.T(0, 0) = eG11*eG11;
        rMetric.T(0, 1) = eG21*eG21;
        rMetric.T(0, 2) = 2.0*eG11*eG21;
        rMetric.T(1, 0) = eG12*eG12;
        rMetric.T(1, 1) = eG22*eG22;
        rMetric.T(1, 2) = 2.0*eG12*eG22;
        rMetric.T(2, 0) = eG11*eG12;
        rMetric.T(2, 1) = eG21*eG22;
        rMetric.T(2, 2) = eG11*eG22 + eG12*eG21;
    }
    //************************************************************************************
    //************************************************************************************
    void MembraneDiscreteElement::CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        Vector strain_vector = ZeroVector(3);

        CalculateStrain(strain_vector, rActualMetric.gab, mInitialMetric.gab);
        rThisConstitutiveVariables.E = prod(mInitialMetric.Q, strain_vector);

        //Constitive Matrices DMembrane and DCurvature
        rValues.SetStrainVector(rThisConstitutiveVariables.E); //this is the input parameter
        rValues.SetStressVector(rThisConstitutiveVariables.S);    //this is an ouput parameter
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //this is an ouput parameter

        mConstitutiveLawVector[0]->CalculateMaterialResponse(rValues, ThisStressMeasure);

        //Local Cartesian Forces and Moments
        rThisConstitutiveVariables.S = prod(
            trans(rThisConstitutiveVariables.D), rThisConstitutiveVariables.E);
    }
} // Namespace Kratos


