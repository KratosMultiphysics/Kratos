//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Alexander Müller
//                   Tobias Teschemacher
//

// System includes
// External includes
// Project includes
// Application includes
#include "custom_elements/shell_5p_element.h"
#include <numeric>
namespace Kratos
{
    ///@name Initialize Functions
    ///@{

    void Shell5pElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        const GeometryType& r_geometry = GetGeometry();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        // Prepare memory
        if (reference_Curvature.size() != r_number_of_integration_points)
            reference_Curvature.resize(r_number_of_integration_points);
        if (reference_TransShear.size() != r_number_of_integration_points)
            reference_TransShear.resize(r_number_of_integration_points);
        if (m_dA_vector.size() != r_number_of_integration_points)
            m_dA_vector.resize(r_number_of_integration_points);
        if (m_cart_deriv.size() != r_number_of_integration_points)
            m_cart_deriv.resize(r_number_of_integration_points);

        KinematicVariables kinematic_variables;

        m_N = GetGeometry().ShapeFunctionsValues(); //get shape functions

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
        {
            m_cart_deriv[point_number] = CalculateCartesianDerivatives(point_number);

            std::tie(kinematic_variables, std::ignore) = CalculateKinematics(point_number);

            reference_Curvature[point_number] = kinematic_variables.curvature;
            reference_TransShear[point_number] = kinematic_variables.transShear;
        }
        InitializeMaterial();
        KRATOS_CATCH("")
    }

    void Shell5pElement::InitializeMaterial()
    {
        KRATOS_TRY
            const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        mConstitutiveLawVector.resize(r_number_of_integration_points);

        std::fill(mConstitutiveLawVector.begin(), mConstitutiveLawVector.end(), GetProperties()[CONSTITUTIVE_LAW]->Clone());

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number)
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(m_N, point_number));

        CalculateSVKMaterialTangent();
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Assembly
    ///@{

    void Shell5pElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        Matrix BOperator(8, number_of_nodes * 5);

        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {

            KinematicVariables kinematic_variables;
            VariationVariables variation_variables;
            std::tie(kinematic_variables, variation_variables) = CalculateKinematics(point_number);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);

            ConstitutiveVariables constitutive_variables(8); //0..2 membrane, 3..5 curvature, 6..7 transverse shear
            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);

            BOperator = CalculateStrainDisplacementOperator(point_number, kinematic_variables, variation_variables);

            double integration_weight = r_integration_points[point_number].Weight() * m_dA_vector[point_number];

            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag == true)
            {
                const Matrix Kg =
                    CalculateGeometricStiffness(
                        point_number,
                        kinematic_variables,
                        variation_variables,
                        constitutive_variables);
                noalias(rLeftHandSideMatrix) += integration_weight * (prod(prod<MatrixType>(trans(BOperator), mC), BOperator) + Kg);
            }
            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag == true)
                noalias(rRightHandSideVector) -= integration_weight * prod(trans(BOperator), constitutive_variables.StressVector);
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Kinematics
    ///@{

    std::pair< Shell5pElement::KinematicVariables, Shell5pElement::VariationVariables> Shell5pElement::CalculateKinematics(
        const IndexType IntegrationPointIndex
    ) const
    {
        KinematicVariables rKin;
        VariationVariables rVar;
        rKin.t = InterpolateNodalVariable(row(m_N, IntegrationPointIndex), m_GetValueFunctor, DIRECTOR);
        rKin.dtd1 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 0), m_GetValueFunctor, DIRECTOR);
        rKin.dtd2 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 1), m_GetValueFunctor, DIRECTOR);
        rKin.a1 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 0), m_GetCoordinatesFunctor);
        rKin.a2 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 1), m_GetCoordinatesFunctor);
        rKin.A1 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 0), m_GetInitialPositionFunctor);
        rKin.A2 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 1), m_GetInitialPositionFunctor);
        rKin.dud1 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 0), m_FastGetSolutionStepValueFunctor, DISPLACEMENT);
        rKin.dud2 = InterpolateNodalVariable(row(m_cart_deriv[IntegrationPointIndex], 1), m_FastGetSolutionStepValueFunctor, DISPLACEMENT);

        const double invL_t = 1.0 / norm_2(rKin.t);
        rKin.t *= invL_t;

        rKin.transShear[0] = inner_prod(rKin.t, rKin.a1);
        rKin.transShear[1] = inner_prod(rKin.t, rKin.a2);

        const BoundedMatrix<double, 3, 3> tdyadt = outer_prod(rKin.t, rKin.t);
        rVar.P = (IdentityMatrix(3) - tdyadt)* invL_t;

        const double txdtd1 = inner_prod(rKin.t, rKin.dtd1);
        const double txdtd2 = inner_prod(rKin.t, rKin.dtd2);
        const double a1xdtd1 = inner_prod(rKin.a1, rKin.dtd1);
        const double a2xdtd2 = inner_prod(rKin.a2, rKin.dtd2);
        const double a1xdtd2 = inner_prod(rKin.a1, rKin.dtd2);
        const double a2xdtd1 = inner_prod(rKin.a2, rKin.dtd1);
        const Matrix3d a1dyadt = outer_prod(rKin.a1, rKin.t);
        const Matrix3d a2dyadt = outer_prod(rKin.a2, rKin.t);
        const Matrix3d dtd1dyadt = outer_prod(rKin.dtd1, rKin.t);
        const Matrix3d dtd2dyadt = outer_prod(rKin.dtd2, rKin.t);
        const Matrix3d a1dyaddtd1 = outer_prod(rKin.a1, rKin.dtd1);
        const Matrix3d a2dyaddtd1 = outer_prod(rKin.a2, rKin.dtd1);
        const Matrix3d a1dyaddtd2 = outer_prod(rKin.a1, rKin.dtd2);
        const Matrix3d a2dyaddtd2 = outer_prod(rKin.a2, rKin.dtd2);

        const double normwquadinv = invL_t * invL_t;
        rVar.Q1 = normwquadinv * (txdtd1 * (3.0 * tdyadt - IdentityMatrix(3)) - dtd1dyadt - trans(dtd1dyadt));
        rVar.Q2 = normwquadinv * (txdtd2 * (3.0 * tdyadt - IdentityMatrix(3)) - dtd2dyadt - trans(dtd2dyadt));
        rVar.S1 = normwquadinv * (rKin.transShear[0] * (3.0 * tdyadt - IdentityMatrix(3)) - a1dyadt - trans(a1dyadt));
        rVar.S2 = normwquadinv * (rKin.transShear[1] * (3.0 * tdyadt - IdentityMatrix(3)) - a2dyadt - trans(a2dyadt));

        const double normwcubinv = normwquadinv * invL_t;
        rVar.Chi11 = normwcubinv * (3.0 * txdtd1 * (a1dyadt + 0.5 * rKin.transShear[0] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * a1xdtd1 * tdyadt + rKin.transShear[0] * dtd1dyadt) - a1dyaddtd1 - a1xdtd1 * 0.5 * IdentityMatrix(3));
        rVar.Chi12Chi21 = normwcubinv * (3.0 * txdtd1 * (a2dyadt + 0.5 * rKin.transShear[1] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * a2xdtd1 * tdyadt + rKin.transShear[1] * dtd1dyadt) - a2dyaddtd1 - a2xdtd1 * 0.5 * IdentityMatrix(3));
        rVar.Chi12Chi21 += normwcubinv * (3.0 * txdtd2 * (a1dyadt + 0.5 * rKin.transShear[0] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * a1xdtd2 * tdyadt + rKin.transShear[0] * dtd2dyadt) - a1dyaddtd2 - a1xdtd2 * 0.5 * IdentityMatrix(3));
        rVar.Chi22 = normwcubinv * (3.0 * txdtd2 * (a2dyadt + 0.5 * rKin.transShear[1] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * a2xdtd2 * tdyadt + rKin.transShear[1] * dtd2dyadt) - a2dyaddtd2 - a2xdtd2 * 0.5 * IdentityMatrix(3));

        rVar.Chi11 = trans(rVar.Chi11) + rVar.Chi11;
        rVar.Chi12Chi21 = trans(rVar.Chi12Chi21) + rVar.Chi12Chi21;
        rVar.Chi22 = trans(rVar.Chi22) + rVar.Chi22;
        //up to there dtd_al corresponds to w_{,al}
        rKin.dtd1 = prod(rVar.P, rKin.dtd1);
        rKin.dtd2 = prod(rVar.P, rKin.dtd2);

        rKin.metricChange[0] = inner_prod(rKin.A1, rKin.dud1) + 0.5 * norm_2_square(rKin.dud1);
        rKin.metricChange[1] = inner_prod(rKin.A2, rKin.dud2) + 0.5 * norm_2_square(rKin.dud2);
        rKin.metricChange[2] = inner_prod(rKin.a1, rKin.a2);

        rKin.curvature[0] = inner_prod(rKin.a1, rKin.dtd1);
        rKin.curvature[1] = inner_prod(rKin.a2, rKin.dtd2);
        rKin.curvature[2] = inner_prod(rKin.a1, rKin.dtd2) + inner_prod(rKin.a2, rKin.dtd1);

        return std::make_pair(rKin, rVar);
    }

    /* Transforms derivatives to obtain cartesian quantities  */
    Matrix Shell5pElement::CalculateCartesianDerivatives(
        const IndexType iP
    )
    {
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(iP);
        array_1d<double, 3> a3;
        Matrix32d J0Cart;
        Matrix J0;
        GetGeometry().Jacobian(J0, iP);

        const array_1d<double, 3> a1 = column(J0, 0);
        const array_1d<double, 3> a2 = column(J0, 1);

        m_dA_vector[iP] = norm_2(MathUtils<double>::CrossProduct(a1, a2));

        column(J0Cart, 0) = a1 / norm_2(a1); //Gram-Schmidt-Orthonormalization
        column(J0Cart, 1) = a2 - inner_prod(column(J0Cart, 0), a2) * column(J0Cart, 0);
        column(J0Cart, 1) = column(J0Cart, 1) / norm_2(column(J0Cart, 1));

        Matrix2d J = prod(trans(J0), J0Cart);

        double detJ;
        Matrix2d invJ;
        MathUtils<double>::InvertMatrix2(J, invJ, detJ);

        return prod(invJ, trans(r_DN_De));
    }

    void Shell5pElement::CalculateConstitutiveVariables(
        const IndexType iP,
        const KinematicVariables& rActualKinematic,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        subrange(rThisConstitutiveVariables.StrainVector, 0, 3) = rActualKinematic.metricChange;
        subrange(rThisConstitutiveVariables.StrainVector, 3, 6) = rActualKinematic.curvature - reference_Curvature[iP];
        subrange(rThisConstitutiveVariables.StrainVector, 6, 8) = rActualKinematic.transShear - reference_TransShear[iP];
        noalias(rThisConstitutiveVariables.StressVector) = prod(mC, rThisConstitutiveVariables.StrainVector);
    }


    Matrix Shell5pElement::CalculateStrainDisplacementOperator(
        const IndexType iP,
        const KinematicVariables& rActualKinematic,
        const VariationVariables& rVariations) const
    {
        const SizeType number_of_nodes = GetGeometry().size();

        Matrix rB{ ZeroMatrix(8, number_of_nodes * 5) };
        for (SizeType r = 0; r < number_of_nodes; r++)
        {
            const SizeType kr = 5 * r;

            const Matrix32d BLAI = GetGeometry()[r].GetValue(DIRECTORTANGENTSPACE);
            const Matrix3d WI1 = rVariations.Q1 * m_N(iP, r) + rVariations.P * m_cart_deriv[iP](0, r);
            const Matrix3d WI2 = rVariations.Q2 * m_N(iP, r) + rVariations.P * m_cart_deriv[iP](1, r);

            for (int s = 0; s < 3; s++)
            {
                rB(0, kr + s) = m_cart_deriv[iP](0, r) * rActualKinematic.a1[s]; //membrane_{,disp}
                rB(1, kr + s) = m_cart_deriv[iP](1, r) * rActualKinematic.a2[s];
                rB(2, kr + s) = m_cart_deriv[iP](0, r) * rActualKinematic.a2[s] + m_cart_deriv[iP](1, r) * rActualKinematic.a1[s];

                rB(3, kr + s) = m_cart_deriv[iP](0, r) * rActualKinematic.dtd1[s]; //bending_{,disp}
                rB(4, kr + s) = m_cart_deriv[iP](1, r) * rActualKinematic.dtd2[s];
                rB(5, kr + s) = m_cart_deriv[iP](0, r) * rActualKinematic.dtd2[s] + m_cart_deriv[iP](1, r) * rActualKinematic.dtd1[s];

                rB(6, kr + s) = m_cart_deriv[iP](0, r) * rActualKinematic.t[s];  //trans_shear_{,disp}
                rB(7, kr + s) = m_cart_deriv[iP](1, r) * rActualKinematic.t[s];
            }

            const Vector Temp = prod(prod<Vector>(trans(rActualKinematic.a1), WI1), BLAI); //bending_{,dir}
            const Vector Temp1 = prod(prod<Vector>(trans(rActualKinematic.a2), WI2), BLAI);
            const Vector Temp2 = prod(prod<Vector>(trans(rActualKinematic.a2), WI1) + prod(trans(rActualKinematic.a1), WI2), BLAI);
            const Vector Temp3 = prod(prod<Vector>(trans(rActualKinematic.a1), rVariations.P) * m_N(iP, r), BLAI); //shear_{,dir}
            const Vector Temp4 = prod(prod<Vector>(trans(rActualKinematic.a2), rVariations.P) * m_N(iP, r), BLAI);

            for (int s = 0; s < 2; s++)
            {
                rB(3, kr + 3 + s) = Temp(s);
                rB(4, kr + 3 + s) = Temp1(s);
                rB(5, kr + 3 + s) = Temp2(s);
                rB(6, kr + 3 + s) = Temp3(s);
                rB(7, kr + 3 + s) = Temp4(s);
            }
        }
        return rB;
    }

    Matrix Shell5pElement::CalculateGeometricStiffness(
        const IndexType iP,  //Integration Point
        const KinematicVariables& rActKin,
        const VariationVariables& ractVar,
        const ConstitutiveVariables& rConstitutive) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        Matrix Kg{ ZeroMatrix(number_of_control_points * 5, number_of_control_points * 5) };
        const array_1d<double, 8>& S = rConstitutive.StressVector;

        const Matrix3d chiAndSfac = ractVar.Chi11 * S[3] + ractVar.Chi22 * S[4] + ractVar.Chi12Chi21 * S[5]  //S[3..5] are moments
            + ractVar.S1 * S[6] + ractVar.S2 * S[7]; //S[6..7] is transverse shear
        for (SizeType i = 0; i < number_of_control_points; i++)
        {
            const SizeType i1 = 5 * i;
            const SizeType i4 = i1 + 3;

            const double Ni = m_N(iP, i);
            const double dN1i = m_cart_deriv[iP](0, i);
            const double dN2i = m_cart_deriv[iP](1, i);

            const Matrix3d WI1 = ractVar.Q1 * Ni + ractVar.P * dN1i;
            const Matrix3d WI2 = ractVar.Q2 * Ni + ractVar.P * dN2i;
            const Matrix23d BLAI_T = trans(r_geometry[i].GetValue(DIRECTORTANGENTSPACE));
            for (SizeType j = i; j < number_of_control_points; j++)
            {
                const SizeType j1 = 5 * j;
                const SizeType j4 = j1 + 3;

                const double Nj = m_N(iP, j);
                const double dN1j = m_cart_deriv[iP](0, j);
                const double dN2j = m_cart_deriv[iP](1, j);

                const Matrix3d WJ1 = ractVar.Q1 * Nj + ractVar.P * dN1j;
                const Matrix3d WJ2 = ractVar.Q2 * Nj + ractVar.P * dN2j;
                const Matrix32d BLAJ = r_geometry[j].GetValue(DIRECTORTANGENTSPACE);

                const double NS = dN1i * dN1j * S[0] + dN2i * dN2j * S[1] + (dN1i * dN2j + dN2i * dN1j) * S[2];
                Kg(i1, j1) = Kg(i1 + 1, j1 + 1) = Kg(i1 + 2, j1 + 2) = NS; // membrane_{,disp,disp}*N

                Matrix3d Temp = S[3] * dN1i * WJ1 + S[4] * dN2i * WJ2 + S[5] * (dN1i * WJ2 + dN2i * WJ1); // bending_{,dir,disp}*M
                Temp += ractVar.P * Nj * (dN1i * S[6] + dN2i * S[7]);  // shear_{,dir,disp}*Q

                noalias(subrange(Kg, i1, i1 + 3, j4, j4 + 2)) = prod(Temp, BLAJ);
                Temp = S[3] * dN1j * WI1 + S[4] * dN2j * WI2 + S[5] * (dN1j * WI2 + dN2j * WI1); // bending_{,disp,dir}*M

                Temp += ractVar.P * Ni * (dN1j * S[6] + dN2j * S[7]);// shear_{,disp,dir}*Q

                noalias(subrange(Kg, i4, i4 + 2, j1, j1 + 3)) = prod(BLAI_T, Temp);

                const double NdN1 = dN1j * Ni + Nj * dN1i;
                const double NdN2 = dN2j * Ni + Nj * dN2i;
                Temp = Ni * Nj * chiAndSfac;  // shear_{,dir,dir}*Q + bending_{,dir,dir}*M
                Temp += ractVar.S1 * NdN1 * S[3] + ractVar.S2 * NdN2 * S[4] + (ractVar.S1 * NdN2 + ractVar.S2 * NdN1) * S[5];  // bending_{,dir,dir}*M

                const Matrix23d Temp2 = prod(BLAI_T, Temp); //useless temp due to nonworking ublas prod(prod())
                noalias(subrange(Kg, i4, i4 + 2, j4, j4 + 2)) = prod(Temp2, BLAJ);
            }
            const auto& r_director = r_geometry[i].GetValue(DIRECTOR);
            const double kgT = -inner_prod(r_director /norm_2(r_director),
                prod(WI1, rActKin.a1) * S[3] +
                prod(WI2, rActKin.a2) * S[4] +
                (prod(WI2, rActKin.a1) + prod(WI1, rActKin.a2)) * S[5] +
                prod(ractVar.P, rActKin.a1) * Ni * S[6] +
                prod(ractVar.P, rActKin.a2) * Ni * S[7]); //P’_{,dir}*F_{int}
            Kg(i1 + 3, i1 + 3) += kgT;
            Kg(i1 + 4, i1 + 4) += kgT;
        }
        for (size_t i = 0; i < Kg.size1(); i++)
            for (size_t j = i + 1; j < Kg.size2(); j++)
                Kg(j, i) = Kg(i, j);

        return Kg;
    }

    ///@}
    ///@name Dynamic Functions
    ///@{

    void Shell5pElement::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 5;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (SizeType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const SizeType index = i * 3;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
            rValues[index + 2] = displacement[2];
        }
    }

    void Shell5pElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (SizeType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const SizeType index = i * 3;

            rValues[index] = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    void Shell5pElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (SizeType i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const SizeType index = i * 3;

            rValues[index] = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    void Shell5pElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;
        const SizeType number_of_control_points = GetGeometry().size();

        if (rResult.size() != 5 * number_of_control_points)
            rResult.resize(5 * number_of_control_points, false);

        const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const SizeType index = i * 5;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof(DIRECTORINC_X).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(DIRECTORINC_Y).EquationId();
        }
        KRATOS_CATCH("")
    }

    void Shell5pElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;
        const SizeType number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(5 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DIRECTORINC_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DIRECTORINC_Y));
        }

        KRATOS_CATCH("")
    };

    template< typename ContainerType, typename NodeFunctor, typename ...Args>
    BoundedVector<double, 3> Shell5pElement::InterpolateNodalVariable(const ContainerType& vec, const NodeFunctor& funct, const Args&... args) const
    {
        auto nodeValuesTimesAnsatzFunction = [&](double nodalVar, const NodeType& node)
        { return nodalVar * (node.*funct)(args...); };
        BoundedVector<double, 3> nullVec;
        nullVec = ZeroVector(3);
        return std::inner_product(vec.begin(), vec.end(), GetGeometry().begin(), nullVec, std::plus<BoundedVector<double, 3>>(), nodeValuesTimesAnsatzFunction);
    }

    void Shell5pElement::CalculateSVKMaterialTangent()
    {
        const double nu = this->GetProperties()[POISSON_RATIO];
        const double Emodul = this->GetProperties()[YOUNG_MODULUS];
        const double thickness = this->GetProperties().GetValue(THICKNESS);
        mC = ZeroMatrix(8, 8);

        const double fac1 = thickness * Emodul / (1 - nu * nu); //membrane
        mC(0, 0) = mC(1, 1) = fac1;
        mC(2, 2) = fac1 * (1 - nu) * 0.5;
        mC(1, 0) = mC(0, 1) = fac1 * nu;

        const double fac2 = thickness * thickness * fac1 / 12; // bending
        mC(3, 3) = mC(4, 4) = fac2;
        mC(5, 5) = fac2 * (1 - nu) * 0.5;
        mC(3, 4) = mC(4, 3) = fac2 * nu;

        const double fac3 = thickness * Emodul * 0.5 / (1 + nu); //trans shear
        mC(6, 6) = mC(7, 7) = fac3;
    }

    BoundedMatrix<double, 3, 2> Shell5pElement::TangentSpaceFromStereographicProjection(const array_1d<double, 3 >& director)
    {
        double st = (director[2] > 0) ? 1.0 : -1.0;
        double s = 1 / (1 + fabs(director[2]));

        const array_1d<double, 2 > y{ director[0] * s, director[1] * s };
        const double ys1 = y[0] * y[0];
        const double ys2 = y[1] * y[1];
        const double s2 = 2 * (1 + ys1 + ys2);

        BoundedMatrix<double, 3, 2> BLA;

        BLA(0, 0) = s2 - 4 * ys1;     BLA(0, 1) = -4 * y[0] * y[1];
        BLA(1, 0) = -4 * y[0] * y[1]; BLA(1, 1) = s2 - 4 * ys2;
        BLA(2, 0) = -st * 4 * y[0];   BLA(2, 1) = -st * 4 * y[1];

        const double normcol0 = norm_2(column(BLA, 0));
        const double normcol1 = norm_2(column(BLA, 1));
        column(BLA, 0) /= normcol0;
        column(BLA, 1) /= normcol1;
        return BLA;
    }

    void Shell5pElement::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
    {
        #pragma omp critical
        GetGeometry().GetGeometryParent(0).SetValue(DIRECTOR_COMPUTED, false);
    }

    void Shell5pElement::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
    {
        bool compute_director = false;

        #pragma omp critical
        if (!GetGeometry().GetGeometryParent(0).GetValue(DIRECTOR_COMPUTED))
        {
            GetGeometry().GetGeometryParent(0).SetValue(DIRECTOR_COMPUTED, true);
            compute_director = true;
        }
        if (compute_director) {
            auto& points = GetGeometry().GetGeometryParent(0).pGetGeometryPart(-1)->Points();
            for (auto& node : points)
            {
                const auto& r_director = node.GetValue(DIRECTOR);
                const double norm_t = norm_2(r_director);
                array_1d<double, 3 > director_new = r_director / norm_2(r_director);
                array_1d<double, 2 > inc2d;
                inc2d[0] = node.FastGetSolutionStepValue(DIRECTORINC_X) / norm_t;
                inc2d[1] = node.FastGetSolutionStepValue(DIRECTORINC_Y) / norm_t;

                const Matrix32d BLA = node.GetValue(DIRECTORTANGENTSPACE);
                const array_1d<double, 3 > inc3d = prod(BLA, inc2d);

                director_new += inc3d;
                director_new *= norm_t / norm_2(director_new);
                node.SetValue(DIRECTOR, director_new);
                node.FastGetSolutionStepValue(DIRECTORINC_X) = 0.0;
                node.FastGetSolutionStepValue(DIRECTORINC_Y) = 0.0;
                node.SetValue(DIRECTORTANGENTSPACE, TangentSpaceFromStereographicProjection(director_new));
            }
        }
    }
} // Namespace Kratos
