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
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/shell_5p_element.h"



namespace Kratos
{
	///@name Initialize Functions
	///@{

	void Shell5pElement::Initialize()
	{
		KRATOS_TRY

			const GeometryType& r_geometry = GetGeometry();

		const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

		// Prepare memory
		if (m_B_ab_covariant_vector.size() != r_number_of_integration_points)
			m_B_ab_covariant_vector.resize(r_number_of_integration_points);
		if (m_dA_vector.size() != r_number_of_integration_points)
			m_dA_vector.resize(r_number_of_integration_points);
		if (m_cart_deriv.size() != r_number_of_integration_points)
			m_cart_deriv.resize(r_number_of_integration_points);

		KinematicVariables kinematic_variables(WorkingSpaceDimension());
		VariationVariables variation_variables;

		m_N = GetGeometry().ShapeFunctionsValues(); //get shape functions

		for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
		{
			CalculateCartesianDerivatives(point_number, kinematic_variables, m_cart_deriv[point_number]);

			CalculateKinematics(
				point_number,
				kinematic_variables, variation_variables);

			m_B_ab_covariant_vector[point_number] = kinematic_variables.b_ab_covariant;
			m_S_a_covariant_vector[point_number] = kinematic_variables.s_a_covariant;
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
		if (mConstitutiveLawVector.size() != r_number_of_integration_points)
			mConstitutiveLawVector.resize(r_number_of_integration_points);


		for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
			mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
			mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(m_N, point_number));
		}

		KRATOS_CATCH("");
	}

	///@}
	///@name Assembly
	///@{

	void Shell5pElement::CalculateAll(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool CalculateStiffnessMatrixFlag,
		const bool CalculateResidualVectorFlag
	)
	{
		KRATOS_TRY

			const auto& r_geometry = GetGeometry();

		// definition of problem size
		const SizeType number_of_nodes = r_geometry.size();
		const SizeType mat_size = number_of_nodes * 5;

		const auto& r_integration_points = r_geometry.IntegrationPoints();

		for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
			// Compute Kinematics and Metric
			KinematicVariables kinematic_variables(
				WorkingSpaceDimension());
			VariationVariables variation_variables;

			CalculateKinematics(
				point_number,
				kinematic_variables,
				variation_variables);

			// Create constitutive law parameters:
			ConstitutiveLaw::Parameters constitutive_law_parameters(
				GetGeometry(), GetProperties(), rCurrentProcessInfo);

			ConstitutiveVariables constitutive_variables(8); //0..2 membrane, 3..5 curvature, 6..7 transverse shear
			CalculateConstitutiveVariables(
				point_number,
				kinematic_variables,
				constitutive_variables_membrane,
				constitutive_variables_curvature,
				constitutive_variables_shear,
				constitutive_law_parameters,
				ConstitutiveLaw::StressMeasure_PK2);

			// calculate B operator
			Matrix BOperator = ZeroMatrix(8, mat_size);
			CalculateBOperator(
				point_number,
				BOperator,
				kinematic_variables,
				variation_variables);


			// Geometric stiffness matrix
			Matrix Kg(mat_size, mat_size);
			CalculateGeometricStiffness(
				point_number,
				Kg,
				kinematic_variables,
				variation_variables,
				constitutive_variables);

			double integration_weight =
				r_integration_points[point_number].Weight()
				* m_dA_vector[point_number]
				* GetProperties()[THICKNESS];

			// LEFT HAND SIDE MATRIX
			if (CalculateStiffnessMatrixFlag == true)
			{
				//adding membrane contributions to the stiffness matrix
				CalculateAndAddKm(
					rLeftHandSideMatrix,
					BMembrane,
					constitutive_variables_membrane.ConstitutiveMatrix,
					integration_weight);
				//adding curvature contributions to the stiffness matrix
				CalculateAndAddKm(
					rLeftHandSideMatrix,
					BCurvature,
					constitutive_variables_curvature.ConstitutiveMatrix,
					integration_weight);

				// adding  non-linear-contribution to Stiffness-Matrix
				CalculateAndAddNonlinearKm(
					rLeftHandSideMatrix,
					second_variations_strain,
					constitutive_variables_membrane.StressVector,
					integration_weight);

				CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
					second_variations_curvature,
					constitutive_variables_curvature.StressVector,
					integration_weight);
			}
			// RIGHT HAND SIDE VECTOR
			if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
			{
				// operation performed: rRightHandSideVector -= Weight*IntForce
				noalias(rRightHandSideVector) -= integration_weight * prod(trans(BMembrane), constitutive_variables_membrane.StressVector);
				noalias(rRightHandSideVector) -= integration_weight * prod(trans(BCurvature), constitutive_variables_curvature.StressVector);
			}
		}
		KRATOS_CATCH("");
	}

	///@}
	///@name Kinematics
	///@{

	void Shell5pElement::CalculateKinematics(
		IndexType IntegrationPointIndex,
		KinematicVariables& rKin,
		VariationVariables& rVar
	)
	{
		Matrix J;
		GetGeometry().Jacobian(J, IntegrationPointIndex);

		rKin.a1 = column(J, 0);
		rKin.a2 = column(J, 1);

		//GetCovariantMetric
		rKin.a_ab_covariant[0] = norm_2_square(rKin.a1); //TODO calculate metric from displacements and not from geometry. This is more accurate for small strains due cancelation with reference metric since 1.0-1.0 ~ loses digits 
		rKin.a_ab_covariant[1] = norm_2_square(rKin.a2);
		rKin.a_ab_covariant[2] = inner_prod(rKin.a1, rKin.a2);

		rKin.t = ZeroVector(3);
		rKin.dtd1 = ZeroVector(3);
		rKin.dtd2 = ZeroVector(3);

		const SizeType number_of_nodes = GetGeometry().size();
		for (size_t i = 0; i < number_of_nodes; i++)
		{
			rKin.t += m_N(IntegrationPointIndex, i) * GetGeometry()[i].FastGetSolutionStepValue(DIRECTOR); //TODO does this really do what i want? t=sum_I N_I *t_I ? and is this the best way to obtain this?
			rKin.dtd1 += m_cart_deriv[IntegrationPointIndex](0, i) * GetGeometry()[i].FastGetSolutionStepValue(DIRECTOR);
			rKin.dtd2 += m_cart_deriv[IntegrationPointIndex](1, i) * GetGeometry()[i].FastGetSolutionStepValue(DIRECTOR);
		}

		double invL_t = 1.0 / norm_2(rKin.t);
		rKin.t *= invL_t;

		array_1d<double, 3>& t = rKin.t; //define alias for cleaner code
		array_1d<double, 3>& dtd1 = rKin.dtd1; //define alias for cleaner code
		array_1d<double, 3>& dtd2 = rKin.dtd2; //define alias for cleaner code
		array_1d<double, 3>& a1 = rKin.a1; //define alias for cleaner code
		array_1d<double, 3>& a2 = rKin.a2; //define alias for cleaner code

		rVar.P = outer_prod(t, t) - IdentityMatrix(3);
		rVar.P *= invL_t;

		invL_t *= invL_t;
		const Matrix tdyadt = outer_prod(t, t);

		rKin.s_a_covariant[0] = inner_prod(t, a1);
		rKin.s_a_covariant[1] = inner_prod(t, a2);


		rVar.Q1 = invL_t * (inner_prod(t, dtd1) * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(dtd1, t) - outer_prod(t, dtd1));
		rVar.Q2 = invL_t * (inner_prod(t, dtd2) * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(dtd2, t) - outer_prod(t, dtd2));
		rVar.S1 = invL_t * (rKin.s_a_covariant[0] * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(a1, t) - outer_prod(t, a1));
		rVar.S2 = invL_t * (rKin.s_a_covariant[1] * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(a2, t) - outer_prod(t, a2));

		invL_t *= invL_t;
		rVar.Chi11 = invL_t * (3.0 * inner_prod(t, dtd1) * (outer_prod(a1, t) + 0.5 * rKin.s_a_covariant[0] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a1, dtd1) * tdyadt + rKin.s_a_covariant[0] * outer_prod(dtd1, t)) - outer_prod(a1, dtd1) - inner_prod(a1, dtd1) * 0.5 * IdentityMatrix(3));
		rVar.Chi21 = invL_t * (3.0 * inner_prod(t, dtd1) * (outer_prod(a2, t) + 0.5 * rKin.s_a_covariant[1] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a2, dtd1) * tdyadt + rKin.s_a_covariant[1] * outer_prod(dtd1, t)) - outer_prod(a2, dtd1) - inner_prod(a2, dtd1) * 0.5 * IdentityMatrix(3));
		rVar.Chi12 = invL_t * (3.0 * inner_prod(t, dtd2) * (outer_prod(a1, t) + 0.5 * rKin.s_a_covariant[0] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a1, dtd2) * tdyadt + rKin.s_a_covariant[0] * outer_prod(dtd2, t)) - outer_prod(a1, dtd2) - inner_prod(a1, dtd2) * 0.5 * IdentityMatrix(3));
		rVar.Chi22 = invL_t * (3.0 * inner_prod(t, dtd2) * (outer_prod(a2, t) + 0.5 * rKin.s_a_covariant[1] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a2, dtd2) * tdyadt + rKin.s_a_covariant[1] * outer_prod(dtd2, t)) - outer_prod(a2, dtd2) - inner_prod(a2, dtd2) * 0.5 * IdentityMatrix(3));

		//up to there dtd_al corresponds to w_al
		rKin.dtd1 = prod(rVar.P, dtd1);
		rKin.dtd2 = prod(rVar.P, dtd2);

		rKin.b_ab_covariant[0] = inner_prod(rKin.a1, dtd1);
		rKin.b_ab_covariant[1] = inner_prod(rKin.a2, dtd2);
		rKin.b_ab_covariant[2] = inner_prod(rKin.a1, dtd2) + inner_prod(a2, dtd1);
	}

	/* Transforms derivatives to obtain cartesian quantities  */
	void Shell5pElement::CalculateCartesianDerivatives(
		IndexType IntegrationPointIndex,
		const KinematicVariables& rKinematicVariables,
		Matrix& mCartD
	)
	{
		const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(IntegrationPointIndex);

		array_1d<double, 3> a3;
		MathUtils<double>::CrossProduct(a3, rKinematicVariables.a1, rKinematicVariables.a2);
		m_dA_vector[IntegrationPointIndex] = norm_2(a3);

		array_1d<double, 3> axi1 = rKinematicVariables.a1 / norm_2(rKinematicVariables.a1);
		array_1d<double, 3> axi2 = rKinematicVariables.a2 / norm_2(rKinematicVariables.a2);

		MathUtils<double>::CrossProduct(a3, axi1, axi2);

		array_1d<double, 3> axi1bar = 0.5 * (axi1 + axi2);
		axi1bar = axi1bar / norm_2(axi1bar);
		array_1d<double, 3> axi2bar;
		MathUtils<double>::CrossProduct(axi2bar, a3, axi1bar);

		array_1d<double, 3> a1Cart = .70710678118654752440 * (axi1bar - axi2bar);
		array_1d<double, 3> a2Cart = .70710678118654752440 * (axi1bar + axi2bar);

		Matrix J = ZeroMatrix(2, 2);

		J(0, 0) = inner_prod(rKinematicVariables.a1, a1Cart);
		J(1, 0) = inner_prod(rKinematicVariables.a2, a1Cart);
		J(0, 1) = inner_prod(rKinematicVariables.a1, a2Cart);
		J(1, 1) = inner_prod(rKinematicVariables.a2, a2Cart);

		double detJ;
		Matrix invJ;
		MathUtils<double>::InvertMatrix2(J, invJ, detJ);

		mCartD = prod(invJ, r_DN_De);
	}

	void Shell5pElement::CalculateConstitutiveVariables(
		IndexType IntegrationPointIndex,
		KinematicVariables& rActualKinematic,
		ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
		ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
		ConstitutiveVariables& rThisConstitutiveVariablesShear,
		ConstitutiveLaw::Parameters& rValues,
		const ConstitutiveLaw::StressMeasure ThisStressMeasure
	)
	{
		rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
		rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
		rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		array_1d<double, 3> strain_vector = (rActualKinematic.a_ab_covariant);
		strain_vector(0) -= 1.0; //TODO this is unnecessary if TODO from callculatekinematics is done
		strain_vector(1) -= 1.0;

		strain_vector(0) *= 0.5;
		strain_vector(1) *= 0.5;

		noalias(rThisConstitutiveVariablesMembrane.StrainVector) = strain_vector;

		array_1d<double, 3> curvature_vector = rActualKinematic.b_ab_covariant - m_B_ab_covariant_vector[IntegrationPointIndex];
		noalias(rThisConstitutiveVariablesCurvature.StrainVector) = curvature_vector;

		// Constitive Matrices DMembrane and DCurvature
		rValues.SetStrainVector(rThisConstitutiveVariablesMembrane.StrainVector); //this is the input parameter
		rValues.SetStressVector(rThisConstitutiveVariablesMembrane.StressVector);    //this is an ouput parameter
		rValues.SetConstitutiveMatrix(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix); //this is an ouput parameter

		mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);

		double thickness = this->GetProperties().GetValue(THICKNESS);
		noalias(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix) = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix * (pow(thickness, 2) / 12);   //TODO this does not work for general material laws especially including shear

		//Local Cartesian Forces and Moments
		noalias(rThisConstitutiveVariablesMembrane.StressVector) = prod(
			trans(rThisConstitutiveVariablesMembrane.ConstitutiveMatrix), rThisConstitutiveVariablesMembrane.StrainVector);
		noalias(rThisConstitutiveVariablesCurvature.StressVector) = prod(
			trans(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix), rThisConstitutiveVariablesCurvature.StrainVector);

		//TODO add shear constiutive relatons
	}

	void Shell5pElement::CalculateBOperator(
		IndexType IntegrationPointIndex,
		Matrix& rB,
		const KinematicVariables& rActualKinematic,
		const VariationVariables& rVariations)
	{
		const SizeType number_of_nodes = GetGeometry().size();
		const SizeType num_dof = number_of_nodes * 5;

		noalias(rB) = ZeroMatrix(8, num_dof);

		Matrix WI1(3, 3);
		Matrix WI2(3, 3);

		for (IndexType r = 0; r < number_of_nodes; r++)
		{
			// local node number kr and dof direction dirr
			IndexType kr = 5 * r;

			WI1 = rVariations.Q1 * m_N(IntegrationPointIndex, r) + rVariations.P * m_cart_deriv[IntegrationPointIndex](0, r);
			WI2 = rVariations.Q2 * m_N(IntegrationPointIndex, r) + rVariations.P * m_cart_deriv[IntegrationPointIndex](1, r);

			//membrane
			subrange(rB, 0, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.a1); //TODO WHY DOES THIS NOT WORK is  the equal sign not overloaded for array_1d?
			subrange(rB, 1, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.a2);
			subrange(rB, 2, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.a2) + m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.a1);

			//bending
			subrange(rB, 3, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.dtd1); //TODO WHY DOES THIS NOT WORK
			subrange(rB, 4, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.dtd2);
			subrange(rB, 5, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.dtd2) + m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.dtd1);

			subrange(rB, 3, 1, kr + 3, 2) = prod(trans(rActualKinematic.a1), WI1); //TODO WHY DOES THIS NOT WORK
			subrange(rB, 4, 1, kr + 3, 2) = prod(trans(rActualKinematic.a2), WI2); //TODO WHY DOES THIS NOT WORK
			subrange(rB, 5, 1, kr + 3, 2) = prod(trans(rActualKinematic.a2), WI1) + prod(trans(rActualKinematic.a1), WI2); //TODO WHY DOES THIS NOT WORK

			//shear
			subrange(rB, 6, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](0, r) * trans(rActualKinematic.t); //TODO WHY DOES THIS NOT WORK
			subrange(rB, 7, 1, kr, 3) = m_cart_deriv[IntegrationPointIndex](1, r) * trans(rActualKinematic.t);

			subrange(rB, 6, 1, kr + 3, 2) = prod(trans(rActualKinematic.a1), rVariations.P) * m_N(IntegrationPointIndex, r);
			subrange(rB, 7, 1, kr + 3, 2) = prod(trans(rActualKinematic.a2), rVariations.P) * m_N(IntegrationPointIndex, r);
		}
	}


	void Shell5pElement::CalculateGeometricStiffness(
		IndexType IntegrationPointIndex,
		Matrix& Kg,
		const KinematicVariables& rActKin,
		const VariationVariables& ractVar,
		const ConstitutiveVariables& rConstitutive)
	{
		const auto& r_geometry = GetGeometry();
		int ii1, ii2, ii3, ii4, ii5, jj1, jj2, jj3, jj4, iimult3, jjmult3;

		const SizeType number_of_control_points = GetGeometry().size();

		const Vector& StressResultants = rConstitutive.StressVector;

		const Matrix chifac = (ractVar.Chi11 * StressResultants[3] + ractVar.Chi22 * StressResultants[4] + (ractVar.Chi12 + ractVar.Chi21) * StressResultants[5]); //S[3..5] are moments
		Matrix kg2directorshear = ractVar.S1 * StressResultants[6] + ractVar.S2 * StressResultants[7]; //S[6..7] is  shear
		Matrix WI1(3, 3);
		Matrix WI2(3, 3);
		Matrix WJ1(3, 3);
		Matrix WJ2(3, 3);

		Matrix BLAI_T(2, 3); //Nodal tangent space representations
		Matrix BLAJ(3, 2);

		Matrix Temp(3, 3);


		double  Nii, dN1ii, dN2ii, Njj, dN1jj, dN2jj, dN1ii_dN2jj_p_dN2ii_dN1jj, NS, NdN1, NdN2; //the first 6 are only for cleaner coding/debugging, can be removed and m_cart_deriv[IntegrationPointIndex](0, ii) directly called

		for (int ii = 0; ii < number_of_control_points; ii++)
		{
			ii1 = 5 * ii;
			ii2 = ii1 + 1;
			ii3 = ii1 + 2;
			ii4 = ii1 + 3;
			ii5 = ii1 + 4;

			Nii = m_N(IntegrationPointIndex, ii);
			dN1ii = m_cart_deriv[IntegrationPointIndex](0, ii);
			dN2ii = m_cart_deriv[IntegrationPointIndex](1, ii);


			noalias(WI1) = ractVar.Q1 * Nii + ractVar.P * dN1ii;
			noalias(WI2) = ractVar.Q2 * Nii + ractVar.P * dN2ii;

			iimult3 = 3 * ii;

			BLAI_T = BLAA.block<3, 2>(iimult3, 0).transpose();
			for (int jj = ii; jj < number_of_control_points; jj++)
			{
				jj1 = 5 * jj;
				jj2 = jj1 + 1;
				jj3 = jj1 + 2;
				jj4 = jj1 + 3;

				Njj = m_N(IntegrationPointIndex, jj);
				dN1jj = m_cart_deriv[IntegrationPointIndex](0, jj);
				dN2jj = m_cart_deriv[IntegrationPointIndex](1, jj);

				jjmult3 = 3 * jj;

				BLAJ.noalias() = BLAA.block<3, 2>(jjmult3, 0);

				noalias(WJ1) = ractVar.Q1 * Njj + ractVar.P * dN1jj;
				noalias(WJ2) = ractVar.Q2 * Njj + ractVar.P * dN2jj;

				dN1ii_dN2jj_p_dN2ii_dN1jj = dN1ii * dN2jj + dN2ii * dN1jj;
				// MEMBRANE CONTRIBUTION
				NS = dN1ii * dN1jj * S[0] + dN2ii * dN2jj * S[1] + dN1ii_dN2jj_p_dN2ii_dN1jj * S[2];
				Kg(ii1, jj1) = Kg(ii2, jj2) = Kg(ii3, jj3) = NS;

				// BENDING CONTRIBUTION ROTATIONAL IJ
				//bending part
				Temp = StressResultants[3] * dN1ii * WJ1 + StressResultants[4] * dN2ii * WJ2 +
					StressResultants[5] * (dN1ii * WJ2 + dN2ii * WJ1);
				//TRANSVERSAL SHEAR ROTATIONAL IJ
				Temp += ractVar.P * Njj * (dN1ii * StressResultants[6] + dN2ii * StressResultants[7]);

				//Kg.block<3, 2>(ii1, jj4).noalias() = prod(Temp , BLAJ);
				noalias(subrange(Kg, ii, 3, jj4, 2)) = prod(Temp, BLAJ);

				// BENDING CONTRIBUTION ROTATIONAL JI
				Temp = S[3] * dN1jj * WI1 + S[4] * dN2jj * WI2 +
					S[5] * (dN1jj * WI2 + dN2jj * WI1);
				//TRANSVERSAL SHEAR ROTATIONAL JI
				Temp += ractVar.P * Nii * (dN1jj * StressResultants[6] + dN2jj * StressResultants[7]);

				//Kg.block<2, 3>(ii4, jj1).noalias() = prod(BLAI_T , Temp);
				noalias(subrange(Kg, ii4, 2, jj1, 3)) = prod(BLAI_T, Temp);

				// LINEARIZATION OF VARIATION OF DIRECTOR
				// linearization of variation of director shear part
				Temp = Nii * Njj * kg2directorshear;
				// linearization of variation of director bending part
				NdN1 = dN1jj * Nii + Njj * dN1ii;
				NdN2 = dN2jj * Nii + Njj * dN2ii;
				Temp += Nii * Njj * chifac;
				Temp += ractVar.S1 * NdN1 * StressResultants[3] + ractVar.S2 * NdN2 * StressResultants[4] + (ractVar.S1 * NdN2 + ractVar.S2 * NdN1) * StressResultants[5];

				//Kg.block<2, 2>(ii4, jj4).noalias() = prod(prod(BLAI_T , Temp) , BLAJ);
				noalias(subrange(Kg, ii4, 2, jj4, 2)) = prod(prod(BLAI_T, Temp), BLAJ);

			}
			///  Residual Part second part of Linearisation
			///   differencce to numerical Hessian
			/// LINEARIZATION OF PROJECTOR
			kgT = -(node[ii]->get_director_x() * fint_euk[6 * ii + 3]
				+ node[ii]->get_director_y() * fint_euk[6 * ii + 4]
				+ node[ii]->get_director_z() * fint_euk[6 * ii + 5]);
			IKg(ii4, ii4) += kgT;
			IKg(ii5, ii5) += kgT;
		}
	}

	///@}
	///@name Stiffness matrix assembly
	///@{

	inline void Shell5pElement::CalculateAndAddKm(
		MatrixType& rLeftHandSideMatrix,
		const Matrix& rB,
		const Matrix& rD,
		const double IntegrationWeight
	)
	{
		noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(rB), Matrix(prod(rD, rB)));
	}

	inline void Shell5pElement::CalculateAndAddNonlinearKm(
		Matrix& rLeftHandSideMatrix,
		const SecondVariations& rSecondVariationsStrain,
		const Vector& rSD,
		const double IntegrationWeight)
	{
		const int number_of_control_points = GetGeometry().size();
		const int mat_size = number_of_control_points * 3;

		for (int n = 0; n < mat_size; n++)
		{
			for (int m = 0; m <= n; m++)
			{
				double nm = (rSD[0] * rSecondVariationsStrain.B11(n, m)
					+ rSD[1] * rSecondVariationsStrain.B22(n, m)
					+ rSD[2] * rSecondVariationsStrain.B12(n, m)) * IntegrationWeight;

				rLeftHandSideMatrix(n, m) += nm;
				if (n != m)
					rLeftHandSideMatrix(m, n) += nm;
			}
		}
	}

	///@}
	///@name Dynamic Functions
	///@{

	void Shell5pElement::GetValuesVector(
		Vector& rValues,
		int Step)
	{
		const unsigned int number_of_control_points = GetGeometry().size();
		const unsigned int mat_size = number_of_control_points * 5;

		if (rValues.size() != mat_size)
			rValues.resize(mat_size, false);

		for (unsigned int i = 0; i < number_of_control_points; ++i)
		{
			const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
			const int index = i * 3; //TODO Geometry are 6 parameters but degrees of freedeom are 5 how to deal with this?

			rValues[index] = displacement[0];
			rValues[index + 1] = displacement[1];
			rValues[index + 2] = displacement[2];
		}
	}

	void Shell5pElement::GetFirstDerivativesVector(
		Vector& rValues,
		int Step)
	{
		const unsigned int number_of_control_points = GetGeometry().size();
		const unsigned int mat_size = number_of_control_points * 3;

		if (rValues.size() != mat_size)
			rValues.resize(mat_size, false);

		for (unsigned int i = 0; i < number_of_control_points; ++i) {
			const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
			const unsigned int index = i * 3;

			rValues[index] = velocity[0];
			rValues[index + 1] = velocity[1];
			rValues[index + 2] = velocity[2];
		}
	}

	void Shell5pElement::GetSecondDerivativesVector(
		Vector& rValues,
		int Step)
	{
		const unsigned int number_of_control_points = GetGeometry().size();
		const unsigned int mat_size = number_of_control_points * 3;

		if (rValues.size() != mat_size)
			rValues.resize(mat_size, false);

		for (unsigned int i = 0; i < number_of_control_points; ++i) {
			const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
			const unsigned int index = i * 3;

			rValues[index] = acceleration[0];
			rValues[index + 1] = acceleration[1];
			rValues[index + 2] = acceleration[2];
		}
	}

	void Shell5pElement::EquationIdVector(
		EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo
	)
	{
		KRATOS_TRY;

		const SizeType number_of_control_points = GetGeometry().size();

		if (rResult.size() != 3 * number_of_control_points)
			rResult.resize(3 * number_of_control_points, false);

		const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

		for (IndexType i = 0; i < number_of_control_points; ++i) {
			const SizeType index = i * 3;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
		}

		KRATOS_CATCH("")
	};

	void Shell5pElement::GetDofList(
		DofsVectorType& rElementalDofList,
		ProcessInfo& rCurrentProcessInfo
	)
	{
		KRATOS_TRY;

		const SizeType number_of_control_points = GetGeometry().size();

		rElementalDofList.resize(0);
		rElementalDofList.reserve(3 * number_of_control_points);

		for (IndexType i = 0; i < number_of_control_points; ++i) {
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}

		KRATOS_CATCH("")
	};

	///@}
	///@name Check
	///@{

	int Shell5pElement::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

			// Verify that the constitutive law exists
			if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
			{
				KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
			}
			else
			{
				// Verify that the constitutive law has the correct dimension
				KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
					<< "THICKNESS not provided for element " << this->Id() << std::endl;

				// Check strain size
				KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
					<< "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
					<< this->Id() << std::endl;
			}

		return 0;
	}

	///@}

} // Namespace Kratos


