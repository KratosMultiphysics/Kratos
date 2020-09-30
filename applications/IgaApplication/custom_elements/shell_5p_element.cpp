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
		if (reference_Curvature.size() != r_number_of_integration_points)
			reference_Curvature.resize(r_number_of_integration_points);
		if (reference_TransShear.size() != r_number_of_integration_points)
			reference_TransShear.resize(r_number_of_integration_points);
		if (m_dA_vector.size() != r_number_of_integration_points)
			m_dA_vector.resize(r_number_of_integration_points);
		if (m_cart_deriv.size() != r_number_of_integration_points)
			m_cart_deriv.resize(r_number_of_integration_points);

		KinematicVariables kinematic_variables;
		VariationVariables variation_variables;

		m_N = GetGeometry().ShapeFunctionsValues(); //get shape functions

		for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
		{
			m_cart_deriv[point_number]= CalculateCartesianDerivatives(point_number, kinematic_variables );

			CalculateKinematics(
				point_number,
				kinematic_variables, variation_variables);

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
			KinematicVariables kinematic_variables;
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
				constitutive_variables,
				constitutive_law_parameters,
				ConstitutiveLaw::StressMeasure_PK2);

			// calculate B operator
			Matrix BOperator = CalculateStrainDisplacementOperator(point_number,
				kinematic_variables,
				variation_variables);

			double integration_weight =
				r_integration_points[point_number].Weight()
				* m_dA_vector[point_number]
				* GetProperties()[THICKNESS]; // divide by onehalf?

			// LEFT HAND SIDE MATRIX
			if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
			{
				// Geometric stiffness matrix
				Matrix Kg =
				CalculateGeometricStiffness(
					point_number,
					kinematic_variables,
					variation_variables,
					constitutive_variables);
				//Matrix temp = prod(trans(BOperator), constitutive_variables.ConstitutiveMatrix);  //defined useless temporary due to ublas prod(prod( not working
				rLeftHandSideMatrix = rLeftHandSideMatrix + integration_weight * (prod(prod<MatrixType>(trans(BOperator), constitutive_variables.ConstitutiveMatrix), BOperator)+ Kg);
			}
			// RIGHT HAND SIDE VECTOR
			if (CalculateResidualVectorFlag == true) 
			{
				// operation performed: rRightHandSideVector -= Weight*IntForce
				noalias(rRightHandSideVector) -= integration_weight * prod(trans(BOperator), constitutive_variables.StressVector);
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
		rKin.metric[0] = norm_2_square(rKin.a1); //TODO calculate metric from displacements and not from geometry. This is more accurate for small strains due cancelation with reference metric since 1.0-1.0 ~ loses digits 
		rKin.metric[1] = norm_2_square(rKin.a2);
		rKin.metric[2] = inner_prod(rKin.a1, rKin.a2);

		rKin.t = ZeroVector(3);
		rKin.dtd1 = ZeroVector(3);
		rKin.dtd2 = ZeroVector(3);


		const SizeType number_of_nodes = GetGeometry().size();
		for (size_t i = 0; i < number_of_nodes; i++)
		{
			rKin.t += m_N(IntegrationPointIndex, i) * GetGeometry()[i].GetValue(DIRECTOR); //TODO does this really do what i want? t=sum_I N_I *t_I ? and is this the best way to obtain this?
			rKin.dtd1 += m_cart_deriv[IntegrationPointIndex](0, i) * GetGeometry()[i].GetValue(DIRECTOR);
			rKin.dtd2 += m_cart_deriv[IntegrationPointIndex](1, i) * GetGeometry()[i].GetValue(DIRECTOR);
		}

		double invL_t = 1.0 / norm_2(rKin.t);
		rKin.t *= invL_t;

		const array_1d<double, 3>& t = rKin.t; //define alias for cleaner code
		const array_1d<double, 3>& dtd1 = rKin.dtd1; //define alias for cleaner code
		const array_1d<double, 3>& dtd2 = rKin.dtd2; //define alias for cleaner code
		const array_1d<double, 3>& a1 = rKin.a1; //define alias for cleaner code
		const array_1d<double, 3>& a2 = rKin.a2; //define alias for cleaner code

		rVar.P = outer_prod(t, t) - IdentityMatrix(3);
		rVar.P *= invL_t;

		invL_t *= invL_t;
		const Matrix tdyadt = outer_prod(t, t);

		rKin.transShear[0] = inner_prod(t, a1);
		rKin.transShear[1] = inner_prod(t, a2);


		rVar.Q1 = invL_t * (inner_prod(t, dtd1) * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(dtd1, t) - outer_prod(t, dtd1));
		rVar.Q2 = invL_t * (inner_prod(t, dtd2) * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(dtd2, t) - outer_prod(t, dtd2));
		rVar.S1 = invL_t * (rKin.transShear[0] * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(a1, t) - outer_prod(t, a1));
		rVar.S2 = invL_t * (rKin.transShear[1] * (3.0 * tdyadt - IdentityMatrix(3)) - outer_prod(a2, t) - outer_prod(t, a2));

		invL_t *= invL_t;
		rVar.Chi11 = invL_t * (3.0 * inner_prod(t, dtd1) * (outer_prod(a1, t) + 0.5 * rKin.transShear[0] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a1, dtd1) * tdyadt + rKin.transShear[0] * outer_prod(dtd1, t)) - outer_prod(a1, dtd1) - inner_prod(a1, dtd1) * 0.5 * IdentityMatrix(3));
		rVar.Chi21 = invL_t * (3.0 * inner_prod(t, dtd1) * (outer_prod(a2, t) + 0.5 * rKin.transShear[1] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a2, dtd1) * tdyadt + rKin.transShear[1] * outer_prod(dtd1, t)) - outer_prod(a2, dtd1) - inner_prod(a2, dtd1) * 0.5 * IdentityMatrix(3));
		rVar.Chi12 = invL_t * (3.0 * inner_prod(t, dtd2) * (outer_prod(a1, t) + 0.5 * rKin.transShear[0] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a1, dtd2) * tdyadt + rKin.transShear[0] * outer_prod(dtd2, t)) - outer_prod(a1, dtd2) - inner_prod(a1, dtd2) * 0.5 * IdentityMatrix(3));
		rVar.Chi22 = invL_t * (3.0 * inner_prod(t, dtd2) * (outer_prod(a2, t) + 0.5 * rKin.transShear[1] * (IdentityMatrix(3) - 5.0 * tdyadt)) + 3.0 * (0.5 * inner_prod(a2, dtd2) * tdyadt + rKin.transShear[1] * outer_prod(dtd2, t)) - outer_prod(a2, dtd2) - inner_prod(a2, dtd2) * 0.5 * IdentityMatrix(3));

		//up to there dtd_al corresponds to w_al
		rKin.dtd1 = prod(rVar.P, dtd1);
		rKin.dtd2 = prod(rVar.P, dtd2);

		rKin.curvature[0] = inner_prod(rKin.a1, dtd1);
		rKin.curvature[1] = inner_prod(rKin.a2, dtd2);
		rKin.curvature[2] = inner_prod(rKin.a1, dtd2) + inner_prod(a2, dtd1);

	}

	/* Transforms derivatives to obtain cartesian quantities  */
	Matrix Shell5pElement::CalculateCartesianDerivatives(
		IndexType IntegrationPointIndex,
		const KinematicVariables& rKinematicVariables
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

		return prod(invJ, r_DN_De);
	}

	void Shell5pElement::CalculateConstitutiveVariables(
		IndexType IntegrationPointIndex,
		KinematicVariables& rActualKinematic,
		ConstitutiveVariables& rThisConstitutiveVariables,
		ConstitutiveLaw::Parameters& rValues,
		const ConstitutiveLaw::StressMeasure ThisStressMeasure
	)
	{
		rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
		rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
		rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		array_1d<double, 8> strain_vector;
		subrange(strain_vector,0,2) = rActualKinematic.metric;
		strain_vector(0) -= 1.0; //TODO this is unnecessary if TODO from calculatekinematics is done
		strain_vector(1) -= 1.0;

		strain_vector(0) *= 0.5;
		strain_vector(1) *= 0.5;

		subrange(strain_vector,3,5) = rActualKinematic.curvature - reference_Curvature[IntegrationPointIndex];

		subrange(strain_vector, 6, 7) = rActualKinematic.transShear - reference_TransShear[IntegrationPointIndex];
		noalias(rThisConstitutiveVariables.StrainVector) = strain_vector;

		rValues.SetStrainVector(rThisConstitutiveVariables.StrainVector); //this is the input parameter
		rValues.SetStressVector(rThisConstitutiveVariables.StressVector);    //this is an ouput parameter
		rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix); //this is an ouput parameter

		mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);

		double thickness = this->GetProperties().GetValue(THICKNESS);
	    //	noalias(rThisConstitutiveVariablesCurvature.ConstitutiveMatrix) = rThisConstitutiveVariablesMembrane.ConstitutiveMatrix * (pow(thickness, 2) / 12);   //TODO this does not work for general material laws especially including shear

		//Local Cartesian Forces and Moments and Shear Forces
		noalias(rThisConstitutiveVariables.StressVector) = prod(
			rThisConstitutiveVariables.ConstitutiveMatrix, rThisConstitutiveVariables.StrainVector);
	}

	Matrix Shell5pElement::CalculateStrainDisplacementOperator(
		IndexType IntPoint,
		const KinematicVariables& rActualKinematic,
		const VariationVariables& rVariations)
	{
		const SizeType number_of_nodes = GetGeometry().size();
		const SizeType num_dof = number_of_nodes * 5;

		Matrix rB{ ZeroMatrix(8, num_dof) };

		for (SizeType r = 0; r < number_of_nodes; r++)
		{
			SizeType kr = 5 * r;

			const Matrix32d BLAI = GetGeometry()[r].GetValue(DIRECTORTANGENTSPACE);
			const Matrix3d WI1 = rVariations.Q1 * m_N(IntPoint, r) + rVariations.P * m_cart_deriv[IntPoint](0, r);
			const Matrix3d WI2 = rVariations.Q2 * m_N(IntPoint, r) + rVariations.P * m_cart_deriv[IntPoint](1, r);

			Vector Temp(1, 2), Temp1(1, 2), Temp2(1, 2), Temp3(1, 2), Temp4(1, 2);

			for (int s = 0; s < 3; s++)
			{
				rB(0, kr + s) = m_cart_deriv[IntPoint](0, r) * rActualKinematic.a1[s]; //membrane_{,disp}
				rB(1, kr + s) = m_cart_deriv[IntPoint](1, r) * rActualKinematic.a2[s];
				rB(2, kr + s) = m_cart_deriv[IntPoint](0, r) * rActualKinematic.a2[s] + m_cart_deriv[IntPoint](1, r) * rActualKinematic.a1[s];

				rB(3, kr + s) = m_cart_deriv[IntPoint](0, r) * rActualKinematic.dtd1[s]; //bending_{,disp}
				rB(4, kr + s) = m_cart_deriv[IntPoint](1, r) * rActualKinematic.dtd2[s];
				rB(5, kr + s) = m_cart_deriv[IntPoint](0, r) * rActualKinematic.dtd2[s] + m_cart_deriv[IntPoint](1, r) * rActualKinematic.dtd1[s];

				rB(6, kr + s) = m_cart_deriv[IntPoint](0, r) * rActualKinematic.t[s];  //trans_shear_{,disp}
				rB(7, kr + s) = m_cart_deriv[IntPoint](1, r) * rActualKinematic.t[s];
			}

			 Temp = prod(prod<Vector>(trans(rActualKinematic.a1), WI1), BLAI); //bending_{,dir}
			 Temp1 = prod(prod<Vector>(trans(rActualKinematic.a2), WI2), BLAI);
			 Temp2 = prod(prod<Vector>(trans(rActualKinematic.a2), WI1) + prod(trans(rActualKinematic.a1), WI2), BLAI);

			 Temp3 = prod(prod<Vector>(trans(rActualKinematic.a1), rVariations.P) * m_N(IntPoint, r), BLAI); //shear_{,dir}
			 Temp4 = prod(prod<Vector>(trans(rActualKinematic.a2), rVariations.P) * m_N(IntPoint, r), BLAI);

			for (int s = 0; s < 2; s++)
			{
				rB(3, kr + 3 + s) = Temp(s);
				rB(4, kr + 3 + s) = Temp1(s);
				rB(5, kr + 3 + s) = Temp2(s);
				rB(6, kr + 3 + s) = Temp3(s);
				rB(7, kr + 3 + s) = Temp4(s);
			}///=============== end of stupid version ===============


			////membrane
			//noalias(subrange(rB, 0, 1, kr, 3)) = m_cart_deriv[IntPoint](0, r) * trans(rActualKinematic.a1);
			//noalias(subrange(rB, 1, 1, kr, 3)) = m_cart_deriv[IntPoint](1, r) * trans(rActualKinematic.a2);

			//////inplane shear
			//noalias(subrange(rB, 2, 1, kr, 3)) = m_cart_deriv[IntPoint](0, r) * trans(rActualKinematic.a2) + m_cart_deriv[IntPoint](1, r) * trans(rActualKinematic.a1);

			//////BENDING PART displacement variation
			//noalias(subrange(rB, 3, 1, kr, 3)) = m_cart_deriv[IntPoint](0, r) * trans(rActualKinematic.dtd1);
			//noalias(subrange(rB, 4, 1, kr, 3)) = m_cart_deriv[IntPoint](1, r) * trans(rActualKinematic.dtd2);
			//noalias(subrange(rB, 5, 1, kr, 3)) = m_cart_deriv[IntPoint](0, r) * trans(rActualKinematic.dtd2); + m_cart_deriv[IntPoint](1, r) * trans(rActualKinematic.dtd1);

			//////BENDING PART director variation
			//noalias(subrange(rB, 3, 1, kr + 3, 2)) = prod(prod<Vector>(trans(rActualKinematic.a1), WI1),BLAI);
			//noalias(subrange(rB, 4, 1, kr + 3, 2)) = prod(prod<Vector>(trans(rActualKinematic.a2), WI2),BLAI);
			//noalias(subrange(rB, 5, 1, kr + 3, 2)) = prod(prod<Vector>(trans(rActualKinematic.a2), WI1) + prod(trans(rActualKinematic.a1), WI2), BLAI);

			//////TransverseShear displacement variation
			//noalias(subrange(rB, 6, 1, kr, 3)) = m_cart_deriv[IntPoint](0, r) * trans(rActualKinematic.t);
			//noalias(subrange(rB, 7, 1, kr, 3)) = m_cart_deriv[IntPoint](1, r) * trans(rActualKinematic.t);

			//////TransverseShear director variation
			//noalias(subrange(rB, 6, 1, kr + 3, 2)) = prod(prod<Vector>(trans(rActualKinematic.a1), rVariations.P) * m_N(IntPoint, r), BLAI);
			//noalias(subrange(rB, 7, 1, kr + 3, 2)) = prod(prod<Vector>(trans(rActualKinematic.a2), rVariations.P) * m_N(IntPoint, r), BLAI);
		}
		return rB;
	}


	Matrix Shell5pElement::CalculateGeometricStiffness(
		IndexType iP,  //Integration Point
		const KinematicVariables& rActKin,
		const VariationVariables& ractVar,
		const ConstitutiveVariables& rConstitutive)
	{
		const auto& r_geometry = GetGeometry();
		const SizeType number_of_control_points = r_geometry.size();
		const SizeType num_dof = number_of_control_points * 5;

		Matrix Kg{ ZeroMatrix(num_dof, num_dof) };
		const array_1d<double, 8>& S = rConstitutive.StressVector;

		const Matrix3d chifac = (ractVar.Chi11 * S[3] + ractVar.Chi22 * S[4] + (ractVar.Chi12 + ractVar.Chi21) * S[5]); //S[3..5] are moments
		const Matrix3d kg2directorshear = ractVar.S1 * S[6] + ractVar.S2 * S[7]; //S[6..7] is  shear

		for (SizeType i = 0; i < number_of_control_points; i++)
		{
			const int i1 = 5 * i;
			const int i4 = i1 + 3;

			const double Ni = m_N(iP, i);
			const double dN1i = m_cart_deriv[iP](0, i);
			const double dN2i = m_cart_deriv[iP](1, i);

			const Matrix3d WI1 = ractVar.Q1 * Ni + ractVar.P * dN1i;
			const Matrix3d WI2 = ractVar.Q2 * Ni + ractVar.P * dN2i;
			const Matrix23d BLAI_T = trans(r_geometry[i].GetValue(DIRECTORTANGENTSPACE));
			for (SizeType j = i; j < number_of_control_points; j++)
			{
				const int j1 = 5 * j;
				const int j4 = j1 + 3;

				const double Nj = m_N(iP, j);
				const double dN1j = m_cart_deriv[iP](0, j);
				const double dN2j = m_cart_deriv[iP](1, j);

				const Matrix3d WJ1 = ractVar.Q1 * Nj + ractVar.P * dN1j;
				const Matrix3d WJ2 = ractVar.Q2 * Nj + ractVar.P * dN2j;
				const Matrix32d BLAJ = r_geometry[j].GetValue(DIRECTORTANGENTSPACE);

				const double NS = dN1i * dN1j * S[0] + dN2i * dN2j * S[1] + (dN1i * dN2j + dN2i * dN1j) * S[2];
				Kg(i1, j1) = Kg(i1+1, j1+1) = Kg(i1+2, j1+2) = NS; // membrane_{,disp,disp}*N

				Matrix3d Temp = S[3] * dN1i * WJ1 + S[4] * dN2i * WJ2 + S[5] * (dN1i * WJ2 + dN2i * WJ1); // bending_{,dir,disp}*M
				Temp += ractVar.P * Nj * (dN1i * S[6] + dN2i * S[7]);  // shear_{,dir,disp}*Q

				noalias(subrange(Kg, i, 3, j4, 2)) = prod(Temp, BLAJ);

				Temp = S[3] * dN1j * WI1 + S[4] * dN2j * WI2 + S[5] * (dN1j * WI2 + dN2j * WI1); // bending_{,disp,dir}*M
	
				Temp += ractVar.P * Ni * (dN1j * S[6] + dN2j * S[7]);// shear_{,disp,dir}*Q

				noalias(subrange(Kg, i4, 2, j1, 3)) = prod(BLAI_T, Temp);

				Temp = Ni * Nj * kg2directorshear;  // shear_{,dir,dir}*Q

				const double NdN1 = dN1j * Ni + Nj * dN1i;
				const double NdN2 = dN2j * Ni + Nj * dN2i;
				Temp += Ni * Nj * chifac;  // bending_{,dir,dir}*M
				Temp += ractVar.S1 * NdN1 * S[3] + ractVar.S2 * NdN2 * S[4] + (ractVar.S1 * NdN2 + ractVar.S2 * NdN1) * S[5];  // bending_{,dir,dir}*M

				const Matrix23d Temp2 = prod(BLAI_T, Temp); //useless temp due to ublas prodprod
				noalias(subrange(Kg, i4, 2, j4, 2)) = prod(Temp2, BLAJ);
			}
			const double kgT = -inner_prod(r_geometry[i].GetValue(DIRECTOR), 
				                     prod(WI1, rActKin.a1) * S[3] +
				                     prod(WI2, rActKin.a2) * S[4] + 
				                     prod(WI2, rActKin.a1) + prod(WI1, rActKin.a2) * S[5] +
			                         prod(ractVar.P, rActKin.a1) * Ni * S[6] +
								     prod(ractVar.P, rActKin.a2) * Ni * S[7]); //P’_{,dir}*F_{int}
			Kg(i1 + 3, i1 + 3) += kgT;
			Kg(i1 + 4, i1 + 4) += kgT;
		}
		return Kg;
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

		for (int i = 0; i < number_of_control_points; ++i)
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

		if (rResult.size() != 5 * number_of_control_points)
			rResult.resize(5 * number_of_control_points, false);

		const SizeType pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

		for (IndexType i = 0; i < number_of_control_points; ++i) {
			const SizeType index = i * 5;
			rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos    ).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
			rResult[index + 3] = GetGeometry()[i].GetDof(DIRECTORINC_X).EquationId();
			rResult[index + 4] = GetGeometry()[i].GetDof(DIRECTORINC_Y).EquationId();
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

	///@}
	///@name Check
	///@{

	int Shell5pElement::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		

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



	//void Shell5pElement::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) //update before next iteration
	//{
	//	for (int i=1 ; i< GetGeometry().size(); i++)  //update for each node after every solution step
	//	{ 

	//		array_1d<double, 3> director = GetGeometry()[i].GetValue(DIRECTOR);
	//		array_1d<double, 2> inc2d = GetGeometry()[i].GetValue(DIRECTORINC);
	//		Matrix BLA(3, 2);
	//		BLA= GetGeometry()[i].GetValue(DIRECTORTANGENTSPACE);

	//		array_1d<double, 3> inc3d = prod(BLA,inc2d);

	//		//projection-based update
	//		director = director+inc3d;
	//		director = director / sqrt(inner_prod(director, director));

	//		GetGeometry()[i].SetValue(DIRECTOR, director);

	//		//TODO: update tangentspace
	//	}
	//	
	//}




	//void Shell5pElement::constructReferenceDirectorL2FitSystem(MatrixType& rLeftHandSideMatrix, MatrixType& rRightHandSideMatrix) 
	//{
	//	    //size of rLeftHandSideMatrix(num_node, num_node);
	//	    //size of rRightHandSideMatrix(num_node, 3);

	//		const auto& r_geometry = GetGeometry();
	//		const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();
	//		const SizeType r_number_of_nodes = r_geometry.size();

	//		array_1d<double, 3 > A3;

	//		for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number)
	//		{
	//			A3 = GetGeometry().Normal(point_number); //this makes only sense if the geometry is undeformed
	//			A3 = A3 / sqrt(inner_prod(A3, A3));
	//			//Vector m_Nvec = trans(row(m_N, point_number));
	//				rRightHandSideMatrix = rRightHandSideMatrix + outer_prod(trans(row(m_N, point_number)),  trans(A3)     ) ;

	//			rLeftHandSideMatrix = rLeftHandSideMatrix + outer_prod(row(m_N,point_number), row(m_N, point_number));
	//		}
	//}
	///@}

} // Namespace Kratos


