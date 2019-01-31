// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Michael Loibl
//					 Tobias Teschemacher
//					 Riccardo Rossi
//
//  Based on work of Tesser and Talledo (University of Padua)		 

//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2

// System includes
#include <fstream>						// check later if really needed (ML) !!!!!!
#include <iostream>

// External includes

// Project includes
#include "tc_plastic_damage_3d_law.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
	double& TCPlasticDamage3DLaw::GetValue(
		const Variable<double>& rThisVariable,
		double& rValue)
	{
		rValue = 0.0;
		//if (rThisVariable == DAMAGE_T)
		//    rValue = m_damage_t;
		//else if (rThisVariable == DAMAGE_C)
		//    rValue = m_damage_c;
		return rValue;
	}

	Vector& TCPlasticDamage3DLaw::GetValue(
		const Variable<Vector>& rThisVariable,
		Vector& rValue)
	{
		//if (rThisVariable == EIGENVALUE_VECTOR)
		//    rValue = m_eigen_values;
		return rValue;
	}

	Matrix& TCPlasticDamage3DLaw::GetValue(
		const Variable<Matrix>& rThisVariable,
		Matrix& rValue)
	{
		//if (rThisVariable == EIGENVECTOR_MATRIX)
		//    rValue = m_eigen_vectors;
		return rValue;
	}

	void TCPlasticDamage3DLaw::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		m_compressive_strength = rMaterialProperties[UNIAXIAL_STRESS_COMPRESSION];
		m_tensile_strength = rMaterialProperties[UNIAXIAL_STRESS_TENSION];

		m_elastic_strain = ZeroVector(6);
		m_plastic_strain = ZeroVector(6);

		// CalculateElasticityMatrix(m_D0);

		/**

		m_rate_biaxial_uniaxial = rMaterialProperties[RATE_BIAXIAL_UNIAXIAL];

		m_beta = rMaterialProperties[BETA];

		m_compression_parameter_A = rMaterialProperties[COMPRESSION_PARAMETER_A];
		m_compression_parameter_B = rMaterialProperties[COMPRESSION_PARAMETER_B];


		m_E = rMaterialProperties[YOUNG_MODULUS];
		m_nu = rMaterialProperties[POISSON_RATIO];

		m_Gf_t = rMaterialProperties[FRACTURE_ENERGY_TENSION];
		m_Gf_c = rMaterialProperties[FRACTURE_ENERGY_COMPRESSION];

		double l_c = 0.5;

		m_tension_parameter_A = 1 / ((1 - m_beta)*((m_Gf_t*m_E / (l_c*m_tensile_strength*m_tensile_strength)) - 0.5));//rMaterialProperties[TENSION_PARAMETER_A];

		KRATOS_WATCH(m_tension_parameter_A)

		m_Di = m_D0;

		m_K = std::sqrt(2) * ((1- m_rate_biaxial_uniaxial) / (1 - 2 * m_rate_biaxial_uniaxial));

		model = 2;

		Vector d_tension = ZeroVector(2);
		d_tension(0) = m_tensile_strength;

		CalculateTresholdTension(d_tension, m_treshold_tension_initial);

		if (model == 1)
		{
			m_treshold_compression_initial = std::sqrt(std::sqrt(3.0)*(m_K - sqrt(2.0))*m_compressive_strength / 3);
			//m_treshold_tension_initial = sqrt(m_tensile_strength / sqrt(m_E));
		}
		if (model == 2)
		{
			m_treshold_compression_initial = std::sqrt(std::sqrt(3.0)*(m_K - sqrt(2.0))*m_compressive_strength / 3);
			//m_treshold_tension_initial = std::sqrt(m_tensile_strength);
		}
		if (model == 3)
		{
			m_treshold_compression_initial = std::sqrt(3.0)*(m_K - std::sqrt(2.0))*m_compressive_strength / 3;
			//m_treshold_tension_initial = m_tensile_strength;
		}

		m_treshold_tension = m_treshold_tension_initial;
		m_treshold_compression = m_treshold_compression_initial;

		m_damage_t = 0.0;
		m_damage_c = 0.0;

		m_gamma_C = 0.0;
		*/
	}

	void TCPlasticDamage3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
	{
		const ProcessInfo&  pinfo = rValues.GetProcessInfo();
		const GeometryType& geom = rValues.GetElementGeometry();
		const Properties&   props = rValues.GetMaterialProperties();

		const Vector& rStrainVector = rValues.GetStrainVector();
		Vector&       rStressVector = rValues.GetStressVector();
		Matrix& rConstitutiveLaw = rValues.GetConstitutiveMatrix();

		this->CalculateMaterialResponseInternal(rStrainVector, rStressVector, rConstitutiveLaw);
	}

	void TCPlasticDamage3DLaw::CalculateMaterialResponseInternal(
		const Vector& rStrainVector,
		Vector& rStressVector,
		Matrix& rConstitutiveLaw)
	{
		if (rStressVector.size() != 6)
		{
			rStressVector.resize(6, false);
		}
		const double tolerance = (1.0e-14)* m_compressive_strength;

		// 1.step: elastic stress tensor
		noalias(rStressVector) = prod(trans(m_D0), (rStrainVector - m_plastic_strain));

		// 2.step: spectral decomposition
		Vector StressVectorTension = ZeroVector(6);
		Vector StressVectorCompression = ZeroVector(6);
		Matrix PMatrixTension = ZeroMatrix(6, 6);
		Matrix PMatrixCompression = ZeroMatrix(6, 6);

		SpectralDecomposition(rStressVector,
			StressVectorTension,
			StressVectorCompression,
			PMatrixTension,
			PMatrixCompression);


		/**

		if (m_beta > 0.0)
		{
			Matrix V = ZeroMatrix(2, 2);
			Vector d = ZeroVector(2);
			solve_eig(rStressVector(0), rStressVector(2), rStressVector(2), rStressVector(1), V, d);

			Vector d_compression = ZeroVector(2);
			d_compression(0) = (d(0) - std::abs(d(0))) / 2;
			d_compression(1) = (d(1) - std::abs(d(1))) / 2;

			Vector d_tension = ZeroVector(2);
			d_tension(0) = (d(0) + std::abs(d(0))) / 2;
			d_tension(1) = (d(1) + std::abs(d(1))) / 2;

			// compute the stress measures
			double treshold_tension = 0.0;
			double treshold_compression = 0.0;

			CalculateTresholdTension(
				d_tension,
				treshold_tension);
			CalculateTresholdCompression(
				d_compression,
				treshold_compression);

			CalculatePlasticStrain(treshold_tension, treshold_compression, rStressVector, m_plastic_strain);
		}


		Matrix V = ZeroMatrix(2, 2);
		Vector d = ZeroVector(2);
		solve_eig(rStressVector(0), rStressVector(2), rStressVector(2), rStressVector(1), V, d);
		m_eigen_values = d;
		m_eigen_vectors = V;
		//KRATOS_WATCH(rStressVector)
		//KRATOS_WATCH(d)
		Vector d_compression = ZeroVector(2);
		d_compression(0) = (d(0) - std::abs(d(0))) / 2;
		d_compression(1) = (d(1) - std::abs(d(1))) / 2;

		if (d_compression(0)>1e-1)
			KRATOS_WATCH(d_compression(0))
		if (d_compression(1)>1e-1)
			KRATOS_WATCH(d_compression(1))


		Vector d_tension = ZeroVector(2);
		d_tension(0) = (d(0) + std::abs(d(0))) / 2;
		d_tension(1) = (d(1) + std::abs(d(1))) / 2;

		// compute the stress measures
		double treshold_tension = 0.0;
		double treshold_compression = 0.0;

		CalculateTresholdTension(
			d_tension,
			treshold_tension);
		CalculateTresholdCompression(
			d_compression,
			treshold_compression);

		// compute the equivalent stress measures
		double damage_treshold_tension = 0.0;
		double damage_treshold_compression = 0.0;

		CalculateUniqueDamageCriterion(treshold_tension, treshold_compression,
			damage_treshold_tension, damage_treshold_compression);

		//std::cout << "we are here: UniqueDamageCriterion: damage_tr_tension: " << damage_treshold_tension << ", damage_treshold_compression: " << damage_treshold_compression << std::endl;
		double damage_t = 0.0;
		double damage_c = 0.0;

		// damage update
		CalculateDamageTension(
			damage_treshold_tension,
			damage_t);

		CalculateDamageCompression(
			damage_treshold_compression,
			damage_c);
		//std::cout << "we are here: CalculateDamageCompression" << std::endl;

		m_damage_c = std::max(damage_c, m_damage_c);
		m_damage_t = std::max(damage_t, m_damage_t);

		double shear_retention_factor = 0.0;
		//if (m_gamma_C > 0.0)
		//{
		//    shear_retention_factor = 1 - std::abs(rStrainVector(2)) / (2 * m_gamma_C);
		//    shear_retention_factor = std::max(shear_retention_factor, 0.0);
		//}

		Matrix supp1 = ZeroMatrix(2, 2);
		Matrix supp2 = ZeroMatrix(2, 2);

		//KRATOS_WATCH(damage_treshold_compression)
		//KRATOS_WATCH(m_treshold_compression_initial)

		//KRATOS_WATCH(damage_treshold_tension)
		//KRATOS_WATCH(m_treshold_tension_initial)

		//KRATOS_WATCH(d_compression)
		//KRATOS_WATCH(d_tension)
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				//supp1[i][j] = 0.0;
				//supp2[i][j] = 0.0;
				supp1(i, j) = supp1(i, j) + V(i, j) * d_compression[j];
				supp2(i, j) = supp2(i, j) + V(i, j) * d_tension[j];
			}
		}
		Matrix D_compression = ZeroMatrix(2, 2);
		Matrix D_tension = ZeroMatrix(2, 2);
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					D_compression(i, j) = D_compression(i, j) + supp1(i, k) * V(j, k);
					D_tension(i, j) = D_tension(i, j) + supp2(i, k) * V(j, k);
				}
			}
		}

		Matrix Stress2d = ZeroMatrix(2, 2);
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				//KRATOS_WATCH(D_compression)
				//KRATOS_WATCH(D_tension)

				//Stress2d(i,j)=(1-dn)*Dn[i][j]+(1-dp)*Dp[i][j];
				// 11/03/2013 Diego Talledo: Added Environmental Chemical Damage
				Stress2d(i, j) = (1 - m_damage_c)*D_compression(i,j) + (1 - m_damage_t)*D_tension(i, j);
				// 13/01/2013 Diego Talledo: Added Shear Retention Factor
				if (((i == 0) && (j == 1)) || ((i == 1) && (j == 0)))
				{
					//if (!srfCompr) // 11/03/2013 Diego Talledo: Apply SRF also to compression.
					//    Stress2d(i, j) = (1 - dnstar)*D_compression[i][j] + (1 - (1 - SRF12)*dpstar)*Dp[i][j];
					//else
					Stress2d(i, j) = (1 - (1 - shear_retention_factor)*m_damage_c)*D_compression(i, j)
						+ (1 - (1 - shear_retention_factor)*m_damage_t)*D_tension(i, j);
				}
			}
		}

		//KRATOS_WATCH(m_damage_c)
		//KRATOS_WATCH(m_damage_t)

		// calculation of stress tensor
		noalias(rStressVector) = (1.0 - m_damage_t)*stress_vector_tension
			+ (1.0 - m_damage_c)*stress_vector_compression;
		rStressVector(2) = (1.0 - (1 - shear_retention_factor)*m_damage_t)*stress_vector_tension(2)
			+ (1.0 - (1 - shear_retention_factor)*m_damage_c)*stress_vector_compression(2);

		//KRATOS_WATCH(rStressVector)
		//KRATOS_WATCH(Stress2d)

		//KRATOS_WATCH(D_compression)
		//KRATOS_WATCH(D_tension)

		//KRATOS_WATCH(d)
		//KRATOS_WATCH(V)
		//KRATOS_WATCH(p_matrix_compression)
		//KRATOS_WATCH(p_matrix_tension)

		Matrix Damage = (1-m_damage_t) * p_matrix_tension + (1-m_damage_c) * p_matrix_compression;
		Damage(0, 1) = (1 - (1 - shear_retention_factor) * m_damage_t) * p_matrix_tension(0, 1)
			+ (1 - (1 - shear_retention_factor) * m_damage_c) * p_matrix_compression(0, 1);
		Damage(1, 0) = (1 - (1 - shear_retention_factor) * m_damage_t) * p_matrix_tension(1, 0)
			+ (1 - (1 - shear_retention_factor) * m_damage_c) * p_matrix_compression(1, 0);

		//KRATOS_WATCH(Damage)
		//KRATOS_WATCH(m_D0)

		noalias(rConstitutiveLaw) = prod(Damage, m_D0);
		//// Second Secant Operator
		//Matrix D_S2 = 0.5*(rConstitutiveLaw + trans(rConstitutiveLaw));
		////noalias(rConstitutiveLaw) = D_S2;

		////KRATOS_WATCH(rStrainVector)
		////KRATOS_WATCH(m_D0)
		////KRATOS_WATCH(rConstitutiveLaw)

		//Vector lalala = prod(trans(rConstitutiveLaw), rStrainVector);
		//Vector error = lalala - rStressVector;
		////KRATOS_WATCH(lalala)
		////KRATOS_WATCH(error)

		//Vector lalala2 = prod(trans(D_S2), rStrainVector);
		//Vector error2 = lalala2 - rStressVector;
		////KRATOS_WATCH(lalala2)
		////KRATOS_WATCH(error2)

		//Matrix Q_CW_Tension = ZeroMatrix(3, 3);
		//SpectralDecompositionStrain(rStrainVector, Q_CW_Tension);

		//Matrix A = std::sqrt(1 - damage_t)*Q_CW_Tension + std::sqrt(1 - damage_c)*(IdentityMatrix(3,3) - Q_CW_Tension);
		//Matrix A_C0 = prod(A, m_D0);
		//Matrix D_E = prod(A_C0, A);
		////noalias(rConstitutiveLaw) = D_E;
		////m_Di = rConstitutiveLaw;
		//Vector lalala3 = prod(trans(D_E), rStrainVector);
		//Vector error3 = lalala3 - rStressVector;
		//KRATOS_WATCH(D_E)
		//KRATOS_WATCH(lalala3)
		//KRATOS_WATCH(error3)
		//noalias(rConstitutiveLaw) = m_D0;
		m_treshold_tension = damage_treshold_tension;
		m_treshold_compression = damage_treshold_compression;
		*/
	}

	void TCPlasticDamage3DLaw::SpectralDecomposition(
		const Vector& rStressVector,
		Vector& rStressVectorTension,
		Vector& rStressVectorCompression,
		Matrix& PMatrixTension,
		Matrix& PMatrixCompression)
	{
		Matrix helpmatrix = ZeroMatrix(6, 6);

		rStressVectorTension = ZeroVector(6);
		rStressVectorCompression = ZeroVector(6);
		PMatrixTension = ZeroMatrix(6, 6);
		PMatrixCompression = ZeroMatrix(6, 6);

		BoundedMatrix<double, 3, 3> Stress33;
		BoundedMatrix<double, 3, 3> EigenvalueStress33;
		BoundedMatrix<double, 3, 3> EigenvectorStress33;

		Stress33 = MathUtils<double>::StressVectorToTensor(rStressVector);

		bool converged = MathUtils<double>::EigenSystem<3>(Stress33, EigenvectorStress33, EigenvalueStress33);

		std::vector<Vector> p_i(3);
		for (IndexType i = 0; i < 3; ++i) {
			p_i[i] = ZeroVector(3);
			for (IndexType j = 0; j < 3; ++j) {
				p_i[i](j) = EigenvectorStress33(j, i);
			}
		}

		for (IndexType i = 0; i < 3; ++i) {
			if (EigenvalueStress33(i, i) > 0.0) {
				helpmatrix  = outer_prod(p_i[i], p_i[i]); // p_i x p_i
				rStressVectorTension += MathUtils<double>::StressTensorToVector(EigenvalueStress33(i, i) * helpmatrix);
			}
		}
		rStressVectorCompression = rStressVector - rStressVectorTension;

		PMatrixTension = ZeroMatrix(6, 6);
		for (SizeType i = 0.0; i < 3; ++i)
		{
			if (EigenvalueStress33(i, i) > 0.0)
			{
				Vector helpvector = MathUtils<double>::StressTensorToVector(helpmatrix);
				PMatrixTension += outer_prod(helpvector, helpvector);
			}
		}
		Matrix I = IdentityMatrix(6, 6);
		PMatrixCompression = I - PMatrixTension;
	};
}
