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
		m_f_01cc = rMaterialProperties[UNIAXIAL_STRESS_COMPRESSION];
		m_f_ct = rMaterialProperties[UNIAXIAL_STRESS_TENSION];
		m_f_02cc = rMaterialProperties[BIAXIAL_STRESS_COMPRESSION];
		m_E = rMaterialProperties[YOUNG_MODULUS];
		m_nu = rMaterialProperties[POISSON_RATIO];

		m_elastic_strain = ZeroVector(6);
		m_plastic_strain = ZeroVector(6);

		CalculateElasticityMatrix(m_D0);

		K = sqrt(2.0) * (m_f_02cc - m_f_01cc)/(2 * m_f_02cc - m_f_01cc);
		usedEquivalentTensionDefinition = COMPDYN;

		/**

		m_rate_biaxial_uniaxial = rMaterialProperties[RATE_BIAXIAL_UNIAXIAL];

		m_beta = rMaterialProperties[BETA];

		m_compression_parameter_A = rMaterialProperties[COMPRESSION_PARAMETER_A];
		m_compression_parameter_B = rMaterialProperties[COMPRESSION_PARAMETER_B];




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
		const double tolerance = (1.0e-14)* m_f_01cc;

		// 1.step: elastic stress tensor
		noalias(rStressVector) = prod(trans(m_D0), (rStrainVector - m_plastic_strain));

		// 2.step: spectral decomposition
		Vector StressVectorTension = ZeroVector(6);
		Vector StressVectorCompression = ZeroVector(6);
		Vector StressEigenvalues = ZeroVector(3);
		Matrix PMatrixTension = ZeroMatrix(6, 6);
		Matrix PMatrixCompression = ZeroMatrix(6, 6);

		SpectralDecomposition(rStressVector,
			StressVectorTension,
			StressVectorCompression,
			StressEigenvalues,
			PMatrixTension,
			PMatrixCompression);

		// 3.step: equivalent effective stress tau
		double tau_n;
		double tau_p;
		
		
		




		/**

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

	void TCPlasticDamage3DLaw::CalculateElasticityMatrix(
    	Matrix& rElasticityMatrix)
	{
    	const double lambda = m_E * m_nu / ((1.0 + m_nu) * (1.0 - 2 * m_nu));
    	const double mu = lambda / (2 * (1 + m_nu));

		if (rElasticityMatrix.size1() != 6 || rElasticityMatrix.size2() != 6)
			rElasticityMatrix.resize(6, 6, false);
		rElasticityMatrix = ZeroMatrix(6,6);

		rElasticityMatrix(0, 0) = 2 * mu + lambda;
		rElasticityMatrix(0, 1) = mu;
		rElasticityMatrix(0, 2) = mu;

		rElasticityMatrix(1, 0) = mu;
		rElasticityMatrix(1, 1) = 2 * mu + lambda;
		rElasticityMatrix(1, 2) = mu;

		rElasticityMatrix(2, 0) = mu;
		rElasticityMatrix(2, 1) = mu;
		rElasticityMatrix(2, 2) = 2 * mu + lambda;

		rElasticityMatrix(3, 3) = mu;

		rElasticityMatrix(4, 4) = mu;

		rElasticityMatrix(5, 5) = mu;
	};

	void TCPlasticDamage3DLaw::SpectralDecomposition(
		const Vector& rStressVector,
		Vector& rStressVectorTension,
		Vector& rStressVectorCompression,
		Vector& rStressEigenvalues,
		Matrix& PMatrixTension,
		Matrix& PMatrixCompression)
	{
		Matrix helpmatrix = ZeroMatrix(6, 6);

		/** necessary if SpectralDecomposition becomes a static member function
		rStressVectorTension = ZeroVector(6);
		rStressVectorCompression = ZeroVector(6);
		rStressEigenvalues = ZeroVector(3);
		PMatrixTension = ZeroMatrix(6, 6);
		PMatrixCompression = ZeroMatrix(6, 6);
		*/

		BoundedMatrix<double, 3, 3> Stress33;
		BoundedMatrix<double, 3, 3> EigenvalueStress33;
		BoundedMatrix<double, 3, 3> EigenvectorStress33;

		Stress33 = MathUtils<double>::StressVectorToTensor(rStressVector);

		bool converged = MathUtils<double>::EigenSystem<3>(Stress33, EigenvectorStress33, EigenvalueStress33);

		for (IndexType i = 0; i < 3; ++i){
			rStressEigenvalues(i) = EigenvalueStress33(i, i);
		}

		vector<Vector> p_i(3);
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

	void TCPlasticDamage3DLaw::ComputeTau(
		const Vector& rStressEigenvalues, 
		double& tau_n, 
		double& tau_p)
	{		
		Vector dgn(3);
		Vector dgp(3);
		for (IndexType i = 0; i < 3; i++)
			dgn(i) = (rStressEigenvalues(i) - fabs(rStressEigenvalues(i)))/2.0;	
		for (IndexType i = 0; i < 3; i++)
			dgp(i) = (rStressEigenvalues(i) + fabs(rStressEigenvalues(i)))/2.0;

		double sigoct = (dgn(0) + dgn(1) + dgn(2))/3.0;
		double tauoct = sqrt((dgn(0)-dgn(1)) * (dgn(0) - dgn(1)) 
		+ (dgn(0) - dgn(2)) * (dgn(0) - dgn(2)) 
		+ (dgn(1) - dgn(2)) * (dgn(1) - dgn(2)))/3.0;
	
		if ((usedEquivalentTensionDefinition == ORIGINAL)  || (usedEquivalentTensionDefinition == HOMOGENEOUS))
			{
			tau_n = sqrt(3.0) * (K * sigoct + tauoct);
			if (tau_n >= 0)
				tau_n = sqrt(tau_n);
			else
				tau_n = 0;
			}
		else if (usedEquivalentTensionDefinition == COMPDYN)
			{
			tau_n = sqrt(3.0) * (K * sigoct + tauoct);
			}

		double diag = (dgp(0) + dgp(1) + dgp(2)) * m_nu/(-m_E);
		Vector elastic_strain_diag_positive(3);
		for (IndexType i = 0; i < 3; i++)
			{
				elastic_strain_diag_positive(i) = dgp(i) * (1 + m_nu)/m_E + diag;
			}
	
		tau_p = 0.0;
		for (IndexType i = 0; i < 3; i++)
			{
				tau_p += elastic_strain_diag_positive(i) * dgp(i);
			}
	
		if (usedEquivalentTensionDefinition == ORIGINAL)
			{
			tau_p = sqrt(tau_p);
			}
		else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
			{
			tau_p = sqrt(sqrt(tau_p * m_E));
			}
		else if (usedEquivalentTensionDefinition == COMPDYN)
			{
			tau_p = sqrt(tau_p * m_E);
			}
	};
}
