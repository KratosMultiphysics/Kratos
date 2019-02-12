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

// System includes

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
		const double blabla=1;         // (ML)
    	KRATOS_WATCH(blabla)        // (ML)
		m_uniaxial_compressive_strength = rMaterialProperties[UNIAXIAL_STRESS_COMPRESSION];
		m_tensile_strength = rMaterialProperties[UNIAXIAL_STRESS_TENSION];
		m_biaxial_compressive_strength = rMaterialProperties[BIAXIAL_STRESS_COMPRESSION];
		m_E = rMaterialProperties[YOUNG_MODULUS];
		m_nu = rMaterialProperties[POISSON_RATIO];
		m_Gf = rMaterialProperties[FRACTURE_ENERGY_TENSION];

		CalculateElasticityMatrix(m_D0);

		m_elastic_strain = ZeroVector(6);
		m_plastic_strain = ZeroVector(6);
		
		m_beta = 0;			// method so far implemented neglected plastic strain (ML)

		m_K = sqrt(2.0) * (m_biaxial_compressive_strength - m_uniaxial_compressive_strength)/(2 * m_biaxial_compressive_strength - m_uniaxial_compressive_strength);
		usedEquivalentEffectiveStressDefinition = 2;

		// D+D- Damage Constitutive laws variables
		if (usedEquivalentEffectiveStressDefinition == 1)
			{
				m_initial_damage_threshold_compression = sqrt(sqrt(3.0)*(m_K - sqrt(2.0))*m_uniaxial_compressive_strength / 3);
				m_initial_damage_threshold_tension = sqrt(m_tensile_strength / sqrt(m_E));
			}
		else if (usedEquivalentEffectiveStressDefinition == 3)
			{
				m_initial_damage_threshold_compression = sqrt(sqrt(3.0)*(m_K - sqrt(2.0))*m_uniaxial_compressive_strength / 3);
				m_initial_damage_threshold_tension = sqrt(m_tensile_strength);
			}
		else if (usedEquivalentEffectiveStressDefinition == 2)
			{
				m_initial_damage_threshold_compression = sqrt(3.0)*(m_K - sqrt(2.0))*m_uniaxial_compressive_strength / 3;
				m_initial_damage_threshold_tension = m_tensile_strength;
			}
		m_damage_compression = 0.0;
		m_damage_tension = 0.0;
		m_damage_threshold_compression = m_initial_damage_threshold_compression;
		m_damage_threshold_compression1 = m_damage_threshold_compression;
		m_damage_threshold_tension = m_initial_damage_threshold_tension;
		m_damage_threshold_tension1 = m_damage_threshold_tension;
		m_compression_parameter_A = rMaterialProperties[COMPRESSION_PARAMETER_A];
		m_compression_parameter_B = rMaterialProperties[COMPRESSION_PARAMETER_B];
		/** be aware: area-function generally not perfectly implemented, returns volume for 3D-elements 
		and is not defined in IGA-Application (ML) */
		double l_c = rElementGeometry.Area();
		KRATOS_WATCH(l_c)
		m_tension_parameter_A = 1 / ((1 - m_beta)*((m_Gf * m_E / (l_c * m_tensile_strength * m_tensile_strength)) - 0.5));

		m_strain_ref = rMaterialProperties[STRAIN_REFERENCE];
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void TCPlasticDamage3DLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
	{
		this->CalculateMaterialResponseCauchy(rValues);
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void TCPlasticDamage3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
	{
		this->CalculateMaterialResponseCauchy(rValues);
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void TCPlasticDamage3DLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
	{
		this->CalculateMaterialResponseCauchy(rValues);
	}
	
	/***********************************************************************************/
	/***********************************************************************************/
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

/***********************************************************************************/
	/***********************************************************************************/
	void TCPlasticDamage3DLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
	{
		BaseType::FinalizeMaterialResponsePK1(rValues);
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void TCPlasticDamage3DLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
	{
		BaseType::FinalizeMaterialResponsePK2(rValues);
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void TCPlasticDamage3DLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
	{
		BaseType::FinalizeMaterialResponseKirchhoff(rValues);
	}

	/***********************************************************************************/
	/***********************************************************************************/
	void TCPlasticDamage3DLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
	{
		BaseType::FinalizeMaterialResponseCauchy(rValues);
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
		const double tolerance = (1.0e-14)* m_uniaxial_compressive_strength;		//move to function "DamageCriterion" if only used once (ML)
		KRATOS_WATCH(tolerance)
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
		double EquEffStressCompression;
		double EquEffStressTension;
		
		ComputeTau(StressEigenvalues,
			EquEffStressCompression, 
			EquEffStressTension);

		// 4.step: check damage criterion and update damage threshold
		DamageCriterion(EquEffStressCompression, 
			EquEffStressTension,
			tolerance);
		
		// 5.step: compute damage variables
		ComputeDamageCompression();

		ComputeDamageTension();

		// 6.step: compute shear retention factors
		ComputeSRF(rStrainVector);

		// 7.step: update stiffness matrix (damaged system)
		rConstitutiveLaw = prod(((1 - m_damage_tension) * PMatrixTension + (1 - m_damage_compression) * PMatrixCompression), m_D0);

		// 8.step: compute Cauchy stress (damaged system)
		rStressVector = (1 - m_damage_tension) * StressVectorTension + (1 - m_damage_compression) * StressVectorCompression;
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
		Matrix& rPMatrixTension,
		Matrix& rPMatrixCompression)
	{
		Matrix helpmatrix = ZeroMatrix(6, 6);

		/** necessary if SpectralDecomposition becomes a static member function
		rStressVectorTension = ZeroVector(6);
		rStressVectorCompression = ZeroVector(6);
		rStressEigenvalues = ZeroVector(3);
		rPMatrixTension = ZeroMatrix(6, 6);
		rPMatrixCompression = ZeroMatrix(6, 6);
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

		rPMatrixTension = ZeroMatrix(6, 6);
		for (SizeType i = 0.0; i < 3; ++i)
		{
			if (EigenvalueStress33(i, i) > 0.0)
			{
				Vector helpvector = MathUtils<double>::StressTensorToVector(helpmatrix);
				rPMatrixTension += outer_prod(helpvector, helpvector);
			}
		}
		Matrix I = IdentityMatrix(6, 6);
		rPMatrixCompression = I - rPMatrixTension;
	};

	void TCPlasticDamage3DLaw::ComputeTau(
		const Vector& rStressEigenvalues, 
		double& rEquEffStressCompression, 
		double& rEquEffStressTension)
	{		
		Vector eigenvalues_compression(3);
		Vector eigenvalues_tension(3);
		for (IndexType i = 0; i < 3; i++)
			eigenvalues_compression(i) = (rStressEigenvalues(i) - fabs(rStressEigenvalues(i)))/2.0;	
		for (IndexType i = 0; i < 3; i++)
			eigenvalues_tension(i) = (rStressEigenvalues(i) + fabs(rStressEigenvalues(i)))/2.0;

		double sigoct = (eigenvalues_compression(0) + eigenvalues_compression(1) + eigenvalues_compression(2))/3.0;
		double tauoct = sqrt((eigenvalues_compression(0)-eigenvalues_compression(1)) * (eigenvalues_compression(0) - eigenvalues_compression(1)) 
		+ (eigenvalues_compression(0) - eigenvalues_compression(2)) * (eigenvalues_compression(0) - eigenvalues_compression(2)) 
		+ (eigenvalues_compression(1) - eigenvalues_compression(2)) * (eigenvalues_compression(1) - eigenvalues_compression(2)))/3.0;
	
		if ((usedEquivalentEffectiveStressDefinition == 1)  || (usedEquivalentEffectiveStressDefinition == 3))
			{
			rEquEffStressCompression = sqrt(3.0) * (m_K * sigoct + tauoct);
			if (rEquEffStressCompression >= 0)
				rEquEffStressCompression = sqrt(rEquEffStressCompression);
			else
				rEquEffStressCompression = 0;
			}
		else if (usedEquivalentEffectiveStressDefinition == 2)
			{
			rEquEffStressCompression = sqrt(3.0) * (m_K * sigoct + tauoct);
			}

		double diag = (eigenvalues_tension(0) + eigenvalues_tension(1) + eigenvalues_tension(2)) * m_nu/(-m_E);
		Vector elastic_strain_diag_positive(3);
		for (IndexType i = 0; i < 3; i++)
			{
				elastic_strain_diag_positive(i) = eigenvalues_tension(i) * (1 + m_nu)/m_E + diag;
			}
	
		rEquEffStressTension = 0.0;
		for (IndexType i = 0; i < 3; i++)
			{
				rEquEffStressTension += elastic_strain_diag_positive(i) * eigenvalues_tension(i);
			}
	
		if (usedEquivalentEffectiveStressDefinition == 1)
			{
			rEquEffStressTension = sqrt(rEquEffStressTension);
			}
		else if (usedEquivalentEffectiveStressDefinition == 3)
			{
			rEquEffStressTension = sqrt(sqrt(rEquEffStressTension * m_E));
			}
		else if (usedEquivalentEffectiveStressDefinition == 2)
			{
			rEquEffStressTension = sqrt(rEquEffStressTension * m_E);
			}
	};

	void TCPlasticDamage3DLaw::DamageCriterion(
		const double& rEquEffStressCompression, 
		const double& rEquEffStressTension,
		const double& rtolerance)
	{
		double g = (rEquEffStressTension/m_damage_threshold_tension)*(rEquEffStressTension/m_damage_threshold_tension)+(rEquEffStressCompression/m_damage_threshold_compression)*(rEquEffStressCompression/m_damage_threshold_compression)-1;
		
		if (g > rtolerance)
			{
			double rhoQ = sqrt(rEquEffStressTension*rEquEffStressTension+rEquEffStressCompression*rEquEffStressCompression);
			double thetaQ;
			if (rEquEffStressCompression > 1e-15) 
				{
				thetaQ=atan(rEquEffStressTension/rEquEffStressCompression);
				} 
			else 
				{
				thetaQ=std::atan(1)*4/2.0;		//atan(1)*4=Pi -> replace maybe (ML)
				}
			double rhoP = m_damage_threshold_tension*m_damage_threshold_compression*sqrt((rEquEffStressCompression*rEquEffStressCompression+rEquEffStressTension*rEquEffStressTension)/((rEquEffStressCompression*m_damage_threshold_tension)*(rEquEffStressCompression*m_damage_threshold_tension)+(rEquEffStressTension*m_damage_threshold_compression)*(rEquEffStressTension*m_damage_threshold_compression)));
			if (m_damage_threshold_compression >= m_damage_threshold_tension)
				{
				if (rhoP<m_damage_threshold_tension)
					rhoP =m_damage_threshold_tension;
				if (rhoP>m_damage_threshold_compression)
					rhoP = m_damage_threshold_compression;
				}
			else if (m_damage_threshold_compression < m_damage_threshold_tension)
				{
				if (rhoP>m_damage_threshold_tension)
					rhoP =m_damage_threshold_tension;
				if (rhoP<m_damage_threshold_compression)
					rhoP = m_damage_threshold_compression;
				}
			double alfa=rhoQ/rhoP;
			double thetaL = atan((m_damage_threshold_tension*m_damage_threshold_tension)/(m_damage_threshold_compression*m_damage_threshold_compression));
			double rhoL=sqrt((m_damage_threshold_tension*m_damage_threshold_tension*m_damage_threshold_compression*m_damage_threshold_compression)/(m_damage_threshold_compression*m_damage_threshold_compression*sin(thetaL)*sin(thetaL)+m_damage_threshold_tension*m_damage_threshold_tension*cos(thetaL)*cos(thetaL)));
			double alfasn;
			if (((rhoP>rhoL) && (rhoP<=m_damage_threshold_compression)) || ((rhoP>=m_damage_threshold_compression) && (rhoP<rhoL))) 
				{
					double alfasp=1+(alfa-1)*(m_damage_threshold_compression-rhoP)/(m_damage_threshold_compression-rhoL);
					m_damage_threshold_tension1=m_damage_threshold_tension*alfasp;
					m_damage_threshold_compression1=sqrt((m_damage_threshold_tension1*m_damage_threshold_tension1*rEquEffStressCompression*rEquEffStressCompression)/(m_damage_threshold_tension1*m_damage_threshold_tension1-rEquEffStressTension*rEquEffStressTension));
				} 
			else if (((rhoP>rhoL) && (rhoP<=m_damage_threshold_tension)) || ((rhoP>=m_damage_threshold_tension) && (rhoP<rhoL))) 
				{
				alfasn=1+(alfa-1)*(rhoP-m_damage_threshold_tension)/(rhoL-m_damage_threshold_tension);
				m_damage_threshold_compression1=m_damage_threshold_compression*alfasn;
				m_damage_threshold_tension1=sqrt((m_damage_threshold_compression1*m_damage_threshold_compression1*rEquEffStressTension*rEquEffStressTension)/(m_damage_threshold_compression1*m_damage_threshold_compression1-rEquEffStressCompression*rEquEffStressCompression));
				} 				
			else 
				{
				alfasn=1+(alfa-1)*(rhoP-m_damage_threshold_tension)/(rhoL-m_damage_threshold_tension);
				m_damage_threshold_compression1=m_damage_threshold_compression*alfasn;
				m_damage_threshold_tension1=sqrt((m_damage_threshold_compression1*m_damage_threshold_compression1*rEquEffStressTension*rEquEffStressTension)/(m_damage_threshold_compression1*m_damage_threshold_compression1-rEquEffStressCompression*rEquEffStressCompression));
				}
			} 
		else
			{
			m_damage_threshold_compression1=m_damage_threshold_compression;
			m_damage_threshold_tension1=m_damage_threshold_tension;
			};
	};

	void TCPlasticDamage3DLaw::ComputeDamageCompression() 
	{
		double rDamageCompression;
		if (m_damage_threshold_compression1 < 1e-7)
        	rDamageCompression = 0.0;
		else
			{
				if ((usedEquivalentEffectiveStressDefinition == 1)  || (usedEquivalentEffectiveStressDefinition == 3))
					rDamageCompression = 1 - m_initial_damage_threshold_compression/m_damage_threshold_compression1 * (1 - m_compression_parameter_A) - m_compression_parameter_A * exp(m_compression_parameter_B * (1 - m_damage_threshold_compression1/m_initial_damage_threshold_compression));
				else if (usedEquivalentEffectiveStressDefinition == 2)
					rDamageCompression = 1 - (sqrt(m_initial_damage_threshold_compression))/(sqrt(m_damage_threshold_compression1)) * (1 - m_compression_parameter_A) - m_compression_parameter_A * exp(m_compression_parameter_B * (1 - (sqrt(m_damage_threshold_compression1)/(sqrt(m_initial_damage_threshold_compression)))));
    		    // limiting damage variable
				rDamageCompression = std::min(rDamageCompression, 1.0 - 1e-7);		// maximum not equal to 1.0 since then the stiffness matrix would be equal to 0.0
        		rDamageCompression = std::max(rDamageCompression, 0.0);
				// update damage variable
				m_damage_compression = std::max(m_damage_compression,rDamageCompression);
			}	
	}
	
	void TCPlasticDamage3DLaw::ComputeDamageTension()
	{	
		double rDamageTension;	
		if (m_damage_threshold_compression1 < 1e-7)
        	rDamageTension = 0.0;
		else
			{
				if ((usedEquivalentEffectiveStressDefinition == 1)  || (usedEquivalentEffectiveStressDefinition == 2))
					rDamageTension = 1 - ((m_initial_damage_threshold_tension * m_initial_damage_threshold_tension)/(m_damage_threshold_tension1 * m_damage_threshold_tension1)) * exp(m_tension_parameter_A * (1 - (m_damage_threshold_tension1 * m_damage_threshold_tension1)/(m_initial_damage_threshold_tension * m_initial_damage_threshold_tension)));
				else if (usedEquivalentEffectiveStressDefinition == 3)
					rDamageTension = 1 - ((m_initial_damage_threshold_tension)/(m_damage_threshold_tension1)) * exp(m_tension_parameter_A * (1 - (m_damage_threshold_tension1)/(m_initial_damage_threshold_tension)));
				// limiting damage variable
				rDamageTension = std::min(rDamageTension, 1.0 - 1e-7);
        		rDamageTension = std::max(rDamageTension, 0.0);
				// update damage variable
				m_damage_tension = std::max(m_damage_tension,rDamageTension);
			}
	}

	void TCPlasticDamage3DLaw::ComputeSRF(const Vector& rStrainVector)
	{
		if (m_strain_ref <= 0.0) 
		{
			m_SRF12 = 0.0;
			m_SRF13 = 0.0;
			m_SRF23 = 0.0;
		}
		else 
		{
			// evolution law for SRF according to Scotta(2001)
			m_SRF12 = 1-abs(rStrainVector(3))/m_strain_ref;
			m_SRF13 = 1-abs(rStrainVector(4))/m_strain_ref;
			m_SRF23 = 1-abs(rStrainVector(5))/m_strain_ref;
			if (m_SRF12 <= 0.0)
				m_SRF12 = 0.0;
			if (m_SRF13 <= 0.0)
				m_SRF13 = 0.0;
			if (m_SRF23 <= 0.0)
				m_SRF23 = 0.0;
		}
	}

	//add formula such as "commitState" in Code from Padua (ML)

} // namespace Kratos
