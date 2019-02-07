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
		m_f_01cc = rMaterialProperties[UNIAXIAL_STRESS_COMPRESSION];
		m_f_ct = rMaterialProperties[UNIAXIAL_STRESS_TENSION];
		m_f_02cc = rMaterialProperties[BIAXIAL_STRESS_COMPRESSION];
		m_E = rMaterialProperties[YOUNG_MODULUS];
		m_nu = rMaterialProperties[POISSON_RATIO];
		m_Gf = rMaterialProperties[FRACTURE_ENERGY_TENSION];

		CalculateElasticityMatrix(m_D0);

		m_elastic_strain = ZeroVector(6);
		m_plastic_strain = ZeroVector(6);
		
		m_beta = 0;			// method so far implemented neglected plastic strain (ML)

		m_K = sqrt(2.0) * (m_f_02cc - m_f_01cc)/(2 * m_f_02cc - m_f_01cc);
		usedEquivalentEffectiveStressDefinition = 2;

		// D+D- Damage Constitutive laws variables
		if (usedEquivalentEffectiveStressDefinition == 1)
			{
				m_r_0n = sqrt(sqrt(3.0)*(m_K - sqrt(2.0))*m_f_01cc / 3);
				m_r_0p = sqrt(m_f_ct / sqrt(m_E));
			}
		else if (usedEquivalentEffectiveStressDefinition == 3)
			{
				m_r_0n = sqrt(sqrt(3.0)*(m_K - sqrt(2.0))*m_f_01cc / 3);
				m_r_0p = sqrt(m_f_ct);
			}
		else if (usedEquivalentEffectiveStressDefinition == 2)
			{
				m_r_0n = sqrt(3.0)*(m_K - sqrt(2.0))*m_f_01cc / 3;
				m_r_0p = m_f_ct;
			}
		m_d_n = 0.0;
		m_d_p = 0.0;
		m_r_n = m_r_0n;
		m_r_n1 = m_r_n;
		m_r_p = m_r_0p;
		m_r_p1 = m_r_p;
		m_A_n = rMaterialProperties[COMPRESSION_PARAMETER_A];
		m_B_n = rMaterialProperties[COMPRESSION_PARAMETER_B];
		/** be aware: area-function generally not perfectly implemented, returns volume for 3D-elements 
		and is not defined in IGA-Application (ML) */
		double l_c = rElementGeometry.Area();
		KRATOS_WATCH(l_c)
		m_A_p = 1 / ((1 - m_beta)*((m_Gf * m_E / (l_c * m_f_ct * m_f_ct)) - 0.5));

		m_strain_ref = rMaterialProperties[STRAIN_REFERENCE];
		/**


		m_Di = m_D0;

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
		const double tolerance = (1.0e-14)* m_f_01cc;		//move to function "DamageCriterion" if only used once (ML)
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
		double tau_n;
		double tau_p;
		
		ComputeTau(StressEigenvalues,
			tau_n, 
			tau_p);

		// 4.step: check damage criterion and update damage threshold
		DamageCriterion(tau_n, 
			tau_p,
			tolerance);
		
		// 5.step: compute damage variables
		ComputeDamageCompression();

		ComputeDamageTension();

		// 6.step: compute shear retention factors
		ComputeSRF(rStrainVector);

		// 7.step: update stiffness matrix (damaged system)
		rConstitutiveLaw = prod(((1 - m_d_p) * PMatrixTension + (1 - m_d_n) * PMatrixCompression), m_D0);

		// 8.step: compute Cauchy stress (damaged system)
		rStressVector = (1 - m_d_p) * StressVectorTension + (1 - m_d_n) * StressVectorCompression;
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
		double& rtau_n, 
		double& rtau_p)
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
	
		if ((usedEquivalentEffectiveStressDefinition == 1)  || (usedEquivalentEffectiveStressDefinition == 3))
			{
			rtau_n = sqrt(3.0) * (m_K * sigoct + tauoct);
			if (rtau_n >= 0)
				rtau_n = sqrt(rtau_n);
			else
				rtau_n = 0;
			}
		else if (usedEquivalentEffectiveStressDefinition == 2)
			{
			rtau_n = sqrt(3.0) * (m_K * sigoct + tauoct);
			}

		double diag = (dgp(0) + dgp(1) + dgp(2)) * m_nu/(-m_E);
		Vector elastic_strain_diag_positive(3);
		for (IndexType i = 0; i < 3; i++)
			{
				elastic_strain_diag_positive(i) = dgp(i) * (1 + m_nu)/m_E + diag;
			}
	
		rtau_p = 0.0;
		for (IndexType i = 0; i < 3; i++)
			{
				rtau_p += elastic_strain_diag_positive(i) * dgp(i);
			}
	
		if (usedEquivalentEffectiveStressDefinition == 1)
			{
			rtau_p = sqrt(rtau_p);
			}
		else if (usedEquivalentEffectiveStressDefinition == 3)
			{
			rtau_p = sqrt(sqrt(rtau_p * m_E));
			}
		else if (usedEquivalentEffectiveStressDefinition == 2)
			{
			rtau_p = sqrt(rtau_p * m_E);
			}
	};

	void TCPlasticDamage3DLaw::DamageCriterion(
		const double& rtau_n, 
		const double& rtau_p,
		const double& rtolerance)
	{
		double g = (rtau_p/m_r_p)*(rtau_p/m_r_p)+(rtau_n/m_r_n)*(rtau_n/m_r_n)-1;
		
		if (g > rtolerance)
			{
			double rhoQ = sqrt(rtau_p*rtau_p+rtau_n*rtau_n);
			double thetaQ;
			if (rtau_n > 1e-15) 
				{
				thetaQ=atan(rtau_p/rtau_n);
				} 
			else 
				{
				thetaQ=std::atan(1)*4/2.0;		//atan(1)*4=Pi -> replace maybe (ML)
				}
			double rhoP = m_r_p*m_r_n*sqrt((rtau_n*rtau_n+rtau_p*rtau_p)/((rtau_n*m_r_p)*(rtau_n*m_r_p)+(rtau_p*m_r_n)*(rtau_p*m_r_n)));
			if (m_r_n >= m_r_p)
				{
				if (rhoP<m_r_p)
					rhoP =m_r_p;
				if (rhoP>m_r_n)
					rhoP = m_r_n;
				}
			else if (m_r_n < m_r_p)
				{
				if (rhoP>m_r_p)
					rhoP =m_r_p;
				if (rhoP<m_r_n)
					rhoP = m_r_n;
				}
			double alfa=rhoQ/rhoP;
			double thetaL = atan((m_r_p*m_r_p)/(m_r_n*m_r_n));
			double rhoL=sqrt((m_r_p*m_r_p*m_r_n*m_r_n)/(m_r_n*m_r_n*sin(thetaL)*sin(thetaL)+m_r_p*m_r_p*cos(thetaL)*cos(thetaL)));
			double alfasn;
			if (((rhoP>rhoL) && (rhoP<=m_r_n)) || ((rhoP>=m_r_n) && (rhoP<rhoL))) 
				{
					double alfasp=1+(alfa-1)*(m_r_n-rhoP)/(m_r_n-rhoL);
					m_r_p1=m_r_p*alfasp;
					m_r_n1=sqrt((m_r_p1*m_r_p1*rtau_n*rtau_n)/(m_r_p1*m_r_p1-rtau_p*rtau_p));
				} 
			else if (((rhoP>rhoL) && (rhoP<=m_r_p)) || ((rhoP>=m_r_p) && (rhoP<rhoL))) 
				{
				alfasn=1+(alfa-1)*(rhoP-m_r_p)/(rhoL-m_r_p);
				m_r_n1=m_r_n*alfasn;
				m_r_p1=sqrt((m_r_n1*m_r_n1*rtau_p*rtau_p)/(m_r_n1*m_r_n1-rtau_n*rtau_n));
				} 				
			else 
				{
				alfasn=1+(alfa-1)*(rhoP-m_r_p)/(rhoL-m_r_p);
				m_r_n1=m_r_n*alfasn;
				m_r_p1=sqrt((m_r_n1*m_r_n1*rtau_p*rtau_p)/(m_r_n1*m_r_n1-rtau_n*rtau_n));
				}
			} 
		else
			{
			m_r_n1=m_r_n;
			m_r_p1=m_r_p;
			};
	};

	void TCPlasticDamage3DLaw::ComputeDamageCompression() 
	{
		double rd_n;
		if (m_r_n1 < 1e-7)
        	rd_n = 0.0;
		else
			{
				if ((usedEquivalentEffectiveStressDefinition == 1)  || (usedEquivalentEffectiveStressDefinition == 3))
					rd_n = 1 - m_r_0n/m_r_n1 * (1 - m_A_n) - m_A_n * exp(m_B_n * (1 - m_r_n1/m_r_0n));
				else if (usedEquivalentEffectiveStressDefinition == 2)
					rd_n = 1 - (sqrt(m_r_0n))/(sqrt(m_r_n1)) * (1 - m_A_n) - m_A_n * exp(m_B_n * (1 - (sqrt(m_r_n1)/(sqrt(m_r_0n)))));
    		    // limiting damage variable
				rd_n = std::min(rd_n, 1.0 - 1e-7);		// maximum not equal to 1.0 since then the stiffness matrix would be equal to 0.0
        		rd_n = std::max(rd_n, 0.0);
				// update damage variable
				m_d_n = std::max(m_d_n,rd_n);
			}	
	}
	
	void TCPlasticDamage3DLaw::ComputeDamageTension()
	{	
		double rd_p;	
		if (m_r_n1 < 1e-7)
        	rd_p = 0.0;
		else
			{
				if ((usedEquivalentEffectiveStressDefinition == 1)  || (usedEquivalentEffectiveStressDefinition == 2))
					rd_p = 1 - ((m_r_0p * m_r_0p)/(m_r_p1 * m_r_p1)) * exp(m_A_p * (1 - (m_r_p1 * m_r_p1)/(m_r_0p * m_r_0p)));
				else if (usedEquivalentEffectiveStressDefinition == 3)
					rd_p = 1 - ((m_r_0p)/(m_r_p1)) * exp(m_A_p * (1 - (m_r_p1)/(m_r_0p)));
				// limiting damage variable
				rd_p = std::min(rd_p, 1.0 - 1e-7);
        		rd_p = std::max(rd_p, 0.0);
				// update damage variable
				m_d_p = std::max(m_d_p,rd_p);
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
