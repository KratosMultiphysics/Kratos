/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/* *********************************************************
*
*   Last Modified by:    $Author:   Massimo Petracca$
*   Date:                $Date:     30-10-2013$
*   Revision:            $Revision: 1.0$
*
* ***********************************************************/

#include "j2_constitutive_law_3d.h"
#include "multiscale_application_variables.h"

namespace Kratos
{

	J2ConstitutiveLaw3D::J2ConstitutiveLaw3D()
		: ConstitutiveLaw()
		, m_initialized(false)
		, m_error_code(0.0)
		, m_lch(0.0)
		, m_lch_multiplier(1.0)
		, m_init_strain()
	{
	}

	ConstitutiveLaw::Pointer J2ConstitutiveLaw3D::Clone() const
	{
		return ConstitutiveLaw::Pointer(new J2ConstitutiveLaw3D());
	}

	J2ConstitutiveLaw3D::SizeType J2ConstitutiveLaw3D::WorkingSpaceDimension()
	{
		return 3;
	}

	J2ConstitutiveLaw3D::SizeType J2ConstitutiveLaw3D::GetStrainSize()
	{
		return 6;
	}

	bool J2ConstitutiveLaw3D::Has(const Variable<double>& rThisVariable)
	{
		if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN) return true;
		return false;
	}

	bool J2ConstitutiveLaw3D::Has(const Variable<Vector>& rThisVariable)
	{
		if(rThisVariable == INITIAL_STRAIN) return true;
		return false;
	}

	bool J2ConstitutiveLaw3D::Has(const Variable<Matrix>& rThisVariable)
	{
		if(rThisVariable == PLASTIC_STRAIN_TENSOR) return true;
		return false;
	}

	bool J2ConstitutiveLaw3D::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	bool J2ConstitutiveLaw3D::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

	double& J2ConstitutiveLaw3D::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		rValue = 0.0;
		if (rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = m_error_code;
		if(rThisVariable == EQUIVALENT_PLASTIC_STRAIN)
			rValue = m_xi;
		return rValue;
	}

	Vector& J2ConstitutiveLaw3D::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		if (rThisVariable == INITIAL_STRAIN) {
			if (rValue.size() != m_init_strain.size())
				rValue.resize(m_init_strain.size());
			noalias(rValue) = m_init_strain;
		}
		return rValue;
	}

	Matrix& J2ConstitutiveLaw3D::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		if(rThisVariable == PLASTIC_STRAIN_TENSOR) {
			if(rValue.size1() != 1 || rValue.size2() != 6) rValue.resize(1, 6, false);
			for(int i = 0; i < 6; i++)
				rValue(0,i) = m_eps_pl(i);
		}
		return rValue;
	}

	array_1d<double, 3 > & J2ConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & J2ConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void J2ConstitutiveLaw3D::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			m_lch_multiplier = rValue;
	}

	void J2ConstitutiveLaw3D::SetValue(const Variable<Vector >& rVariable,
		const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == INITIAL_STRAIN) {
			if (rValue.size() == m_init_strain.size())
				noalias(m_init_strain) = rValue;
		}
	}

	void J2ConstitutiveLaw3D::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void J2ConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void J2ConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool J2ConstitutiveLaw3D::ValidateInput(const Properties& rMaterialProperties)
	{
		if(!rMaterialProperties.Has(YOUNG_MODULUS)) return false;
		if(!rMaterialProperties.Has(POISSON_RATIO)) return false;
		if(!rMaterialProperties.Has(DENSITY)) return false;
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING)) return false;
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING_EXPONENT)) return false;
		if(!rMaterialProperties.Has(KINEMATIC_HARDENING)) return false;
		if(!rMaterialProperties.Has(YIELD_STRESS)) return false;
		if(!rMaterialProperties.Has(YIELD_STRESS_INFINITY)) return false;
		if(!rMaterialProperties.Has(FRACTURE_ENERGY_T)) return false;
		if(!rMaterialProperties.Has(VISCOSITY)) return false;
		return true;
	}

	J2ConstitutiveLaw3D::StrainMeasure J2ConstitutiveLaw3D::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	J2ConstitutiveLaw3D::StressMeasure J2ConstitutiveLaw3D::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool J2ConstitutiveLaw3D::IsIncremental()
	{
		return false;
	}

	void J2ConstitutiveLaw3D::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if(!m_initialized)
		{
			m_error_code = 0.0;
			m_xi = 0.0;
			m_xi_converged = 0.0;
			m_eps_pl.clear();
			m_eps_pl_converged.clear();
			m_lch = rElementGeometry.Length();
			m_lch_multiplier = 1.0;
			m_init_strain = ZeroVector(this->GetStrainSize());
			m_initialized = true;
		}
	}

	void J2ConstitutiveLaw3D::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		m_xi = m_xi_converged;
		m_eps_pl = m_eps_pl_converged;
	}

	void J2ConstitutiveLaw3D::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		m_xi_converged = m_xi;
		m_eps_pl_converged = m_eps_pl;
	}

	void J2ConstitutiveLaw3D::InitializeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void J2ConstitutiveLaw3D::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void J2ConstitutiveLaw3D::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void J2ConstitutiveLaw3D::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void J2ConstitutiveLaw3D::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void J2ConstitutiveLaw3D::CalculateMaterialResponseCauchy (Parameters& rValues)
	{
		m_error_code = 0.0;
		// get some references
		const Properties& props = rValues.GetMaterialProperties();
		//const Vector& eps = rValues.GetStrainVector();
		Vector eps = rValues.GetStrainVector();

		// Subtract Initial_Strain
		noalias(eps) -= m_init_strain;

		// Subtract Temperature_Strain
		double DeltaTemp = 0.0;
		Vector temp_strain(eps.size(), 0.0);
		CalculateStrainTemperature(rValues, props, DeltaTemp, temp_strain);
		noalias(eps) -= temp_strain;

		// get a reference of the output stress vector and tangent operator
		// and initialize them
		Vector& stressVector = rValues.GetStressVector();
		Matrix& tangentMatrix = rValues.GetConstitutiveMatrix();
		if(stressVector.size() != 6)
			stressVector.resize(6, false);
		noalias(stressVector) = ZeroVector(6);
		if(tangentMatrix.size1() != 6 || tangentMatrix.size2() != 6)
			tangentMatrix.resize(6, 6, false);
		noalias(tangentMatrix) = ZeroMatrix(6, 6);

		// get material parameters
		double E    = props[YOUNG_MODULUS];
		double ni   = props[POISSON_RATIO];
		double Hi   = props[ISOTROPIC_HARDENING];
		double Hk   = props[KINEMATIC_HARDENING];
		double sy   = props[YIELD_STRESS];
		double sinf = props[YIELD_STRESS_INFINITY];
		double Gf   = props[FRACTURE_ENERGY_T];
		double eta  = props[VISCOSITY];
		// elastic constants
		double G    = E/(2.0*(1.0+ni));
		double K    = E/(3.0*(1.0-2.0*ni));
		// trial state
		m_xi = m_xi_converged;
		m_eps_pl = m_eps_pl_converged;

		// compute the viscosity parameter eta/dT
		// avoid numerical errors
		double dt = rValues.GetProcessInfo()[DELTA_TIME];
		double eta_over_dt = dt > 0.0 ? eta/dt : 0.0;

		// Compute volumetric and deviatoric strain split
		double trace_eps = eps(0) + eps(1) + eps(2);

		array_1d<double, 6> eps_dev;
		eps_dev(0) = eps(0) - trace_eps / 3.0;
		eps_dev(1) = eps(1) - trace_eps / 3.0;
		eps_dev(2) = eps(2) - trace_eps / 3.0;
		eps_dev(3) = eps(3)/2.0;
		eps_dev(4) = eps(4)/2.0;
		eps_dev(5) = eps(5)/2.0;

		// trial deviatoric stress
		array_1d<double, 6> sig_dev;
		noalias( sig_dev ) = 2.0 * G * (eps_dev  - m_eps_pl);

		// trial relative stress and its norm
		array_1d<double, 6> sig_rel;
		noalias( sig_rel ) = sig_dev - 2.0/3.0*Hk * m_eps_pl;
		double norm_sig_rel = std::sqrt(
			sig_rel(0)*sig_rel(0) +
			sig_rel(1)*sig_rel(1) +
			sig_rel(2)*sig_rel(2) +
			2.0*sig_rel(3)*sig_rel(3) +
			2.0*sig_rel(4)*sig_rel(4) +
			2.0*sig_rel(5)*sig_rel(5)
			);

		// Initialize the tangent operator to Zero
		//
		// Note: initialize delta1 to 1 so that in case of no plastic corrector
		// the finalization of the tangent operator will lead to the elastic
		// constitutive matrix
		double delta1 = 1.0;

		// Evaluate the Yield Function

		double F = norm_sig_rel - std::sqrt(2.0/3.0)*HardeningFunction(m_xi, sy, sinf, Hi, Gf);

		if(F > 0.0) // Plastic corrector
		{
			const double tolerance = 1.0E-4;
			const int max_iter = 1000;
			double gamma = 0.0;
			double resid = 1.0;
			int iter = 0;

			for(iter = 0; iter < max_iter; iter++)
			{
				resid = norm_sig_rel
					-(2.0*G + 2.0/3.0*Hk)*gamma
					- std::sqrt(2.0/3.0)*HardeningFunction(m_xi + std::sqrt(2.0/3.0)*gamma, sy, sinf, Hi, Gf)
					- eta_over_dt*gamma;

				if(resid < tolerance)
					break;

				double tangent = -2.0*G
					- 2.0/3.0 * Hk
					- 2.0/3.0 * HardeningDerivative(m_xi + std::sqrt(2.0/3.0)*gamma, sy, sinf, Hi, Gf)
					- eta_over_dt;

				double dgamma = -resid/tangent;
				gamma += dgamma;
			}

			if(iter == max_iter) {
				m_error_code = -1.0;
				std::cout << "J2 Constitutive Law 3D - Integration algorithm did not converge\n";
			}

			m_xi += gamma*sqrt(2.0/3.0);  // accumulated equivalent plastic strain
			array_1d<double, 6> N = sig_rel/norm_sig_rel;   // normal to the yield surface

			m_eps_pl += gamma*N;       // updated plastic strains
			sig_dev -= gamma*2.0*G*N;  // updated deviatoric stresses

			// compute the parameters for the tangent operator
			// NOTE: the elastic part will be added later on
			delta1 = 1.0-2.0*G*gamma/norm_sig_rel;
			double delta2 = 1.0/(1.0
				+1.0/(3.0*G)*(Hk+HardeningDerivative(m_xi, sy, sinf, Hi, Gf))
				+eta_over_dt/(2.0*G)
				)-(1.0-delta1);
			tangentMatrix -= 2.0*G*delta2*outer_prod(N,N);
		}

		// Finalize the computation of the total stress vector at t(n+1)
		double pressure = K * trace_eps;
		stressVector(0) = sig_dev(0) + pressure;
		stressVector(1) = sig_dev(1) + pressure;
		stressVector(2) = sig_dev(2) + pressure;
		stressVector(3) = sig_dev(3);
		stressVector(4) = sig_dev(4);
		stressVector(5) = sig_dev(5);

		// Finalize the computation of the tangent operator at t(n+1)
		array_1d<double, 6> Id2; // 2nd-order identity tensor
		Id2(0) = 1.0; Id2(1) = 1.0; Id2(2) = 1.0; Id2(3) = 0.0; Id2(4) = 0.0; Id2(5) = 0.0;
		Matrix Id4(IdentityMatrix(6,6)); // 4th-order identity tensor
		Id4(3,3) /= 2.0; Id4(4,4) /= 2.0; Id4(5,5) /= 2.0;
		Matrix ItensI( outer_prod(Id2, Id2) );
		noalias(tangentMatrix) += K*ItensI;
		noalias(tangentMatrix) += 2.0*G*delta1*(Id4-ItensI/3.0);
	}

	void J2ConstitutiveLaw3D::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void J2ConstitutiveLaw3D::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void J2ConstitutiveLaw3D::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void J2ConstitutiveLaw3D::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{
	}

	void J2ConstitutiveLaw3D::ResetMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if(m_initialized)
		{
			m_error_code = 0.0;
			m_xi = 0.0;
			m_xi_converged = 0.0;
			m_eps_pl.clear();
			m_eps_pl_converged.clear();
			m_lch = 0.0;
			m_lch_multiplier = 1.0;
			m_initialized = false;
		}
	}

	//******************************* COMPUTE DOMAIN TEMPERATURE  ************************
	//************************************************************************************


	double &  J2ConstitutiveLaw3D::CalculateDomainTemperature(Parameters& rValues,
		double& rDeltaTemperature)
	{

		//1.-Temperature from nodes
		const GeometryType& DomainGeometry = rValues.GetElementGeometry();
		//GeometryType aaa = DomainGeometry;
		const Vector& ShapeFunctionsValues = rValues.GetShapeFunctionsValues();
		const unsigned int number_of_nodes = DomainGeometry.size();
		double ambient_T = 0.0; // Set to Zero by default

		double iterpolated_temp = 0.0;

		for (unsigned int j = 0; j < number_of_nodes; j++)
		{
			if (DomainGeometry[j].SolutionStepsDataHas(TEMPERATURE))
			{
				if (DomainGeometry[j].SolutionStepsDataHas(RVE_FULL_TEMPERATURE))
				{
					//std::cout << "Is RVE_FULL_TEMPERATURE" << std::endl;
					iterpolated_temp += ShapeFunctionsValues[j] * DomainGeometry[j].FastGetSolutionStepValue(RVE_FULL_TEMPERATURE);
				}
				else
				{
					//std::cout << "Is TEMPERATURE" << std::endl;
					iterpolated_temp += ShapeFunctionsValues[j] * DomainGeometry[j].FastGetSolutionStepValue(TEMPERATURE);
				}
			}
			else
			{
				//std::cout << "IsNotTemperature" << std::endl;
				iterpolated_temp = 0.0;
			}
		}

		rDeltaTemperature = iterpolated_temp - ambient_T;

		//std::cout << TEMPERATURE << std::endl;
		//std::cout << "------------MEC-------------";
		//aaa[0].SolutionStepData().PrintData(std::cout);
		//std::cout << std::endl;
		//std::cout << "rDeltaTemperature = " << rDeltaTemperature << std::endl;

		return rDeltaTemperature;
	}

	//******************************* COMPUTE STRAIN TEMPERATURE  ************************
	//************************************************************************************


	Vector &  J2ConstitutiveLaw3D::CalculateStrainTemperature(Parameters& rValues,
		const Properties& rMaterialProperties,
		double& rDeltaTemperature,
		Vector & rStrainTemp)
	{
		double DeltaT = CalculateDomainTemperature(rValues, rDeltaTemperature);
		double strain_size = this->GetStrainSize();

		double alpha = rMaterialProperties.Has(COEFFICIENT_THERMAL_EXPANSION) ? rMaterialProperties[COEFFICIENT_THERMAL_EXPANSION] : 0.0;

		// calculate Emod = -(- alpha * Delta_Temp)
		Vector Emod(strain_size, 0.0);
		Emod[0] = alpha * DeltaT;
		Emod[1] = alpha * DeltaT;
		Emod[2] = alpha * DeltaT;
		Emod[3] = 0.0;
		Emod[4] = 0.0;
		Emod[5] = 0.0;

		noalias(rStrainTemp) = Emod;

		return rStrainTemp;
	}

	void J2ConstitutiveLaw3D::GetLawFeatures(Features& rFeatures)
	{
		rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( ISOTROPIC );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
		rFeatures.mStrainSize = GetStrainSize();
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    }

    int J2ConstitutiveLaw3D::Check(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

		if(!rMaterialProperties.Has(YOUNG_MODULUS)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing YOUNG_MODULUS", "");
		}
		if(!rMaterialProperties.Has(POISSON_RATIO)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing POISSON_RATIO", "");
		}
		if(!rMaterialProperties.Has(DENSITY)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing DENSITY", "");
		}
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing ISOTROPIC_HARDENING", "");
		}
		if(!rMaterialProperties.Has(KINEMATIC_HARDENING)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing KINEMATIC_HARDENING", "");
		}
		if(!rMaterialProperties.Has(ISOTROPIC_HARDENING_EXPONENT)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing ISOTROPIC_HARDENING_EXPONENT", "");
		}
		if(!rMaterialProperties.Has(YIELD_STRESS)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing YIELD_STRESS", "");
		}
		if(!rMaterialProperties.Has(YIELD_STRESS_INFINITY)) {
			KRATOS_THROW_ERROR( std::logic_error, "J2ConstitutiveLaw3D - missing YIELD_STRESS_INFINITY", "");
		}
		return 0;
        KRATOS_CATCH("");
    }

} /* namespace Kratos.*/
