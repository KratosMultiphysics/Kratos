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
*   Date:                $Date:     19-09-2014$
*   Revision:            $Revision: 1.0$
*
* ***********************************************************/

#include "damage_iso_plane_stress_2d_law.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"
#include "custom_utilities/imperfection_utilities.h"
#include "includes/variables.h"

#define DAM_ISO_PREC 1.0E-12

#define DAM_ISO_SIGN(X) (X == 0.0 ? 0.0 : ( X > 0.0 ? 1.0 : -1.0 ))

//#define DAM_ISO_USE_SMOOTH_TENSILE_LAW
#define DAM_ISO_STLAW_A1 0.9

namespace Kratos
{

	DamageIsoPlaneStress2DLaw::DamageIsoPlaneStress2DLaw()
		: ConstitutiveLaw()
		, m_initialized(false)
		, m_r(0.0)
		, m_r_converged(0.0)
		, m_damage(0.0)
		, m_lch(0.0)
		, m_lch_multiplier(1.0)
		, m_initial_strain()
#ifdef DAM_ISO_2D_IMPLEX
		, m_r_converged_old(0.0)
		, m_strain()
		, m_dTime_n(0.0)
		, m_dTime_n_converged(0.0)
		, m_r_impl_temp(0.0)
#endif // DAM_ISO_2D_IMPLEX
		, m_suggested_time_step(0.0)
		, m_error_code(0.0)
	{
	}

	ConstitutiveLaw::Pointer DamageIsoPlaneStress2DLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer( new DamageIsoPlaneStress2DLaw() );
	}

	DamageIsoPlaneStress2DLaw::SizeType DamageIsoPlaneStress2DLaw::WorkingSpaceDimension()
	{
		return 2;
	}

	DamageIsoPlaneStress2DLaw::SizeType DamageIsoPlaneStress2DLaw::GetStrainSize()
	{
		return 3;
	}

	bool DamageIsoPlaneStress2DLaw::Has(const Variable<double>& rThisVariable)
	{
		if(rThisVariable == DAMAGE_T)
			return true;
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			return true;
		if(rThisVariable == SUGGESTED_TIME_STEP)
			return true;
		return false;
	}

	bool DamageIsoPlaneStress2DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		if(rThisVariable == INITIAL_STRAIN)
			return true;
		return false;
	}

	bool DamageIsoPlaneStress2DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool DamageIsoPlaneStress2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	bool DamageIsoPlaneStress2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

	double& DamageIsoPlaneStress2DLaw::GetValue(
		const Variable<double>& rThisVariable,
		double& rValue)
	{
		rValue = 0.0;
		if(rThisVariable == DAMAGE_T)
			rValue = m_damage;
		if(rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = m_error_code;
		if(rThisVariable == SUGGESTED_TIME_STEP)
			rValue = m_suggested_time_step;
		return rValue;
	}

	Vector& DamageIsoPlaneStress2DLaw::GetValue(
		const Variable<Vector>& rThisVariable,
		Vector& rValue)
	{
		if(rThisVariable == INITIAL_STRAIN) {
			if(rValue.size() != m_initial_strain.size())
				rValue.resize(m_initial_strain.size());
			noalias(rValue) = m_initial_strain;
		}
		return rValue;
	}

	Matrix& DamageIsoPlaneStress2DLaw::GetValue(
		const Variable<Matrix>& rThisVariable,
		Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & DamageIsoPlaneStress2DLaw::GetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & DamageIsoPlaneStress2DLaw::GetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void DamageIsoPlaneStress2DLaw::SetValue(
		const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable == DAMAGE_T)
			m_damage = rValue;
		else if(rVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			m_lch_multiplier = rValue;
	}

	void DamageIsoPlaneStress2DLaw::SetValue(
		const Variable<Vector >& rVariable,
		const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable == INITIAL_STRAIN) {
			if(rValue.size() == m_initial_strain.size())
				noalias(m_initial_strain) = rValue;
		}
	}

	void DamageIsoPlaneStress2DLaw::SetValue(
		const Variable<Matrix >& rVariable,
		const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageIsoPlaneStress2DLaw::SetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageIsoPlaneStress2DLaw::SetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool DamageIsoPlaneStress2DLaw::ValidateInput(const Properties& rMaterialProperties)
	{
		if( !rMaterialProperties.Has(YOUNG_MODULUS) ) return false;
		if( !rMaterialProperties.Has(POISSON_RATIO) ) return false;
		if( !rMaterialProperties.Has(DAMAGE_STRESS_T_0) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_T) ) return false;
		if( !rMaterialProperties.Has(VISCOSITY) ) return false;
		return true;
	}

	DamageIsoPlaneStress2DLaw::StrainMeasure DamageIsoPlaneStress2DLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	DamageIsoPlaneStress2DLaw::StressMeasure DamageIsoPlaneStress2DLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool DamageIsoPlaneStress2DLaw::IsIncremental()
	{
		return false;
	}

	void DamageIsoPlaneStress2DLaw::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if(!m_initialized)
		{

			// random imperfections
			double impf = ImperfectionUtilties::CalculateRandomImperfectionScaleFactor(
				rElementGeometry,rShapeFunctionsValues);

#ifndef DAM_ISO_USE_SMOOTH_TENSILE_LAW
			m_r              = rMaterialProperties[DAMAGE_STRESS_T_0]*impf;
#else
			m_r              = rMaterialProperties[DAMAGE_STRESS_T_0]*impf*DAM_ISO_STLAW_A1;
#endif // !DAM_ISO_USE_SMOOTH_TENSILE_LAW

			m_r_converged    = m_r;
			m_damage         = 0.0;
			m_lch            = rElementGeometry.Length();
			m_lch_multiplier = 1.0;
			m_initial_strain = ZeroVector(this->GetStrainSize());
			m_temp_strain = ZeroVector(this->GetStrainSize());
			m_initialized    = true;

#ifdef DAM_ISO_2D_IMPLEX
			m_r_converged_old = m_r;
			m_strain = ZeroVector(this->GetStrainSize());
			m_dTime_n = 0.0;
			m_dTime_n_converged = 0.0;
#endif // DAM_ISO_2D_IMPLEX

			m_error_code = 0.0;

		}
	}

	void DamageIsoPlaneStress2DLaw::InitializeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageIsoPlaneStress2DLaw::FinalizeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
#ifdef DAM_ISO_2D_IMPLEX

		m_r = m_r_impl_temp;

		// move from n to n-1
		m_r_converged_old  = m_r_converged;
		m_dTime_n_converged = m_dTime_n;

#endif // DAM_ISO_2D_IMPLEX

		m_r_converged = m_r;
	}

	void DamageIsoPlaneStress2DLaw::InitializeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageIsoPlaneStress2DLaw::FinalizeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void DamageIsoPlaneStress2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageIsoPlaneStress2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageIsoPlaneStress2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void DamageIsoPlaneStress2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
	{
		const ProcessInfo&  pinfo = rValues.GetProcessInfo();
		const GeometryType& geom  = rValues.GetElementGeometry();
		const Properties&   props = rValues.GetMaterialProperties();

#ifdef DAM_ISO_2D_IMPLEX
		this->m_strain = rValues.GetStrainVector();
#endif // DAM_ISO_2D_IMPLEX

		Vector& strain_vector = rValues.GetStrainVector();
		Vector&       stress_vector       = rValues.GetStressVector();

		CalculationData data;
		this->InitializeCalculationData(props, geom, rValues.GetShapeFunctionsValues(), pinfo, data);

		//-1.- Initial & Thermal Strain
		// Subtract Initial_Strain
		noalias(strain_vector) -= m_initial_strain;

		// Subtract Temperature_Strain
		double DeltaTemp = 0.0;
		Vector temp_strain(strain_vector.size(), 0.0);
		CalculateStrainTemperature(rValues, props, DeltaTemp, temp_strain);
		noalias(strain_vector) -= temp_strain;

		this->CalculateMaterialResponseInternal(strain_vector, stress_vector, data);

#ifdef DAM_ISO_2D_IMPLEX
		if(rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR))
		{
			size_t n = GetStrainSize();
			Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
			if(constitutive_matrix.size1() != n || constitutive_matrix.size2() != n)
				constitutive_matrix.resize(n, n);
			noalias(constitutive_matrix) = (1.0-m_damage)*data.C0;
		}
		return;
#endif // DAM_ISO_2D_IMPLEX

		if (m_error_code != 0.0) return;

		if(rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR))
		{
			size_t n = GetStrainSize();

			// prepare constitutive matrix
			Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
			if(constitutive_matrix.size1() != n || constitutive_matrix.size2() != n)
				constitutive_matrix.resize(n, n);

			if(props.Has(DAMAGE_SECANT_MATRIX))
			{
				if(props[DAMAGE_SECANT_MATRIX] > 0.0)
				{
					noalias(constitutive_matrix) = (1.0-m_damage)*data.C0;
					return;
				}
			}

			// save internal variables
			double save_r = m_r;
			double save_d = m_damage;

			// perturbation parameter
			double h = 1.0E-10;

			// perturbed vectors
			Vector strain_bar(n);
			Vector S1(n);
			Vector S2(n);

			// apply perturbation to each strain component...
			for(size_t j = 0; j < n; j++)
			{
				noalias(strain_bar) = strain_vector;

				strain_bar(j) = strain_vector(j) - h;
				this->CalculateMaterialResponseInternal(strain_bar, S1, data);

				strain_bar(j) = strain_vector(j) + h;
				this->CalculateMaterialResponseInternal(strain_bar, S2, data);

				for(size_t i = 0; i < n; i++)
					constitutive_matrix(i,j) = (S2(i) - S1(i))/(2.0*h);

				/*strain_bar(j) = strain_vector(j) + h;
				this->CalculateMaterialResponseInternal(strain_bar, S2, data);

				for(size_t i = 0; i < n; i++)
					constitutive_matrix(i,j) = (S2(i) - stress_vector(i))/h;*/
			}

			// restore internal variables
			m_r	     = save_r;
			m_damage = save_d;
		}
	}

	void DamageIsoPlaneStress2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageIsoPlaneStress2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageIsoPlaneStress2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void DamageIsoPlaneStress2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{

	}

	void DamageIsoPlaneStress2DLaw::ResetMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		m_error_code = 0.0;
		m_r = 0.0;
		m_r_converged = 0.0;
		m_damage = 0.0;
		m_lch = 0.0;
		m_lch_multiplier = 1.0;
		m_initialized = false;
	}

	void DamageIsoPlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set( PLANE_STRESS_LAW );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mOptions.Set( ISOTROPIC );

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}

	int DamageIsoPlaneStress2DLaw::Check(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if( !rMaterialProperties.Has(YOUNG_MODULUS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS", "");

		if( !rMaterialProperties.Has(POISSON_RATIO) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: POISSON_RATIO", "");

		if( !rMaterialProperties.Has(DAMAGE_STRESS_T_0) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_STRESS_T_0", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_T) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_T", "");

		if( !rMaterialProperties.Has(VISCOSITY) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: VISCOSITY", "");

		return 0;

		KRATOS_CATCH("");
	}

	void DamageIsoPlaneStress2DLaw::CalculateMaterialResponse(const Vector& StrainVector,
			const Matrix& DeformationGradient,
			Vector& StressVector,
			Matrix& AlgorithmicTangent,
			const ProcessInfo& rCurrentProcessInfo,
			const Properties& rMaterialProperties,
			const GeometryType& rElementGeometry,
			const Vector& rShapeFunctionsValues,
			bool CalculateStresses,
			int CalculateTangent,
			bool SaveInternalVariables)
	{
		ConstitutiveLaw::Parameters parameters(rElementGeometry, rMaterialProperties, rCurrentProcessInfo);
		Vector E(StrainVector);
		parameters.SetStrainVector( E );
		parameters.SetStressVector( StressVector );
		parameters.SetConstitutiveMatrix( AlgorithmicTangent );
		Flags& options = parameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, CalculateStresses);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, CalculateTangent);
		double detF = 1.0;
		Matrix F(IdentityMatrix(2,2));
		parameters.SetDeterminantF(detF);
		parameters.SetDeformationGradientF(F);
		parameters.SetShapeFunctionsValues(rShapeFunctionsValues);
		this->CalculateMaterialResponseCauchy(parameters);
	}

	void DamageIsoPlaneStress2DLaw::InitializeCalculationData(const Properties& props,
															 const GeometryType& geom,
															 const Vector& N,
															 const ProcessInfo& pinfo,
															 CalculationData& data)
	{
		// random imperfections
		double impf = ImperfectionUtilties::CalculateRandomImperfectionScaleFactor(geom, N);

		// elasticity
		data.E   = props[YOUNG_MODULUS];
		data.nu  = props[POISSON_RATIO];
		this->CalculateElasticityMatrix(data);

		// tension
		data.ft = props[DAMAGE_STRESS_T_0]*impf;
		data.Gt = props[FRACTURE_ENERGY_T]*impf*impf;

		// effective stress data
		data.S.resize(3,false);
		data.Si.resize(2,false);

		// misc
		data.eta = props[VISCOSITY];
		data.lch = m_lch * m_lch_multiplier;
		data.dTime = pinfo[DELTA_TIME];
		if(data.dTime > 0.0 && data.eta > 0.0)
		{
			data.rate_coeff_1 = data.eta/(data.eta+data.dTime);
			data.rate_coeff_2 = data.dTime/(data.eta+data.dTime);
		}
		else
		{
			data.rate_coeff_1 = 0.0;
			data.rate_coeff_2 = 1.0;
		}
	}

	void DamageIsoPlaneStress2DLaw::CalculateElasticityMatrix(CalculationData& data)
	{
		if(data.C0.size1() != 3 || data.C0.size2() != 3)
			data.C0.resize(3,3,false);

		double c1 = data.E / (1.0 - data.nu*data.nu);
		double c2 = c1 * data.nu;
		double c3 = c1 * (1.0 - data.nu) / 2.0;

		data.C0(0,0) = c1;		data.C0(0,1) = c2;		data.C0(0,2) = 0.0;
		data.C0(1,0) = c2;		data.C0(1,1) = c1;		data.C0(1,2) = 0.0;
		data.C0(2,0) = 0.0;		data.C0(2,1) = 0.0;		data.C0(2,2) = c3;
	}

	void DamageIsoPlaneStress2DLaw::StressDecomposition(CalculationData& data)
	{
		const Vector& S  = data.S;
		Vector&       Si = data.Si;
		if(Si.size() != 2) Si.resize(2,false);

		if(std::abs(S(2)) < DAM_ISO_PREC)
		{
			if(S(0) > S(1))
			{
				Si(0) = S(0);
				Si(1) = S(1);
			}
			else
			{
				Si(0) = S(1);
				Si(1) = S(0);
			}
		}
		else
		{
			double Tr = S(0)+S(1);
			double Dt = S(0)*S(1) - S(2)*S(2);

			Si(0) = Tr/2.0 + std::sqrt(Tr*Tr/4.0 - Dt);
			Si(1) = Tr/2.0 - std::sqrt(Tr*Tr/4.0 - Dt);

			if(std::abs(Si(0)) < DAM_ISO_PREC) Si(0) = 0.0;
			if(std::abs(Si(1)) < DAM_ISO_PREC) Si(1) = 0.0;
		}
	}

	void DamageIsoPlaneStress2DLaw::TrialEquivalentStress(CalculationData& data, double& r_trial)
	{
		r_trial = std::max(std::max(data.Si(0), data.Si(1)),0.0);
	}






	inline double b3_eval_bezier(double xi,
								 double x0, double x1, double x2,
								 double y0, double y1, double y2)
	{
		double A = x0-2.0*x1+x2;
		double B = 2.0*(x1-x0);
		double C = x0-xi;
		if(std::abs(A)<1.0E-12)
		{
			x1 = x1+1.0E-6*(x2-x0);
			A = x0-2.0*x1+x2;
			B = 2.0*(x1-x0);
			C = x0-xi;
		}
		double D = B*B-4.0*A*C;
		double t = (-B+std::sqrt(D))/(2.0*A);
		return (y0-2.0*y1+y2)*t*t+2.0*(y1-y0)*t+y0;
	}

	inline double b3_eval_area(double x1,double x2,double x3,double y1,double y2,double y3)
	{
		return x2*y1/3.0 + x3*y1/6.0 - x2*y3/3.0 + x3*y2/3.0 + x3*y3/2.0
			  -x1*(y1/2.0 + y2/3.0 + y3/6.0);
	}

	inline double b3_calc_G(double sp, double sk, double sr,
							double ep, double ej, double ek, double er, double eu)
	{
		double G1 = ep*sp/2.0;
		double G2 = b3_eval_area(ep,ej,ek,sp,sp,sk);
		double G3 = b3_eval_area(ek,er,eu,sk,sr,sr);
		return G1+G2+G3;
	}

	inline void b3_stretch(double S, double ep, double &ej, double &ek, double &er, double &eu)
	{
		ej = ej+(ej-ep)*S;
		ek = ek+(ek-ep)*S;
		er = er+(er-ep)*S;
		eu = eu+(eu-ep)*S;
	}


	void DamageIsoPlaneStress2DLaw::CalculateDamage(CalculationData& data, double r, double& d)
	{
#ifndef DAM_ISO_USE_SMOOTH_TENSILE_LAW

		if(r <= data.ft)
		{
			d = 0.0;
		}
		else
		{
			double lch = data.lch;
			double E   = data.E;
			double ft  = data.ft;
			double G   = data.Gt;
			double r0  = ft;
			double lt  = 2.0*E*G/ft/ft;

			if(lch >= lt)
			{
				std::stringstream ss;
				ss << "FRACTURE_ENERGY_T is to low:  2*E*Gt/(ft*ft) = " << lt
					<< ",   Characteristic Length = " << lch << std::endl;
				std::cout << ss.str();
				exit(-1);
			}

			double Hs  = lch/(lt-lch);
			double A   = 2.0*Hs;
			d          = 1.0-r0/r*std::exp(A*(1.0-r/r0));
		}

#else

		if(r <= data.ft*DAM_ISO_STLAW_A1)
		{
			d = 0.0;
		}
		else
		{
			double lch = data.lch;
			double E   = data.E;
			double ft  = data.ft;
			double G   = data.Gt/lch;

			// extract material parameters
			double s0 = ft*DAM_ISO_STLAW_A1;
			double sp = ft;
			double sr = ft/1000.0;
			double ep = ft/E*1.3;
			double c1 = 0.9;
			double c2 = 0.4;
			double c3 = 4.0;

			// auto-computation of other parameters
			double sk    = sr + (sp-sr)*c1;
			double e0    = s0/E;
			double alpha = 2.0*(ep-sp/E);
			double ej    = ep + alpha*c2;
			double ek    = ej + alpha*(1-c2);
			double er    = (ek-ej)/(sp-sk)*(sp-sr)+ej;
			double eu    = er*c3;

			// regularization
			double G_bar   = b3_calc_G( sp,sk,sr,ep,ej,ek,er,eu );
			double G1      = sp*ep/2.0;
			double stretch = (G-G1)/(G_bar-G1)-1.0;
			if(stretch <= -1.0)
			{
				std::stringstream ss;
				ss << "Damage Iso Error: fracture energy is too low" << std::endl;
				ss << "Input G/lch = " << G << std::endl;
				ss << "Minimum G to avoid constitutive snap-back = " << G1 << std::endl;
				std::cout << ss.str();
				exit(-1);
			}
			b3_stretch( stretch,ep,ej,ek,er,eu );

			// current abscissa
			double xi = r/E;

			// compute damage
			double s = r;
			if(xi <= ep)
				s = b3_eval_bezier(xi,e0,sp/E,ep,s0,sp,sp);
			else if(xi <= ek)
				s = b3_eval_bezier(xi,ep,ej,ek,sp,sp,sk);
			else if(xi <= eu)
				s = b3_eval_bezier(xi,ek,er,eu,sk,sr,sr);
			else
				s = sr;
			d = 1.0-s/r;
		}

#endif // !DAM_ISO_USE_SMOOTH_TENSILE_LAW

	}

	void DamageIsoPlaneStress2DLaw::CalculateMaterialResponseInternal(const Vector& strain_vector,
																	 Vector& stress_vector,
																	 CalculationData& data)
	{
		m_error_code = 0.0;
		m_suggested_time_step = data.dTime;

		size_t strain_size = this->GetStrainSize();

		if(stress_vector.size() != strain_size)
			stress_vector.resize(strain_size,false);

		// set up coefficients for the rate-dependent model

		double rate_coeff_1 = data.rate_coeff_1;
		double rate_coeff_2 = data.rate_coeff_2;

		// elastic predictor + stress decomposition

		m_r = m_r_converged;

		noalias(data.S) = prod(data.C0, strain_vector);

		if(std::abs(data.S(0)) < DAM_ISO_PREC) data.S(0) = 0.0;
		if(std::abs(data.S(1)) < DAM_ISO_PREC) data.S(1) = 0.0;
		if(std::abs(data.S(2)) < DAM_ISO_PREC) data.S(2) = 0.0;

		this->StressDecomposition(data);

		// compute the equivalent stress measure

		double r_trial;
		this->TrialEquivalentStress(data, r_trial);

		// damage update

#ifdef DAM_ISO_2D_IMPLEX

		double dtime_n = m_dTime_n; // save for next time step estimation

		double time_factor = 0.0;
		if(m_dTime_n_converged>0.0) time_factor = data.dTime/m_dTime_n_converged;
		m_dTime_n = data.dTime;
		m_r = m_r_converged + time_factor * (m_r_converged-m_r_converged_old);

		// new damage (implicit)
		double r_impl=m_r_converged;
		double d_impl=0.0;
		if(r_trial > r_impl)
			r_impl = rate_coeff_1*r_impl + rate_coeff_2*r_trial;
		m_r_impl_temp = r_impl;
		this->CalculateDamage(data, r_impl, d_impl);
		// new damage explicit
		this->CalculateDamage(data, m_r, m_damage);
		// check error
		//double max_damage_diff = 0.01; // damage goes from 0 to 1
		//double damage_diff = std::abs(m_damage-d_impl);
		//if(damage_diff > max_damage_diff) {
		//	m_error_code = -1.0;
		//	m_suggested_time_step = max_damage_diff/damage_diff*data.dTime;
		//}

		/*
		r(n+1)
		r_converged (n)
		r_converged_old (n-1)
		*/
		double r_error = m_r - r_impl;
		double tolerance = 0.5 * data.ft;
		if(r_error > tolerance)
		{
			m_error_code = -1.0;
			m_suggested_time_step = std::sqrt( tolerance*dtime_n*dtime_n*1.0/std::abs(2.0*(m_r_converged-m_r_converged_old)) );
		}

#else

		if(r_trial > m_r)
			m_r = rate_coeff_1*m_r + rate_coeff_2*r_trial;
		this->CalculateDamage(data, m_r, m_damage);

#endif // DAM_ISO_2D_IMPLEX

		// calculation of stress tensor

		noalias(stress_vector)  = (1.0 - m_damage)*data.S;
	}

	double &  DamageIsoPlaneStress2DLaw::CalculateDomainTemperature(Parameters& rValues,
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


	Vector &  DamageIsoPlaneStress2DLaw::CalculateStrainTemperature(Parameters& rValues,
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
		Emod[2] = 0.0;

		noalias(rStrainTemp) = Emod;

		return rStrainTemp;
		// END MOD STEFANO
	}


} /* namespace Kratos.*/
