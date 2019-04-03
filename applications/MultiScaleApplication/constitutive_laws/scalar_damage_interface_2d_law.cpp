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

#include "scalar_damage_interface_2d_law.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"

//#define USE_AS_BRICK_INTERFACE

namespace Kratos
{

    ScalarDamageInterface2DLaw::ScalarDamageInterface2DLaw()
        : ConstitutiveLaw()
		, mInitialized(false)
		, m_initial_strain()
#ifdef INTERF_DAM_2D_IMPLEX
		, mK1_converged_old(0.0)
		, mK2_converged_old(0.0)
		, m_strain()
		, m_dTime_n(0.0)
		, m_dTime_n_converged(0.0)
#endif // INTERF_DAM_2D_IMPLEX
	{
    }

    ConstitutiveLaw::Pointer ScalarDamageInterface2DLaw::Clone() const
    {
        return ConstitutiveLaw::Pointer( new ScalarDamageInterface2DLaw() );
    }

    ScalarDamageInterface2DLaw::SizeType ScalarDamageInterface2DLaw::WorkingSpaceDimension()
    {
        return 2;
    }

    ScalarDamageInterface2DLaw::SizeType ScalarDamageInterface2DLaw::GetStrainSize()
    {
        return 2;
    }

    bool ScalarDamageInterface2DLaw::Has(const Variable<double>& rThisVariable)
    {
		if(rThisVariable == DAMAGE_T)
			return true;
		if(rThisVariable == DAMAGE_C)
			return true;
		if(rThisVariable == YIELD_FUNCTION_VALUE)
			return true;
        return false;
    }

    bool ScalarDamageInterface2DLaw::Has(const Variable<Vector>& rThisVariable)
    {
		if(rThisVariable == YIELD_SURFACE_DATA_2D_X || rThisVariable == YIELD_SURFACE_DATA_2D_Y)
			return true;
		if (rThisVariable == INITIAL_STRAIN)
			return true;
		return false;
    }

    bool ScalarDamageInterface2DLaw::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    bool ScalarDamageInterface2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }

    bool ScalarDamageInterface2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    double& ScalarDamageInterface2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
		rValue = 0.0;
		if(rThisVariable == DAMAGE_T)
			rValue = mD1;
		if(rThisVariable == DAMAGE_C)
			rValue = mD2;
		if(rThisVariable == YIELD_FUNCTION_VALUE)
			rValue = mYieldValue;
        return rValue;
    }

    Vector& ScalarDamageInterface2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
		if (rThisVariable == INITIAL_STRAIN) {
			if (rValue.size() != m_initial_strain.size())
				rValue.resize(m_initial_strain.size());
			noalias(rValue) = m_initial_strain;
		}
		if(rThisVariable == YIELD_SURFACE_DATA_2D_X || rThisVariable == YIELD_SURFACE_DATA_2D_Y)
		{
			int ijob = rThisVariable == YIELD_SURFACE_DATA_2D_X ? 1 : 2;
			double Ft = 1.0;
			double C0 = 2.0;
			double Fs = std::tan(Globals::Pi / 180 * 30.0);

			SizeType nn = 10;
			if(rValue.size() != nn)
				rValue.resize(nn, false);

			double Ft_d = (1.0 - mD1)*Ft;
			double C0_d = (1.0 - mD2)*C0;

			double x = Ft_d;
			double x_incr = (4.0*Ft+Ft_d) / (double)(nn-1);

			rValue(0) = ijob == 1 ? x : 0.0;

			for(SizeType i = 1; i < nn; i++)
			{
				double y = -x*Fs+C0_d;
				rValue(i) = ijob == 1 ? x : y;
				x -= x_incr;
			}
		}
        return rValue;
    }

    Matrix& ScalarDamageInterface2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
        return rValue;
    }

    array_1d<double, 3 > & ScalarDamageInterface2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    array_1d<double, 6 > & ScalarDamageInterface2DLaw::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
    {
        return rValue;
    }

    void ScalarDamageInterface2DLaw::SetValue(const Variable<double>& rVariable,
                                              const double& rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
		if(rVariable == DAMAGE_T)
			mD1 = rValue;
		if(rVariable == DAMAGE_C)
			mD2 = rValue;
    }

    void ScalarDamageInterface2DLaw::SetValue(const Variable<Vector >& rVariable,
                                              const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
		if (rVariable == INITIAL_STRAIN) {
			if (rValue.size() == m_initial_strain.size())
				noalias(m_initial_strain) = rValue;
		}
    }

    void ScalarDamageInterface2DLaw::SetValue(const Variable<Matrix >& rVariable,
                                              const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                              const array_1d<double, 3 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface2DLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                              const array_1d<double, 6 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
    }

    bool ScalarDamageInterface2DLaw::ValidateInput(const Properties& rMaterialProperties)
    {
		if( !rMaterialProperties.Has(NORMAL_STIFFNESS) ) return false;
		if( !rMaterialProperties.Has(TANGENTIAL_STIFFNESS) ) return false;
		if( !rMaterialProperties.Has(INTERFACE_TENSILE_LAW_S0) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_I) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_II) ) return false;
		if( !rMaterialProperties.Has(INITIAL_COHESION) ) return false;
		if( !rMaterialProperties.Has(INITIAL_FRICTION_ANGLE) ) return false;
        return true;
    }

    ScalarDamageInterface2DLaw::StrainMeasure ScalarDamageInterface2DLaw::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }

    ScalarDamageInterface2DLaw::StressMeasure ScalarDamageInterface2DLaw::GetStressMeasure()
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    bool ScalarDamageInterface2DLaw::IsIncremental()
    {
        return false;
    }

    void ScalarDamageInterface2DLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const Vector& rShapeFunctionsValues)
    {
		if(!mInitialized)
		{
			mK1 = 0.0;
			mK1_converged = 0.0;
			mK2 = 0.0;
			mK2_converged = 0.0;
			mD1 = 0.0;
			mD2 = 0.0;
			mD2_bar = 0.0;
			mD2_bar_converged = 0.0;
			mYieldValue = 0.0;
			m_initial_strain = ZeroVector(this->GetStrainSize());
			mInitialized = true;
#ifdef INTERF_DAM_2D_IMPLEX
			mK1_converged_old = mK1;
			mK2_converged_old = mK2;
			m_strain = ZeroVector(this->GetStrainSize());
			m_dTime_n = 0.0;
			m_dTime_n_converged = 0.0;
#endif // INTERF_DAM_2D_IMPLEX
		}
    }

    void ScalarDamageInterface2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
                                                            const GeometryType& rElementGeometry,
                                                            const Vector& rShapeFunctionsValues,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
#ifdef INTERF_DAM_2D_IMPLEX

		// implicit step
		// create dummy material parameters
		Vector dummy_stress(this->GetStrainSize());
		Matrix dummy_tangent(this->GetStrainSize(), this->GetStrainSize());
		ConstitutiveLaw::Parameters parameters(rElementGeometry, rMaterialProperties, rCurrentProcessInfo);
		parameters.SetStrainVector( m_strain );
		parameters.SetStressVector( dummy_stress );
		parameters.SetConstitutiveMatrix( dummy_tangent );
		Flags& options = parameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
		options.Set(ConstitutiveLaw::INITIAL_CONFIGURATION);
		double detF = 1.0;
		double detF0 = 1.0;
		Matrix F(IdentityMatrix(2,2));
		Matrix F0(IdentityMatrix(2,2));
		parameters.SetDeterminantF(detF);
		parameters.SetDeterminantF0(detF0);
		parameters.SetDeformationGradientF(F);
		parameters.SetDeformationGradientF0(F0);
		// initialize calculation data
		CalculationData data;
		InitializeCalculationData(rMaterialProperties, rElementGeometry, m_strain, rCurrentProcessInfo, data);
		CalculateElasticStressVector( data, m_strain );
		// calculate internal variables implicitly
		CalculateEquivalentMeasure( data );
		UpdateDamage( data );
		mD1 = data.D1;
		mD2 = data.D2;

		// move from n to n-1
		mK1_converged_old  = mK1_converged;
		mK2_converged_old  = mK2_converged;
		m_dTime_n_converged = m_dTime_n;

#endif // INTERF_DAM_2D_IMPLEX

		// save converged values
		mK1_converged = mK1;
		mK2_converged = mK2;
		mD2_bar_converged = mD2_bar;
    }

    void ScalarDamageInterface2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                                  const GeometryType& rElementGeometry,
                                                                  const Vector& rShapeFunctionsValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                                const GeometryType& rElementGeometry,
                                                                const Vector& rShapeFunctionsValues,
                                                                const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
        const Properties& props = rValues.GetMaterialProperties();
		const Vector& strainVector = rValues.GetStrainVector();
		Vector& stressVector = rValues.GetStressVector();
		Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
		Flags& Options = rValues.GetOptions();
		bool compute_constitutive_tensor = Options.Is(COMPUTE_CONSTITUTIVE_TENSOR);
		bool compute_stress = Options.Is(COMPUTE_STRESS) || compute_constitutive_tensor;

#ifdef INTERF_DAM_2D_IMPLEX
		this->m_strain = rValues.GetStrainVector();
#endif // INTERF_DAM_2D_IMPLEX

		SizeType size = GetStrainSize();
		if(compute_stress)
			if(stressVector.size() != size)
				stressVector.resize(size, false);
		if(compute_constitutive_tensor)
			if(constitutiveMatrix.size1() != size || constitutiveMatrix.size2() != size)
				constitutiveMatrix.resize(size, size, false);

		CalculationData data;
		InitializeCalculationData( props, rValues.GetElementGeometry(), strainVector, rValues.GetProcessInfo(), data );
		CalculateElasticStressVector( data, strainVector );

#ifdef INTERF_DAM_2D_IMPLEX

		double time_factor = 0.0;
		if(m_dTime_n_converged>0.0) time_factor = data.dTime/m_dTime_n_converged;
		m_dTime_n = data.dTime;
		mK1 = mK1_converged + time_factor * (mK1_converged-mK1_converged_old);
		mK2 = mK2_converged + time_factor * (mK2_converged-mK2_converged_old);
		if(mK1 > 0.0)
		{
			data.D1 = 1.0 - data.Ft/(mK1+data.Ft) * std::exp( -data.Ft/(data.GI*data.Kn) * mK1 );
			data.D1 = std::max( std::min( data.D1, 1.0 ), 0.0 );
		}
		if(mK2 > 0.0)
		{
			data.D2 = 1.0 - data.C0/(mK2+data.C0) * std::exp( -data.C0/data.GII/data.Kt * mK2 );
			data.D2 = std::max( std::min( data.D2, 1.0 ), 0.0 );
		}
#ifdef USE_AS_BRICK_INTERFACE
		data.D2 = 0.0;
		if(data.D1 > 0.99)
			data.D2 = 1.0;
#endif // USE_AS_BRICK_INTERFACE
#else

		CalculateEquivalentMeasure( data );
		UpdateDamage( data );

#endif // INTERF_DAM_2D_IMPLEX

		mD1 = data.D1;
		mD2 = data.D2;

		if( compute_stress )
			CalculateStress( data, stressVector );

		//**********************************************
		double sig_n = stressVector(1);
		double sig_t = std::abs(stressVector(0));
		double C0_d = (1.0 - mD2)*data.C0;
		mYieldValue = sig_n*data.Fs + sig_t - C0_d;
		//**********************************************

		if( compute_constitutive_tensor )
		{
			if(data.ForceSecant) {
				constitutiveMatrix.clear();
				constitutiveMatrix(0,0) = data.Kt*(1.0-mD2);
				constitutiveMatrix(1,1) = data.Kn;
				if(stressVector(1) > 0.0) {
					constitutiveMatrix(1,1) *= (1.0-mD1);
				}
			}
			else {
				CalculateConstitutiveMatrix( data, strainVector, stressVector, constitutiveMatrix );
			}
		}
    }

    void ScalarDamageInterface2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {

    }

    void ScalarDamageInterface2DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                   const GeometryType& rElementGeometry,
                                                   const Vector& rShapeFunctionsValues)
    {
        mInitialized = false;
		mK1 = 0.0;
		mK1_converged = 0.0;
		mK2 = 0.0;
		mK2_converged = 0.0;
		mD1 = 0.0;
		mD2 = 0.0;
		mYieldValue = 0.0;
    }

    void ScalarDamageInterface2DLaw::GetLawFeatures(Features& rFeatures)
    {
        //Set the type of law
		rFeatures.mOptions.Set( PLANE_STRESS_LAW ); // TODO: INTERFACE 2D LAW
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mOptions.Set( ISOTROPIC );

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
    }

    int ScalarDamageInterface2DLaw::Check(const Properties& rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

		if( !rMaterialProperties.Has(NORMAL_STIFFNESS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: NORMAL_STIFFNESS", "");

		if( !rMaterialProperties.Has(TANGENTIAL_STIFFNESS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TANGENTIAL_STIFFNESS", "");

		if( !rMaterialProperties.Has(INTERFACE_TENSILE_LAW_S0) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INTERFACE_TENSILE_LAW_S0", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_I) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_MODE_I", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_II) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_MODE_II", "");

		if( !rMaterialProperties.Has(INITIAL_COHESION) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INITIAL_COHESION", "");

		if( !rMaterialProperties.Has(INITIAL_FRICTION_ANGLE) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INITIAL_FRICTION_ANGLE", "");

        return 0;

        KRATOS_CATCH("");
    }

	void ScalarDamageInterface2DLaw::InitializeCalculationData(const Properties& props,
		                                                       const GeometryType& geom,
															   const Vector& strainVector,
															   const ProcessInfo& pinfo,
															   CalculationData& data)
	{
		data.Kn  = props[NORMAL_STIFFNESS];
		data.Kt  = props[TANGENTIAL_STIFFNESS];

		data.Kn_compression_multiplier = 1.0;
		if(props.Has(NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER)) {
			data.Kn_compression_multiplier = props[NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER];
			if(data.Kn_compression_multiplier < 1.0)
				data.Kn_compression_multiplier = 1.0;
			if(strainVector(1) <= 0.0)
				data.Kn *= data.Kn_compression_multiplier;
		}

		data.Ft  = props[INTERFACE_TENSILE_LAW_S0];

		data.C0 = props[INITIAL_COHESION];
		data.Fs = std::tan( props[INITIAL_FRICTION_ANGLE] * Globals::Pi / 180.0 );

		data.GI  = props[FRACTURE_ENERGY_MODE_I];
		data.GII = props[FRACTURE_ENERGY_MODE_II];

		// check
		if(data.C0 <= 0.0) {
			data.C0 = 0.0;
			data.Ft = 0.0;
			data.MRatio = 1.0;
		}
		else {
			if(data.Ft <= 0.0) {
				data.Ft = 0.0;
				data.MRatio = 1.0;
			}
			else {
				double Ft_max = data.C0 / data.Fs;
				if(data.Ft > Ft_max) {
					data.Ft = Ft_max;
				}
				data.MRatio = data.Ft / data.C0 * data.GII / data.GI;
			}
		}

		data.D1 = 0.0;
		data.D2 = 0.0;
		data.K1 = 0.0;
		data.K2 = 0.0;

		data.ElasticStressVector.clear();

		data.ForceSecant = false;
		if(props.Has(DAMAGE_SECANT_MATRIX))
			data.ForceSecant = props[DAMAGE_SECANT_MATRIX] != 0;

		data.lch = 1.0; // * m_lch_multiplier;
		data.eta = 0.0;
		if(props.Has(VISCOSITY)) data.eta = props[VISCOSITY];
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

	void ScalarDamageInterface2DLaw::CalculateElasticStressVector(CalculationData& data,
		                                                          const Vector& strainVector)
	{
		data.ElasticStressVector(0) = data.Kt * strainVector(0);
		data.ElasticStressVector(1) = data.Kn * strainVector(1);
	}

	void ScalarDamageInterface2DLaw::CalculateEquivalentMeasure(CalculationData& data)
	{
		//double sigma_t = std::abs(data.ElasticStressVector(0));
		//double sigma_n = data.ElasticStressVector(1);

		//double lambda1 = sigma_n - data.Ft;
		//double lambda2 = sigma_n * data.Fs + sigma_t - data.C0;

		//data.K1 = 0.0;
		//data.K2 = 0.0;

		//if(lambda1 > 0.0)
		//{
		//	if(lambda2 > 0.0)
		//	{
		//		// both active
		//		data.K1 = std::sqrt( lambda1*lambda1 + (lambda2/data.MRatio)*(lambda2/data.MRatio) );
		//		data.K2 = std::sqrt( lambda2*lambda2 + (lambda1*data.MRatio)*(lambda1*data.MRatio) );
		//	}
		//	else
		//	{
		//		// only tension cut-off
		//		data.K1 = lambda1;
		//		data.K2 = lambda1 * data.MRatio;
		//	}
		//}
		//else
		//{
		//	if(lambda2 > 0.0)
		//	{
		//		// only coulomb friction
		//		data.K1 = lambda2 / data.MRatio;
		//		data.K2 = lambda2;
		//	}
		//}

		double sigma_t = std::abs(data.ElasticStressVector(0));
		double sigma_n = data.ElasticStressVector(1);

		data.K1 = 0.0;
		data.K2 = 0.0;

		if(sigma_n <= 0.0)
		{
			double lambda2 = sigma_n * data.Fs + sigma_t - data.C0;
			if(lambda2 > 0.0)
			{
				// only coulomb friction
				//data.K1 = lambda2 / data.MRatio;
				data.K2 = lambda2;
			}
		}
		else
		{
			/*double lambda1 = sigma_n - data.Ft;
			double lambda2 = sigma_n * data.Fs + sigma_t - data.C0;*/
			double lambda = std::sqrt(std::pow(sigma_n,2) + std::pow(data.Ft/data.C0*sigma_t,2)) - data.Ft;
			if(lambda>0.0)
			{
				data.K1 = lambda;
				data.K2 = lambda;
			}
			//if(lambda1 > 0.0)
			//{
			//	if(lambda2 > 0.0)
			//	{
			//		// both active
			//		data.K1 = std::sqrt( lambda1*lambda1 + (lambda2/data.MRatio)*(lambda2/data.MRatio) );
			//		data.K2 = std::sqrt( lambda2*lambda2 + (lambda1*data.MRatio)*(lambda1*data.MRatio) );
			//	}
			//	else
			//	{
			//		// only tension cut-off
			//		data.K1 = lambda1;
			//		data.K2 = lambda1 * data.MRatio;
			//	}
			//}
			//else
			//{
			//	if(lambda2 > 0.0)
			//	{
			//		// only coulomb friction
			//		data.K1 = lambda2 / data.MRatio;
			//		data.K2 = lambda2;
			//	}
			//}
		}
	}

	void ScalarDamageInterface2DLaw::UpdateDamageIndicator(CalculationData& data)
	{
		mK1 = std::max( data.K1, mK1_converged );
		mK2 = std::max( data.K2, mK2_converged );
	}

	void ScalarDamageInterface2DLaw::CalculateDamage(CalculationData& data)
	{
		if(mK1 > 0.0)
		{
			data.D1 = 1.0 - data.Ft/(mK1+data.Ft) * std::exp( -data.Ft/(data.GI*data.Kn) * mK1 );
			data.D1 = std::max( std::min( data.D1, 1.0 ), 0.0 );
		}
		if(mK2 > 0.0)
		{
			data.D2 = 1.0 - data.C0/(mK2+data.C0) * std::exp( -data.C0/data.GII/data.Kt * mK2 );
			data.D2 = std::max( std::min( data.D2, 1.0 ), 0.0 );
		}
	}

	void ScalarDamageInterface2DLaw::UpdateDamage(CalculationData& data)
	{
		bool update_equ_shear_damage = false;

		mK1 = mK1_converged;
		mK2 = mK2_converged;
		mD2_bar = mD2_bar_converged;

		if(data.K1 > mK1)
		{
			//mK1 = data.K1;
			mK1 = data.rate_coeff_1*mK1 + data.rate_coeff_2*data.K1;
		}
		if(data.K2 > mK2)
		{
			//mK2 = data.K2;
			mK2 = data.rate_coeff_1*mK2 + data.rate_coeff_2*data.K2;
			update_equ_shear_damage = true;
		}



		if(mK1 > 0.0)
		{
			//data.D1 = 1.0 - data.Ft/(mK1+data.Ft) * std::exp( -data.Ft/(data.GI*data.Kn) * mK1 );
			double lt  = 2.0*data.Kn*data.GI/data.Ft/data.Ft;
			double Hs  = data.lch/(lt-data.lch);
			double A   = 2.0*Hs;
			data.D1 = 1.0-data.Ft/(mK1+data.Ft)*std::exp(A*(1.0-(mK1+data.Ft)/data.Ft));
			data.D1 = std::max( std::min( data.D1, 1.0 ), 0.0 );
		}
		if(mK2 > 0.0)
		{
			//data.D2 = 1.0 - data.C0/(mK2+data.C0) * std::exp( -data.C0/data.GII/data.Kt * mK2 );
			double lt  = 2.0*data.Kt*data.GII/data.C0/data.C0;
			double Hs  = data.lch/(lt-data.lch);
			double A   = 2.0*Hs;
			data.D2 = 1.0-data.C0/(mK2+data.C0)*std::exp(A*(1.0-(mK2+data.C0)/data.C0));
			data.D2 = std::max( std::min( data.D2, 1.0 ), 0.0 );
		}

#ifdef USE_AS_BRICK_INTERFACE
		data.D2 = 0.0;
		if(data.D1 > 0.25)
			data.D2 = 1.0;
#endif // USE_AS_BRICK_INTERFACE

		if(update_equ_shear_damage)
		{
			double sigma_n = data.ElasticStressVector(1);
			double sigma_t = data.ElasticStressVector(0);
			if(sigma_t != 0.0)
			{
				if(sigma_n < 0.0) {
					double C0_d = (1.0 - data.D2)*data.C0;
					double tau = std::max(C0_d - sigma_n*data.Fs, 0.0);

					if(sigma_t < 0.0) tau = -tau;

					if(std::abs(sigma_t) > std::abs(tau))
					{
						mD2_bar = 1.0 - tau / sigma_t;
						mD2_bar = std::max(std::min(mD2_bar,1.0),0.0);
					}
				}
				else {
					mD2_bar = data.D2;
				}
			}
		}
	}

	void ScalarDamageInterface2DLaw::CalculateStress(CalculationData& data,
		                                             Vector& stressVector)
	{
		double sigma_n = data.ElasticStressVector(1);
		double sigma_t = data.ElasticStressVector(0);

		if(sigma_n > 0.0) sigma_n *= (1.0 - data.D1);
		stressVector(1) = sigma_n;

		stressVector(0) = (1.0 - mD2_bar) * sigma_t;
		//stressVector(0) = (1.0 - data.D2) * sigma_t;
	}

	void ScalarDamageInterface2DLaw::CalculateConstitutiveMatrix(CalculationData& data,
		                                                         const Vector& strainVector,
		                                                         const Vector& stressVector,
																 Matrix& constitutiveMatrix)
	{
		// elastic case
		constitutiveMatrix(0, 0) = data.Kt;
		constitutiveMatrix(1, 1) = data.Kn;
		constitutiveMatrix(0, 1) = constitutiveMatrix(1, 0) = 0.0;

		double perturbation = 1.0E-9;
		Vector perturbedStrainVector(2);
		Vector stressPerturbation(2);

		for(int j = 0; j < 2; j++)
		{
			// save internal variables
			double save_k1 = mK1;
			double save_k2 = mK2;
			double save_d2_bar = mD2_bar;

			// FORWARD difference
			noalias( perturbedStrainVector ) = strainVector;

			//double delta_strain = perturbedStrainVector(j) < 0.0 ? -perturbation : perturbation;
			double delta_strain = perturbation;
			perturbedStrainVector(j) += delta_strain;

			CalculateElasticStressVector( data, perturbedStrainVector );
			CalculateEquivalentMeasure( data );
			UpdateDamage( data );
			CalculateStress( data, stressPerturbation );

			// fill the numerical tangent operator
			noalias( stressPerturbation ) -= stressVector;
			constitutiveMatrix(0, j) = stressPerturbation(0) / delta_strain;
			constitutiveMatrix(1, j) = stressPerturbation(1) / delta_strain;

			// restore internal variables
			mK1 = save_k1;
			mK2 = save_k2;
			mD2_bar = save_d2_bar;
		}
	}

} /* namespace Kratos.*/
