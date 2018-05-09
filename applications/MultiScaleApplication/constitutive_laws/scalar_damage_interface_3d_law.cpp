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

#include "scalar_damage_interface_3d_law.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"

namespace Kratos
{

    ScalarDamageInterface3DLaw::ScalarDamageInterface3DLaw()
        : ConstitutiveLaw()
		, mInitialized(false)
		, m_initial_strain()
    {
    }

    ConstitutiveLaw::Pointer ScalarDamageInterface3DLaw::Clone() const
    {
        return ConstitutiveLaw::Pointer( new ScalarDamageInterface3DLaw() );
    }

    ScalarDamageInterface3DLaw::SizeType ScalarDamageInterface3DLaw::WorkingSpaceDimension()
    {
        return 3;
    }

    ScalarDamageInterface3DLaw::SizeType ScalarDamageInterface3DLaw::GetStrainSize()
    {
        return 3;
    }

    bool ScalarDamageInterface3DLaw::Has(const Variable<double>& rThisVariable)
    {
		if(rThisVariable == DAMAGE_T)
			return true;
		if(rThisVariable == DAMAGE_C)
			return true;
		if(rThisVariable == YIELD_FUNCTION_VALUE)
			return true;
        return false;
    }

    bool ScalarDamageInterface3DLaw::Has(const Variable<Vector>& rThisVariable)
    {
		if(rThisVariable == YIELD_SURFACE_DATA_3D_X || rThisVariable == YIELD_SURFACE_DATA_3D_Y || rThisVariable == YIELD_SURFACE_DATA_3D_Z)
			return true;
		if (rThisVariable == INITIAL_STRAIN)
			return true;
        return false;
    }

    bool ScalarDamageInterface3DLaw::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    bool ScalarDamageInterface3DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }

    bool ScalarDamageInterface3DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    double& ScalarDamageInterface3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
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

    Vector& ScalarDamageInterface3DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
		if (rThisVariable == INITIAL_STRAIN) {
			if (rValue.size() != m_initial_strain.size())
				rValue.resize(m_initial_strain.size());
			noalias(rValue) = m_initial_strain;
		}

		if(rThisVariable == YIELD_SURFACE_DATA_3D_X || rThisVariable == YIELD_SURFACE_DATA_3D_Y || rThisVariable == YIELD_SURFACE_DATA_3D_Z)
		{
			int ijob = rThisVariable == YIELD_SURFACE_DATA_3D_X ? 1 : 2;
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

    Matrix& ScalarDamageInterface3DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
        return rValue;
    }

    array_1d<double, 3 > & ScalarDamageInterface3DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    array_1d<double, 6 > & ScalarDamageInterface3DLaw::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
    {
        return rValue;
    }

    void ScalarDamageInterface3DLaw::SetValue(const Variable<double>& rVariable,
                                              const double& rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
		if(rVariable == DAMAGE_T)
			mD1 = rValue;
		if(rVariable == DAMAGE_C)
			mD2 = rValue;
    }

    void ScalarDamageInterface3DLaw::SetValue(const Variable<Vector >& rVariable,
                                              const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
		if (rVariable == INITIAL_STRAIN) {
			if (rValue.size() == m_initial_strain.size())
				noalias(m_initial_strain) = rValue;
		}
    }

    void ScalarDamageInterface3DLaw::SetValue(const Variable<Matrix >& rVariable,
                                              const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface3DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                              const array_1d<double, 3 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface3DLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                              const array_1d<double, 6 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
    }

    bool ScalarDamageInterface3DLaw::ValidateInput(const Properties& rMaterialProperties)
    {
		if( !rMaterialProperties.Has(NORMAL_STIFFNESS) ) return false;
		if( !rMaterialProperties.Has(TANGENTIAL_STIFFNESS) ) return false;
		if( !rMaterialProperties.Has(YIELD_STRESS_T) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_I) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_II) ) return false;
		if( !rMaterialProperties.Has(INITIAL_COHESION) ) return false;
		if( !rMaterialProperties.Has(INTERNAL_FRICTION_ANGLE) ) return false;
        return true;
    }

    ScalarDamageInterface3DLaw::StrainMeasure ScalarDamageInterface3DLaw::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }

    ScalarDamageInterface3DLaw::StressMeasure ScalarDamageInterface3DLaw::GetStressMeasure()
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    bool ScalarDamageInterface3DLaw::IsIncremental()
    {
        return false;
    }

    void ScalarDamageInterface3DLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const Vector& rShapeFunctionsValues)
    {
		if(!mInitialized)
		{
			mK1 = 0.0;
			mK1_converged = 0.0;
			mK2 = 0.0;
			mK2_converged = 0.0;
			mK3 = 0.0;
			mK3_converged = 0.0;
			mD1 = 0.0;
			mD2 = 0.0;
			mD2_bar = 0.0;
			mD2_bar_converged = 0.0;
			mD3 = 0.0;
			mD3_bar = 0.0;
			mD3_bar_converged = 0.0;
			mYieldValue = 0.0;
			m_initial_strain = ZeroVector(this->GetStrainSize());
			mInitialized = true;
		}
    }

    void ScalarDamageInterface3DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
                                                            const GeometryType& rElementGeometry,
                                                            const Vector& rShapeFunctionsValues,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
		mK1 = mK1_converged;
		mK2 = mK2_converged;
		mK3 = mK3_converged;
		mD2_bar = mD2_bar_converged;
		mD3_bar = mD3_bar_converged;
    }

    void ScalarDamageInterface3DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
		mK1_converged = mK1;
		mK2_converged = mK2;
		mK3_converged = mK3;
		mD2_bar_converged = mD2_bar;
		mD3_bar_converged = mD3_bar;
    }

    void ScalarDamageInterface3DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                                  const GeometryType& rElementGeometry,
                                                                  const Vector& rShapeFunctionsValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface3DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                                const GeometryType& rElementGeometry,
                                                                const Vector& rShapeFunctionsValues,
                                                                const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ScalarDamageInterface3DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
        const Properties& props = rValues.GetMaterialProperties();
		const Vector& strainVector = rValues.GetStrainVector();
		Vector& stressVector = rValues.GetStressVector();
		Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
		Flags& Options = rValues.GetOptions();
		bool compute_constitutive_tensor = Options.Is(COMPUTE_CONSTITUTIVE_TENSOR);
		bool compute_stress = Options.Is(COMPUTE_STRESS) || compute_constitutive_tensor;

		SizeType size = GetStrainSize();
		if(compute_stress)
			if(stressVector.size() != size)
				stressVector.resize(size, false);
		if(compute_constitutive_tensor)
			if(constitutiveMatrix.size1() != size || constitutiveMatrix.size2() != size)
				constitutiveMatrix.resize(size, size, false);

		CalculationData data;

		InitializeCalculationData( props, rValues.GetElementGeometry(), strainVector, data );

		CalculateElasticStressVector( data, strainVector );

		stressVector = data.ElasticStressVector;
		if(compute_constitutive_tensor)
		{
			// elastic case
			constitutiveMatrix(0, 0) = data.Kt1;
			constitutiveMatrix(1, 1) = data.Kt2;
			constitutiveMatrix(2, 2) = data.Kn;
			constitutiveMatrix(0, 1) = constitutiveMatrix(1, 0) = 0.0;
			constitutiveMatrix(0, 2) = constitutiveMatrix(2, 0) = 0.0;
			constitutiveMatrix(1, 2) = constitutiveMatrix(2, 1) = 0.0;
		}

		return;

		CalculateEquivalentMeasure( data );

		UpdateDamage( data );
		mD1 = data.D1;
		mD2 = data.D2;
		mD3 = data.D3;

		if( compute_stress )
		{
			CalculateStress( data, stressVector );
		}

		//**********************************************
		double sig_n = stressVector(2);
		double normtau = std::sqrt(stressVector(0)*stressVector(0) + stressVector(1)*stressVector(1));
		double C0_d = (1.0 - mD2)*data.C0;
		mYieldValue = sig_n*data.Fs + normtau - C0_d;
		//**********************************************

		if( compute_constitutive_tensor )
		{
			CalculateConstitutiveMatrix( data, strainVector, stressVector, constitutiveMatrix );
		}
    }

    void ScalarDamageInterface3DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface3DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface3DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ScalarDamageInterface3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {

    }

    void ScalarDamageInterface3DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                   const GeometryType& rElementGeometry,
                                                   const Vector& rShapeFunctionsValues)
    {
        mInitialized = false;
		mK1 = 0.0;
		mK1_converged = 0.0;
		mK2 = 0.0;
		mK2_converged = 0.0;
		mK3 = 0.0;
		mK3_converged = 0.0;
		mD1 = 0.0;
		mD2 = 0.0;
		mD3 = 0.0;
		mYieldValue = 0.0;
    }

    void ScalarDamageInterface3DLaw::GetLawFeatures(Features& rFeatures)
    {
        //Set the type of law
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mOptions.Set( ISOTROPIC );

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
    }

    int ScalarDamageInterface3DLaw::Check(const Properties& rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

		if( !rMaterialProperties.Has(NORMAL_STIFFNESS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: NORMAL_STIFFNESS", "");

		if( !rMaterialProperties.Has(TANGENTIAL_STIFFNESS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TANGENTIAL_STIFFNESS", "");

		if( !rMaterialProperties.Has(YIELD_STRESS_T) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_T", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_I) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_MODE_I", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_MODE_II) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_MODE_II", "");

		if( !rMaterialProperties.Has(INITIAL_COHESION) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INITIAL_COHESION", "");

		if( !rMaterialProperties.Has(INTERNAL_FRICTION_ANGLE) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: INTERNAL_FRICTION_ANGLE", "");

        return 0;

        KRATOS_CATCH("");
    }

	void ScalarDamageInterface3DLaw::InitializeCalculationData(const Properties& props,
		                                                       const GeometryType& geom,
															   const Vector& strainVector,
															   CalculationData& data)
	{
		data.Kn  = props[NORMAL_STIFFNESS];
		data.Kt1  = props[TANGENTIAL_STIFFNESS];
		data.Kt2  = props[TANGENTIAL_STIFFNESS];

		data.Kn_compression_multiplier = 1.0;
		/*if(props.Has(NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER)) {
			data.Kn_compression_multiplier = props[NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER];
			if(data.Kn_compression_multiplier < 1.0)
				data.Kn_compression_multiplier = 1.0;
			if(strainVector(2) <= 0.0)
				data.Kn *= data.Kn_compression_multiplier;
		}*/

		data.Ft  = props[YIELD_STRESS_T];

		data.C0 = props[INITIAL_COHESION];
		data.Fs = std::tan( props[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0 );

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
		data.D3 = 0.0;
		data.K1 = 0.0;
		data.K2 = 0.0;
		data.K3 = 0.0;

		data.ElasticStressVector.clear();

		data.ForceSecant = false;
		if(props.Has(DAMAGE_SECANT_MATRIX))
			data.ForceSecant = props[DAMAGE_SECANT_MATRIX] != 0;
	}

	void ScalarDamageInterface3DLaw::CalculateElasticStressVector(CalculationData& data,
		                                                          const Vector& strainVector)
	{
		data.ElasticStressVector(0) = data.Kt1 * strainVector(0);
		data.ElasticStressVector(1) = data.Kt2 * strainVector(1);
		data.ElasticStressVector(2) = data.Kn * strainVector(2);
	}

	void ScalarDamageInterface3DLaw::CalculateEquivalentMeasure(CalculationData& data)
	{
		double norm_t = std::sqrt(data.ElasticStressVector(0)*data.ElasticStressVector(0)+
			                       data.ElasticStressVector(1)*data.ElasticStressVector(1));
		double sigma_n = data.ElasticStressVector(2);

		double lambda1 = sigma_n - data.Ft;
		double lambda2 = sigma_n * data.Fs + norm_t - data.C0;

		data.K1 = 0.0;
		data.K2 = 0.0;
		//lambda2 = 0.0; // only for test!!!

		// elastic test
		lambda1 = 0.0;
		lambda2 = 0.0;

		if(lambda1 > 0.0)
		{
			if(lambda2 > 0.0)
			{
				// both active
				data.K1 = std::sqrt( lambda1*lambda1 + (lambda2/data.MRatio)*(lambda2/data.MRatio) );
				data.K2 = std::sqrt( lambda2*lambda2 + (lambda1*data.MRatio)*(lambda1*data.MRatio) );
			}
			else
			{
				// only tension cut-off
				data.K1 = lambda1;
				data.K2 = lambda1 * data.MRatio;
			}
		}
		else
		{
			if(lambda2 > 0.0)
			{
				// only coulomb friction
				data.K1 = lambda2 / data.MRatio;
				data.K2 = lambda2;
			}
		}
		//data.K2 = 0.0; // test
		//data.K3 = 0.0; // test
	}

	void ScalarDamageInterface3DLaw::UpdateDamageIndicator(CalculationData& data)
	{
		mK1 = std::max( data.K1, mK1_converged );
		mK2 = std::max( data.K2, mK2_converged );
		mK3 = std::max( data.K3, mK3_converged );
	}

	void ScalarDamageInterface3DLaw::CalculateDamage(CalculationData& data)
	{
		if(mK1 > 0.0)
		{
			data.D1 = 1.0 - data.Ft/(mK1+data.Ft) * std::exp( -data.Ft/(data.GI*data.Kn) * mK1 );
			data.D1 = std::max( std::min( data.D1, 1.0 ), 0.0 );
		}
		if(mK2 > 0.0)
		{
			data.D2 = 1.0 - data.C0/(mK2+data.C0) * std::exp( -data.C0/data.GII/data.Kt1 * mK2 );
			data.D2 = std::max( std::min( data.D2, 1.0 ), 0.0 );
		}
		if(mK3 > 0.0)
		{
			data.D3 = 1.0 - data.C0/(mK3+data.C0) * std::exp( -data.C0/data.GII/data.Kt2 * mK3 );
			data.D3 = std::max( std::min( data.D3, 1.0 ), 0.0 );
		}
	}

	void ScalarDamageInterface3DLaw::UpdateDamage(CalculationData& data)
	{
		bool update_equ_shear_damage = false;

		mK1 = mK1_converged;
		mK2 = mK2_converged;
		mK3 = mK3_converged;
		mD2_bar = mD2_bar_converged;
		mD3_bar = mD3_bar_converged;

		if(data.K1 > mK1)
		{
			mK1 = data.K1;
		}
		if(data.K2 > mK2)
		{
			mK2 = data.K2;
			update_equ_shear_damage = true;
		}
		if(data.K3 > mK3)
		{
			mK3 = data.K3;
			update_equ_shear_damage = true;
		}

		if(mK1 > 0.0)
		{
			data.D1 = 1.0 - data.Ft/(mK1+data.Ft) * std::exp( -data.Ft/(data.GI*data.Kn) * mK1 );
			data.D1 = std::max( std::min( data.D1, 1.0 ), 0.0 );
		}
		if(mK2 > 0.0)
		{
			data.D2 = 1.0 - data.C0/(mK2+data.C0) * std::exp( -data.C0/data.GII/data.Kt1 * mK2 );
			data.D2 = std::max( std::min( data.D2, 1.0 ), 0.0 );
		}
		if(mK3 > 0.0)
		{
			data.D3 = 1.0 - data.C0/(mK3+data.C0) * std::exp( -data.C0/data.GII/data.Kt2 * mK3 );
			data.D3 = std::max( std::min( data.D3, 1.0 ), 0.0 );
		}

		if(update_equ_shear_damage)
		{
			double sigma_n = data.ElasticStressVector(2);
			double sigma_t1 = data.ElasticStressVector(0);
			double sigma_t2 = data.ElasticStressVector(1);
            if(sigma_t1 != 0.0)
			{
				double C0_d = (1.0 - data.D2)*data.C0;
				double tau = std::max(C0_d - sigma_n*data.Fs, 0.0);

                if(sigma_t1 < 0.0) tau = -tau;

				if(std::abs(sigma_t1) > std::abs(tau))
				{
					mD2_bar = 1.0 - tau / sigma_t1;
					mD2_bar = std::max(std::min(mD2_bar,1.0),0.0);
                }
			}
            if(sigma_t2 != 0.0)
            {
                double C0_d = (1.0 - data.D2)*data.C0;
                double tau = std::max(C0_d - sigma_n*data.Fs, 0.0);

                if(sigma_t2 < 0.0) tau = -tau;

                if(std::abs(sigma_t2) > std::abs(tau))
                {
                    mD3_bar = 1.0 - tau / sigma_t2;
                    mD3_bar = std::max(std::min(mD3_bar,1.0),0.0);
                }
            }
		}
	}

	void ScalarDamageInterface3DLaw::CalculateStress(CalculationData& data,
		                                             Vector& stressVector)
	{
		double sigma_n = data.ElasticStressVector(2);
		double sigma_t1 = data.ElasticStressVector(0);
		double sigma_t2 = data.ElasticStressVector(1);
		//sigma_t1 = 0.0; // TEST ONLY TENSION
		//sigma_t2 = 0.0; // TEST ONLY TENSION

		if(sigma_n > 0.0) sigma_n *= (1.0 - data.D1);
		stressVector(2) = sigma_n;

		stressVector(0) = (1.0 - mD2_bar) * sigma_t1;
		stressVector(1) = (1.0 - mD3_bar) * sigma_t2;
	}

	void ScalarDamageInterface3DLaw::CalculateConstitutiveMatrix(CalculationData& data,
		                                                         const Vector& strainVector,
		                                                         const Vector& stressVector,
																 Matrix& constitutiveMatrix)
	{
		// elastic case
		constitutiveMatrix(0, 0) = data.Kt1;
		constitutiveMatrix(1, 1) = data.Kt2;
		constitutiveMatrix(2, 2) = data.Kn;
		constitutiveMatrix(0, 1) = constitutiveMatrix(1, 0) = 0.0;
		constitutiveMatrix(0, 2) = constitutiveMatrix(2, 0) = 0.0;
		constitutiveMatrix(1, 2) = constitutiveMatrix(2, 1) = 0.0;

		double perturbation = 1.0E-8;
		Vector perturbedStrainVector(3);
		Vector stressPerturbation(3);

		for(int j = 0; j < 3; j++)
		{
			// save internal variables
			double save_k1 = mK1;
			double save_k2 = mK2;
			double save_k3 = mK3;
			double save_d2_bar = mD2_bar;
			double save_d3_bar = mD3_bar;

			// FORWARD difference
			noalias( perturbedStrainVector ) = strainVector;

			//double delta_strain = perturbedStrainVector(j) < 0.0 ? -perturbation : perturbation;
			//perturbedStrainVector(j) += delta_strain;
			perturbedStrainVector(j) += perturbation;

			CalculateElasticStressVector( data, perturbedStrainVector );
			CalculateEquivalentMeasure( data );
			UpdateDamage( data );
			CalculateStress( data, stressPerturbation );

			// fill the numerical tangent operator
			noalias( stressPerturbation ) -= stressVector;
			constitutiveMatrix(0, j) = stressPerturbation(0) / perturbation;
			constitutiveMatrix(1, j) = stressPerturbation(1) / perturbation;
			constitutiveMatrix(2, j) = stressPerturbation(2) / perturbation;

			// restore internal variables
			mK1 = save_k1;
			mK2 = save_k2;
			mK3 = save_k3;
			mD2_bar = save_d2_bar;
			mD3_bar = save_d3_bar;
		}
	}

} /* namespace Kratos.*/
