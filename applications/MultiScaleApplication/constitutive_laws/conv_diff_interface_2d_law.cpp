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

#include "conv_diff_interface_2d_law.h"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"

//#define EXPONENTIAL_DAMAGE

namespace Kratos
{

    ConvDiffInterface2DLaw::ConvDiffInterface2DLaw()
        : ConstitutiveLaw()
		, mInitialized(false)
		, m_init_gradT()
		, mStressVector()
		, m_gap_interface()
	{
    }

    ConstitutiveLaw::Pointer ConvDiffInterface2DLaw::Clone() const
    {
        return ConstitutiveLaw::Pointer( new ConvDiffInterface2DLaw() );
    }

    ConvDiffInterface2DLaw::SizeType ConvDiffInterface2DLaw::WorkingSpaceDimension()
    {
        return 2;
    }

    ConvDiffInterface2DLaw::SizeType ConvDiffInterface2DLaw::GetStrainSize()
    {
        return 2;
    }

    bool ConvDiffInterface2DLaw::Has(const Variable<double>& rThisVariable)
    {
		if (rThisVariable == CONVECTION_DEGRADATION)
			return true;
		if(rThisVariable == GAP_INTERFACE)
			return true;
		if(rThisVariable == YIELD_FUNCTION_VALUE)
			return true;
        return false;
    }

    bool ConvDiffInterface2DLaw::Has(const Variable<Vector>& rThisVariable)
    {
		if(rThisVariable == YIELD_SURFACE_DATA_2D_X || rThisVariable == YIELD_SURFACE_DATA_2D_Y)
			return true;
		if (rThisVariable == INITIAL_STRAIN)
			return true;
		return false;
    }

    bool ConvDiffInterface2DLaw::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    bool ConvDiffInterface2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }

    bool ConvDiffInterface2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    double& ConvDiffInterface2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
		rValue = 0.0;
		if (rThisVariable == CONVECTION_DEGRADATION)
			rValue = mD1;
		return rValue;
    }

    Vector& ConvDiffInterface2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
		if (rThisVariable == INITIAL_TEMP_GRAD) {
			if (rValue.size() != m_init_gradT.size())
				rValue.resize(m_init_gradT.size());
			noalias(rValue) = m_init_gradT;
		}
		if (rThisVariable == FLUX_RVE) {
			if (rValue.size() != mStressVector.size())
				rValue.resize(mStressVector.size());
			noalias(rValue) = mStressVector;
		}
		if (rThisVariable == HEAT_FLUX_RVE) { //For Output
			if (rValue.size() != 3)
				rValue.resize(3);
			rValue(0) = mStressVector(0); // 1.0e6;//[W/mm^2]
			rValue(1) = mStressVector(1); // 1.0e6;
			rValue(2) = 0.0;
		}
		if (rThisVariable == GAP_INTERFACE) {
			if (rValue.size() != m_gap_interface.size())
				rValue.resize(m_gap_interface.size());
			noalias(rValue) = m_gap_interface;
		}
        return rValue;
    }

    Matrix& ConvDiffInterface2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
        return rValue;
    }

    array_1d<double, 2 > & ConvDiffInterface2DLaw::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
    {
        return rValue;
    }

    array_1d<double, 3 > & ConvDiffInterface2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    void ConvDiffInterface2DLaw::SetValue(const Variable<double>& rVariable,
                                              const double& rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == CONVECTION_DEGRADATION)
			mD1 = rValue;
    }

    void ConvDiffInterface2DLaw::SetValue(const Variable<Vector >& rVariable,
                                              const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
		if (rVariable == FLUX_RVE)
		{
			if (rValue.size() == mStressVector.size())
				noalias(mStressVector) = rValue;
		}
		if (rVariable == GAP_INTERFACE)
		{
			if (rValue.size() == m_gap_interface.size())
				noalias(m_gap_interface) = rValue;
		}
    }

    void ConvDiffInterface2DLaw::SetValue(const Variable<Matrix >& rVariable,
                                              const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ConvDiffInterface2DLaw::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
                                              const array_1d<double, 2 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {

    }

    void ConvDiffInterface2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                              const array_1d<double, 3 > & rValue,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
    }

    bool ConvDiffInterface2DLaw::ValidateInput(const Properties& rMaterialProperties)
    {
		//if (!rMaterialProperties.Has(CONDUCTIVITY)) return false;
		//if (!rMaterialProperties.Has(NORMAL_STIFFNESS)) return false; //Same of the solid phase
		//if (!rMaterialProperties.Has(YIELD_STRESS_T)) return false; //Same of the solid phase
		//if (!rMaterialProperties.Has(FRACTURE_ENERGY_MODE_I)) return false; //Same of the solid phase
        return true;
    }

    ConvDiffInterface2DLaw::StrainMeasure ConvDiffInterface2DLaw::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }

    ConvDiffInterface2DLaw::StressMeasure ConvDiffInterface2DLaw::GetStressMeasure()
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    bool ConvDiffInterface2DLaw::IsIncremental()
    {
        return false;
    }

    void ConvDiffInterface2DLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const Vector& rShapeFunctionsValues)
    {
		if(!mInitialized)
		{
			mK1 = 0.0;
			mK1_converged = 0.0;
			mD1 = 0.0;
			mYieldValue = 0.0;
			m_init_gradT = ZeroVector(this->GetStrainSize());
			mStressVector = ZeroVector(this->GetStrainSize());
			m_gap_interface = ZeroVector(this->GetStrainSize());
			mInitialized = true;
		}
    }

    void ConvDiffInterface2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
                                                            const GeometryType& rElementGeometry,
                                                            const Vector& rShapeFunctionsValues,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
		mK1 = mK1_converged;
    }

    void ConvDiffInterface2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
		mK1_converged = mK1;
    }

    void ConvDiffInterface2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                                  const GeometryType& rElementGeometry,
                                                                  const Vector& rShapeFunctionsValues,
                                                                  const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ConvDiffInterface2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                                const GeometryType& rElementGeometry,
                                                                const Vector& rShapeFunctionsValues,
                                                                const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void ConvDiffInterface2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ConvDiffInterface2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ConvDiffInterface2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void ConvDiffInterface2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
        const Properties& props = rValues.GetMaterialProperties();
		Vector& strainVector = rValues.GetStrainVector();
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

		InitializeCalculationData( props, rValues.GetElementGeometry(), strainVector, rValues.GetProcessInfo(), data );

		std::stringstream ss;
		//std::cout << "CalculateMaterialResponseCauchy -  strainVector = " << strainVector << std::endl;
		//std::cout << "CalculateMaterialResponseCauchy -  constitutiveMatrix = " << constitutiveMatrix << std::endl;

		if (data.ExpCurveTypeFlag)
		{
			CalculateElasticStressVector(data, strainVector);
			//std::cout << "CalculateMaterialResponseCauchy -  ElasticStressVector = " << data.ElasticStressVector << std::endl;
			CalculateEquivalentMeasure(data);
			UpdateDamage(data);
			CalculateContactConductivity(data);
		}
		else
		{
			CalculateEffectiveOpening(data);
			UpdateDamage(data);
			CalculateContactConductivity(data);
		}

		mD1 = data.D1; // Update the mK1 to the current value

		if( compute_stress )
		{
			CalculateStress(data, strainVector, stressVector);
		}
		//std::cout << "CalculateMaterialResponseCauchy -  CalculateStress = " << stressVector << std::endl;

		if( compute_constitutive_tensor )
		{
			CalculateConstitutiveMatrix( data, strainVector, stressVector, constitutiveMatrix );
		}

		std::cout << ss.str();
    }

    void ConvDiffInterface2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ConvDiffInterface2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ConvDiffInterface2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void ConvDiffInterface2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {

    }

    void ConvDiffInterface2DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                                   const GeometryType& rElementGeometry,
                                                   const Vector& rShapeFunctionsValues)
    {
        mInitialized = false;
		mK1 = 0.0;
		mK1_converged = 0.0;
		mD1 = 0.0;
		mYieldValue = 0.0;
    }

    void ConvDiffInterface2DLaw::GetLawFeatures(Features& rFeatures)
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

    int ConvDiffInterface2DLaw::Check(const Properties& rMaterialProperties,
                                          const GeometryType& rElementGeometry,
                                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

			if (!rMaterialProperties.Has(CONDUCTIVITY))
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: CONDUCTIVITY", "");

		if (!rMaterialProperties.Has(AIR_CONDUCTIVITY))
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: AIR_CONDUCTIVITY", "");

		if (!rMaterialProperties.Has(ASPERITIES_CONSTANT))
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: ASPERITIES_CONSTANT", "");

        return 0;

        KRATOS_CATCH("");
    }

	void ConvDiffInterface2DLaw::InitializeCalculationData(const Properties& props,
														   const GeometryType& geom,
														   const Vector& strainVector,
														   const ProcessInfo& pinfo,
														   CalculationData& data)
	{
		// Check EXPONENTIAL\LINEAR DAMAGE
		data.ExpCurveTypeFlag = false;
		if (props.Has(EXPONENTIAL_DAMAGE))
			data.ExpCurveTypeFlag = props[EXPONENTIAL_DAMAGE] != 0;

		if (data.ExpCurveTypeFlag)
		{
			data.Kn = props[NORMAL_STIFFNESS];
			data.Kt = props[TANGENTIAL_STIFFNESS];

			data.Ks = props[CONDUCTIVITY];
			data.Kg = props[AIR_CONDUCTIVITY];
			//data.Kt = 0.0; // TEST ONLY TENSION

			data.Kn_compression_multiplier = 1.0;
			if (props.Has(NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER))
				data.Kn_compression_multiplier = props[NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER];

			// boundary materials properties
			data.m_exp = props[ASPERITIES_CONSTANT];
			data.E_1 = props[YOUNG_MAT1];
			data.K_1 = props[CONDUCT_MAT1];
			data.nu_1 = props[POISSON_MAT1];
			data.E_2 = props[YOUNG_MAT2];
			data.K_2 = props[CONDUCT_MAT2];
			data.nu_2 = props[POISSON_MAT2];

			//double Ft_mec = props[YIELD_STRESS_T];
			//double C0_mec = props[INITIAL_COHESION];
			//double normal_scaling_factor = Ft_mec / data.Ks;
			//double tangential_scaling_factor = C0_mec / data.Kg;

			data.Ft = props[YIELD_STRESS_T]; //Ft_mec * normal_scaling_factor;

			data.GI = props[FRACTURE_ENERGY_MODE_I]; // *normal_scaling_factor;

			data.D1 = 0.0;
			data.K1 = 0.0;

			data.ElasticStressVector.clear();
		}
		else
		{
			// Leggo le prop meccaniche
			data.Kn_compression_multiplier = props[NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER];
			data.m_exp = props[ASPERITIES_CONSTANT];
			data.beta = props[BETA_THERMAL_COEFF];
			data.norm_stiff = props[NORMAL_STIFFNESS];
			data.tang_stiff = props[TANGENTIAL_STIFFNESS];
			data.Ft = props[YIELD_STRESS_T];
			data.C0 = props[INITIAL_COHESION];
			// boundary materials properties
			data.E_1 = props[YOUNG_MAT1];
			data.K_1 = props[CONDUCT_MAT1];
			data.nu_1 = props[POISSON_MAT1];
			data.E_2 = props[YOUNG_MAT2];
			data.K_2 = props[CONDUCT_MAT2];
			data.nu_2 = props[POISSON_MAT2];

			// Modifico Proporzionalmente le proprieta per degradare allo stesso modo del meccanico
			data.Ks = props[CONDUCTIVITY];
			data.Kg = props[AIR_CONDUCTIVITY];

			data.D1 = 0.0;
			data.K1 = 0.0;

			CalculateCriticalOpening(data);
		}

		data.ForceSecant = false;
		if(props.Has(DAMAGE_SECANT_MATRIX))
			data.ForceSecant = props[DAMAGE_SECANT_MATRIX] != 0;
	}

	void ConvDiffInterface2DLaw::CalculateContactConductivity(CalculationData& data)
	{
		double ro_1 = 1 / data.K_1;
		double ro_2 = 1 / data.K_2;
		double aux1 = (2 / data.E_1)*(1.0 - pow(data.nu_1,2));
		double aux2 = (2 / data.E_2)*(1.0 - pow(data.nu_2,2));
		double gn = std::abs(m_gap_interface(1))*(m_gap_interface(1)); // ? 0 : 1);
		//std::cout << MathHelpers::VectorToString(m_gap_interface, 4, std::scientific);
		data.Kc = (aux1 + aux2)*(1 / (ro_1 + ro_2))*data.Kn_compression_multiplier*data.m_exp*pow(gn,data.m_exp - 1);

		//std::stringstream ss;
		//ss << "--------  ro_1 = " << ro_1 << ", " << std::endl;
		//ss << "--------  ro_2 = " << ro_2 << ", " << std::endl;
		//ss << "--------  aux1 = " << aux1 << ", " << std::endl;
		//ss << "--------  aux2 = " << aux2 << ", " << std::endl;
		//if (gn > 0)
		//	ss << "--------  gn = " << gn << ", " << std::endl;
			//ss << "--------  data.Kc = " << data.Kc << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffInterface2DLaw::CalculateCriticalOpening(CalculationData& data)
	{
		double critical_normal_gap = data.Ft / data.norm_stiff;
		double critical_tangential_gap = data.C0 / data.tang_stiff;
		data.crticalDispl = std::sqrt(pow(data.beta, 2)*pow(critical_tangential_gap, 2) + pow(critical_normal_gap, 2));

		//std::stringstream ss;
		//ss << "--------  critical_normal_gap = " << critical_normal_gap << ", " << std::endl;
		//ss << "--------  crticalDispl = " << data.crticalDispl << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffInterface2DLaw::CalculateEffectiveOpening(CalculationData& data)
	{
		data.effDispl = std::sqrt(pow(data.beta, 2)*pow(m_gap_interface(0), 2) + pow(m_gap_interface(1), 2));

		//std::stringstream ss;
		//ss << "--------  effDispl = " << data.effDispl << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffInterface2DLaw::UpdateDamage(CalculationData& data)
	{
		if (data.ExpCurveTypeFlag)
		{
			mK1 = mK1_converged;

			if (data.K1 > mK1)
			{
				mK1 = data.K1; // Update the mK1 to the current value
			}
			if (mK1 > 0.0)
			{
				data.D1 = 1.0 - data.Ft / (mK1 + data.Ft) * std::exp(-data.Ft / (data.GI*data.Ks) * mK1);
				data.D1 = std::max(std::min(data.D1, 1.0), 0.0);
			}
		}
		else
		{
			data.D1 = data.effDispl / data.crticalDispl;
			data.D1 = std::max(std::min(data.D1, 1.0), 0.0);
		}

		//std::stringstream ss;
		//if (data.D1 > 1.0e-8)
		//	ss << "-------- Therm D1 = " << data.D1 << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffInterface2DLaw::CalculateElasticStressVector(CalculationData& data,
		const Vector& strainVector)
	{
		data.ElasticStressVector(0) = data.Ks * strainVector(0);
		data.ElasticStressVector(1) = data.Ks * strainVector(1);
	}

	void ConvDiffInterface2DLaw::CalculateEquivalentMeasure(CalculationData& data)
	{
		double sigma_n = data.ElasticStressVector(1);

		double lambda1 = sigma_n - data.Ft;

		data.K1 = 0.0;
		if (lambda1 > 0.0)
			data.K1 = lambda1;
	}

	void ConvDiffInterface2DLaw::UpdateDamageIndicator(CalculationData& data)
	{
		mK1 = std::max(data.K1, mK1_converged);
	}

	void ConvDiffInterface2DLaw::CalculateStress(CalculationData& data,
													 const Vector& strainVector,
		                                             Vector& stressVector)
	{
		double sigma_n;
		double sigma_t;

		if (data.ExpCurveTypeFlag)
		{
			sigma_t = (1.0 - data.D1) * data.ElasticStressVector(0);
			sigma_n = data.ElasticStressVector(1);

			if (m_gap_interface(1) /*sigma_n*/ <= 0.0) {
				sigma_n = ((1.0 - data.D1)*data.Ks + data.Ks /*Kc*/) * strainVector(1);
			}
			else {
				sigma_n *= (1.0 - data.D1);
			}
		}
		else
		{
			sigma_n = ((1.0 - data.D1)*data.Ks + data.D1*data.Kc + data.Kg)*strainVector(1);
			sigma_t = (1.0 - data.D1)*data.Ks*strainVector(0);
		}
		mStressVector(0) = sigma_n;
		mStressVector(1) = sigma_t;
		noalias(stressVector) = mStressVector;

		/*std::stringstream ss;
		ss << "--------  data.Kc = " << data.Kc << ", " << std::endl;
		ss << "--------  data.D1 = " << data.D1 << ", " << std::endl;
		ss << "--------  m_gap_interface = " << m_gap_interface << ", " << std::endl;
		ss << "PlaneStressIterf - stressVector = " << stressVector << ", " << std::endl;
		std::cout << ss.str();*/
	}

	void ConvDiffInterface2DLaw::CalculateConstitutiveMatrix(CalculationData& data,
		                                                         const Vector& strainVector,
		                                                         const Vector& stressVector,
																 Matrix& constitutiveMatrix)
	{
		//std::stringstream ss;

		size_t n = GetStrainSize();

		double perturbation = 1.0E-6;
		Vector perturbedStrainVector(n);
		Vector stressPerturbation(n);

		for (size_t j = 0; j < n; j++)
		{
			// save internal variables
			double save_k1 = mK1;

			// FORWARD difference
			noalias(perturbedStrainVector) = strainVector;

			perturbedStrainVector(j) += perturbation;

			if (data.ExpCurveTypeFlag)
			{
				CalculateElasticStressVector(data, perturbedStrainVector);
				CalculateEquivalentMeasure(data);
				UpdateDamage(data);
				CalculateContactConductivity(data);
			}
			else
			{
				CalculateEffectiveOpening(data);
				UpdateDamage(data);
				CalculateContactConductivity(data);
			}

			CalculateStress(data, perturbedStrainVector, stressPerturbation);

			// fill the numerical tangent operator
			noalias(stressPerturbation) -= stressVector;
			//std::cout << "--------  stressPerturbation = " << stressPerturbation << std::endl;
			constitutiveMatrix(0, j) = stressPerturbation(0) / perturbation;
			constitutiveMatrix(1, j) = stressPerturbation(1) / perturbation;

			// restore internal variables
			mK1 = save_k1;
		}
		//std::cout << "--------  strainVector = " << strainVector << std::endl;
		//std::cout << "--------  constitutiveMatrix = " << constitutiveMatrix << std::endl;
		//std::cout << ss.str();
	}

} /* namespace Kratos.*/
