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

#include "isotropic_damage_plane_stress_2d_law.h"
#include "multiscale_application.h"
#include "custom_utilities/math_helpers.h"
#include "includes/variables.h"

namespace Kratos
{

    IsotropicDamagePlaneStress2DLaw::IsotropicDamagePlaneStress2DLaw() 
        : ConstitutiveLaw()
		, mInitialized(false)
		, mK(0.0)
		, mK_converged(0.0)
		, mDamageT(0.0)
		, mDamageC(0.0)
		, mClen0(0.0)
		, mClen0_multiplier(1.0)
    {
		
    }
 
    ConstitutiveLaw::Pointer IsotropicDamagePlaneStress2DLaw::Clone() const
    {
        return ConstitutiveLaw::Pointer( new IsotropicDamagePlaneStress2DLaw() );
    }

    IsotropicDamagePlaneStress2DLaw::SizeType IsotropicDamagePlaneStress2DLaw::WorkingSpaceDimension()
    {
        return 2;
    }

    IsotropicDamagePlaneStress2DLaw::SizeType IsotropicDamagePlaneStress2DLaw::GetStrainSize()
    {
        return 3;
    }

    bool IsotropicDamagePlaneStress2DLaw::Has(const Variable<double>& rThisVariable)
    {
		if(rThisVariable == DAMAGE_T)
			return true;
		if(rThisVariable == DAMAGE_C)
			return true;
        return false;
    }
    
    bool IsotropicDamagePlaneStress2DLaw::Has(const Variable<Vector>& rThisVariable)
    {
        return false;
    }

    bool IsotropicDamagePlaneStress2DLaw::Has(const Variable<Matrix>& rThisVariable)
    {
        return false;
    }

    bool IsotropicDamagePlaneStress2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
    {
        return false;
    }
    
    bool IsotropicDamagePlaneStress2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
    {
        return false;
    }

    double& IsotropicDamagePlaneStress2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
		rValue = 0.0;
		if(rThisVariable == DAMAGE_T)
			rValue = mDamageT;
		else if(rThisVariable == DAMAGE_C)
			rValue = mDamageC;
        return rValue;
    }

    Vector& IsotropicDamagePlaneStress2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
    {
        return rValue;
    }

    Matrix& IsotropicDamagePlaneStress2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
    {
        return rValue;
    }

    array_1d<double, 3 > & IsotropicDamagePlaneStress2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
    {
        return rValue;
    }

    array_1d<double, 6 > & IsotropicDamagePlaneStress2DLaw::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
    {
        return rValue;
    }

     void IsotropicDamagePlaneStress2DLaw::SetValue(const Variable<double>& rVariable,
                                         const double& rValue,
                                         const ProcessInfo& rCurrentProcessInfo)
    {
		if(rVariable == DAMAGE_T)
			mDamageT = rValue;
		else if(rVariable == CHARACTERISTIC_LENGTH_MULTIPLIER)
			mClen0_multiplier = rValue;
    }

    void IsotropicDamagePlaneStress2DLaw::SetValue(const Variable<Vector >& rVariable,
                                        const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicDamagePlaneStress2DLaw::SetValue(const Variable<Matrix >& rVariable,
                                        const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicDamagePlaneStress2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                                         const array_1d<double, 3 > & rValue,
                                         const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicDamagePlaneStress2DLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                                        const array_1d<double, 6 > & rValue,
                                        const ProcessInfo& rCurrentProcessInfo)
    {
    }

    bool IsotropicDamagePlaneStress2DLaw::ValidateInput(const Properties& rMaterialProperties)
    {
		if( !rMaterialProperties.Has(YOUNG_MODULUS) ) return false;
		if( !rMaterialProperties.Has(POISSON_RATIO) ) return false;
		if( !rMaterialProperties.Has(YIELD_STRESS_T) ) return false;
		if( !rMaterialProperties.Has(YIELD_STRESS_C) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_T) ) return false;
		if( !rMaterialProperties.Has(FRACTURE_ENERGY_C) ) return false;
        return true;
    }

    IsotropicDamagePlaneStress2DLaw::StrainMeasure IsotropicDamagePlaneStress2DLaw::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }
    
    IsotropicDamagePlaneStress2DLaw::StressMeasure IsotropicDamagePlaneStress2DLaw::GetStressMeasure()
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }

    bool IsotropicDamagePlaneStress2DLaw::IsIncremental()
    {
        return false;
    }

    void IsotropicDamagePlaneStress2DLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                  const GeometryType& rElementGeometry,
                                                  const Vector& rShapeFunctionsValues)
    {
		if(!mInitialized)
		{
			mK = 0.0;
			mDamageT = 0.0;
			mDamageC = 0.0;
			mK_converged = 0.0;
			
			mClen0 = rElementGeometry.Length();

			mInitialized = true;
		}
    }

    void IsotropicDamagePlaneStress2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
                                                      const GeometryType& rElementGeometry,
                                                      const Vector& rShapeFunctionsValues,
                                                      const ProcessInfo& rCurrentProcessInfo)
    {
		mK = mK_converged;
    }

    void IsotropicDamagePlaneStress2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
                                                    const GeometryType& rElementGeometry,
                                                    const Vector& rShapeFunctionsValues,
                                                    const ProcessInfo& rCurrentProcessInfo)
    {
		mK_converged = mK;
    }

    void IsotropicDamagePlaneStress2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
                                                            const GeometryType& rElementGeometry,
                                                            const Vector& rShapeFunctionsValues,
                                                            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicDamagePlaneStress2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
                                                          const GeometryType& rElementGeometry,
                                                          const Vector& rShapeFunctionsValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IsotropicDamagePlaneStress2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void IsotropicDamagePlaneStress2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void IsotropicDamagePlaneStress2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }

    void IsotropicDamagePlaneStress2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
    {
        const Properties& props = rValues.GetMaterialProperties();
		const Vector& strainVector = rValues.GetStrainVector();
		Vector& stressVector = rValues.GetStressVector();
		Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
		Flags& Options = rValues.GetOptions();
		bool compute_constitutive_tensor = Options.Is(COMPUTE_CONSTITUTIVE_TENSOR);
		bool compute_stress = Options.Is(COMPUTE_STRESS) || compute_constitutive_tensor;

		CalculationData data;

		InitializeCalculationData( props, rValues.GetElementGeometry(), data );

		CalculateElasticConstitutiveMatrix( data );

		noalias( data.ElasticStressVector ) = prod( data.ElasticConstitutiveMatrix, strainVector );

		CalculateStressPrincipalValues( data.ElasticStressVector, data.ElasticStressVectorPrincipalValues );

		CalculateEquivalentMeasure( data );

		UpdateDamageIndicator( data );

		CalculateDamage( data );
		mDamageT = data.DamageT;
		mDamageC = data.DamageC;

		if( compute_stress )
		{
			CalculateStress( data, stressVector );
		}

		if( compute_constitutive_tensor )
		{
			CalculateConstitutiveMatrix( data, strainVector, stressVector, constitutiveMatrix );
		}
    }

    void IsotropicDamagePlaneStress2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void IsotropicDamagePlaneStress2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void IsotropicDamagePlaneStress2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
    {
        FinalizeMaterialResponseCauchy(rValues);
    }

    void IsotropicDamagePlaneStress2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
    {
        
    }

    void IsotropicDamagePlaneStress2DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const Vector& rShapeFunctionsValues)
    {
        mInitialized = false;
		mK = 0.0;
		mK_converged = 0.0;
		mDamageT = 0.0;
		mDamageC = 0.0;
    }

    void IsotropicDamagePlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
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

    int IsotropicDamagePlaneStress2DLaw::Check(const Properties& rMaterialProperties,
                                            const GeometryType& rElementGeometry,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

		if( !rMaterialProperties.Has(YOUNG_MODULUS) )
			KRATOS_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS", "");

		if( !rMaterialProperties.Has(POISSON_RATIO) )
			KRATOS_ERROR(std::logic_error, "Missing variable: POISSON_RATIO", "");

		if( !rMaterialProperties.Has(YIELD_STRESS_T) )
			KRATOS_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_T", "");

		if( !rMaterialProperties.Has(YIELD_STRESS_C) )
			KRATOS_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_C", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_T) )
			KRATOS_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_T", "");

		if( !rMaterialProperties.Has(FRACTURE_ENERGY_C) )
			KRATOS_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_C", "");
		
		mClen0 = rElementGeometry.Length();
		CalculationData data;
		InitializeCalculationData(rMaterialProperties, rElementGeometry, data);

		if(data.Ft < 0.0)
			KRATOS_ERROR(std::logic_error, "YIELD_STRESS_T should be a positive real number", "");
		if(data.Gt < 0.0)
			KRATOS_ERROR(std::logic_error, "FRACTURE_ENERGY_T should be a positive real number", "");

		if(data.Ft > 0.0 && data.Gt > 0.0)
		{
			double lt = 2.0 * data.E * data.Gt / (data.Ft * data.Ft);
			if(data.CLen >= lt) 
			{
				std::stringstream ss;
				ss << "FRACTURE_ENERGY_T is to low:  2*E*Gt/(Fy*Fy) = " << lt 
				   << ",   Characteristic Length = " << data.CLen << std::endl;
				KRATOS_ERROR(std::logic_error, ss.str(), "");
			}
		}

		if(data.Fc > 0.0 && data.Gc > 0.0)
		{
			double lc = 2.0 * data.E * data.Gc / (data.Fc * data.Fc);
			if(data.CLen >= lc) 
			{
				std::stringstream ss;
				ss << "FRACTURE_ENERGY_C is to low:  2*E*Gc/(Fy*Fy) = " << lc 
				   << ",   Characteristic Length = " << data.CLen << std::endl;
				KRATOS_ERROR(std::logic_error, ss.str(), "");
			}
		}

        return 0;

        KRATOS_CATCH("");
    }

	void IsotropicDamagePlaneStress2DLaw::InitializeCalculationData(const Properties& props, 
		                                                            const GeometryType& geom, 
																    CalculationData& data)
	{
		data.Formulation = DamageFormulationType_Rankine;
		if(props.Has(DAMAGE_MODEL))
		{
			int theFormulation = props[DAMAGE_MODEL];
			if(theFormulation >= 0 || theFormulation < DamageFormulationType_End)
				data.Formulation = (DamageFormulationType)theFormulation;
		}
		data.E  = props[YOUNG_MODULUS];
		data.nu = props[POISSON_RATIO];
		data.Ft = props[YIELD_STRESS_T];
		data.Fc = props[YIELD_STRESS_C];
		data.Gt = props[FRACTURE_ENERGY_T];
		data.Gc = props[FRACTURE_ENERGY_C];

		data.DamageT = 0.0;
		data.DamageC = 0.0;
		data.EquivalentMeasure = 0.0;

		data.ElasticConstitutiveMatrix = ZeroMatrix(3, 3);
		data.ElasticStressVector = ZeroVector(3);
		data.ElasticStressVectorPrincipalValues = ZeroVector(2);

		data.CLen = mClen0;
		data.CLen *= mClen0_multiplier;

		data.ForceSecant = false;
		if(props.Has(DAMAGE_SECANT_MATRIX))
			data.ForceSecant = props[DAMAGE_SECANT_MATRIX] != 0;

		switch (data.Formulation)
		{

		case DamageFormulationType_MohrCoulomb:
			data.R0 = data.Fc;
			data.G = data.Gc;
			break;

		default: // DamageFormulationType_Rankine
			data.R0 = data.Ft;
			data.G = data.Gt;
			break;
		}
	}

	void IsotropicDamagePlaneStress2DLaw::CalculateElasticConstitutiveMatrix(CalculationData& data)
	{
		double c1 = data.E / (1.0 - data.nu * data.nu);
		double c2 = c1 * data.nu;
		double c3 = c1 * (1.0 - data.nu) / 2.0;

		data.ElasticConstitutiveMatrix(0, 0) =  c1;
		data.ElasticConstitutiveMatrix(0, 1) =  c2;
		data.ElasticConstitutiveMatrix(0, 2) = 0.0;

		data.ElasticConstitutiveMatrix(1, 0) =  c2;
		data.ElasticConstitutiveMatrix(1, 1) =  c1;
		data.ElasticConstitutiveMatrix(1, 2) = 0.0;

		data.ElasticConstitutiveMatrix(2, 0) = 0.0;
		data.ElasticConstitutiveMatrix(2, 1) = 0.0;
		data.ElasticConstitutiveMatrix(2, 2) =  c3;
	}

	void IsotropicDamagePlaneStress2DLaw::CalculateStressPrincipalValues(const Vector& S, 
		                                                              Vector& Sp)
	{
		double sig_m = 0.5 * (S(0) + S(1));
		double sig_r = std::sqrt( S(2) * S(2) + 0.25 * ( S(0) - S(1) ) * ( S(0) - S(1)) );
		Sp(0) = sig_m + sig_r;
		Sp(1) = sig_m - sig_r;
	}

	void IsotropicDamagePlaneStress2DLaw::CalculateEquivalentMeasure(CalculationData& data)
	{
		if(data.Formulation == DamageFormulationType_MohrCoulomb)
		{
			double m = data.Fc / data.Ft;
			double K = (m-1.0)/(m+1.0);
			double s1 = data.ElasticStressVectorPrincipalValues(0);
			double s2 = data.ElasticStressVectorPrincipalValues(1);
			double s3 = 0.0;
			data.EquivalentMeasure = (m+1.0)/2.0*std::max(std::abs(s1-s2)+K*(s1+s2), 
												 std::max(std::abs(s1-s3)+K*(s1+s3), 
														  std::abs(s2-s3)+K*(s2+s3)));
		}
		else // DamageFormulationType_Rankine
		{
			double s1 = data.ElasticStressVectorPrincipalValues(0);
			data.EquivalentMeasure = std::max(s1, 0.0);
		}
	}

	void IsotropicDamagePlaneStress2DLaw::UpdateDamageIndicator(CalculationData& data)
	{
		mK = mK_converged;
		if(data.EquivalentMeasure > mK)
			mK = data.EquivalentMeasure;
	}

	void IsotropicDamagePlaneStress2DLaw::CalculateDamage(CalculationData& data)
	{
		if(mK > data.R0)
		{
			double lt = 2.0*data.E * data.G / (data.R0 * data.R0);
			double Hs = data.CLen / ( lt - data.CLen );

			data.DamageT = 1.0 - data.R0 / mK * std::exp( - 2.0*Hs * (mK - data.R0) / data.R0 );

			data.DamageT = std::max(std::min(data.DamageT,1.0),0.0);
		}
	}

	void IsotropicDamagePlaneStress2DLaw::CalculateStress(CalculationData& data, 
		                                               Vector& stressVector)
	{
		noalias( stressVector ) = (1.0 - data.DamageT) * data.ElasticStressVector;
	}

	void IsotropicDamagePlaneStress2DLaw::CalculateConstitutiveMatrix(CalculationData& data, 
		                                                           const Vector& strainVector,
		                                                           const Vector& stressVector, 
																   Matrix& constitutiveMatrix)
	{
		if(data.DamageT == 0.0) { // elastic
			noalias( constitutiveMatrix ) = data.ElasticConstitutiveMatrix;
			return;
		}
		else if( mK < mK_converged ) { // unloading
			noalias( constitutiveMatrix ) = (1.0 - data.DamageT) * data.ElasticConstitutiveMatrix;
			return;
		}
		if(data.ForceSecant) {
			noalias( constitutiveMatrix ) = (1.0 - data.DamageT) * data.ElasticConstitutiveMatrix;
			return;
		}

		double perturbation = 1.0E-8;

		Vector perturbedStrainVector(3);
		Vector stressPerturbation(3);

		for(int j = 0; j < 3; j++)
		{
			// save internal variables
			double save_k = mK; 

			// FORWARD difference
			noalias( perturbedStrainVector ) = strainVector;

			double delta_strain = perturbation;
			if(perturbedStrainVector(j) < 0.0)
				delta_strain = -delta_strain;
			perturbedStrainVector(j) += delta_strain;

			noalias( data.ElasticStressVector ) = prod( data.ElasticConstitutiveMatrix, perturbedStrainVector );
			CalculateStressPrincipalValues( data.ElasticStressVector, data.ElasticStressVectorPrincipalValues );
			CalculateEquivalentMeasure( data );
			UpdateDamageIndicator( data );
			CalculateDamage( data );
			CalculateStress( data, stressPerturbation );

			// fill the numerical tangent operator
			noalias( stressPerturbation ) -= stressVector;
			for(int i = 0; i < 3; i++)
				constitutiveMatrix(i, j) = stressPerturbation(i) / delta_strain;

			// restore internal variables
			mK = save_k; 
		}
	}

} /* namespace Kratos.*/
