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

#include "linear_elastic_no_tension_plane_stress_2d_law.h"
#include "multiscale_application.h"
#include "custom_utilities/math_helpers.h"
#include "includes/variables.h"

namespace Kratos
{

	LinearElasticNoTensionPlaneStress2DLaw::LinearElasticNoTensionPlaneStress2DLaw()
		: ConstitutiveLaw()
	{
	}

	ConstitutiveLaw::Pointer LinearElasticNoTensionPlaneStress2DLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer( new LinearElasticNoTensionPlaneStress2DLaw() );
	}

	LinearElasticNoTensionPlaneStress2DLaw::SizeType LinearElasticNoTensionPlaneStress2DLaw::WorkingSpaceDimension()
	{
		return 2;
	}

	LinearElasticNoTensionPlaneStress2DLaw::SizeType LinearElasticNoTensionPlaneStress2DLaw::GetStrainSize()
	{
		return 3;
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::Has(const Variable<int>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

	int& LinearElasticNoTensionPlaneStress2DLaw::GetValue(
		const Variable<int>& rThisVariable,
		int& rValue)
	{
		rValue = 0;
		return rValue;
	}

	double& LinearElasticNoTensionPlaneStress2DLaw::GetValue(
		const Variable<double>& rThisVariable,
		double& rValue)
	{
		rValue = 0.0;
		return rValue;
	}

	Vector& LinearElasticNoTensionPlaneStress2DLaw::GetValue(
		const Variable<Vector>& rThisVariable,
		Vector& rValue)
	{
		return rValue;
	}

	Matrix& LinearElasticNoTensionPlaneStress2DLaw::GetValue(
		const Variable<Matrix>& rThisVariable,
		Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & LinearElasticNoTensionPlaneStress2DLaw::GetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & LinearElasticNoTensionPlaneStress2DLaw::GetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void LinearElasticNoTensionPlaneStress2DLaw::SetValue(
		const Variable<int>& rVariable,
		const int& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::SetValue(
		const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::SetValue(
		const Variable<Vector >& rVariable,
		const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::SetValue(
		const Variable<Matrix >& rVariable,
		const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::SetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::SetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::ValidateInput(const Properties& rMaterialProperties)
	{
		if( !rMaterialProperties.Has(YOUNG_MODULUS) ) return false;
		return true;
	}

	LinearElasticNoTensionPlaneStress2DLaw::StrainMeasure LinearElasticNoTensionPlaneStress2DLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	LinearElasticNoTensionPlaneStress2DLaw::StressMeasure LinearElasticNoTensionPlaneStress2DLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool LinearElasticNoTensionPlaneStress2DLaw::IsIncremental()
	{
		return false;
	}

	void LinearElasticNoTensionPlaneStress2DLaw::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::InitializeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::FinalizeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::InitializeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::FinalizeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void LinearElasticNoTensionPlaneStress2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void LinearElasticNoTensionPlaneStress2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void LinearElasticNoTensionPlaneStress2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
	{
		const Properties&   props = rValues.GetMaterialProperties();

		const Vector& strain_vector       = rValues.GetStrainVector();
		Vector&       stress_vector       = rValues.GetStressVector();
		Matrix&       constitutive_matrix  = rValues.GetConstitutiveMatrix();

		if(stress_vector.size() != 3)
			stress_vector.resize(3, false);
		stress_vector.clear();

		if(constitutive_matrix.size1() != 3 || constitutive_matrix.size2() != 3)
			constitutive_matrix.resize(3,3,false);
		constitutive_matrix.clear();

		constitutive_matrix.clear();

		double E  = props[YOUNG_MODULUS];

		if(strain_vector(0) <= 0.0)
		{
			stress_vector(0) = E*strain_vector(0);
			constitutive_matrix(0,0) = E;
		}
		else
		{
			stress_vector(0) = 1.0e-5*E*strain_vector(0);
			constitutive_matrix(0,0) = 1.0e-5*E;
		}
		stress_vector(1) = E*strain_vector(1);
		constitutive_matrix(1,1) = E;
		stress_vector(2) = E/2.0*strain_vector(2);
		constitutive_matrix(2,2) = E/2.0;
	}

	void LinearElasticNoTensionPlaneStress2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void LinearElasticNoTensionPlaneStress2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void LinearElasticNoTensionPlaneStress2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void LinearElasticNoTensionPlaneStress2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::ResetMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
	}

	void LinearElasticNoTensionPlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
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

	int LinearElasticNoTensionPlaneStress2DLaw::Check(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if( !rMaterialProperties.Has(YOUNG_MODULUS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS", "");

		return 0;

		KRATOS_CATCH("");
	}

	void LinearElasticNoTensionPlaneStress2DLaw::CalculateMaterialResponse(const Vector& StrainVector,
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
		Matrix F(IdentityMatrix(3,3));
		parameters.SetDeterminantF(detF);
		parameters.SetDeformationGradientF(F);
		parameters.SetShapeFunctionsValues(rShapeFunctionsValues);
		this->CalculateMaterialResponseCauchy(parameters);
	}

} /* namespace Kratos.*/
