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

#include "linear_elastic_thick_shell_law.h"
#include "multiscale_application.h"
#include "custom_utilities/math_helpers.h"
#include "includes/variables.h"

namespace Kratos
{

	LinearElasticThickShellLaw::LinearElasticThickShellLaw()
		: ConstitutiveLaw()
		, mDensity(0.0)
		, m_initial_strain(8,0.0)
	{
	}

	ConstitutiveLaw::Pointer LinearElasticThickShellLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer( new LinearElasticThickShellLaw() );
	}

	LinearElasticThickShellLaw::SizeType LinearElasticThickShellLaw::WorkingSpaceDimension()
	{
		return 3;
	}

	LinearElasticThickShellLaw::SizeType LinearElasticThickShellLaw::GetStrainSize()
	{
		return 8;
	}

	bool LinearElasticThickShellLaw::Has(const Variable<int>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticThickShellLaw::Has(const Variable<double>& rThisVariable)
	{
		if(rThisVariable == DENSITY)
			return true;
		return false;
	}

	bool LinearElasticThickShellLaw::Has(const Variable<Vector>& rThisVariable)
	{
		if(rThisVariable == INITIAL_STRAIN)
			return true;
		return false;
	}

	bool LinearElasticThickShellLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool LinearElasticThickShellLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	bool LinearElasticThickShellLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

	int& LinearElasticThickShellLaw::GetValue(
		const Variable<int>& rThisVariable,
		int& rValue)
	{
		rValue = 0;
		return rValue;
	}

	double& LinearElasticThickShellLaw::GetValue(
		const Variable<double>& rThisVariable,
		double& rValue)
	{
		rValue = 0.0;
		if(rThisVariable == DENSITY)
			rValue = mDensity;
		return rValue;
	}

	Vector& LinearElasticThickShellLaw::GetValue(
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

	Matrix& LinearElasticThickShellLaw::GetValue(
		const Variable<Matrix>& rThisVariable,
		Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & LinearElasticThickShellLaw::GetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & LinearElasticThickShellLaw::GetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void LinearElasticThickShellLaw::SetValue(
		const Variable<int>& rVariable,
		const int& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::SetValue(
		const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::SetValue(
		const Variable<Vector >& rVariable,
		const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable == INITIAL_STRAIN) {
			if(rValue.size() == m_initial_strain.size())
				noalias(m_initial_strain) = rValue;
		}
	}

	void LinearElasticThickShellLaw::SetValue(
		const Variable<Matrix >& rVariable,
		const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::SetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::SetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool LinearElasticThickShellLaw::ValidateInput(const Properties& rMaterialProperties)
	{
		if( !rMaterialProperties.Has(YOUNG_MODULUS) ) return false;
		if( !rMaterialProperties.Has(POISSON_RATIO) ) return false;
		if( !rMaterialProperties.Has(THICKNESS) ) return false;
		return true;
	}

	LinearElasticThickShellLaw::StrainMeasure LinearElasticThickShellLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	LinearElasticThickShellLaw::StressMeasure LinearElasticThickShellLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool LinearElasticThickShellLaw::IsIncremental()
	{
		return false;
	}

	void LinearElasticThickShellLaw::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		mDensity = 0.0;
		if(rMaterialProperties.Has(DENSITY))
			mDensity = rMaterialProperties[DENSITY];
		m_initial_strain.clear();
	}

	void LinearElasticThickShellLaw::InitializeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::FinalizeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::InitializeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::FinalizeNonLinearIteration(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void LinearElasticThickShellLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void LinearElasticThickShellLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void LinearElasticThickShellLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void LinearElasticThickShellLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
	{
		const Properties&   props = rValues.GetMaterialProperties();

		const Vector& strain_vector       = rValues.GetStrainVector();
		Vector&       stress_vector       = rValues.GetStressVector();
		Matrix&       constitutive_matrix  = rValues.GetConstitutiveMatrix();

		if(stress_vector.size() != 8)
			stress_vector.resize(8, false);

		if(constitutive_matrix.size1() != 8 || constitutive_matrix.size2() != 8)
			constitutive_matrix.resize(8,8,false);

		constitutive_matrix.clear();

		double E  = props[YOUNG_MODULUS];
		double nu = props[POISSON_RATIO];
		double H  = props[THICKNESS];
		double K  = 5.0/6.0;

		double c1 = E / (1.0 - nu*nu);
		double c2 = c1 * nu;
		double G  = E/(2.0*(nu+1.0));

		double I   = H*H*H/12.0;
		double GKH = G*K*H;

		constitutive_matrix(0,0) = c1*H;
		constitutive_matrix(0,1) = c2*H;
		constitutive_matrix(1,0) = constitutive_matrix(0,1);
		constitutive_matrix(1,1) = c1*H;
		constitutive_matrix(2,2) = G*H;
		constitutive_matrix(3,3) = c1*I;
		constitutive_matrix(3,4) = c2*I;
		constitutive_matrix(4,3) = constitutive_matrix(3,4);
		constitutive_matrix(4,4) = c1*I;
		constitutive_matrix(5,5) = G*I;
		constitutive_matrix(6,6) = GKH;
		constitutive_matrix(7,7) = GKH;

		stress_vector(0) = c1*H*(strain_vector(0)-m_initial_strain(0)) + c2*H*(strain_vector(1)-m_initial_strain(1));
		stress_vector(1) = c2*H*(strain_vector(0)-m_initial_strain(0)) + c1*H*(strain_vector(1)-m_initial_strain(1));
		stress_vector(2) = G*H*(strain_vector(2)-m_initial_strain(2));
		stress_vector(3) = c1*I*(strain_vector(3)-m_initial_strain(3)) + c2*I*(strain_vector(4)-m_initial_strain(4));
		stress_vector(4) = c2*I*(strain_vector(3)-m_initial_strain(3)) + c1*I*(strain_vector(4)-m_initial_strain(4));
		stress_vector(5) = G*I*(strain_vector(5)-m_initial_strain(5));
		stress_vector(6) = GKH*(strain_vector(6)-m_initial_strain(6));
		stress_vector(7) = GKH*(strain_vector(7)-m_initial_strain(7));
	}

	void LinearElasticThickShellLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void LinearElasticThickShellLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void LinearElasticThickShellLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void LinearElasticThickShellLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{
	}

	void LinearElasticThickShellLaw::ResetMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
	}

	void LinearElasticThickShellLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mOptions.Set( ISOTROPIC );

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}

	int LinearElasticThickShellLaw::Check(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if( !rMaterialProperties.Has(YOUNG_MODULUS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS", "");

		if( !rMaterialProperties.Has(POISSON_RATIO) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: POISSON_RATIO", "");

		if( !rMaterialProperties.Has(THICKNESS) )
			KRATOS_THROW_ERROR(std::logic_error, "Missing variable: THICKNESS", "");

		return 0;

		KRATOS_CATCH("");
	}

	void LinearElasticThickShellLaw::CalculateMaterialResponse(const Vector& StrainVector,
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
