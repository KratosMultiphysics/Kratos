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

#include "conv_diff_anisotropic_3d_law.h"
#include "multiscale_application_variables.h"

namespace Kratos
{

	ConvDiffAnisotropic3DLaw::ConvDiffAnisotropic3DLaw()
		: ConstitutiveLaw()
		, m_initialized(false)
		, m_init_gradT()
	{
	}

	ConstitutiveLaw::Pointer ConvDiffAnisotropic3DLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer(new ConvDiffAnisotropic3DLaw());
	}

	ConvDiffAnisotropic3DLaw::SizeType ConvDiffAnisotropic3DLaw::WorkingSpaceDimension()
	{
		return 3;
	}

	ConvDiffAnisotropic3DLaw::SizeType ConvDiffAnisotropic3DLaw::GetStrainSize()
	{
		return 3;
	}

	bool ConvDiffAnisotropic3DLaw::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffAnisotropic3DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		if (rThisVariable == INITIAL_TEMP_GRAD) return true;
		return false;
	}

	bool ConvDiffAnisotropic3DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffAnisotropic3DLaw::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
	{
		return false;
	}

	bool ConvDiffAnisotropic3DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	double& ConvDiffAnisotropic3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		return rValue;
	}

	Vector& ConvDiffAnisotropic3DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
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
			if (rValue.size() != 6)
				rValue.resize(6);
			rValue(0) = mStressVector(0); // / 1.0e6; //[W/mm^2]
			rValue(1) = mStressVector(1); // / 1.0e6;
			rValue(2) = mStressVector(2); // / 1.0e6;
			rValue(3) = 0.0;
			rValue(4) = 0.0;
			rValue(5) = 0.0;
		}
		return rValue;
	}

	Matrix& ConvDiffAnisotropic3DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 2 > & ConvDiffAnisotropic3DLaw::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & ConvDiffAnisotropic3DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	void ConvDiffAnisotropic3DLaw::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic3DLaw::SetValue(const Variable<Vector >& rVariable,
		const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == INITIAL_TEMP_GRAD) {
			if (rValue.size() == m_init_gradT.size())
				noalias(m_init_gradT) = rValue;
		}
		if (rVariable == FLUX_RVE)
		{
			if (rValue.size() == mStressVector.size())
				noalias(mStressVector) = rValue;
		}
	}

	void ConvDiffAnisotropic3DLaw::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic3DLaw::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
		const array_1d<double, 2 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic3DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool ConvDiffAnisotropic3DLaw::ValidateInput(const Properties& rMaterialProperties)
	{
		if (!rMaterialProperties.Has(CONDUCTIVITY)) return false;
		return true;
	}

	ConvDiffAnisotropic3DLaw::StrainMeasure ConvDiffAnisotropic3DLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	ConvDiffAnisotropic3DLaw::StressMeasure ConvDiffAnisotropic3DLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool ConvDiffAnisotropic3DLaw::IsIncremental()
	{
		return false;
	}

	void ConvDiffAnisotropic3DLaw::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (!m_initialized)
		{
			m_init_gradT = ZeroVector(this->GetStrainSize());
			mStressVector = ZeroVector(this->GetStrainSize());
			m_initialized = true;
		}
	}

	void ConvDiffAnisotropic3DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic3DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic3DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic3DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic3DLaw::CalculateMaterialResponsePK1(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void ConvDiffAnisotropic3DLaw::CalculateMaterialResponsePK2(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void ConvDiffAnisotropic3DLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void ConvDiffAnisotropic3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
	{
		// get some references
		const Properties& props = rValues.GetMaterialProperties();
		Vector& strainVector = rValues.GetStrainVector();
		Vector& stressVector = rValues.GetStressVector();
		Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
		Flags& Options = rValues.GetOptions();
		bool compute_constitutive_tensor = Options.Is(COMPUTE_CONSTITUTIVE_TENSOR);
		bool compute_stress = Options.Is(COMPUTE_STRESS) || compute_constitutive_tensor;

		SizeType size = GetStrainSize();
		if (compute_stress)
			if (stressVector.size() != size)
				stressVector.resize(size, false);
		if (compute_constitutive_tensor)
			if (constitutiveMatrix.size1() != size || constitutiveMatrix.size2() != size)
				constitutiveMatrix.resize(size, size, false);

		if (compute_stress)
		{
			CalculateStress(props, strainVector, stressVector);
		}

		if (compute_constitutive_tensor)
		{
			CalculateConstitutiveMatrix(props, strainVector, stressVector, constitutiveMatrix);
		}
	}

	void ConvDiffAnisotropic3DLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void ConvDiffAnisotropic3DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void ConvDiffAnisotropic3DLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void ConvDiffAnisotropic3DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
	{
	}

	void ConvDiffAnisotropic3DLaw::ResetMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (m_initialized)
		{
			m_init_gradT = ZeroVector(this->GetStrainSize());
		}
	}

	void ConvDiffAnisotropic3DLaw::CalculateStress(const Properties& props,
		const Vector& strainVector,
		Vector& stressVector)
	{
		double CONDUCTIVITY_11_11 = props[CONDUCTIVITY_1111];
		double CONDUCTIVITY_11_22 = props[CONDUCTIVITY_1122];
		double CONDUCTIVITY_11_33 = props[CONDUCTIVITY_1133];
		double CONDUCTIVITY_22_11 = props[CONDUCTIVITY_2211];
		double CONDUCTIVITY_22_22 = props[CONDUCTIVITY_2222];
		double CONDUCTIVITY_22_33 = props[CONDUCTIVITY_2233];
		double CONDUCTIVITY_33_11 = props[CONDUCTIVITY_3311];
		double CONDUCTIVITY_33_22 = props[CONDUCTIVITY_3322];
		double CONDUCTIVITY_33_33 = props[CONDUCTIVITY_3333];

		SizeType strain_size = strainVector.size();

		Matrix constitutiveMatrix(strain_size, strain_size, false);
		constitutiveMatrix(0, 0) = CONDUCTIVITY_11_11;
		constitutiveMatrix(0, 1) = CONDUCTIVITY_11_22;
		constitutiveMatrix(0, 2) = CONDUCTIVITY_11_33;
		constitutiveMatrix(1, 0) = CONDUCTIVITY_22_11;
		constitutiveMatrix(1, 1) = CONDUCTIVITY_22_22;
		constitutiveMatrix(1, 2) = CONDUCTIVITY_22_33;
		constitutiveMatrix(2, 0) = CONDUCTIVITY_33_11;
		constitutiveMatrix(2, 1) = CONDUCTIVITY_33_22;
		constitutiveMatrix(2, 2) = CONDUCTIVITY_33_33;

		mStressVector.clear();
		mStressVector = prod(constitutiveMatrix, strainVector - m_init_gradT);

		// mStressVector(0) = conductivity*(strainVector(0) - m_init_gradT(0));
		// mStressVector(1) = conductivity*(strainVector(1) - m_init_gradT(1));
		// mStressVector(2) = conductivity*(strainVector(2) - m_init_gradT(2));

		noalias(stressVector) = mStressVector;
	}

	void ConvDiffAnisotropic3DLaw::CalculateConstitutiveMatrix(const Properties& props,
		const Vector& strainVector,
		const Vector& stressVector,
		Matrix& constitutiveMatrix)
	{
		//size_t n = GetStrainSize();
		//double perturbation = 1.0E-6;
		//Vector perturbedStrainVector(n);
		//Vector stressPerturbation(n);

		//for (int j = 0; j < 3; j++)
		//{
		//	// FORWARD difference
		//	noalias(perturbedStrainVector) = strainVector;
		//	perturbedStrainVector(j) += perturbation;
		//	CalculateStress(props, perturbedStrainVector, stressPerturbation);

		//	// fill the numerical tangent operator
		//	noalias(stressPerturbation) -= stressVector;
		//	constitutiveMatrix(0, j) = stressPerturbation(0) / perturbation;
		//	constitutiveMatrix(1, j) = stressPerturbation(1) / perturbation;
		//	constitutiveMatrix(2, j) = stressPerturbation(2) / perturbation;
		//}

		double CONDUCTIVITY_11_11 = props[CONDUCTIVITY_1111];
		double CONDUCTIVITY_11_22 = props[CONDUCTIVITY_1122];
		double CONDUCTIVITY_11_33 = props[CONDUCTIVITY_1133];
		double CONDUCTIVITY_22_11 = props[CONDUCTIVITY_2211];
		double CONDUCTIVITY_22_22 = props[CONDUCTIVITY_2222];
		double CONDUCTIVITY_22_33 = props[CONDUCTIVITY_2233];
		double CONDUCTIVITY_33_11 = props[CONDUCTIVITY_3311];
		double CONDUCTIVITY_33_22 = props[CONDUCTIVITY_3322];
		double CONDUCTIVITY_33_33 = props[CONDUCTIVITY_3333];

		constitutiveMatrix.clear();
		constitutiveMatrix(0, 0) = CONDUCTIVITY_11_11;
		constitutiveMatrix(0, 1) = CONDUCTIVITY_11_22;
		constitutiveMatrix(0, 2) = CONDUCTIVITY_11_33;
		constitutiveMatrix(1, 0) = CONDUCTIVITY_22_11;
		constitutiveMatrix(1, 1) = CONDUCTIVITY_22_22;
		constitutiveMatrix(1, 2) = CONDUCTIVITY_22_33;
		constitutiveMatrix(2, 0) = CONDUCTIVITY_33_11;
		constitutiveMatrix(2, 1) = CONDUCTIVITY_33_22;
		constitutiveMatrix(2, 2) = CONDUCTIVITY_33_33;
	}

	void ConvDiffAnisotropic3DLaw::GetLawFeatures(Features& rFeatures)
	{
		rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( ISOTROPIC );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
		rFeatures.mStrainSize = GetStrainSize();
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	}

	int ConvDiffAnisotropic3DLaw::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			if(!rMaterialProperties.Has(CONDUCTIVITY)) {
				KRATOS_THROW_ERROR( std::logic_error, "ConvDiffAnisotropic3DLaw - missing CONDUCTIVITY", "");
			}
			return 0;

		KRATOS_CATCH("");
	}

} /* namespace Kratos.*/
