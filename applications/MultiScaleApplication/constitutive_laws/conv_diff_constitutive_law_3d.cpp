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

#include "conv_diff_constitutive_law_3d.h"
#include "multiscale_application_variables.h"

namespace Kratos
{

	ConvDiffConstitutiveLaw3D::ConvDiffConstitutiveLaw3D()
		: ConstitutiveLaw()
		, m_initialized(false)
		, m_init_gradT()
	{
	}

	ConstitutiveLaw::Pointer ConvDiffConstitutiveLaw3D::Clone() const
	{
		return ConstitutiveLaw::Pointer(new ConvDiffConstitutiveLaw3D());
	}

	ConvDiffConstitutiveLaw3D::SizeType ConvDiffConstitutiveLaw3D::WorkingSpaceDimension()
	{
		return 3;
	}

	ConvDiffConstitutiveLaw3D::SizeType ConvDiffConstitutiveLaw3D::GetStrainSize()
	{
		return 3;
	}

	bool ConvDiffConstitutiveLaw3D::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffConstitutiveLaw3D::Has(const Variable<Vector>& rThisVariable)
	{
		if (rThisVariable == INITIAL_TEMP_GRAD) return true;
		return false;
	}

	bool ConvDiffConstitutiveLaw3D::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffConstitutiveLaw3D::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
	{
		return false;
	}

	bool ConvDiffConstitutiveLaw3D::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	double& ConvDiffConstitutiveLaw3D::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		return rValue;
	}

	Vector& ConvDiffConstitutiveLaw3D::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		//std::stringstream ss;
		//ss << "ConvDiffConstitutiveLaw3D::GetValue" << std::endl;
		if (rThisVariable == INITIAL_TEMP_GRAD) {
			if (rValue.size() != m_init_gradT.size())
				rValue.resize(m_init_gradT.size());
			noalias(rValue) = m_init_gradT;
		}
		if (rThisVariable == FLUX_RVE || rThisVariable == HEAT_FLUX_RVE) {
			if (rValue.size() != mStressVector.size())
				rValue.resize(mStressVector.size());
			noalias(rValue) = mStressVector;
		}
		//if (rThisVariable == HEAT_FLUX_RVE) { //For Output
		//	if (rValue.size() != 6)
		//		rValue.resize(6);
		//	rValue(0) = mStressVector(0); // / 1.0e6; //[W/mm^2]
		//	rValue(1) = mStressVector(1); // / 1.0e6;
		//	rValue(2) = mStressVector(2); // / 1.0e6;
		//	rValue(3) = 0.0;
		//	rValue(4) = 0.0;
		//	rValue(5) = 0.0;

		//	//ss << "HEAT_FLUX_RVE = " << rValue << ", " << std::endl;
		//	//std::cout << ss.str();
		//}
		return rValue;
	}

	Matrix& ConvDiffConstitutiveLaw3D::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 2 > & ConvDiffConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & ConvDiffConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	void ConvDiffConstitutiveLaw3D::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffConstitutiveLaw3D::SetValue(const Variable<Vector >& rVariable,
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

	void ConvDiffConstitutiveLaw3D::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
		const array_1d<double, 2 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool ConvDiffConstitutiveLaw3D::ValidateInput(const Properties& rMaterialProperties)
	{
		if (!rMaterialProperties.Has(CONDUCTIVITY)) return false;
		return true;
	}

	ConvDiffConstitutiveLaw3D::StrainMeasure ConvDiffConstitutiveLaw3D::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	ConvDiffConstitutiveLaw3D::StressMeasure ConvDiffConstitutiveLaw3D::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool ConvDiffConstitutiveLaw3D::IsIncremental()
	{
		return false;
	}

	void ConvDiffConstitutiveLaw3D::InitializeMaterial(
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

	void ConvDiffConstitutiveLaw3D::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffConstitutiveLaw3D::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffConstitutiveLaw3D::InitializeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffConstitutiveLaw3D::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffConstitutiveLaw3D::CalculateMaterialResponsePK1(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void ConvDiffConstitutiveLaw3D::CalculateMaterialResponsePK2(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void ConvDiffConstitutiveLaw3D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void ConvDiffConstitutiveLaw3D::CalculateMaterialResponseCauchy(Parameters& rValues)
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

	void ConvDiffConstitutiveLaw3D::FinalizeMaterialResponsePK1(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void ConvDiffConstitutiveLaw3D::FinalizeMaterialResponsePK2(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void ConvDiffConstitutiveLaw3D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void ConvDiffConstitutiveLaw3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
	{
	}

	void ConvDiffConstitutiveLaw3D::ResetMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (m_initialized)
		{
			m_init_gradT = ZeroVector(this->GetStrainSize());
		}
	}

	void ConvDiffConstitutiveLaw3D::CalculateStress(const Properties& props,
		const Vector& strainVector,
		Vector& stressVector)
	{
		double conductivity = props[CONDUCTIVITY];
		mStressVector(0) = conductivity*(strainVector(0) - m_init_gradT(0));
		mStressVector(1) = conductivity*(strainVector(1) - m_init_gradT(1));
		mStressVector(2) = conductivity*(strainVector(2) - m_init_gradT(2));
		noalias(stressVector) = mStressVector;
		//std::stringstream ss;
		//ss << "stressVector = " << stressVector << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffConstitutiveLaw3D::CalculateConstitutiveMatrix(const Properties& props,
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

		double conductivity = props[CONDUCTIVITY];
		constitutiveMatrix.clear();
		constitutiveMatrix(0, 0) = conductivity;
		constitutiveMatrix(1, 1) = conductivity;
		constitutiveMatrix(2, 2) = conductivity;
	}

	void ConvDiffConstitutiveLaw3D::GetLawFeatures(Features& rFeatures)
	{
		rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
		rFeatures.mOptions.Set( ISOTROPIC );
		rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
		rFeatures.mStrainSize = GetStrainSize();
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	}

	int ConvDiffConstitutiveLaw3D::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			if(!rMaterialProperties.Has(CONDUCTIVITY)) {
				KRATOS_THROW_ERROR( std::logic_error, "ConvDiffConstitutiveLaw3D - missing CONDUCTIVITY", "");
			}
			return 0;

		KRATOS_CATCH("");
	}

} /* namespace Kratos.*/
