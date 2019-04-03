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

#include "conv_diff_anisotropic_2d_law.h"
#include "multiscale_application_variables.h"

namespace Kratos
{

	ConvDiffAnisotropic2DLaw::ConvDiffAnisotropic2DLaw()
		: ConstitutiveLaw()
		, m_initialized(false)
		, m_init_gradT()
		, mStressVector()
	{
	}

	ConstitutiveLaw::Pointer ConvDiffAnisotropic2DLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer(new ConvDiffAnisotropic2DLaw());
	}

	ConvDiffAnisotropic2DLaw::SizeType ConvDiffAnisotropic2DLaw::WorkingSpaceDimension()
	{
		return 2;
	}

	ConvDiffAnisotropic2DLaw::SizeType ConvDiffAnisotropic2DLaw::GetStrainSize()
	{
		return 2;
	}

	bool ConvDiffAnisotropic2DLaw::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffAnisotropic2DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		if(rThisVariable == INITIAL_TEMP_GRAD) return true;
		return false;
	}

	bool ConvDiffAnisotropic2DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffAnisotropic2DLaw::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
	{
		return false;
	}

	bool ConvDiffAnisotropic2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	double& ConvDiffAnisotropic2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		rValue = 0.0;
		return rValue;
	}

	Vector& ConvDiffAnisotropic2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
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
				rValue.resize(3, false);
			rValue(0) = mStressVector(0); // 1.0e6; //  //[W/mm^2]
			rValue(1) = mStressVector(1); // 1.0e6; //
			rValue(2) = 0.0;
			//std::stringstream ss;
			//ss << "HEAT_FLUX_RVE = " << rValue << ", " << std::endl;
			//std::cout << ss.str();
		}
		return rValue;
	}

	Matrix& ConvDiffAnisotropic2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 2 > & ConvDiffAnisotropic2DLaw::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & ConvDiffAnisotropic2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	void ConvDiffAnisotropic2DLaw::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic2DLaw::SetValue(const Variable<Vector >& rVariable,
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

	void ConvDiffAnisotropic2DLaw::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic2DLaw::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
		const array_1d<double, 2 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool ConvDiffAnisotropic2DLaw::ValidateInput(const Properties& rMaterialProperties)
	{
		if(!rMaterialProperties.Has(CONDUCTIVITY)) return false;
		if(!rMaterialProperties.Has(THICKNESS)) return false;
		return true;
	}

	ConvDiffAnisotropic2DLaw::StrainMeasure ConvDiffAnisotropic2DLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	ConvDiffAnisotropic2DLaw::StressMeasure ConvDiffAnisotropic2DLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool ConvDiffAnisotropic2DLaw::IsIncremental()
	{
		return false;
	}

	void ConvDiffAnisotropic2DLaw::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if(!m_initialized)
		{
			m_init_gradT = ZeroVector(this->GetStrainSize());
			mStressVector = ZeroVector(this->GetStrainSize());
			m_initialized = true;
		}
	}

	void ConvDiffAnisotropic2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffAnisotropic2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void ConvDiffAnisotropic2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void ConvDiffAnisotropic2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void ConvDiffAnisotropic2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
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

	void ConvDiffAnisotropic2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void ConvDiffAnisotropic2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void ConvDiffAnisotropic2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void ConvDiffAnisotropic2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{
	}

	void ConvDiffAnisotropic2DLaw::ResetMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (m_initialized)
		{
			m_init_gradT = ZeroVector(this->GetStrainSize());
			m_initialized = false;
		}
	}

	void ConvDiffAnisotropic2DLaw::CalculateStress(const Properties& props,
												   const Vector& strainVector,
												   Vector& stressVector)
	{
		double CONDUCTIVITY_11_11 = props[CONDUCTIVITY_1111];
		double CONDUCTIVITY_11_22 = props[CONDUCTIVITY_1122];
		double CONDUCTIVITY_22_11 = props[CONDUCTIVITY_2211];
		double CONDUCTIVITY_22_22 = props[CONDUCTIVITY_2222];

		SizeType strain_size = strainVector.size();

		Matrix constitutiveMatrix(strain_size, strain_size, false);
		constitutiveMatrix(0, 0) = CONDUCTIVITY_11_11;
		constitutiveMatrix(0, 1) = CONDUCTIVITY_11_22;
		constitutiveMatrix(1, 0) = CONDUCTIVITY_22_11;
		constitutiveMatrix(1, 1) = CONDUCTIVITY_22_22;

		mStressVector.clear();
		mStressVector = prod(constitutiveMatrix, strainVector - m_init_gradT);

		// mStressVector(0) = conductivity*(strainVector(0) - m_init_gradT(0));
		// mStressVector(1) = conductivity*(strainVector(1) - m_init_gradT(1));

		noalias(stressVector) = mStressVector;
	}

	void ConvDiffAnisotropic2DLaw::CalculateConstitutiveMatrix(const Properties& props,
															   const Vector& strainVector,
															   const Vector& stressVector,
															   Matrix& constitutiveMatrix)
	{
		//size_t n = GetStrainSize();
		//double perturbation = 1.0E-10;
		//Vector perturbedStrainVector(n);
		//Vector stressPerturbation(n);

		//for (int j = 0; j < 2; j++)
		//{
		//	// FORWARD difference
		//	noalias(perturbedStrainVector) = strainVector;
		//	perturbedStrainVector(j) += perturbation;
		//	CalculateStress(props, perturbedStrainVector, stressPerturbation);

		//	// fill the numerical tangent operator
		//	noalias(stressPerturbation) -= stressVector;
		//	constitutiveMatrix(0, j) = stressPerturbation(0) / perturbation;
		//	constitutiveMatrix(1, j) = stressPerturbation(1) / perturbation;
		//}
		double CONDUCTIVITY_11_11 = props[CONDUCTIVITY_1111];
		double CONDUCTIVITY_11_22 = props[CONDUCTIVITY_1122];
		double CONDUCTIVITY_22_11 = props[CONDUCTIVITY_2211];
		double CONDUCTIVITY_22_22 = props[CONDUCTIVITY_2222];

		constitutiveMatrix.clear();
		constitutiveMatrix(0, 0) = CONDUCTIVITY_11_11;
		constitutiveMatrix(0, 1) = CONDUCTIVITY_11_22;
		constitutiveMatrix(1, 0) = CONDUCTIVITY_22_11;
		constitutiveMatrix(1, 1) = CONDUCTIVITY_22_22;
		//std::stringstream ss;
		//ss << "PlaneStress - constitutiveMatrix = " << constitutiveMatrix << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffAnisotropic2DLaw::GetLawFeatures(Features& rFeatures)
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

	int ConvDiffAnisotropic2DLaw::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			if(!rMaterialProperties.Has(CONDUCTIVITY)) {
				KRATOS_THROW_ERROR( std::logic_error, "ConvDiffAnisotropic2DLaw - missing CONDUCTIVITY", "");
			}
			return 0;

		KRATOS_CATCH("");
	}

} /* namespace Kratos.*/
