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

#include "conv_diff_plane_stress_2d_law.h"
#include "multiscale_application_variables.h"

namespace Kratos
{

	ConvDiffPlaneStress2DLaw::ConvDiffPlaneStress2DLaw()
		: ConstitutiveLaw()
		, m_initialized(false)
		, m_init_gradT()
		, mStressVector()
	{
	}

	ConstitutiveLaw::Pointer ConvDiffPlaneStress2DLaw::Clone() const
	{
		return ConstitutiveLaw::Pointer(new ConvDiffPlaneStress2DLaw());
	}

	ConvDiffPlaneStress2DLaw::SizeType ConvDiffPlaneStress2DLaw::WorkingSpaceDimension()
	{
		return 2;
	}

	ConvDiffPlaneStress2DLaw::SizeType ConvDiffPlaneStress2DLaw::GetStrainSize()
	{
		return 2;
	}

	bool ConvDiffPlaneStress2DLaw::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffPlaneStress2DLaw::Has(const Variable<Vector>& rThisVariable)
	{
		if(rThisVariable == INITIAL_TEMP_GRAD) return true;
		return false;
	}

	bool ConvDiffPlaneStress2DLaw::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool ConvDiffPlaneStress2DLaw::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
	{
		return false;
	}

	bool ConvDiffPlaneStress2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	double& ConvDiffPlaneStress2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		return rValue;
	}

	Vector& ConvDiffPlaneStress2DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
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
			/*std::stringstream ss;
			ss << "HEAT_FLUX_RVE = " << rValue << ", " << std::endl;
			std::cout << ss.str();*/
		}
		return rValue;
	}

	Matrix& ConvDiffPlaneStress2DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 2 > & ConvDiffPlaneStress2DLaw::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & ConvDiffPlaneStress2DLaw::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	void ConvDiffPlaneStress2DLaw::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffPlaneStress2DLaw::SetValue(const Variable<Vector >& rVariable,
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

	void ConvDiffPlaneStress2DLaw::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffPlaneStress2DLaw::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
		const array_1d<double, 2 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffPlaneStress2DLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool ConvDiffPlaneStress2DLaw::ValidateInput(const Properties& rMaterialProperties)
	{
		if(!rMaterialProperties.Has(CONDUCTIVITY)) return false;
		if(!rMaterialProperties.Has(THICKNESS)) return false;
		return true;
	}

	ConvDiffPlaneStress2DLaw::StrainMeasure ConvDiffPlaneStress2DLaw::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	ConvDiffPlaneStress2DLaw::StressMeasure ConvDiffPlaneStress2DLaw::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool ConvDiffPlaneStress2DLaw::IsIncremental()
	{
		return false;
	}

	void ConvDiffPlaneStress2DLaw::InitializeMaterial(
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

	void ConvDiffPlaneStress2DLaw::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffPlaneStress2DLaw::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffPlaneStress2DLaw::InitializeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffPlaneStress2DLaw::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffPlaneStress2DLaw::CalculateMaterialResponsePK1 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void ConvDiffPlaneStress2DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void ConvDiffPlaneStress2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
	{
		CalculateMaterialResponseCauchy( rValues );
	}

	void ConvDiffPlaneStress2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
	{
		// get some references
		const Properties& props = rValues.GetMaterialProperties();
		Vector& strainVector = rValues.GetStrainVector();
		Vector& stressVector = rValues.GetStressVector();
		Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
		Flags& Options = rValues.GetOptions();
		bool compute_constitutive_tensor = Options.Is(COMPUTE_CONSTITUTIVE_TENSOR);
		bool compute_stress = Options.Is(COMPUTE_STRESS) || compute_constitutive_tensor;

		noalias(strainVector) -= m_init_gradT;

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

	void ConvDiffPlaneStress2DLaw::FinalizeMaterialResponsePK1 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void ConvDiffPlaneStress2DLaw::FinalizeMaterialResponsePK2 (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void ConvDiffPlaneStress2DLaw::FinalizeMaterialResponseKirchhoff (Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy( rValues );
	}

	void ConvDiffPlaneStress2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
	{
	}

	void ConvDiffPlaneStress2DLaw::ResetMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (m_initialized)
		{
			m_init_gradT = ZeroVector(this->GetStrainSize());
			m_initialized = false;
		}
	}

	void ConvDiffPlaneStress2DLaw::CalculateStress(const Properties& props,
												   const Vector& strainVector,
												   Vector& stressVector)
	{
		mStressVector.clear();
		double conductivity = props[CONDUCTIVITY];
		mStressVector(0) = conductivity*(strainVector(0));
		mStressVector(1) = conductivity*(strainVector(1));
		noalias(stressVector) = mStressVector;

		//std::stringstream ss;
		//ss << "PlaneStress - m_init_gradT = " << m_init_gradT << ", " << std::endl;
		//if (strainVector(1) > 0.01)
		//	ss << "PlaneStress - strainVector = " << strainVector << ", " << std::endl;
		////if (stressVector(1) > 0.01)
		////	ss << "PlaneStress - stressVector = " << stressVector << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffPlaneStress2DLaw::CalculateConstitutiveMatrix(const Properties& props,
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
		//std::stringstream ss;
		//ss << "Perturbed - constitutiveMatrix = " << constitutiveMatrix << ", " << std::endl;
		//std::cout << ss.str();

		double conductivity = props[CONDUCTIVITY];
		constitutiveMatrix.clear();
		constitutiveMatrix(0, 0) = conductivity;
		constitutiveMatrix(1, 1) = conductivity;
		//std::stringstream ss;
		//ss << "NO Perturbed - constitutiveMatrix = " << constitutiveMatrix << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffPlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
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

	int ConvDiffPlaneStress2DLaw::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			if(!rMaterialProperties.Has(CONDUCTIVITY)) {
				KRATOS_THROW_ERROR( std::logic_error, "ConvDiffPlaneStress2DLaw - missing CONDUCTIVITY", "");
			}
			return 0;

		KRATOS_CATCH("");
	}

} /* namespace Kratos.*/
