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
*   Last Modified by:    $Author:   Stefano Zaghi$
*   Date:                $Date:     30-10-2016$
*   Revision:            $Revision: 1.0$
*
* ***********************************************************/

#include <utility>

#include "includes/kratos_parameters.h"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include "interpolated_constitutive_law_3d.h"
#include "multiscale_application_variables.h"

#define SIGMA_SIGN(X) (X == 0.0 ? 1.0 : ( X > 0.0 ? 1.0 : -1.0 ))

namespace Kratos
{
	template <typename T>
	vector<size_t> sort_indexes(const vector<T> &v) { //only c++11

		// initialize original index locations
		vector<size_t> idx(v.size());
		for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

		// sort indexes based on comparing values in v
		sort(idx.begin(), idx.end(),
			[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

		return idx;
	}

	InterpolatedConstitutiveLaw3D::InterpolatedConstitutiveLaw3D(const RveMaterialDatabase::Pointer& pNewRveMaterialDatabase)
		: ConstitutiveLaw()
		, mpRveMaterialDatabase(RveMaterialDatabase::Pointer())
		, mInitialized(false)
		, mInitStrain()
		, mStressVector()
		, mC0()
		, mStrainMap()
		, mTagDamageMap()
		, mDamage(0.0)
		, m_r(0.0)
		, m_r_converged(0.0)
		, mDamageConverged(0.0)
		, mIsElastic(false)
		, mlch(0.0)
		, mlchRef(0.0)
		, mG0(0.0)
		, mGf_bar(0.0)
		, mAlpha(0.0)
	{
		//std::cout << "READY TO READ MATERIAL DATABASE\n";
		if (pNewRveMaterialDatabase == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - The input RveMaterialDatabase is NULL", "");
		mpRveMaterialDatabase = pNewRveMaterialDatabase;
	}

	ConstitutiveLaw::Pointer InterpolatedConstitutiveLaw3D::Clone() const
	{
		return ConstitutiveLaw::Pointer(new InterpolatedConstitutiveLaw3D(mpRveMaterialDatabase));
	}

	InterpolatedConstitutiveLaw3D::SizeType InterpolatedConstitutiveLaw3D::WorkingSpaceDimension()
	{
		return 3;
	}

	InterpolatedConstitutiveLaw3D::SizeType InterpolatedConstitutiveLaw3D::GetStrainSize()
	{
		return 6;
	}

	bool InterpolatedConstitutiveLaw3D::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool InterpolatedConstitutiveLaw3D::Has(const Variable<Vector>& rThisVariable)
	{
		if (rThisVariable == INITIAL_STRAIN) return true;
		return false;
	}

	bool InterpolatedConstitutiveLaw3D::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool InterpolatedConstitutiveLaw3D::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	bool InterpolatedConstitutiveLaw3D::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return false;
	}

	double& InterpolatedConstitutiveLaw3D::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == LCH_REF_RVE){
			rValue = mlchRef;
		}
		else if (rThisVariable == EQUIVALENT_DAMAGE)
		{
			rValue = mDamage;
		}
		return rValue;
	}

	Vector& InterpolatedConstitutiveLaw3D::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		if (rThisVariable == INITIAL_STRAIN) {
			if (rValue.size() != mInitStrain.size())
				rValue.resize(mInitStrain.size());
			noalias(rValue) = mInitStrain;
		}
		return rValue;
	}

	Matrix& InterpolatedConstitutiveLaw3D::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		if (rThisVariable == PREDICTED_STRESS_TENSOR)
			StressVectorToTensor(mStressVector, rValue);
		return rValue;
	}

	array_1d<double, 3 > & InterpolatedConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & InterpolatedConstitutiveLaw3D::GetValue(const Variable<array_1d<double, 6 > >& rVariable, array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void InterpolatedConstitutiveLaw3D::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == LCH_REF_RVE) {
			mlchRef = rValue;
		}
	}

	void InterpolatedConstitutiveLaw3D::SetValue(const Variable<Vector >& rVariable,
		const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == INITIAL_STRAIN) {
			if (rValue.size() == mInitStrain.size())
				noalias(mInitStrain) = rValue;
		}
	}

	void InterpolatedConstitutiveLaw3D::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw3D::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool InterpolatedConstitutiveLaw3D::ValidateInput(const Properties& rMaterialProperties)
	{
		if (!rMaterialProperties.Has(LCH_REF_RVE)) return false;
		return true;
	}

	InterpolatedConstitutiveLaw3D::StrainMeasure InterpolatedConstitutiveLaw3D::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	InterpolatedConstitutiveLaw3D::StressMeasure InterpolatedConstitutiveLaw3D::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool InterpolatedConstitutiveLaw3D::IsIncremental()
	{
		return false;
	}

	void InterpolatedConstitutiveLaw3D::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (!mInitialized)
		{
			m_r = 0.0;
			m_r_converged = m_r;
			mInitStrain = ZeroVector(this->GetStrainSize());
			mStressVector = ZeroVector(this->GetStrainSize());
			mDamage = 0.0;
			mDamageConverged = mDamage;
			mlch = rElementGeometry.Length();
			mlchRef = rMaterialProperties[LCH_REF_RVE];
			mC0 = mpRveMaterialDatabase->GetConstitutiveTensor();
			//mStrainMap = mpRveMaterialDatabase->RveElasticMap();
			//if (mStrainMap.empty()) //Check if strain_map it is not empty
			//{
			//	std::stringstream ss;
			//	ss << "RvePredictorCalculator: STRAIN MAP IS EMPTY" << std::endl;
			//	ss << "RvePredictorCalculator: Cannot calculate the material response without a micromodel history." << std::endl;
			//	std::cout << ss.str();
			//	KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - Cannot calculate the material response without a micromodel history.", "");
			//}
			//mTagDamageMap = mpRveMaterialDatabase->RveMaterialMap();
			//if (mTagDamageMap.empty()) //Check if strain_map it is not empty
			//{
			//	std::stringstream ss;
			//	ss << "RvePredictorCalculator: STRAIN MAP IS EMPTY" << std::endl;
			//	ss << "RvePredictorCalculator: Cannot calculate the material response without a micromodel history." << std::endl;
			//	std::cout << ss.str();
			//	KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - Cannot calculate the material response without a micromodel history.", "");
			//}
			mInitialized = true;
		}
	}

	void InterpolatedConstitutiveLaw3D::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw3D::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		m_r_converged = m_r;
		mDamageConverged = mDamage;
	}

	void InterpolatedConstitutiveLaw3D::InitializeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw3D::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw3D::CalculateMaterialResponsePK1(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw3D::CalculateMaterialResponsePK2(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw3D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw3D::CalculateMaterialResponseCauchy(Parameters& rValues)
	{
		// get some references
		const ProcessInfo& pinfo = rValues.GetProcessInfo();
		const Properties& props = rValues.GetMaterialProperties();
		Vector& strainVector = rValues.GetStrainVector();
		Vector& stressVector = rValues.GetStressVector();
		Matrix& constitutiveMatrix = rValues.GetConstitutiveMatrix();
		Flags& Options = rValues.GetOptions();
		bool compute_constitutive_tensor = Options.Is(COMPUTE_CONSTITUTIVE_TENSOR);
		bool compute_stress = Options.Is(COMPUTE_STRESS) || compute_constitutive_tensor;

		// Viscosity
		rValues.GetProcessInfo().Has(RVE_DAMAGE_SURFACE_FLAG) ? rValues.GetProcessInfo()[RVE_DAMAGE_SURFACE_FLAG] : 0;
		m_eta = props.Has(VISCOSITY) ? props[VISCOSITY] : 0.0;
		m_dTime = pinfo[DELTA_TIME];
		if (m_dTime > 0.0 && m_eta > 0.0)
		{
			m_rate_coeff_1 = m_eta / (m_eta + m_dTime);
			m_rate_coeff_2 = m_dTime / (m_eta + m_dTime);
		}
		else
		{
			m_rate_coeff_1 = 0.0;
			m_rate_coeff_2 = 1.0;
		}

		noalias(strainVector) -= mInitStrain;

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

	void InterpolatedConstitutiveLaw3D::FinalizeMaterialResponsePK1(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw3D::FinalizeMaterialResponsePK2(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw3D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw3D::FinalizeMaterialResponseCauchy(Parameters& rValues)
	{
	}

	void InterpolatedConstitutiveLaw3D::ResetMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (mInitialized)
		{
			m_r = 0.0;
			m_r_converged = 0.0;
			mInitStrain = ZeroVector(this->GetStrainSize());
			mDamage = 0.0;
			mDamageConverged = 0.0;
			mIsElastic = false;
			mlch = 0.0;
			mInitialized = false;
		}
	}

	void InterpolatedConstitutiveLaw3D::CalculateStress(const Properties& props,
		const Vector& strainVector,
		Vector& stressVector)
	{
		m_r = m_r_converged;

		size_t strain_size = strainVector.size();
		mStressVector.clear();

		//Initialize data for the interpolation
		Vector real_tag(strain_size - 1);
		Vector Theta(strain_size - 1);
		Vector eps(strain_size);
		double radius = 0.0;

		//Calculate phi & number of division
		const size_t m = mpRveMaterialDatabase->NumberOfDivision();
		//std::cout << "m: " << m << std::endl;
		double phi = Globals::Pi / m;
		//std::cout << "phi: " << phi << " rad" << std::endl;

		//1. Calculate the normalized Strain
		double real_radius = norm_2(strainVector);
		//std::cout << "real_radius: " << real_radius << std::endl;

		Theta = this->CalculateTheta(strainVector);
		//std::cout << "Theta: " << Theta << std::endl;

		real_tag = Theta / phi;
		//std::cout << "PredictElasticity -> real_tag: " << real_tag << std::endl;

		for (size_t i_tag = 0; i_tag < real_tag.size(); i_tag++)
		{
			if (abs(real_tag[i_tag]) > m)
			{
				double mult_0 = abs(floor(real_tag[i_tag] / m));
				if (real_tag[i_tag] > m)
					real_tag[i_tag] -= m*mult_0;
				else
					real_tag[i_tag] += m*mult_0;
			}
		}

		//LINEAR INTERPOLATION OVER TAG
		//1. Found the minimum value of the tag
		vector<int> min_tag(5);	vector<int> max_tag(5);
		for (size_t i_tag = 0; i_tag < real_tag.size(); i_tag++)
		{
			min_tag[i_tag] = floor(real_tag[i_tag]);
			max_tag[i_tag] = ceil(real_tag[i_tag]);
		}

		this->TagsLinearInterpolation(real_tag, min_tag, radius);

		//std::cout << "real_radius VS radius: " << real_radius << " VS " << radius << std::endl;
		//std::cout << "difference: " << real_radius - radius << std::endl;
		//KRATOS_WATCH(mIsElastic)
		if (real_radius < radius)
		{
			mIsElastic = true;
			mStressVector = prod(mC0, strainVector);
		}
		else
		{
			mIsElastic = false;
			//m_r = m_rate_coeff_1*m_r + m_rate_coeff_2*real_radius;
			this->StressPrediction(real_tag, min_tag, real_radius, Theta, strainVector, mStressVector);
		}

		noalias(stressVector) = mStressVector;
	}

	void InterpolatedConstitutiveLaw3D::TagsLinearInterpolation(const Vector& xy, const vector<int>& min_xy, double& z)
	{
		double m = mpRveMaterialDatabase->NumberOfDivision();

		if (xy[0] == min_xy[0] && xy[1] == min_xy[1] && xy[2] == min_xy[2] && xy[3] == min_xy[3] && xy[4] == min_xy[4])
		{
			mpRveMaterialDatabase->GetElasticRadius(xy,z);
		}
		else
		{
			// local_x and y defines where to estimate the value on the interpolated line,
			// it is 0 at the first point and 1 and the second point
			double local_x = (xy[0] - min_xy[0]);
			double local_y = (xy[1] - min_xy[1]);
			double local_z = (xy[2] - min_xy[2]);
			double local_k = (xy[3] - min_xy[3]);
			double local_l = (xy[4] - min_xy[4]);
			//std::cout << "local_x, local_y: " << local_x << ", " << local_y << std::endl;

			array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> p;

			size_t k = 0;
			for (size_t kk = 0; kk < 2; kk++)
			{
				size_t l = 0;
				for (size_t ll = 0; ll < 2; ll++)
				{
					size_t m = 0;
					for (size_t mm = 0; mm < 2; ll++)
					{
						size_t n = 0;
						for (size_t nn = 0; nn < 2; ll++)
						{
							size_t o = 0;
							for (size_t oo = 0; oo < 2; ll++)
							{
								vector<int> tag_i(5);
								tag_i[0] = min_xy[0] + kk;
								tag_i[1] = min_xy[1] + ll;
								tag_i[2] = min_xy[2] + mm;
								tag_i[3] = min_xy[3] + nn;
								tag_i[4] = min_xy[4] + oo;

								for (size_t iter = 0; iter < tag_i.size(); iter++)
								{
									if (abs(tag_i[iter]) > m)
									{
										double mult_0 = abs(floor(tag_i[iter] / m));
										if (tag_i[iter] > m)
											tag_i[iter] -= m*mult_0;
										else
											tag_i[iter] += m*mult_0;
									}
								}

								//std::cout << "tag_i: " << tag_i << std::endl;

								double i_z = 0.0;
								mpRveMaterialDatabase->GetElasticRadius(tag_i, i_z);

								p[k][l][m][n][o] = i_z;
								oo++;
							}
							n++;
						}
						m++;
					}
					l++;
				}
				k++;
			}
			z = this->pentaDimInterpolate(p, local_x, local_y, local_z, local_k, local_l);
		}
	}

	double InterpolatedConstitutiveLaw3D::Interpolate(array_1d<double, 2> p, double x) const {
		double x0 = 0.0;
		double x1 = 1.0;

		return (p[0] * (x1 - x) + p[1] * (x - x0)) / (x1 - x0);
	}

	double InterpolatedConstitutiveLaw3D::biDimInterpolate(array_1d<array_1d<double, 2>, 2> p, double x, double y) const {
		array_1d<double, 2> arr;
		arr[0] = Interpolate(p[0], y);
		arr[1] = Interpolate(p[1], y);
		return Interpolate(arr, x);
	}

	double InterpolatedConstitutiveLaw3D::triDimInterpolate(array_1d<array_1d<array_1d<double, 2>, 2>, 2> p, double x, double y, double z) const {
		array_1d<double, 2> arr;
		arr[0] = biDimInterpolate(p[0], y, z);
		arr[1] = biDimInterpolate(p[1], y, z);
		return Interpolate(arr, x);
	}

	double InterpolatedConstitutiveLaw3D::quadDimInterpolate(array_1d<array_1d<array_1d<array_1d<double, 2>, 2>, 2>, 2> p, double x, double y, double z, double h) const {
		array_1d<double, 2> arr;
		arr[0] = triDimInterpolate(p[0], y, z, h);
		arr[1] = triDimInterpolate(p[1], y, z, h);
		return Interpolate(arr, x);
	}

	double InterpolatedConstitutiveLaw3D::pentaDimInterpolate(array_1d<array_1d<array_1d<array_1d<array_1d<double, 2>, 2>, 2>, 2>, 2> p, double x, double y, double z, double h, double k) const {
		array_1d<double, 2> arr;
		arr[0] = quadDimInterpolate(p[0], y, z, h, k);
		arr[1] = quadDimInterpolate(p[1], y, z, h, k);
		return Interpolate(arr, x);
	}

	void InterpolatedConstitutiveLaw3D::StressPrediction(const Vector& real_tag, const vector<int>& min_tag, const double& real_radius, const Vector& Theta, const Vector& strain_vector, Vector& stress_vector)
	{
		size_t strain_size = strain_vector.size();
		double m = mpRveMaterialDatabase->NumberOfDivision();
		const size_t n_damage = mpRveMaterialDatabase->DamageSize();

		Vector InterpDamageVector(n_damage, 0.0);
		Vector InterpRadiusVector(n_damage, 0.0);
		Vector InterpStressVectorXX(n_damage, 0.0);
		Vector InterpStressVectorYY(n_damage, 0.0);
		Vector InterpStressVectorZZ(n_damage, 0.0);
		Vector InterpStressVectorXY(n_damage, 0.0);
		Vector InterpStressVectorYZ(n_damage, 0.0);
		Vector InterpStressVectorXZ(n_damage, 0.0);

		if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1] && real_tag[2] == min_tag[2] && real_tag[3] == min_tag[3] && real_tag[4] == min_tag[4])
		{
			Vector DamageVector(n_damage, 0.0);
			Vector RadiusVector(n_damage, 0.0);
			Matrix StressMatrix(strain_size, n_damage, 0.0);

			mpRveMaterialDatabase->GetMaterialValues(min_tag, DamageVector, RadiusVector, StressMatrix);
			noalias(InterpDamageVector) = DamageVector;
			noalias(InterpRadiusVector) = RadiusVector;
			for (size_t index = 0; index < n_damage; index++)
			{
				InterpStressVectorXX(index) = StressMatrix(0, index);
				InterpStressVectorYY(index) = StressMatrix(1, index);
				InterpStressVectorZZ(index) = StressMatrix(2, index);
				InterpStressVectorXY(index) = StressMatrix(3, index);
				InterpStressVectorYZ(index) = StressMatrix(4, index);
				InterpStressVectorXZ(index) = StressMatrix(5, index);
			}
		}
		else
		{
			// local_x and y defines where to estimate the value on the interpolated line,
			// it is 0 at the first point and 1 and the second point
			double local_x = (real_tag[0] - min_tag[0]);
			double local_y = (real_tag[1] - min_tag[1]);
			double local_z = (real_tag[2] - min_tag[2]);
			double local_k = (real_tag[3] - min_tag[3]);
			double local_l = (real_tag[4] - min_tag[4]);

			// std::cout << "**************************************************************" << std::endl;
			for (size_t i_dam = 0; i_dam < n_damage; i_dam++)
			{
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> p;
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> r;
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> xx;
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> yy;
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> zz;
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> xy;
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> yz;
				array_1d< array_1d< array_1d< array_1d< array_1d<double, 2>, 2>, 2>, 2>, 2> xz;

				size_t k = 0;
				for (size_t kk = 0; kk < 2; kk++)
				{
					size_t l = 0;
					for (size_t ll = 0; ll < 2; ll++)
					{
						size_t m = 0;
						for (size_t mm = 0; mm < 2; ll++)
						{
							size_t n = 0;
							for (size_t nn = 0; nn < 2; ll++)
							{
								size_t o = 0;
								for (size_t oo = 0; oo < 2; ll++)
								{
									vector<int> tag_i(5);
									tag_i[0] = min_tag[0] + kk;
									tag_i[1] = min_tag[1] + ll;
									tag_i[2] = min_tag[2] + mm;
									tag_i[3] = min_tag[3] + nn;
									tag_i[4] = min_tag[4] + oo;

									for (size_t iter = 0; iter < tag_i.size(); iter++)
									{
										if (abs(tag_i[iter]) > m)
										{
											double mult_0 = abs(floor(tag_i[iter] / m));
											if (tag_i[iter] > m)
												tag_i[iter] -= m*mult_0;
											else
												tag_i[iter] += m*mult_0;
										}
									}

									Vector DamageVector(n_damage, 0.0);
									Vector RadiusVector(n_damage, 0.0);
									Matrix StressMatrix(strain_size, n_damage, 0.0);

									mpRveMaterialDatabase->GetMaterialValues(tag_i, DamageVector, RadiusVector, StressMatrix);

									p[k][l][m][n][o] = DamageVector[i_dam];
									r[k][l][m][n][o] = RadiusVector[i_dam];
									xx[k][l][m][n][o] = StressMatrix(0, i_dam);
									yy[k][l][m][n][o] = StressMatrix(1, i_dam);
									zz[k][l][m][n][o] = StressMatrix(2, i_dam);
									xy[k][l][m][n][o] = StressMatrix(3, i_dam);
									yz[k][l][m][n][o] = StressMatrix(4, i_dam);
									xz[k][l][m][n][o] = StressMatrix(5, i_dam);
									oo++;
								}
								n++;
							}
							m++;
						}
						l++;
					}
					k++;
				}
				InterpDamageVector[i_dam] = this->pentaDimInterpolate(p, local_x, local_y, local_z, local_k, local_l);
				InterpRadiusVector[i_dam] = this->pentaDimInterpolate(r, local_x, local_y, local_z, local_k, local_l);
				InterpStressVectorXX[i_dam] = this->pentaDimInterpolate(xx, local_x, local_y, local_z, local_k, local_l);
				InterpStressVectorYY[i_dam] = this->pentaDimInterpolate(yy, local_x, local_y, local_z, local_k, local_l);
				InterpStressVectorZZ[i_dam] = this->pentaDimInterpolate(zz, local_x, local_y, local_z, local_k, local_l);
				InterpStressVectorXY[i_dam] = this->pentaDimInterpolate(xy, local_x, local_y, local_z, local_k, local_l);
				InterpStressVectorYZ[i_dam] = this->pentaDimInterpolate(yz, local_x, local_y, local_z, local_k, local_l);
				InterpStressVectorXZ[i_dam] = this->pentaDimInterpolate(xz, local_x, local_y, local_z, local_k, local_l);
			}
			// std::cout << "**************************************************************" << std::endl;
		}

		//Regularization of Fracture Energy with lf
		//std::cout << "damages: " << damages << std::endl;
		this->CalculateFractureEnergy(InterpRadiusVector, Theta, InterpStressVectorXX, InterpStressVectorYY, InterpStressVectorZZ, InterpStressVectorXY, InterpStressVectorYZ, InterpStressVectorXZ);
		this->CalculateAlpha(mlch);
		//
		//std::cout << "InterpRadiusVector BEFORE REG: " << InterpRadiusVector << std::endl;
		//std::cout << "InterpStressVector BEFORE REG: " << InterpStressVectorXY << std::endl;
		this->RegularizeInterpRadius(InterpRadiusVector, Theta);
		//std::cout << "InterpRadiusVector AFTER REG: " << InterpRadiusVector << std::endl;

		//double gf = 0.0;
		//this->CheckFractureEnergy(InterpRadiusVector, Theta, InterpStressVectorXX, InterpStressVectorYY, InterpStressVectorXY, gf);
		//if ((abs(mGf_bar - gf) / gf) > 1.1)
		//{
		//	std::wcout << "real_tag: " << real_tag << std::endl;
		//	std::cout << "gf_ref*lch_ref: " << mGf_bar*1.0e3 << " ; gf*lch_macro: " << gf*lch_macro << std::endl;
		//}

		if (real_radius < InterpRadiusVector[0])
		{
			//std::cout << "Is_Elastic = true\n";
			double d = 0.0;
			Matrix const_tensor = (1 - d)*mC0;
			noalias(stress_vector) = prod(const_tensor, strain_vector);
			//return;
		}
		else if (real_radius > InterpRadiusVector[n_damage - 1]) // || mDamage == 1.0)
		{
			if (mDamage >= 1.0)
			{
				double d = 1.0;
				Matrix const_tensor = (1 - d)*mC0;
				noalias(stress_vector) = prod(const_tensor, strain_vector);
				//return;
			}
			else
			{
				double x = real_radius;
				double x_1 = InterpRadiusVector[n_damage - 2];
				double x0 = InterpRadiusVector[n_damage - 1];
				double diff = (x0 - x_1);
				if (diff < 1.0e-12)
					diff = 1.0;

				stress_vector[0] = InterpStressVectorXX[n_damage - 1] + (InterpStressVectorXX[n_damage - 1] - InterpStressVectorXX[n_damage - 2])*(x - x0) / diff;
				if (SIGMA_SIGN(stress_vector[0]) != SIGMA_SIGN(InterpStressVectorXX[n_damage - 2]))
					stress_vector[0] = 1.0e-3*InterpStressVectorXX[0];

				stress_vector[1] = InterpStressVectorYY[n_damage - 1] + (InterpStressVectorYY[n_damage - 1] - InterpStressVectorYY[n_damage - 2])*(x - x0) / diff;
				if (SIGMA_SIGN(stress_vector[1]) != SIGMA_SIGN(InterpStressVectorYY[n_damage - 2]))
					stress_vector[1] = 1.0e-3*InterpStressVectorYY[0];

				stress_vector[2] = InterpStressVectorZZ[n_damage - 1] + (InterpStressVectorZZ[n_damage - 1] - InterpStressVectorZZ[n_damage - 2])*(x - x0) / diff;
				if (SIGMA_SIGN(stress_vector[2]) != SIGMA_SIGN(InterpStressVectorZZ[n_damage - 2]))
					stress_vector[2] = 1.0e-3*InterpStressVectorZZ[0];

				stress_vector[2] = InterpStressVectorXY[n_damage - 1] + (InterpStressVectorXY[n_damage - 1] - InterpStressVectorXY[n_damage - 2])*(x - x0) / diff;
				if (SIGMA_SIGN(stress_vector[2]) != SIGMA_SIGN(InterpStressVectorXY[n_damage - 2]))
					stress_vector[2] = 1.0e-3*InterpStressVectorXY[0];

				stress_vector[2] = InterpStressVectorYZ[n_damage - 1] + (InterpStressVectorYZ[n_damage - 1] - InterpStressVectorYZ[n_damage - 2])*(x - x0) / diff;
				if (SIGMA_SIGN(stress_vector[2]) != SIGMA_SIGN(InterpStressVectorYZ[n_damage - 2]))
					stress_vector[2] = 1.0e-3*InterpStressVectorYZ[0];

				stress_vector[2] = InterpStressVectorXZ[n_damage - 1] + (InterpStressVectorXZ[n_damage - 1] - InterpStressVectorXZ[n_damage - 2])*(x - x0) / diff;
				if (SIGMA_SIGN(stress_vector[2]) != SIGMA_SIGN(InterpStressVectorXZ[n_damage - 2]))
					stress_vector[2] = 1.0e-3*InterpStressVectorXZ[0];
			}
		}
		else
		{
			for (size_t jj = 1; jj < n_damage; jj++)
			{
				if (real_radius <= InterpRadiusVector[jj])
				{
					double x = real_radius;
					double x0 = InterpRadiusVector[jj - 1];
					double x1 = InterpRadiusVector[jj];

					double diff = (x1 - x0);
					if (diff < 1.0e-12)
						diff = 1.0;

					double y0 = InterpDamageVector[jj - 1];
					double y1 = InterpDamageVector[jj];

					//mDamage = y0 + (y1 - y0)*(x - x0) / diff;

					stress_vector[0] = InterpStressVectorXX[jj - 1] + (InterpStressVectorXX[jj] - InterpStressVectorXX[jj - 1])*(x - x0) / diff;
					stress_vector[1] = InterpStressVectorYY[jj - 1] + (InterpStressVectorYY[jj] - InterpStressVectorYY[jj - 1])*(x - x0) / diff;
					stress_vector[2] = InterpStressVectorZZ[jj - 1] + (InterpStressVectorZZ[jj] - InterpStressVectorZZ[jj - 1])*(x - x0) / diff;
					stress_vector[3] = InterpStressVectorXY[jj - 1] + (InterpStressVectorXY[jj] - InterpStressVectorXY[jj - 1])*(x - x0) / diff;
					stress_vector[4] = InterpStressVectorYZ[jj - 1] + (InterpStressVectorYZ[jj] - InterpStressVectorYZ[jj - 1])*(x - x0) / diff;
					stress_vector[5] = InterpStressVectorXZ[jj - 1] + (InterpStressVectorXZ[jj] - InterpStressVectorXZ[jj - 1])*(x - x0) / diff;

					break;
				}
			}
		}

		//std::cout << "stress_vector: " << stress_vector << std::endl;
		//Calculate damage
		Vector elastic_stress(strain_vector.size(), 0.0);
		noalias(elastic_stress) = prod(mC0, strain_vector);
		// Calculate norm
		double norm_elastic_stress = norm_2(elastic_stress);
		//double norm_stress_difference = norm_2(stress_vector - elastic_stress);
		double norm_stress_difference = norm_2(elastic_stress - stress_vector);

		// Calculate equivalent damage
		if (norm_elastic_stress < 1.0e-12)
			mDamage = 0.0;
		else
			mDamage = norm_stress_difference / norm_elastic_stress;

		mDamage = std::max(std::min(mDamage, 1.0), 0.0);
		if (mDamage < mDamageConverged)
			mDamage = mDamageConverged;

		//Calculate damaged Constitutive Tensor
		//noalias(const_tensor) = (1 - mDamage)*mC0;
		//noalias(const_tensor) = mC0;

		//Calculate stress
		//stress_vector.clear();
		//noalias(stress_vector) = prod(const_tensor, trial_macro_scale_strain_vector);

		//noalias(const_tensor) = (1 - EquivalentDamageConverged)*mC0;
	}

	void InterpolatedConstitutiveLaw3D::BiLinearInterpolation(const size_t& m, const Vector& xy, const vector<int>& min_xy)
	{
		// local_x and y defines where to estimate the value on the interpolated line,
		// it is 0 at the first point and 1 and the second point
		double local_x = (xy[0] - min_xy[0]); // / (max_tag[0] - min_tag[0]);
		double local_y = (xy[1] - min_xy[1]); // / (max_tag[1] - min_tag[1]);

		vector<double> xi(4);		vector<double> eta(4);
		xi[0] = -1.0;				eta[0] = -1.0;
		xi[1] = -1.0;				eta[1] = 1.0;
		xi[2] = 1.0;				eta[2] = -1.0;
		xi[3] = 1.0;				eta[3] = 1.0;

		size_t i = 0;
		double N_i = 0;
		double N_sum = 0;

		for (size_t kk = 0; kk < 2; kk++)
		{
			for (size_t ll = 0; ll < 2; ll++)
			{
				vector<int> tag_i(2);
				tag_i[0] = min_xy[0] + kk;
				tag_i[1] = min_xy[1] + ll;
				if (abs(tag_i[0]) > m)
				{
					double mult_0 = abs(floor(tag_i[0] / m));
					if (tag_i[0] > m)
						tag_i[0] -= m*mult_0;
					else
						tag_i[0] += m*mult_0;
				}

				if (abs(tag_i[1]) > m)
				{
					double mult_1 = abs(floor(tag_i[1] / m));
					if (tag_i[1] > m)
						tag_i[1] -= m*mult_1;
					else
						tag_i[1] += m*mult_1;
				}

				N_i = 1.0 / 4.0 * (1.0 + xi[i] * (local_x * 2.0 - 1.0))*(1.0 + eta[i] * (local_y * 2.0 - 1.0));
				N_sum += N_i;
				i++;

				// DO SOMETHING
			}
		}
		if (abs(1 - N_sum) > 1.0e-15)
		{
			std::cout << "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION ON PredictElasticity - N_sum = " << N_sum << " ; shape_index = " << i << std::endl;
			KRATOS_THROW_ERROR(std::logic_error, "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
		}
	}

	void InterpolatedConstitutiveLaw3D::ReconstructStrain(Vector& eps, const double& radius, const Vector& Theta)
	{
		eps[0] = radius*cos(Theta[0]);
		eps[1] = radius*sin(Theta[0])*cos(Theta[1]);
		eps[2] = radius*sin(Theta[0])*sin(Theta[1])*cos(Theta[2]);
		eps[3] = radius*sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*cos(Theta[3]);
		eps[4] = radius*sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*cos(Theta[4]);
		eps[5] = radius*sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*sin(Theta[4]);
	}

	void InterpolatedConstitutiveLaw3D::CalculateFractureEnergy(const Vector& radius, const Vector& Theta, Vector& sigma_xx, Vector& sigma_yy, Vector& sigma_zz, Vector& sigma_xy, Vector& sigma_yz, Vector& sigma_xz)
	{
		const size_t n_damage = radius.size();

		Vector norm_s(n_damage, 0.0);
		for (size_t i = 0; i < n_damage; i++)
		{
			norm_s[i] = sqrt(sigma_xx[i] * sigma_xx[i] + sigma_yy[i] * sigma_yy[i] + sigma_zz[i] * sigma_zz[i] + 2.0*sigma_xy[i] * sigma_xy[i] + 2.0*sigma_yz[i] * sigma_yz[i] + 2.0*sigma_xz[i] * sigma_xz[i]);
		}

		auto biggest_ft = std::max_element(std::begin(norm_s), std::end(norm_s));

		m_position_biggest_ft = std::distance(std::begin(norm_s), biggest_ft);

		//Calculate Area -> Fracture Energy
		Vector e0(3, 0.0);
		this->ReconstructStrain(e0, radius[0], Theta);
		double norm_e0 = sqrt(e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2] + 2.0 * ((0.5*e0[3]) * (0.5*e0[3])) + 2.0 * ((0.5*e0[4]) * (0.5*e0[4])) + 2.0 * ((0.5*e0[5]) * (0.5*e0[5])));
		mG0 = 0.5 * norm_e0 * norm_s[0];
		mGf_bar = mG0;
		Vector ei1(6, 0.0);
		Vector ei2(6, 0.0);
		double norm_ei1 = 0.0;
		double norm_ei2 = 0.0;
		for (size_t i = m_position_biggest_ft - 1; i < n_damage - 1; i++)
			//for (size_t i = 0; i < n_damage - 1; i++)
		{
			this->ReconstructStrain(ei1, radius[i], Theta);
			this->ReconstructStrain(ei2, radius[i + 1], Theta);
			norm_ei1 = sqrt(ei1[0] * ei1[0] + ei1[1] * ei1[1] + ei1[2] * ei1[2] + 2.0 * ((0.5*ei1[3]) * (0.5*ei1[3])) + 2.0 * ((0.5*ei1[4]) * (0.5*ei1[4])) + 2.0 * ((0.5*ei1[5]) * (0.5*ei1[5])));
			norm_ei2 = sqrt(ei2[0] * ei2[0] + ei2[1] * ei2[1] + ei2[2] * ei2[2] + 2.0 * ((0.5*ei2[3]) * (0.5*ei2[3])) + 2.0 * ((0.5*ei2[4]) * (0.5*ei2[4])) + 2.0 * ((0.5*ei2[5]) * (0.5*ei2[5])));
			mGf_bar += std::abs(0.5 * (norm_ei2 - norm_ei1)*(norm_s[i] + norm_s[i + 1]));
		}
	}

	void InterpolatedConstitutiveLaw3D::CalculateAlpha(const double& lch_macro)
	{
		// TODO: DEFINE lch_ref AS INPUT
		double lch_ref = mlchRef;
		double Gf = 0.0;

		//lch_ref = 1.1283791670955 * sqrt(fabs(Area));;
		//mAlpha = (Gc-G1)/(G_bar-G1)-1.0;
		double mult = lch_ref / lch_macro;
		//std::cout << "lch_ref: " << lch_ref << std::endl;
		//std::cout << "lch_macro: " << lch_macro << std::endl;
		//xx
		Gf = mGf_bar * mult;
		//mAlpha = lch_ref / lch_macro - 1.0;
		mAlpha = (Gf - mG0) / (mGf_bar - mG0) - 1.0;
		//std::cout << "mAlphaxx: " << mAlphaxx << std::endl;
		if (mAlpha <= -1.0)
		{
			std::stringstream ss;
			ss << "RvePredictorCalculator Error: mAlphaxx <= -1.0" << std::endl;
			ss << "mAlpha = " << mAlpha << std::endl;
			ss << "mG0 = " << mG0 << std::endl;
			ss << "mGf_bar = " << mGf_bar << std::endl;
			std::cout << ss.str();
			exit(-1);
		}
		//mAlpha = -0.80;
	}

	void InterpolatedConstitutiveLaw3D::RegularizeInterpRadius(Vector& InterpRadius, const Vector& Theta)
	{
		size_t n_radius = InterpRadius.size();
		for (size_t i = m_position_biggest_ft - 1; i < n_radius; i++)
			//for (size_t i = 0; i < n_radius; i++)
		{
			// eps[i] = eps[i] + mAlpha * (eps[i] - eps[0]);
			// ej = ej+(ej-ep)*S;

			Vector e(3, false);
			e[0] = InterpRadius[i] * cos(Theta[0]) + mAlpha * (InterpRadius[i] * cos(Theta[0]) - InterpRadius[m_position_biggest_ft - 1] * cos(Theta[0]));
			e[1] = InterpRadius[i] * sin(Theta[0])*cos(Theta[1]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*cos(Theta[1]) - InterpRadius[m_position_biggest_ft - 1] * sin(Theta[0])*cos(Theta[1]));
			e[2] = InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*cos(Theta[2]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*cos(Theta[2]) - InterpRadius[m_position_biggest_ft - 1] * sin(Theta[0])*sin(Theta[1])*cos(Theta[2]));
			e[3] = InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*cos(Theta[3]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*cos(Theta[3]) - InterpRadius[m_position_biggest_ft - 1] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*cos(Theta[3]));
			e[4] = InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*cos(Theta[4]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*cos(Theta[4]) - InterpRadius[m_position_biggest_ft - 1] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*cos(Theta[4]));
			e[5] = InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*sin(Theta[4]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*sin(Theta[4]) - InterpRadius[m_position_biggest_ft - 1] * sin(Theta[0])*sin(Theta[1])*sin(Theta[2])*sin(Theta[3])*sin(Theta[4]));

			InterpRadius[i] = norm_2(e);

			//InterpRadius[i] = InterpRadius[i] + mAlpha*(InterpRadius[i] - InterpRadius[m_position_biggest_ft]);
		}
	}

	void InterpolatedConstitutiveLaw3D::CalculateConstitutiveMatrix(const Properties& props,
		const Vector& strainVector,
		const Vector& stressVector,
		Matrix& constitutiveMatrix)
	{
		size_t strain_size = GetStrainSize();
		Matrix aux_matrix(strain_size, strain_size, 0.0);
		double aux_damage = 0.0;
		double aux_damage2 = 0.0;
		// perturbation parameter
		double h = 0.0;
		h = std::max(1.0e-8, 1.0e-6*norm_2(strainVector));
		// perturbed vectors
		Vector strain_bar(strain_size);
		Vector S1(strain_size);
		Vector S2(strain_size);
		// apply perturbation to each strain component...
		for (size_t j = 0; j < strain_size; j++)
		{
			noalias(strain_bar) = strainVector;

			strain_bar(j) = strainVector(j) - h;
			CalculateStress(props, strain_bar, S1);

			strain_bar(j) = strainVector(j) + h;
			CalculateStress(props, strain_bar, S2);

			for (size_t i = 0; i < strain_size; i++)
				constitutiveMatrix(i, j) = (S2(i) - S1(i)) / (2.0*h);
		}
		//std::stringstream ss;
		//ss << "Perturbed - constitutiveMatrix = " << constitutiveMatrix << ", " << std::endl;
		//std::cout << ss.str();
	}

	Vector InterpolatedConstitutiveLaw3D::CalculateTheta(const Vector& StrainVector) const {
		//Calculate Theta (for the real strain)
		size_t theta_size = StrainVector.size() - 1;
		//std::cout << "PredictStress2D -> StrainVector: " << StrainVector << std::endl;
		Vector Theta(theta_size, 0.0);

		//Calculate Theta1
		double theta_1 = acos(StrainVector[0] / sqrt(pow(StrainVector[5], 2) + pow(StrainVector[4], 2) + pow(StrainVector[3], 2) + pow(StrainVector[2], 2) + pow(StrainVector[1], 2) + pow(StrainVector[0], 2)));
		if (abs(theta_1) > Globals::Pi)
		{
			double mult = abs(floor(theta_1 / Globals::Pi));
			if (theta_1 > Globals::Pi)
				theta_1 = theta_1 - 2.0 * Globals::Pi*mult;
			else
				theta_1 = theta_1 + 2.0 * Globals::Pi*mult;
		}

		//Calculate Theta2
		double theta_2 = acos(StrainVector[1] / sqrt(pow(StrainVector[5], 2) + pow(StrainVector[4], 2) + pow(StrainVector[3], 2) + pow(StrainVector[2], 2) + pow(StrainVector[1], 2)));
		if (abs(theta_2) > Globals::Pi)
		{
			double mult = abs(floor(theta_2 / Globals::Pi));
			if (theta_2 > Globals::Pi)
				theta_2 = theta_2 - 2.0 * Globals::Pi*mult;
			else
				theta_2 = theta_2 + 2.0 * Globals::Pi*mult;
		}

		//Calculate Theta3
		double theta_3 = acos(StrainVector[2] / sqrt(pow(StrainVector[5], 2) + pow(StrainVector[4], 2) + pow(StrainVector[3], 2) + pow(StrainVector[2], 2)));
		if (abs(theta_3) > Globals::Pi)
		{
			double mult = abs(floor(theta_3 / Globals::Pi));
			if (theta_3 > Globals::Pi)
				theta_3 = theta_3 - 2.0 * Globals::Pi*mult;
			else
				theta_3 = theta_3 + 2.0 * Globals::Pi*mult;
		}

		//Calculate Theta4
		double theta_4 = acos(StrainVector[2] / sqrt(pow(StrainVector[5], 2) + pow(StrainVector[4], 2) + pow(StrainVector[3], 2)));
		if (abs(theta_4) > Globals::Pi)
		{
			double mult = abs(floor(theta_4 / Globals::Pi));
			if (theta_4 > Globals::Pi)
				theta_4 = theta_4 - 2.0 * Globals::Pi*mult;
			else
				theta_4 = theta_4 + 2.0 * Globals::Pi*mult;
		}

		//Calculate Theta5
		double theta_5;
		double denom = sqrt(pow(StrainVector[5], 2) + pow(StrainVector[4], 2));
		if (denom == 0)
		{
			theta_5 = 0;
		}
		else
		{
			if (StrainVector[5] < 0.0)
			{
				theta_5 = 2 * Globals::Pi - acos(StrainVector[4] / denom);
			}
			else
			{
				theta_5 = acos(StrainVector[4] / denom);
			}
		}
		if (abs(theta_5) > Globals::Pi)
		{
			double mult = abs(floor(theta_5 / Globals::Pi));
			if (theta_5 > Globals::Pi)
				theta_5 = theta_5 - 2.0 * Globals::Pi*mult;
			else
				theta_5 = theta_5 + 2.0 * Globals::Pi*mult;
		}

		//Assemble Theta
		Theta[0] = theta_1;
		Theta[1] = theta_2;
		Theta[2] = theta_3;
		Theta[3] = theta_4;
		Theta[4] = theta_5;

		return Theta;
	}

	void InterpolatedConstitutiveLaw3D::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
		rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
		rFeatures.mOptions.Set(ISOTROPIC);

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}

	int InterpolatedConstitutiveLaw3D::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
			if (!rMaterialProperties.Has(LCH_REF_RVE)) {
				KRATOS_THROW_ERROR(std::logic_error, "InterpolatedConstitutiveLaw3D - missing LCH_REF_RVE", "");
			}
		return 0;

		KRATOS_CATCH("");
	}

} /* namespace Kratos.*/
