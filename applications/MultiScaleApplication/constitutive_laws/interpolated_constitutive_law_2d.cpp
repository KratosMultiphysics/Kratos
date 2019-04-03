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

#include "interpolated_constitutive_law_2d.h"
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

	InterpolatedConstitutiveLaw2D::InterpolatedConstitutiveLaw2D(const RveMaterialDatabase::Pointer& pNewRveMaterialDatabase)
		: ConstitutiveLaw()
		, mpRveMaterialDatabase(RveMaterialDatabase::Pointer())
		, mInitialized(false)
		, mInitStrain()
		, mErrorCode(0.0)
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

	ConstitutiveLaw::Pointer InterpolatedConstitutiveLaw2D::Clone() const
	{
		return ConstitutiveLaw::Pointer(new InterpolatedConstitutiveLaw2D(mpRveMaterialDatabase));
	}

	InterpolatedConstitutiveLaw2D::SizeType InterpolatedConstitutiveLaw2D::WorkingSpaceDimension()
	{
		return 2;
	}

	InterpolatedConstitutiveLaw2D::SizeType InterpolatedConstitutiveLaw2D::GetStrainSize()
	{
		return 3;
	}

	bool InterpolatedConstitutiveLaw2D::Has(const Variable<double>& rThisVariable)
	{
		return false;
	}

	bool InterpolatedConstitutiveLaw2D::Has(const Variable<Vector>& rThisVariable)
	{
		if (rThisVariable == INITIAL_STRAIN) return true;
		return false;
	}

	bool InterpolatedConstitutiveLaw2D::Has(const Variable<Matrix>& rThisVariable)
	{
		return false;
	}

	bool InterpolatedConstitutiveLaw2D::Has(const Variable<array_1d<double, 2 > >& rThisVariable)
	{
		return false;
	}

	bool InterpolatedConstitutiveLaw2D::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return false;
	}

	double& InterpolatedConstitutiveLaw2D::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == LCH_REF_RVE)
			rValue = mlchRef;
		else if (rThisVariable == EQUIVALENT_DAMAGE)
			rValue = mDamage;
		else if (rThisVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			rValue = mErrorCode;
		return rValue;
	}

	Vector& InterpolatedConstitutiveLaw2D::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		if (rThisVariable == INITIAL_STRAIN) {
			if (rValue.size() != mInitStrain.size())
				rValue.resize(mInitStrain.size());
			noalias(rValue) = mInitStrain;
		}
		return rValue;
	}

	Matrix& InterpolatedConstitutiveLaw2D::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		if (rThisVariable == PREDICTED_STRESS_TENSOR)
			StressVectorToTensor(mStressVector, rValue);
		return rValue;
	}

	array_1d<double, 2 > & InterpolatedConstitutiveLaw2D::GetValue(const Variable<array_1d<double, 2 > >& rVariable, array_1d<double, 2 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & InterpolatedConstitutiveLaw2D::GetValue(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	void InterpolatedConstitutiveLaw2D::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == LCH_REF_RVE) {
			mlchRef = rValue;
		}
	}

	void InterpolatedConstitutiveLaw2D::SetValue(const Variable<Vector >& rVariable,
		const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == INITIAL_STRAIN) {
			if (rValue.size() == mInitStrain.size())
				noalias(mInitStrain) = rValue;
		}
	}

	void InterpolatedConstitutiveLaw2D::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw2D::SetValue(const Variable<array_1d<double, 2 > >& rVariable,
		const array_1d<double, 2 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw2D::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	bool InterpolatedConstitutiveLaw2D::ValidateInput(const Properties& rMaterialProperties)
	{
		if (!rMaterialProperties.Has(LCH_REF_RVE)) return false;
		if (!rMaterialProperties.Has(THICKNESS)) return false;
		return true;
	}

	InterpolatedConstitutiveLaw2D::StrainMeasure InterpolatedConstitutiveLaw2D::GetStrainMeasure()
	{
		return ConstitutiveLaw::StrainMeasure_Infinitesimal;
	}

	InterpolatedConstitutiveLaw2D::StressMeasure InterpolatedConstitutiveLaw2D::GetStressMeasure()
	{
		return ConstitutiveLaw::StressMeasure_Cauchy;
	}

	bool InterpolatedConstitutiveLaw2D::IsIncremental()
	{
		return false;
	}

	void InterpolatedConstitutiveLaw2D::InitializeMaterial(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (!mInitialized)
		{
			mErrorCode = 0.0;
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

	void InterpolatedConstitutiveLaw2D::InitializeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw2D::FinalizeSolutionStep(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		m_r_converged = m_r;
		mDamageConverged = mDamage;
	}

	void InterpolatedConstitutiveLaw2D::InitializeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw2D::FinalizeNonLinearIteration(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void InterpolatedConstitutiveLaw2D::CalculateMaterialResponsePK1(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw2D::CalculateMaterialResponsePK2(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw2D::CalculateMaterialResponseKirchhoff(Parameters& rValues)
	{
		CalculateMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw2D::CalculateMaterialResponseCauchy(Parameters& rValues)
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

	void InterpolatedConstitutiveLaw2D::FinalizeMaterialResponsePK1(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw2D::FinalizeMaterialResponsePK2(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw2D::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
	{
		FinalizeMaterialResponseCauchy(rValues);
	}

	void InterpolatedConstitutiveLaw2D::FinalizeMaterialResponseCauchy(Parameters& rValues)
	{
	}

	void InterpolatedConstitutiveLaw2D::ResetMaterial(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues)
	{
		if (mInitialized)
		{
			m_r = 0.0;
			m_r_converged = 0.0;
			mErrorCode = 0.0;
			mInitStrain = ZeroVector(this->GetStrainSize());
			mDamage = 0.0;
			mDamageConverged = 0.0;
			mIsElastic = false;
			mlch = 0.0;
			mInitialized = false;
		}
	}

	void InterpolatedConstitutiveLaw2D::CalculateStress(const Properties& props,
		const Vector& strainVector,
		Vector& stressVector)
	{
		m_r = m_r_converged;

		size_t strain_size = strainVector.size();
		mStressVector.clear();

		//Initialize data for the interpolation
		Vector Theta(strain_size - 1);
		Vector eps(strain_size);
		double radius = 0.0;

		//Calculate phi & number of division
		const size_t m = mpRveMaterialDatabase->NumberOfDivision();
		//std::cout << "m: " << m << std::endl;
		double phi = Globals::Pi / m;
		//std::cout << "phi: " << phi << " rad" << std::endl;

		//1. Calculate the normalized Strain
		Matrix S(3, 3, false);
		S(0, 0) = 1.0;
		S(1, 1) = 1.0;
		S(2, 2) = 1.0;
		Vector Se = prod(S, strainVector);
		double eTSe = strainVector(0)*Se(0) + strainVector(1)*Se(1) + strainVector(2)*Se(2);
		double real_radius = sqrt(eTSe);
		//std::cout << "real_radius: " << real_radius << std::endl;

		Theta = this->CalculateTheta(strainVector);
		//std::cout << "Theta: " << Theta << std::endl;

		Vector real_tag = Theta / phi;
		//std::cout << "PredictElasticity -> real_tag: " << real_tag << std::endl;
		if (abs(real_tag[0]) > m)
		{
			double mult_0 = abs(floor(real_tag[0] / m));
			if (real_tag[0] > m)
				real_tag[0] -= m*mult_0;
			else
				real_tag[0] += m*mult_0;
		}

		if (abs(real_tag[1]) > m)
		{
			double mult_1 = abs(floor(real_tag[1] / m));
			if (real_tag[1] > m)
				real_tag[1] -= m*mult_1;
			else
				real_tag[1] += m*mult_1;
		}

		//BI-LINEAR INTERPOLATION OVER TAG
		//1. Found the minimum value of the tag
		vector<int> min_tag(2);
		min_tag[0] = floor(real_tag[0]);
		min_tag[1] = floor(real_tag[1]);
		//std::cout << "min_tag: " << min_tag << std::endl;
		vector<int> max_tag(2);
		max_tag[0] = ceil(real_tag[0]);
		max_tag[1] = ceil(real_tag[1]);
		//std::cout << "max_tag: " << max_tag << std::endl;

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

		//noalias(stressVector) = (1.0 - mDamageConverged)*mStressVector;
		noalias(stressVector) = mStressVector;
	}

	void InterpolatedConstitutiveLaw2D::TagsLinearInterpolation(const Vector& xy, const vector<int>& min_xy, double& z)
	{
		double m = mpRveMaterialDatabase->NumberOfDivision();

		if (xy[0] == min_xy[0] && xy[1] == min_xy[1])
		{
			mpRveMaterialDatabase->GetElasticRadius(xy,z);
		}
		else
		{
			// local_x and y defines where to estimate the value on the interpolated line,
			// it is 0 at the first point and 1 and the second point
			double local_x = (xy[0] - min_xy[0]); // / (max_tag[0] - min_tag[0]);
			double local_y = (xy[1] - min_xy[1]); // / (max_tag[1] - min_tag[1]);
			//std::cout << "local_x, local_y: " << local_x << ", " << local_y << std::endl;
			Vector coordinates(2, 0.0);
			coordinates[0] = local_x;
			coordinates[1] = local_y;

			vector<double> xi(4);
			xi[0] = -1.0;
			xi[1] = -1.0;
			xi[2] = 1.0;
			xi[3] = 1.0;

			vector<double> eta(4);
			eta[0] = -1.0;
			eta[1] = 1.0;
			eta[2] = -1.0;
			eta[3] = 1.0;

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

					//std::cout << "tag_i: " << tag_i << std::endl;

					double i_z = 0.0;
					mpRveMaterialDatabase->GetElasticRadius(tag_i, i_z);
					//std::cout << "i_radius: " << i_radius << std::endl;

					N_i = 1.0 / 4.0 * (1.0 + xi[i] * (local_x * 2.0 - 1.0))*(1.0 + eta[i] * (local_y * 2.0 - 1.0));
					//N_i = ShapeFunctionValue4N(i, coordinates);
					N_sum += N_i;

					//std::cout << "N_i: " << N_i << std::endl;
					z += N_i*i_z;

					i++;
				}
			}
			if (abs(1 - N_sum) > 1.0e-15)
			{
				std::cout << "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION ON PredictElasticity - N_sum = " << N_sum << " ; shape_index = " << i << std::endl;
				KRATOS_THROW_ERROR(std::logic_error, "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
			}
		}
	}

	void InterpolatedConstitutiveLaw2D::StressPrediction(const Vector& real_tag, const vector<int>& min_tag, const double& real_radius, const Vector& Theta, const Vector& strain_vector, Vector& stress_vector)
	{
		size_t strain_size = strain_vector.size();
		double m = mpRveMaterialDatabase->NumberOfDivision();
		const size_t n_damage = mpRveMaterialDatabase->DamageSize();

		vector<double> xi(4);		vector<double> eta(4);
		xi[0] = -1.0;				eta[0] = -1.0;
		xi[1] = -1.0;				eta[1] = 1.0;
		xi[2] = 1.0;				eta[2] = -1.0;
		xi[3] = 1.0;				eta[3] = 1.0;

		Vector InterpDamageVector(n_damage, 0.0);
		Vector InterpRadiusVector(n_damage, 0.0);
		Vector InterpStressVectorXX(n_damage, 0.0);
		Vector InterpStressVectorYY(n_damage, 0.0);
		Vector InterpStressVectorXY(n_damage, 0.0);

		if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
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
				InterpStressVectorXY(index) = StressMatrix(2, index);
			}
		}
		else
		{
			// local_x and y defines where to estimate the value on the interpolated line,
			// it is 0 at the first point and 1 and the second point
			double local_x = (real_tag[0] - min_tag[0]); // / (max_tag[0] - min_tag[0]);
			double local_y = (real_tag[1] - min_tag[1]); // / (max_tag[1] - min_tag[1]);

			Vector coordinates(2, 0.0);
			coordinates[0] = local_x;
			coordinates[1] = local_y;
			// std::cout << "**************************************************************" << std::endl;
			for (size_t i_dam = 0; i_dam < n_damage; i_dam++)
			{
				size_t i = 0;
				double N_i = 0;
				double N_sum = 0;

				for (size_t kk = 0; kk < 2; kk++)
				{
					for (size_t ll = 0; ll < 2; ll++)
					{
						vector<int> tag_i(2);
						tag_i[0] = min_tag[0] + kk;
						tag_i[1] = min_tag[1] + ll;
						// std::cout << "PredictStress2D -> tag_i: " << tag_i << std::endl;
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

						Vector DamageVector(n_damage, 0.0);
						Vector RadiusVector(n_damage, 0.0);
						Matrix StressMatrix(strain_size, n_damage, 0.0);

						mpRveMaterialDatabase->GetMaterialValues(tag_i, DamageVector, RadiusVector, StressMatrix);
						InterpDamageVector[i_dam] += N_i*DamageVector[i_dam]; // Check it!!!
						InterpRadiusVector[i_dam] += N_i*RadiusVector[i_dam];
						InterpStressVectorXX[i_dam] += N_i*StressMatrix(0, i_dam);
						InterpStressVectorYY[i_dam] += N_i*StressMatrix(1, i_dam);
						InterpStressVectorXY[i_dam] += N_i*StressMatrix(2, i_dam);
					}
				}
				if (abs(1 - N_sum) > 1.0e-15)
				{
					std::cout << "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION ON PredictElasticity - N_sum = " << N_sum << " ; shape_index = " << i << std::endl;
					KRATOS_THROW_ERROR(std::logic_error, "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
				}
			}
			// std::cout << "**************************************************************" << std::endl;
		}

		//Regularization of Fracture Energy with lf
		//std::cout << "damages: " << damages << std::endl;
		this->CalculateFractureEnergy(InterpRadiusVector, Theta, InterpStressVectorXX, InterpStressVectorYY, InterpStressVectorXY);
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

				// stress_vector[0] = InterpStressVectorXX[n_damage - 1] + (InterpStressVectorXX[n_damage - 1] - InterpStressVectorXX[n_damage - 2])*(x - x0) / diff;
				// if (SIGMA_SIGN(stress_vector[0]) != SIGMA_SIGN(InterpStressVectorXX[n_damage - 2]))
				// 	stress_vector[0] = 1.0e-3*InterpStressVectorXX[0];
				//
				// stress_vector[1] = InterpStressVectorYY[n_damage - 1] + (InterpStressVectorYY[n_damage - 1] - InterpStressVectorYY[n_damage - 2])*(x - x0) / diff;
				// if (SIGMA_SIGN(stress_vector[1]) != SIGMA_SIGN(InterpStressVectorYY[n_damage - 2]))
				// 	stress_vector[1] = 1.0e-3*InterpStressVectorYY[0];
				//
				// stress_vector[2] = InterpStressVectorXY[n_damage - 1] + (InterpStressVectorXY[n_damage - 1] - InterpStressVectorXY[n_damage - 2])*(x - x0) / diff;
				// if (SIGMA_SIGN(stress_vector[2]) != SIGMA_SIGN(InterpStressVectorXY[n_damage - 2]))
				// 	stress_vector[2] = 1.0e-3*InterpStressVectorXY[0];

				stress_vector[0] = InterpStressVectorXX[n_damage - 1] + (InterpStressVectorXX[n_damage - 1] - InterpStressVectorXX[n_damage - 2])*(x - x0) / diff;
				stress_vector[1] = InterpStressVectorYY[n_damage - 1] + (InterpStressVectorYY[n_damage - 1] - InterpStressVectorYY[n_damage - 2])*(x - x0) / diff;
				stress_vector[2] = InterpStressVectorXY[n_damage - 1] + (InterpStressVectorXY[n_damage - 1] - InterpStressVectorXY[n_damage - 2])*(x - x0) / diff;

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
					stress_vector[2] = InterpStressVectorXY[jj - 1] + (InterpStressVectorXY[jj] - InterpStressVectorXY[jj - 1])*(x - x0) / diff;

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

	void InterpolatedConstitutiveLaw2D::BiLinearInterpolation(const size_t& m, const Vector& xy, const vector<int>& min_xy)
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

	void InterpolatedConstitutiveLaw2D::ReconstructStrain(Vector& eps, const double& radius, const Vector& Theta)
	{
		eps[0] = radius*cos(Theta[0]);
		eps[1] = radius*sin(Theta[0])*cos(Theta[1]);
		eps[2] = radius*sin(Theta[0])*sin(Theta[1]);
	}

	void InterpolatedConstitutiveLaw2D::CalculateFractureEnergy(const Vector& radius, const Vector& Theta, Vector& sigma_xx, Vector& sigma_yy, Vector& sigma_xy)
	{
		const size_t n_damage = radius.size();

		Vector norm_s(n_damage, 0.0);

		size_t count = 0;
		for (size_t i = 0; i < n_damage; i++)
		{
			norm_s[i] = sqrt(sigma_xx[i] * sigma_xx[i] + sigma_yy[i] * sigma_yy[i] + 2.0*sigma_xy[i] * sigma_xy[i]);
			//norm_s[i] = sqrt(sigma_xx[i] * sigma_xx[i] + sigma_yy[i] * sigma_yy[i] + sigma_xy[i] * sigma_xy[i]);
		}

		//std::cout << "radius: " << radius << std::endl;
		//std::cout << "norm_s: " << norm_s << std::endl;
		auto biggest_ft = std::max_element(std::begin(norm_s), std::end(norm_s));

		m_position_biggest_ft = std::distance(std::begin(norm_s), biggest_ft);
		if (m_position_biggest_ft == n_damage)
		{
			//std::cout << "Real Tag: " << Theta / (3.1415926535897932384626433832795 / 8) << std::endl;
			//std::cout << "m_position_biggest_ft: " << m_position_biggest_ft << std::endl;

			//Print to file
			std::string dName = "normS_vs_normE.txt";
			std::ofstream mdI(dName.c_str(), std::ios_base::app | std::ios_base::out);
			for (size_t i = 0; i < n_damage; i++)
			{
				//std::cout << radius[i] << "	" << norm_s[i] << std::endl;
				mdI << radius[i] << "	" << norm_s[i] << std::endl;
				//std::cout << radius[i] << "	" << norm_s[i] << std::endl;
			}
			//std::cout << "	" << std::endl;
			mdI << "	" << std::endl;
			mdI.flush();
			mdI.close();
		}

		//std::cout << ".....................................................\n";

		//Calculate Area -> Fracture Energy
		mG0 = 0.0;
		if (m_position_biggest_ft == 0)
		{
			Vector e0(3, 0.0);
			this->ReconstructStrain(e0, radius[0], Theta);
			double norm_e0 = sqrt(e0[0] * e0[0] + e0[1] * e0[1] + 2.0 * ((0.5*e0[2]) * (0.5*e0[2])));
			mG0 = 0.5 * norm_e0 * norm_s[0];
		}
		else
		{
			Vector ei1(3, 0.0);
			Vector ei2(3, 0.0);
			double norm_ei1 = 0.0;
			double norm_ei2 = 0.0;
			for (size_t i = 0; i < m_position_biggest_ft - 1; i++)
			{
				this->ReconstructStrain(ei1, radius[i], Theta);
				this->ReconstructStrain(ei2, radius[i + 1], Theta);
				norm_ei1 = sqrt(ei1[0] * ei1[0] + ei1[1] * ei1[1] + 2.0 * ((0.5*ei1[2]) * (0.5*ei1[2])));
				norm_ei2 = sqrt(ei2[0] * ei2[0] + ei2[1] * ei2[1] + 2.0 * ((0.5*ei2[2]) * (0.5*ei2[2])));
				mG0 += std::abs(0.5 * (norm_ei2 - norm_ei1)*(norm_s[i] + norm_s[i + 1]));
				//std::cout << "------------------------> mGf_bar: " << mGf_bar << std::endl;
			}
		}

		mGf_bar = mG0;
		Vector ei1(3, 0.0);
		Vector ei2(3, 0.0);
		double norm_ei1 = 0.0;
		double norm_ei2 = 0.0;
		for (size_t i = m_position_biggest_ft; i < n_damage - 1; i++)
		{
			this->ReconstructStrain(ei1, radius[i], Theta);
			this->ReconstructStrain(ei2, radius[i + 1], Theta);
			norm_ei1 = sqrt(ei1[0] * ei1[0] + ei1[1] * ei1[1] + 2.0 * ((0.5*ei1[2]) * (0.5*ei1[2])));
			norm_ei2 = sqrt(ei2[0] * ei2[0] + ei2[1] * ei2[1] + 2.0 * ((0.5*ei2[2]) * (0.5*ei2[2])));
			mGf_bar += std::abs(0.5 * (norm_ei2 - norm_ei1)*(norm_s[i] + norm_s[i + 1]));
			//std::cout << "------------------------> mGf_bar: " << mGf_bar << std::endl;
		}
	}

	void InterpolatedConstitutiveLaw2D::CheckFractureEnergy(const Vector& radius, const Vector& Theta, Vector& sigma_xx, Vector& sigma_yy, Vector& sigma_xy, double& gf)
	{
		const size_t n_damage = radius.size();

		Vector norm_s(n_damage, 0.0);
		for (size_t i = 0; i < n_damage; i++)
		{
			norm_s[i] = sqrt(sigma_xx[i] * sigma_xx[i] + sigma_yy[i] * sigma_yy[i] + sigma_xy[i] * sigma_xy[i]);
		}

		//Calculate Area -> Fracture Energy
		Vector e0(3, 0.0);
		this->ReconstructStrain(e0, radius[0], Theta);
		double norm_e0 = sqrt(e0[0] * e0[0] + e0[1] * e0[1] + 2.0 * ((0.5*e0[2]) * (0.5*e0[2])));
		double g0 = 0.5 * norm_e0 * norm_s[0];
		gf = g0;
		Vector ei1(3, 0.0);
		Vector ei2(3, 0.0);
		double norm_ei1 = 0.0;
		double norm_ei2 = 0.0;
		for (size_t i = m_position_biggest_ft; i < n_damage - 1; i++)
		{
			this->ReconstructStrain(ei1, radius[i], Theta);
			this->ReconstructStrain(ei2, radius[i + 1], Theta);
			norm_ei1 = sqrt(ei1[0] * ei1[0] + ei1[1] * ei1[1] + 2.0 * ((0.5*ei1[2]) * (0.5*ei1[2])));
			norm_ei2 = sqrt(ei2[0] * ei2[0] + ei2[1] * ei2[1] + 2.0 * ((0.5*ei2[2]) * (0.5*ei2[2])));
			gf += std::abs(0.5 * (norm_ei2 - norm_ei1)*(norm_s[i] + norm_s[i + 1]));
		}
	}

	void InterpolatedConstitutiveLaw2D::CalculateAlpha(const double& lch_macro)
	{
		mErrorCode = 0.0;
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
		mAlpha = ((Gf - mG0) / (mGf_bar - mG0)) - 1.0;
		//std::cout << "mult: " << mult << std::endl;
		//std::cout << "mG0: " << mG0 << std::endl;
		//std::cout << "mGf_bar: " << mGf_bar << std::endl;
		//std::cout << "Gf: " << Gf << std::endl;
		//std::cout << "mAlpha: " << mAlpha << std::endl;
		if (mAlpha <= -1.0)
		{
			std::stringstream ss;
			ss << "RvePredictorCalculator Error: mAlphaxx <= -1.0" << std::endl;
			ss << "mAlpha = " << mAlpha << std::endl;
			ss << "mG0 = " << mG0 << std::endl;
			ss << "mGf_bar = " << mGf_bar << std::endl;
			ss << "m_position_biggest_ft = " << m_position_biggest_ft << std::endl;
			std::cout << ss.str();
			mErrorCode = -1.0;
			//exit(-1);
		}
	}

	void InterpolatedConstitutiveLaw2D::RegularizeInterpRadius(Vector& InterpRadius, const Vector& Theta)
	{
		size_t n_radius = InterpRadius.size();
		for (size_t i = m_position_biggest_ft; i < n_radius; i++)
			//for (size_t i = 0; i < n_radius; i++)
		{
			// eps[i] = eps[i] + mAlpha * (eps[i] - eps[0]);
			// ej = ej+(ej-ep)*S;

			Vector e(3, false);
			e(0) = InterpRadius[i] * cos(Theta[0]) + mAlpha * (InterpRadius[i] * cos(Theta[0]) - InterpRadius[m_position_biggest_ft] * cos(Theta[0]));
			e(1) = InterpRadius[i] * sin(Theta[0])*cos(Theta[1]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*cos(Theta[1]) - InterpRadius[m_position_biggest_ft] * sin(Theta[0])*cos(Theta[1]));
			e(2) = InterpRadius[i] * sin(Theta[0])*sin(Theta[1]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*sin(Theta[1]) - InterpRadius[m_position_biggest_ft] * sin(Theta[0])*sin(Theta[1]));

			Matrix S(3, 3, false);
			S(0, 0) = 1.0;
			S(1, 1) = 1.0;
			S(2, 2) = 1.0;
			Vector Se = prod(S, e);
			double eTSe = e(0)*Se(0) + e(1)*Se(1) + e(2)*Se(2);
			InterpRadius[i] = sqrt(eTSe);

			//InterpRadius[i] = InterpRadius[i] + mAlpha*(InterpRadius[i] - InterpRadius[m_position_biggest_ft]);
		}
	}

	void InterpolatedConstitutiveLaw2D::CalculateConstitutiveMatrix(const Properties& props,
		const Vector& strainVector,
		const Vector& stressVector,
		Matrix& constitutiveMatrix)
	{
		if (props.Has(DAMAGE_SECANT_MATRIX)) {
			if (props[DAMAGE_SECANT_MATRIX] == 1) {
				noalias(constitutiveMatrix) = (1 - mDamage)*mC0;
				return;
			}
			else
			{
				noalias(constitutiveMatrix) = (1 - mDamageConverged)*mC0;
				return;
			}
		}

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
			//h = std::max(1.0e-9, 1.0e-6*std::abs(strainVector(j)));

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

	Vector InterpolatedConstitutiveLaw2D::CalculateTheta(const Vector& StrainVector) const {
		//Calculate Theta (for the real strain)
		size_t theta_size = StrainVector.size() - 1;
		//std::cout << "PredictStress2D -> StrainVector: " << StrainVector << std::endl;
		Vector Theta(theta_size, 0.0);
		//Calculate Theta1
		double theta_1;
		double denom_theta_1 = sqrt(pow(StrainVector[2], 2) + pow(StrainVector[1], 2) + pow(StrainVector[0], 2));
		//std::cout << "PredictStress2D -> denom_theta_1: " << denom_theta_1 << std::endl;
		if (denom_theta_1 == 0.0)
		{
			theta_1 = 0.0;
		}
		else
		{
			theta_1 = acos(StrainVector[0] / denom_theta_1);
		}
		if (abs(theta_1) > Globals::Pi)
		{
			double mult_1 = abs(floor(theta_1 / Globals::Pi));
			if (theta_1 > Globals::Pi)
				theta_1 = theta_1 - 2.0 * Globals::Pi*mult_1;
			else
				theta_1 = theta_1 + 2.0 * Globals::Pi*mult_1;
		}

		//Calculate Theta2
		double theta_2;
		double denom_theta_2 = sqrt(pow(StrainVector[2], 2) + pow(StrainVector[1], 2));
		//std::cout << "PredictStress2D -> denom_theta_2: " << denom_theta_2 << std::endl;
		if (denom_theta_2 == 0.0)
		{
			theta_2 = 0.0;
		}
		else
		{
			// theta_2 = atan(StrainVector[2] / StrainVector[1]);
			if (StrainVector[2] < 0.0)
			{
				//std::cout << "PredictStress2D -> StrainVector[2] < 0.0" << std::endl;
				//std::cout << "PredictStress2D -> 2.0 * Globals::Pi" << 2.0 * Globals::Pi << std::endl;
				//std::cout << "PredictStress2D -> StrainVector[1] / denom_theta_2" << StrainVector[1] / denom_theta_2 << std::endl;
				//std::cout << "PredictStress2D -> acos(StrainVector[1] / denom_theta_2)" << acos(StrainVector[1] / denom_theta_2) << std::endl;
				theta_2 = 2.0 * Globals::Pi - acos(StrainVector[1] / denom_theta_2);
			}
			else
			{
				//std::cout << "PredictStress2D -> StrainVector[2] > 0.0" << std::endl;
				theta_2 = acos(StrainVector[1] / denom_theta_2);
			}
		}
		if (abs(theta_2) > Globals::Pi)
		{
			double mult_1 = abs(floor(theta_2 / Globals::Pi));
			if (theta_2 > Globals::Pi)
				theta_2 = theta_2 - 2.0 * Globals::Pi*mult_1;
			else
				theta_2 = theta_2 + 2.0 * Globals::Pi*mult_1;
		}

		//Assemble Theta
		Theta[0] = theta_1;
		Theta[1] = theta_2;

		return Theta;
	}

	void InterpolatedConstitutiveLaw2D::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set(PLANE_STRESS_LAW);
		rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
		rFeatures.mOptions.Set(ISOTROPIC);

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the space dimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}

	int InterpolatedConstitutiveLaw2D::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
			if (!rMaterialProperties.Has(LCH_REF_RVE)) {
				KRATOS_THROW_ERROR(std::logic_error, "InterpolatedConstitutiveLaw2D - missing LCH_REF_RVE", "");
			}
			if (!rMaterialProperties.Has(THICKNESS)) {
				KRATOS_THROW_ERROR(std::logic_error, "InterpolatedConstitutiveLaw2D - missing THICKNESS", "");
			}
		return 0;

		KRATOS_CATCH("");
	}

} /* namespace Kratos.*/
