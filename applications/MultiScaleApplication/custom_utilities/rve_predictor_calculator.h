/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Stefano Zaghi $
//   Date:                $Date: 2013-11-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_PREDICTOR_CALCULATOR_H_INCLUDED)
#define RVE_PREDICTOR_CALCULATOR_H_INCLUDED

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <boost/unordered_map.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <utility>

#include "rapidjson/document.h"

namespace Kratos
{

	class RvePredictorCalculator
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION(RvePredictorCalculator);

		struct KeyComparor
		{
			bool operator()(const vector<double>& lhs, const vector<double>& rhs) const
			{
				if (lhs.size() != rhs.size())
					return false;

				for (unsigned int i = 0; i<lhs.size(); i++)
				{
					if (lhs[i] != rhs[i]) return false;
				}

				return true;
			}
		};

		struct KeyHasher
		{
			std::size_t operator()(const vector<int>& k) const
			{
				return boost::hash_range(k.begin(), k.end());
			}
		};

		/*Unordered maps are associative containers that store elements
		* formed by the combination of a key value and a mapped value,
		* and which allows for fast retrieval of individual elements based on their keys.
		*
		* In an unordered_map, the key value is generally used to uniquely identify the element,
		* while the mapped value is an object with the content associated to this key.
		* Types of key and mapped value may differ.
		*
		* Internally, the elements in the unordered_map are not sorted in any particular order
		* with respect to either their key or mapped values, but organized into buckets depending
		* on their hash values to allow for fast access to individual elements directly by their
		* key values (with a constant average time complexity on average).
		*
		* unordered_map containers are faster than map containers to access individual elements
		* by their key, although they are generally less efficient for range iteration through
		* a subset of their elements.
		*
		* Unordered maps implement the direct access operator (operator[]) which allows for
		* direct access of the mapped value using its key value as argument.
		* Ref.: http://www.cplusplus.com/reference/unordered_map/unordered_map/
		*/
		typedef boost::unordered_map<vector<int>, Vector, KeyHasher, KeyComparor > TagStrainVectorMap;

		typedef std::map<double, vector<double> > DamageMap;

		typedef boost::unordered_map<vector<int>, DamageMap, KeyHasher, KeyComparor > TagNestedMap;

		typedef boost::unordered_map<vector<int>, Matrix, KeyHasher, KeyComparor> TagDamageMap;


	public:

		RvePredictorCalculator(std::string C0_file_name, std::string hash_file_name, std::string json_file_name);

		virtual ~RvePredictorCalculator();

	private:

		void ReadAndCreatePredictorData(std::string file_name);

		void ReadAndCreateConstTensorData(std::string file_name);

		void ReadAndCreatePredictorDataDamageHistory(std::string file_name);

	public:

		Matrix& GetTangentMatrix(Matrix& tangent_matrix);

		void PredictElasticity(const Vector& trial_macro_scale_strain_vector, bool& Is_Elastic, Vector& stress_vector) const;

		void PredictStress2D(const Vector& trial_macro_scale_strain_vector, Vector& stress_vector, double& lch_macro, Matrix& const_tens, double& EquivalentDamage, double& EquivalentDamageConverged);

		void ReconstructStrain(Vector& eps, const double& radius, const Vector& Theta);

		void CalculateFractureEnergy(const Vector& radius, const Vector& Theta, Vector& sigma_xx, Vector& sigma_yy, Vector& sigma_xy);

		void CalculateAlpha(double& lch_macro);

		void RegularizeInterpRadius(Vector& InterpRadius, Vector& Theta);

		void PredictStress3D(const Vector& trial_macro_scale_strain_vector, Vector& stress_vector, Matrix& const_tens, double& EquivalentDamage, double& EquivalentDamageConverged);

		double ShapeFunctionValue4N(size_t ShapeFunctionIndex, const Vector& rPoint) const;

		double ShapeFunctionValue8N(size_t ShapeFunctionIndex, const Vector& rPoint) const;

		double ShapeFunctionValue27N(size_t ShapeFunctionIndex, const Vector& rPoint) const;

		Vector CalculateTheta(Vector NormalizedStrainVector) const;

		Matrix SelectInterpolationType(size_t IType) const;

		//double cubicInterpolate(array_1d<double, 4> p, double x, Matrix& F) const;
		double Interpolate(array_1d<double, 4> p, double x, Matrix& F) const;

		double biDimInterpolate(array_1d<array_1d<double, 4>, 4> p, double x, double y, Matrix& F) const;

		double triDimInterpolate(array_1d<array_1d<array_1d<double, 4>, 4>, 4> p, double x, double y, double z, Matrix F) const;

		double quadDimInterpolate(array_1d<array_1d<array_1d<array_1d<double, 4>, 4>, 4>, 4> p, double x, double y, double z, double h, Matrix F) const;

		double pentaDimInterpolate(array_1d<array_1d<array_1d<array_1d<array_1d<double, 4>, 4>, 4>, 4>, 4> p, double x, double y, double z, double h, double k, Matrix F) const;

		double nDimInterpolate(size_t n, double* p, double coordinates[], Matrix& F) const;

	public:

		void FindDimension(ModelPart& modelPart);

		void FindDomainSize(ModelPart& modelPart);

		bool Has(const Variable<double>& rThisVariable);

		bool Has(const Variable<Vector>& rThisVariable);

		bool Has(const Variable<Matrix>& rThisVariable);

		bool Has(const Variable<array_1d<double, 3 > >& rThisVariable);

		bool Has(const Variable<array_1d<double, 6 > >& rThisVariable);

		double& GetValue(const Variable<double>& rThisVariable, double& rValue);

		Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue);

		Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue);

		array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rVariable,
			array_1d<double, 3 > & rValue);

		array_1d<double, 6 > & GetValue(const Variable<array_1d<double, 6 > >& rVariable,
			array_1d<double, 6 > & rValue);

		void SetValue(const Variable<double>& rVariable,
			const double& rValue,
			const ProcessInfo& rCurrentProcessInfo);

		void SetValue(const Variable<Vector >& rVariable,
			const Vector& rValue,
			const ProcessInfo& rCurrentProcessInfo);

		void SetValue(const Variable<Matrix >& rVariable,
			const Matrix& rValue,
			const ProcessInfo& rCurrentProcessInfo);

		void SetValue(const Variable<array_1d<double, 3 > >& rVariable,
			const array_1d<double, 3 > & rValue,
			const ProcessInfo& rCurrentProcessInfo);

		void SetValue(const Variable<array_1d<double, 6 > >& rVariable,
			const array_1d<double, 6 > & rValue,
			const ProcessInfo& rCurrentProcessInfo);

	private:

		std::string mC0FileName;
		std::string mHashFileName;
		std::string mJsonFileName;
		Matrix mC0;
		TagStrainVectorMap mStrainMap;
		TagNestedMap mTagNestedMap;
		DamageMap mDamageMap;
		TagDamageMap mTagDamageMap;
		rapidjson::Value* mValue;

		size_t m_position_biggest_ft;

		double mG0;
		double mGf_bar;
		double mAlpha;

	public:

		inline const size_t Dimension()const { return m_dimension; }

		inline const double DomainSize()const { return m_domain_size; }

	protected:

		size_t m_dimension;
		double m_domain_size;

	}; // class RvePredictor

}; // namespace Kratos

#endif // RVE_PREDICTOR_CALCULATOR_H_INCLUDED
