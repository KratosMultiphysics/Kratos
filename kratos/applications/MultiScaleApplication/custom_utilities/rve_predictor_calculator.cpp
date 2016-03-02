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

#include <list>
#include <limits>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "multiscale_application.h"
#include "rve_predictor_calculator.h"
#include "rve_utilities.h"
#include "math_helpers.h"

#include "includes/kratos_parameters.h"

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <boost/unordered_map.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <boost/lexical_cast.hpp>
#include <utility>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // !M_PI

namespace Kratos
{

	inline double round(double d)
	{
	  return floor(d + 0.5);
	}
	
	RvePredictorCalculator::RvePredictorCalculator(std::string& file_name, ModelPart& modelPart)
		: mFileName(file_name)
		, mC0()
		, mStrainMap()
		, m_dimension(0)
		, m_domain_size(0.0)
	{
		std::cout << "Costruttore RvePredictorCalculator" << std::endl;
		mStrainMap.reserve(sizeof(double)*5*3000);
		std::string hash_file_name = mFileName; // "E:/PC_UNI/Paper1and2/TestSurf/RveSurf/Restored_d0.1_25.txt";
		std::cout << "hash_file_name: " << hash_file_name << std::endl;
		this->ReadAndCreatePredictorData(hash_file_name, modelPart);

		//std::cout << "mStrainMap contains:" << std::endl;
		//for (auto& x : mStrainMap) //only for msvc 11 and newest
		//	std::cout << x.first << ": " << x.second << std::endl;
		//std::cout << std::endl;
	}

	RvePredictorCalculator::~RvePredictorCalculator()
	{
	}

	void RvePredictorCalculator::FindDimension(ModelPart& modelPart)
	{
		double zmin = std::numeric_limits<double>::max();
		double zmax = -zmin;
		for (ModelPart::NodeIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
		{
			double iz = it->Z0();
			zmin = std::min(zmin, iz);
			zmax = std::max(zmax, iz);
		}
		double tolerance = RveUtilities::Precision()*std::max(zmin, zmax);
		if (std::abs(zmax - zmin) <= tolerance)
			m_dimension = 2;
		else
			m_dimension = 3;
	}

	void RvePredictorCalculator::FindDomainSize(ModelPart& modelPart)
	{
		m_domain_size = 0.0;
		if (m_dimension == 2)
		{
			for (ModelPart::ElementIterator it = modelPart.ElementsBegin(); it != modelPart.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				m_domain_size += igeom.Area();
			}
		}
		else if (m_dimension == 3)
		{
			for (ModelPart::ElementIterator it = modelPart.ElementsBegin(); it != modelPart.ElementsEnd(); ++it)
			{
				Element& ielem = *it;
				Element::GeometryType& igeom = ielem.GetGeometry();
				m_domain_size += igeom.Volume();
			}
		}
	}

	bool RvePredictorCalculator::Has(const Variable<double>& rThisVariable)
	{
		return true;
	}

	bool RvePredictorCalculator::Has(const Variable<Vector>& rThisVariable)
	{
		return true;
	}

	bool RvePredictorCalculator::Has(const Variable<Matrix>& rThisVariable)
	{
		return true;
	}

	bool RvePredictorCalculator::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
	{
		return true;
	}

	bool RvePredictorCalculator::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
	{
		return true;
	}

	double& RvePredictorCalculator::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		return rValue;
	}

	Vector& RvePredictorCalculator::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
	{
		return rValue;
	}

	Matrix& RvePredictorCalculator::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
	{
		return rValue;
	}

	array_1d<double, 3 > & RvePredictorCalculator::GetValue(const Variable<array_1d<double, 3 > >& rVariable,
		array_1d<double, 3 > & rValue)
	{
		return rValue;
	}

	array_1d<double, 6 > & RvePredictorCalculator::GetValue(const Variable<array_1d<double, 6 > >& rVariable,
		array_1d<double, 6 > & rValue)
	{
		return rValue;
	}

	void RvePredictorCalculator::SetValue(const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void RvePredictorCalculator::SetValue(const Variable<Vector >& rVariable,
		const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void RvePredictorCalculator::SetValue(const Variable<Matrix >& rVariable,
		const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void RvePredictorCalculator::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void RvePredictorCalculator::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void RvePredictorCalculator::ReadAndCreatePredictorData(std::string file_name, ModelPart& modelPart)
	{
		size_t strain_size = 3; // Modify it with FindDimension()

		std::cout << "ReadAndCreatePredictorData: " << file_name << std::endl;

		if (strain_size == 3) {
			if (mC0.size1() != strain_size || mC0.size2() != strain_size)
				mC0.resize(strain_size, strain_size, false);
			noalias(mC0) = ZeroMatrix(strain_size, strain_size);

			bool UseJson = false;
			if (UseJson)
			{
				//Read the Tags and the Strains From JsonFile
				/* JSON PRETTY FORMAT
				/ {
				/ 	"[-23, -13]": {									/ 1.		TAG
				/ 		"[0.5]": {									/ 1.1		DAMAGE
				/			"lambda": 0.02066451529793182,			/ 1.1.1		LAMBDA -> STRAIN INTENSITY
				/ 			"stress": [-9.191925896294595,			/ 1.1.2		STRESS
				/ 			-1.233581037907645,
				/ 			2.0069675053082023],
				/ 		}
				/		add mode damage ...
				/ 	}
				/	add more tags ...
				/ }
				*/
				std::cout << "Reading Json File" << std::endl;
				KratosParameters kp(file_name.c_str());

				// Fill HashMap
				int tag_1 = -23;	int tag_2 = -13;
				std::string str_tag_1 = boost::lexical_cast<std::string>(tag_1);
				std::string str_tag_2 = boost::lexical_cast<std::string>(tag_2);
				std::string tag_key = "[" + str_tag_1 + ", " + str_tag_2 + "]";

				system("pause");

				Parameters a = kp[tag_key];
				kp[tag_key]["stress"][2].GetDouble();
			}
			else
			{
				std::ifstream hashfile(file_name.c_str());
				std::string line;
				int count_line = 0;
				while (std::getline(hashfile, line)) //loop over lines of the .hash
				{
					std::istringstream iss(line);
					//std::cout << "IN-Line: " << line << std::endl;
					// Define variables
					int tag_1, tag_2;
					char comma;
					char bracket;
					double damage;
					double e_xx, e_yy, e_xy;
					double s_xx, s_yy, s_xy;
					double C11, C12, C13;
					double C21, C22, C23;
					double C31, C32, C33;

					// Read the File
					if (!(iss >> bracket >> damage >> bracket >> tag_1 >> comma >> tag_2 >> e_xx >> e_yy >> e_xy >> s_xx >> s_yy >> s_xy >> bracket >> C11 >> comma >> C12 >> comma >> C13 >> comma >> C21 >> comma >> C22 >> comma >> C23 >> comma >> C31 >> comma >> C32 >> comma >> C33 >> bracket))
					{
						std::cout << "damage: " << damage << std::endl;
						std::cout << "tag_1,tag_2: " << tag_1 << ", " << tag_2 << std::endl;
						std::cout << "e_xx, e_yy, e_xy: " << e_xx << ", " << e_yy << ", " << e_xy << std::endl;
						std::cout << "s_xx, s_yy, s_xy: " << s_xx << ", " << s_yy << ", " << s_xy << std::endl;

						std::cout << "C11, C12, C13: " << C11 << ", " << C12 << ", " << C13 << std::endl;
						std::cout << "C21, C22, C23: " << C21 << ", " << C22 << ", " << C23 << std::endl;
						std::cout << "C31, C32, C33: " << C31 << ", " << C32 << ", " << C33 << std::endl;
						std::cout << "Error during Reading the Hash File" << std::endl;
						break; // error
					}

					if (count_line == 0)
					{
						mC0(0, 0) = C11; mC0(0, 1) = C12; mC0(0, 2) = C13;
						mC0(1, 0) = C21; mC0(1, 1) = C22; mC0(1, 2) = C23;
						mC0(2, 0) = C31; mC0(2, 1) = C32; mC0(2, 2) = C33;
						count_line++;
					}

					// Fill variables
					vector<int> i_tag(2);
					i_tag[0] = tag_1;
					i_tag[1] = tag_2;

					vector<double> r_and_strain(4);
					vector<double> i_strain(3);
					i_strain[0] = e_xx;
					i_strain[1] = e_yy;
					i_strain[2] = e_xy;

					double radius = norm_2(i_strain);
					r_and_strain[0] = radius;
					r_and_strain[1] = e_xx;
					r_and_strain[2] = e_yy;
					r_and_strain[3] = e_xy;

					mStrainMap.insert(TagStrainVectorMap::value_type(i_tag, r_and_strain));
				}
			}
		}
		else
		{
			if (file_name.substr(file_name.find_last_of(".") + 1) == "txt")
			{
				if (mC0.size1() != strain_size || mC0.size2() != strain_size)
					mC0.resize(strain_size, strain_size, false);
				noalias(mC0) = ZeroMatrix(strain_size, strain_size);

				//Read the ConstitutiveMatrix From file
				std::ifstream mC0file(file_name.c_str());
				std::string line;
				while (std::getline(mC0file, line)) //loop over lines of the .hash
				{
					std::istringstream iss(line);
					char comma;
					char bracket;
					double C11, C12, C13, C14, C15, C16;
					double C21, C22, C23, C24, C25, C26;
					double C31, C32, C33, C34, C35, C36;
					double C41, C42, C43, C44, C45, C46;
					double C51, C52, C53, C54, C55, C56;
					double C61, C62, C63, C64, C65, C66;
					if (!(iss >> bracket >> C11 >> comma >> C12 >> comma >> C13 >> comma >> C14 >> comma >> C15 >> comma >> C16 >>
						comma >> C21 >> comma >> C22 >> comma >> C23 >> comma >> C24 >> comma >> C25 >> comma >> C26 >>
						comma >> C31 >> comma >> C32 >> comma >> C33 >> comma >> C34 >> comma >> C35 >> comma >> C36 >>
						comma >> C41 >> comma >> C42 >> comma >> C43 >> comma >> C44 >> comma >> C45 >> comma >> C46 >>
						comma >> C51 >> comma >> C52 >> comma >> C53 >> comma >> C54 >> comma >> C55 >> comma >> C56 >>
						comma >> C61 >> comma >> C62 >> comma >> C63 >> comma >> C64 >> comma >> C65 >> comma >> C66 >> bracket))
					{
						std::cout << "Error during Reading the C0 File" << std::endl;
						break; // error
					}
					mC0(0, 0) = C11; mC0(0, 1) = C12; mC0(0, 2) = C13; mC0(0, 3) = C14; mC0(0, 4) = C15; mC0(0, 5) = C16;
					mC0(1, 0) = C21; mC0(1, 1) = C22; mC0(1, 2) = C23; mC0(1, 3) = C24; mC0(1, 4) = C25; mC0(1, 5) = C26;
					mC0(2, 0) = C31; mC0(2, 1) = C32; mC0(2, 2) = C33; mC0(2, 3) = C34; mC0(2, 4) = C35; mC0(2, 5) = C36;
					mC0(3, 0) = C41; mC0(3, 1) = C42; mC0(3, 2) = C43; mC0(3, 3) = C44; mC0(3, 4) = C45; mC0(3, 5) = C46;
					mC0(4, 0) = C51; mC0(4, 1) = C52; mC0(4, 2) = C53; mC0(4, 3) = C54; mC0(4, 4) = C55; mC0(4, 5) = C56;
					mC0(5, 0) = C61; mC0(5, 1) = C62; mC0(5, 2) = C63; mC0(5, 3) = C64; mC0(5, 4) = C65; mC0(5, 5) = C66;
				}
			}
			else
			{
				//Read the Tags and the Strains From file
				std::ifstream hashfile(file_name.c_str());
				std::string line;
				while (std::getline(hashfile, line)) //loop over lines of the .hash
				{
					std::istringstream iss(line);

					// Define variables
					int tag_1, tag_2, tag_3, tag_4, tag_5;
					char comma;
					double e_xx, e_yy, e_zz, e_xy, e_yz, e_xz;

					// Read the File
					if (!(iss >> tag_1 >> comma >> tag_2 >> comma >> tag_3 >> comma >> tag_4 >> comma >> tag_5 >> e_xx >> e_yy >> e_zz >> e_xy >> e_yz >> e_xz))
					{
						std::cout << "Error during Reading the Hash File" << std::endl;
						break; // error
					}

					// Fill variables
					vector<int> i_tag(2);
					i_tag[0] = tag_1;
					i_tag[1] = tag_2;
					i_tag[2] = tag_3;
					i_tag[3] = tag_4;
					i_tag[4] = tag_5;

					vector<double> r_and_strain(4);
					vector<double> i_strain(3);
					i_strain[0] = e_xx;
					i_strain[1] = e_yy;
					i_strain[2] = e_zz;
					i_strain[3] = e_xy;
					i_strain[4] = e_yz;
					i_strain[5] = e_xz;

					double radius = norm_2(i_strain);
					r_and_strain[0] = radius;
					r_and_strain[1] = e_xx;
					r_and_strain[2] = e_yy;
					r_and_strain[3] = e_zz;
					r_and_strain[4] = e_xy;
					r_and_strain[5] = e_yz;
					r_and_strain[6] = e_xz;

					// Fill HashMap
					mStrainMap.insert(TagStrainVectorMap::value_type(i_tag, r_and_strain));
				}
			}
		}
		//modelPart.GetProcessInfo()[RVE_CLAW_MAP] = mStrainMap;
	}

	//bool RvePredictorCalculator::Predict(const TagStrainVectorMap&, const Vector& trial_macro_scale_strain_vector) const
	bool RvePredictorCalculator::Predict(const Vector& trial_macro_scale_strain_vector) const
	{
		std::cout << "IN - Predict\n";
		bool Is_Elastic = false;
		double strain_size = trial_macro_scale_strain_vector.size();

		// Check if the analysis is Thermal or Mechanical
		if (strain_size == 3) //Only 2D
		{
			std::cout << "if (strain_size == 3) //Only 2D mStrainMap.size()-> " << mStrainMap.size() << std::endl;
			//The Analysis is Mechanical
			if (mStrainMap.empty()) //Check if strain_map it is not empty
			{
				std::stringstream ss;
				ss << "RvePredictorCalculator: STRAIN MAP IS EMPTY" << std::endl;
				ss << "RvePredictorCalculator: Cannot calculate the material response without a micromodel history" << std::endl;
				ss << "RvePredictorCalculator: The Rve Will be Generated" << std::endl;
				std::cout << ss.str();
				//KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - Cannot calculate the material response without a micromodel history", "");
				Is_Elastic = false;
			}
			else
			{
				//Initialize data for the interpolation
				Vector Theta(strain_size - 1);
				Vector eps(strain_size);
				double radius = 0.0;

				//Calculate phi & number of division
				double m = (sqrtl(mStrainMap.size()) - 1) / 2; // m = max of tag
				//std::cout << "m: " << m << std::endl;
				double phi = M_PI / m;
				//std::cout << "phi: " << phi << " rad" << std::endl;

				//1. Calculate the normalized Strain
				double norm_real_eps = norm_2(trial_macro_scale_strain_vector);
				//std::cout << "norm_real_eps: " << norm_real_eps << std::endl;
				Vector normalized_eps = norm_real_eps > 0 ? trial_macro_scale_strain_vector / norm_real_eps : trial_macro_scale_strain_vector;//This is the value to compare with the RADIUS
				//std::cout << "normalized_eps: " << normalized_eps << std::endl;

				Theta = this->CalculateTheta(normalized_eps);
				//std::cout << "Theta: " << Theta << std::endl;

				Vector real_tag = Theta / phi;
				//std::cout << "real_tag: " << real_tag << std::endl;

				bool first_interpolation = false;
				bool second_interpolation = false;
				bool third_interpolation = false;
				bool BiCubicInterpolation = true;
				bool RadialBasis = false;
				bool nCubicInterpolation = false;
				if (first_interpolation)
				{
					//FIRST INTERPOLATION - NEAREST POINT
					vector<int> round_tag(2);
					round_tag[0] = round(real_tag[0]);
					round_tag[1] = round(real_tag[1]);

					//2. Read the corresponding radius and strain associated to the tag
					TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(round_tag);
					if (it_strain == mStrainMap.end())
					{
						std::cout << "no strain in map" << std::endl;
						exit(-1);
					}
					const Vector & i_r_and_strain = it_strain->second;
					radius = i_r_and_strain[0];
					eps[0] = i_r_and_strain[1];
					eps[1] = i_r_and_strain[2];
					eps[2] = i_r_and_strain[3];

					std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
					if (norm_real_eps < radius)
						Is_Elastic = true;
					else
						Is_Elastic = false;
				}
				else if (second_interpolation)
				{
					//SECOND INTERPOLATION - BI-LINEAR OVER TAG
					//1. Found the minimum value of the tag
					vector<int> min_tag(2);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "no strain in map" << std::endl;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
						eps[0] = i_r_and_strain[1];
						eps[1] = i_r_and_strain[2];
						eps[2] = i_r_and_strain[3];
					}
					else
					{
						double local_x = real_tag[0] - min_tag[0];
						double local_y = real_tag[1] - min_tag[1];

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
								N_i = 1.0 / 4.0 * (1.0 + xi[i] * (local_x * 2.0 - 1.0))*(1.0 + eta[i] * (local_y * 2.0 - 1.0));
								N_sum += N_i;
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk;
								tag_i[1] = min_tag[1] + ll;
								if (abs(tag_i(0)) > m)
								{
									if (tag_i(0) > m)
										tag_i(0) = tag_i(0) - m;
									else
										tag_i(0) = tag_i(0) + m;
								}
								else if (abs(tag_i(1)) > m)
								{
									if (tag_i(1) > m)
										tag_i(1) = tag_i(1) - m;
									else
										tag_i(1) = tag_i(1) + m;
								}

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "no strain in map" << std::endl;
									exit(-1);
								}
								const Vector & i_r_and_strain = it_strain->second;
								Vector i_eps(3);
								double i_radius = i_r_and_strain[0];

								radius += N_i*i_radius;

								i++;
							}
						}
						if (abs(1 - N_sum) > 1.0e-15)
						{
							KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
							Is_Elastic = false;
						}
					}

					std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
					if (norm_real_eps < radius)
						Is_Elastic = true;
					else
						Is_Elastic = false;
				}
				else if (third_interpolation)
				{
					std::cout << "THIRD INTERPOLATION - BI-LINEAR OVER EPS << BEST FITTING" << std::endl;
					//THIRD INTERPOLATION - BI-LINEAR OVER EPS << BEST FITTING
					//1. Found the minimum value of the tag
					vector<int> min_tag(2);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						std::cout << "real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "no strain in map" << std::endl;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
						eps[0] = i_r_and_strain[1];
						eps[1] = i_r_and_strain[2];
						eps[2] = i_r_and_strain[3];
					}
					else
					{
						//std::cout << "real_tag[0] != min_tag[0] && real_tag[1] != min_tag[1]" << std::endl;
						double local_x = real_tag[0] - min_tag[0];
						double local_y = real_tag[1] - min_tag[1];

						//std::cout << "real_tag: " << real_tag << std::endl;
						//std::cout << "min_tag: " << min_tag << std::endl;
						//std::cout << "local_x, local_y: " << local_x << ", " << local_y << std::endl;

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
								N_i = 1.0 / 4.0 * (1.0 + xi[i] * (local_x * 2.0 - 1.0))*(1.0 + eta[i] * (local_y * 2.0 - 1.0));
								N_sum += N_i;
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk;
								tag_i[1] = min_tag[1] + ll;
								if (abs(tag_i(0)) > m)
								{
									if (tag_i(0) > m)
										tag_i(0) = tag_i(0) - m;
									else
										tag_i(0) = tag_i(0) + m;
								}
								else if (abs(tag_i(1)) > m)
								{
									if (tag_i(1) > m)
										tag_i(1) = tag_i(1) - m;
									else
										tag_i(1) = tag_i(1) + m;
								}

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "no strain in map" << std::endl;
									exit(-1);
								}
								const Vector & i_r_and_strain = it_strain->second;
								Vector i_eps(3);
								double i_radius = i_r_and_strain[0];
								i_eps[0] = i_r_and_strain[1];
								i_eps[1] = i_r_and_strain[2];
								i_eps[2] = i_r_and_strain[3];

								double aux = N_i*i_radius;
								noalias(eps) += aux*i_eps;

								i++;
							}
						}
						if (abs(1 - N_sum) > 1.0e-15)
						{
							KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
							Is_Elastic = false;
						}
						radius = norm_2(eps);
					}

					std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
					if (norm_real_eps < radius)
						Is_Elastic = true;
					else
						Is_Elastic = false;
				}
				else if (BiCubicInterpolation)
				{
					//std::cout << "FOURTH INTERPOLATION - BI-CUBIC" << std::endl;
					//1. Found the minimum value of the tag
					vector<int> min_tag(2);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						//std::cout << "real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "no strain in map -> real_tag: " << real_tag << " min_tag: " << min_tag << std::endl;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
						eps[0] = i_r_and_strain[1];
						eps[1] = i_r_and_strain[2];
						eps[2] = i_r_and_strain[3];
					}
					else
					{
						// Create array
						double p[4][4];
						//array_1d< array_1d<double, 4>, 4 > p;

						double local_x = real_tag[0] - min_tag[0];
						double local_y = real_tag[1] - min_tag[1];

						size_t ii = 0;
						for (size_t kk = 0; kk < 4; kk++)
						{
							size_t jj = 0;
							for (size_t ll = 0; ll < 4; ll++)
							{
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk - 1;
								tag_i[1] = min_tag[1] + ll - 1;
								if (abs(tag_i(0)) > m)
								{
									if (tag_i(0) > m)
										tag_i(0) = tag_i(0) - m;
									else
										tag_i(0) = tag_i(0) + m;
								}
								else if (abs(tag_i(1)) > m)
								{
									if (tag_i(1) > m)
										tag_i(1) = tag_i(1) - m;
									else
										tag_i(1) = tag_i(1) + m;
								}

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									//std::cout << "no strain in map -> real_tag: " << real_tag << " min_tag: " << min_tag << std::endl;
									exit(-1);
								}
								const Vector & i_r_and_strain = it_strain->second;
								double i_radius = i_r_and_strain[0];
								p[ii][jj] = i_radius;
								//std::cout << "Vector of radius:" << p[ii][jj] << std::endl;
								jj++;
							}
							ii++;
						}
						size_t IType = 1;
						Matrix F = this->SelectInterpolationType(IType);
						radius = this->bicubicInterpolate(p, local_x, local_y, F);
						/*int n = strain_size - 1;
						double coordinates[2] = { local_x, local_y };
						std::cout << "n: " << n << " - coordinates: " << coordinates[0] << ", " << coordinates[1] << std::endl;
						radius = this->nCubicInterpolate(n, (double*)p, coordinates, F);*/
					}
					//std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
					if (norm_real_eps < radius)
						Is_Elastic = true;
					else
						Is_Elastic = false;
				}
				else if (RadialBasis)
				{
					std::cout << "Radial basis function - weight = norm(strain)" << std::endl;
					//1. Found the minimum value of the tag
					vector<int> min_tag(2);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						std::cout << "real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "no strain in map" << std::endl;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
						eps[0] = i_r_and_strain[1];
						eps[1] = i_r_and_strain[2];
						eps[2] = i_r_and_strain[3];
					}
					else
					{
						double summ_w_i = 0.0;
						double numeratore = 0.0;

						double local_x = real_tag[0] - min_tag[0];
						double local_y = real_tag[1] - min_tag[1];

						for (size_t kk = 0; kk < 4; kk++)
						{
							for (size_t ll = 0; ll < 4; ll++)
							{
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk - 1;
								tag_i[1] = min_tag[1] + ll - 1;
								if (abs(tag_i(0)) > m)
								{
									if (tag_i(0) > m)
										tag_i(0) = tag_i(0) - m;
									else
										tag_i(0) = tag_i(0) + m;
								}
								else if (abs(tag_i(1)) > m)
								{
									if (tag_i(1) > m)
										tag_i(1) = tag_i(1) - m;
									else
										tag_i(1) = tag_i(1) + m;
								}

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "no strain in map" << std::endl;
									exit(-1);
								}
								const Vector & i_r_and_strain = it_strain->second;
								double i_radius = i_r_and_strain[0];
								Vector i_eps(3);
								i_eps[0] = i_r_and_strain[1];
								i_eps[1] = i_r_and_strain[2];
								i_eps[2] = i_r_and_strain[3];
								Vector normalized_i_eps = i_eps / i_radius;
								double w_i = inner_prod(i_eps, normalized_eps);
								summ_w_i += w_i;
								numeratore += w_i * i_radius;
							}
						}
						radius = numeratore / summ_w_i;
					}
					std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
					if (norm_real_eps < radius)
						Is_Elastic = true;
					else
						Is_Elastic = false;
				}
				else if (nCubicInterpolation)
				{
					//1. Found the minimum value of the tag
					vector<int> min_tag(2);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						std::cout << "real_tag[:] == min_tag[:]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "no strain in map" << std::endl;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
						eps[0] = i_r_and_strain[1];
						eps[1] = i_r_and_strain[2];
					}
					else
					{
						// Create array
						double p[4][4];
						//array_1d< array_1d<double, 4>, 4 > p;

						double local_x = real_tag[0] - min_tag[0];
						double local_y = real_tag[1] - min_tag[1];

						size_t ii = 0;
						for (size_t kk = 0; kk < 4; kk++)
						{
							size_t jj = 0;
							for (size_t ll = 0; ll < 4; ll++)
							{
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk - 1;
								tag_i[1] = min_tag[1] + ll - 1;

								for (size_t i = 0; i < tag_i.size(); i++)
								{
									if (abs(tag_i(i)) > m)
									{
										if (tag_i(i) > m)
											tag_i(i) = tag_i(i) - m;
										else
											tag_i(i) = tag_i(i) + m;
									}
								}

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "no strain in map" << std::endl;
									exit(-1);
								}
								const Vector & i_r_and_strain = it_strain->second;
								double i_radius = i_r_and_strain[0];
								p[ii][jj] = i_radius;
								jj++;
							}
							ii++;
						}
						size_t IType = 1;
						Matrix F = this->SelectInterpolationType(IType);
						int n = strain_size - 1;
						double coordinates[2] = { local_x, local_y };
						radius = this->nCubicInterpolate(n, (double*)p, coordinates, F);
					}
					if (norm_real_eps < radius)
						Is_Elastic = true;
					else
						Is_Elastic = false;
				}
			}
		}
		else if (strain_size == 6)
		{
			//The Analysis is Mechanical
			if (mStrainMap.empty()) //Check if strain_map it is not empty
			{
				std::stringstream ss;
				ss << "RvePredictorCalculator: Cannot calculate the material response without a micromodel history" << std::endl;
				ss << "RvePredictorCalculator: The Rve Will be Generated" << std::endl;
				std::cout << ss.str();
				//KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - Cannot calculate the material response without a micromodel history", "");
				Is_Elastic = false;
			}
			else
			{
				//Initialize data for the interpolation
				Vector Theta(strain_size - 1);
				Vector eps(strain_size);
				double radius;

				//Calculate phi & number of division
				double m = (sqrtl(mStrainMap.size()) - 1) / 5;
				//std::cout << "m: " << m << std::endl;
				double phi = M_PI / m;

				//1. Calculate the normalized Strain
				double norm_real_eps = norm_2(trial_macro_scale_strain_vector);
				Vector normalized_eps = norm_real_eps > 0 ? trial_macro_scale_strain_vector / norm_real_eps : trial_macro_scale_strain_vector;//This is the value to compare with the RADIUS

				Theta = this->CalculateTheta(normalized_eps);


				//std::cout << "Theta: " << Theta << std::endl;

				Vector real_tag = Theta / phi;

				bool nCubicInterpolation = true;
				if (nCubicInterpolation)
				{
					//1. Found the minimum value of the tag
					vector<int> min_tag(5);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);
					min_tag[2] = floor(real_tag[2]);
					min_tag[3] = floor(real_tag[3]);
					min_tag[4] = floor(real_tag[4]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1] && real_tag[2] == min_tag[2] && real_tag[3] == min_tag[3] && real_tag[4] == min_tag[4])
					{
						std::cout << "real_tag[:] == min_tag[:]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "no strain in map" << std::endl;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
						eps[0] = i_r_and_strain[1];
						eps[1] = i_r_and_strain[2];
						eps[2] = i_r_and_strain[3];
						eps[3] = i_r_and_strain[4];
						eps[4] = i_r_and_strain[5];
						eps[5] = i_r_and_strain[6];
					}
					else
					{
						// Create array
						double p[4][4][4][4][4];
						//array_1d< array_1d< array_1d< array_1d< array_1d<double, 4>, 4>, 4>, 4>, 4> p;
						double local_x = real_tag[0] - min_tag[0];
						double local_y = real_tag[1] - min_tag[1];
						double local_z = real_tag[2] - min_tag[2];
						double local_l = real_tag[3] - min_tag[3];
						double local_m = real_tag[4] - min_tag[4];

						size_t ff = 0;
						for (size_t vv = 0; vv < 4; vv++)
						{
							size_t gg = 0;
							for (size_t mm = 0; mm < 4; mm++)
							{
								size_t hh = 0;
								for (size_t nn = 0; nn < 4; nn++)
								{
									size_t ii = 0;
									for (size_t kk = 0; kk < 4; kk++)
									{
										size_t jj = 0;
										for (size_t ll = 0; ll < 4; ll++)
										{
											vector<int> tag_i(5);
											tag_i[0] = min_tag[0] + vv - 1;
											tag_i[1] = min_tag[1] + mm - 1;
											tag_i[2] = min_tag[2] + nn - 1;
											tag_i[3] = min_tag[3] + kk - 1;
											tag_i[4] = min_tag[4] + ll - 1;

											for (size_t iii = 0; iii < tag_i.size(); iii++)
											{
												if (abs(tag_i(iii)) > m)
												{
													if (tag_i(iii) > m)
														tag_i(iii) = tag_i(iii) - m;
													else
														tag_i(iii) = tag_i(iii) + m;
												}
											}

											//2. Read the corresponding radius and strain associated to the tag
											TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
											if (it_strain == mStrainMap.end())
											{
												std::cout << "no strain in map" << std::endl;
												exit(-1);
											}
											const Vector & i_r_and_strain = it_strain->second;
											double i_radius = i_r_and_strain[0];
											p[ff][gg][hh][ii][jj] = i_radius;
											jj++;
										}
										ii++;
									}
									hh++;
								}
								gg++;
							}
							ff++;
						}
						size_t IType = 1;
						Matrix F = this->SelectInterpolationType(IType);
						int n = strain_size - 1;
						double coordinates[2] = { local_x, local_y };
						radius = this->nCubicInterpolate(n, (double*)p, coordinates, F);
					}
					if (norm_real_eps < radius)
						Is_Elastic = true;
					else
						Is_Elastic = false;
				}
			}
		}
		else
		{
			// The analysis is Thermal and always Elastic
			Is_Elastic = true;
		}
		return Is_Elastic;
	}

	Vector RvePredictorCalculator::CalculateTheta(Vector NormalizedStrainVector) const {
		//Calculate Theta (for the real strain)
		size_t theta_size = NormalizedStrainVector.size() - 1;
		Vector Theta(theta_size, 0.0);
		if (theta_size == 2)
		{
			//Calculate Theta1
			double theta_1;
			double denom_theta_1 = sqrt(pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2) + pow(NormalizedStrainVector[0], 2));
			if (denom_theta_1 == 0)
			{
				theta_1 = 0;
			}
			else
			{
				theta_1 = acos(NormalizedStrainVector[0] / denom_theta_1);
			}
			if (theta_1 > M_PI)
				theta_1 = theta_1 - 2 * M_PI;
			if (theta_1 < -M_PI)
				theta_1 = theta_1 + 2 * M_PI;

			//Calculate Theta2
			double theta_2;
			double denom_theta_2 = sqrt(pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2));
			if (denom_theta_2 == 0)
			{
				theta_2 = 0;
			}
			else
			{
				if (NormalizedStrainVector[2] < 0.0)
				{
					theta_2 = 2 * M_PI - acos(NormalizedStrainVector[1] / denom_theta_2);
				}
				else
				{
					theta_2 = acos(NormalizedStrainVector[1] / denom_theta_2);
				}
			}
			if (theta_2 > M_PI)
				theta_2 = theta_2 - 2 * M_PI;
			if (theta_2 < -M_PI)
				theta_2 = theta_2 + 2 * M_PI;

			//Assemble Theta
			Theta[0] = theta_1;
			Theta[1] = theta_2;
		}
		else
		{
			//Calculate Theta1
			double theta_1 = acos(NormalizedStrainVector[0] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2) + pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2) + pow(NormalizedStrainVector[0], 2)));
			if (theta_1 > M_PI)
				theta_1 = theta_1 - 2 * M_PI;
			if (theta_1 < -M_PI)
				theta_1 = theta_1 + 2 * M_PI;

			//Calculate Theta2
			double theta_2 = acos(NormalizedStrainVector[1] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2) + pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2)));
			if (theta_2 > M_PI)
				theta_2 = theta_2 - 2 * M_PI;
			if (theta_2 < -M_PI)
				theta_2 = theta_2 + 2 * M_PI;

			//Calculate Theta3
			double theta_3 = acos(NormalizedStrainVector[2] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2) + pow(NormalizedStrainVector[2], 2)));
			if (theta_3 > M_PI)
				theta_3 = theta_3 - 2 * M_PI;
			if (theta_3 < -M_PI)
				theta_3 = theta_3 + 2 * M_PI;

			//Calculate Theta4
			double theta_4 = acos(NormalizedStrainVector[2] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2)));
			if (theta_4 > M_PI)
				theta_4 = theta_4 - 2 * M_PI;
			if (theta_4 < -M_PI)
				theta_4 = theta_4 + 2 * M_PI;

			//Calculate Theta5
			double theta_5;
			double denom = sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2));
			if (denom == 0)
			{
				theta_5 = 0;
			}
			else
			{
				if (NormalizedStrainVector[5] < 0.0)
				{
					theta_5 = 2 * M_PI - acos(NormalizedStrainVector[4] / denom);
				}
				else
				{
					theta_5 = acos(NormalizedStrainVector[4] / denom);
				}
			}
			if (theta_5 > M_PI)
				theta_5 = theta_5 - 2 * M_PI;
			if (theta_5 < -M_PI)
				theta_5 = theta_5 + 2 * M_PI;

			//Assemble Theta
			Theta[0] = theta_1;
			Theta[1] = theta_2;
			Theta[2] = theta_3;
			Theta[3] = theta_4;
			Theta[4] = theta_5;
		}

		return Theta;
	}

	Matrix RvePredictorCalculator::SelectInterpolationType(size_t IType) const {
		// Now only Kochanek-Bartels was implemented for 2D and 3D!!

		/*% ' Basis matrix for global Tension, Continuity  & Bias
		% ' As T goes from 0 to 1 you go from n to n+1
		% '        n-1      n        n+1       n+2
		% '   T^3  -A    4+A-B-C   -4+B+C-D     D     /
		% '   T^2 +2A  -6-2A+2B+C  6-2B-C+D    -D    /
		% '   T^1  -A     A-B         B         0   /  2
		% '   T^0   0      2          0         0  /
		% '*/
		double FT = 0;  //%' Tension       T=+1-->Tight             T=-1--> Round
		double FB = 0;  //%' Bias          B=+1-->Post Shoot        B=-1--> Pre shoot
		double FC = -1; //%' Continuity    C=+1-->Inverted corners  C=-1--> Box corners
		/*% '
		% ' When T=B=C=0 this is the Catmul-Rom.
		% ' When T=1 & B=C=0 this is the Simple Cubic.
		% ' When T=B=0 & C=-1 this is the linear interp.
		% '
		% ' Here, the three parameters are folded into the basis matrix.  If you want
		% ' independent control of the segment start and end, make two T, C & Bs.
		% ' One for A & B (beginning) and one for C & D (end).  For local control of
		% ' each point, you'll need an array of T, C & Bs for each individual point.
		% ' If you want the local control as shown on the video or in the paper, you use
		% ' the "A" & "B" for the current segment and the "C" & "D" from the previous
		% ' segment.
		% '*/
		double FFA = (1 - FT)*(1 + FC)*(1 + FB);
		double FFB = (1 - FT)*(1 - FC)*(1 - FB);
		double FFC = (1 - FT)*(1 - FC)*(1 + FB);
		double FFD = (1 - FT)*(1 + FC)*(1 - FB);

		//% '  Here, the basis matrix coefficients are defined
		double DD = 2; //%' divisor for K-B
		double Ft1 = -FFA / DD, Ft2 = (+4 + FFA - FFB - FFC) / DD, Ft3 = (-4 + FFB + FFC - FFD) / DD, Ft4 = FFD / DD;
		double Fu1 = +2 * FFA / DD, Fu2 = (-6 - 2 * FFA + 2 * FFB + FFC) / DD, Fu3 = (+6 - 2 * FFB - FFC + FFD) / DD, Fu4 = -FFD / DD;
		double Fv1 = -FFA / DD, Fv2 = (FFA - FFB) / DD, Fv3 = FFB / DD, Fv4 = 0;
		double Fw1 = 0, Fw2 = +2 / DD, Fw3 = 0, Fw4 = 0;

		Matrix F(4, 4, 0.0);

		F(0, 0) = Ft1; F(0, 1) = Ft2; F(0, 2) = Ft3; F(0, 3) = Ft4;
		F(1, 0) = Fu1; F(1, 1) = Fu2; F(1, 2) = Fu3; F(1, 3) = Fu4;
		F(2, 0) = Fv1; F(2, 1) = Fv2; F(2, 2) = Fv3; F(2, 3) = Fv4;
		F(3, 0) = Fw1; F(3, 1) = Fw2; F(3, 2) = Fw3; F(3, 3) = Fw4;

		return F;
	}

	double RvePredictorCalculator::cubicInterpolate(double p[4], double x, Matrix& F) const {
		/*std::cout << ">>> cubicInterpolate \n";
		std::cout << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << '\n';*/
		//return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

		double FAX = F(0, 0)*p[0] + F(0, 1)*p[1] + F(0, 2)*p[2] + F(0, 3)*p[3];
		double FBX = F(1, 0)*p[0] + F(1, 1)*p[1] + F(1, 2)*p[2] + F(1, 3)*p[3];
		double FCX = F(2, 0)*p[0] + F(2, 1)*p[1] + F(2, 2)*p[2] + F(2, 3)*p[3];
		double FDX = F(3, 0)*p[0] + F(3, 1)*p[1] + F(3, 2)*p[2] + F(3, 3)*p[3];
		return ((FAX*x + FBX)*x + FCX)*x + FDX;
	}

	double RvePredictorCalculator::nCubicInterpolate(int n, double* p, double coordinates[], Matrix& F) const {
		assert(n > 0);
		if (n == 1) {
			return cubicInterpolate(p, *coordinates, F);
		}
		else {
			double arr[4];
			int skip = 1 << (n - 1) * 2;
			arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1, F);
			arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1, F);
			arr[2] = nCubicInterpolate(n - 1, p + 2 * skip, coordinates + 1, F);
			arr[3] = nCubicInterpolate(n - 1, p + 3 * skip, coordinates + 1, F);
			return cubicInterpolate(arr, *coordinates, F);
		}
	}

	double RvePredictorCalculator::bicubicInterpolate(double p[4][4], double x, double y, Matrix& F) const {
		/*std::cout << ">>> bicubicInterpolate \n";
		std::cout << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << '\n';*/
		double arr[4];
		arr[0] = cubicInterpolate(p[0], y, F);
		arr[1] = cubicInterpolate(p[1], y, F);
		arr[2] = cubicInterpolate(p[2], y, F);
		arr[3] = cubicInterpolate(p[3], y, F);
		return cubicInterpolate(arr, x, F);
	}

} // namespace Kratos
