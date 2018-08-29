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

#include "multiscale_application_variables.h"
#include "rve_predictor_calculator.h"
#include "rve_utilities.h"
#include "math_helpers.h"

#include "includes/kratos_parameters.h"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <boost/unordered_map.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <boost/lexical_cast.hpp>
#include <utility>

namespace Kratos
{

#define SIGMA_SIGN(X) (X == 0.0 ? 1.0 : ( X > 0.0 ? 1.0 : -1.0 ))

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

	inline double round(double d)
	{
		return floor(d + 0.5);
	}

	inline double b3_eval_bezier(double xi,
		double x0, double x1, double x2,
		double y0, double y1, double y2)
	{
		double A = x0 - 2.0*x1 + x2;
		double B = 2.0*(x1 - x0);
		double C = x0 - xi;
		if (std::abs(A)<1.0E-12)
		{
			x1 = x1 + 1.0E-6*(x2 - x0);
			A = x0 - 2.0*x1 + x2;
			B = 2.0*(x1 - x0);
			C = x0 - xi;
		}
		double D = B*B - 4.0*A*C;
		double t = (-B + std::sqrt(D)) / (2.0*A);
		return (y0 - 2.0*y1 + y2)*t*t + 2.0*(y1 - y0)*t + y0;
	}

	RvePredictorCalculator::RvePredictorCalculator(std::string C0_file_name, std::string hash_file_name, std::string json_file_name)
		: mC0FileName(C0_file_name)
		, mHashFileName(hash_file_name)
		, mJsonFileName(json_file_name)
		, mC0()
		, mStrainMap()
		, mDamageMap()
		, mTagNestedMap()
		, mTagDamageMap()
		, m_dimension(0)
		, m_domain_size(0.0)
		, mG0(0.0)
		, mGf_bar(0.0)
		, mAlpha(0.0)
	{

		std::string mCo_file = mC0FileName; // "E:\\PC_UNI\\Paper1and2\\TestSurf\\RveSurf\\C0_Matrix_Damage2.txt";
		std::cout << "mCo_file: " << mCo_file << std::endl;
		this->ReadAndCreateConstTensorData(mCo_file);

		std::string hash_file = mHashFileName; // "E:/PC_UNI/Paper1and2/TestSurf/RveSurf/Restored_d0.1_25.txt"; //
		std::cout << "hash_file_name: " << hash_file << std::endl;
		this->ReadAndCreatePredictorData(hash_file);

		std::string json_file = mJsonFileName; // "C:\\Users\\zstefano\\Desktop\\json_surf_damage_m25.json";
		std::cout << "json_file: " << json_file << std::endl;
		this->ReadAndCreatePredictorDataDamageHistory(json_file);
	}

	RvePredictorCalculator::~RvePredictorCalculator()
	{
	}

	void RvePredictorCalculator::ReadAndCreatePredictorData(std::string file_name)
	{
		std::ifstream hashfile(file_name.c_str());
		std::string line;
		int count_line = 0;
		while (std::getline(hashfile, line)) //loop over lines of the .hash
		{
			std::istringstream iss(line);
			//std::cout << "IN-Line: " << line << std::endl;
			// Define variables
			int tag_1, tag_2, tag_3, tag_4, tag_5;
			char comma;
			char bracket;
			double damage;
			double e_xx, e_yy, e_zz, e_xy, e_yz, e_xz;
			double s_xx, s_yy, s_zz, s_xy, s_yz, s_xz;

			// Read the File
			m_dimension = 2;
			if (!(iss >> bracket >> damage >> bracket >> tag_1 >> comma >> tag_2 >> e_xx >> e_yy >> e_xy >> s_xx >> s_yy >> s_xy)) //2D
			{
				m_dimension = 3;
				if (!(iss >> bracket >> damage >> bracket >> tag_1 >> comma >> tag_2 >> comma >> tag_3 >> comma >> tag_4 >> comma >> tag_5 >> e_xx >> e_yy >> e_zz >> e_xy >> e_yz >> e_xz >> s_xx >> s_yy >> s_zz >> s_xy >> s_yz >> s_xz)) //3D
				{
					std::cout << "Error during Reading the Hash File" << std::endl;
					break; // error
				}
			}

			if (m_dimension == 2)
			{
				// Fill variables
				vector<int> i_tag(2);
				i_tag[0] = tag_1;
				i_tag[1] = tag_2;
				//std::cout << "i_tag: " << i_tag << std::endl;

				vector<double> r_and_strain_stress(7);

				vector<double> i_strain(3);
				i_strain[0] = e_xx;
				i_strain[1] = e_yy;
				i_strain[2] = 2.0*e_xy;

				vector<double> i_stress(3);
				i_stress[0] = s_xx;
				i_stress[1] = s_yy;
				i_stress[2] = s_xy;

				double radius = norm_2(i_strain);
				//double radius = sqrt(e_xx*e_xx + e_yy*e_yy + 2.0*(e_xy*e_xy));
				Matrix S(3, 3, false);
				S(0, 0) = 1.0;
				S(1, 1) = 1.0;
				S(2, 2) = 1.0;
				Vector Se = prod(S, i_strain);
				double eTSe = i_strain(0)*Se(0) + i_strain(1)*Se(1) + i_strain(2)*Se(2);
				radius = sqrt(eTSe);

				////Calculate phi & number of division
				//double m = 15.0; // m = max of tag
				////std::cout << "m: " << m << std::endl;
				//double phi = Globals::Pi / m;
				////std::cout << "phi: " << phi << " rad" << std::endl;
				//
				////std::cout << "norm_real_eps: " << norm_real_eps << std::endl;
				//Vector normalized_eps = radius > 0 ? i_strain / radius : i_strain;
				////std::cout << "normalized_eps: " << normalized_eps << std::endl;
				//Vector Theta(2, 0.0);
				//Theta = this->CalculateTheta(normalized_eps);
				////std::cout << "Theta: " << Theta << std::endl;
				//
				//Vector real_tag = Theta / phi;
				//std::cout << "real_tag: " << real_tag << std::endl;

				r_and_strain_stress[0] = radius;//This is the value to compare with the REAL RADIUS
				r_and_strain_stress[1] = e_xx;
				r_and_strain_stress[2] = e_yy;
				r_and_strain_stress[3] = 2.0 * e_xy; //multiply by 2.0 in order to obtain gamma xy
				r_and_strain_stress[4] = s_xx;
				r_and_strain_stress[5] = s_yy;
				r_and_strain_stress[6] = s_xy;

				mStrainMap.insert(TagStrainVectorMap::value_type(i_tag, r_and_strain_stress));
			}
			else
			{
				// Fill variables
				vector<int> i_tag(5);
				i_tag[0] = tag_1;
				i_tag[1] = tag_2;
				i_tag[2] = tag_3;
				i_tag[3] = tag_4;
				i_tag[4] = tag_5;

				vector<double> r_and_strain_stress(13);

				vector<double> i_strain(6);
				i_strain[0] = e_xx;
				i_strain[1] = e_yy;
				i_strain[2] = e_zz;
				i_strain[3] = 2.0*e_xy;
				i_strain[4] = 2.0*e_yz;
				i_strain[5] = 2.0*e_xz;

				vector<double> i_stress(6);
				i_stress[0] = s_xx;
				i_stress[1] = s_yy;
				i_stress[2] = s_zz;
				i_stress[3] = s_xy;
				i_stress[4] = s_yz;
				i_stress[5] = s_xz;

				double radius = norm_2(i_strain);
				r_and_strain_stress[0] = radius;
				r_and_strain_stress[1] = e_xx;
				r_and_strain_stress[2] = e_yy;
				r_and_strain_stress[3] = e_zz;
				r_and_strain_stress[4] = e_xy;
				r_and_strain_stress[5] = e_yz;
				r_and_strain_stress[6] = e_xz;
				r_and_strain_stress[7] = s_xx;
				r_and_strain_stress[8] = s_yy;
				r_and_strain_stress[9] = s_zz;
				r_and_strain_stress[10] = s_xy;
				r_and_strain_stress[11] = s_yz;
				r_and_strain_stress[12] = s_xz;

				mStrainMap.insert(TagStrainVectorMap::value_type(i_tag, r_and_strain_stress));
			}
		}
	}

	void RvePredictorCalculator::ReadAndCreateConstTensorData(std::string file_name)
	{

		std::ifstream hashfile(file_name.c_str());
		std::string line;
		int count_line = 0;
		while (std::getline(hashfile, line)) //loop over lines of the .hash
		{
			std::istringstream iss(line);
			//std::cout << "IN-Line: " << line << std::endl;
			// Define variables
			char comma;
			char bracket;
			double C11, C12, C13, C14, C15, C16;
			double C21, C22, C23, C24, C25, C26;
			double C31, C32, C33, C34, C35, C36;
			double C41, C42, C43, C44, C45, C46;
			double C51, C52, C53, C54, C55, C56;
			double C61, C62, C63, C64, C65, C66;

			// Read the File
			size_t strain_size = 3;
			if (!(iss >> bracket >> C11 >> comma >> C12 >> comma >> C13 >> comma >> C21 >> comma >> C22 >> comma >> C23 >> comma >> C31 >> comma >> C32 >> comma >> C33 >> bracket))
			{
				size_t strain_size = 6;
				if (!(iss >> bracket >> C11 >> comma >> C12 >> comma >> C13 >> comma >> C14 >> comma >> C15 >> comma >> C16 >> comma >> C21 >> comma >> C22 >> comma >> C23 >> comma >> C24 >> comma >> C25 >> comma >> C26 >> comma >> C31 >> comma >> C32 >> comma >> C33 >> comma >> C34 >> comma >> C35 >> comma >> C36 >> comma >> C41 >> comma >> C42 >> comma >> C43 >> comma >> C44 >> comma >> C45 >> comma >> C46 >> comma >> C51 >> comma >> C52 >> comma >> C53 >> comma >> C54 >> comma >> C55 >> comma >> C56 >> comma >> C61 >> comma >> C62 >> comma >> C63 >> comma >> C64 >> comma >> C65 >> comma >> C66 >> bracket))
				{
					std::cout << "Error during Reading the C0 File" << std::endl;
					break; // error
				}
			}

			if (mC0.size1() != strain_size || mC0.size2() != strain_size)
				mC0.resize(strain_size, strain_size, false);
			noalias(mC0) = ZeroMatrix(strain_size, strain_size);


			if (strain_size == 3)
			{
				mC0(0, 0) = C11; mC0(0, 1) = C12; mC0(0, 2) = C13;
				mC0(1, 0) = C21; mC0(1, 1) = C22; mC0(1, 2) = C23;
				mC0(2, 0) = C31; mC0(2, 1) = C32; mC0(2, 2) = C33;
			}
			else
			{
				mC0(0, 0) = C11; mC0(0, 1) = C12; mC0(0, 2) = C13; mC0(0, 3) = C14; mC0(0, 4) = C15; mC0(0, 5) = C16;
				mC0(1, 0) = C21; mC0(1, 1) = C22; mC0(1, 2) = C23; mC0(1, 3) = C24; mC0(1, 4) = C25; mC0(1, 5) = C26;
				mC0(2, 0) = C31; mC0(2, 1) = C32; mC0(2, 2) = C33; mC0(2, 3) = C34; mC0(2, 4) = C35; mC0(2, 5) = C36;
				mC0(3, 0) = C41; mC0(3, 1) = C42; mC0(3, 2) = C43; mC0(3, 3) = C44; mC0(3, 4) = C45; mC0(3, 5) = C46;
				mC0(4, 0) = C51; mC0(4, 1) = C52; mC0(4, 2) = C53; mC0(4, 3) = C54; mC0(4, 4) = C55; mC0(4, 5) = C56;
				mC0(5, 0) = C61; mC0(5, 1) = C62; mC0(5, 2) = C63; mC0(5, 3) = C64; mC0(5, 4) = C65; mC0(5, 5) = C66;
			}
		}

	}

	void RvePredictorCalculator::ReadAndCreatePredictorDataDamageHistory(std::string file_name)
	{
		//Read the Tags and the Strains From JsonFile
		/* JSON PRETTY FORMAT
		/ {
		/ 	"-23,-13": {									/ 1.		TAG
		/ 		"0.5": {									/ 1.1		DAMAGE
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

		/*
		std::cout << "Reading Json File" << std::endl;
		Parameters jsonFile("C:\\Users\\zstefano\\Desktop\\json_test.json");

		mValue = jsonFile.GetUnderlyingStorage(); // TODO: Change when iterator process is iplemented

		for (rapidjson::Value::ConstMemberIterator itr1 = mValue->MemberBegin(); itr1 != mValue->MemberEnd(); ++itr1) // Loop over TAGS
		{
		std::string i_tags_name = itr1->name.GetString();
		std::istringstream iss(i_tags_name);
		vector<int> i_tags;
		char comma;
		if (!(iss >> i_tags[0] >> comma >> i_tags[1]))
		std::cout << "ERROR DURING READING TAGS\n";

		rapidjson::Value* dValue = (jsonFile[i_tags_name.c_str()]).GetUnderlyingStorage();
		for (rapidjson::Value::ConstMemberIterator itr2 = dValue->MemberBegin(); itr2 != dValue->MemberEnd(); ++itr2) // Loop over DAMAGE
		{
		std::string i_damage_name = itr2->name.GetString();
		double i_damage = stod(i_damage_name); //convert string to double -> only c++11

		mTagNestedMap[i_tags][i_damage][0] = jsonFile[i_tags_name][i_damage_name]["lambda"].GetDouble();
		mTagNestedMap[i_tags][i_damage][1] = jsonFile[i_tags_name][i_damage_name]["stress"][0].GetDouble();
		mTagNestedMap[i_tags][i_damage][2] = jsonFile[i_tags_name][i_damage_name]["stress"][1].GetDouble();
		mTagNestedMap[i_tags][i_damage][3] = jsonFile[i_tags_name][i_damage_name]["stress"][2].GetDouble();
		}
		}
		*/

		FILE* pFile = fopen(file_name.c_str(), "rb");
		char buffer[65536];
		rapidjson::FileReadStream is(pFile, buffer, sizeof(buffer));
		rapidjson::Document document;
		rapidjson::ParseResult ok = document.ParseStream<0>(is);

		if (!ok)
		{
			std::stringstream msg;
			msg << rapidjson::GetParseError_En(ok.Code()) << " offset of the error from the beginning of the string = " << ok.Offset() << std::endl;
			msg << "a much more explicative error message can be obtained by analysing the input string with an online analyzer such for example json lint" << std::endl;
			msg << "the value of the string that was attempted to parse is :" << std::endl << std::endl;
			msg << file_name;
			KRATOS_THROW_ERROR(std::invalid_argument, "error found in parsing the json_string, the value of the json string was: \n", msg.str());
		}

		if (m_dimension == 2)
		{
			vector<int> i_tags(2);
			for (rapidjson::Value::ConstMemberIterator itr1 = document.MemberBegin(); itr1 != document.MemberEnd(); ++itr1) // Loop over TAGS
			{
				std::string i_tags_name = itr1->name.GetString();
				//std::cout << "i_tags_name: " << i_tags_name << std::endl;

				//convert string to int -> only c++11
				size_t sz;     // alias of size_t
				i_tags[0] = std::stoi(i_tags_name, &sz);
				i_tags[1] = std::stoi(i_tags_name.substr(sz));

				size_t i_D = 0;
				size_t itr_dim = document[i_tags_name.c_str()].Size();
				matrix<double> d_lambda_sigma(5, itr_dim, 0.0);
				for (rapidjson::Value::ConstMemberIterator itr2 = document[i_tags_name.c_str()].MemberBegin(); itr2 != document[i_tags_name.c_str()].MemberEnd(); ++itr2) // Loop over DAMAGE
				{
					std::string i_damage_name = itr2->name.GetString();
					//std::cout << "i_damage_name: " << i_damage_name << std::endl;
					double i_damage = stod(i_damage_name); //convert string to double -> only c++11

					vector<double> lambda_sigma(4, 0.0);
					lambda_sigma[0] = document[i_tags_name.c_str()][i_damage_name.c_str()]["lambda"].GetDouble();
					//std::cout << "lambda_sigma[0]: " << lambda_sigma[0] << std::endl;
					lambda_sigma[1] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][0].GetDouble();
					lambda_sigma[2] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][1].GetDouble();
					lambda_sigma[3] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][2].GetDouble();

					//mDamageMap.insert(DamageMap::value_type(i_damage, lambda_sigma));

					d_lambda_sigma(0, i_D) = i_damage;
					d_lambda_sigma(1, i_D) = lambda_sigma[0];
					d_lambda_sigma(2, i_D) = lambda_sigma[1];
					d_lambda_sigma(3, i_D) = lambda_sigma[2];
					d_lambda_sigma(4, i_D) = lambda_sigma[3];
					i_D += 1;
				}

				size_t n_damage = i_D;
				//std::cout << "n_damage: " << n_damage << std::endl;
				Vector damagesVector_aux(n_damage, 0.0);
				for (size_t ii = 0; ii < n_damage; ii++)
					damagesVector_aux[ii] = d_lambda_sigma(0, ii);
				size_t idx_final = 0;
				Matrix dls(5, itr_dim, 0.0);
				for (auto i_sort : sort_indexes(damagesVector_aux)) { // only c++11
					dls(0, idx_final) = d_lambda_sigma(0, i_sort);
					dls(1, idx_final) = d_lambda_sigma(1, i_sort);
					dls(2, idx_final) = d_lambda_sigma(2, i_sort);
					dls(3, idx_final) = d_lambda_sigma(3, i_sort);
					dls(4, idx_final) = d_lambda_sigma(4, i_sort);
					idx_final += 1;
				}

				mTagDamageMap.insert(TagDamageMap::value_type(i_tags, dls));
				//mTagNestedMap.insert(TagNestedMap::value_type(i_tags, mDamageMap));
			}
		}
		else
		{
			vector<int> i_tags(5);
			for (rapidjson::Value::ConstMemberIterator itr1 = document.MemberBegin(); itr1 != document.MemberEnd(); ++itr1) // Loop over TAGS
			{
				std::string i_tags_name = itr1->name.GetString();
				std::cout << "i_tags_name: " << i_tags_name << std::endl;

				//convert string to int -> only c++11
				size_t sz;
				i_tags[0] = std::stoi(i_tags_name, &sz);
				size_t pos0 = i_tags_name.find(i_tags[0]);

				i_tags[1] = std::stoi(i_tags_name.substr(pos0 + sz, sz));
				size_t pos1 = i_tags_name.find(i_tags[1]);

				i_tags[2] = std::stoi(i_tags_name.substr(pos1 + sz, sz));
				size_t pos2 = i_tags_name.find(i_tags[2]);

				i_tags[3] = std::stoi(i_tags_name.substr(pos2 + sz, sz));
				size_t pos3 = i_tags_name.find(i_tags[3]);

				i_tags[4] = std::stoi(i_tags_name.substr(pos3 + sz, sz));
				size_t pos4 = i_tags_name.find(i_tags[4]);

				std::cout << "i_tags: " << i_tags << std::endl;

				size_t i_D = 0;
				for (rapidjson::Value::ConstMemberIterator itr2 = document[i_tags_name.c_str()].MemberBegin(); itr2 != document[i_tags_name.c_str()].MemberEnd(); ++itr2) // Loop over DAMAGE
				{
					std::string i_damage_name = itr2->name.GetString();
					//std::cout << "i_damage_name: " << i_damage_name << std::endl;
					double i_damage = stod(i_damage_name); //convert string to double -> only c++11

					vector<double> lambda_sigma(7, 0.0);
					lambda_sigma[0] = document[i_tags_name.c_str()][i_damage_name.c_str()]["lambda"].GetDouble();
					//std::cout << "lambda_sigma[0]: " << lambda_sigma[0] << std::endl;
					lambda_sigma[1] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][0].GetDouble();
					lambda_sigma[2] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][1].GetDouble();
					lambda_sigma[3] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][2].GetDouble();
					lambda_sigma[4] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][3].GetDouble();
					lambda_sigma[5] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][4].GetDouble();
					lambda_sigma[6] = document[i_tags_name.c_str()][i_damage_name.c_str()]["stress"][5].GetDouble();

					mDamageMap.insert(DamageMap::value_type(i_damage, lambda_sigma));

					size_t itr_dim = document[i_tags_name.c_str()].Size();
					matrix<double> d_lambda_sigma(8, itr_dim);
					d_lambda_sigma(0, i_D) = i_damage;
					d_lambda_sigma(1, i_D) = lambda_sigma[0];
					d_lambda_sigma(2, i_D) = lambda_sigma[1];
					d_lambda_sigma(3, i_D) = lambda_sigma[2];
					d_lambda_sigma(4, i_D) = lambda_sigma[3];
					d_lambda_sigma(5, i_D) = lambda_sigma[4];
					d_lambda_sigma(6, i_D) = lambda_sigma[5];
					d_lambda_sigma(7, i_D) = lambda_sigma[6];

					mTagDamageMap.insert(TagDamageMap::value_type(i_tags, d_lambda_sigma));
					i_D += 1;
				}
				mTagNestedMap.insert(TagNestedMap::value_type(i_tags, mDamageMap));
			}
		}
	}

	Matrix& RvePredictorCalculator::GetTangentMatrix(Matrix& tangent_matrix)
	{
		noalias(tangent_matrix) = mC0;
		return tangent_matrix;
	}

	void RvePredictorCalculator::PredictElasticity(const Vector& trial_macro_scale_strain_vector, bool& Is_Elastic, Vector& stress_vector) const
	{
		//std::cout << "IN - PredictElasticity\n";
		double strain_size = trial_macro_scale_strain_vector.size();
		Vector InStrain = trial_macro_scale_strain_vector;

		// Check if the analysis is Thermal or Mechanical
		if (strain_size == 3) //Only 2D
		{
			//std::cout << "if (strain_size == 3) //Only 2D mStrainMap.size()-> " << mStrainMap.size() << std::endl;
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
				double phi = Globals::Pi / m;
				//std::cout << "phi: " << phi << " rad" << std::endl;

				//1. Calculate the normalized Strain
				double norm_real_eps = norm_2(InStrain);

				Matrix S(3, 3, false);
				S(0, 0) = 1.0;
				S(1, 1) = 1.0;
				S(2, 2) = 1.0;
				Vector Se = prod(S, InStrain);
				double eTSe = InStrain(0)*Se(0) + InStrain(1)*Se(1) + InStrain(2)*Se(2);
				norm_real_eps = sqrt(eTSe);

				//std::cout << "norm_real_eps2: " << norm_real_eps2 << std::endl;
				//double norm_real_eps = sqrt(InStrain[0] * InStrain[0] + InStrain[1] * InStrain[1] + 2.0*((InStrain[2] / 2.0) * (InStrain[2] / 2.0)));
				//std::cout << "norm_real_eps: " << norm_real_eps << std::endl;
				Vector normalized_eps = norm_real_eps > 0 ? InStrain / norm_real_eps : InStrain;//This is the value to compare with the RADIUS
																								//std::cout << "normalized_eps: " << normalized_eps << std::endl;

																								//Theta = this->CalculateTheta(normalized_eps);
				Theta = this->CalculateTheta(InStrain);
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

				bool first_interpolation = false;
				bool second_interpolation = true;
				bool third_interpolation = false;
				bool BiCubicInterpolation = false;
				bool RadialBasis = false;
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
						std::cout << "1. no strain in map, round_tag: " << round_tag << std::endl;
						exit(-1);
					}
					const Vector & i_r_and_strain = it_strain->second;
					radius = i_r_and_strain[0];

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
					//std::cout << "min_tag: " << min_tag << std::endl;
					vector<int> max_tag(2);
					max_tag[0] = ceil(real_tag[0]);
					max_tag[1] = ceil(real_tag[1]);
					//std::cout << "max_tag: " << max_tag << std::endl;

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "2. no strain in map, real_tag: " << real_tag << std::endl;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
					}
					else
					{
						// local_x and y defines where to estimate the value on the interpolated line,
						// it is 0 at the first point and 1 and the second point
						double local_x = (real_tag[0] - min_tag[0]); // / (max_tag[0] - min_tag[0]);
						double local_y = (real_tag[1] - min_tag[1]); // / (max_tag[1] - min_tag[1]);
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
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk;
								tag_i[1] = min_tag[1] + ll;
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

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "3. no strain in map, tag_i: " << tag_i << std::endl;
									exit(-1);
								}
								const Vector & i_r_and_strain = it_strain->second;
								double i_radius = i_r_and_strain[0];
								//std::cout << "i_radius: " << i_radius << std::endl;

								N_i = 1.0 / 4.0 * (1.0 + xi[i] * (local_x * 2.0 - 1.0))*(1.0 + eta[i] * (local_y * 2.0 - 1.0));
								N_sum += N_i;

								//std::cout << "N_i: " << N_i << std::endl;
								radius += N_i*i_radius;

								i++;
							}
						}
						if (abs(1 - N_sum) > 1.0e-15)
						{
							std::cout << "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION ON PredictElasticity - N_sum = " << N_sum << " ; shape_index = " << i << std::endl;
							KRATOS_THROW_ERROR(std::logic_error, "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
							Is_Elastic = false;
						}
					}

					//std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
					//std::cout << "difference: " << norm_real_eps - radius << std::endl;
					if (norm_real_eps < radius)
					{
						Is_Elastic = true;
						//std::cout << "Is_Elastic = true\n";
						noalias(stress_vector) = prod(mC0, trial_macro_scale_strain_vector);
					}
					else
					{
						Is_Elastic = false;
						//std::cout << "Is_Elastic = false\n";
						//std::cout << "AMMOUNT OF EXTRA STRAIN: " << norm_real_eps - radius << std::endl;
						noalias(stress_vector) = ZeroVector(strain_size);
					}
				}
				else if (third_interpolation)
				{
					//std::cout << "THIRD INTERPOLATION - BI-LINEAR OVER EPS << BEST FITTING" << std::endl;
					//THIRD INTERPOLATION - BI-LINEAR OVER EPS << BEST FITTING
					//1. Found the minimum value of the tag
					vector<int> min_tag(2);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);
					vector<int> max_tag(2);
					max_tag[0] = ceil(real_tag[0]);
					max_tag[1] = ceil(real_tag[1]);

					// local_x and y defines where to estimate the value on the interpolated line,
					// it is 0 at the first point and 1 and the second point
					double local_x = (real_tag[0] - min_tag[0]); // / (max_tag[0] - min_tag[0]);
					double local_y = (real_tag[1] - min_tag[1]); // / (max_tag[1] - min_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						//std::cout << "real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "4. no strain in map, real_tag: " << real_tag << std::endl;
							Is_Elastic = false;
							//return;
							exit(-1);
						}
						const Vector & i_r_and_strain = it_strain->second;
						radius = i_r_and_strain[0];
						eps[0] = i_r_and_strain[1];
						eps[1] = i_r_and_strain[2];
						eps[2] = i_r_and_strain[3];
						stress_vector[0] = i_r_and_strain[4];
						stress_vector[1] = i_r_and_strain[5];
						stress_vector[2] = i_r_and_strain[6];
					}
					else
					{
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

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "5. no strain in map, tag_i: " << tag_i << std::endl;
									Is_Elastic = false;
									//return;
									exit(-1);
								}
								const Vector & i_r_and_strain = it_strain->second;
								Vector i_eps(3);
								Vector i_sigma(3);
								double i_radius = i_r_and_strain[0];
								i_eps[0] = i_r_and_strain[1];
								i_eps[1] = i_r_and_strain[2];
								i_eps[2] = i_r_and_strain[3];

								i_sigma[0] = i_r_and_strain[4];
								i_sigma[1] = i_r_and_strain[5];
								i_sigma[2] = i_r_and_strain[6];

								double aux = N_i*i_radius;
								noalias(eps) += aux*i_eps;
								double aux2 = N_i*norm_2(i_sigma);
								noalias(stress_vector) += aux2*i_sigma;

								i++;
							}
						}
						if (abs(1 - N_sum) > 1.0e-15)
						{
							KRATOS_THROW_ERROR(std::logic_error, "RveAdapterV2 - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
							Is_Elastic = false;
						}
						//radius = norm_2(eps);
						radius = sqrt(eps[0] * eps[0] + eps[1] * eps[1] + 2.0*((eps[2] / 2.0)*(eps[2] / 2.0)));
					}

					//std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
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
					vector<int> max_tag(2);
					max_tag[0] = ceil(real_tag[0]);
					max_tag[1] = ceil(real_tag[1]);

					// local_x and y defines where to estimate the value on the interpolated line,
					// it is 0 at the first point and 1 and the second point
					double local_x = (real_tag[0] - min_tag[0]); // / (max_tag[0] - min_tag[0]);
					double local_y = (real_tag[1] - min_tag[1]); // / (max_tag[1] - min_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						//std::cout << "real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "6. no strain in map -> real_tag: " << real_tag << " min_tag: " << min_tag << std::endl;
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
						array_1d<array_1d<double, 4>, 4> p;

						size_t ii = 0;
						for (size_t kk = 0; kk < 4; kk++)
						{
							size_t jj = 0;
							for (size_t ll = 0; ll < 4; ll++)
							{
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk - 1;
								tag_i[1] = min_tag[1] + ll - 1;
								//std::cout << "real_tag: " << real_tag << std::endl;
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

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "7. no strain in map -> real_tag: " << real_tag << " tag_i: " << tag_i << std::endl;
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
						//size_t IType = 0; // Linear
						size_t IType = 1; // Cubic
						Matrix F = this->SelectInterpolationType(IType);
						radius = this->biDimInterpolate(p, local_x, local_y, F);
						/*int n = strain_size - 1;
						double coordinates[2] = { local_x, local_y };
						std::cout << "n: " << n << " - coordinates: " << coordinates[0] << ", " << coordinates[1] << std::endl;
						radius = this->nDimInterpolate(n, (double*)p, coordinates, F);*/
					}
					//std::cout << "norm_real_eps VS radius: " << norm_real_eps << " VS " << radius << std::endl;
					if (norm_real_eps < radius)
					{
						Is_Elastic = true;
						//std::cout << "Is_Elastic = true\n";
						stress_vector = ZeroVector(strain_size);
					}
					else
					{
						Is_Elastic = false;
						//std::cout << "Is_Elastic = false\n";
						noalias(stress_vector) = ZeroVector(strain_size);
					}
				}
				else if (RadialBasis)
				{
					std::cout << "Radial basis function - weight = norm(strain)" << std::endl;
					//1. Found the minimum value of the tag
					vector<int> min_tag(2);
					min_tag[0] = floor(real_tag[0]);
					min_tag[1] = floor(real_tag[1]);
					vector<int> max_tag(2);
					max_tag[0] = ceil(real_tag[0]);
					max_tag[1] = ceil(real_tag[1]);

					// local_x and y defines where to estimate the value on the interpolated line,
					// it is 0 at the first point and 1 and the second point
					double local_x = (real_tag[0] - min_tag[0]); // / (max_tag[0] - min_tag[0]);
					double local_y = (real_tag[1] - min_tag[1]); // / (max_tag[1] - min_tag[1]);

					if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
					{
						std::cout << "real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1]" << std::endl;
						//2. Read the corresponding radius and strain associated to the tag
						TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(real_tag);
						if (it_strain == mStrainMap.end())
						{
							std::cout << "8. no strain in map, real_tag: " << real_tag << std::endl;
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

						for (size_t kk = 0; kk < 4; kk++)
						{
							for (size_t ll = 0; ll < 4; ll++)
							{
								vector<int> tag_i(2);
								tag_i[0] = min_tag[0] + kk - 1;
								tag_i[1] = min_tag[1] + ll - 1;
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

								//2. Read the corresponding radius and strain associated to the tag
								TagStrainVectorMap::const_iterator it_strain = mStrainMap.find(tag_i);
								if (it_strain == mStrainMap.end())
								{
									std::cout << "9. no strain in map, tag_i: " << tag_i << std::endl;
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
				double phi = Globals::Pi / m;

				//1. Calculate the normalized Strain
				double norm_real_eps = norm_2(InStrain);
				Vector normalized_eps = norm_real_eps > 0 ? InStrain / norm_real_eps : InStrain;//This is the value to compare with the RADIUS

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
							std::cout << "10. no strain in map, real_tag: " << real_tag << std::endl;
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
						//double p[4][4][4][4][4];
						array_1d< array_1d< array_1d< array_1d< array_1d<double, 4>, 4>, 4>, 4>, 4> p;
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
												std::cout << "11. no strain in map, tag_i: " << tag_i << std::endl;
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
						//size_t IType = 0; // Linear
						size_t IType = 1; // Cubic
						Matrix F = this->SelectInterpolationType(IType);
						//size_t n = strain_size - 1;
						//double coordinates[5] = { local_x, local_y, local_z, local_l, local_m };
						//radius = this->nDimInterpolate(n, (double*)p, coordinates, F);
						radius = this->pentaDimInterpolate(p, local_x, local_y, local_z, local_l, local_m, F);
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
	}

	void RvePredictorCalculator::PredictStress2D(const Vector& trial_macro_scale_strain_vector, Vector& stress_vector, double& lch_macro, Matrix& const_tensor, double& EquivalentDamage, double& EquivalentDamageConverged)
	{
		//Initialize data for the interpolation
		double strain_size = trial_macro_scale_strain_vector.size();
		// std::cout << "trial_macro_scale_strain_vector: " << trial_macro_scale_strain_vector << std::endl;

		Vector Theta(strain_size - 1);
		Vector eps(strain_size);
		EquivalentDamage = 0.0;
		stress_vector.clear();

		//Calculate phi & number of division
		double m = (sqrtl(mStrainMap.size()) - 1) / 2; // m = max of tag
													   //std::cout << "PredictStress2D -> m: " << m << std::endl;
		double phi = Globals::Pi / m;
		//std::cout << "PredictStress2D -> Globals::Pi: " << Globals::Pi << std::endl;
		//std::cout << "PredictStress2D -> phi: " << phi << std::endl;

		//1. Calculate the normalized Strain
		// in JSON: norm_radius = sqrt(e_xx[0]**2 + e_yy[1]**2 + (2.0*e_xy[2])**2)
		double norm_real_eps = norm_2(trial_macro_scale_strain_vector);//This IS the value to compare with the RADIUS sored in JSON File
																	   //double norm_real_eps = sqrt(trial_macro_scale_strain_vector[0] * trial_macro_scale_strain_vector[0] + trial_macro_scale_strain_vector[1] * trial_macro_scale_strain_vector[1] + 2.0*((trial_macro_scale_strain_vector[2] / 2.0) * (trial_macro_scale_strain_vector[2] / 2.0)));
		Vector strain_direction = norm_real_eps > 0 ? trial_macro_scale_strain_vector / norm_real_eps : trial_macro_scale_strain_vector;//This is the direction
																																		//std::cout << "PredictStress2D -> strain_direction: " << strain_direction << std::endl;

		Theta = this->CalculateTheta(strain_direction);
		// std::cout << "PredictStress2D -> Theta: " << Theta << std::endl;
		Vector real_tag = Theta / phi;
		//std::cout << "PredictStress2D -> real_tag: " << real_tag << std::endl;
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

		vector<double> xi(4);		vector<double> eta(4);
		xi[0] = -1.0;				eta[0] = -1.0;
		xi[1] = -1.0;				eta[1] = 1.0;
		xi[2] = 1.0;				eta[2] = -1.0;
		xi[3] = 1.0;				eta[3] = 1.0;

		//START INTERPOLATION
		//1. Found the minimum value of the tag
		vector<int> min_tag(2);
		min_tag[0] = floor(real_tag[0]);
		min_tag[1] = floor(real_tag[1]);
		vector<int> max_tag(2);
		max_tag[0] = ceil(real_tag[0]);
		max_tag[1] = ceil(real_tag[1]);

		//Find n_damage
		vector<int> zero_tag(2);
		zero_tag[0] = 0;
		zero_tag[1] = 0;
		TagDamageMap::const_iterator it_n_dam = mTagDamageMap.find(zero_tag);
		if (it_n_dam == mTagDamageMap.end())
		{
			std::cout << "12. no strain in map, zero_tag: " << zero_tag << std::endl;
			exit(-1);
		}
		const Matrix & i_n_dam = it_n_dam->second;
		const size_t n_damage = i_n_dam.size2();

		Vector InterpRadiusVector(n_damage, 0.0);
		Vector InterpStressVectorXX(n_damage, 0.0);
		Vector InterpStressVectorYY(n_damage, 0.0);
		Vector InterpStressVectorXY(n_damage, 0.0);

		Vector InterpRadiusVectorC(n_damage, 0.0);
		Vector InterpStressVectorXXC(n_damage, 0.0);
		Vector InterpStressVectorYYC(n_damage, 0.0);
		Vector InterpStressVectorXYC(n_damage, 0.0);

		Vector damages(n_damage, 0.0);
		Vector i_radiusVector_aux(n_damage, 0.0);

		Vector i_stressVectorXX_aux(n_damage, 0.0);
		Vector i_stressVectorYY_aux(n_damage, 0.0);
		Vector i_stressVectorXY_aux(n_damage, 0.0);

		Vector damagesVector_aux(n_damage, 0.0);

		bool LinearInterp = true;
		bool CubicInterp = false; //Better then Linear Interpolation

		if (real_tag[0] == min_tag[0] && real_tag[1] == min_tag[1])
		{
			TagDamageMap::const_iterator it_n_damage = mTagDamageMap.find(min_tag);
			if (it_n_damage == mTagDamageMap.end())
			{
				std::cout << "13. no strain in map, min_tag: " << min_tag << std::endl;
				exit(-1);
			}
			const Matrix & i_n_damage = it_n_damage->second;
			for (size_t j_dls = 0; j_dls < n_damage; j_dls++)
			{
				damagesVector_aux[j_dls] = i_n_damage(0, j_dls);
				i_radiusVector_aux[j_dls] = i_n_damage(1, j_dls);
				//sigma
				i_stressVectorXX_aux[j_dls] = i_n_damage(2, j_dls);
				i_stressVectorYY_aux[j_dls] = i_n_damage(3, j_dls);
				i_stressVectorXY_aux[j_dls] = i_n_damage(4, j_dls);
			}

			size_t idx_final = 0;
			for (auto i_sort : sort_indexes(damagesVector_aux)) { // only c++11
				damages[idx_final] = damagesVector_aux[i_sort];
				InterpRadiusVector[idx_final] = i_radiusVector_aux[i_sort];
				//sigma
				InterpStressVectorXX[idx_final] = i_stressVectorXX_aux[i_sort];
				InterpStressVectorYY[idx_final] = i_stressVectorYY_aux[i_sort];
				InterpStressVectorXY[idx_final] = i_stressVectorXY_aux[i_sort];
				idx_final += 1;
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
				if (LinearInterp)
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

							TagDamageMap::const_iterator it_dls = mTagDamageMap.find(tag_i);
							if (it_dls == mTagDamageMap.end())
							{
								std::cout << "14. no strain in map for tag: " << tag_i << std::endl;
								exit(-1);
							}
							const Matrix & i_dls = it_dls->second;

							damages[i_dam] = i_dls(0, i_dam); //Used for check the order

							double i_radius = i_dls(1, i_dam);
							double i_stressXX = i_dls(2, i_dam);
							double i_stressYY = i_dls(3, i_dam);
							double i_stressXY = i_dls(4, i_dam);

							N_i = 1.0 / 4.0 * (1.0 + xi[i] * (local_x * 2.0 - 1.0))*(1.0 + eta[i] * (local_y * 2.0 - 1.0));
							// std::cout << "PredictStress2D -> N_i: " << N_i << std::endl;
							//N_i = ShapeFunctionValue4N(i, coordinates);
							N_sum += N_i;
							i++;
							// std::cout << "PredictStress2D -> i_radius: " << i_radius << std::endl;
							// std::cout << "PredictStress2D -> i_stressXX: " << i_stressXX << std::endl;
							// std::cout << "PredictStress2D -> i_stressYY: " << i_stressYY << std::endl;
							// std::cout << "PredictStress2D -> i_stressXY: " << i_stressXY << std::endl;

							InterpRadiusVector[i_dam] += N_i*i_radius;
							// std::cout << "PredictStress2D -> InterpRadiusVector[i_dam]: " << InterpRadiusVector[i_dam] << std::endl;

							InterpStressVectorXX[i_dam] += N_i*i_stressXX;
							InterpStressVectorYY[i_dam] += N_i*i_stressYY;
							InterpStressVectorXY[i_dam] += N_i*i_stressXY;
						}
					}
					if (abs(1 - N_sum) > 1.0e-15)
					{
						std::cout << "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION ON PredictElasticity - N_sum = " << N_sum << " ; shape_index = " << i << std::endl;
						KRATOS_THROW_ERROR(std::logic_error, "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
					}
					// Vector i_eps(3, 0.0);
					// this->ReconstructStrain(i_eps, InterpRadiusVector[i_dam], Theta);
					// std::cout << "PredictStress2D -> ReconstructStrain: " << i_eps << std::endl;

					// double norm_i_eps = norm_2(i_eps);//This IS the value to compare with the RADIUS sored in JSON File
					// Vector i_eps_direction = norm_i_eps > 0 ? i_eps / norm_i_eps : i_eps;//This is the direction
					// std::cout << "PredictStress2D -> strain_direction: " << strain_direction << std::endl;
					// Vector i_eps_Theta(strain_size - 1);
					// i_eps_Theta = this->CalculateTheta(i_eps_direction);
					// std::cout << "PredictStress2D -> Theta: " << Theta << std::endl;
					// Vector real_i_eps_tag = i_eps_Theta / phi;
					// std::cout << "PredictStress2D -> real_i_eps_tag: " << real_i_eps_tag << std::endl;


					// Vector i_stress_el = prod(mC0, i_eps);
					// std::cout << "PredictStress2D -> mC0: " << mC0 << std::endl;
					// std::cout << "PredictStress2D -> i_stress_el: " << i_stress_el << std::endl;
					// Vector final_stress = (1 - damages[i_dam])*i_stress_el;
					// std::cout << "PredictStress2D -> damages[i_dam]: " << damages[i_dam] << std::endl;
					// std::cout << InterpRadiusVector[i_dam] << i_eps << "   " << final_stress << std::endl;
					// std::cout << "INTERPOLATED: " << std::endl;
					// std::cout << InterpRadiusVector[i_dam] << i_eps << "   (" << InterpStressVectorXX[i_dam] << "," << InterpStressVectorYY[i_dam] << "," << InterpStressVectorXY[i_dam] << ")" << std::endl;

					// InterpStressVectorXX[i_dam] = final_stress[0];
					// InterpStressVectorYY[i_dam] = final_stress[1];
					// InterpStressVectorXY[i_dam] = final_stress[2];
					// std::cout << "(1-D)*SIGMA_EL: " << std::endl;
					// std::cout << InterpRadiusVector[i_dam] << i_eps << "   (" << final_stress[0] << "," << final_stress[1] << "," << final_stress[2] << ")" << std::endl;

				}
				else if (CubicInterp)
				{
					//std::cout << "CubicInterp" << std::endl;
					array_1d<array_1d<double, 4>, 4> p;
					array_1d<array_1d<double, 4>, 4> q;
					array_1d<array_1d<double, 4>, 4> r;
					array_1d<array_1d<double, 4>, 4> s;

					size_t ii = 0;
					for (size_t kk = 0; kk < 4; kk++)
					{
						size_t jj = 0;
						for (size_t ll = 0; ll < 4; ll++)
						{
							vector<int> tag_i(2);
							tag_i[0] = min_tag[0] + kk - 1;
							tag_i[1] = min_tag[1] + ll - 1;
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

							TagDamageMap::const_iterator it_dls = mTagDamageMap.find(tag_i);
							if (it_dls == mTagDamageMap.end())
							{
								std::cout << "15. no strain in map, tag_i: " << tag_i << std::endl;
								exit(-1);
							}
							const Matrix & i_dls = it_dls->second;

							damages[i_dam] = i_dls(0, i_dam); //Used for check the order

							double i_radius = i_dls(1, i_dam);
							double i_stressXX = i_dls(2, i_dam);
							double i_stressYY = i_dls(3, i_dam);
							double i_stressXY = i_dls(4, i_dam);

							p[ii][jj] = i_radius;
							q[ii][jj] = i_stressXX;
							r[ii][jj] = i_stressYY;
							s[ii][jj] = i_stressXY;
							jj++;
						}
						ii++;
					}
					//size_t IType = 0; // Linear
					size_t IType = 1; // Cubic
					Matrix F = this->SelectInterpolationType(IType);
					InterpRadiusVector[i_dam] = this->biDimInterpolate(p, local_x, local_y, F);

					Vector i_eps(3, 0.0);
					this->ReconstructStrain(i_eps, InterpRadiusVector[i_dam], Theta);
					Vector i_stress_el = prod(mC0, i_eps);
					Vector final_stress = (1 - damages[i_dam])*i_stress_el;

					InterpStressVectorXX[i_dam] = final_stress[0];
					InterpStressVectorYY[i_dam] = final_stress[1];
					InterpStressVectorXY[i_dam] = final_stress[2];
				}
			}
			// std::cout << "**************************************************************" << std::endl;
		}

		//Regularization of Fracture Energy with lf
		//std::cout << "damages: " << damages << std::endl;

		this->CalculateFractureEnergy(InterpRadiusVector, Theta, InterpStressVectorXX, InterpStressVectorYY, InterpStressVectorXY);
		this->CalculateAlpha(lch_macro);

		Vector InterpRadiusVectorOld(n_damage, 0.0);
		noalias(InterpRadiusVectorOld) = InterpRadiusVector;
		//std::cout << "InterpRadiusVector BEFORE REG: " << InterpRadiusVector << std::endl;
		this->RegularizeInterpRadius(InterpRadiusVector, Theta);
		//std::cout << "InterpRadiusVector AFTER REG: " << InterpRadiusVector << std::endl;

		LinearInterp = true; //Better then Cubic Interpolation
		CubicInterp = false;

		if (norm_real_eps < InterpRadiusVector[0])
		{
			//std::cout << "Is_Elastic = true\n";
			EquivalentDamage = 0.0;
			noalias(const_tensor) = (1 - EquivalentDamage)*mC0;
			noalias(stress_vector) = prod(const_tensor, trial_macro_scale_strain_vector);
			return;
		}
		else if (norm_real_eps > InterpRadiusVector[n_damage - 1])
		{
			//std::cout << "Is_Elastic = false -> OUT OF BOUND\n";
			double x = norm_real_eps;
			double x_1 = InterpRadiusVector[n_damage - 2];
			double x0 = InterpRadiusVector[n_damage - 1];
			double diff = (x0 - x_1);
			if (diff < 1.0e-12)
				diff = 1.0;

			double y_1 = damages[n_damage - 2];
			double y0 = damages[n_damage - 1];

			EquivalentDamage = y0 + (y0 - y_1)*(x - x0) / diff;

			//sigma
			stress_vector[0] = InterpStressVectorXX[n_damage - 1] + (InterpStressVectorXX[n_damage - 1] - InterpStressVectorXX[n_damage - 2])*(x - x0) / diff;
			stress_vector[1] = InterpStressVectorYY[n_damage - 1] + (InterpStressVectorYY[n_damage - 1] - InterpStressVectorYY[n_damage - 2])*(x - x0) / diff;
			stress_vector[2] = InterpStressVectorXY[n_damage - 1] + (InterpStressVectorXY[n_damage - 1] - InterpStressVectorXY[n_damage - 2])*(x - x0) / diff;

			// //std::cout << "Is_Elastic = false -> OUT OF BOUND\n";
			// double x = norm_real_eps;
			// Vector x_eps(3, 0.0);
			// this->ReconstructStrain(x_eps, x, Theta);
			// double x0 = InterpRadiusVector[n_damage - 1];
			// Vector x0_eps(3, 0.0);
			// this->ReconstructStrain(x0_eps, x0, Theta);
			// double x_1 = InterpRadiusVector[n_damage - 2];
			// Vector x_1_eps(3, 0.0);
			// this->ReconstructStrain(x_1_eps, x_1, Theta);
			//
			// //sigma
			// stress_vector[0] = (x0_eps[0] - x_1_eps[0]) ? InterpStressVectorXX[n_damage - 1] + (InterpStressVectorXX[n_damage - 1] - InterpStressVectorXX[n_damage - 2])*(x_eps[0] - x0_eps[0]) / (x0_eps[0] - x_1_eps[0]) : InterpStressVectorXX[n_damage - 1] + (InterpStressVectorXX[n_damage - 1] - InterpStressVectorXX[n_damage - 2])*(x_eps[0] - x0_eps[0]);
			// stress_vector[1] = (x0_eps[1] - x_1_eps[1]) ? InterpStressVectorYY[n_damage - 1] + (InterpStressVectorYY[n_damage - 1] - InterpStressVectorYY[n_damage - 2])*(x_eps[1] - x0_eps[1]) / (x0_eps[1] - x_1_eps[1]) : InterpStressVectorYY[n_damage - 1] + (InterpStressVectorYY[n_damage - 1] - InterpStressVectorYY[n_damage - 2])*(x_eps[1] - x0_eps[1]);
			// stress_vector[2] = (x0_eps[2] - x_1_eps[2]) ? InterpStressVectorXY[n_damage - 1] + (InterpStressVectorXY[n_damage - 1] - InterpStressVectorXY[n_damage - 2])*(x_eps[2] - x0_eps[2]) / (x0_eps[2] - x_1_eps[2]) : InterpStressVectorXY[n_damage - 1] + (InterpStressVectorXY[n_damage - 1] - InterpStressVectorXY[n_damage - 2])*(x_eps[2] - x0_eps[2]);
		}
		else
		{
			if (LinearInterp)
			{
				for (size_t jj = 1; jj < n_damage; jj++)
				{
					if (norm_real_eps <= InterpRadiusVector[jj])
					{
						//std::cout << "Is_Elastic = false -> INSIDE OF BOUND\n";
						double x = norm_real_eps;
						double x0 = InterpRadiusVector[jj - 1];
						double x1 = InterpRadiusVector[jj];

						double diff = (x1 - x0);
						if (diff < 1.0e-12)
							diff = 1.0;

						double y0 = damages[jj - 1];
						double y1 = damages[jj];

						EquivalentDamage = y0 + (y1 - y0)*(x - x0) / diff;

						stress_vector[0] = InterpStressVectorXX[jj - 1] + (InterpStressVectorXX[jj] - InterpStressVectorXX[jj - 1])*(x - x0) / diff;
						stress_vector[1] = InterpStressVectorYY[jj - 1] + (InterpStressVectorYY[jj] - InterpStressVectorYY[jj - 1])*(x - x0) / diff;
						stress_vector[2] = InterpStressVectorXY[jj - 1] + (InterpStressVectorXY[jj] - InterpStressVectorXY[jj - 1])*(x - x0) / diff;

						break;

						// //std::cout << "Is_Elastic = false -> INSIDE OF BOUND\n";
						// double x = norm_real_eps;
						// Vector x_eps(3, 0.0);
						// //std::cout << "norm_real_eps: " << norm_real_eps << std::endl;
						// this->ReconstructStrain(x_eps, x, Theta);
						// //std::cout << "x_eps: " << x_eps << std::endl;
						// double x0 = InterpRadiusVector[jj - 1];
						// Vector x0_eps(3, 0.0);
						// this->ReconstructStrain(x0_eps, x0, Theta);
						// //std::cout << "damages[jj - 1]: " << damages[jj - 1] << std::endl;
						// //std::cout << "x0_eps: " << x0_eps << std::endl;
						// double x1 = InterpRadiusVector[jj];
						// Vector x1_eps(3, 0.0);
						// this->ReconstructStrain(x1_eps, x1, Theta);
						// //std::cout << "damages[jj]: " << damages[jj] << std::endl;
						// //std::cout << "x1_eps: " << x1_eps << std::endl;
						//
						// stress_vector[0] = (x1_eps[0] - x0_eps[0]) > 1.0e-12 ? InterpStressVectorXX[jj - 1] + (InterpStressVectorXX[jj] - InterpStressVectorXX[jj - 1])*(trial_macro_scale_strain_vector[0] - x0_eps[0]) / (x1_eps[0] - x0_eps[0]) : InterpStressVectorXX[jj - 1] + (InterpStressVectorXX[jj] - InterpStressVectorXX[jj - 1])*(x_eps[0] - x0_eps[0]);
						// stress_vector[1] = (x1_eps[1] - x0_eps[1]) > 1.0e-12 ? InterpStressVectorYY[jj - 1] + (InterpStressVectorYY[jj] - InterpStressVectorYY[jj - 1])*(trial_macro_scale_strain_vector[1] - x0_eps[1]) / (x1_eps[1] - x0_eps[1]) : InterpStressVectorYY[jj - 1] + (InterpStressVectorYY[jj] - InterpStressVectorYY[jj - 1])*(x_eps[1] - x0_eps[1]);
						// stress_vector[2] = (x1_eps[2] - x0_eps[2]) > 1.0e-12 ? InterpStressVectorXY[jj - 1] + (InterpStressVectorXY[jj] - InterpStressVectorXY[jj - 1])*(trial_macro_scale_strain_vector[2] - x0_eps[2]) / (x1_eps[2] - x0_eps[2]) : InterpStressVectorXY[jj - 1] + (InterpStressVectorXY[jj] - InterpStressVectorXY[jj - 1])*(x_eps[2] - x0_eps[2]);

						break;
					}
					//else if (norm_real_eps == InterpRadiusVector[jj])
					//{
					//	EquivalentDamage = damages[jj];
					//	stress_vector[0] = InterpStressVectorXX[jj];
					//	stress_vector[1] = InterpStressVectorYY[jj];
					//	stress_vector[2] = InterpStressVectorXY[jj];
					//
					//	break;
					//}
				}
			}
			else if (CubicInterp)
			{
				array_1d<double, 4> t;
				array_1d<double, 4> s_xx;
				array_1d<double, 4> s_yy;
				array_1d<double, 4> s_xy;
				double local_coord;

				for (size_t jj = 1; jj < n_damage; jj++)
				{
					if (norm_real_eps < InterpRadiusVector[jj])
					{
						//std::cout << "Is_Elastic = false -> INSIDE OF BOUND\n";
						double x = norm_real_eps;

						double x0 = InterpRadiusVector[jj - 1];
						double x1 = InterpRadiusVector[jj];

						local_coord = (x - x0) / (x1 - x0);

						double y0 = damages[jj - 1];
						double y1 = damages[jj];

						for (size_t ind = 0; ind < 4; ind++)
						{
							size_t a_ind = jj - 1 + ind - 1;
							if (a_ind > n_damage - 1)
							{
								a_ind = n_damage - 1;
							}
							else if (a_ind < 0)
							{
								a_ind = 0;
							}
							t[ind] = InterpRadiusVector[a_ind];
							s_xx[ind] = InterpStressVectorXX[a_ind];
							s_yy[ind] = InterpStressVectorYY[a_ind];
							s_xy[ind] = InterpStressVectorXY[a_ind];
						}
						break;
					}
				}
				//size_t IType = 0; // Linear
				size_t IType = 1; // Cubic
				Matrix F = this->SelectInterpolationType(IType);
				EquivalentDamage = this->Interpolate(t, local_coord, F);
				stress_vector[0] = this->Interpolate(s_xx, local_coord, F);
				stress_vector[1] = this->Interpolate(s_yy, local_coord, F);
				stress_vector[2] = this->Interpolate(s_xy, local_coord, F);

				//if (SIGMA_SIGN(stress_vector[0]) != SIGMA_SIGN(InterpStressVectorXX[0]))
				//{
				//	//std::cout << "XX InterpRadiusVector(" << InterpRadiusVector << ") != InterpRadiusVectorOld(" << InterpRadiusVectorOld << ")" << std::endl;
				//	//std::cout << "InterpStressVectorXX(" << InterpStressVectorXX << ")" << std::endl;
				//	//std::cout << "XX SIGMA_SIGN(" << stress_vector[0] << ") != SIGMA_SIGN(" << InterpStressVectorXX[0] << ")" << std::endl;
				//	//std::cout << "SIGMA_SIGN(stress_vector[0]): " << SIGMA_SIGN(stress_vector[0]) << " != SIGMA_SIGN(InterpStressVectorXX[0]): " << SIGMA_SIGN(InterpStressVectorXX[n_damage - 1]) << std::endl;
				//	stress_vector[0] = 0.0; // InterpStressVectorXX[0] * 1.0e-4;
				//}
				//if (SIGMA_SIGN(stress_vector[1]) != SIGMA_SIGN(InterpStressVectorYY[0]))
				//{
				//	//std::cout << "YY InterpRadiusVector(" << InterpRadiusVector << ") != InterpRadiusVectorOld(" << InterpRadiusVectorOld << ")" << std::endl;
				//	//std::cout << "InterpStressVectorYY(" << InterpStressVectorYY << ")" << std::endl;
				//	//std::cout << "YY SIGMA_SIGN(" << stress_vector[1] << ") != SIGMA_SIGN(" << InterpStressVectorYY[0] << ")" << std::endl;
				//	//std::cout << "SIGMA_SIGN(stress_vector[1]): " << SIGMA_SIGN(stress_vector[0]) << " != SIGMA_SIGN(InterpStressVectorYY[n_damage - 1]): " << SIGMA_SIGN(InterpStressVectorYY[n_damage - 1]) << std::endl;
				//	stress_vector[1] = 0.0; // InterpStressVectorYY[0] * 1.0e-4;
				//}
				//if (SIGMA_SIGN(stress_vector[2]) != SIGMA_SIGN(InterpStressVectorXY[0]))
				//{
				//	//std::cout << "XY InterpRadiusVector(" << InterpRadiusVector << ") != InterpRadiusVectorOld(" << InterpRadiusVectorOld << ")" << std::endl;
				//	//std::cout << "InterpStressVectorXY(" << InterpStressVectorXY << ")" << std::endl;
				//	//std::cout << "XY SIGMA_SIGN(" << stress_vector[2] << ") != SIGMA_SIGN(" << InterpStressVectorXY[0] << ")" << std::endl;
				//	//std::cout << "SIGMA_SIGN(stress_vector[2]): " << SIGMA_SIGN(stress_vector[2]) << " != SIGMA_SIGN(InterpStressVectorXY[n_damage - 1]): " << SIGMA_SIGN(InterpStressVectorXY[n_damage - 1]) << std::endl;
				//	stress_vector[2] = 0.0; // InterpStressVectorXY[0] * 1.0e-4;
				//}
			}
		}

		//std::cout << "stress_vector: " << stress_vector << std::endl;
		//Calculate damage
		EquivalentDamage = 0.0;
		Vector elastic_stress(strain_size, 0.0);
		noalias(elastic_stress) = prod(mC0, trial_macro_scale_strain_vector);
		double norm_stress_difference = norm_2(elastic_stress - stress_vector);
		EquivalentDamage = norm_stress_difference / norm_2(elastic_stress);
		//double d_sigma_xx = elastic_stress[0] - stress_vector[0];
		//double d_sigma_yy = elastic_stress[1] - stress_vector[1];
		//double d_sigma_xy = elastic_stress[2] - stress_vector[2];
		//double norm_stress_difference = sqrt(d_sigma_xx*d_sigma_xx + d_sigma_yy*d_sigma_yy + 2.0 * d_sigma_xy*d_sigma_xy);
		//EquivalentDamage = norm_stress_difference / sqrt(elastic_stress[0] * elastic_stress[0] + elastic_stress[1] * elastic_stress[1] + 2.0 * elastic_stress[1] * elastic_stress[1]);

		EquivalentDamage = std::max(std::min(EquivalentDamage, 1.0), 0.0);
		if (EquivalentDamage < EquivalentDamageConverged)
			EquivalentDamage = EquivalentDamageConverged;

		//Calculate damaged Constitutive Tensor
		const_tensor.clear();
		noalias(const_tensor) = (1 - EquivalentDamage)*mC0;
		//noalias(const_tensor) = mC0;

		//Calculate stress
		//stress_vector.clear();
		//noalias(stress_vector) = prod(const_tensor, trial_macro_scale_strain_vector);

		//noalias(const_tensor) = (1 - EquivalentDamageConverged)*mC0;
	}

	void RvePredictorCalculator::ReconstructStrain(Vector& eps, const double& radius, const Vector& Theta)
	{
		eps[0] = radius*cos(Theta[0]);
		eps[1] = radius*sin(Theta[0])*cos(Theta[1]);
		eps[2] = radius*sin(Theta[0])*sin(Theta[1]);
	}

	void RvePredictorCalculator::CalculateFractureEnergy(const Vector& radius, const Vector& Theta, Vector& sigma_xx, Vector& sigma_yy, Vector& sigma_xy)
	{
		const size_t n_damage = radius.size();

		Vector norm_s(n_damage, 0.0);
		for (size_t i = 0; i < n_damage; i++)
		{
			norm_s[i] = sqrt(sigma_xx[i] * sigma_xx[i] + sigma_yy[i] * sigma_yy[i] + 2.0*sigma_xy[i] * sigma_xy[i]);
		}

		auto biggest_ft = std::max_element(std::begin(norm_s), std::end(norm_s));

		m_position_biggest_ft = std::distance(std::begin(norm_s), biggest_ft);

		//Calculate Area -> Fracture Energy
		Vector e0(3, 0.0);
		this->ReconstructStrain(e0, radius[0], Theta);
		double norm_e0 = sqrt(e0[0] * e0[0] + e0[1] * e0[1] + 2.0 * ((0.5*e0[2]) * (0.5*e0[2])));
		mG0 = 0.5 * norm_e0 * norm_s[0];
		mGf_bar = mG0;
		Vector ei1(3, 0.0);
		Vector ei2(3, 0.0);
		double norm_ei1 = 0.0;
		double norm_ei2 = 0.0;
		for (size_t i = m_position_biggest_ft - 1; i < n_damage - 1; i++)
			//for (size_t i = 0; i < n_damage - 1; i++)
		{
			this->ReconstructStrain(ei1, radius[i], Theta);
			this->ReconstructStrain(ei2, radius[i + 1], Theta);
			norm_ei1 = sqrt(ei1[0] * ei1[0] + ei1[1] * ei1[1] + 2.0 * ((0.5*ei1[2]) * (0.5*ei1[2])));
			norm_ei2 = sqrt(ei2[0] * ei2[0] + ei2[1] * ei2[1] + 2.0 * ((0.5*ei2[2]) * (0.5*ei2[2])));
			mGf_bar += std::abs(0.5 * (norm_ei2 - norm_ei1)*(norm_s[i] + norm_s[i + 1]));
		}
	}

	void RvePredictorCalculator::CalculateAlpha(double& lch_macro)
	{
		double Area = 0.5; //Reference Macro Triangle 2D3N
		double lch_ref = 0.000;
		lch_ref = std::sqrt(2.0*std::abs(Area));
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
	}

	void RvePredictorCalculator::RegularizeInterpRadius(Vector& InterpRadius, Vector& Theta)
	{
		size_t n_radius = InterpRadius.size();
		for (size_t i = m_position_biggest_ft - 1; i < n_radius; i++)
			//for (size_t i = 0; i < n_radius; i++)
		{
			// eps[i] = eps[i] + mAlpha * (eps[i] - eps[0]);
			// ej = ej+(ej-ep)*S;

			Vector e(3, false);
			e(0) = InterpRadius[i] * cos(Theta[0]) + mAlpha * (InterpRadius[i] * cos(Theta[0]) - InterpRadius[m_position_biggest_ft - 1] * cos(Theta[0]));
			e(1) = InterpRadius[i] * sin(Theta[0])*cos(Theta[1]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*cos(Theta[1]) - InterpRadius[m_position_biggest_ft - 1] * sin(Theta[0])*cos(Theta[1]));
			e(2) = InterpRadius[i] * sin(Theta[0])*sin(Theta[1]) + mAlpha * (InterpRadius[i] * sin(Theta[0])*sin(Theta[1]) - InterpRadius[m_position_biggest_ft - 1] * sin(Theta[0])*sin(Theta[1]));

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

	void RvePredictorCalculator::PredictStress3D(const Vector& trial_macro_scale_strain_vector, Vector& stress_vector, Matrix& const_tensor, double& EquivalentDamage, double& EquivalentDamageConverged)
	{
		//Initialize data for the interpolation
		double strain_size = trial_macro_scale_strain_vector.size();
		Vector InStrain = trial_macro_scale_strain_vector;
		//InStrain[2] = 0.5 * InStrain[2];
		Vector Theta(strain_size - 1);
		Vector eps(strain_size);
		EquivalentDamage = 0.0;

		//Calculate phi & number of division
		double m = (sqrtl(mStrainMap.size()) - 1) / 5; // m = max of tag
													   //std::cout << "m: " << m << std::endl;
		double phi = Globals::Pi / m;
		//std::cout << "phi: " << phi << " rad" << std::endl;

		//1. Calculate the normalized Strain
		double norm_real_eps = norm_2(InStrain);//This is the value to compare with the RADIUS
												//std::cout << "norm_real_eps: " << norm_real_eps << std::endl;
												//std::cout << "InStrain: " << InStrain << std::endl;
		Vector strain_direction = norm_real_eps > 0 ? InStrain / norm_real_eps : InStrain;//This is the direction
																						  //std::cout << "strain_direction: " << strain_direction << std::endl;

		Theta = this->CalculateTheta(strain_direction);
		Vector real_tag = Theta / phi;
		//std::cout << "real_tag: " << real_tag << std::endl;

		vector<double> xi(4);		vector<double> eta(4);
		xi[0] = -1.0;				eta[0] = -1.0;
		xi[1] = -1.0;				eta[1] = 1.0;
		xi[2] = 1.0;				eta[2] = -1.0;
		xi[3] = 1.0;				eta[3] = 1.0;

		//START INTERPOLATION
		//1. Found the minimum value of the tag
		vector<int> min_tag(2);
		min_tag[0] = floor(real_tag[0]);
		min_tag[1] = floor(real_tag[1]);
		//std::cout << "min_tag: " << min_tag << std::endl;
		vector<int> max_tag(2);
		max_tag[0] = ceil(real_tag[0]);
		max_tag[1] = ceil(real_tag[1]);
		//std::cout << "max_tag: " << max_tag << std::endl;

		// local_x and y defines where to estimate the value on the interpolated line,
		// it is 0 at the first point and 1 and the second point
		double local_x = (real_tag[0] - min_tag[0]) / (max_tag[0] - min_tag[0]);
		double local_y = (real_tag[1] - min_tag[1]) / (max_tag[1] - min_tag[1]);

		Vector coordinates(2, 0.0);
		coordinates[0] = local_x;
		coordinates[1] = local_y;

		//Find n_damage
		vector<int> zero_tag(2);
		zero_tag[0] = 0;
		zero_tag[1] = 0;
		TagDamageMap::const_iterator it_n_damage = mTagDamageMap.find(zero_tag);
		if (it_n_damage == mTagDamageMap.end())
		{
			std::cout << "16. no strain in map, zero_tag: " << zero_tag << std::endl;
			exit(-1);
		}
		const Matrix & i_n_damage = it_n_damage->second;
		const size_t n_damage = i_n_damage.size2();
		Vector InterpRadiusVector(n_damage, 0.0);
		//sigma
		Vector InterpStressVectorXX(n_damage, 0.0);
		Vector InterpStressVectorYY(n_damage, 0.0);
		Vector InterpStressVectorXY(n_damage, 0.0);


		Vector damages(n_damage, 0.0);

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
				//N_i = 1.0 / 4.0 * (1.0 + xi[i] * (local_x * 2.0 - 1.0))*(1.0 + eta[i] * (local_y * 2.0 - 1.0));
				N_i = ShapeFunctionValue4N(i, coordinates);
				N_sum += N_i;
				i++;

				//std::cout << "tag_i: " << tag_i << std::endl;
				//2. Read the corresponding radius and strain associated to the tag
				TagDamageMap::const_iterator it_dls = mTagDamageMap.find(tag_i);
				if (it_dls == mTagDamageMap.end())
				{
					std::cout << "17. no strain in map, tag_i: " << tag_i << std::endl;
					exit(-1);
				}
				const Matrix & i_dls = it_dls->second;
				Vector damagesVector_aux(n_damage, 0.0);
				Vector i_radiusVector_aux(n_damage, 0.0);

				Vector i_stressVectorXX_aux(n_damage, 0.0);
				Vector i_stressVectorYY_aux(n_damage, 0.0);
				Vector i_stressVectorXY_aux(n_damage, 0.0);
				for (size_t j_dls = 0; j_dls < n_damage; j_dls++)
				{
					damagesVector_aux[j_dls] = i_dls(0, j_dls);
					i_radiusVector_aux[j_dls] = i_dls(1, j_dls);
					//sigma
					i_stressVectorXX_aux[j_dls] = i_dls(2, j_dls);
					i_stressVectorYY_aux[j_dls] = i_dls(3, j_dls);
					i_stressVectorXY_aux[j_dls] = i_dls(4, j_dls);
				}

				Vector i_radiusVector(n_damage, 0.0);
				Vector i_stressVectorXX(n_damage, 0.0);
				Vector i_stressVectorYY(n_damage, 0.0);
				Vector i_stressVectorXY(n_damage, 0.0);
				size_t idx_final = 0;
				for (auto i_sort : sort_indexes(damagesVector_aux)) { // only c++11
					damages[idx_final] = damagesVector_aux[i_sort];
					i_radiusVector[idx_final] = i_radiusVector_aux[i_sort];
					//sigma
					i_stressVectorXX[idx_final] = i_stressVectorXX_aux[i_sort];
					i_stressVectorYY[idx_final] = i_stressVectorYY_aux[i_sort];
					i_stressVectorXY[idx_final] = i_stressVectorXY_aux[i_sort];
					idx_final += 1;
				}
				for (size_t ii = 0; ii < n_damage; ii++)
				{
					InterpRadiusVector[ii] += N_i*i_radiusVector[ii];
					//sigma
					InterpStressVectorXX[ii] += N_i*i_stressVectorXX[ii];
					InterpStressVectorYY[ii] += N_i*i_stressVectorYY[ii];
					InterpStressVectorXY[ii] += N_i*i_stressVectorXY[ii];
				}
			}
		}
		if (abs(1 - N_sum) > 1.0e-15)
		{
			std::cout << "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION ON PredictElasticity - N_sum = " << N_sum << " ; shape_index = " << i << std::endl;
			KRATOS_THROW_ERROR(std::logic_error, "RvePredictorCalculator - ERROR DURING THE SHAPE FUNCTION COMPUTATION", "");
		}

		//std::cout << "norm_real_eps: " << norm_real_eps << std::endl;
		//std::cout << "InterpRadiusVector: " << InterpRadiusVector << std::endl;
		//std::cout << "damages: " << damages << std::endl;
		if (norm_real_eps < InterpRadiusVector[0])
		{
			//std::cout << "Is_Elastic = true\n";
			EquivalentDamage = 0.0;

			noalias(stress_vector) = prod(mC0, InStrain);
		}
		else if (norm_real_eps > InterpRadiusVector[n_damage - 1])
		{
			double x = norm_real_eps;
			double x1 = InterpRadiusVector[n_damage - 1];
			double y1 = damages[n_damage - 1];
			//double x2 = 1.0;
			double x2 = InterpRadiusVector[n_damage - 1] + 0.5 * (InterpRadiusVector[n_damage - 1] - InterpRadiusVector[n_damage - 2]);
			double y2 = 1.0;

			double N1 = (x - x2) / (x1 - x2);
			double N2 = (x - x1) / (x1 - x2);

			EquivalentDamage = N1*y1 - N2*y2;

			//sigma
			double xx1 = InterpStressVectorXX[n_damage - 1];
			double xx2 = 0.0;
			double yy1 = InterpStressVectorYY[n_damage - 1];
			double yy2 = 0.0;
			double xy1 = InterpStressVectorXY[n_damage - 1];
			double xy2 = 0.0;

			stress_vector[0] = N1*xx1 - N2*xx2;
			stress_vector[1] = N1*yy1 - N2*yy2;
			stress_vector[2] = N1*xy1 - N2*xy2;

		}
		else
		{
			for (size_t jj = 1; jj < n_damage; jj++)
			{
				if (norm_real_eps < InterpRadiusVector[jj])
				{
					double x = norm_real_eps;
					double x1 = InterpRadiusVector[jj - 1];
					double x2 = InterpRadiusVector[jj];

					double N1 = (x - x2) / (x1 - x2);
					double N2 = (x - x1) / (x1 - x2);

					//Linear Interpolation of damage
					EquivalentDamage = N1*damages[jj - 1] - N2*damages[jj];

					stress_vector[0] = N1*InterpStressVectorXX[jj - 1] - N2*InterpStressVectorXX[jj];
					stress_vector[1] = N1*InterpStressVectorYY[jj - 1] - N2*InterpStressVectorYY[jj];
					stress_vector[2] = N1*InterpStressVectorXY[jj - 1] - N2*InterpStressVectorXY[jj];

					break;
				}
			}
		}

		//Calculate damaged Constitutive Tensor
		if (EquivalentDamage < 1.0e-12)
			EquivalentDamage = 0.0;
		if (EquivalentDamage > 1.0)
			EquivalentDamage = 1.0;

		if (EquivalentDamage < EquivalentDamageConverged)
			EquivalentDamage = EquivalentDamageConverged;

		//std::cout << "EquivalentDamage: " << EquivalentDamage << std::endl;
		//std::cout << "EquivalentDamageConverged: " << EquivalentDamageConverged << std::endl;

		noalias(const_tensor) = (1 - EquivalentDamage)*mC0;
		//std::cout << "const_tensor: " << const_tensor << std::endl;
		//noalias(stress_vector) = prod(const_tensor, InStrain);
		//std::cout << "stress_vector: " << stress_vector << std::endl;
		noalias(const_tensor) = (1 - EquivalentDamageConverged)*mC0;

	}

	double RvePredictorCalculator::ShapeFunctionValue4N(size_t ShapeFunctionIndex,
		const Vector& rPoint) const
	{
		//if (ShapeFunctionIndex == 0)
		//	return(0.25*(1.0 - rPoint[0])*(1.0 - rPoint[1]));
		//else if (ShapeFunctionIndex == 1)
		//	return(0.25*(1.0 + rPoint[0])*(1.0 - rPoint[1]));
		//else if (ShapeFunctionIndex == 2)
		//	return(0.25*(1.0 + rPoint[0])*(1.0 + rPoint[1]));
		//else if (ShapeFunctionIndex == 3)
		//	return(0.25*(1.0 - rPoint[0])*(1.0 + rPoint[1]));
		//else
		//	KRATOS_THROW_ERROR(std::logic_error,
		//	"Wrong index of shape function!", "");
		//return 0;
		if (ShapeFunctionIndex == 0)
			return(0.25*(1.0 - rPoint[0])*(1.0 - rPoint[1]));
		else if (ShapeFunctionIndex == 1)
			return(0.25*(1.0 - rPoint[0])*(1.0 + rPoint[1]));
		else if (ShapeFunctionIndex == 2)
			return(0.25*(1.0 + rPoint[0])*(1.0 - rPoint[1]));
		else if (ShapeFunctionIndex == 3)
			return(0.25*(1.0 + rPoint[0])*(1.0 + rPoint[1]));
		else
			KRATOS_THROW_ERROR(std::logic_error,
				"Wrong index of shape function!", "");
		return 0;
	}

	double RvePredictorCalculator::ShapeFunctionValue8N(size_t ShapeFunctionIndex,
		const Vector& rPoint) const
	{
		/*switch (ShapeFunctionIndex)
		{
		case 0:
		return(0.125*(1.0 - rPoint[0])*(1.0 - rPoint[1])*(1.0 - rPoint[2]));
		case 1:
		return(0.125*(1.0 + rPoint[0])*(1.0 - rPoint[1])*(1.0 - rPoint[2]));
		case 2:
		return(0.125*(1.0 + rPoint[0])*(1.0 + rPoint[1])*(1.0 - rPoint[2]));
		case 3:
		return(0.125*(1.0 - rPoint[0])*(1.0 + rPoint[1])*(1.0 - rPoint[2]));
		case 4:
		return(0.125*(1.0 - rPoint[0])*(1.0 - rPoint[1])*(1.0 + rPoint[2]));
		case 5:
		return(0.125*(1.0 + rPoint[0])*(1.0 - rPoint[1])*(1.0 + rPoint[2]));
		case 6:
		return(0.125*(1.0 + rPoint[0])*(1.0 + rPoint[1])*(1.0 + rPoint[2]));
		case 7:
		return(0.125*(1.0 - rPoint[0])*(1.0 + rPoint[1])*(1.0 + rPoint[2]));
		default:
		KRATOS_THROW_ERROR(std::logic_error,
		"Wrong index of shape function!", "");*/
		if (ShapeFunctionIndex == 0)
			return(0.125*(1.0 - rPoint[0])*(1.0 - rPoint[1])*(1.0 - rPoint[2]));
		else if (ShapeFunctionIndex == 1)
			return(0.125*(1.0 + rPoint[0])*(1.0 - rPoint[1])*(1.0 - rPoint[2]));
		else if (ShapeFunctionIndex == 2)
			return(0.125*(1.0 + rPoint[0])*(1.0 + rPoint[1])*(1.0 - rPoint[2]));
		else if (ShapeFunctionIndex == 3)
			return(0.125*(1.0 - rPoint[0])*(1.0 + rPoint[1])*(1.0 - rPoint[2]));
		else if (ShapeFunctionIndex == 4)
			return(0.125*(1.0 - rPoint[0])*(1.0 - rPoint[1])*(1.0 + rPoint[2]));
		else if (ShapeFunctionIndex == 5)
			return(0.125*(1.0 + rPoint[0])*(1.0 - rPoint[1])*(1.0 + rPoint[2]));
		else if (ShapeFunctionIndex == 6)
			return(0.125*(1.0 + rPoint[0])*(1.0 + rPoint[1])*(1.0 + rPoint[2]));
		else if (ShapeFunctionIndex == 7)
			return(0.125*(1.0 - rPoint[0])*(1.0 + rPoint[1])*(1.0 + rPoint[2]));
		else
			KRATOS_THROW_ERROR(std::logic_error,
				"Wrong index of shape function!", "");
		return 0;
	}

	double RvePredictorCalculator::ShapeFunctionValue27N(size_t ShapeFunctionIndex,
		const Vector& rPoint) const
	{
		double fx1 = 0.5 * (rPoint[0] - 1.0) * (rPoint[0]);
		double fx2 = 0.5 * (rPoint[0] + 1.0) * (rPoint[0]);
		double fx3 = 1.0 - (rPoint[0] * rPoint[0]);
		double fy1 = 0.5 * (rPoint[1] - 1.0) * (rPoint[1]);
		double fy2 = 0.5 * (rPoint[1] + 1.0) * (rPoint[1]);
		double fy3 = 1.0 - (rPoint[1] * rPoint[1]);
		double fz1 = 0.5 * (rPoint[2] - 1.0) * (rPoint[2]);
		double fz2 = 0.5 * (rPoint[2] + 1.0) * (rPoint[2]);
		double fz3 = 1.0 - (rPoint[2] * rPoint[2]);

		switch (ShapeFunctionIndex)
		{
		case 0:
			return(fx1*fy1*fz1);
		case 1:
			return(fx2*fy1*fz1);
		case 2:
			return(fx2*fy2*fz1);
		case 3:
			return(fx1*fy2*fz1);
		case 4:
			return(fx1*fy1*fz2);
		case 5:
			return(fx2*fy1*fz2);
		case 6:
			return(fx2*fy2*fz2);
		case 7:
			return(fx1*fy2*fz2);
		case 8:
			return(fx3*fy1*fz1);
		case 9:
			return(fx2*fy3*fz1);
		case 10:
			return(fx3*fy2*fz1);
		case 11:
			return(fx1*fy3*fz1);
		case 12:
			return(fx1*fy1*fz3);
		case 13:
			return(fx2*fy1*fz3);
		case 14:
			return(fx2*fy2*fz3);
		case 15:
			return(fx1*fy2*fz3);
		case 16:
			return(fx3*fy1*fz2);
		case 17:
			return(fx2*fy3*fz2);
		case 18:
			return(fx3*fy2*fz2);
		case 19:
			return(fx1*fy3*fz2);
		case 20:
			return(fx3*fy3*fz1);
		case 21:
			return(fx3*fy1*fz3);
		case 22:
			return(fx2*fy3*fz3);
		case 23:
			return(fx3*fy2*fz3);
		case 24:
			return(fx1*fy3*fz3);
		case 25:
			return(fx3*fy3*fz2);
		case 26:
			return(fx3*fy3*fz3);

		default:
			KRATOS_THROW_ERROR(std::logic_error, "Wrong index of shape function!", "");
		}

		return 0;
	}

	Vector RvePredictorCalculator::CalculateTheta(Vector NormalizedStrainVector) const {
		//Calculate Theta (for the real strain)
		size_t theta_size = NormalizedStrainVector.size() - 1;
		//std::cout << "PredictStress2D -> NormalizedStrainVector: " << NormalizedStrainVector << std::endl;
		Vector Theta(theta_size, 0.0);
		if (theta_size == 2)
		{
			//Calculate Theta1
			double theta_1;
			double denom_theta_1 = sqrt(pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2) + pow(NormalizedStrainVector[0], 2));
			//std::cout << "PredictStress2D -> denom_theta_1: " << denom_theta_1 << std::endl;
			if (denom_theta_1 == 0.0)
			{
				theta_1 = 0;
			}
			else
			{
				theta_1 = acos(NormalizedStrainVector[0] / denom_theta_1);
			}
			if (theta_1 > Globals::Pi)
				theta_1 = theta_1 - 2.0 * Globals::Pi;
			if (theta_1 < -Globals::Pi)
				theta_1 = theta_1 + 2.0 * Globals::Pi;

			//Calculate Theta2
			double theta_2;
			double denom_theta_2 = sqrt(pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2));
			//std::cout << "PredictStress2D -> denom_theta_2: " << denom_theta_2 << std::endl;
			if (denom_theta_2 == 0.0)
			{
				theta_2 = 0;
			}
			else
			{
				// theta_2 = atan(NormalizedStrainVector[2] / NormalizedStrainVector[1]);
				if (NormalizedStrainVector[2] < 0.0)
				{
					//std::cout << "PredictStress2D -> NormalizedStrainVector[2] < 0.0" << std::endl;
					//std::cout << "PredictStress2D -> 2.0 * Globals::Pi" << 2.0 * Globals::Pi << std::endl;
					//std::cout << "PredictStress2D -> NormalizedStrainVector[1] / denom_theta_2" << NormalizedStrainVector[1] / denom_theta_2 << std::endl;
					//std::cout << "PredictStress2D -> acos(NormalizedStrainVector[1] / denom_theta_2)" << acos(NormalizedStrainVector[1] / denom_theta_2) << std::endl;
					theta_2 = 2.0 * Globals::Pi - acos(NormalizedStrainVector[1] / denom_theta_2);
				}
				else
				{
					//std::cout << "PredictStress2D -> NormalizedStrainVector[2] > 0.0" << std::endl;
					theta_2 = acos(NormalizedStrainVector[1] / denom_theta_2);
				}
			}
			if (theta_2 > Globals::Pi)
				theta_2 = theta_2 - 2.0 * Globals::Pi;
			if (theta_2 < -Globals::Pi)
				theta_2 = theta_2 + 2.0 * Globals::Pi;

			//Assemble Theta
			Theta[0] = theta_1;
			Theta[1] = theta_2;
		}
		else
		{
			//Calculate Theta1
			double theta_1 = acos(NormalizedStrainVector[0] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2) + pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2) + pow(NormalizedStrainVector[0], 2)));
			if (theta_1 > Globals::Pi)
				theta_1 = theta_1 - 2 * Globals::Pi;
			if (theta_1 < -Globals::Pi)
				theta_1 = theta_1 + 2 * Globals::Pi;

			//Calculate Theta2
			double theta_2 = acos(NormalizedStrainVector[1] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2) + pow(NormalizedStrainVector[2], 2) + pow(NormalizedStrainVector[1], 2)));
			if (theta_2 > Globals::Pi)
				theta_2 = theta_2 - 2 * Globals::Pi;
			if (theta_2 < -Globals::Pi)
				theta_2 = theta_2 + 2 * Globals::Pi;

			//Calculate Theta3
			double theta_3 = acos(NormalizedStrainVector[2] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2) + pow(NormalizedStrainVector[2], 2)));
			if (theta_3 > Globals::Pi)
				theta_3 = theta_3 - 2 * Globals::Pi;
			if (theta_3 < -Globals::Pi)
				theta_3 = theta_3 + 2 * Globals::Pi;

			//Calculate Theta4
			double theta_4 = acos(NormalizedStrainVector[2] / sqrt(pow(NormalizedStrainVector[5], 2) + pow(NormalizedStrainVector[4], 2) + pow(NormalizedStrainVector[3], 2)));
			if (theta_4 > Globals::Pi)
				theta_4 = theta_4 - 2 * Globals::Pi;
			if (theta_4 < -Globals::Pi)
				theta_4 = theta_4 + 2 * Globals::Pi;

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
					theta_5 = 2 * Globals::Pi - acos(NormalizedStrainVector[4] / denom);
				}
				else
				{
					theta_5 = acos(NormalizedStrainVector[4] / denom);
				}
			}
			if (theta_5 > Globals::Pi)
				theta_5 = theta_5 - 2 * Globals::Pi;
			if (theta_5 < -Globals::Pi)
				theta_5 = theta_5 + 2 * Globals::Pi;

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
		// https://www.engr.colostate.edu/ECE481A2/Recommended_Readings/TCB.pdf

		/*% ' Basis matrix for global Tension, Continuity  & Bias
		% ' As T goes from 0 to 1 you go from n to n+1
		% '        n-1      n        n+1       n+2
		% '   T^3  -A    4+A-B-C   -4+B+C-D     D     /
		% '   T^2 +2A  -6-2A+2B+C  6-2B-C+D    -D    /
		% '   T^1  -A     A-B         B         0   /  2
		% '   T^0   0      2          0         0  /
		% '*/
		double FT = 0.0;  //%' Tension       T=+1-->Tight             T=-1--> Round
		double FB = 0.0;  //%' Bias          B=+1-->Post Shoot        B=-1--> Pre shoot
		double FC = 0.0;  //%' Continuity    C=+1-->Inverted corners  C=-1--> Box corners
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

		if (IType == 0) // Linear
		{
			FT = 0.0;  //%' Tension       T=+1-->Tight             T=-1--> Round
			FB = 0.0;  //%' Bias          B=+1-->Post Shoot        B=-1--> Pre shoot
			FC = -1.0;  //%' Continuity    C=+1-->Inverted corners  C=-1--> Box corners
		}
		else if (IType == 1) // Cubic
		{
			FT = 1.0;  //%' Tension       T=+1-->Tight             T=-1--> Round
			FB = 0.0;  //%' Bias          B=+1-->Post Shoot        B=-1--> Pre shoot
			FC = 0.0;  //%' Continuity    C=+1-->Inverted corners  C=-1--> Box corners
		}
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

	double RvePredictorCalculator::Interpolate(array_1d<double, 4> p, double x, Matrix& F) const {
		/*std::cout << ">>> cubicInterpolate \n";
		std::cout << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << '\n';*/
		return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
		//std::cout << "F: " << F << std::endl;
		//double FAX = F(0, 0)*p[0] + F(0, 1)*p[1] + F(0, 2)*p[2] + F(0, 3)*p[3];
		//double FBX = F(1, 0)*p[0] + F(1, 1)*p[1] + F(1, 2)*p[2] + F(1, 3)*p[3];
		//double FCX = F(2, 0)*p[0] + F(2, 1)*p[1] + F(2, 2)*p[2] + F(2, 3)*p[3];
		//double FDX = F(3, 0)*p[0] + F(3, 1)*p[1] + F(3, 2)*p[2] + F(3, 3)*p[3];
		//return ((FAX*x + FBX)*x + FCX)*x + FDX;
	}

	double RvePredictorCalculator::biDimInterpolate(array_1d<array_1d<double, 4>, 4> p, double x, double y, Matrix& F) const {
		/*std::cout << ">>> bicubicInterpolate \n";
		std::cout << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << '\n';*/
		array_1d<double, 4> arr;
		arr[0] = Interpolate(p[0], y, F);
		arr[1] = Interpolate(p[1], y, F);
		arr[2] = Interpolate(p[2], y, F);
		arr[3] = Interpolate(p[3], y, F);
		return Interpolate(arr, x, F);
	}

	double RvePredictorCalculator::triDimInterpolate(array_1d<array_1d<array_1d<double, 4>, 4>, 4> p, double x, double y, double z, Matrix F) const {
		array_1d<double, 4> arr;
		arr[0] = biDimInterpolate(p[0], y, z, F);
		arr[1] = biDimInterpolate(p[1], y, z, F);
		arr[2] = biDimInterpolate(p[2], y, z, F);
		arr[3] = biDimInterpolate(p[3], y, z, F);
		return Interpolate(arr, x, F);
	}

	double RvePredictorCalculator::quadDimInterpolate(array_1d<array_1d<array_1d<array_1d<double, 4>, 4>, 4>, 4> p, double x, double y, double z, double h, Matrix F) const {
		array_1d<double, 4> arr;
		arr[0] = triDimInterpolate(p[0], y, z, h, F);
		arr[1] = triDimInterpolate(p[1], y, z, h, F);
		arr[2] = triDimInterpolate(p[2], y, z, h, F);
		arr[3] = triDimInterpolate(p[3], y, z, h, F);
		return Interpolate(arr, x, F);
	}

	double RvePredictorCalculator::pentaDimInterpolate(array_1d<array_1d<array_1d<array_1d<array_1d<double, 4>, 4>, 4>, 4>, 4> p, double x, double y, double z, double h, double k, Matrix F) const {
		array_1d<double, 4> arr;
		arr[0] = quadDimInterpolate(p[0], y, z, h, k, F);
		arr[1] = quadDimInterpolate(p[1], y, z, h, k, F);
		arr[2] = quadDimInterpolate(p[2], y, z, h, k, F);
		arr[3] = quadDimInterpolate(p[3], y, z, h, k, F);
		return Interpolate(arr, x, F);
	}

	double RvePredictorCalculator::nDimInterpolate(size_t n, double* p, double coordinates[], Matrix& F) const {
		//assert(n > 0);
		//if (n == 1) {
		//	return Interpolate(p, *coordinates, F);
		//}
		//else {
		//	array_1d<double, 4> arr;
		//	int skip = 1 << (n - 1) * 2;
		//	arr[0] = nDimInterpolate(n - 1, p, coordinates + 1, F);
		//	arr[1] = nDimInterpolate(n - 1, p + skip, coordinates + 1, F);
		//	arr[2] = nDimInterpolate(n - 1, p + 2 * skip, coordinates + 1, F);
		//	arr[3] = nDimInterpolate(n - 1, p + 3 * skip, coordinates + 1, F);
		//	return Interpolate(arr, *coordinates, F);
		//}
		return 0.0;
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
		/*for(ModelPart::ElementIterator it = modelPart.ElementsBegin(); it != modelPart.ElementsEnd(); ++it)
		{
		Element& ielem = *it;
		Element::GeometryType& igeom = ielem.GetGeometry();
		m_domain_size += igeom.DomainSize();
		}*/
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

} // namespace Kratos
