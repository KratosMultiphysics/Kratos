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

#if !defined(RVE_MATERIAL_DATABASE_3D_H_INCLUDED)
#define RVE_MATERIAL_DATABASE_3D_H_INCLUDED

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include "includes/kratos_parameters.h"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include <algorithm>

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "multiscale_application.h"

namespace Kratos
{
    class RveMaterialDatabase3D
    {

    public:

		KRATOS_CLASS_POINTER_DEFINITION(RveMaterialDatabase3D);
        typedef double RealType;
        typedef Vector VectorType;
        typedef Matrix MatrixType;
        typedef size_t IndexType;

		typedef Kratos::KratosMultiScaleApplication::UnorderedVectorVectorMap VectorVectorMap;
		typedef Kratos::KratosMultiScaleApplication::UnorderedVectorMatrixMap VectorMatrixMap;

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

    public:

		RveMaterialDatabase3D(std::string C0_file_name, std::string hash_file_name, std::string json_file_name)
			: mC0FileName(C0_file_name)
			, mHashFileName(hash_file_name) // TODO: Remove it
			, mJsonFileName(json_file_name)
			, mConstitutiveTensor()
			, mRveElasticMap()
			, mRveMaterialMap()
			, mDamageSize(0)
			, mNumberOfDivision(0)
		{
			std::string mCo_file = mC0FileName;
			std::cout << "mCo_file: " << mCo_file << std::endl;
			this->ReadAndCreateConstTensorData(mCo_file);

			std::string hash_file = mHashFileName; // TODO: Remove it
			std::cout << "hash_file_name: " << hash_file << std::endl;
			this->ReadAndCreatePredictorData(hash_file); // TODO: Remove it

			mNumberOfDivision = (sqrtl(mRveElasticMap.size()) - 1) / 5;

			std::string json_file = mJsonFileName;
			std::cout << "json_file: " << json_file << std::endl;
			this->ReadAndCreatePredictorDataDamageHistory(json_file);
        }

		virtual ~RveMaterialDatabase3D()
        {
        }

    public:

		virtual void GetMaterialMapSize(size_t& RveMaterialMapSize)
		{
			RveMaterialMapSize = mRveMaterialMap.size();
		}

		virtual void GetElasticValues(const vector<int>& Tag,
									  double& Radius,
									  Vector& ElasticStrain,
									  Vector& ElasticStress)
		{
			VectorVectorMap::const_iterator it = mRveElasticMap.find(Tag);
			if (it == mRveElasticMap.end())
			{
				std::cout << "Radius-Strain-Stress Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Vector & i = it->second;
			Radius = i[0];
			ElasticStrain[0] = i[1];
			ElasticStrain[1] = i[2];
			ElasticStrain[2] = i[3];
			ElasticStrain[3] = i[4];
			ElasticStrain[4] = i[5];
			ElasticStrain[5] = i[6];
			ElasticStress[0] = i[7];
			ElasticStress[1] = i[8];
			ElasticStress[2] = i[9];
			ElasticStress[3] = i[10];
			ElasticStress[4] = i[11];
			ElasticStress[5] = i[12];
		}

		virtual void GetElasticRadius(const vector<int>& Tag, double& Radius)
		{
			VectorVectorMap::const_iterator it = mRveElasticMap.find(Tag);
			if (it == mRveElasticMap.end())
			{
				std::cout << "Radius Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Vector & i = it->second;
			Radius = i[0];
		}

		virtual void GetElasticStrain(const vector<int>& Tag, Vector& ElasticStrain)
		{
			VectorVectorMap::const_iterator it = mRveElasticMap.find(Tag);
			if (it == mRveElasticMap.end())
			{
				std::cout << "Strain Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Vector & i = it->second;
			ElasticStrain[0] = i[1];
			ElasticStrain[1] = i[2];
			ElasticStrain[2] = i[3];
			ElasticStrain[3] = i[4];
			ElasticStrain[4] = i[5];
			ElasticStrain[5] = i[6];
		}

		virtual void GetElasticStress(const vector<int>& Tag, Vector& ElasticStress)
		{
			VectorVectorMap::const_iterator it = mRveElasticMap.find(Tag);
			if (it == mRveElasticMap.end())
			{
				std::cout << "Stress Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Vector & i = it->second;
			ElasticStress[0] = i[7];
			ElasticStress[1] = i[8];
			ElasticStress[2] = i[9];
			ElasticStress[3] = i[10];
			ElasticStress[4] = i[11];
			ElasticStress[5] = i[12];
		}

		virtual void GetMaterialValues(const vector<int>& Tag,
									   Vector& DamageVector,
									   Vector& RadiusVector,
									   Matrix& StressMatrix)
		{
			size_t dam_dim = DamageSize();
			Vector DamageVector_aux(dam_dim,0.0);
			Vector RadiusVector_aux(dam_dim,0.0);
			Matrix StressMatrix_aux(6,dam_dim,0.0);

			VectorMatrixMap::const_iterator it = mRveMaterialMap.find(Tag);
			if (it == mRveMaterialMap.end())
			{
				std::cout << "Damage Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Matrix & i = it->second;
			for (size_t j = 0; j < mDamageSize; j++)
			{
				DamageVector_aux[j] = i(0, j);
				RadiusVector_aux[j] = i(1, j);
				StressMatrix_aux(0, j) = i(2, j);
				StressMatrix_aux(1, j) = i(3, j);
				StressMatrix_aux(2, j) = i(4, j);
				StressMatrix_aux(3, j) = i(5, j);
				StressMatrix_aux(4, j) = i(6, j);
				StressMatrix_aux(5, j) = i(7, j);
			}
			size_t idx_final = 0;
			for (auto i_sort : sort_indexes(DamageVector_aux)) { // only c++11
				DamageVector[idx_final] = DamageVector_aux[i_sort];
				RadiusVector[idx_final] = RadiusVector_aux[i_sort];
				//sigma
				StressMatrix(0, idx_final) = StressMatrix_aux(0, i_sort);
				StressMatrix(1, idx_final) = StressMatrix_aux(1, i_sort);
				StressMatrix(2, idx_final) = StressMatrix_aux(2, i_sort);
				StressMatrix(3, idx_final) = StressMatrix_aux(3, i_sort);
				StressMatrix(4, idx_final) = StressMatrix_aux(4, i_sort);
				StressMatrix(5, idx_final) = StressMatrix_aux(5, i_sort);
				idx_final += 1;
			}
		}

		virtual void GetDamageVector(const vector<int>& Tag, Vector& DamageVector)
		{
			VectorMatrixMap::const_iterator it = mRveMaterialMap.find(Tag);
			if (it == mRveMaterialMap.end())
			{
				std::cout << "Damage Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Matrix & i = it->second;
			for (size_t j = 0; j < mDamageSize; j++)
			{
				DamageVector[j] = i(0, j);
			}
		}

		virtual void GetRadiusVector(const vector<int>& Tag, Vector& RadiusVector)
		{
			VectorMatrixMap::const_iterator it = mRveMaterialMap.find(Tag);
			if (it == mRveMaterialMap.end())
			{
				std::cout << "Damage Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Matrix & i = it->second;
			for (size_t j = 0; j < mDamageSize; j++)
			{
				RadiusVector[j] = i(1, j);
			}
		}

		virtual void GetStressMatrix(const vector<int>& Tag, Matrix& StressMatrix)
		{
			VectorMatrixMap::const_iterator it = mRveMaterialMap.find(Tag);
			if (it == mRveMaterialMap.end())
			{
				std::cout << "Damage Data not present in map, Tag: " << Tag << std::endl;
				exit(-1);
			}
			const Matrix & i = it->second;
			for (size_t j = 0; j < mDamageSize; j++)
			{
				StressMatrix(0, j) = i(2, j);
				StressMatrix(1, j) = i(3, j);
				StressMatrix(2, j) = i(4, j);
				StressMatrix(3, j) = i(5, j);
				StressMatrix(4, j) = i(6, j);
				StressMatrix(5, j) = i(7, j);
			}
		}

		virtual void PrintData()
		{
			std::stringstream ss;
			ss << "RVE Constitutive Tensor:" << std::endl;

			ss << "Size: (" << mConstitutiveTensor.size1() << ", " <<mConstitutiveTensor.size2() << ") \n";

			ss << "Tensor: [" << mConstitutiveTensor << "]" << std::endl;

			std::cout << ss.str();
		}

		virtual void SetData()
		{
		}

		inline const IndexType& NumberOfDivision()const{ return mNumberOfDivision; }
		//inline IndexType&       NumberOfDivision()     { return mNumberOfDivision; }

		inline const IndexType& DamageSize()const{ return mDamageSize; }
		//inline IndexType&       DamageSize()     { return mDamageSize; }

		inline Matrix& GetConstitutiveTensor() { return mConstitutiveTensor; }
		inline const VectorVectorMap& RveElasticMap() { return mRveElasticMap; }
		inline const VectorMatrixMap& RveMaterialMap() { return mRveMaterialMap; }

		std::string GetInfo()const
        {
			/* JSON PRETTY FORMAT
			/ {
			/ 	"-23,-13": {									/ 1.		TAG
			/ 		"0.5": {									/ 1.1		DAMAGE
			/			"lambda": radius,						/ 1.1.1		LAMBDA -> STRAIN INTENSITY
			/ 			"stress": [Sxx,							/ 1.1.2		STRESS
			/ 					   Syy,
			/ 					   Szz,
			/ 					   Sxy,
			/ 					   Syz,
			/ 					   Sxz],
			/ 		}
			/		add mode damage ...
			/ 	}
			/	add more tags ...
			/ }
			*/

            std::stringstream ss;
            ss << "RVE Material Database:" << std::endl;

			ss << "Size: " << mRveMaterialMap.size() << " \n";

			ss << "Values: \n";
			for (VectorMatrixMap::const_iterator i = mRveMaterialMap.begin(); i != mRveMaterialMap.end(); ++i)
				ss << " " << i->first << ":" << i->second << std::endl; // TODO: Improve GetInfo() in order to return JSON PRETTY FORMAT
            ss << std::endl;
            return ss.str();
        }

	private:

		void ReadAndCreatePredictorData(const std::string& file_name)
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
				if (!(iss >> bracket >> damage >> bracket >> tag_1 >> comma >> tag_2 >> comma >> tag_3 >> comma >> tag_4 >> comma >> tag_5 >> e_xx >> e_yy >> e_zz >> e_xy >> e_yz >> e_xz >> s_xx >> s_yy >> s_zz >> s_xy >> s_yz >> s_xz)) //3D
				{
					std::cout << "Error during Reading the Hash File" << std::endl;
					break; // error
				}
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

				// Matrix S(3, 3, false);
				// S(0, 0) = 1.0;
				// S(1, 1) = 1.0;
				// S(2, 2) = 1.0;
				// S(3, 3) = 1.0;
				// Vector e = i_strain;
				// Vector Se = prod(S, e);
				// double eTSe = e(0)*Se(0) + e(1)*Se(1) + e(2)*Se(2) + e(3)*Se(3);
				//
				// double radius = sqrt(eTSe);
				double radius = norm_2(i_strain);
				// //Calculate phi & number of division
				// double m = 10.0; // m = max of tag
				// //std::cout << "m: " << m << std::endl;
				// double phi = Globals::Pi / m;
				// //std::cout << "phi: " << phi << " rad" << std::endl;
				//
				// Vector Theta(2, 0.0);
				// Theta = this->CalculateTheta(i_strain);
				// //std::cout << "Theta: " << Theta << std::endl;
				//
				// Vector real_tag = Theta / phi;
				// std::cout << "real_tag: " << real_tag << std::endl;

				r_and_strain_stress[0] = radius;
				r_and_strain_stress[1] = i_strain[0];
				r_and_strain_stress[2] = i_strain[1];
				r_and_strain_stress[3] = i_strain[2];
				r_and_strain_stress[4] = i_strain[3];
				r_and_strain_stress[5] = i_strain[4];
				r_and_strain_stress[6] = i_strain[5];
				r_and_strain_stress[7] = i_stress[0];
				r_and_strain_stress[8] = i_stress[1];
				r_and_strain_stress[9] = i_stress[2];
				r_and_strain_stress[10] = i_stress[3];
				r_and_strain_stress[11] = i_stress[4];
				r_and_strain_stress[12] = i_stress[5];

				mRveElasticMap.insert(VectorVectorMap::value_type(i_tag, r_and_strain_stress));
			}

			hashfile.close();
		}

		void ReadAndCreateConstTensorData(const std::string& file_name)
		{
			Matrix mC0;
			std::ifstream C0file(file_name.c_str());
			std::string line;
			int count_line = 0;
			while (std::getline(C0file, line)) //loop over lines of the .hash
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
				size_t strain_size = 6;
				if (!(iss >> bracket >> C11 >> comma >> C12 >> comma >> C13 >> comma >> C14 >> comma >> C15 >> comma >> C16 >> comma >> C21 >> comma >> C22 >> comma >> C23 >> comma >> C24 >> comma >> C25 >> comma >> C26 >> comma >> C31 >> comma >> C32 >> comma >> C33 >> comma >> C34 >> comma >> C35 >> comma >> C36 >> comma >> C41 >> comma >> C42 >> comma >> C43 >> comma >> C44 >> comma >> C45 >> comma >> C46 >> comma >> C51 >> comma >> C52 >> comma >> C53 >> comma >> C54 >> comma >> C55 >> comma >> C56 >> comma >> C61 >> comma >> C62 >> comma >> C63 >> comma >> C64 >> comma >> C65 >> comma >> C66 >> bracket))
				{
					std::cout << "Error during Reading the C0 File" << std::endl;
					break; // error
				}

				if (mC0.size1() != strain_size || mC0.size2() != strain_size)
					mC0.resize(strain_size, strain_size, false);
				noalias(mC0) = ZeroMatrix(strain_size, strain_size);

				mC0(0, 0) = C11; mC0(0, 1) = C12; mC0(0, 2) = C13; mC0(0, 3) = C14; mC0(0, 4) = C15; mC0(0, 5) = C16;
				mC0(1, 0) = C21; mC0(1, 1) = C22; mC0(1, 2) = C23; mC0(1, 3) = C24; mC0(1, 4) = C25; mC0(1, 5) = C26;
				mC0(2, 0) = C31; mC0(2, 1) = C32; mC0(2, 2) = C33; mC0(2, 3) = C34; mC0(2, 4) = C35; mC0(2, 5) = C36;
				mC0(3, 0) = C41; mC0(3, 1) = C42; mC0(3, 2) = C43; mC0(3, 3) = C44; mC0(3, 4) = C45; mC0(3, 5) = C46;
				mC0(4, 0) = C51; mC0(4, 1) = C52; mC0(4, 2) = C53; mC0(4, 3) = C54; mC0(4, 4) = C55; mC0(4, 5) = C56;
				mC0(5, 0) = C61; mC0(5, 1) = C62; mC0(5, 2) = C63; mC0(5, 3) = C64; mC0(5, 4) = C65; mC0(5, 5) = C66;
			}
			if (mConstitutiveTensor.size1() != mC0.size1() || mConstitutiveTensor.size2() != mC0.size2())
				mConstitutiveTensor.resize(mC0.size1(), mC0.size2(), false);
			noalias(mConstitutiveTensor) = mC0;

			C0file.close();
		}

		void ReadAndCreatePredictorDataDamageHistory(const std::string& file_name)
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

			size_t damage_size = 0;

			if (!ok)
			{
				std::stringstream msg;
				msg << rapidjson::GetParseError_En(ok.Code()) << " offset of the error from the beginning of the string = " << ok.Offset() << std::endl;
				msg << "a much more explicative error message can be obtained by analysing the input string with an online analyzer such for example json lint" << std::endl;
				msg << "the value of the string that was attempted to parse is :" << std::endl << std::endl;
				msg << file_name;
				KRATOS_THROW_ERROR(std::invalid_argument, "error found in parsing the json_string, the value of the json string was: \n", msg.str());
			}

			vector<int> i_tags(2);
			for (rapidjson::Value::ConstMemberIterator itr1 = document.MemberBegin(); itr1 != document.MemberEnd(); ++itr1) // Loop over TAGS
			{
				std::string i_tags_name = itr1->name.GetString();
				//std::cout << "i_tags_name: " << i_tags_name << std::endl;

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

				size_t i_D = 0;
				size_t itr_dim = document[i_tags_name.c_str()].Size();
				matrix<double> d_lambda_sigma(8, itr_dim, 0.0);
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

					//mDamageMap.insert(DamageMap::value_type(i_damage, lambda_sigma));

					d_lambda_sigma(0, i_D) = i_damage;
					d_lambda_sigma(1, i_D) = lambda_sigma[0];
					d_lambda_sigma(2, i_D) = lambda_sigma[1];
					d_lambda_sigma(3, i_D) = lambda_sigma[2];
					d_lambda_sigma(4, i_D) = lambda_sigma[3];
					d_lambda_sigma(5, i_D) = lambda_sigma[4];
					d_lambda_sigma(6, i_D) = lambda_sigma[5];
					d_lambda_sigma(7, i_D) = lambda_sigma[6];
					i_D += 1;
				}

				damage_size = i_D;
				//std::cout << "n_damage: " << n_damage << std::endl;
				Vector damagesVector_aux(damage_size, 0.0);
				for (size_t ii = 0; ii < damage_size; ii++)
					damagesVector_aux[ii] = d_lambda_sigma(0, ii);
				size_t idx_final = 0;
				Matrix dls(8, itr_dim, 0.0);
				for (auto i_sort : sort_indexes(damagesVector_aux)) { // only c++11
					dls(0, idx_final) = d_lambda_sigma(0, i_sort);
					dls(1, idx_final) = d_lambda_sigma(1, i_sort);
					dls(2, idx_final) = d_lambda_sigma(2, i_sort);
					dls(3, idx_final) = d_lambda_sigma(3, i_sort);
					dls(4, idx_final) = d_lambda_sigma(4, i_sort);
					dls(5, idx_final) = d_lambda_sigma(5, i_sort);
					dls(6, idx_final) = d_lambda_sigma(6, i_sort);
					dls(7, idx_final) = d_lambda_sigma(7, i_sort);
					idx_final += 1;
				}

				mRveMaterialMap.insert(VectorMatrixMap::value_type(i_tags, dls));
			}

			mDamageSize = damage_size;

			fclose(pFile);
		}


    protected:

		Matrix mConstitutiveTensor;
		VectorVectorMap mRveElasticMap; // TODO: Remove it
		VectorMatrixMap mRveMaterialMap;

		std::string mC0FileName;
		std::string mHashFileName; // TODO: Remove it
		std::string mJsonFileName;

		IndexType mNumberOfDivision;
		IndexType mDamageSize;

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const
		{
			rSerializer.save("strain", mConstitutiveTensor);
			//Serialize mRveElasticMap & mRveMaterialMap
		}

		virtual void load(Serializer& rSerializer)
		{
			rSerializer.load("strain", mConstitutiveTensor);
			//Serialize mRveElasticMap & mRveMaterialMap
		}

    };

	inline std::ostream & operator << (std::ostream& rOStream, const RveMaterialDatabase3D& rThis)
    {
		return rOStream << rThis.GetInfo();
    }

} // namespace Kratos

#endif // RVE_MATERIAL_DATABASE_3D_H_INCLUDED
