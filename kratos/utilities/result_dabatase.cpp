//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/result_dabatase.h"

namespace Kratos
{

const EntityDatabase::GPDatabaseType& EntityDatabase::GetResultaData(const SizeType GPIndex) const
{
    KRATOS_DEBUG_ERROR_IF(GPIndex == mData.size()) << "Incompatible size. GPIndex: " << GPIndex << ". Size: " << mData.size() << std::endl;
    return mData[GPIndex];
}

/***********************************************************************************/
/***********************************************************************************/

double EntityDatabase::GetValue(
    const double Time,
    const SizeType ComponentIndex,
    const SizeType GPIndex
    ) const
{
    KRATOS_DEBUG_ERROR_IF(ComponentIndex == GetResultaData(GPIndex).size()) << "Incompatible size. ComponentIndex: " << ComponentIndex << ". Size: " << GetResultaData(GPIndex).size() << std::endl;
    return mData[GPIndex][ComponentIndex]->GetValue(Time);
}

/***********************************************************************************/
/***********************************************************************************/

void EntityDatabase::SetValues(
    const Vector& rValuesX,
    const Vector& rValuesY,
    const SizeType ComponentIndex,
    const SizeType GPIndex
    )
{
    auto& p_table = GetResultaData(GPIndex)[ComponentIndex];

    KRATOS_DEBUG_ERROR_IF(p_table == nullptr) << "No table defined for ComponentIndex: " << ComponentIndex << " GPIndex: " << GPIndex << std::endl;

    if (p_table->Data().size() > 0) {
        p_table->Clear(); // We clear to avoid reassign
    }

    KRATOS_ERROR_IF_NOT(rValuesX.size() == rValuesY.size()) << "The input vectors don't have the same size" << std::endl;

    // Push values
    for (IndexType i = 0; i < rValuesX.size(); ++i) {
        p_table->PushBack(rValuesX[i], rValuesY[i]);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void EntityDatabase::Clear()
{
    // Clear stored databases
    for (auto& r_database : mData) {
        r_database.clear();
    }

    // Clear the container
    mData.clear();
}

/***********************************************************************************/
/***********************************************************************************/

const EntityDatabase& VariableDatabase::GetEntityData(const IndexType EntityIndex) const
{
    KRATOS_DEBUG_ERROR_IF(EntityIndex == mData.size()) << "Incompatible size. EntityIndex: " << EntityIndex << ". Size: " << mData.size() << std::endl;
    return mData[EntityIndex];
}

/***********************************************************************************/
/***********************************************************************************/

double VariableDatabase::GetValue(
    const IndexType EntityIndex,
    const double Time,
    const SizeType ComponentIndex,
    const SizeType GPIndex
    ) const
{
    KRATOS_DEBUG_ERROR_IF(EntityIndex == mData.size()) << "Incompatible size. EntityIndex: " << EntityIndex << ". Size: " << mData.size() << std::endl;
    return mData[EntityIndex].GetValue(Time, ComponentIndex, GPIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void VariableDatabase::SetValues(
    const Vector& rValuesX,
    const Vector& rValuesY,
    const IndexType EntityIndex,
    const SizeType ComponentIndex,
    const SizeType GPIndex
    )
{
    auto& r_database = mData[EntityIndex];
    r_database.SetValues(rValuesX, rValuesY, ComponentIndex, GPIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void VariableDatabase::Clear()
{
    // Clear stored databases
    for (auto& r_database : mData) {
        r_database.Clear();
    }

    // Clear the container
    mData.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void ResultDatabase::Initialize(
    const std::vector<IndexType>& rVariablesIndexes,
    const std::vector<IndexType>& rValuesSizes,
    const SizeType NumberOfEntites,
    const SizeType NumberOfGP
    )
{
    // If no variables we skip
    if (rVariablesIndexes.size() == 0) {
        return void();
    }

    KRATOS_ERROR_IF_NOT(rVariablesIndexes.size() == rValuesSizes.size()) << "Inconsistent sizes in the values sizes and the variable indexes" << std::endl;

    // Auxiliar lambda to generate vectors of tables
    auto table_generator =[](const SizeType NumberOfEntites, const SizeType NumberOfComponents, const SizeType NumberOfGP){
        EntityDatabase::GPDatabaseType aux_1(NumberOfComponents, nullptr);
        EntityDatabase aux_2(NumberOfGP, aux_1);
        VariableDatabase data(NumberOfEntites, aux_2);
        for (IndexType k = 0; k < NumberOfEntites; ++k){
            for (IndexType j = 0; j < NumberOfGP; ++j){
                for (IndexType i = 0; i < NumberOfComponents; ++i){
                    data[k][j][i] = new Table<double, double>();
                }
            }
        };
        return data;
    };

    // Fill the inner map of tables
    for (IndexType i = 0; i < rVariablesIndexes.size(); ++i) {
        const IndexType index = rVariablesIndexes[i];
        const SizeType size = rValuesSizes[i];
        mData.insert(std::pair<IndexType, VariableDatabase>(index, table_generator(NumberOfEntites, size, NumberOfGP)));
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ResultDatabase::Clear()
{
    // Clear stored databases
    for (auto& r_data : mData) {
        for (auto& r_sub_data : r_data.second) {
            for (auto& r_sub_sub_data : r_sub_data) {
                std::for_each( r_sub_sub_data.begin(), r_sub_sub_data.end(), [](Table<double, double>* p){delete p;});
                r_sub_sub_data.clear();
            }
        }
        (r_data.second).Clear();
    }

    // Clear the container
    mData.clear();
    mCommonColumn.clear();
}

/***********************************************************************************/
/***********************************************************************************/

int ResultDatabase::Check()
{
    // Doing check in the table size
    for (const auto& r_pair : mData) {
        const auto& r_table_vector_vector_vector = r_pair.second;
        for (const auto& r_table_vector_vector : r_table_vector_vector_vector) {
            for (const auto& r_table_vector : r_table_vector_vector) {
                for (const auto& p_table : r_table_vector) {
                    KRATOS_ERROR_IF(p_table == nullptr) << "Table not created" << std::endl;
                    KRATOS_ERROR_IF_NOT(p_table->Data().size() == mCommonColumn.size()) << "Inconsistent size of the tables. Time vector size: " << mCommonColumn.size() << " vs table size " << p_table->Data().size() << std::endl;
                }
            }
        }
    }
    return 0;
}

}  // namespace Kratos.
