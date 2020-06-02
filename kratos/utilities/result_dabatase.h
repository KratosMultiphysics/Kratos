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

#if !defined(KRATOS_RESULT_DATABASE_H_INCLUDED )
#define  KRATOS_RESULT_DATABASE_H_INCLUDED

// System includes

// External includes
#include <unordered_map>

// Project includes
#include "includes/table.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class EntityDatabase
 * @ingroup KratosCore
 * @brief This class stores the results of a entity
 * @author Vicente Mataix Ferrandiz
*/
class EntityDatabase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResultDatabase
    KRATOS_CLASS_POINTER_DEFINITION(EntityDatabase);

    /// GP database definition
    typedef std::vector<Table<double, double>*> GPDatabaseType;

    /// Base type definition
    typedef std::vector<GPDatabaseType> DataType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    EntityDatabase(
        const SizeType SizeVector,
        const GPDatabaseType& rBaseData
        ) : mData(SizeVector, rBaseData)
    {
    }

    /// Destructor.
    virtual ~EntityDatabase()
    {
        this->Clear();
    }

    /// Copy constructor.
    EntityDatabase(EntityDatabase const& rOther)
        : mData(rOther.mData)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief The [] operator
     * @param i The index of the value to access
     */
    GPDatabaseType& operator[](const std::size_t i)
    {
        return mData[i];
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief The begin iterator
     */
    DataType::iterator begin()
    {
        return mData.begin();
    }

    /**
     * @brief The const begin iterator
     */
    DataType::const_iterator begin() const
    {
        return mData.begin();
    }

    /**
     * @brief The end iterator
     */
    DataType::iterator end()
    {
        return mData.end();
    }

    /**
     * @brief The const end iterator
     */
    DataType::const_iterator end() const
    {
        return mData.end();
    }

    /**
     * @brief This method retrieves the entity database
     * @param EntityIndex The index of the entity (not the entity Id, the index in the database)
     */
    const GPDatabaseType& GetResultaData(const SizeType GPIndex = 0) const
    {
        KRATOS_DEBUG_ERROR_IF(GPIndex == mData.size()) << "Incompatible size. GPIndex: " << GPIndex << ". Size: " << mData.size() << std::endl;
        return mData[GPIndex];
    }

    /**
     * @brief This method retrieves the interpolated value from the database
     * @param Time The time value to be retrieved
     * @param ComponentIndex The component index of the vector/array
     * @param GPIndex The Gauss point index
     */
    double GetValue(
        const double Time,
        const SizeType ComponentIndex = 0,
        const SizeType GPIndex = 0
        ) const
    {
        KRATOS_DEBUG_ERROR_IF(ComponentIndex == GetResultaData(GPIndex).size()) << "Incompatible size. ComponentIndex: " << ComponentIndex << ". Size: " << GetResultaData(GPIndex).size() << std::endl;
        return mData[GPIndex][ComponentIndex]->GetValue(Time);
    }

    /**
     * @brief This method set the values into the tables
     * @param rValuesX The values of the X axis
     * @param rValuesY The values of the Y axis
     * @param ComponentIndex The component index of the vector/array
     * @param GPIndex The GP index
     */
    void SetValues(
        const Vector& rValuesX,
        const Vector& rValuesY,
        const SizeType ComponentIndex = 0,
        const SizeType GPIndex = 0
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
//             p_table->insert(rValuesX[i], rValuesY[i]);
        }
    }

    /**
     * @brief This function is designed for clear all the databases
     */
    void Clear()
    {
        // Clear stored databases
        for (auto& r_database : mData) {
            r_database.clear();
        }

        // Clear the container
        mData.clear();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "EntityDatabase";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "EntityDatabase";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    DataType mData; // The database storing the values

    ///@}

}; // Class EntityDatabase

/**
 * @class VariableDatabase
 * @ingroup KratosCore
 * @brief This class stores the results of a variable
 * @author Vicente Mataix Ferrandiz
*/
class VariableDatabase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResultDatabase
    KRATOS_CLASS_POINTER_DEFINITION(VariableDatabase);

    /// Base type definition
    typedef std::vector<EntityDatabase> DataType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    VariableDatabase(
        const SizeType SizeVector,
        const EntityDatabase& rBaseData
        ) : mData(SizeVector, rBaseData)
    {
    }

    /// Destructor.
    virtual ~VariableDatabase()
    {
        this->Clear();
    }

    /// Copy constructor.
    VariableDatabase(VariableDatabase const& rOther)
        : mData(rOther.mData)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief The [] operator
     * @param i The index of the value to access
     */
    EntityDatabase& operator[](const std::size_t i)
    {
        return mData[i];
    }

    ///@}
    ///@name Operations
    ///@{

    DataType::iterator begin()
    {
        return mData.begin();
    }

    DataType::const_iterator begin() const
    {
        return mData.begin();
    }

    DataType::iterator end()
    {
        return mData.end();
    }

    DataType::const_iterator end() const
    {
        return mData.end();
    }

    /**
     * @brief This method retrieves the entity database
     * @param EntityIndex The index of the entity (not the entity Id, the index in the database)
     */
    const EntityDatabase& GetEntityData(const IndexType EntityIndex) const
    {
        KRATOS_DEBUG_ERROR_IF(EntityIndex == mData.size()) << "Incompatible size. EntityIndex: " << EntityIndex << ". Size: " << mData.size() << std::endl;
        return mData[EntityIndex];
    }

    /**
     * @brief This method retrieves the interpolated value from the database
     * @param EntityIndex The index of the entity (not the entity Id, the index in the database)
     * @param Time The time value to be retrieved
     * @param ComponentIndex The component index of the vector/array
     * @param GPIndex The Gauss point index
     */
    double GetValue(
        const IndexType EntityIndex,
        const double Time,
        const SizeType ComponentIndex = 0,
        const SizeType GPIndex = 0
        ) const
    {
        KRATOS_DEBUG_ERROR_IF(EntityIndex == mData.size()) << "Incompatible size. EntityIndex: " << EntityIndex << ". Size: " << mData.size() << std::endl;
        return mData[EntityIndex].GetValue(Time, ComponentIndex, GPIndex);
    }

    /**
     * @brief This method set the values into the tables
     * @param rValuesX The values of the X axis
     * @param rValuesY The values of the Y axis
     * @param EntityIndex The index of the entity (not the entity Id, the index in the database)
     * @param ComponentIndex The component index of the vector/array
     * @param GPIndex The GP index
     */
    void SetValues(
        const Vector& rValuesX,
        const Vector& rValuesY,
        const IndexType EntityIndex,
        const SizeType ComponentIndex = 0,
        const SizeType GPIndex = 0
        )
    {
        auto& r_database = mData[EntityIndex];
        r_database.SetValues(rValuesX, rValuesY, ComponentIndex, GPIndex);
    }

    /**
     * @brief This function is designed for clear all the databases
     */
    void Clear()
    {
        // Clear stored databases
        for (auto& r_database : mData) {
            r_database.Clear();
        }

        // Clear the container
        mData.clear();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "VariableDatabase";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "VariableDatabase";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    DataType mData; // The database storing the values

    ///@}
}; // Class VariableDatabase

/**
 * @class ResultDatabase
 * @ingroup KratosCore
 * @brief This class stores the results of a simulation for a later comparison
 * @details The results are stored in a map which stores the results in tables. The following structure is followed:
 *  - The key is the Id of the variable considered (so it does not depend if a double or component)
 *  - A vector containing a set of tables. This tables are which actually contain the results
 * @note If the table could store N columns without requiring a template argument the vector would not be required, saving NvariablesxNColumnsxNsteps doubles of memory
 * @author Vicente Mataix Ferrandiz
*/
class ResultDatabase
    : public std::unordered_map<IndexType, VariableDatabase>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResultDatabase
    KRATOS_CLASS_POINTER_DEFINITION(ResultDatabase);

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ResultDatabase(){}

    /// Destructor.
    virtual ~ResultDatabase() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method initializes the database
     * @param rVariablesIndexes The list of indexes of the variables
     * @param NumberOfEntites The number of entities considered
     */
    void Initialize(
        const std::vector<IndexType>& rVariablesIndexes,
        const std::vector<IndexType>& rValuesSizes,
        const SizeType NumberOfEntites,
        const SizeType NumberOfGP = 1
        )
    {
        // If no variables we skip
        if (rVariablesIndexes.size() == 0) {
            return void();
        }

        KRATOS_ERROR_IF_NOT(rVariablesIndexes.size() == rValuesSizes.size()) << "Inconsistent sizes in the values sizes and the variable indexes" << std::endl;

        // Auxiliar lambda to generate vectors of tables
        auto table_generator =[](const SizeType NumberOfEntites, const SizeType NumberOfComponents, const SizeType NumberOfGP){std::vector<Table<double, double>*> aux_1(NumberOfComponents, nullptr); EntityDatabase aux_2(NumberOfGP, aux_1); VariableDatabase data(NumberOfEntites, aux_2); for (IndexType k = 0; k < NumberOfEntites; ++k){for (IndexType j = 0; j < NumberOfGP; ++j){ for (IndexType i = 0; i < NumberOfComponents; ++i){ data[k][j][i] = new Table<double, double>();}}}; return data;};

        // Fill the inner map of tables
        for (IndexType i = 0; i < rVariablesIndexes.size(); ++i) {
            const IndexType index = rVariablesIndexes[i];
            const SizeType size = rValuesSizes[i];
            this->insert(std::pair<IndexType, VariableDatabase>(index, table_generator(NumberOfEntites, size, NumberOfGP)));
        }
    }

    /**
     * @brief This method retrieves the variable database
     * @param rVariable The variable to be retrieved
     * @tparam TVariableType The variable type considered
     */
    template<class TVariableType>
    VariableDatabase& GetVariableData(const TVariableType& rVariable)
    {
        const auto it = find(rVariable.Key());
        if (it != end()) {
            return it->second;
        } else {
            KRATOS_ERROR << "Not allocated Variable: " << rVariable.Name() << std::endl;
        }
    }

    /**
     * @brief This method retrieves the variable database
     * @param rVariable The variable to be retrieved
     * @tparam TVariableType The variable type considered
     */
    template<class TVariableType>
    const VariableDatabase& GetVariableData(const TVariableType& rVariable) const
    {
        const auto it = find(rVariable.Key());
        if (it != end()) {
            return it->second;
        } else {
            KRATOS_ERROR << "Not allocated Variable: " << rVariable.Name() << std::endl;
        }
    }

    /**
     * @brief This method set the values into the tables
     * @param rValuesX The values of the X axis
     * @param rValuesY The values of the Y axis
     * @param rVariable The variable to be retrieved
     * @param EntityIndex The index of the entity (not the entity Id, the index in the database)
     * @param ComponentIndex The component index of the vector/array
     * @param GPIndex The GP index
     * @tparam TVariableType The variable type considered
     */
    template<class TVariableType>
    void SetValues(
        const Vector& rValuesX,
        const Vector& rValuesY,
        const TVariableType& rVariable,
        const IndexType EntityIndex,
        const SizeType ComponentIndex = 0,
        const SizeType GPIndex = 0
        )
    {
        auto it = find(rVariable.Key());
        if (it != end()) {
            auto& r_database = it->second;
            r_database.SetValues(rValuesX, rValuesY, EntityIndex, ComponentIndex, GPIndex);
        } else {
            KRATOS_ERROR << "Not allocated Variable: " << rVariable.Name() << std::endl;
        }
    }

    /**
     * @brief This function is designed for clear all the databases
     */
    void Clear()
    {
        this->clear();
        mCommonColumn.clear();
    }

    /**
     * @brief This function is designed for being called after ExecuteInitialize ONCE to verify that the input is correct.
     */
    int Check()
    {
        // Doing check in the table size
        for (const auto& r_pair : *this) {
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

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the common column vector
     * @return The common column vector
     */
    Vector& GetCommonColumn()
    {
        return mCommonColumn;
    }

    /**
     * @brief This method sets the common column
     * @param rCommonColumn The common column (time) vector
     */
    void SetCommonColumn(const Vector& rCommonColumn)
    {
        mCommonColumn = rCommonColumn;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ResultDatabase";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ResultDatabase";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    Vector mCommonColumn; /// This vector stores the common column (usually TIME), so it needs to be initialized at the begining

    ///@}

}; // Class ResultDatabase

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_RESULT_DATABASE_H_INCLUDED  defined
