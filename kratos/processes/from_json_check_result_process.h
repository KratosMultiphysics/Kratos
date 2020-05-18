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

#if !defined(KRATOS_FROM_JSON_CHECK_RESULT_PROCESS_H_INCLUDED )
#define  KRATOS_FROM_JSON_CHECK_RESULT_PROCESS_H_INCLUDED

// System includes

// External includes
#include <unordered_map>

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/table.h"
#include "includes/kratos_parameters.h"

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
    : public std::vector<std::vector<Table<double, double>*>>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResultDatabase
    KRATOS_CLASS_POINTER_DEFINITION(EntityDatabase);

    /// GP database definition
    typedef std::vector<Table<double, double>*> GPDatabaseType;

    /// Base type definition
    typedef std::vector<GPDatabaseType> BaseType;

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
        ) : BaseType(SizeVector, rBaseData)
    {
    }

    /// Destructor.
    virtual ~EntityDatabase() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method retrieves the entity database
     * @param EntityIndex The index of the entity (not the entity Id, the index in the database)
     */
    const GPDatabaseType& GetResultaData(const SizeType GPIndex = 0) const
    {
        KRATOS_DEBUG_ERROR_IF(GPIndex == this->size()) << "Incompatible size. GPIndex: " << GPIndex << ". Size: " << this->size() << std::endl;
        return operator[](GPIndex);
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
        return operator[](GPIndex)[ComponentIndex]->GetValue(Time);
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

}; // Class EntityDatabase

/**
 * @class VariableDatabase
 * @ingroup KratosCore
 * @brief This class stores the results of a variable
 * @author Vicente Mataix Ferrandiz
*/
class VariableDatabase
    : public std::vector<EntityDatabase>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResultDatabase
    KRATOS_CLASS_POINTER_DEFINITION(VariableDatabase);

    /// Base type definition
    typedef std::vector<EntityDatabase> BaseType;

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
        ) : BaseType(SizeVector, rBaseData)
    {
    }

    /// Destructor.
    virtual ~VariableDatabase() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method retrieves the entity database
     * @param EntityIndex The index of the entity (not the entity Id, the index in the database)
     */
    const EntityDatabase& GetEntityData(const IndexType EntityIndex) const
    {
        KRATOS_DEBUG_ERROR_IF(EntityIndex == this->size()) << "Incompatible size. EntityIndex: " << EntityIndex << ". Size: " << this->size() << std::endl;
        return operator[](EntityIndex);
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
        KRATOS_DEBUG_ERROR_IF(EntityIndex == this->size()) << "Incompatible size. EntityIndex: " << EntityIndex << ". Size: " << this->size() << std::endl;
        return operator[](EntityIndex).GetValue(Time, ComponentIndex, GPIndex);
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
        auto& r_database = operator[](EntityIndex);
        r_database.SetValues(rValuesX, rValuesY, ComponentIndex, GPIndex);
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

/**
 * @class FromJSONCheckResultProcess
 * @ingroup KratosCore
 * @brief A
 * @details F
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) FromJSONCheckResultProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FromJSONCheckResultProcess
    KRATOS_CLASS_POINTER_DEFINITION(FromJSONCheckResultProcess);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( CORRECT_RESULT );                 /// This flag is used in order to check that the result is correct
    KRATOS_DEFINE_LOCAL_FLAG( HISTORICAL_VALUE );               /// This flag is used in order to check if the values are historical
    KRATOS_DEFINE_LOCAL_FLAG( CHECK_ONLY_LOCAL_ENTITIES );      /// This flag is used in order to check only local entities
    KRATOS_DEFINE_LOCAL_FLAG( NODES_CONTAINER_INITIALIZED );    /// This flag is used in order to check that nodes container are initialized
    KRATOS_DEFINE_LOCAL_FLAG( ELEMENTS_CONTAINER_INITIALIZED ); /// This flag is used in order to check that elements container are initialized
    KRATOS_DEFINE_LOCAL_FLAG( NODES_DATABASE_INITIALIZED );     /// This flag is used in order to check that nodes database are initialized
    KRATOS_DEFINE_LOCAL_FLAG( ELEMENTS_DATABASE_INITIALIZED );  /// This flag is used in order to check that elements database are initialized

    /// Containers definition
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ElementsContainerType        ElementsArrayType;

    /// The node type definiton
    typedef Node<3> NodeType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModel The model where the where the simulation is performed
     * @param ThisParameters The parameters of configuration
     */
    FromJSONCheckResultProcess(
        Model& rModel,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief Default constructor.
     * @param rModelPart The model part where the simulation is performed
     * @param ThisParameters The parameters of configuration
     */
    FromJSONCheckResultProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    virtual ~FromJSONCheckResultProcess() {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This function is designed for being called at the end of the computations
     */
    void ExecuteFinalize() override;

    /**
     * @brief This function is designed for being called after ExecuteInitialize ONCE to verify that the input is correct.
     */
    int Check() override;

    /**
     * @brief This function returns if the result is correct
     * @return If the result is correct
     */
    bool IsCorrectResult()
    {
        return this->Is(CORRECT_RESULT);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FromJSONCheckResultProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FromJSONCheckResultProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    /**
     * @brief This initializes the databases
     */
    void InitializeDatabases();

    /**
     * @brief This method fills the list of variables
     * @param rNodalVariablesNames The names of the nodal variables
     * @param rGPVariablesNames The names of the GP variables
     */
    void FillVariablesList(
        const std::vector<std::string>& rNodalVariablesNames,
        const std::vector<std::string>& rGPVariablesNames
        );

    /**
     * @brief This method checks if a flag is active in a given entity
     * @param rEntity The entity to check
     * @param pFLag The pointer to the flag to check
     */
    template<class TEntity>
    bool CheckFlag(
        const TEntity& rEntity,
        const Flags* pFlag
        )
    {
        if (pFlag != nullptr) {
            if (rEntity.IsNot(*pFlag)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief This method checks the results
     * @param ValueEntity The value on the entity
     * @param ValueJSON The reference value from the JSON
     */
    bool CheckValues(
        const double ValueEntity,
        const double ValueJSON
        );

    /**
     * @brief This returns a message in case of fail
     * @param EntityId The Kratos node or element to check
     * @param rEntityType The type of the entity
     * @param ValueEntity The value on the entity
     * @param ValueJSON The reference value from the json
     * @param rVariableName The name of the variable
     * @param ComponentIndex The component index
     * @param GPIndex The GP index
     */
    void FailMessage(
        const IndexType EntityId,
        const std::string& rEntityType,
        const double ValueEntity,
        const double ValueJSON,
        const std::string& rVariableName,
        const int ComponentIndex = -1,
        const int GPIndex = -1
        );

    /**
     * @brief This method check the nodal values
     * @param rCheckCounter The check counter
     * @tparam THistorical If the value is historical or not
     */
    template<bool THistorical>
    void CheckNodeValues(IndexType& rCheckCounter)
    {
        // Get time
        const double time = mrModelPart.GetProcessInfo().GetValue(TIME);

        // Node database
        const auto& r_node_database = GetNodeDatabase();

        // Iterate over nodes
        const auto& r_nodes_array = GetNodes();
        const auto it_node_begin = r_nodes_array.begin();

        // Auxiliar check counter (MVSC does not accept references)
        IndexType check_counter = rCheckCounter;

        for (auto& p_var_double : mpNodalVariableDoubleList) {
            const auto& r_var_database = r_node_database.GetVariableData(*p_var_double);

            #pragma omp parallel for reduction(+:check_counter)
            for (int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;

                const double result = GetValue<THistorical>(it_node, p_var_double);
                const double reference = r_var_database.GetValue(i, time);
                if (!CheckValues(result, reference)) {
                    FailMessage(it_node->Id(), "Node", result, reference, p_var_double->Name());
                    check_counter += 1;
                }
            }
        }
        for (auto& p_var_array : mpNodalVariableArrayList) {
            const auto& r_var_database = r_node_database.GetVariableData(*p_var_array);

            #pragma omp parallel for reduction(+:check_counter)
            for (int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;

                const auto& r_entity_database = r_var_database.GetEntityData(i);
                const array_1d<double, 3>& r_result = GetValue<THistorical>(it_node, p_var_array);
                for (IndexType i_comp = 0; i_comp < 3; ++i_comp) {
                    const double reference = r_entity_database.GetValue(time, i_comp);
                    if (!CheckValues(r_result[i_comp], reference)) {
                        FailMessage(it_node->Id(), "Node", r_result[i_comp], reference, p_var_array->Name());
                        check_counter += 1;
                    }
                }
            }
        }
        for (auto& p_var_vector : mpNodalVariableVectorList) {
            const auto& r_var_database = r_node_database.GetVariableData(*p_var_vector);

            #pragma omp parallel for reduction(+:check_counter)
            for (int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
                auto it_node = it_node_begin + i;

                const auto& r_entity_database = r_var_database.GetEntityData(i);
                const Vector& r_result = GetValue<THistorical>(it_node, p_var_vector);
                for (IndexType i_comp = 0; i_comp < r_result.size(); ++i_comp) {
                    const double reference = r_entity_database.GetValue(time, i_comp);
                    if (!CheckValues(r_result[i_comp], reference)) {
                        FailMessage(it_node->Id(), "Node", r_result[i_comp], reference, p_var_vector->Name());
                        check_counter += 1;
                    }
                }
            }
        }

        // Save the reference
        rCheckCounter = check_counter;
    }

    /**
     * @brief This method check the GP values
     * @param rCheckCounter The check counter
     */
    void CheckGPValues(IndexType& rCheckCounter);

    /**
     * @brief
     */
    SizeType SizeDatabase(
        const Parameters& rResults,
        const NodesArrayType& rNodesArray,
        const ElementsArrayType& rElementsArray
        );

    /**
     * @brief
     */
    void FillDatabase(
        const Parameters& rResults,
        const NodesArrayType& rNodesArray,
        const ElementsArrayType& rElementsArray,
        const SizeType NumberOfGP
        );

    /**
     * @brief Returns the identifier/key for saving nodal results in the json this can be either the node Id or its coordinates
     * @details The coordinates can be used to check the nodal results in MPI
     * @param rNode The Kratos node to get the identifier for
     */
    std::string GetNodeIdentifier(NodeType& rNode);

    /**
     * @brief This method returns the nodes of the model part
     * @return The nodes of the model part
     */
    NodesArrayType& GetNodes(const Flags* pFlag = nullptr);

    /**
     * @brief This method returns the elements of the model part
     * @return The elements of the model part
     */
    ElementsArrayType& GetElements(const Flags* pFlag = nullptr);

    /**
     * @brief This method computes the relevant digits to take into account
     */
    std::size_t ComputeRelevantDigits(const double Value);

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Protected  Access
    ///@{

    /**
     * @brief This method returns the model part
     * @return The model part of the problem
     */
    const ModelPart& GetModelPart() const;

    /**
     * @brief This method returns the settings
     * @return The settings of the problem
     */
    const Parameters GetSettings() const;

    /**
     * @brief This method returns the Nodes database. If not initialized it will try initialize again
     * @return The nodes database
     */
    const ResultDatabase& GetNodeDatabase();

    /**
     * @brief This method returns the GP database. If not initialized it will try initialize again
     * @return The GP database
     */
    const ResultDatabase& GetGPDatabase();

    ///@}
    ///@name Protected LifeCycle
    ///@{

    /// Protected constructor with modified default settings to be defined by derived class.
    FromJSONCheckResultProcess(ModelPart& rModelPart, Parameters Settings, Parameters DefaultSettings);

    ///@}

private:
    ///@name Member Variables
    ///@{

    /* Model part and different settings */
    ModelPart& mrModelPart;         /// The main model part
    Parameters mThisParameters;     /// The parameters (can be used for general pourposes)

    /* Additional values */
    double mFrequency;              /// The check frequency
    double mRelativeTolerance;      /// The relative tolerance
    double mAbsoluteTolerance;      /// The absolute tolerance
    SizeType mRelevantDigits;       /// This is the number of relevant digits

    /* Counters */
    double mTimeCounter = 0.0;      /// A time counter

    /* The entities of the containers */
    NodesArrayType mNodesArray;       /// The nodes of study
    ElementsArrayType mElementsArray; /// The elements of study

    /* The vectors storing the variables of study */
    std::vector<const Variable<double>*> mpNodalVariableDoubleList;                     /// The scalar variable list to compute
    std::vector<const Variable<array_1d<double,3>>*> mpNodalVariableArrayList;          /// The array variable list to compute
    std::vector<const Variable<Vector>*> mpNodalVariableVectorList;                     /// The vector variable list to compute

    std::vector<const Variable<double>*> mpGPVariableDoubleList;                        /// The scalar variable list to compute
    std::vector<const Variable<array_1d<double,3>>*> mpGPVariableArrayList;             /// The array variable list to compute
    std::vector<const Variable<Vector>*> mpGPVariableVectorList;                        /// The vector variable list to compute

    /* The databases which store the values */
    ResultDatabase mDatabaseNodes;  /// The database containing the information to compare the results for the nodes
    ResultDatabase mDatabaseGP;     /// The database containing the information to compare the results for the Gauss Points

    ///@name Private Operations
    ///@{

    /**
     * @brief This gets the double value
     * @param itNode Node iterator
     * @param pVariable The double variable
     * @tparam THistorical If the value is historical or not
     */
    template<bool THistorical>
    double GetValue(
        NodesArrayType::const_iterator& itNode,
        const Variable<double>* pVariable
        );

    /**
     * @brief This gets the array value
     * @param itNode Node iterator
     * @param pVariable The array variable
     * @tparam THistorical If the value is historical or not
     */
    template<bool THistorical>
    const array_1d<double, 3>& GetValue(
        NodesArrayType::const_iterator& itNode,
        const Variable<array_1d<double, 3>>* pVariable
        );

    /**
     * @brief This gets the vector value
     * @param itNode Node iterator
     * @param pVariable The vector variable
     * @tparam THistorical If the value is historical or not
     */
    template<bool THistorical>
    const Vector& GetValue(
        NodesArrayType::const_iterator& itNode,
        const Variable<Vector>* pVariable
        );

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FromJSONCheckResultProcess& operator=(FromJSONCheckResultProcess const& rOther);

    ///@}

}; // Class FromJSONCheckResultProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FromJSONCheckResultProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FromJSONCheckResultProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FROM_JSON_CHECK_RESULT_PROCESS_H_INCLUDED  defined
