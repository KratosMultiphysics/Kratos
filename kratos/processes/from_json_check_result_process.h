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

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/result_dabatase.h"

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
 * @class FromJSONCheckResultProcess
 * @ingroup KratosCore
 * @brief This class is used in order to check results using a json file containing the solution a given model part with a certain frequency
 * @details This stores the dababase in a class denominated ResultDatabase which considers Table to store the information, therefore being able to interpolate results
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
    typedef Node NodeType;

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

    /**
     * @brief This function returns the error message
     * @return The error message
     */
    std::string GetErrorMessage()
    {
        return mErrorMessage;
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
     */
    void CheckNodeValues(IndexType& rCheckCounter);

    /**
     * @brief This method check the nodal historical values
     * @param rCheckCounter The check counter
     */
    void CheckNodeHistoricalValues(IndexType& rCheckCounter);

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

    /* Error message */
    std::string mErrorMessage = "";

    ///@name Private Operations
    ///@{

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
