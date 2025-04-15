//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ASSIGN_SCALAR_INPUT_TO_ENTITIES_PROCESS_H_INCLUDED )
#define  KRATOS_ASSIGN_SCALAR_INPUT_TO_ENTITIES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/result_dabatase.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @brief This struct is used in order to identify when using the hitorical and non historical variables
 */
struct AssignScalarInputToEntitiesProcessSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @class AssignScalarInputToEntitiesProcess
 * @ingroup KratosCore
 * @brief This function assigns a value from an input to a variable belonging to all of the entities in a given mesh
 * @details Can be used to any entities
 * @tparam TEntity The entity type
 * @author Vicente Mataix Ferrandiz
*/
template<class TEntity, bool THistorical = false>
class KRATOS_API(KRATOS_CORE) AssignScalarInputToEntitiesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Node type definition
    typedef Node NodeType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// The container of the entities
    typedef PointerVectorSet<TEntity, IndexedObject> EntityContainerType;

    /// Pointer definition of AssignScalarInputToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarInputToEntitiesProcess);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( GEOMETRIC_DEFINITION ); /// This flag is used in order to check if the definition is defined in the geometry or in the entities

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief This enum helps us to identify the algorithm considered
     */
    enum class Algorithm {
        NEAREST_NEIGHBOUR = 0
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModelPart The model part to be set
     * @param rParameters The configuration parameters
     */
    AssignScalarInputToEntitiesProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        );

    /// Destructor.
    ~AssignScalarInputToEntitiesProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AssignScalarInputToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarInputToEntitiesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;                                                  /// The model part where to assign the values

    const Variable<double>* mpVariable = nullptr;                            /// The pointer of the variable

    ResultDatabase mDatabase;                                                /// The database containing the information to assign the values

    std::vector<std::unordered_map<IndexType, double>> mWeightExtrapolation; /// This maps specifies the extrapolation weights for each entity

    std::vector<array_1d<double, 3>> mCoordinates;                           /// The coordinates of the database

    Algorithm mAlgorithm = Algorithm::NEAREST_NEIGHBOUR;                     /// The algorithm considered

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method identifies the database from a TXT
     * @param rFileName The name of the TXT file name
     */
    void IdentifyDataTXT(const std::string& rFileName);

    /**
     * @brief This method identifies the database from a JSON
     * @param rFileName The name of the JSON file name
     */
    void IdentifyDataJSON(const std::string& rFileName);

    /**
     * @brief This method reads the database from a TXT
     * @param rFileName The name of the TXT file name
     */
    void ReadDataTXT(const std::string& rFileName);

    /**
     * @brief This method reads the database from a JSON
     * @param rFileName The name of the JSON file name
     */
    void ReadDataJSON(const std::string& rFileName);

    /**
     * @brief It computes the extrapolation weights
     */
    void ComputeExtrapolationWeight();

    /**
     * @brief This method assigns the value
     * @param rVariable The variable to be assigned
     * @param Value The value to assign
     */
    void InternalAssignValue(
        const Variable<double>& rVariable,
        const double Value
        );

    /**
     * @brief This converts the algorithm string to an enum
     * @param Str The string that you want to convert in the equivalent enum
     * @return Algorithm: The equivalent enum (this requires less memory and is eassier to compare than a std::string)
     */
    Algorithm ConvertAlgorithmString(const std::string& Str)
    {
        if(Str == "NEAREST_NEIGHBOUR" || Str == "nearest_neighbour")
            return Algorithm::NEAREST_NEIGHBOUR;
        else
            return Algorithm::NEAREST_NEIGHBOUR;
    }

    /**
     * @brief This method returns the current entity label
     * @param Id The id of the entity
     * @return The current entity label
     */
    array_1d<double, 3> GetCoordinatesEntity(const IndexType Id);

    /**
     * @brief This method returns the current entity container
     * @return The current entity container
     */
    EntityContainerType& GetEntitiesContainer();

    /**
     * @brief This method resets values in the entities
     */
    void ResetValues();

    /**
     * @brief This method sets values in a entity
     * @param rEntity The entity reference
     * @param rVariable The variable to set
     * @param Value The value to set
     */
    void SetValue(
        TEntity& rEntity,
        const Variable<double>& rVariable,
        const double Value
        );

    /**
     * @brief This method gets values in a entity
     * @param rEntity The entity reference
     * @param rVariable The variable to set
     * @return The value to get
     */
    double& GetValue(
        TEntity& rEntity,
        const Variable<double>& rVariable
        );

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarInputToEntitiesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TEntity>
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarInputToEntitiesProcess<TEntity>& rThis);

/// output stream function
template<class TEntity>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarInputToEntitiesProcess<TEntity>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_INPUT_TO_ENTITIES_PROCESS_H_INCLUDED  defined
