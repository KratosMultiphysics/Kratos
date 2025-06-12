//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief This struct is used in order to identify when using the historical and non historical variables
 */
struct AssignScalarVariableToEntitiesProcessSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @class AssignScalarVariableToEntitiesProcess
 * @ingroup KratosCore
 * @brief This function assigns a value to a variable belonging to all of the entities in a given mesh
 * @details Can be used to any entities
 * @tparam TEntity The entity type
 * @author Josep Maria Carbonell
 * @author Vicente Mataix Ferrandiz
*/
template<class TEntity, bool THistorical = AssignScalarVariableToEntitiesProcessSettings::SaveAsNonHistoricalVariable>
class KRATOS_API(KRATOS_CORE) AssignScalarVariableToEntitiesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The container of the entities
    using EntityContainerType = PointerVectorSet<TEntity, IndexedObject>;

    /// Pointer definition of AssignScalarVariableToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarVariableToEntitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModel The model part to be set
     * @param rParameters The configuration parameters
     */
    AssignScalarVariableToEntitiesProcess(
        Model& rModel,
        Parameters rParameters
        );

    /**
     * @brief Default constructor
     * @param rModelPart The model part to be set
     * @param rParameters The configuration parameters
     */
    AssignScalarVariableToEntitiesProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        );

    /// Destructor.
    ~AssignScalarVariableToEntitiesProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the AssignScalarVariableToEntitiesProcess algorithms.
     */
    void Execute() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override
    {
        Execute();
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AssignScalarVariableToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarVariableToEntitiesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;     /// The model part where to assign the values
    std::string mVariableName;  /// The name of the variable
    double mDoubleValue;        /// The double value to assign
    int mIntValue;              /// The integer value to assign
    bool mBoolValue;            /// The boolean value to assign

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method assigns the value (with OMP)
     * @param rVar The variable to be assigned
     * @param Value The value to assign
     * @tparam TVarType The variable data type
     */
    template<class TVarType>
    void InternalAssignValue(
        TVarType& rVar,
        const typename TVarType::Type Value
        );

    /**
     * @brief This method returns the current entity container
     * @return The current entity container
     */
    EntityContainerType& GetEntitiesContainer();

    ///@}
}; // Class AssignScalarVariableToEntitiesProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TEntity, bool THistorical>
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarVariableToEntitiesProcess<TEntity, THistorical>& rThis);

/// output stream function
template<class TEntity, bool THistorical>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarVariableToEntitiesProcess<TEntity, THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.