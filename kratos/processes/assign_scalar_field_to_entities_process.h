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
//  Main authors:    Riccardo Rossi
//                   Josep Maria Carbonell
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/function_parser_utility.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief This struct is used in order to identify when using the historical and non historical variables
 */
struct AssignScalarFieldToEntitiesProcessSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @class AssignScalarFieldToEntitiesProcess
 * @ingroup KratosCore
 * @brief This process is used in order to assign a function to a entity
 * @details This function assigns a value depending on a function to a variable in all of the conditions in a given mesh.
 * The behaviour is the following:
 * - Option 1 - Variable<double>: It is evaluated in the center of the entities
 * - Option 2 - array_1d_component_type: The same as Variable evaluated in the center of the element
 * - Option 3 - Variable< Vector > : The vector has to have a size equal to the number of nodes, and its values are computed per each entry of the vector using the coordinates of the nodes which occupies the same position in the geometry
 * @author Riccardo Rossi
 * @author Josep Maria Carbonell
 * @author Vicente Mataix Ferrandiz
*/
template<class TEntity, bool THistorical = AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>
class KRATOS_API(KRATOS_CORE) AssignScalarFieldToEntitiesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the geometry
    using GeometryType = Geometry<Node>;

    /// The IndexType definition
    using IndexType = std::size_t;

    /// The SizeType definition
    using SizeType = std::size_t;

    /// The container of the entities
    using EntityContainerType = PointerVectorSet<TEntity, IndexedObject>;

    /// Pointer definition of AssignScalarFieldToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarFieldToEntitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The default constructor
     * @param rModel The model where the scalar field will be applied
     * @param rParameters The configuration parameters
     */
    AssignScalarFieldToEntitiesProcess(
        Model& rModel,
        Parameters rParameters
        );

    /**
     * @brief The default constructor
     * @param rModelPart The model part where the scalar field will be applied
     * @param rParameters The configuration parameters
     */
    AssignScalarFieldToEntitiesProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        );

    /// Destructor.
    ~AssignScalarFieldToEntitiesProcess() override {}

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
    ///@{.

    /**
     * @brief Execute method is used to execute the AssignScalarFieldToEntitiesProcess algorithms.
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
        return "AssignScalarFieldToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarFieldToEntitiesProcess";
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

    ModelPart& mrModelPart;                                /// The modelpart where compute
    Kratos::unique_ptr<GenericFunctionUtility> mpFunction; /// The python function used, depends on X, Y, Z, and t (it could also depend on X0, Y0, Z0)
    std::string mVariableName;                             /// The name of the variable to assign

    ///@}
    ///@name Private Operators
    ///@{

    /// Assignment operator.
    AssignScalarFieldToEntitiesProcess& operator=(AssignScalarFieldToEntitiesProcess const& rOther);

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief It calls the function for a vector variable
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     */
    void CallFunction(
        const typename TEntity::Pointer pEntity,
        const double Time,
        Vector& rValue,
        GenericFunctionUtility& rFunction
        );

    /**
     * @brief It calls the function for components
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     * @param rFunction The function to call
     */
    void CallFunctionComponents(
        const typename TEntity::Pointer pEntity,
        const double Time,
        double& rValue,
        GenericFunctionUtility& rFunction
        );

    /**
     * @brief It calls the function (local system)
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     * @param rFunction The function to call
     */
    void CallFunctionLocalSystem(
        const typename TEntity::Pointer pEntity,
        const double Time,
        Vector& rValue,
        GenericFunctionUtility& rFunction
        );

    /**
     * @brief It calls the function (local system) for components
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     * @param rFunction The function to call
     */
    void CallFunctionLocalSystemComponents(
        const typename TEntity::Pointer pEntity,
        const double Time,
        double& rValue,
        GenericFunctionUtility& rFunction
        );

    /**
     * @brief It assigns time dependency
     * @param pEntity The pointer to the condition where set the function
     * @param Time The current time
     * @param rValue The value to set
     * @param Value The dependency value
     */
    void AssignTimeDependentValue(
        const typename TEntity::Pointer pEntity,
        const double Time,
        Vector& rValue,
        const double Value
        );

    /**
     * @brief This is the methods that set the values globally (tries all the possible options)
     * @param rVar The variable to set
     * @param Time The current time
     */
    void InternalAssignValueVector(
        const Variable<Vector>& rVar,
        const double Time
        );

    /**
     * @brief This is the methods that set the values globally (tries all the possible options) for components
     * @param rVar The variable to set
     * @param Time The current time
     */
    void InternalAssignValueScalar(
        const Variable<double>& rVar,
        const double Time
        );

    /**
     * @brief This method returns the current entity container
     * @return The current entity container
     */
    EntityContainerType& GetEntitiesContainer();

    ///@}
}; // Class AssignScalarFieldToEntitiesProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TEntity, bool THistorical>
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarFieldToEntitiesProcess<TEntity, THistorical>& rThis);

/// output stream function
template<class TEntity, bool THistorical>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarFieldToEntitiesProcess<TEntity, THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
