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

#if !defined(KRATOS_SPECIFICATIONS_UTILITIES)
#define KRATOS_SPECIFICATIONS_UTILITIES

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "geometries/geometry_data.h"

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
 * @namespace SpecificationsUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for evaluate specifications
 * @author Vicente Mataix Ferrandiz
 */
namespace SpecificationsUtilities
{
    /// A definition of the component variable
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > Component3VarType;
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 4> > > Component4VarType;
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 6> > > Component6VarType;
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 9> > > Component9VarType;
    
    /**
     * @brief This method adds to the model part the missing variables
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    void KRATOS_API(KRATOS_CORE) AddMissingVariables(ModelPart& rModelPart);
    
    /**
     * @brief This method adds to the model part the missing variables from a given set of specifications
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param SpecificationsParameters The specification parameters
     * @param EntityName The name of the entity considered
     */
    void KRATOS_API(KRATOS_CORE) AddMissingVariablesFromSpecifications(
        ModelPart& rModelPart,
        const Parameters SpecificationsParameters,
        const std::string EntityName = "NOT_DEFINED"
        );

    /**
     * @brief This method adds to the model part the missing dofs
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    void KRATOS_API(KRATOS_CORE) AddMissingDofs(ModelPart& rModelPart);
    
    /**
     * @brief This method adds to the model part the missing dofs from a given set of specifications
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param SpecificationsParameters The specification parameters
     * @param EntityName The name of the entity considered
     */
    void KRATOS_API(KRATOS_CORE) AddMissingDofsFromSpecifications(
        ModelPart& rModelPart,
        const Parameters SpecificationsParameters,
        const std::string EntityName = "NOT_DEFINED"
        );
    
    /**
     * @brief This method determine the flags used on the simulation
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    void KRATOS_API(KRATOS_CORE) DetermineFlagsUsed(ModelPart& rModelPart);
    
    /**
     * @brief This method detects the framework considered. It throws a warning if incompatible framework
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return The framework of the problem
     */
    std::string KRATOS_API(KRATOS_CORE) DetermineFramework(ModelPart& rModelPart);
    
    /**
     * @brief This method detects if the LHS is symmetric. It throws a warning if incompatible
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the LHS are symmetric, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) DetermineSymmetricLHS(ModelPart& rModelPart);
    
    /**
     * @brief This method detects if the LHS is positive definite. It throws a warning if incompatible 
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the LHS are positive definite, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) DeterminePositiveDefiniteLHS(ModelPart& rModelPart);
    
    /**
     * @brief This method detects if the elements/conditions are compatible with its geometry
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all geometries are compatible, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) DetermineIfCompatibleGeometries(ModelPart& rModelPart);
    
    /**
     * @brief This method detects if all elements/conditions are implicit. It throws a warning if incompatible
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the elements/conditions are implicit, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) DetermineIfImplicitSimulation(ModelPart& rModelPart);
    
    /**
     * @brief This method detects if all elements/conditions require time integration. It throws a warning if incompatible
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the elements/conditions require time integration, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) DetermineIfRequiresTimeIntegration(ModelPart& rModelPart);

    /**
     * @brief This method detects if all elements/conditions are considering the proper CL
     * @param rModelPart Reference to the ModelPart containing the problem
     * @return True if all the elements/conditions are considering the adequate CL, false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) CheckCompatibleConstitutiveLaws(ModelPart& rModelPart);
    
}; // namespace SpecificationsUtilities
}  // namespace Kratos
#endif /* KRATOS_SPECIFICATIONS_UTILITIES defined */
