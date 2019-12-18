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

// External includes

// Project includes
#include "includes/model_part.h"
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
     */
    void KRATOS_API(KRATOS_CORE) AddMissingVariablesFromSpecifications(
        ModelPart& rModelPart,
        const Parameters SpecificationsParameters 
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
     */
    void KRATOS_API(KRATOS_CORE) AddMissingDofsFromSpecifications(
        ModelPart& rModelPart,
        const Parameters SpecificationsParameters 
        );

}; // namespace SpecificationsUtilities
}  // namespace Kratos
#endif /* KRATOS_SPECIFICATIONS_UTILITIES defined */
