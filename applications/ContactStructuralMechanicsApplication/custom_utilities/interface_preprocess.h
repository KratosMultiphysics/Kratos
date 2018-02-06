// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 


#if !defined(KRATOS_INTERFACE_PREPROCESS_CONDITION_H_INCLUDED )
#define  KRATOS_INTERFACE_PREPROCESS_CONDITION_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// TODO: Add parallellization!!!

namespace Kratos
{
///@name Kratos Globals
///@{
    
///@}
///@name Type Definitions
///@{
    
    typedef Point                                     PointType;
    typedef Node<3>                                    NodeType;
    typedef Geometry<NodeType>                     GeometryType;
    typedef Geometry<PointType>               GeometryPointType;
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{
  
/** \brief InterfacePreprocessCondition 
 * Creates Model Parts containing the interface
 */
class InterfacePreprocessCondition
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef ModelPart::NodesContainerType                   NodesArrayType;
    typedef ModelPart::ElementsContainerType             ElementsArrayType;
    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;
    
    /// Pointer definition of ExactMortarIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(InterfacePreprocessCondition);
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Constructor
    
    /**
     * This is the default constructor
     * @param rMainModelPrt The model part to consider
     */
    
    InterfacePreprocessCondition(ModelPart& rMainModelPrt)
    :mrMainModelPart(rMainModelPrt)
    {
    }
    
    /// Destructor.
    virtual ~InterfacePreprocessCondition() = default;
    
    ///@}
    ///@name Operators
    ///@{
    
    ///@}
    ///@name Operations
    ///@{

    /**
     * Generate a new ModelPart containing only the interface. It will contain the conditions addressed in the call 
     * @param rOriginPart The original model part
     * @param rInterfacePart The interface model part
     * @param ThisParameters The configuration parameters
     */
    
    template<const unsigned int TDim>
    void GenerateInterfacePart(
            ModelPart& rOriginPart,
            ModelPart& rInterfacePart,
            Parameters ThisParameters =  Parameters(R"({})")
            );
    
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

    ModelPart&  mrMainModelPart;
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Creates a new condition with a giving name
     * @param prThisProperties The pointer to the element
     * @param rGeometry The  geometry considered
     * @param CondId The Id of the condition
     * @param rCondition The base condition
     */

    void CreateNewCondition(
            Properties::Pointer prThisProperties,
            GeometryType& rGeometry,
            const unsigned int CondId,
            Condition const& rCondition
            );
    
    /**
     * It prints the nodes and conditions in the interface, gives an error otherwise there are not
     * @param NodesCounter Number of nodes in the interface
     * @param CondCounter Number of conditions in the interface
     */

    void PrintNodesAndConditions(
            const int NodesCounter,
            const int CondCounter
            );
    
    /**
     * It reorders the Ids of the conditions
     * @return cond_id: The Id from the last condition
     */

    unsigned int ReorderConditions();
    
    /**
     * This method creates the conditions for the edges
     * @param rInterfacePart The model part of the interface
     * @param prThisProperties The properties of the base element
     * @param EdgeGeometry Geometry considered
     * @param SimplestGeometry If consider or not the simplest geometry
     * @param CondCounter The counter of conditions
     * @param CondId The condition id
     */

    inline void GenerateEdgeCondition(
        ModelPart& rInterfacePart,
        Properties::Pointer prThisProperties,
        GeometryType& EdgeGeometry,
        const bool SimplestGeometry,
        unsigned int& CondCounter,
        unsigned int& CondId
        );
    
    /**
     * This method creates the conditions for the faces
     * @param rInterfacePart The model part of the interface
     * @param rThisProperties The properties of the base element
     * @param FaceGeometry Geometry considered
     * @param SimplestGeometry If consider or not the simplest geometry
     * @param CondCounter The counter of conditions
     * @param CondId The condition id
     */

    inline void GenerateFaceCondition(
        ModelPart& rInterfacePart,
        Properties::Pointer prThisProperties,
        GeometryType& FaceGeometry,
        const bool SimplestGeometry,
        unsigned int& CondCounter,
        unsigned int& CondId
        );
    
    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class InterfacePreprocessCondition
}
#endif  /* KRATOS_INTERFACE_PREPROCESS_CONDITION_H_INCLUDED defined */
