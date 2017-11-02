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
    
    typedef Point                                  PointType;
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
     * @param rModelPart: The model part to consider
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
     * @param rOriginPart: The original model part
     * @param condition_name: Name of the condition to be created
     * @return InterfacePart: The interface model part
     */
    
    template< const unsigned int TDim>
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
     * @return InterfacePart: The interface model part
     */

    void CreateNewCondition(
            Element::Pointer rpElem,
            Geometry<Node<3> > & rGeometry,
            const unsigned int CondId,
            const std::string& ConditionName,
            const bool IsMortar
            );
    
    /**
     * It prints the nodes and conditions in the interface, gives an error otherwise there are not
     * @param NodesCounter: Number of nodes in the interface
     * @return CondCounter: Number of conditions in the interface
     */

    void PrintNodesAndConditions(
            const int& NodesCounter,
            const int& CondCounter
            );
    
    /**
     * It reorders the Ids of the conditions
     * @return cond_id: The Id from the last condition
     */

    unsigned int ReorderConditions();
    
    /**
     * This method creates the conditions for the edges
     * @param rInterfacePart: The model part of the interface
     * @param rpElem: Pointer to the element
     * @param EdgeGeometry: Geometry considered
     * @param ConditionName: The name of the condition
     * @param FinalString: The last part added to the name condition
     * @param SimplestGeometry: If consider or not the simplest geometry
     * @param CondCounter: The counter of conditions
     * @param CondId: The condition id
     */

    inline void GenerateEdgeCondition(
        ModelPart& rInterfacePart,
        Element::Pointer rpElem,
        GeometryType& EdgeGeometry,
        const std::string& ConditionName,
        const std::string& FinalString,
        const bool& SimplestGeometry,
        unsigned int& CondCounter,
        unsigned int& CondId,
        const bool IsMortar
        );
    
    /**
     * This method creates the conditions for the faces
     * @param rInterfacePart: The model part of the interface
     * @param rpElem: Pointer to the element
     * @param FaceGeometry: Geometry considered
     * @param ConditionName: The name of the condition
     * @param FinalString: The last part added to the name condition
     * @param SimplestGeometry: If consider or not the simplest geometry
     * @param CondCounter: The counter of conditions
     * @param CondId: The condition id
     */

    inline void GenerateFaceCondition(
        ModelPart& rInterfacePart,
        Element::Pointer rpElem,
        GeometryType& FaceGeometry,
        const std::string& ConditionName,
        const std::string& FinalString,
        const bool& SimplestGeometry,
        unsigned int& CondCounter,
        unsigned int& CondId,
        const bool IsMortar
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
