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
            ModelPart& rInterfacePart,
            Element::Pointer rpElem,
            Geometry<Node<3> > & rGeometry,
            const unsigned int CondId,
            const std::string ConditionName
            );
    
    /**
     * Generate a new ModelPart containing only the interface. It will contain only linear linear conditions 
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateLine2NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            );
    
    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadratic lines conditions 
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateLine3NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            );

    /**
     * Generate a new ModelPart containing only the interface. It will contain only triangular conditions of 3 nodes, regardless of what was used while meshing
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateTriangle3NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            );

    /**
     * Generate a new ModelPart containing only the interface. It will contain only triangular of 6 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateTriangle6NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            );

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 4 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral4NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            );

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 8 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral8NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            );

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 9 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral9NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            );
    
    /**
     * Create a set of 2D linear lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine2D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            );
    
    /**
     * Create a set of 3D linear lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine3D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            );

    /**
     * Create a set of 2D quadratic lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine2D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            );

    /**
     * Create a set of 3D quadratic lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            );
    
    /**
     * Create a set of linear triangular conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rTriConds: The triangular conditions created
     */

    void GenerateTriangular3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            const std::string ConditionName
            );

    /**
     * Create a set of quadratic triangular conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rTriConds: The triangular conditions created
     */

    void GenerateTriangular3D6NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            const std::string ConditionName
            );

    /**
     * Create a set of linear quadratic conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D4NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            );
    
    /**
     * Create a set of quadratic quadratic (8 nodes) conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D8NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            );
    
    /**
     * Create a set of quadratic quadratic (9 nodes) conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D9NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            );
    
    /**
     * It prints the nodes and conditions in the interface, gives an error otherwise there are not
     * @param NodesCounter: Number of nodes in the interface
     * @return CondCounter: Number of conditions in the interface
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
     * It method helps to reduce code duplication due two the use for flat and solids elements
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
