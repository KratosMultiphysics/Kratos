// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
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
#include "contact_structural_mechanics_application_variables.h"
/* Geometries */
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

// TODO: Add parallellization!!!

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
    virtual ~InterfacePreprocessCondition() {}
    
    ///@}
    ///@name Operators
    ///@{
    
    ///@}
    ///@name Operations
    ///@{

    /**
     * Generate a new ModelPart containing only the interface. It will contain the conditions addressed in the call 
     * @param rOriginPart: The original model part
     * @param ConditionName: Name of the condition to be created
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
            )
    {
        KRATOS_TRY;
    
        Condition const & rCondition = KratosComponents<Condition>::Get(ConditionName); 
        Condition::Pointer pCond = Condition::Pointer(rCondition.Create(CondId, rGeometry, rpElem->pGetProperties()));
        rInterfacePart.AddCondition(pCond);
        if (ConditionName.find("Mortar") != std::string::npos)
        {
            pCond->GetValue(ELEMENT_POINTER) = rpElem;
            
            // We set the condition as master or slave (master by default)
            if (ConditionName.find("ALM") != std::string::npos || ConditionName.find("MeshTying") != std::string::npos)
            {
                bool IsSlave = true;
                for (unsigned int iNode = 0; iNode < pCond->GetGeometry().size(); iNode++)
                {
                    if (pCond->GetGeometry()[iNode].Is(SLAVE) == false)
                    {
                        IsSlave = false;
                    }
                }
                if (IsSlave == true)
                {
                    pCond->Set(SLAVE, true);
                }
                else
                {
                    pCond->Set(MASTER, true);
                }
            }
        }
        
        KRATOS_CATCH("");
    }
    
    /**
     * Generate a new ModelPart containing only the interface. It will contain only linear linear conditions 
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateLine2NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        const unsigned int dimension = mrMainModelPart.GetProcessInfo()[DOMAIN_SIZE];
                    
        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator iNode = rOriginPart.NodesBegin(); iNode != rOriginPart.NodesEnd(); iNode++)
        {
            if (iNode->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(iNode.base()) );
                NodesCounter ++;
            }
        }

        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator itCond = rOriginPart.ConditionsBegin(); itCond != rOriginPart.ConditionsEnd(); itCond++)
        {
            if (
                ((*itCond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[1].Is(INTERFACE) == true))
            {
                aux.push_back( *(itCond.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        if (dimension == 2) // By default, but someone can be interested in project values to a BEAM for example
        {
            GenerateLine2D2NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }
        else
        {
            GenerateLine3D2NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }

        KRATOS_CATCH("");
    }
    
    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadratic lines conditions 
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateLine3NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;
        
        const unsigned int dimension = mrMainModelPart.GetProcessInfo()[DOMAIN_SIZE];
        
        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator iNode = rOriginPart.NodesBegin(); iNode != rOriginPart.NodesEnd(); iNode++)
        {
            if (iNode->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(iNode.base()) );
                NodesCounter ++;
            }
        }
        
        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator itCond = rOriginPart.ConditionsBegin(); itCond != rOriginPart.ConditionsEnd(); itCond++)
        {
            if (
                ((*itCond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(itCond.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        if (dimension == 2) // By default, but someone can be interested in project values to a BEAM for example
        {
            GenerateLine2D3NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }
        else
        {
            GenerateLine3D3NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Generate a new ModelPart containing only the interface. It will contain only triangular conditions of 3 nodes, regardless of what was used while meshing
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateTriangle3NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;
        
        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator iNode = rOriginPart.NodesBegin(); iNode != rOriginPart.NodesEnd(); iNode++)
        {
            if (iNode->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(iNode.base()) );
                NodesCounter ++;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator itCond = rOriginPart.ConditionsBegin(); itCond != rOriginPart.ConditionsEnd(); itCond++)
        {
            if (
                ((*itCond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(itCond.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateTriangular3D3NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /**
     * Generate a new ModelPart containing only the interface. It will contain only triangular of 6 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateTriangle6NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator iNode = rOriginPart.NodesBegin(); iNode != rOriginPart.NodesEnd(); iNode++)
        {
            if (iNode->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(iNode.base()) );
                NodesCounter ++;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator itCond = rOriginPart.ConditionsBegin(); itCond != rOriginPart.ConditionsEnd(); itCond++)
        {
            if (
                ((*itCond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[5].Is(INTERFACE) == true))
            {
                aux.push_back( *(itCond.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateTriangular3D6NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 4 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral4NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator iNode = rOriginPart.NodesBegin(); iNode != rOriginPart.NodesEnd(); iNode++)
        {
            if (iNode->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(iNode.base()) );
                NodesCounter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator itCond = rOriginPart.ConditionsBegin(); itCond != rOriginPart.ConditionsEnd(); itCond++)
        {
            if (
                ((*itCond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[3].Is(INTERFACE) == true))
            {
                aux.push_back( *(itCond.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateQuadrilateral3D4NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 8 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral8NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator iNode = rOriginPart.NodesBegin(); iNode != rOriginPart.NodesEnd(); iNode++)
        {
            if (iNode->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(iNode.base()) );
                NodesCounter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator itCond = rOriginPart.ConditionsBegin(); itCond != rOriginPart.ConditionsEnd(); itCond++)
        {
            if (
                ((*itCond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[5].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[6].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[7].Is(INTERFACE) == true))
            {
                aux.push_back( *(itCond.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateQuadrilateral3D8NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 9 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral9NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator iNode = rOriginPart.NodesBegin(); iNode != rOriginPart.NodesEnd(); iNode++)
        {
            if (iNode->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(iNode.base()) );
                NodesCounter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator itCond = rOriginPart.ConditionsBegin(); itCond != rOriginPart.ConditionsEnd(); itCond++)
        {
            if (
                ((*itCond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[5].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[6].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[7].Is(INTERFACE) == true) &&
                ((*itCond).GetGeometry()[8].Is(INTERFACE) == true))
            {
                aux.push_back( *(itCond.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateQuadrilateral3D9NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }
    
    /**
     * Create a set of 2D linear lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine2D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition2D2N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType LinId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 2)
            {
                rLinConds.push_back( rCondition.Create(LinId++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 3)
            {
                Line2D2< Node<3> > Lin1(it->GetGeometry()(0), it->GetGeometry()(1));
                Line2D2< Node<3> > Lin2(it->GetGeometry()(1), it->GetGeometry()(2));

                rLinConds.push_back(rCondition.Create(LinId++, Lin1, it->pGetProperties()));
                rLinConds.push_back(rCondition.Create(LinId++, Lin2, it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }
    
    /**
     * Create a set of 3D linear lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine3D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D2N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType LinId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 2)
            {
                rLinConds.push_back( rCondition.Create(LinId++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 3)
            {
                Line3D2< Node<3> > Lin1(it->GetGeometry()(0), it->GetGeometry()(1));
                Line3D2< Node<3> > Lin2(it->GetGeometry()(1), it->GetGeometry()(2));

                rLinConds.push_back(rCondition.Create(LinId++, Lin1, it->pGetProperties()));
                rLinConds.push_back(rCondition.Create(LinId++, Lin2, it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }

    /**
     * Create a set of 2D quadratic lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine2D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition2D3N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType LinId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 3)
            {
                rLinConds.push_back( rCondition.Create(LinId++, it->GetGeometry(), it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }

    /**
     * Create a set of 3D quadratic lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D3N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType LinId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 3)
            {
                rLinConds.push_back( rCondition.Create(LinId++, it->GetGeometry(), it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }
    
    /**
     * Create a set of linear triangular conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rTriConds: The triangular conditions created
     */

    void GenerateTriangular3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType TriId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 3)
            {
                rTriConds.push_back( rCondition.Create(TriId++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 4)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(2));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(2), it->GetGeometry()(3), it->GetGeometry()(0));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 6)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri4(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 8)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(7));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(5), it->GetGeometry()(7));
                Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri4(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > Tri5(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri5, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri6, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 9)
            {
                Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(8));
                Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(8));
                Triangle3D3< Node<3> > Tri4(it->GetGeometry()(8), it->GetGeometry()(3), it->GetGeometry()(4));
                Triangle3D3< Node<3> > Tri5(it->GetGeometry()(8), it->GetGeometry()(4), it->GetGeometry()(5));
                Triangle3D3< Node<3> > Tri6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));
                Triangle3D3< Node<3> > Tri7(it->GetGeometry()(5), it->GetGeometry()(7), it->GetGeometry()(8));
                Triangle3D3< Node<3> > Tri8(it->GetGeometry()(0), it->GetGeometry()(8), it->GetGeometry()(7));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri5, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri6, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri7, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri8, it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear triangles " << std::endl;
            }
        }
    }

    /**
     * Create a set of quadratic triangular conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rTriConds: The triangular conditions created
     */

    void GenerateTriangular3D6NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D6N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType TriId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 6)
            {
                rTriConds.push_back( rCondition.Create(TriId++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 9)
            {
                Triangle3D6< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(8), it->GetGeometry()(7));
                Triangle3D6< Node<3> > Tri2(it->GetGeometry()(2), it->GetGeometry()(3), it->GetGeometry()(1), it->GetGeometry()(6), it->GetGeometry()(8), it->GetGeometry()(5));

                rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using quadratic triangles " << std::endl;
            }
        }
    }

    /**
     * Create a set of linear quadratic conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D4NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D4N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType QuadId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {

            if (it->GetGeometry().PointsNumber() == 4)
            {
                rQuadConds.push_back( rCondition.Create(QuadId++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 9)
            {
                Quadrilateral3D4< Node<3> > Quad1(it->GetGeometry()(0), it->GetGeometry()(4), it->GetGeometry()(8), it->GetGeometry()(7));
                Quadrilateral3D4< Node<3> > Quad2(it->GetGeometry()(4), it->GetGeometry()(1), it->GetGeometry()(5), it->GetGeometry()(8));
                Quadrilateral3D4< Node<3> > Quad3(it->GetGeometry()(8), it->GetGeometry()(5), it->GetGeometry()(2), it->GetGeometry()(6));
                Quadrilateral3D4< Node<3> > Quad4(it->GetGeometry()(7), it->GetGeometry()(8), it->GetGeometry()(6), it->GetGeometry()(3));

                rQuadConds.push_back(rCondition.Create(QuadId++, Quad1, it->pGetProperties()));
                rQuadConds.push_back(rCondition.Create(QuadId++, Quad2, it->pGetProperties()));
                rQuadConds.push_back(rCondition.Create(QuadId++, Quad3, it->pGetProperties()));
                rQuadConds.push_back(rCondition.Create(QuadId++, Quad4, it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear quadrilaterals " << std::endl;
            }
        }
    }
    
    /**
     * Create a set of quadratic quadratic (8 nodes) conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D8NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D8N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType QuadId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 8)
            {
                rQuadConds.push_back( rCondition.Create(QuadId++, it->GetGeometry(), it->pGetProperties()));
            }
            else
	    {
                KRATOS_ERROR << "The geometry can not be divided using quadratic (8 nodes) quadrilaterals " << std::endl;
	    }
        }
    }
    
    /**
     * Create a set of quadratic quadratic (9 nodes) conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D9NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D9N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType QuadId = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {

            if (it->GetGeometry().PointsNumber() == 9)
            {
                rQuadConds.push_back( rCondition.Create(QuadId++, it->GetGeometry(), it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using quadratic (9 nodes) quadrilaterals " << std::endl;
            }
        }
    }
    
    /**
     * It prints the nodes and conditions in the interface, gives an error otherwise there are not
     * @param NodesCounter: Number of nodes in the interface
     * @return CondCounter: Number of conditions in the interface
     */

    void PrintNodesAndConditions(
            const int NodesCounter,
            const int CondCounter
            )
    {
        std::cout << "    " << NodesCounter << " nodes ";
        std::cout << "and " << CondCounter <<  " conditions found." << std::endl;

        // Check that we actually found something
        if( NodesCounter == 0)
        {
            KRATOS_ERROR << "No interface nodes found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true." << std::endl;
        }
        if( CondCounter == 0)
        {
            KRATOS_ERROR << "No interface conditions found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true and that the contact surfaces have been assigned conditions." << std::endl;
        }
    }
    
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

///@name Explicit Specializations
///@{

    template<>
    void InterfacePreprocessCondition::GenerateInterfacePart<2>(
        ModelPart& rOriginPart,
        ModelPart& rInterfacePart,
        Parameters ThisParameters
        )
    {
        KRATOS_TRY;
        
        Parameters DefaultParameters = Parameters(R"(
        {
            "condition_name"                       : "", 
            "final_string"                         : "", 
            "simplify_geometry"                    : false
        })" );
        
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        const std::string ConditionName = ThisParameters["condition_name"].GetString();
        const std::string FinalString = ThisParameters["final_string"].GetString();
        const bool SimplestGeometry = ThisParameters["simplify_geometry"].GetBool();
        
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        unsigned int CondCounter = 0;
        unsigned int CondId = 0;
        
        // We reorder the conditions
        ConditionsArrayType& pCondition = mrMainModelPart.Conditions();
        auto numConditions = pCondition.end() - pCondition.begin();
        
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCondition = pCondition.begin() + i;
            CondId += 1;
            itCondition->SetId(CondId);
        }
        
        // Generate Conditions from original the edges that can be considered interface
        for (ModelPart::ElementsContainerType::const_iterator itElem = rOriginPart.ElementsBegin(); itElem != rOriginPart.ElementsEnd(); itElem++)
        {
            for (unsigned int iEdge = 0; iEdge < (*itElem).GetGeometry().EdgesNumber(); iEdge++)
            {
                unsigned int count = 0;
                const unsigned int NumberOfPoints = (*itElem).GetGeometry().Edges()[iEdge].PointsNumber();
                for (unsigned int iNode = 0; iNode < NumberOfPoints; iNode++)
                {
                    if ((*itElem).GetGeometry().Edges()[iEdge][iNode].IsDefined(INTERFACE) == true)  
                    {
                        if ((*itElem).GetGeometry().Edges()[iEdge][iNode].Is(INTERFACE) == true)  
                        {
                            count++;
                        }
                    }
                }
                
                if (count == NumberOfPoints)
                {
                    std::string EdgeConditionName = ConditionName;
                    if (NumberOfPoints == 2)
                    {
                        EdgeConditionName.append("Condition2D2N");
                        EdgeConditionName.append(FinalString);
                        
                        CondId += 1; // NOTE: To paralellize be careful with this ID
                        CreateNewCondition(rInterfacePart, *(itElem.base()), (*itElem).GetGeometry().Edges()[iEdge], CondId, EdgeConditionName);
                        CondCounter ++;
                    }
                    else
                    {                            
                        if (SimplestGeometry == false)
                        {
                            EdgeConditionName.append("Condition2D3N"); 
                            EdgeConditionName.append(FinalString); 
                            
                            CondId += 1; 
                            CreateNewCondition(rInterfacePart, *(itElem.base()), (*itElem).GetGeometry().Edges()[iEdge], CondId, EdgeConditionName);
                            CondCounter ++;
                        }
                        else
                        {
                            EdgeConditionName.append("Condition2D2N"); 
                            EdgeConditionName.append(FinalString); 
                            
                            CondId += 1;
                            Line2D2< Node<3> > Lin1((*itElem).GetGeometry().Edges()[iEdge](0), (*itElem).GetGeometry().Edges()[iEdge](1));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Lin1, CondId, EdgeConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Line2D2< Node<3> > Lin2((*itElem).GetGeometry().Edges()[iEdge](1), (*itElem).GetGeometry().Edges()[iEdge](2));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Lin2, CondId, EdgeConditionName);
                            CondCounter ++;
                        }
                    }
                }
            }
        }
       
      
        // NOTE: Reorder ID if parallellization
      
        PrintNodesAndConditions(numNodes, CondCounter);
      
        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>
    void InterfacePreprocessCondition::GenerateInterfacePart<3>(
        ModelPart& rOriginPart,
        ModelPart& rInterfacePart,
        Parameters ThisParameters
        )
    {
        KRATOS_TRY;
        
        Parameters DefaultParameters = Parameters(R"(
        {
            "condition_name"                       : "", 
            "final_string"                         : "", 
            "simplify_geometry"                    : false
        })" );
        
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        const std::string ConditionName = ThisParameters["condition_name"].GetString();
        const std::string FinalString = ThisParameters["final_string"].GetString();
        const bool SimplestGeometry = ThisParameters["simplify_geometry"].GetBool();
        
        NodesArrayType& pNode = mrMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        unsigned int CondCounter = 0;
        unsigned int CondId = 0;
        
        // We reorder the conditions
        ConditionsArrayType& pCondition = mrMainModelPart.Conditions();
        auto numConditions = pCondition.end() - pCondition.begin();
        
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCondition = pCondition.begin() + i;
            CondId += 1;
            itCondition->SetId(CondId);
        }
        
        // Generate Conditions from original the faces that can be considered interface
        for (ModelPart::ElementsContainerType::const_iterator itElem = rOriginPart.ElementsBegin(); itElem != rOriginPart.ElementsEnd(); itElem++)
        {                
            for (unsigned int iFace = 0; iFace < (*itElem).GetGeometry().FacesNumber(); iFace++)
            {
                unsigned int count = 0;
                const unsigned int NumberOfPoints = (*itElem).GetGeometry().Faces()[iFace].PointsNumber();
                for (unsigned int iNode = 0; iNode < NumberOfPoints; iNode++)
                {
                    if ((*itElem).GetGeometry().Faces()[iFace][iNode].IsDefined(INTERFACE) == true)  
                    {
                        if ((*itElem).GetGeometry().Faces()[iFace][iNode].Is(INTERFACE) == true)  
                        {
                            count++;
                        }
                    }
                }
                
                if (count == NumberOfPoints)
                {
                    std::string FaceConditionName = ConditionName;
                    if (NumberOfPoints == 3)
                    {
                        FaceConditionName.append("Condition3D3N");
                        FaceConditionName.append(FinalString);
                        
                        CondId += 1;
                        CreateNewCondition(rInterfacePart, *(itElem.base()), (*itElem).GetGeometry().Faces()[iFace], CondId, FaceConditionName);
                        CondCounter ++;
                    }
                    else if (NumberOfPoints == 4)
                    {
                        if (SimplestGeometry == false)
                        {
                            FaceConditionName.append("Condition3D4N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            CreateNewCondition(rInterfacePart, *(itElem.base()), (*itElem).GetGeometry().Faces()[iFace], CondId, FaceConditionName);
                            CondCounter ++;
                        }
                        else
                        {
                            FaceConditionName.append("Condition3D3N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri1((*itElem).GetGeometry().Faces()[iFace](0), (*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](2));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri1, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri2((*itElem).GetGeometry().Faces()[iFace](2), (*itElem).GetGeometry().Faces()[iFace](3), (*itElem).GetGeometry().Faces()[iFace](0));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri2, CondId, FaceConditionName);
                            CondCounter ++;
                        }
                    }
                    else if (NumberOfPoints == 6)
                    {
                            if (SimplestGeometry == false)
                        {
                            FaceConditionName.append("Condition3D6N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            CreateNewCondition(rInterfacePart, *(itElem.base()), (*itElem).GetGeometry().Faces()[iFace], CondId, FaceConditionName);
                            CondCounter ++;
                        }
                        else
                        {
                            FaceConditionName.append("Condition3D3N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri1((*itElem).GetGeometry().Faces()[iFace](0), (*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](5));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri1, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri2((*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](2), (*itElem).GetGeometry().Faces()[iFace](3));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri2, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri3((*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](3), (*itElem).GetGeometry().Faces()[iFace](5));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri3, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri4((*itElem).GetGeometry().Faces()[iFace](3), (*itElem).GetGeometry().Faces()[iFace](4), (*itElem).GetGeometry().Faces()[iFace](5));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri4, CondId, FaceConditionName);
                            CondCounter ++;
                        }
                    }
                    else if (NumberOfPoints == 8)
                    {
                        if (SimplestGeometry == false)
                        {
                            FaceConditionName.append("Condition3D8N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            CreateNewCondition(rInterfacePart, *(itElem.base()), (*itElem).GetGeometry().Faces()[iFace], CondId, FaceConditionName);
                            CondCounter ++;
                        }
                        else
                        {
                            FaceConditionName.append("Condition3D3N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri1((*itElem).GetGeometry().Faces()[iFace](0), (*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](7));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri1, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri2((*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](5), (*itElem).GetGeometry().Faces()[iFace](7));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri2, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri3((*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](3), (*itElem).GetGeometry().Faces()[iFace](5));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri3, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri4((*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](2), (*itElem).GetGeometry().Faces()[iFace](3));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri4, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri5((*itElem).GetGeometry().Faces()[iFace](3), (*itElem).GetGeometry().Faces()[iFace](4), (*itElem).GetGeometry().Faces()[iFace](5));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri5, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri6((*itElem).GetGeometry().Faces()[iFace](5), (*itElem).GetGeometry().Faces()[iFace](6), (*itElem).GetGeometry().Faces()[iFace](7));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri6, CondId, FaceConditionName);
                            CondCounter ++;
                        }
                    }
                    else // Assuming it will not be a very weird geometry
                    {
                        if (SimplestGeometry == false)
                        {
                            FaceConditionName.append("Condition3D4N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            CreateNewCondition(rInterfacePart, *(itElem.base()), (*itElem).GetGeometry().Faces()[iFace], CondId, FaceConditionName);
                            CondCounter ++;
                        }
                        else
                        {
                            FaceConditionName.append("Condition3D3N");
                            FaceConditionName.append(FinalString);
                            
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri1((*itElem).GetGeometry().Faces()[iFace](0), (*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](8));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri1, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri2((*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](2), (*itElem).GetGeometry().Faces()[iFace](3));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri2, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri3((*itElem).GetGeometry().Faces()[iFace](1), (*itElem).GetGeometry().Faces()[iFace](3), (*itElem).GetGeometry().Faces()[iFace](8));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri3, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri4((*itElem).GetGeometry().Faces()[iFace](8), (*itElem).GetGeometry().Faces()[iFace](3), (*itElem).GetGeometry().Faces()[iFace](4));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri4, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri5((*itElem).GetGeometry().Faces()[iFace](8), (*itElem).GetGeometry().Faces()[iFace](4), (*itElem).GetGeometry().Faces()[iFace](5));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri5, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri6((*itElem).GetGeometry().Faces()[iFace](5), (*itElem).GetGeometry().Faces()[iFace](6), (*itElem).GetGeometry().Faces()[iFace](7));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri6, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri7((*itElem).GetGeometry().Faces()[iFace](5), (*itElem).GetGeometry().Faces()[iFace](7), (*itElem).GetGeometry().Faces()[iFace](8));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri7, CondId, FaceConditionName);
                            CondCounter ++;
                            CondId += 1;
                            Triangle3D3< Node<3> > Tri8((*itElem).GetGeometry().Faces()[iFace](0), (*itElem).GetGeometry().Faces()[iFace](8), (*itElem).GetGeometry().Faces()[iFace](7));
                            CreateNewCondition(rInterfacePart, *(itElem.base()), Tri8, CondId, FaceConditionName);
                            CondCounter ++;
                        }
                    }
                }
            }
        }
      
        // NOTE: Reorder ID if parallellization
      
        PrintNodesAndConditions(numNodes, CondCounter);
      
        KRATOS_CATCH("");
    }

}

#endif  /* KRATOS_INTERFACE_PREPROCESS_CONDITION_H_INCLUDED defined */
