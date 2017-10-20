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
            )
    {
        KRATOS_TRY;
    
        Condition const & rCondition = KratosComponents<Condition>::Get(ConditionName); 
        Condition::Pointer p_cond = Condition::Pointer(rCondition.Create(CondId, rGeometry, rpElem->pGetProperties()));
        rInterfacePart.AddCondition(p_cond);
        if (ConditionName.find("Mortar") != std::string::npos)
        {
            p_cond->GetValue(ELEMENT_POINTER) = rpElem;
            
            // We set the condition as master or slave (master by default)
            if (ConditionName.find("ALM") != std::string::npos || ConditionName.find("MeshTying") != std::string::npos)
            {
                bool IsSlave = true;
                for (unsigned int i_node = 0; i_node < p_cond->GetGeometry().size(); i_node++)
                {
                    if (p_cond->GetGeometry()[i_node].Is(SLAVE) == false)
                    {
                        IsSlave = false;
                    }
                }
                if (IsSlave == true)
                {
                    p_cond->Set(SLAVE, true);
                }
                else
                {
                    p_cond->Set(MASTER, true);
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
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); i_node++)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                nodes_counter ++;
            }
        }

        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); it_cond++)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                cond_counter ++;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

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
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); i_node++)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                nodes_counter ++;
            }
        }
        
        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); it_cond++)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                cond_counter ++;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

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
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); i_node++)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                nodes_counter ++;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); it_cond++)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                cond_counter ++;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

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
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); i_node++)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                nodes_counter ++;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); it_cond++)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[5].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                cond_counter ++;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

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
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); i_node++)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                nodes_counter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); it_cond++)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[3].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                cond_counter ++;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

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
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); i_node++)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                nodes_counter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); it_cond++)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[5].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[6].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[7].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                cond_counter ++;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

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
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); i_node++)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                nodes_counter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); it_cond++)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[5].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[6].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[7].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[8].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                cond_counter ++;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

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
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 2)
            {
                rLinConds.push_back( rCondition.Create(lin_id++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 3)
            {
                Line2D2< Node<3> > lin_1(it->GetGeometry()(0), it->GetGeometry()(1));
                Line2D2< Node<3> > lin_2(it->GetGeometry()(1), it->GetGeometry()(2));

                rLinConds.push_back(rCondition.Create(lin_id++, lin_1, it->pGetProperties()));
                rLinConds.push_back(rCondition.Create(lin_id++, lin_2, it->pGetProperties()));
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
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 2)
            {
                rLinConds.push_back( rCondition.Create(lin_id++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 3)
            {
                Line3D2< Node<3> > lin_1(it->GetGeometry()(0), it->GetGeometry()(1));
                Line3D2< Node<3> > lin_2(it->GetGeometry()(1), it->GetGeometry()(2));

                rLinConds.push_back(rCondition.Create(lin_id++, lin_1, it->pGetProperties()));
                rLinConds.push_back(rCondition.Create(lin_id++, lin_2, it->pGetProperties()));
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
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 3)
            {
                rLinConds.push_back( rCondition.Create(lin_id++, it->GetGeometry(), it->pGetProperties()));
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
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 3)
            {
                rLinConds.push_back( rCondition.Create(lin_id++, it->GetGeometry(), it->pGetProperties()));
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
        Condition::IndexType tri_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 3)
            {
                rTriConds.push_back( rCondition.Create(tri_id++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 4)
            {
                Triangle3D3< Node<3> > tri_1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(2));
                Triangle3D3< Node<3> > tri_2(it->GetGeometry()(2), it->GetGeometry()(3), it->GetGeometry()(0));

                rTriConds.push_back(rCondition.Create(tri_id++, tri_1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_2, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 6)
            {
                Triangle3D3< Node<3> > tri_1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(5));
                Triangle3D3< Node<3> > tri_2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > tri_3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                Triangle3D3< Node<3> > tri_4(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));

                rTriConds.push_back(rCondition.Create(tri_id++, tri_1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_4, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 8)
            {
                Triangle3D3< Node<3> > tri_1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(7));
                Triangle3D3< Node<3> > tri_2(it->GetGeometry()(1), it->GetGeometry()(5), it->GetGeometry()(7));
                Triangle3D3< Node<3> > tri_3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                Triangle3D3< Node<3> > tri_4(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > tri_5(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));
                Triangle3D3< Node<3> > tri_6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));

                rTriConds.push_back(rCondition.Create(tri_id++, tri_1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_4, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_5, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_6, it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 9)
            {
                Triangle3D3< Node<3> > tri_1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(8));
                Triangle3D3< Node<3> > tri_2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                Triangle3D3< Node<3> > tri_3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(8));
                Triangle3D3< Node<3> > tri_4(it->GetGeometry()(8), it->GetGeometry()(3), it->GetGeometry()(4));
                Triangle3D3< Node<3> > tri_5(it->GetGeometry()(8), it->GetGeometry()(4), it->GetGeometry()(5));
                Triangle3D3< Node<3> > tri_6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));
                Triangle3D3< Node<3> > tri_7(it->GetGeometry()(5), it->GetGeometry()(7), it->GetGeometry()(8));
                Triangle3D3< Node<3> > tri_8(it->GetGeometry()(0), it->GetGeometry()(8), it->GetGeometry()(7));

                rTriConds.push_back(rCondition.Create(tri_id++, tri_1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_2, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_3, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_4, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_5, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_6, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_7, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_8, it->pGetProperties()));
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
        Condition::IndexType tri_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 6)
            {
                rTriConds.push_back( rCondition.Create(tri_id++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 9)
            {
                Triangle3D6< Node<3> > tri_1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(8), it->GetGeometry()(7));
                Triangle3D6< Node<3> > tri_2(it->GetGeometry()(2), it->GetGeometry()(3), it->GetGeometry()(1), it->GetGeometry()(6), it->GetGeometry()(8), it->GetGeometry()(5));

                rTriConds.push_back(rCondition.Create(tri_id++, tri_1, it->pGetProperties()));
                rTriConds.push_back(rCondition.Create(tri_id++, tri_2, it->pGetProperties()));
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
        Condition::IndexType quad_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {

            if (it->GetGeometry().PointsNumber() == 4)
            {
                rQuadConds.push_back( rCondition.Create(quad_id++, it->GetGeometry(), it->pGetProperties()));
            }
            else if (it->GetGeometry().PointsNumber() == 9)
            {
                Quadrilateral3D4< Node<3> > quad_1(it->GetGeometry()(0), it->GetGeometry()(4), it->GetGeometry()(8), it->GetGeometry()(7));
                Quadrilateral3D4< Node<3> > quad_2(it->GetGeometry()(4), it->GetGeometry()(1), it->GetGeometry()(5), it->GetGeometry()(8));
                Quadrilateral3D4< Node<3> > quad_3(it->GetGeometry()(8), it->GetGeometry()(5), it->GetGeometry()(2), it->GetGeometry()(6));
                Quadrilateral3D4< Node<3> > quad_4(it->GetGeometry()(7), it->GetGeometry()(8), it->GetGeometry()(6), it->GetGeometry()(3));

                rQuadConds.push_back(rCondition.Create(quad_id++, quad_1, it->pGetProperties()));
                rQuadConds.push_back(rCondition.Create(quad_id++, quad_2, it->pGetProperties()));
                rQuadConds.push_back(rCondition.Create(quad_id++, quad_3, it->pGetProperties()));
                rQuadConds.push_back(rCondition.Create(quad_id++, quad_4, it->pGetProperties()));
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
        Condition::IndexType quad_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {
            if (it->GetGeometry().PointsNumber() == 8)
            {
                rQuadConds.push_back( rCondition.Create(quad_id++, it->GetGeometry(), it->pGetProperties()));
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
        Condition::IndexType quad_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); it++)
        {

            if (it->GetGeometry().PointsNumber() == 9)
            {
                rQuadConds.push_back( rCondition.Create(quad_id++, it->GetGeometry(), it->pGetProperties()));
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
    
    /**
     * It reorders the Ids of the conditions
     * @return cond_id: The Id from the last condition
     */

    unsigned int ReorderConditions()
    {
        // We reorder the conditions
        ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
        const unsigned int num_conditions = static_cast<int>(conditions_array.size());
        
        for(unsigned int i = 0; i < num_conditions; i++) 
        {
            auto it_condition = conditions_array.begin() + i;
            it_condition->SetId(i + 1);
        }
        
        return num_conditions;
    }
    
    /**
     * It method helps to reduce code duplication due two the use for flat and solids elements
     * @return rInterfacePart: The model part of the interface
     * @return rpElem: Pointer to the element
     * @return FaceGeometry: Geometry considered
     * @return ConditionName: The name of the condition
     * @return FinalString: The last part added to the name condition
     * @return SimplestGeometry: If consider or not the simplest geometry
     * @return CondCounter: The counter of conditions
     * @return CondId: The condition id
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
        )
    {
        unsigned int count = 0;
        const unsigned int number_of_points = FaceGeometry.PointsNumber();
        for (unsigned int i_node = 0; i_node < number_of_points; i_node++)
        {
            if (FaceGeometry[i_node].IsDefined(INTERFACE) == true)  
            {
                if (FaceGeometry[i_node].Is(INTERFACE) == true)  
                {
                    count++;
                }
            }
        }
        
        if (count == number_of_points)
        {
            std::string face_condition_name = ConditionName;
            CondId += 1;
            if (number_of_points == 3)
            {
                face_condition_name.append("Condition3D3N");
                face_condition_name.append(FinalString);
                
                CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                CondCounter ++;
            }
            else if (number_of_points == 4)
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D4N");
                    face_condition_name.append(FinalString);
                
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    CondCounter ++;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);
                    
                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(2));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(2), FaceGeometry(3), FaceGeometry(0));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    CondCounter ++;
                }
            }
            else if (number_of_points == 6)
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D6N");
                    face_condition_name.append(FinalString);
                    
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    CondCounter ++;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);
                    
                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_3, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_4(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_4, CondId, face_condition_name);
                    CondCounter ++;
                }
            }
            else if (number_of_points == 8)
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D8N");
                    face_condition_name.append(FinalString);
                    
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    CondCounter ++;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);

                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(5), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_3, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_4(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                    CreateNewCondition(rInterfacePart, rpElem, tri_4, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_5(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_5, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_6, CondId, face_condition_name);
                    CondCounter ++;
                }
            }
            else // Assuming it will not be a very weird geometry
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D4N");
                    face_condition_name.append(FinalString);
                    
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    CondCounter ++;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);

                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(8));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(8));
                    CreateNewCondition(rInterfacePart, rpElem, tri_3, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_4(FaceGeometry(8), FaceGeometry(3), FaceGeometry(4));
                    CreateNewCondition(rInterfacePart, rpElem, tri_4, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_5(FaceGeometry(8), FaceGeometry(4), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_5, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_6, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_7(FaceGeometry(5), FaceGeometry(7), FaceGeometry(8));
                    CreateNewCondition(rInterfacePart, rpElem, tri_7, CondId, face_condition_name);
                    CondCounter ++;
                    CondId += 1;
                    Triangle3D3< Node<3> > tri_8(FaceGeometry(0), FaceGeometry(8), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_8, CondId, face_condition_name);
                    CondCounter ++;
                }
            }
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
        
        Parameters default_parameters = Parameters(R"(
        {
            "condition_name"                       : "", 
            "final_string"                         : "", 
            "simplify_geometry"                    : false
        })" );
        
        ThisParameters.ValidateAndAssignDefaults(default_parameters);
        
        const std::string condition_name = ThisParameters["condition_name"].GetString();
        const std::string final_string = ThisParameters["final_string"].GetString();
        const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
        
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const unsigned int num_nodes = static_cast<int>(nodes_array.size());
        
        unsigned int cond_counter = 0;
        
        // We reorder the conditions
        unsigned int cond_id = ReorderConditions();
        
        // Generate Conditions from original the edges that can be considered interface
        for (ModelPart::ElementsContainerType::const_iterator it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); it_elem++)
        {
            for (unsigned int i_edge = 0; i_edge < (*it_elem).GetGeometry().EdgesNumber(); i_edge++)
            {
                unsigned int count = 0;
                const unsigned int number_of_points = (*it_elem).GetGeometry().Edges()[i_edge].PointsNumber();
                for (unsigned int i_node = 0; i_node < number_of_points; i_node++)
                {
                    if ((*it_elem).GetGeometry().Edges()[i_edge][i_node].IsDefined(INTERFACE) == true)  
                    {
                        if ((*it_elem).GetGeometry().Edges()[i_edge][i_node].Is(INTERFACE) == true)  
                        {
                            count++;
                        }
                    }
                }
                
                if (count == number_of_points)
                {
                    std::string Edgecondition_name = condition_name;
                    cond_id += 1; // NOTE: To paralellize be careful with this ID
                    if (number_of_points == 2)
                    {
                        Edgecondition_name.append("Condition2D2N");
                        Edgecondition_name.append(final_string);
                        
                        CreateNewCondition(rInterfacePart, *(it_elem.base()), (*it_elem).GetGeometry().Edges()[i_edge], cond_id, Edgecondition_name);
                        cond_counter ++;
                    }
                    else
                    {                            
                        if (simplest_geometry == false)
                        {
                            Edgecondition_name.append("Condition2D3N"); 
                            Edgecondition_name.append(final_string); 
                            
                            CreateNewCondition(rInterfacePart, *(it_elem.base()), (*it_elem).GetGeometry().Edges()[i_edge], cond_id, Edgecondition_name);
                            cond_counter ++;
                        }
                        else
                        {
                            Edgecondition_name.append("Condition2D2N"); 
                            Edgecondition_name.append(final_string); 

                            Line2D2< Node<3> > lin_1((*it_elem).GetGeometry().Edges()[i_edge](0), (*it_elem).GetGeometry().Edges()[i_edge](1));
                            CreateNewCondition(rInterfacePart, *(it_elem.base()), lin_1, cond_id, Edgecondition_name);
                            cond_counter ++;
                            cond_id += 1;
                            Line2D2< Node<3> > lin_2((*it_elem).GetGeometry().Edges()[i_edge](1), (*it_elem).GetGeometry().Edges()[i_edge](2));
                            CreateNewCondition(rInterfacePart, *(it_elem.base()), lin_2, cond_id, Edgecondition_name);
                            cond_counter ++;
                        }
                    }
                }
            }
        }
       
      
        // NOTE: Reorder ID if parallellization
      
        PrintNodesAndConditions(num_nodes, cond_counter);
      
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
        
        Parameters default_parameters = Parameters(R"(
        {
            "condition_name"                       : "", 
            "final_string"                         : "", 
            "simplify_geometry"                    : false
        })" );
        
        ThisParameters.ValidateAndAssignDefaults(default_parameters);
        
        const std::string& condition_name = ThisParameters["condition_name"].GetString();
        const std::string& final_string = ThisParameters["final_string"].GetString();
        const bool& simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
        
        unsigned int cond_counter = 0;
        
        // We reorder the conditions
        unsigned int cond_id = ReorderConditions();

        // Generate Conditions from original the faces that can be considered interface
        for (ModelPart::ElementsContainerType::const_iterator it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); it_elem++)
        {                
            if ((*it_elem).GetGeometry().LocalSpaceDimension() == 3)
            {
                for (unsigned int i_face = 0; i_face < (*it_elem).GetGeometry().FacesNumber(); i_face++)
                {
                    GenerateFaceCondition(rInterfacePart, *(it_elem.base()), (*it_elem).GetGeometry().Faces()[i_face], condition_name, final_string, simplest_geometry, cond_counter, cond_id);
                }
            }
            else
            {
                GenerateFaceCondition(rInterfacePart, *(it_elem.base()), (*it_elem).GetGeometry(), condition_name, final_string, simplest_geometry, cond_counter, cond_id);
            }
        }
      
        // NOTE: Reorder ID if parallellization
      
        const unsigned int num_nodes = static_cast<int>(mrMainModelPart.Nodes().size());
        PrintNodesAndConditions(num_nodes, cond_counter);
      
        KRATOS_CATCH("");
    }
}
#endif  /* KRATOS_INTERFACE_PREPROCESS_CONDITION_H_INCLUDED defined */
