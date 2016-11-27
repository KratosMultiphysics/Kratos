// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
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
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "contact_structural_mechanics_application_variables.h"

// TODO: Add parallellization!!!
// TODO: Check geometry creation

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
    
    ///@}
    ///@name Life Cycle
    ///@{
    
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
    
    void GenerateInterfacePart(
            ModelPart& rOriginPart,
            ModelPart& rInterfacePart,
            std::string ConditionName,
            unsigned int CondId,
            std::string final_string
            )
    {
        KRATOS_TRY;
        
        const unsigned int dimension = rOriginPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();   
        
        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rInterfacePart.NodesBegin(); node_it != rInterfacePart.NodesEnd(); node_it++)
        {
            NodesCounter++;
        }
        
        unsigned int CondCounter = 0;
        
//         const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
//         OpenMPUtils::PartitionVector ElementPartition;
//         OpenMPUtils::DivideInPartitions(rOriginPart.Elements().size(), NumThreads, ElementPartition);
// 
//         const unsigned int nelem = static_cast<int>( rOriginPart.Elements().size() );
//         ElementsArrayType::iterator ElemBegin = rOriginPart.Elements().begin();
        
        if (dimension == 2)
        {
            // Generate Conditions from original the edges that can be considered interface
            for (ModelPart::ElementsContainerType::const_iterator elem_it = rOriginPart.ElementsBegin(); elem_it != rOriginPart.ElementsEnd(); elem_it++)
            {
//             #pragma omp parallel for
//             for(unsigned int i = 0;  i < nelem; i++)
//             {
//                 ElementsArrayType::iterator elem_it = ElemBegin + i; 
                
                for (unsigned int it_edge = 0; it_edge < (*elem_it).GetGeometry().EdgesNumber(); it_edge++)
                {
                    unsigned int count = 0;
                    unsigned int number_points = (*elem_it).GetGeometry().Edges()[it_edge].PointsNumber();
                    for (unsigned int node_it = 0; node_it < number_points; node_it++)
                    {
                        if ((*elem_it).GetGeometry().Edges()[it_edge][node_it].IsDefined(INTERFACE) == true)  
                        {
                            if ((*elem_it).GetGeometry().Edges()[it_edge][node_it].Is(INTERFACE) == true)  
                            {
                                count++;
                            }
                        }
                    }
                    
                    if (count == number_points)
                    {
                        CondId += 1; // NOTE: Para paralelizar cuidado con esta ID
                        std::string EdgeConditionName = ConditionName;
                        if (number_points == 2)
                        {
                            EdgeConditionName.append("Condition2D2N");
                            EdgeConditionName.append(final_string);
                        }
                        else
                        {
                            EdgeConditionName.append("Condition2D3N"); 
                            EdgeConditionName.append(final_string); 
                        }
                        
                        Condition const & rCondition = KratosComponents<Condition>::Get(EdgeConditionName);
                        Condition::Pointer pCond = Condition::Pointer(rCondition.Create(CondId, (*elem_it).GetGeometry().Edges()[it_edge], (*elem_it).pGetProperties()));
                        rInterfacePart.AddCondition(pCond);
                        if (ConditionName.find("Mortar") != std::string::npos)
                        {
                             Element::Pointer & pElem = const_cast<Condition &>(rCondition).GetValue(ELEMENT_POINTER);
                             pElem = *(elem_it.base());
                             // KRATOS_WATCH(pElem->Id());
                        }
                        CondCounter ++;
                    }
                }
            }
        }
        else
        {
            // Generate Conditions from original the faces that can be considered interface
            for (ModelPart::ElementsContainerType::const_iterator elem_it = rOriginPart.ElementsBegin(); elem_it != rOriginPart.ElementsEnd(); elem_it++)
            {
//             #pragma omp parallel for
//             for(unsigned int i = 0;  i < nelem; i++)
//             {
//                 ElementsArrayType::iterator elem_it = ElemBegin + i; 
                
                for (unsigned int it_face = 0; it_face < (*elem_it).GetGeometry().FacesNumber(); it_face++)
                {
                    unsigned int count = 0;
                    unsigned int number_points = (*elem_it).GetGeometry().Faces()[it_face].PointsNumber();
                    for (unsigned int node_it = 0; node_it < number_points; node_it++)
                    {
                        if ((*elem_it).GetGeometry().Faces()[it_face][node_it].IsDefined(INTERFACE) == true)  
                        {
                            if ((*elem_it).GetGeometry().Faces()[it_face][node_it].Is(INTERFACE) == true)  
                            {
                                count++;
                            }
                        }
                    }
                    
                    if (count == number_points)
                    {
                        CondId += 1;
                        std::string FaceConditionName = ConditionName;
                        if (number_points == 3)
                        {
                            FaceConditionName.append("Condition3D3N");
                            FaceConditionName.append(final_string);
                        }
                        else if (number_points == 4)
                        {
                            FaceConditionName.append("Condition3D4N");
                            FaceConditionName.append(final_string);
                        }
                        else if (number_points == 6)
                        {
                            FaceConditionName.append("Condition3D6N");
                            FaceConditionName.append(final_string);
                        }
                        else if (number_points == 8)
                        {
                            FaceConditionName.append("Condition3D8N");
                            FaceConditionName.append(final_string);
                        }
                        else // Assuming it will not be a very weird geometry
                        {
                            FaceConditionName.append("Condition3D9N");
                            FaceConditionName.append(final_string);
                        }
  
                        Condition const & rCondition = KratosComponents<Condition>::Get(FaceConditionName); 
                        Condition::Pointer pCond = Condition::Pointer(rCondition.Create(CondId, (*elem_it).GetGeometry().Faces()[it_face], (*elem_it).pGetProperties()));
                        rInterfacePart.AddCondition(pCond);
                        if (ConditionName.find("Mortar") != std::string::npos)
                        {
                             Element::Pointer & pElem = const_cast<Condition &>(rCondition).GetValue(ELEMENT_POINTER);
                             pElem = *(elem_it.base());
                             // KRATOS_WATCH(pElem->Id());
                        }
                        CondCounter ++;
                    }
                }
            }
        }
      
        // NOTE: Reorder ID if parallellization
      
        PrintNodesAndConditions(NodesCounter, CondCounter);
      
        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Generate a new ModelPart containing only the interface. It will contain only linear linear conditions 
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateLine2NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            std::string ConditionName
            )
    {
        KRATOS_TRY;

        const unsigned int dimension = rOriginPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
                    
        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
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

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadratic lines conditions 
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateLine3NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            std::string ConditionName
            )
    {
        KRATOS_TRY;
        
        const unsigned int dimension = rOriginPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
        
        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }
        
        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
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
            std::string ConditionName
            )
    {
        KRATOS_TRY;
        
        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateTriangular3D3NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Generate a new ModelPart containing only the interface. It will contain only triangular of 6 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateTriangle6NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[5].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateTriangular3D6NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 4 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral4NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[3].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateQuadrilateral3D4NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 8 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral8NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[5].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[6].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[7].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateQuadrilateral3D8NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Generate a new ModelPart containing only the interface. It will contain only quadrialterals of 9 nodes
     * @param rOriginPart: The original model part
     * @return InterfacePart: The interface model part
     */

    void GenerateQuadrilateral9NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int NodesCounter = 0;
        for (ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin(); node_it != rOriginPart.NodesEnd(); node_it++)
        {
            if (node_it->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(node_it.base()) );
                NodesCounter ++;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int CondCounter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin(); cond_it != rOriginPart.ConditionsEnd(); cond_it++)
        {
            if (
                ((*cond_it).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[1].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[2].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[3].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[4].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[5].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[6].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[7].Is(INTERFACE) == true) &&
                ((*cond_it).GetGeometry()[8].Is(INTERFACE) == true))
            {
                aux.push_back( *(cond_it.base()) );
                CondCounter ++;
            }
        }

        PrintNodesAndConditions(NodesCounter, CondCounter);

        GenerateQuadrilateral3D9NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of 2D linear lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine2D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using linear lines ", "");
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of 3D linear lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine3D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using linear lines ", "");
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of 2D quadratic lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine2D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using linear lines ", "");
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of 3D quadratic lines conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rLinConds: The linear conditions created
     */

    void GenerateLine3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using linear lines ", "");
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of linear triangular conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rTriConds: The triangular conditions created
     */

    void GenerateTriangular3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using linear triangles ", "");
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of quadratic triangular conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rTriConds: The triangular conditions created
     */

    void GenerateTriangular3D6NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using quadratic triangles ", "");
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Create a set of linear quadratic conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D4NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using linear quadrilaterals ", "");
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Create a set of quadratic quadratic (8 nodes) conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D8NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            std::string ConditionName
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
	      KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using quadratic (8 nodes) quadrilaterals ", "");
	    }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Create a set of quadratic quadratic (9 nodes) conditions from a generic condition set
     * @param rOriginConds: The original conditions
     * @return rQuadConds: The triangular conditions created
     */

    void GenerateQuadrilateral3D9NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            std::string ConditionName
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
                KRATOS_THROW_ERROR( std::logic_error, "The geometry can not be divided using quadratic (9 nodes) quadrilaterals ", "");
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
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
            KRATOS_THROW_ERROR(std::invalid_argument,"No interface nodes found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true.","");
        }
        if( CondCounter == 0)
        {
            KRATOS_THROW_ERROR(std::invalid_argument,"No interface conditions found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true and that the contact surfaces have been assigned conditions.","");
        }
    }
    
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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

#endif
